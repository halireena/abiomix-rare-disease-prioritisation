"""HPO extraction from clinical text — FastHPOCR candidates + (preferred) LLM-aided adjudication.

Two lanes, the RfastHPOCR (github.com/sounkou-bioinfo/RfastHPOCR) hybrid harness modes:

  1. `tool_only` — `extract()`: DETERMINISTIC FastHPOCR concept recognition (pure Python) + a REGEX
     negation/subject-scope layer. This is the fast, reproducible FLOOR, but the regex clause-splitter is
     brittle on real narratives: it over-excludes "personal AND family history of X" (drops the proband's own
     X) and under-scopes "the patient's father has Y" (leaks the relative's Y into observed).
  2. `candidates_model` — `extract_llm()` (PREFERRED for real cases): FastHPOCR proposes candidates, then an
     LLM (pi, default gpt-5.3-codex-spark) does keep/drop adjudication with `patient_context`
     (patient / family_history / negated / uncertain) + an auditable `evidence_span` + `short_reason`. The
     model — not a regex — resolves subject-scope and negation, which is where the deterministic layer fails.
     Same schema/prompt contract as RfastHPOCR::hpo_adjudicate_candidates. The model only adjudicates the
     supplied candidates (it does not invent HPO IDs); the decision is auditable and can be gated (acmg.agent).
"""
from __future__ import annotations
import hashlib
import json
import os
import re
import time
from dataclasses import dataclass, field

_NEG = re.compile(r"\b(no|not|without|denies|denied|absent|negative for|ruled out|non|resolved)\b", re.I)
_FAMILY = re.compile(r"\b(family history|fh of|mother|father|sister|brother|sibling|maternal|paternal|cousin|aunt|uncle|grandmother|grandfather)\b", re.I)
_CLAUSE_SPLIT = re.compile(r"[.,;]|\b(?:but|however|except|although|though)\b", re.I)


@dataclass
class HPOResult:
    observed: list[str] = field(default_factory=list)      # asserted, proband-intrinsic
    excluded: list[str] = field(default_factory=list)      # explicitly negated ("not ataxic")
    family_scope: list[str] = field(default_factory=list)  # belongs to a relative, not the proband
    decisions: list[dict] = field(default_factory=list)    # full per-candidate audit rows (patient_context / evidence_span / short_reason)


def _annotator(cr_index_file: str):
    from FastHPOCR.HPOAnnotator import HPOAnnotator
    return HPOAnnotator(cr_index_file)


def _terms_in(annotator, text: str) -> list[str]:
    """HP ids FastHPOCR finds in a text span (AnnotationObject.hpoUri, e.g. 'HP:0001250')."""
    return [a.hpoUri for a in annotator.annotate(text) if getattr(a, "hpoUri", None)]


def extract(text: str, cr_index_file: str) -> HPOResult:
    """Extract HPO terms from one clinical note with negation + subject-scope handling.

    cr_index_file: the FastHPOCR CR index (build once from hp.obo via FastHPOCR.IndexHPO)."""
    ann = _annotator(cr_index_file)
    res = HPOResult()
    seen = set()
    for clause in _CLAUSE_SPLIT.split(text or ""):
        if not clause or not clause.strip():
            continue
        negated = bool(_NEG.search(clause))
        family = bool(_FAMILY.search(clause))
        for hp in _terms_in(ann, clause):
            if hp in seen:
                continue
            seen.add(hp)
            if family:
                res.family_scope.append(hp)   # relative's phenotype, not the proband's
            elif negated:
                res.excluded.append(hp)        # explicitly absent
            else:
                res.observed.append(hp)
    return res


# --------------------------------------------------------------------------------------------------
# candidates_model: FastHPOCR candidates -> LLM keep/drop adjudication (RfastHPOCR contract)
# --------------------------------------------------------------------------------------------------
def _candidates(annotator, text: str, case_id: str) -> list[dict]:
    """FastHPOCR annotations -> adjudication rows (candidate_id, candidate_span, hpo_id)."""
    cands, seen = [], set()
    for i, a in enumerate(annotator.annotate(text or "")):
        hp = getattr(a, "hpoUri", None)
        if not hp:
            continue
        span = (getattr(a, "text", None) or getattr(a, "textSpan", None) or getattr(a, "matchedText", None))
        if not span:
            s, e = getattr(a, "startOffset", None), getattr(a, "endOffset", None)
            span = text[s:e] if (s is not None and e is not None) else ""
        key = (hp, (span or "").lower())
        if key in seen:
            continue
        seen.add(key)
        cands.append({"candidate_id": f"{case_id}:{i:04d}", "candidate_span": span or "", "hpo_id": hp})
    return cands


def adjudication_prompt(note: str, candidates: list[dict], case_id: str = "case") -> str:
    """The candidates_model prompt (RfastHPOCR contract): keep/drop each FastHPOCR candidate with an auditable
    patient_context, evidence_span and short_reason. No hidden chain-of-thought; do not invent HPO IDs."""
    return "\n".join([
        "You are adjudicating candidate Human Phenotype Ontology (HPO) extractions from de-identified clinical text.",
        "Return ONLY valid JSON (no markdown fences, no prose outside JSON, no hidden chain-of-thought).",
        f"Case id: {case_id}",
        "Decision rules:",
        "- Return exactly one decision for each supplied candidate_id.",
        "- decision='keep' ONLY if the phenotype is present FOR THE PATIENT (proband) in the note.",
        "- decision='drop' if the mention is negated, family-history-only, uncertain/not established, not about the "
        "patient, too generic, or duplicated by a more specific kept candidate.",
        "- 'personal and family history of X' means the patient HAS X -> keep, patient_context='patient'.",
        "- 'the patient's father/mother/sibling has Y' means Y belongs to a relative -> drop, patient_context='family_history'.",
        "- patient_context must be one of: patient, family_history, negated, uncertain.",
        "- evidence_span = the shortest exact quote from the note supporting the decision (include negation/family wording).",
        "- short_reason = one concise sentence. Do not invent new HPO IDs; only adjudicate the supplied candidates.",
        'Output: {"case_id":"...","decisions":[{"candidate_id":"","hpo_id":"HP:.......","decision":"keep|drop",'
        '"patient_context":"patient|family_history|negated|uncertain","evidence_span":"","short_reason":""}]}',
        "",
        "Source note:",
        note or "",
        "",
        "Candidate rows as JSON:",
        json.dumps(candidates, ensure_ascii=False),
    ])


def _parse_decisions(raw) -> list[dict]:
    """Parse the model's JSON reply (str / dict / {text|content|reply}) -> the decisions list."""
    if isinstance(raw, dict) and "decisions" not in raw:
        raw = raw.get("text") or raw.get("content") or raw.get("reply") or ""
    if isinstance(raw, dict):
        return raw.get("decisions", [])
    s = re.sub(r"^```(?:json)?|```$", "", str(raw).strip(), flags=re.M).strip()
    m = re.search(r"\{.*\}", s, re.S)
    return (json.loads(m.group(0)).get("decisions", []) if m else [])


def _select(note: str, candidates: list[dict], llm, case_id: str) -> HPOResult:
    """LLM keep/drop selection over `candidates` -> HPOResult (the shared 'model selection at the end' stage)."""
    res = HPOResult()
    if not candidates:
        return res
    decisions = _parse_decisions(llm(adjudication_prompt(note, candidates, case_id)))
    res.decisions = decisions
    by_id = {c["candidate_id"]: c["hpo_id"] for c in candidates}
    seen = set()
    for d in decisions:
        hp = d.get("hpo_id") or by_id.get(d.get("candidate_id"))
        if not hp or hp in seen:
            continue
        seen.add(hp)
        ctx, dec = d.get("patient_context"), d.get("decision")
        if dec == "keep" and ctx == "patient":
            res.observed.append(hp)
        elif ctx == "family_history":
            res.family_scope.append(hp)
        elif ctx == "negated":
            res.excluded.append(hp)
        # 'uncertain' / generic drops -> asserted in neither list (not present for the proband)
    return res


def augment_phrases(note: str, llm) -> list[str]:
    """LLM AUGMENTATION (recall): list the patient's phenotype PHRASES verbatim — including paraphrases the
    FastHPOCR dictionary misses. The model returns PHRASES ONLY; grounding to HP IDs is done deterministically by
    FastHPOCR (_ground), so the model cannot hallucinate ontology codes."""
    prompt = "\n".join([
        "From the clinical note below, list the distinct PHENOTYPE phrases present FOR THE PATIENT (the proband) —",
        "the exact words describing an abnormal clinical feature. Exclude negated findings, a relative's phenotypes,",
        "and normal findings. Return ONLY JSON: {\"phrases\": [\"...\"]}. Do NOT output HPO codes; just the phrases.",
        "", "Note:", note or "",
    ])
    try:
        raw = llm(prompt)
        obj = raw if isinstance(raw, dict) and "phrases" in raw else json.loads(re.search(r"\{.*\}", str(raw), re.S).group(0))
        return [p for p in obj.get("phrases", []) if isinstance(p, str) and p.strip()]
    except Exception:
        return []


def _ground(annotator, phrases: list[str], case_id: str, offset: int = 5000) -> list[dict]:
    """Ontology tool: ground LLM phrases to HP IDs via FastHPOCR -> candidate rows. Only valid HP IDs survive."""
    cands, seen = [], set()
    for i, ph in enumerate(phrases):
        for a in annotator.annotate(ph):
            hp = getattr(a, "hpoUri", None)
            if hp and (hp, ph.lower()) not in seen:
                seen.add((hp, ph.lower()))
                cands.append({"candidate_id": f"{case_id}:aug{offset + i:04d}", "candidate_span": ph, "hpo_id": hp})
    return cands


HPO_MODES = ("tool_only", "candidates_model", "augment_select", "model_only")
HPO_STORE = ".cache/hpo_runs.jsonl"


def _persist(store: str | None, record: dict) -> None:
    """Append an IMMUTABLE HPO adjudication record (append-only JSONL, SQL-queryable via DuckDB read_json_auto).
    The model judgment (candidates, per-candidate decisions with evidence_span/short_reason, augmented phrases,
    model, note digest) is durable + auditable — same discipline as acmg.agent's ledger. store=None disables."""
    if not store:
        return
    os.makedirs(os.path.dirname(store) or ".", exist_ok=True)
    with open(store, "a") as fh:
        fh.write(json.dumps(record, ensure_ascii=False) + "\n")


def extract_hpo(text: str, cr_index_file: str, *, mode: str = "augment_select", llm=None,
                model: str = "gpt-5.3-codex-spark", case_id: str = "case", store: "str | None" = HPO_STORE) -> HPOResult:
    """HPO extraction with a selectable HARNESS MODE (RfastHPOCR ablations) — different configs per task:
      - 'augment_select'   : LLM AUGMENT phrases + FastHPOCR direct candidates, both GROUNDED by FastHPOCR,
                             -> LLM selection (recall + precision; the full augment+ground+select). DEFAULT.
      - 'candidates_model' : FastHPOCR candidates -> LLM keep/drop selection (no augmentation lane).
      - 'model_only'       : LLM augment phrases -> FastHPOCR grounding -> LLM selection (no FastHPOCR-direct lane).
      - 'tool_only'        : FastHPOCR + regex negation/scope (deterministic FLOOR; no LLM).
    LLM augmentation returns PHRASES only; FastHPOCR grounds them to HP ids (codes are never hallucinated).
    Every run is PERSISTED to `store` (append-only JSONL; set store=None to disable). `llm`: callable(prompt)->
    str|dict; default = pi.run (spark model)."""
    if mode not in HPO_MODES:
        raise ValueError(f"mode must be one of {HPO_MODES}")
    note = text or ""
    base = {"ts": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "event": "hpo_extraction",
            "case_id": case_id, "mode": mode, "model": (None if mode == "tool_only" else model),
            "note_sha256": hashlib.sha256(note.encode("utf-8")).hexdigest(), "note_len": len(note)}
    if mode == "tool_only":
        r = extract(text, cr_index_file)
        _persist(store, {**base, "observed": r.observed, "excluded": r.excluded, "family_scope": r.family_scope})
        return r
    if llm is None:
        from . import pi
        llm = lambda p: pi.run(p, model=model, mode="text")
    ann = _annotator(cr_index_file)
    direct = _candidates(ann, note, case_id) if mode in ("candidates_model", "augment_select") else []
    phrases = augment_phrases(note, llm) if mode in ("augment_select", "model_only") else []
    augmented = _ground(ann, phrases, case_id)
    merged, seen = [], set()
    for c in direct + augmented:            # dedup on hpo_id, keep first (direct before augmented)
        if c["hpo_id"] not in seen:
            seen.add(c["hpo_id"]); merged.append(c)
    res = _select(note, merged, llm, case_id)
    _persist(store, {**base, "augmented_phrases": phrases, "candidates": merged, "decisions": res.decisions,
                     "observed": res.observed, "excluded": res.excluded, "family_scope": res.family_scope})
    return res


def extract_llm(text: str, cr_index_file: str, *, llm=None, model: str = "gpt-5.3-codex-spark",
                case_id: str = "case") -> HPOResult:
    """candidates_model: FastHPOCR candidates -> LLM keep/drop adjudication. Thin alias for
    extract_hpo(mode='candidates_model'); kept for direct use."""
    return extract_hpo(text, cr_index_file, mode="candidates_model", llm=llm, model=model, case_id=case_id)
