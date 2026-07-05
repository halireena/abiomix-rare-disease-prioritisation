"""LLM driver over acmg.repl.CaseSession — the AGENTIC layer kept OUT of the deterministic batch on purpose.

Given an open case session + the referral phenotype text, the model (BLIND to any answer) proposes, per round:
  1. additional phenotype PHRASES it reads in the narrative -> grounded to HP ids by FastHPOCR (acmg.hpo.extract,
     never model-invented codes) -> session.refine_hpo -> deterministic re-rank;
  2. PP4 phenotype-specificity on the candidate(s) whose GENE's known disease matches the phenotype -> gated
     session.propose_fact -> approve -> deterministic re-fold.
Bounded rounds; stops when nothing changes. The model sees ONLY the phenotype + the pipeline's own candidate table
(gene/class/scores) — no truth, no solved/unsolved. The re-rank/re-fold stays deterministic; the model only supplies
grounded HPO + gated judgment. This is the fair, no-cheating "agent-augmented" arm to compare against the baseline.

NOTE (eval mode): PP4 proposals are AUTO-APPROVED here so the arm is fully automated. In production the same
propose->digest->approve gate routes approval to a human/critic; auto-approve is the "agent decides" setting, logged.
"""
from __future__ import annotations
import json
from . import hpo as hpo_mod

_SYSTEM = ("You are a clinical-genomics reviewer prioritising variants for a rare-disease re-analysis. You are shown "
           "the referral phenotype and the pipeline's current candidate shortlist. You DO NOT know the diagnosis. "
           "Reason ONLY from the phenotype and public gene-disease knowledge. Be conservative and specific.")


def _prompt(clinical_text: str, hpo_terms, top_df) -> str:
    cands = top_df[["candidate_id", "gene", "acmg_class", "phenotype_score", "genotype_fit"]].head(15).to_dict("records")
    return (
        f"Referral phenotype:\n{clinical_text}\n\n"
        f"HPO already captured: {sorted(hpo_terms)}\n\n"
        f"Candidate shortlist (current ranking, best first):\n{json.dumps(cands, default=str)}\n\n"
        "Do TWO things, from the phenotype ONLY:\n"
        "1. add_phrases: short phenotype noun-phrases present in the narrative that may be missing from the HPO list "
        "(e.g. 'lower limb weakness'). Do NOT invent HP codes.\n"
        "2. pp4: for candidates whose gene's known disease SPECIFICALLY matches this phenotype, award PP4 (+1) for a "
        "plausible match or PP4_Moderate (+2) only for a highly specific match. Omit if nothing matches.\n"
        'Return STRICT JSON only: {"add_phrases": ["..."], "pp4": [{"candidate_id":"...","criterion":"PP4","why":"..."}]}'
    )


def agent_review(session, clinical_text: str, *, cr_index: str, llm=None, model: str = "gpt-5.3-codex-spark",
                 max_rounds: int = 2, approver: str = "agent"):
    """Drive `session` with the LLM for up to max_rounds; return the agent-refined top-10. `llm`: a prompt->text
    callable (defaults to acmg.pi). Deterministic given the model's outputs (grounding + fold are deterministic)."""
    from .case import CASE_POINTS
    if llm is None:
        from . import pi
        llm = lambda p: pi.run(p, model=model, mode="text", system=_SYSTEM)

    for _ in range(max_rounds):
        raw = llm(_prompt(clinical_text, session.hpo, session.top(15)))
        try:
            obj = json.loads(raw[raw.find("{"): raw.rfind("}") + 1])
        except Exception:
            break
        changed = False

        phrases = [p for p in (obj.get("add_phrases") or []) if str(p).strip()]
        if phrases:
            grounded = set(hpo_mod.extract(" . ".join(phrases), cr_index).observed)  # phrases -> HP ids (grounded)
            new = grounded - session.hpo
            if new:
                session.refine_hpo(add=new)
                changed = True

        valid = set(session.top(50)["candidate_id"])
        for p in (obj.get("pp4") or []):
            vk, crit = p.get("candidate_id"), p.get("criterion", "PP4")
            if vk in valid and crit in ("PP4", "PP4_Moderate"):
                dg = session.propose_fact(vk, crit, CASE_POINTS[crit], evidence={"why": p.get("why"), "by": "agent_review"})
                session.approve_fact(dg, approver)                                  # eval mode: auto-approve, logged
                changed = True

        if not changed:
            break
    return session.top(10)
