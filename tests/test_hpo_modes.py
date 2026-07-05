"""HPO harness modes (RfastHPOCR candidates_model / augment_select) — hermetic, no network.

A fake `llm` returns canned JSON so the LLM-aided lanes are testable offline. Covers the two things the regex
layer got wrong: 'personal and family history of X' (proband HAS X -> observed) and "patient's father has Y"
(relative -> family_scope)."""
import json
import os
import pytest
from acmg import hpo

IDX = os.path.join(os.path.dirname(__file__), "..", ".cache", "hp.index")
pytestmark = pytest.mark.skipif(not os.path.exists(IDX), reason="FastHPOCR hp.index not built")


def _fake(decisions, phrases=None):
    def llm(prompt):
        if "PHENOTYPE phrases" in prompt:
            return json.dumps({"phrases": phrases or []})
        return json.dumps({"case_id": "c", "decisions": decisions})
    return llm


def test_personal_and_family_kept_for_proband():
    # regex over-excludes this to family; the model keeps it as the proband's
    note = "Personal and family history of recurrent hip dislocation."
    llm = _fake([{"candidate_id": "case:0000", "hpo_id": "HP:0002827", "decision": "keep", "patient_context": "patient"}])
    r = hpo.extract_hpo(note, IDX, mode="candidates_model", llm=llm, store=None)
    assert "HP:0002827" in r.observed and r.family_scope == []


def test_relative_phenotype_scoped_to_family():
    note = "The patient's father has scoliosis."
    llm = _fake([{"candidate_id": "case:0000", "hpo_id": "HP:0002650", "decision": "drop", "patient_context": "family_history"}])
    r = hpo.extract_hpo(note, IDX, mode="candidates_model", llm=llm, store=None)
    assert "HP:0002650" in r.family_scope and r.observed == []


def test_negation_scoped_to_excluded():
    note = "No seizures were observed."
    llm = _fake([{"candidate_id": "case:0000", "hpo_id": "HP:0001250", "decision": "drop", "patient_context": "negated"}])
    r = hpo.extract_hpo(note, IDX, mode="candidates_model", llm=llm, store=None)
    assert "HP:0001250" in r.excluded and r.observed == []


def test_augment_grounds_via_fasthpocr_not_hallucinated_ids():
    # model augments a PHRASE; FastHPOCR grounds it to a real HP id; then selection keeps it
    note = "The child has developmental delay."
    llm = _fake(
        decisions=[{"candidate_id": "case:0000", "hpo_id": "HP:0001263", "decision": "keep", "patient_context": "patient"}],
        phrases=["developmental delay"],
    )
    r = hpo.extract_hpo(note, IDX, mode="augment_select", llm=llm, store=None)
    assert any(h.startswith("HP:") for h in r.observed)


def test_bad_mode_raises():
    with pytest.raises(ValueError):
        hpo.extract_hpo("x", IDX, mode="not_a_mode", llm=_fake([]), store=None)


def test_parse_decisions_strips_fences():
    fenced = "```json\n{\"decisions\":[{\"candidate_id\":\"x\",\"hpo_id\":\"HP:0000001\",\"decision\":\"drop\",\"patient_context\":\"negated\"}]}\n```"
    assert len(hpo._parse_decisions(fenced)) == 1


def test_run_is_persisted_with_audit_fields(tmp_path):
    import json
    store = str(tmp_path / "hpo_runs.jsonl")
    note = "The patient has seizures."
    llm = _fake([{"candidate_id": "case:0000", "hpo_id": "HP:0001250", "decision": "keep",
                  "patient_context": "patient", "evidence_span": "seizures", "short_reason": "stated for the patient"}])
    hpo.extract_hpo(note, IDX, mode="candidates_model", llm=llm, store=store)
    rec = json.loads(open(store).read().splitlines()[-1])
    assert rec["event"] == "hpo_extraction" and rec["mode"] == "candidates_model"
    assert rec["note_sha256"] and rec["decisions"] and rec["observed"] == ["HP:0001250"]
