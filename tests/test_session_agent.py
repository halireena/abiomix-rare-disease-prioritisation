"""session_agent: the LLM driver refines HPO + proposes GATED PP4, deterministically re-folding. Fake LLM = hermetic."""
import duckdb
import pandas as pd
from acmg.repl import CaseSession
from acmg.session_agent import agent_review


def _session():
    con = duckdb.connect()
    free = pd.DataFrame([
        {"variant_key": "1-1-A-G", "gene": "AAA", "total_points": 3, "criteria": "PM2", "acmg_class": "VUS"},
        {"variant_key": "2-1-C-T", "gene": "BBB", "total_points": 3, "criteria": "PM2", "acmg_class": "VUS"},
    ])
    return CaseSession(con, "CASE0001", free, store=None)


def test_agent_awards_gated_pp4_and_refolds():
    s = _session()
    base = s.table.set_index("candidate_id").loc["1-1-A-G", "case_points"]
    # fake LLM: no phrases (skip grounding), one PP4_Moderate on the AAA candidate
    llm = lambda _p: '{"add_phrases": [], "pp4": [{"candidate_id": "1-1-A-G", "criterion": "PP4_Moderate", "why": "specific match"}]}'
    top = agent_review(s, "some phenotype", cr_index="", llm=llm, max_rounds=1)
    got = s.table.set_index("candidate_id").loc["1-1-A-G", "case_points"]
    assert base == 0 and got == 2                       # PP4_Moderate (+2) gated fact applied via propose->approve
    assert "PP4_Moderate" in s.table.set_index("candidate_id").loc["1-1-A-G", "case_criteria"]
    assert len(top) == 2


def test_agent_ignores_invalid_or_out_of_scope_proposals():
    s = _session()
    # a criterion the agent may NOT assert (PS2 family) + a non-existent candidate -> both ignored
    llm = lambda _p: '{"add_phrases": [], "pp4": [{"candidate_id": "9-9-Z-Z", "criterion": "PP4"}, {"candidate_id": "1-1-A-G", "criterion": "PS2"}]}'
    agent_review(s, "pheno", cr_index="", llm=llm, max_rounds=1)
    assert s.table.set_index("candidate_id").loc["1-1-A-G", "case_points"] == 0   # PS2 not in PP4/PP4_Moderate -> ignored
