"""Persistent case REPL (acmg.repl.CaseSession): SQL over the case table, HPO refinement -> deterministic
re-rank, gated case facts -> re-fold, every action persisted + replayable. Hermetic."""
import json
import duckdb
import pandas as pd
from acmg.repl import CaseSession


def _free():
    # all three band to VUS (same class tier) so phenotype/genotype break the tie in ranking
    return pd.DataFrame([
        ("1-1-A-G", "AAA", 5, "PM2_Supporting"),
        ("2-1-C-T", "BBB", 4, "PM2_Supporting"),
        ("3-1-G-A", "CCC", 3, "PM2_Supporting"),
    ], columns=["variant_key", "gene", "total_points", "criteria"]).assign(acmg_class="VUS")


# test rank_fn: +2 for a gene named in the HPO set (test treats hpo as gene symbols)
def _rank(df, hpo):
    return [2.0 if g in hpo else 0.0 for g in df["gene"]]


def test_session_persists_and_reranks(tmp_path):
    store = str(tmp_path / "sessions.jsonl")
    con = duckdb.connect()
    s = CaseSession(con, "CASE0001", _free(), store=store, rank_fn=_rank)
    # SQL over the case table works + is recorded
    n = s.sql("SELECT count(*) c FROM case_variants").iloc[0]["c"]
    assert n == 3
    # refine HPO -> re-rank: boost gene AAA above the others (within the same class tier)
    s.refine_hpo(add=["AAA"])
    assert s.top().iloc[0]["candidate_id"] == "1-1-A-G"      # AAA now ranks first (deterministic re-rank)
    # gated case fact: pending doesn't move points; approved does
    dg = s.propose_fact("2-1-C-T", "PS3", 4, evidence={"pmid": "1"})
    assert s.table.set_index("candidate_id").loc["2-1-C-T", "case_points"] == 0    # pending
    s.approve_fact(dg, "curator")
    assert s.table.set_index("candidate_id").loc["2-1-C-T", "case_points"] == 4    # approved -> counts
    # persistence: the ledger has one line per action, SQL-queryable
    lines = [json.loads(x) for x in open(store).read().splitlines()]
    actions = [r["action"] for r in lines]
    assert actions == ["open", "sql", "refine_hpo", "propose_fact", "approve_fact"]
    assert all(r["session_id"] == s.session_id and r["case_id"] == "CASE0001" for r in lines)
    assert lines[-1]["state"] != lines[0]["state"]           # state digest changed (hpo + approved fact)


def test_store_none_disables_persistence(tmp_path):
    con = duckdb.connect()
    s = CaseSession(con, "C", _free(), store=None, rank_fn=_rank)
    s.refine_hpo(add=["BBB"])           # no file written, no error
    assert not (tmp_path / "sessions.jsonl").exists()


def test_cnv_competes_in_the_repl(tmp_path):
    # a pathogenic CNV enters the SAME case table and outranks the SNVs (SNV-only top-10 would miss it)
    con = duckdb.connect()
    cnv = pd.DataFrame([{"cnv_key": "X-1-2-DEL", "gene": "DMD", "cnv_class": "Pathogenic", "svtype": "DEL", "cn": 0}])
    s = CaseSession(con, "CASE0002", _free(), cnv=cnv, store=None, rank_fn=_rank)
    top = s.top()
    assert top.iloc[0]["candidate_id"] == "X-1-2-DEL" and top.iloc[0]["kind"] == "cnv"
    assert set(s.sql("SELECT DISTINCT kind FROM case_variants")["kind"]) == {"snv", "cnv"}
