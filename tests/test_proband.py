"""run_proband — end-to-end per-proband composition: phenotype-aware SNV+CNV top-N with IC-weighted rank. Hermetic."""
import duckdb
import pandas as pd
from acmg.proband import run_proband


def _con_with_phenotype():
    con = duckdb.connect()
    # gene_phenotype fixture: AAA is specific to the patient's term, TTN is pleiotropic (many terms -> low IC)
    gp = pd.DataFrame(
        [("AAA", "HP:0001250")] + [("TTN", f"HP:{i:07d}") for i in range(1, 40)] + [("TTN", "HP:0001250")],
        columns=["gene", "hpo_id"])
    con.register("_gp", gp)
    con.execute("CREATE TABLE gene_phenotype AS SELECT * FROM _gp")
    return con


def _snv_free():
    return pd.DataFrame([
        ("1-1-A-G", "AAA", 5, "PM2_Supporting"),
        ("2-1-C-T", "TTN", 5, "PM2_Supporting"),
    ], columns=["variant_key", "gene", "total_points", "criteria"]).assign(acmg_class="VUS")


def test_ic_weighted_rank_prefers_specific_gene():
    con = _con_with_phenotype()
    top, s = run_proband(con, "CASE0001", _snv_free(), hpo=["HP:0001250"], store=None)
    # AAA (specific -> high IC) ranks above TTN (pleiotropic -> low IC) within the same VUS tier
    assert top.iloc[0]["candidate_id"] == "1-1-A-G"
    assert s.session_id  # an open REPL session is returned


def test_pathogenic_cnv_makes_the_top_n():
    con = _con_with_phenotype()
    cnv = pd.DataFrame([{"cnv_key": "X-1-2-DEL", "gene": "DMD", "cnv_class": "Pathogenic", "svtype": "DEL", "cn": 0}])
    top, s = run_proband(con, "CASE0002", _snv_free(), cnv=cnv, hpo=["HP:0001250"], store=None, top_n=10)
    assert top.iloc[0]["candidate_id"] == "X-1-2-DEL" and top.iloc[0]["kind"] == "cnv"   # Pathogenic CNV wins
    assert "snv" in set(top["kind"]) and "cnv" in set(top["kind"])                        # both tracks present


def test_no_phenotype_table_still_runs():
    con = duckdb.connect()   # no gene_phenotype table -> phenotype_score falls back to 0, still ranks by class
    top, s = run_proband(con, "C", _snv_free(), hpo=["HP:0001250"], store=None)
    assert len(top) == 2
