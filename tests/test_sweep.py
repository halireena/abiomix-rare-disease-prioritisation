"""Multiplex sweep gate: shared-homozygous candidates are kept when ClinVar-P/LP OR cohort-rare, and dropped
when common-and-not-pathogenic. This is the 'don't discard a common-but-pathogenic ancestry variant' rule
(the G6PD-rescue). Run: `PYTHONPATH=. python3 tests/test_sweep.py`."""
import os
import duckdb
import pandas as pd
from acmg import family
from scripts.run_bundle import sweep_multiplex


def _con(tmp_cohort):
    # two affected brothers share three homozygous variants
    rows = []
    for case in ("S1", "S2"):
        for key in ("1-100-A-G", "1-200-C-T", "1-300-G-A"):
            chrom, pos, ref, alt = key.split("-")
            rows.append(dict(case_id=case, family_id="FAM_M", member_id=case,
                             proband_status="unclear_multi_affected", relationship="affected_son",
                             family_size=2, variant_kind="snv", variant_key=key, chrom=chrom, pos=int(pos),
                             ref=ref, alt=alt, gene=None, gt="1/1", gq=40, dp=30, phase_set=None))
    con = duckdb.connect()
    con.register("_df", pd.DataFrame(rows)); con.execute("CREATE TABLE sample_call AS SELECT * FROM _df")
    con.unregister("_df")
    family.build_carriage(con)
    # synthetic exact ClinVar: only 1-100 is Pathogenic
    con.execute("""CREATE TABLE clinvar AS SELECT * FROM (VALUES
        ('1', 100, 'A', 'G', 'Pathogenic', 2)) t(chrom, pos, ref, alt, clnsig, review_stars)""")
    # synthetic cohort freq: 1-100 common (0.90), 1-200 rare (0.01), 1-300 common (0.90)
    pd.DataFrame([
        dict(chr="1", pos=100, ref="A", alt="G", cohort_freq=0.90),
        dict(chr="1", pos=200, ref="C", alt="T", cohort_freq=0.01),
        dict(chr="1", pos=300, ref="G", alt="A", cohort_freq=0.90),
    ]).to_parquet(tmp_cohort)
    return con


def test_gate_keeps_rare_or_pathogenic_drops_common_benign(tmp_path=None):
    tmp = str((tmp_path or __import__("pathlib").Path("/tmp")) / "_cohort_freq_test.parquet")
    con = _con(tmp)
    try:
        out = sweep_multiplex(con, rare=0.05, cohort_freq=tmp)
        keys = set(out.variant_key)
        assert "1-100-A-G" in keys, "common-but-ClinVar-Pathogenic must be kept (the G6PD-rescue)"
        assert "1-200-C-T" in keys, "rare variant must be kept"
        assert "1-300-G-A" not in keys, "common + not-pathogenic must be dropped"
        assert bool(out.loc[out.variant_key == "1-100-A-G", "clinvar_plp"].iloc[0]) is True
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


if __name__ == "__main__":
    import pathlib
    test_gate_keeps_rare_or_pathogenic_drops_common_benign(pathlib.Path("/tmp"))
    print("all passed")
