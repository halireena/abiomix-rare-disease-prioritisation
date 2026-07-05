"""Family segregation tests: multi-affected shared-homozygous (affected_shared_candidates), the QC gate,
min_affected, and X-linked hemizygous inclusion. Run: `PYTHONPATH=. python3 tests/test_family.py`."""
import duckdb
import pandas as pd
from acmg import family


def _sample_call(rows):
    """rows: (case_id, relationship, proband_status, variant_key, gt, gq, dp). Builds a minimal canonical
    sample_call for one multiplex family (FAM_M)."""
    recs = []
    for case_id, rel, status, key, gt, gq, dp in rows:
        chrom, pos, ref, alt = key.split("-")
        recs.append(dict(case_id=case_id, family_id="FAM_M", member_id=case_id, proband_status=status,
                         relationship=rel, family_size=3, variant_kind="snv", variant_key=key,
                         chrom=chrom, pos=int(pos), ref=ref, alt=alt, gene=None, gt=gt, gq=gq, dp=dp,
                         phase_set=None))
    con = duckdb.connect()
    con.register("_df", pd.DataFrame(recs))
    con.execute("CREATE TABLE sample_call AS SELECT * FROM _df")
    con.unregister("_df")
    family.build_carriage(con)
    return con


# Two affected brothers + unaffected mother.
ROWS = [
    ("S1", "affected_son", "unclear_multi_affected", "1-100-A-G", "1/1", 40, 30),  # shared hom (autosomal) -> hit
    ("S2", "affected_son", "unclear_multi_affected", "1-100-A-G", "1/1", 40, 30),
    ("MO", "mother", "unclear", "1-100-A-G", "0/1", 40, 30),
    ("S1", "affected_son", "unclear_multi_affected", "1-200-C-T", "1/1", 40, 30),  # hom in ONE affected -> miss
    ("S2", "affected_son", "unclear_multi_affected", "1-200-C-T", "0/1", 40, 30),
    ("S1", "affected_son", "unclear_multi_affected", "1-300-G-A", "0/1", 40, 30),  # het in both -> miss (not hom)
    ("S2", "affected_son", "unclear_multi_affected", "1-300-G-A", "0/1", 40, 30),
    ("S1", "affected_son", "unclear_multi_affected", "X-500-C-T", "1/1", 40, 30),  # shared hom on X (hemizygous) -> hit
    ("S2", "affected_son", "unclear_multi_affected", "X-500-C-T", "1/1", 40, 30),
    ("S1", "affected_son", "unclear_multi_affected", "1-600-T-C", "1/1", 5, 30),   # shared hom but low GQ -> QC drop
    ("S2", "affected_son", "unclear_multi_affected", "1-600-T-C", "1/1", 5, 30),
]


def test_shared_homozygous_surfaces_only_shared_hom():
    con = _sample_call(ROWS)
    got = family.affected_shared_candidates(con).df()
    keys = set(got.variant_key)
    assert "1-100-A-G" in keys, "autosomal shared hom missed"
    assert "X-500-C-T" in keys, "X-linked shared hom (hemizygous males) missed"
    assert "1-200-C-T" not in keys, "hom in only one affected should not surface"
    assert "1-300-G-A" not in keys, "het-in-both is not homozygous, should not surface"
    assert "1-600-T-C" not in keys, "low-GQ call must be dropped by the QC gate"
    assert int(got.loc[got.variant_key == "1-100-A-G", "n_affected_hom"].iloc[0]) == 2


def test_min_affected_threshold():
    con = _sample_call(ROWS)
    # only 2 affected exist; requiring 3 yields nothing
    assert len(family.affected_shared_candidates(con, min_affected=3).df()) == 0


def test_affected_label_vocabulary_is_parameterised():
    con = _sample_call(ROWS)
    # a label not in the affected set (and status not matched) => those members are not 'affected' => no hits
    out = family.affected_shared_candidates(con, affected_labels=("affected_daughter",),
                                            affected_status=("proband",)).df()
    assert len(out) == 0


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name)
    print("all passed")
