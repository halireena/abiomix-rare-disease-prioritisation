"""ClinGen gene_curation ingest: parse the versioned Gene-Disease Validity + Dosage feeds (fixtures, no
network) into the kernel's gene_curation, and confirm the shape the kernel consumes. Run:
`PYTHONPATH=. python3 tests/test_clingen.py`."""
import os
import duckdb
from acmg import clingen

FIX = os.path.join(os.path.dirname(__file__), "fixtures")


def _built():
    con = duckdb.connect()
    clingen.load_gene_curation(con, os.path.join(FIX, "clingen_gv.csv"), os.path.join(FIX, "clingen_dosage.csv"))
    return con


def test_schema_matches_kernel_gene_curation():
    con = _built()
    cols = [c[0] for c in con.execute("DESCRIBE gene_curation").fetchall()]
    assert cols == ["gene", "disease_id", "mode_of_inheritance", "hi_score", "ts_score", "lof_mechanism",
                    "gene_disease_validity", "source"]


def test_preamble_and_separator_rows_filtered():
    con = _built()
    # only the 2 real genes survive; the '++++' separator and preamble never become rows
    assert con.execute("SELECT count(*) FROM gene_curation").fetchone()[0] == 2
    assert set(con.execute("SELECT gene FROM gene_curation").fetchall()) == {("TP53",), ("XYZAR",)}


def test_haploinsufficiency_and_lof_mechanism():
    con = _built()
    tp53 = con.execute("SELECT hi_score, lof_mechanism, mode_of_inheritance, gene_disease_validity "
                       "FROM gene_curation WHERE gene='TP53'").fetchone()
    assert tp53 == (3, True, "AD", "Definitive")  # sufficient HI -> score 3, LoF mechanism true
    ar = con.execute("SELECT hi_score, lof_mechanism, mode_of_inheritance FROM gene_curation WHERE gene='XYZAR'").fetchone()
    assert ar == (30, True, "AR")  # AR-associated dosage (30) -> LoF IS a (recessive) disease mechanism; MOI normalised


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print("ok", name)
    print("all passed")
