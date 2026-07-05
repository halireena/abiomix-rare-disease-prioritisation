"""CNV track tests. The Riggs banding is tested standalone (no ClassifyCNV needed); the full ClassifyCNV
integration runs only if ClassifyCNV is installed. Run: python tests/test_cnv.py"""
import os
import duckdb
import pandas as pd
from pathlib import Path

_MANIFEST = Path(__file__).parent.parent / "acmg" / "manifests" / "cnv_riggs.sql"
_CLASSIFYCNV = os.environ.get("CLASSIFYCNV_DIR", "/root/ClassifyCNV")


def test_riggs_banding_from_scoresheet():
    """The manifest bands total_score into the 5 ClinGen Riggs tiers (no external tool needed)."""
    con = duckdb.connect()
    ev = pd.DataFrame([
        ("dup_path", "chr22", 1, 2, "DUP", 1.00, "MECP2"),
        ("del_lp", "chr1", 1, 2, "DEL", 0.92, "NRXN1"),
        ("vus", "chrX", 1, 2, "DEL", 0.00, ""),
        ("lb", "chr2", 1, 2, "DUP", -0.95, ""),
        ("ben", "chr3", 1, 2, "DUP", -1.00, ""),
    ], columns=["VariantID", "Chromosome", "Start", "End", "Type", "total_score",
                "Known or predicted dosage-sensitive genes"])
    con.register("_ev", ev)
    con.execute("CREATE TABLE cnv_evidence AS SELECT * FROM _ev")
    con.execute(_MANIFEST.read_text())
    got = dict(con.execute("SELECT variant_id, acmg_class FROM cnv_classification").fetchall())
    assert got["dup_path"] == "Pathogenic"
    assert got["del_lp"] == "Likely Pathogenic"
    assert got["vus"] == "VUS"
    assert got["lb"] == "Likely Benign"
    assert got["ben"] == "Benign"


def test_classifycnv_integration_if_available():
    """End-to-end on real CNVs, only when ClassifyCNV is installed (needs bedtools)."""
    if not os.path.isdir(_CLASSIFYCNV):
        print("skip: ClassifyCNV not installed (set CLASSIFYCNV_DIR)")
        return
    from acmg.cnv import classify_cnv
    cnvs = pd.DataFrame([
        dict(chrom="22", start=18893863, end=20307536, svtype="DUP"),  # 22q11.21 dup -> Pathogenic
        dict(chrom="1", start=1000000, end=1050000, svtype="DUP"),      # benign control
    ])
    res = classify_cnv(cnvs).set_index("chromosome")["acmg_class"].to_dict()
    assert res.get("chr22") == "Pathogenic"
    assert res.get("chr1") == "Benign"


if __name__ == "__main__":
    for n, f in sorted(globals().items()):
        if n.startswith("test_") and callable(f):
            f(); print("ok", n)
    print("all passed")
