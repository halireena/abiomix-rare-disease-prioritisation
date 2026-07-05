"""Format-agnostic ingest tests: VCF (duckhts), tidy TSV (custom mapping), Excel, and PED pedigree — each
producing the SAME canonical `sample_call`, verified end-to-end through acmg.family. Run:
`PYTHONPATH=. python3 tests/test_ingest.py`."""
import os
import duckdb
from acmg import ingest
from acmg.ingest import CANONICAL
from acmg import family

FIX = os.path.join(os.path.dirname(__file__), "fixtures")


def _assert_canonical(con):
    cols = [c[0] for c in con.execute("DESCRIBE sample_call").fetchall()]
    assert cols == CANONICAL, f"schema drift: {cols}"


def test_vcf_ingest_unpivots_per_sample():
    con = duckdb.connect()
    ingest.ingest_variants(con, os.path.join(FIX, "trio.vcf"))
    _assert_canonical(con)
    # 3 samples x 3 SNVs = 9 rows; one row per sample x variant.
    assert con.execute("SELECT count(*) FROM sample_call").fetchone()[0] == 9
    assert set(con.execute("SELECT DISTINCT member_id FROM sample_call").fetchall()) == {("PROB",), ("MOM",), ("DAD",)}
    # sample id doubles as case_id; FORMAT parsed.
    r = con.execute("SELECT gt, gq, dp FROM sample_call WHERE member_id='PROB' AND variant_key='1-100-A-G'").fetchone()
    assert r == ("0/1", 40, 30)


def test_vcf_plus_ped_de_novo():
    con = duckdb.connect()
    ingest.ingest_variants(con, os.path.join(FIX, "trio.vcf"))
    ingest.ingest_pedigree(con, os.path.join(FIX, "trio.ped"))
    ped = con.execute("SELECT proband_id, mother_id, father_id FROM pedigree WHERE family_id='FAM1'").fetchone()
    assert ped == ("PROB", "MOM", "DAD")
    assert con.execute("SELECT DISTINCT family_id FROM sample_call").fetchone()[0] == "FAM1"  # backfilled
    family.build_carriage(con)
    dn = family.de_novo_candidates(con).df()
    # only 1-100 is absent in both called parents; 1-200 is inherited from MOM; 2-300 from DAD.
    assert set(dn.variant_key) == {"1-100-A-G"}


def test_tsv_custom_column_mapping():
    con = duckdb.connect()
    cols = {"case_id": "sample_id", "family_id": "fam", "member_id": "sample_id",
            "proband_status": "status", "relationship": "rel", "chrom": "chr", "pos": "position",
            "ref": "ref_allele", "alt": "alt_allele", "gt": "genotype", "gq": "geno_qual", "dp": "depth"}
    ingest.ingest_variants(con, os.path.join(FIX, "tidy.tsv"), columns=cols)
    _assert_canonical(con)
    assert con.execute("SELECT count(*) FROM sample_call").fetchone()[0] == 6
    # 'chr' prefix stripped; variant_kind derived (snv); variant_key built.
    assert con.execute("SELECT DISTINCT chrom FROM sample_call").fetchone()[0] == "1"
    assert con.execute("SELECT DISTINCT variant_kind FROM sample_call").fetchone()[0] == "snv"
    # curated relationship fields -> pedigree, then de novo.
    ingest.ingest_pedigree(con)
    assert con.execute("SELECT proband_id, mother_id, father_id FROM pedigree WHERE family_id='FX'").fetchone() == ("S1", "S2", "S3")
    family.build_carriage(con)
    assert set(family.de_novo_candidates(con).df().variant_key) == {"1-100-A-G"}


def test_excel_ingest_matches_tsv():
    import pandas as pd
    con = duckdb.connect()
    df = pd.read_csv(os.path.join(FIX, "tidy.tsv"), sep="\t")
    xlsx = os.path.join(FIX, "_tidy.xlsx")
    df.to_excel(xlsx, index=False)
    try:
        cols = {"case_id": "sample_id", "family_id": "fam", "member_id": "sample_id",
                "proband_status": "status", "relationship": "rel", "chrom": "chr", "pos": "position",
                "ref": "ref_allele", "alt": "alt_allele", "gt": "genotype", "gq": "geno_qual", "dp": "depth"}
        ingest.ingest_variants(con, xlsx, columns=cols)
        _assert_canonical(con)
        assert con.execute("SELECT count(*) FROM sample_call").fetchone()[0] == 6
        assert con.execute("SELECT DISTINCT chrom FROM sample_call").fetchone()[0] == "1"
    finally:
        os.remove(xlsx)


def test_snv_only_filter_keeps_snvs():
    con = duckdb.connect()
    ingest.ingest_variants(con, os.path.join(FIX, "trio.vcf"), snv_only=False)
    n_all = con.execute("SELECT count(*) FROM sample_call").fetchone()[0]
    ingest.ingest_variants(con, os.path.join(FIX, "trio.vcf"), snv_only=True)
    n_snv = con.execute("SELECT count(*) FROM sample_call").fetchone()[0]
    assert n_all == n_snv == 9  # fixture is all SNVs; filter keeps them, drops nothing


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name)
    print("all passed")
