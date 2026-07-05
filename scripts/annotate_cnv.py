"""CNV incremental annotation datalake — the CNV twin of scripts/annotate_cohort.py, on the SAME key-generic
store (keyed by cnv_key). Distinct CNVs are classified ONCE via ClassifyCNV (ClinGen/ACMG Riggs 2020); the
expensive scoring accumulates in a per-release datalake and only the INCREMENT (new CNVs not yet in the release)
is (re)classified — reused across cohorts. Same store_increment / store_merge / store_read_cohort as the SNV
path, just key='cnv_key'. GRCh37.

Setup once: bash scripts/setup_classifycnv.sh   (clones pinned ClassifyCNV, needs bedtools).
Run:  PYTHONPATH=. python3 scripts/annotate_cnv.py --parquet /root/bioconnect/prepared/cnv.parquet --cache .cache
"""
from __future__ import annotations
import argparse
from pathlib import Path
import duckdb
import pandas as pd
from acmg.cnv import classify_cnv
from scripts.annotate_cohort import store_increment, store_merge, store_read_cohort

CNV_RELEASE = "classifycnv-v1.1.1/riggs2020"      # pinned ClassifyCNV commit + Riggs manifest = the CNV source release


def build_cnv_cohort(con: duckdb.DuckDBPyConnection, cnv_parquet: str) -> int:
    """Distinct CNVs to classify ONCE -> cnv_cohort(cnv_key, chrom, start, end, svtype). This bundle emits a
    generic svtype='CNV', so DEL/DUP is derived from copy number (cn<2 -> DEL, cn>2 -> DUP; cn==2 is neutral,
    dropped). cnv_key = chrom(no 'chr')-start-end-svtype is the CNV IDENTITY (cn is a per-call attribute used only
    to derive direction, not part of the key), so DISTINCT gives one row per CNV."""
    con.execute(f"""
        CREATE OR REPLACE TABLE cnv_cohort AS
        SELECT DISTINCT
            replace(CAST(chrom AS VARCHAR),'chr','')||'-'||CAST(start AS VARCHAR)||'-'||CAST("end" AS VARCHAR)||'-'||sv AS cnv_key,
            replace(CAST(chrom AS VARCHAR),'chr','') AS chrom, CAST(start AS BIGINT) AS start,
            CAST("end" AS BIGINT) AS "end", sv AS svtype
        FROM (
            SELECT chrom, start, "end", CASE WHEN cn < 2 THEN 'DEL' WHEN cn > 2 THEN 'DUP' END AS sv
            FROM read_parquet('{cnv_parquet}')
            WHERE start IS NOT NULL AND "end" IS NOT NULL AND cn IS NOT NULL
        ) WHERE sv IS NOT NULL
    """)
    return con.execute("SELECT count(*) FROM cnv_cohort").fetchone()[0]


def annotate_cnv(con: duckdb.DuckDBPyConnection, inc: pd.DataFrame, *, classifycnv_dir: str | None = None) -> pd.DataFrame:
    """Classify the CNV increment via ClassifyCNV (Riggs 2020) -> (cnv_key, cnv_class, cnv_points). The expensive
    step, cached in the CNV datalake exactly as VEP is for SNVs. Empty frame if the increment is empty."""
    cols = ["cnv_key", "cnv_class", "cnv_points", "dosage_genes"]
    if not len(inc):
        return pd.DataFrame(columns=cols)
    res = classify_cnv(inc, con=duckdb.connect(), build="hg19", classifycnv_dir=classifycnv_dir)
    if not len(res):
        return pd.DataFrame(columns=cols)
    con.register("_res", res)
    return con.execute("""
        SELECT DISTINCT replace(CAST(chromosome AS VARCHAR),'chr','')||'-'||CAST(start_pos AS VARCHAR)||'-'||CAST(end_pos AS VARCHAR)||'-'||type AS cnv_key,
               acmg_class AS cnv_class, total_score AS cnv_points, dosage_genes
        FROM _res
    """).df().drop_duplicates("cnv_key")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--parquet", default="/root/bioconnect/prepared/cnv.parquet")
    ap.add_argument("--cache", default=".cache")
    ap.add_argument("--store", default=".cache/cnv_store.parquet", help="CNV annotation datalake (keyed by cnv_key)")
    ap.add_argument("--out", default=".cache/cnv_classified.parquet")
    ap.add_argument("--classifycnv-dir", default=None, help="ClassifyCNV clone (default $CLASSIFYCNV_DIR or /root/ClassifyCNV)")
    args = ap.parse_args()

    con = duckdb.connect()
    n = build_cnv_cohort(con, args.parquet)
    inc = store_increment(con, args.store, "classifycnv_release", CNV_RELEASE, cohort="cnv_cohort", key="cnv_key")
    print(f"distinct CNVs: {n}; increment (new, not in release {CNV_RELEASE}): {len(inc)}")
    if len(inc):
        got = annotate_cnv(con, inc, classifycnv_dir=args.classifycnv_dir)
        total = store_merge(con, args.store, got, "classifycnv_release", CNV_RELEASE, key="cnv_key")
        print(f"  classified {len(got)} new CNVs -> datalake now {total} ({args.store})")
    else:
        print("  increment empty — every CNV already in the release (full reuse, no ClassifyCNV run).")

    cls = store_read_cohort(con, args.store, cohort="cnv_cohort", key="cnv_key", cols=["cnv_key", "cnv_class", "cnv_points", "dosage_genes"])
    con.register("_cls", cls)
    con.execute(f"COPY (SELECT * FROM _cls) TO '{args.out}' (FORMAT parquet)")
    print(f"wrote {len(cls)} CNV classifications -> {args.out}")
    if len(cls):
        print("  Riggs banding:", cls["cnv_class"].value_counts().to_dict())


if __name__ == "__main__":
    main()
