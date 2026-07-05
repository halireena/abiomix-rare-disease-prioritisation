"""Reshape the monolithic tidy bundle (`dataset.parquet` — one row per sample x variant, 7.7M rows, cryptic
FORMAT_* columns) into a saner, pipeline-ready layout that students run end-to-end from:

IMPORTANT: source the DE-IDENTIFIED student file `dataset.parquet` (its `clinical_indication_text` is the
de-identified text). Do NOT source `dataset_clinical_curated.parquet` — that curated file carries raw PII
(clinician names, DOB, accession IDs) in `clinical_indication_text`. The variant calls are identical in both.


  cases.parquet     one row per case: case_id, family_id, member_id, proband_status, relationship,
                    family_size, clinical_indication_text, ancestry
  pedigree.parquet  one row per family: proband_id / mother_id / father_id + proband_ambiguous (multiplex flag)
  snv.parquet       canonical SNV/indel calls (acmg.ingest schema: case/member/variant_key/gt/gq/dp/...)
  cnv.parquet       canonical CNV calls: case/member/chrom/start/end/svtype/cn

Usage:  python scripts/prepare_bundle.py --source /root/bioconnect/dataset_clinical_curated.parquet --out prepared
"""
from __future__ import annotations
import argparse
from pathlib import Path
import duckdb
from acmg import ingest, family


def prepare(source: str, out: str) -> dict:
    con = duckdb.connect()
    outp = Path(out); outp.mkdir(parents=True, exist_ok=True)
    src = f"read_parquet('{source}')"

    # cases: one clean row per case
    con.execute(f"""COPY (
        SELECT DISTINCT student_case_id AS case_id, student_family_id AS family_id, student_member_id AS member_id,
               family_size, proband_status_curated AS proband_status, family_relationship_curated AS relationship,
               clinical_indication_text, ancestry_or_ethnicity AS ancestry
        FROM {src}
    ) TO '{outp}/cases.parquet' (FORMAT parquet)""")

    # canonical SNV/indel calls via the format-agnostic ingest
    ingest.load_bundle(con, source, snv_only=True)
    con.execute(f"COPY sample_call TO '{outp}/snv.parquet' (FORMAT parquet)")
    family.build_pedigree(con)
    con.execute(f"COPY pedigree TO '{outp}/pedigree.parquet' (FORMAT parquet)")

    # canonical CNV calls (carry the case ids + SV fields the SNV schema drops)
    con.execute(f"""COPY (
        SELECT student_case_id AS case_id, student_family_id AS family_id, student_member_id AS member_id,
               replace(CAST(CHROM AS VARCHAR),'chr','') AS chrom, CAST(POS AS BIGINT) AS start,
               TRY_CAST(INFO_END AS BIGINT) AS "end", CAST(INFO_SVTYPE AS VARCHAR) AS svtype,
               TRY_CAST(FORMAT_CN AS INTEGER) AS cn, FORMAT_GT AS gt
        FROM {src} WHERE variant_kind = 'cnv'
    ) TO '{outp}/cnv.parquet' (FORMAT parquet)""")

    cnt = lambda f: con.execute(f"SELECT count(*) FROM read_parquet('{outp}/{f}')").fetchone()[0]
    stats = {f: cnt(f) for f in ("cases.parquet", "pedigree.parquet", "snv.parquet", "cnv.parquet")}
    print(f"wrote {out}/: " + ", ".join(f"{k}={v}" for k, v in stats.items()))
    return stats


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", default="/root/bioconnect/dataset.parquet")  # de-identified student file
    ap.add_argument("--out", default="prepared")
    a = ap.parse_args()
    prepare(a.source, a.out)


if __name__ == "__main__":
    main()
