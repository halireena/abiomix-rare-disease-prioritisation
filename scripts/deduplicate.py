#!/usr/bin/env python3
"""
deduplicate.py

Reproduces the SNV/CNV deduplication of the clinical variant dataset.

Steps:
  1. Read the curated clinical dataset (parquet).
  2. Filter to variant_kind == 'snv' AND proband_status_curated in
     ('proband', 'single_member_case'), and deduplicate on (CHROM, POS, REF, ALT).
  3. Filter to variant_kind == 'cnv' and deduplicate on
     (CHROM, POS, INFO_END, INFO_SVTYPE).

Note: CHROM is stored in mixed naming conventions in the source (e.g. both '1'
and 'chr1'); this script deduplicates on the raw CHROM string and does NOT
normalize them, matching the agreed 'snv + patient filter' definition.
  4. Write unique_snvs_for_annotation.parquet and unique_cnvs.parquet.
  5. Append a run record to dedup_log.txt (input rows, output rows, timestamp).

Deduplication uses DuckDB's DISTINCT ON so that one full row is kept per key,
made deterministic by ordering on the key columns.
"""

from datetime import datetime, timezone

import duckdb
import pandas as pd

# --- Configuration -----------------------------------------------------------
INPUT_FILE = "dataset_clinical_curated.parquet"
SNV_OUTPUT = "unique_snvs_for_annotation.parquet"
CNV_OUTPUT = "unique_cnvs.parquet"
LOG_FILE = "dedup_log.txt"

SNV_KEYS = ["CHROM", "POS", "REF", "ALT"]
SNV_COLS = ["CHROM", "POS", "REF", "ALT"]
# SNVs are restricted to proband / single-member cases (original request spec).
SNV_WHERE_EXTRA = ("AND proband_status_curated IN "
                   "('proband', 'single_member_case')")

CNV_KEYS = ["CHROM", "POS", "INFO_END", "INFO_SVTYPE"]
CNV_COLS = ["CHROM", "POS", "INFO_END", "INFO_SVTYPE",
            "INFO_SVLEN", "INFO_NEXONS", "FORMAT_CN"]
CNV_WHERE_EXTRA = ""


def dedup(con, kind, keys, cols, output, where_extra=""):
    """Filter to `kind` (plus any `where_extra`), deduplicate on `keys`,
    write `cols` to `output`.

    Returns (input_rows_after_filter, output_rows).
    """
    key_list = ", ".join(keys)
    col_list = ", ".join(cols)
    where = f"WHERE variant_kind = '{kind}' {where_extra}"

    con.sql(f"""
        COPY (
            SELECT DISTINCT ON ({key_list})
                {col_list}
            FROM '{INPUT_FILE}'
            {where}
            ORDER BY {key_list}
        ) TO '{output}' (FORMAT PARQUET)
    """)

    in_rows = con.sql(
        f"SELECT COUNT(*) FROM '{INPUT_FILE}' {where}"
    ).fetchone()[0]
    out_rows = con.sql(f"SELECT COUNT(*) FROM '{output}'").fetchone()[0]
    return in_rows, out_rows


def main():
    run_ts = datetime.now(timezone.utc).astimezone()

    # --- Version / provenance logging (stdout) -------------------------------
    header = [
        "=" * 60,
        "Deduplication run",
        "=" * 60,
        f"Date run     : {run_ts.isoformat()}",
        f"Input file   : {INPUT_FILE}",
        f"duckdb       : {duckdb.__version__}",
        f"pandas       : {pd.__version__}",
        "=" * 60,
    ]
    print("\n".join(header))

    con = duckdb.connect()

    total_rows = con.sql(f"SELECT COUNT(*) FROM '{INPUT_FILE}'").fetchone()[0]
    print(f"Total input rows: {total_rows:,}")

    snv_in, snv_out = dedup(con, "snv", SNV_KEYS, SNV_COLS, SNV_OUTPUT,
                            SNV_WHERE_EXTRA)
    print(f"SNV (proband/single_member_case): {snv_in:,} rows -> "
          f"{snv_out:,} unique  ({SNV_OUTPUT})")

    cnv_in, cnv_out = dedup(con, "cnv", CNV_KEYS, CNV_COLS, CNV_OUTPUT,
                            CNV_WHERE_EXTRA)
    print(f"CNV: {cnv_in:,} rows -> {cnv_out:,} unique  ({CNV_OUTPUT})")

    # --- Summary log file ----------------------------------------------------
    log_lines = [
        "=" * 60,
        f"Deduplication run: {run_ts.isoformat()}",
        "=" * 60,
        f"Input file        : {INPUT_FILE}",
        f"duckdb version    : {duckdb.__version__}",
        f"pandas version    : {pd.__version__}",
        "",
        f"Total input rows  : {total_rows:,}",
        "",
        f"SNV filter        : variant_kind='snv' {SNV_WHERE_EXTRA}",
        f"SNV input rows    : {snv_in:,}  (after filter)",
        f"SNV output rows   : {snv_out:,}  -> {SNV_OUTPUT}",
        "",
        f"CNV input rows    : {cnv_in:,}",
        f"CNV output rows   : {cnv_out:,}  -> {CNV_OUTPUT}",
        "",
    ]
    # Append so repeated runs keep a history; the file records every run.
    with open(LOG_FILE, "a") as fh:
        fh.write("\n".join(log_lines) + "\n")

    print(f"\nLog written to {LOG_FILE}")


if __name__ == "__main__":
    main()
