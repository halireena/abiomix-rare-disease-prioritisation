#!/usr/bin/env bash
# CADD v1.7 GRCh37 PHRED for the EXONIC slice only (protein-coding exons ±500bp), by remote tabix range-read.
# CADD's whole_genome_SNVs is 79GB; we never download it. We enumerate exon±pad regions from the local GENCODE
# GTF, tabix -R only those bytes off the public host, and combine to one parquet (chrom,pos,ref,alt,cadd_phred).
# Supplementary non-missense breadth (REVEL is missense-only, SpliceAI is splice-only); PP3/BP4 stay REVEL/SpliceAI.
#   bash scripts/build_cadd_slice.sh            # default ±500bp, 12-way
#   PAD=50 JOBS=8 bash scripts/build_cadd_slice.sh
set -euo pipefail
cd "$(dirname "$0")/.."
CACHE="${CACHE:-.cache}"; PAD="${PAD:-500}"; JOBS="${JOBS:-12}"
URL="https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh37/whole_genome_SNVs.tsv.gz"
WORK="$CACHE/cadd_work"; mkdir -p "$WORK"
OUT="$CACHE/cadd_exome_slice.parquet"
[[ -s "$OUT" ]] && { echo "ok  $OUT already built ($(du -h "$OUT"|cut -f1))"; exit 0; }

echo "1/3 exon±${PAD}bp regions from GENCODE GTF -> per-chrom BED (merged)"
python - "$CACHE/gencode.lift37.gtf.gz" "$WORK" "$PAD" <<'PY'
import sys, duckdb
gtf, work, pad = sys.argv[1], sys.argv[2], int(sys.argv[3])
con = duckdb.connect(); con.execute("LOAD duckhts")
# protein-coding exons, ±pad, merged per chrom (BED: chrom, start-1, end). read_gtf flattens attributes.
rows = con.execute(f"""
  WITH e AS (
    SELECT replace(CAST(seqname AS VARCHAR),'chr','') AS chrom,
           greatest(CAST(start AS BIGINT)-{pad},0) AS s, CAST("end" AS BIGINT)+{pad} AS e
    FROM read_gtf('{gtf}') WHERE feature='exon'
  )
  SELECT chrom, s, e FROM e WHERE chrom IN ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y')
  ORDER BY chrom, s
""").fetchall()
from collections import defaultdict
by = defaultdict(list)
for c, s, e in rows: by[c].append((s, e))
tot = 0
for c, ivs in by.items():
    merged = []
    for s, e in sorted(ivs):
        if merged and s <= merged[-1][1]: merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else: merged.append((s, e))
    with open(f"{work}/reg_{c}.bed", "w") as fh:
        for s, e in merged: fh.write(f"{c}\t{s}\t{e}\n")
    tot += len(merged)
print(f"    {tot} merged regions across {len(by)} chroms")
PY

echo "2/3 ${JOBS}-way remote tabix range-read (CADD v1.7) -> cadd_<chrom>.tsv.gz"
( cd "$WORK"
  ls reg_*.bed | sed 's/reg_//;s/.bed//' | xargs -P "$JOBS" -I{} bash -c \
    "tabix -R reg_{}.bed '$URL' | awk -F'\t' 'BEGIN{OFS=\"\t\"}{print \$1,\$2,\$3,\$4,\$6}' | gzip > cadd_{}.tsv.gz" )

echo "3/3 combine -> $OUT (chrom,pos,ref,alt,cadd_phred)"
python - "$WORK" "$OUT" <<'PY'
import sys, glob, duckdb
work, out = sys.argv[1], sys.argv[2]
files = sorted(glob.glob(f"{work}/cadd_*.tsv.gz"))
con = duckdb.connect()
con.execute(f"""COPY (
  SELECT column0 AS chrom, CAST(column1 AS BIGINT) AS pos, column2 AS ref, column3 AS alt,
         CAST(column4 AS DOUBLE) AS cadd_phred
  FROM read_csv({files!r}, delim='\t', header=false, columns={{'column0':'VARCHAR','column1':'BIGINT','column2':'VARCHAR','column3':'VARCHAR','column4':'DOUBLE'}})
) TO '{out}' (FORMAT parquet)""")
n = con.execute(f"SELECT count(*) FROM read_parquet('{out}')").fetchone()[0]
print(f"    {n} CADD sites -> {out}")
PY
echo "done."
