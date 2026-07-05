#!/usr/bin/env bash
# Ensembl VEP 116 via Docker (ensemblorg/ensembl-vep:release_116.0) against the LOCAL merged 116 GRCh37 cache.
# Software release == cache release (both 116) — the correct fix for the earlier 105-binary-vs-116-cache mismatch.
# Offline, MANE Select transcript (matches our local NMD/PS1/PM5 joins), tab output.
#   bash scripts/vep116.sh input.vcf [output.tsv]
set -euo pipefail
IN="$(readlink -f "$1")"; OUT="${2:-${IN%.vcf}.vep.tsv}"
IMG="${VEP_IMAGE:-ensemblorg/ensembl-vep:release_116.0}"
CACHE_DIR="${VEP_CACHE:-$(cd "$(dirname "$0")/.." && pwd)/.cache/vep}"
[[ -d "$CACHE_DIR/homo_sapiens_merged/116_GRCh37" ]] || { echo "no 116 GRCh37 cache under $CACHE_DIR (extract merged116.tar.gz first)"; exit 1; }
docker run --rm -u "$(id -u):$(id -g)" \
  -v "$CACHE_DIR":/opt/vep/.vep \
  -v "$(dirname "$IN")":/in -v "$(dirname "$(readlink -f "$OUT")")":/out \
  "$IMG" vep --offline --cache --merged --dir_cache /opt/vep/.vep \
    --assembly GRCh37 --cache_version 116 --use_given_ref \
    --mane_select --pick --symbol --canonical --numbers --force_overwrite --no_stats --tab \
    --input_file "/in/$(basename "$IN")" --output_file "/out/$(basename "$OUT")"
echo "wrote $OUT"
