#!/usr/bin/env bash
# Optional ON-DEVICE SPLICING lane (OpenSpliceAI, SpliceAI-equivalent) — SpliceAI-quality splice scores with no
# auth-gated Illumina precompute. HEAVY: pulls PyTorch (~2GB). Kept OUT of the core requirements.lock so the base
# student install stays lean (same policy as VEP-via-Docker and ClassifyCNV-via-git). Run this only if you want
# the local `spliceai` column (splice PP3 >=0.2 / BP7 <0.2, Walker 2023) filled on-device.
#
#   bash scripts/setup_openspliceai.sh
#
# Fills the `spliceai` source of the annotation datalake: scripts/annotate_cohort.py picks up
# .cache/spliceai_openspliceai.parquet (produced via acmg.annotate.spliceai_local) as its SpliceAI increment.
set -euo pipefail
cd "$(dirname "$0")/.."
CACHE="${CACHE:-.cache}"
MODELS_DIR="$CACHE/openspliceai_models"
OSAI_REPO="${OSAI_REPO:-https://github.com/Kuanhao-Chao/OpenSpliceAI.git}"

echo "== 1/2 install the [splice] extra (PyTorch CPU wheel) =="
# CPU-only torch keeps it to ~200MB of wheels instead of the CUDA stack; drop --index-url for a GPU box.
pip install --index-url https://download.pytorch.org/whl/cpu torch==2.12.1 || pip install torch
pip install -e '.[splice]'

echo "== 2/2 OSAI_MANE model ensembles (400/2000/10000nt) -> $MODELS_DIR =="
# The pretrained human/MANE ensembles ship IN the OpenSpliceAI repo (models/openspliceai-mane/<size>nt/
# model_<size>nt_rs*.pt; ~11MB total for all three sizes). --flanking selects which set: 400=+-200 (~50x
# faster), 2000=+-1000 (~6.5x), 10000=+-5000 (ClinGen/Walker-calibrated). Each size is a DIFFERENT model, so
# flanking must match the model set. Genome-agnostic (score against any FASTA); we use GRCh37 human_g1k_v37.
mkdir -p "$MODELS_DIR"
if ls "$MODELS_DIR"/mane_10000nt/model_10000nt_rs*.pt >/dev/null 2>&1; then
  echo "  ok  models already present ($(ls -d "$MODELS_DIR"/mane_*nt 2>/dev/null | wc -l) sizes)"
else
  TMP="$(mktemp -d)"; echo "  cloning $OSAI_REPO (shallow) ..."
  git clone --depth 1 "$OSAI_REPO" "$TMP/OpenSpliceAI"
  for SZ in 400 2000 10000; do
    src="$TMP/OpenSpliceAI/models/openspliceai-mane/${SZ}nt"
    [ -d "$src" ] && cp -r "$src" "$MODELS_DIR/mane_${SZ}nt" && echo "  installed mane_${SZ}nt ($(ls "$MODELS_DIR/mane_${SZ}nt"/*.pt | wc -l) models)"
  done
  rm -rf "$TMP"
fi
echo "  Reference FASTA (GRCh37): /root/GRCh37/human_g1k_v37.fasta (or pass --fasta / spliceai_local(fasta=...))."

echo "done. Verify:  PYTHONPATH=. python3 -c \"from acmg.annotate import spliceai_local; print(spliceai_local([('11',5248232,'T','A')], model_dir='$MODELS_DIR/mane_400nt', flanking_size=400, threads=2))\""
