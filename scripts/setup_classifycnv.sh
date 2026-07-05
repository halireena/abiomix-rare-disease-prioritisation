#!/usr/bin/env bash
# Reproducible setup for the CNV scorer used by acmg.cnv. ClassifyCNV is NOT vendored (it ships a ~120 MB
# ClinGen dosage map); this clones a PINNED commit so results are reproducible.
#
# bedtools is a compiled C++ binary, NOT a Python package -> uv / requirements.lock CANNOT pin it (that is the
# reproducibility gap). Pin non-Python tools via bioconda instead; the version is pinned here + in
# docs/data_versions.md. Set INSTALL_BEDTOOLS=1 to auto-install the pinned build via micromamba/conda/pixi.
#   usage: bash scripts/setup_classifycnv.sh [target_dir]   (INSTALL_BEDTOOLS=1 to fetch pinned bedtools)
set -euo pipefail
DIR="${1:-ClassifyCNV}"
PIN="${CLASSIFYCNV_COMMIT:-148757c}"        # pinned Genotek/ClassifyCNV commit (v1.1.1)
BEDTOOLS_VERSION="${BEDTOOLS_VERSION:-2.31.1}"  # pinned bioconda build (see docs/data_versions.md)

pinned_bedtools_install() {  # try, in order, whatever conda-family manager exists
  for mgr in micromamba mamba conda; do
    command -v "$mgr" >/dev/null 2>&1 && { "$mgr" install -y -c bioconda -c conda-forge "bedtools=${BEDTOOLS_VERSION}"; return; }
  done
  command -v pixi >/dev/null 2>&1 && { pixi global install "bedtools=${BEDTOOLS_VERSION}"; return; }
  echo "ERROR: no conda-family manager found; install pinned bedtools=${BEDTOOLS_VERSION} via bioconda manually"; exit 1
}

if ! command -v bedtools >/dev/null 2>&1; then
  if [ "${INSTALL_BEDTOOLS:-0}" = "1" ]; then pinned_bedtools_install; else
    echo "ERROR: bedtools not found. Install the PINNED build (uv cannot — it is a system binary):"
    echo "  micromamba install -c bioconda -c conda-forge bedtools=${BEDTOOLS_VERSION}"
    echo "  (or re-run with INSTALL_BEDTOOLS=1 to auto-install)"; exit 1
  fi
fi
have="$(bedtools --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' || echo '?')"
[ "$have" = "$BEDTOOLS_VERSION" ] || echo "WARNING: bedtools $have != pinned $BEDTOOLS_VERSION (results may differ; see docs/data_versions.md)"
if [ ! -d "$DIR/.git" ]; then
  git clone https://github.com/Genotek/ClassifyCNV.git "$DIR"
fi
git -C "$DIR" fetch -q --all 2>/dev/null || true
git -C "$DIR" checkout -q "$PIN"
echo "ClassifyCNV pinned at $(git -C "$DIR" rev-parse --short HEAD) in $DIR"
echo "  export CLASSIFYCNV_DIR=$(realpath "$DIR")"
