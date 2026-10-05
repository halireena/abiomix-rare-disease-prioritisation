"""Build the FastHPOCR concept-recognition index (`.cache/hp.index`) from the HPO ontology file (`hp.obo`).

Stage 2 (HPO extraction), tests/test_hpo_modes.py and the per-proband drivers read `.cache/hp.index`. It is
derived data (git-ignored), so build it once after downloading hp.obo. Needs network only for that download:

    curl -fL -o .cache/hp.obo https://purl.obolibrary.org/obo/hp.obo
    pip install FastHPOCR pronto          # pronto: FastHPOCR's indexer needs it but does not declare it
    python scripts/build_hpo_index.py     # -> .cache/hp.index (~140 MB; about 15 minutes on one CPU)

Record the HPO release (the `data-version:` line at the top of hp.obo) in docs/data_versions.md.
"""
from __future__ import annotations
import argparse
import os
import sys


def main(argv: list[str] | None = None) -> str:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--obo", default=os.path.join(".cache", "hp.obo"), help="path to hp.obo (default: .cache/hp.obo)")
    ap.add_argument("--out-dir", default=".cache", help="folder to write hp.index into (default: .cache)")
    a = ap.parse_args(argv)

    if not os.path.isfile(a.obo):
        sys.exit(f"hp.obo not found at {a.obo!r}. Download it first (needs network):\n"
                 f"  mkdir -p .cache && curl -fL -o .cache/hp.obo https://purl.obolibrary.org/obo/hp.obo")
    try:
        from FastHPOCR.IndexHPO import IndexHPO
    except ImportError as e:
        sys.exit(f"cannot import the FastHPOCR indexer ({e}). Install it with:  pip install FastHPOCR pronto")

    os.makedirs(a.out_dir, exist_ok=True)
    IndexHPO(a.obo, a.out_dir).index()
    out = os.path.join(a.out_dir, "hp.index")
    if not os.path.isfile(out):
        sys.exit(f"FastHPOCR finished but {out} was not written; check the messages above.")
    print(f"wrote {out}")
    return out


if __name__ == "__main__":
    main()
