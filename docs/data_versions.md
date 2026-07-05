# Reference data & tool versions (reproducibility)

`uv` / `requirements.lock` pin the **Python** dependencies. They do **not** pin the reference data or the
external tools — and in genomics that is where reproducibility actually lives (a ClinVar release or a VEP
version changes classifications). Pin these too; record the exact versions used for a run.

**Fetch everything with one script:** `bash scripts/fetch_data.sh` (idempotent — skips present files). It is
the single home for every download that used to be ad hoc, and it writes `.cache/MANIFEST.txt` (sha256 + size
+ date) recording what a run actually used. Heavy/optional slices are opt-in: `WITH_CADD_SLICE=1`
(`scripts/build_cadd_slice.sh`), `WITH_VEP=1`. This table is the human-readable version log for that script.

## Annotation engine

- **Ensembl VEP** — the cohort path uses **offline VEP**, pinned via bioconda: **`ensembl-vep 105.0`**
  (`micromamba -n vep`) + the matching **`homo_sapiens_vep_105_GRCh37`** cache under `/root/.vep`. 105 is a
  deliberate pin, not conda's default accident: **bioconda's `noarch` channel tops out at 106.1**; truly
  current VEP (112+) needs a non-bioconda build. On GRCh37 the 105-vs-106 gap is immaterial because we fix
  the reported transcript to **MANE Select** (`--mane_select --pick`), and our GENCODE **v46lift37** is
  MANE-tagged — so VEP consequence and the local NMD/PS1/PM5 joins align on the same MANE transcript
  regardless of the Ensembl-release gap (MANE ENST/NM share identical CDS by construction). To bump: pin
  `ensembl-vep=106.1` and re-fetch the 106 cache. The REST client (`acmg.vep_rest`) stays a small-residual
  fallback, its rolling release recorded per run.

  **Cache-version — RESOLVED via Docker.** The downloaded merged cache is release **116**
  (`homo_sapiens_merged/116_GRCh37/`, ~26 GB extracted under `$HOME/.vep`). Rather than force it under the
  105 bioconda binary (VEP refuses a version mismatch), we run the **matching software from the official
  image `ensemblorg/ensembl-vep:release_116.0`** — `scripts/vep116.sh` (offline, `--merged --cache_version
  116 --assembly GRCh37 --mane_select --pick`). Software release == cache release. The bioconda 105 env stays
  for quick local calls; the Docker 116 path is the version-matched cohort annotator. MANE Select on GRCh37
  means either release gives the same reported call.

## Reference tables (all GRCh37)

| source | version | used for | notes |
|---|---|---|---|
| ClinVar `variant_summary.txt.gz` | release **2026-06-28** | PS1/PM5 + classification + review stars | record the download date |
| ClinVar `clinvar.vcf.gz` (GRCh37) | rolling | **CLNSIGCONF** conflict-composition for PS1/PM5 (via duckhts read_bcf) | optional; recovers P/LP-vs-VUS conflicts, keeps P/LP-vs-Benign out |
| gnomAD | **v2.1.1** exomes (AF) + by-gene constraint (`mis_z`) | BA1/BS1/PM2, PP2 | GRCh37-native |
| GENCODE | **v46lift37** (basic, MANE-tagged) | NMD / transcript (PVS1) | MANE on GRCh37 via lift37 |
| REVEL | **v1.3** | PP3/BP4 (missense) | local `revel_grch37.parquet` (77.9M loci) |
| CADD | **v1.7** GRCh37, exon±500bp SLICE | supplementary non-missense score (optional) | **not required** — PP3/BP4 come from REVEL (missense) + SpliceAI (splice). `scripts/build_cadd_slice.sh` range-reads only exon±500bp off the 79GB whole-genome file → `cadd_exome_slice.parquet`. Indels need the gnomAD-v4-indel-liftover set (not enumerable) |
| HPO (`hp.obo` + `hp.index`) | download **~2026-07** | FastHPOCR phenotype index | `hp.obo` from OBO PURL; `hp.index` built by `.cache/_build_index.py`. Record the HPO release |
| Monarch KG | `monarch-kg/latest` | phenotype→gene, cross-species | **rolling** — snapshot for a reproducible run |
| ClinGen Gene-Disease Validity | download **2026-07-03** | `gene_curation` (gene×disease MOI + validity) | rolling; record FILE CREATED date. `acmg.clingen` |
| ClinGen Dosage Sensitivity | download **2026-07-03** | `hi_score` / LoF mechanism (PVS1) | rolling; record FILE CREATED date. `acmg.clingen` |
| ClassifyCNV | commit **148757c** (v1.1.1) + bundled ClinGen dosage map | CNV Riggs 2020 | pinned in `scripts/setup_classifycnv.sh` |
| OpenSpliceAI | **0.0.7** + torch 2.12.1 (CPU) + OSAI_MANE 10000nt ensemble | on-device SpliceAI (splice PP3/BP7) | **OPTIONAL, out of the core lock** (torch ~2GB). `pip install -e '.[splice]'` or `scripts/setup_openspliceai.sh` (also fetches models to `.cache/openspliceai_models/`). GRCh37 FASTA `human_g1k_v37.fasta` |
| bedtools | **2.31.1** (bioconda) | interval ops inside ClassifyCNV | **system binary — uv/requirements.lock CANNOT pin it**; pin via bioconda (`scripts/setup_classifycnv.sh`) |

**Non-Python tools are a real reproducibility gap for `uv`.** `requirements.lock` pins only Python packages; compiled binaries (bedtools, and htslib under duckhts) are invisible to it. Pin them via **bioconda** (micromamba/conda/pixi) at the versions above, or containerise. The CNV setup script auto-installs the pinned bedtools with `INSTALL_BEDTOOLS=1`.

## For a fully reproducible run, record

1. `requirements.lock` (Python) — committed.
2. The versions above (this file) — update the dates/releases when you refresh a source.
3. The VEP release (REST header, or the pinned offline cache name).
4. The ACMG rule-set is versioned by the SQL in `acmg/manifests/` (git commit).

Guideline migration (2015 → 2020 → v4) is a change of the manifest tables, not the engine.
