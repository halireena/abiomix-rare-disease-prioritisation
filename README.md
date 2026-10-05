# Re-Analysing the Unsolved

### A transparent, reproducible pipeline for prioritising unsolved rare-disease cases for genomic re-analysis

**BioConnect 2026 Consultancy Sprint — Abiomix Challenge**
Variant Re-Analysis and Case Re-Classification at Scale

Team: Halireena (Rush), Nubla, Rayane
University of Birmingham Dubai · 2-5 July 2026

> [!WARNING]
> **Research use only. Not for clinical diagnosis.** This is a hackathon proof of concept built and evaluated
> on a synthetic dataset. Its ACMG/AMP classifications are automated evidence summaries that cover only the
> criteria a machine can apply (no functional, literature or case-level expert evidence), and its rankings are
> triage suggestions. Nothing it outputs is a diagnosis or a clinical report. Any real-world use needs
> validation and sign-off by qualified clinical scientists under an accredited laboratory process.

---

## Contents

[What this does](#what-this-does) ·
[Quick start (5 minutes, offline)](#quick-start-5-minutes-offline) ·
[What runs where](#what-runs-offline-and-what-needs-data-or-network) ·
[Overview](#overview) ·
[Key results](#key-results) ·
[Pipeline](#pipeline-architecture) ·
[Data schema](#data-schema) ·
[Rubric](#prioritisation-rubric) ·
[Validation](#validation-and-benchmarking) ·
[Where to get the data](#where-to-get-the-data) ·
[Glossary](#glossary) ·
[Installation](#installation) ·
[Tests](#running-the-tests) ·
[Troubleshooting](#troubleshooting) ·
[Repository structure](#repository-structure) ·
[Limitations](#limitations) ·
[Licence and citation](#licence-and-citation)

---

## What this does

Given the variants and clinical notes of a rare-disease case that is still unsolved, the pipeline:

1. classifies each DNA variant with a transparent, rule-based **ACMG/AMP** kernel (written in SQL; every call
   lists the evidence codes that fired and their points),
2. turns the clinical notes into **HPO** phenotype terms,
3. ranks candidate genes by how well they match the phenotype and the family's inheritance pattern, and
4. gives each case a **Reportable / Review / Insufficient** tier with a short rationale, so a laboratory can
   decide which historical cases deserve expert re-analysis first.

**Who it is for:** bioinformaticians, students and reviewers who want to see how variant re-analysis can be
made explainable and auditable. You need some Python; the [glossary](#glossary) covers the genetics terms.
It is not a clinical tool (see the warning above).

The Python package is imported as `acmg`. It installs under the project name `bioconnect-sprint-py`, the
name used during the sprint.

---

## Quick start (5 minutes, offline)

This runs the classification kernel on a small built-in example. It needs only Python 3.10+ and the core
dependencies: no reference data, no network after `pip install`, no patient data.

```bash
git clone https://github.com/halireena/abiomix-rare-disease-prioritisation.git
cd abiomix-rare-disease-prioritisation
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -e '.[dev]'                                 # duckdb, pandas, requests + pytest tooling
python examples/demo.py
```

The demo first prints the four made-up input variants in the kernel's `annotations` format, then the
classification. Expected output (last part):

```
ACMG classification:
        variant_key   gene variant_kind  total_points                         criteria         acmg_class  pm1_hotspot
0  17-43093464-AC-A  BRCA1        indel           9.0  PM2_Supporting, PVS1_VeryStrong  Likely Pathogenic            0
1      1-100000-A-G  SCN1A          snv           5.0       PM2_Supporting, PP3_Strong                VUS            0
2      2-200000-C-T    XYZ          snv           1.0                              PP3                VUS            0
3    6-26090951-C-G    HFE          snv          -8.0                              BA1             Benign            0
```

How to read it:

| Row | What happened |
|-----|---------------|
| BRCA1 frameshift | Null variant in a gene where loss of function causes disease (**PVS1**, +8) and rare in gnomAD (**PM2_Supporting**, +1). 9 points = Likely Pathogenic. |
| SCN1A missense | Rare (+1) and a very high REVEL score (**PP3_Strong**, +4). 5 points = still a VUS: computational evidence alone cannot make it pathogenic. |
| XYZ missense | No frequency data, so PM2 **abstains** (missing data is not treated as "rare"). REVEL 0.70 gives PP3 (+1). VUS. |
| HFE missense | 13.5% allele frequency, far above the 5% stand-alone benign threshold (**BA1**). Benign. |

Points follow Tavtigian et al. 2020: Pathogenic ≥ 10, Likely Pathogenic 6 to 9, VUS 0 to 5, Likely Benign
-1 to -6, Benign ≤ -7, and BA1 on its own means Benign. The demo's coordinates are illustrative, not real loci.

Then check everything works:

```bash
python -m pytest -q      # ~80 tests in a few seconds; ~10 skip without optional resources (expected)
```

---

## What runs offline and what needs data or network

| Command | Offline? | Needs |
|---------|----------|-------|
| `python examples/demo.py` | Yes | Core install only |
| `python -m pytest -q` | Yes | `.[dev]` extras. VCF-ingest tests download the `duckhts` DuckDB extension on first use and skip without network; HPO-mode tests skip until `.cache/hp.index` exists |
| `ruff check .` | Yes | `pip install ruff` |
| `bash scripts/fetch_data.sh` | No | Network; downloads public reference data (about 1 GB) into `.cache/` (see [Where to get the data](#where-to-get-the-data)) |
| `python scripts/build_hpo_index.py` | Yes, after download | `.cache/hp.obo` plus `pip install FastHPOCR pronto`; about 15 minutes, writes a ~140 MB `.cache/hp.index` |
| `python examples/run_case.py CASE0003` | No | The Abiomix challenge bundle (not public), `.cache/` reference data, network (Ensembl VEP REST, Monarch KG, `duckhts`). HPO extraction is local by default; the `pi` CLI is needed only with `--hpo-mode <llm mode> --allow-external-llm` (see [Privacy defaults](#privacy-defaults)) |
| `python examples/literature_arm.py` | No | Network (Europe PMC, MARRVEL, DECIPHER, LitVar2) and the `pi` CLI |
| `PYTHONPATH=. python scripts/<name>.py --help` | Yes | Each driver documents its inputs; most need the challenge bundle and `.cache/` data. Run from the repository root with `PYTHONPATH=.` (several import `scripts.*`) |

The challenge bundle (`dataset.parquet`, synthetic, about 7.7 million variant rows) was provided to sprint
participants by Abiomix and is not redistributed here. Everything else the pipeline uses is public.

---

## Overview

Around half of rare-disease patients remain without a molecular diagnosis after their initial genomic analysis. Clinical genomics, however, is not static: new gene-disease associations, updated variant databases, refined phenotype ontologies, and improved computational tools mean that cases which were unsolvable a few years ago may be solvable today. Published studies show that systematic re-analysis increases diagnostic yield by 10-15%.

The practical bottleneck is scale. When a laboratory holds hundreds or thousands of historical cases, it must decide which ones justify limited expert time. This project delivers a proof-of-concept pipeline that converts historical case information into a structured phenotype x genotype representation and produces an explainable, case-level prioritisation - without bypassing clinical oversight or making unvalidated diagnostic claims.

All work was carried out on a synthetic dataset provided by Abiomix. No real patient data was used at any stage. All outputs are traceable to non-sensitive inputs, and any language-model output is treated strictly as decision support, not clinical decision-making.

---

## Key Results

- 138 de-identified cases across 119 families processed end to end.
- 7.7 million variant rows reduced to 430,284 unique variants after cleaning and deduplication (a 93% reduction).
- A transparent, rule-based ACMG classification kernel achieving 96.3% concordance with ClinVar across more than 84,000 variants.
- Independent cross-validation against Exomiser, including a direct top-ranked match on KMT2D and top-five agreement on candidate genes in four of five phenotyped probands.
- Six probands examined in depth, each assigned a Reportable, Review, or Insufficient tier with a documented rationale.

Representative findings from the deep-dive probands:

| Case | Lead candidate(s) | Tier | Basis |
|------|-------------------|------|-------|
| CASE0133 | CCNO | Reportable | Homozygous proband, both parents carriers; textbook autosomal-recessive segregation confirmed against source data |
| CASE0008 | PACS1 / EBF3 | Reportable | Phenotype-confirmed; EBF3 re-prioritised from a variant of uncertain significance and corroborated by Exomiser (ranked second) |
| CASE0007 | FBN1 / CUL7 / COL27A1 | Reportable | Co-segregation with affected father and son in a consanguineous family; OMIM-confirmed |
| CASE0003 | TRIO / NIPBL / KMT2A | Review | Neurodevelopmental candidates matching phenotype; father-only duo limits segregation |
| CASE0067 | KMT2D | Review | Kabuki syndrome, Pathogenic, Exomiser top-ranked; singleton with weak phenotype |
| CASE0004 | None | Insufficient | Proband phenotype is absent (family-history context); correctly declined to avoid a false lead |

"Reportable" is the sprint's triage tier on synthetic cases ("put this in front of a clinical scientist
first"), not a statement that a finding is ready for a clinical report.

---

## Pipeline Architecture

In one line: **input → annotation → evidence → classification → ranking → case decision.**

```
 case variants (VCF / bundle)     clinical notes          family / pedigree
            │                           │                        │
            ▼                           ▼                        │
  1. Ingest + QC  (acmg.ingest)   2. HPO terms (acmg.hpo)        │
            │                           │                        │
            ▼                           │                        │
  3. Annotate: consequence, gene,       │                        │
     gnomAD AF, REVEL, SpliceAI (VEP)   │                        │
     + ClinVar, NMD, constraint         │                        │
     (acmg.vep_map, acmg.annotate,      │                        │
      acmg.clinvar, acmg.nmd, ...)      │                        │
            │                           │                        │
            ▼                           │                        │
  4. Evidence → ACMG class              │                        │
     SQL rules fire criteria,           │                        │
     points are summed (acmg.kernel,    │                        │
     acmg/manifests/*.sql)              │                        │
            │                           │                        │
            └──────────┬────────────────┘                        │
                       ▼                                         ▼
  5. Rerank: phenotype match (Monarch HPO→gene) + inheritance / segregation (acmg.rank, acmg.family)
                       │
                       ▼
  6. GO decision: Reportable / Review / Insufficient + rationale (acmg.decision)
                       │
                       ▼
  7. Literature: gather evidence, an LLM proposes, a human approves (acmg.evidence, acmg.agent)
```

The pipeline is organised as seven sequential stages, progressively narrowing millions of variants toward a small set of prioritised candidates.

```
1. Ingest        Parse the synthetic case dataset; normalise chromosomes, positions, and alleles.
2. HPO           Extract phenotype terms from clinical notes; apply negation and family-history filtering.
3. Annotate      Assign consequence, allele frequency, and clinical assertions to each variant.
4. ACMG          Classify variants through a transparent, rule-based kernel with ClinGen curation.
5. Rerank        Score candidates by phenotype match (information-content weighted, Monarch HPO to gene).
6. GO decision   Assign a case-level tier with an explainable rationale.
7. Literature    Gather supporting evidence for the leading candidate per case.
```

### Stage 1 - Ingest and quality control

The raw dataset was validated before any analysis. Checks confirmed zero null values across chromosome, position, reference, and alternate allele fields; all positions greater than zero; no records where the reference allele equalled the alternate allele; valid DNA bases only; and consistent chromosome naming (for example `chr1` normalised to `1`, `chrM` to `MT`). Deduplication collapsed 7.7 million variant rows to 430,284 unique variants.

Two dataset characteristics were treated as first-class considerations rather than footnotes. The cohort is predominantly of Middle-Eastern ancestry (Saudi, Egyptian, Iraqi, Yemeni), for which reference frequency resources such as gnomAD are less complete; some variants may therefore appear rarer than they are. Consanguinity was identified in 43 cases and factored into inheritance expectations.

### Stage 2 - Phenotype extraction

Phenotypes were extracted from clinical notes and normalised to Human Phenotype Ontology (HPO) terms. A rule-based extractor (FastHPOCR) was combined with a language-model validation layer that applies negation handling and, critically, family-history filtering.

The value of the validation layer is illustrated by CASE0004, where the clinical notes described findings that belonged to the patient's relatives rather than the proband. The validation layer correctly identified these as family-history context rather than proband-intrinsic phenotype and returned zero proband terms, preventing a false lead downstream. A comparison of rule-based versus language-model-assisted extraction is included as part of the validation work.

`acmg.hpo.extract_hpo(..., mode="tool_only")` is the deterministic, no-LLM path. The library function's default
mode calls an LLM through the `pi` CLI and so sends the note text to that provider; only do that with
de-identified text.

#### Privacy defaults

The per-case drivers (`examples/run_case.py`, `scripts/run_proband.py`, `scripts/run_all_probands.py`) are
safe by default (logic in `acmg/privacy.py`):

- **Input.** They read the de-identified `dataset.parquet`. `run_case.py` refuses the raw curated
  `dataset_clinical_curated.parquet`, which carries PII in `clinical_indication_text`, unless you pass
  `--allow-pii-input`. The variant calls are identical in both files.
- **External LLM.** HPO extraction defaults to `--hpo-mode tool_only` (local FastHPOCR, nothing leaves the
  machine). The LLM modes (`augment_select`, `candidates_model`, `model_only`) and the `run_all_probands.py
  --agent` arm send the clinical-indication text to an LLM provider through `pi`, and need
  `--allow-external-llm`. That flag cannot be combined with `--allow-pii-input`.
- A one-line warning is printed to stderr whenever either opt-in is used.

For example, to reproduce the LLM-assisted HPO extraction on the de-identified bundle:
`python examples/run_case.py CASE0003 --hpo-mode augment_select --allow-external-llm`.

### Stage 3 - Variant annotation

Variants were annotated for consequence, population allele frequency (gnomAD), known clinical assertions (ClinVar), and in-silico pathogenicity predictors (REVEL, SpliceAI). This provides the evidence base consumed by the classification kernel. Offline VEP annotation is identified as a production enhancement in the roadmap.

### Stage 4 - ACMG classification kernel

Variant classification is performed by a transparent, rule-based kernel implemented in SQL. Each classification is expressed as a specification that fires the relevant ACMG/AMP criteria, converts them to Tavtigian points, and derives a final class. Every call records exactly which criteria fired, and the kernel abstains (returning uncertain significance) when the evidence is insufficient, rather than over-calling.

This design was a deliberate choice in favour of traceability: unlike an opaque third-party classifier, every decision can be audited against its inputs. Concordance against ClinVar is reported in the Validation section.

The rules live in `acmg/manifests/01_acmg_rules.sql`, the thresholds and gene curation in `00_spec.sql`, and
the points combination in `02_combine.sql`. The kernel applies the computable criteria (PVS1, PS1, PM2, PM4,
PM5, PP2, PP3, BA1, BS1, BP4, BP7, and PM1 as an opt-in); the literature- and case-dependent ones (PS2, PS3,
PS4, PM3, PP1, PP4, ...) are left to the gated literature arm and the human reviewer.

### Stage 5 - Phenotype x genotype reranking

Candidates are reranked by how well the gene matches the patient's phenotype, using information-content-weighted scoring via Monarch (HPO to gene). Specific phenotype terms are weighted more heavily than generic ones. Where two candidates carry variants of equivalent variant-evidence strength, phenotype relevance and inheritance evidence break the tie.

### Stage 6 - Case-level prioritisation and tiering

The prioritisation score combines the four evidence types named in the challenge brief: phenotype match, inheritance-model compatibility, variant evidence strength, and database signals. Inheritance compatibility is derived from family segregation where parental data is available; database signals are derived from OMIM phenotype-relevance.

Segregation and OMIM function as evidence layers that elevate candidates when present, not as filters that every gene must pass. Each case receives a Reportable, Review, or Insufficient tier and a short rationale explaining why it should or should not be prioritised for expert re-review.

### Stage 7 - Literature support

For the leading candidate in each prioritised case, supporting literature is gathered to accompany the evidence summary. This arm operates under a propose-then-approve model, keeping a human reviewer in the loop.

---

## Data Schema

The core representation is a phenotype x genotype view per case, combining:

- Case metadata: case identifier, family structure, consanguinity status, ancestry, sex.
- Phenotype: HPO term set, extraction source, negation and family-history flags.
- Variant records: genomic coordinates, reference and alternate alleles, gene symbol, consequence, allele frequency, in-silico scores, ClinVar assertion, ACMG classification and fired criteria.
- Segregation: proband, paternal, and maternal genotypes where available; inheritance interpretation.
- Prioritisation: composite score, tier, and rationale.

Two structural details of the source data informed the implementation. The annotated files use the GRCh38 build while the raw genotype data uses GRCh37, requiring liftover for segregation; and the raw data is variant-only, so the absence of a record at a position denotes a homozygous-reference (0/0) genotype rather than missing data.

The kernel's own input is one flat table per variant, the `annotations` schema (`acmg.vep_map.REQUIRED_COLS`):
`variant_key` (`chrom-pos-ref-alt`), `gene`, `consequence`, `variant_kind`, `filtering_af`, `gnomad_mis_z`,
`revel`, `spliceai`, `clinvar_same_aa`, `clinvar_same_codon_lp`, `nmd_escaping`, `pm1`, `loeuf`, `cadd_phred`,
`clinvar_classification`. Any column may be missing or NULL; the rules that need it then abstain.
`acmg.vep_to_annotations` builds this table from a VEP output.

---

## Prioritisation Rubric

Cases are tiered as follows:

- **Reportable** - a strong candidate supported by converging evidence (for example phenotype match with confirming segregation and OMIM relevance, or confirmed recessive segregation).
- **Review** - credible phenotype-matched candidates, but with segregation unresolved (for example a duo lacking one parent, or a singleton without family data).
- **Insufficient** - no reportable candidate; for example when the proband has no recorded phenotype, so no candidate can be tied to the patient.

The rubric is deliberately conservative to avoid unvalidated diagnostic claims. A candidate gene scoring highly does not by itself make a case reportable; case-level confidence also depends on the supporting phenotype and family data. A variant classified as Pathogenic is not a diagnosis unless it also explains the patient's phenotype.

---

## Validation and Benchmarking

The pipeline was evaluated on four independent axes.

**Exomiser cross-validation.** Candidates were compared against Exomiser 14.0.0, an independent phenotype-driven prioritisation tool. Across the five phenotyped probands, our candidate genes appear within Exomiser's top five in four of five cases, including a direct top-ranked match on KMT2D (CASE0067). EBF3 and PACS1 (CASE0008) appear at Exomiser ranks two and five, and FBN1 (CASE0007) at rank two. For CASE0003, our neurodevelopmental candidates (NIPBL, KMT2A, NAXE) also appear in Exomiser's top five, though the leading gene differs. CASE0133 has no recorded phenotype and was therefore not Exomiser-scored; it is validated by segregation instead.

**ClinVar concordance.** The ACMG kernel was benchmarked against ClinVar classifications across more than 84,000 variants, achieving 96.3% overall concordance. The residual discordances were investigated: after excluding conflicting ClinVar entries, they are almost entirely common population variants (allele frequencies up to 65%) that ClinVar labels as risk factors or protective rather than Mendelian-pathogenic. The kernel correctly applies the stand-alone benign frequency rule (BA1) to these, demonstrating that it distinguishes rare-disease-causing variants from common modifiers.

How to read this number: PS1/PM5 draw on ClinVar itself, so the comparison is not fully independent, and an
overall concordance is dominated by whichever class is most common. `scripts/validate_kernel.py` (needs the
ClinVar VCF and `duckhts`) reproduces the comparison and prints its own caveats.

**Rule-based versus language-model HPO extraction.** Rule-based (FastHPOCR) and language-model-assisted extraction were compared for accuracy and completeness across the cohort.

**Reproducibility.** The pipeline is deterministic: the same inputs always produce the same structured outputs. Language-model steps are run at temperature zero with the model version logged.

---

## Where to get the data

None of this is needed for the [quick start](#quick-start-5-minutes-offline). `bash scripts/fetch_data.sh`
downloads the core public resources into `.cache/` (git-ignored), skips files already present, and writes
`.cache/MANIFEST.txt` with a checksum, size and date for each file. `docs/data_versions.md` records the exact
versions the sprint used. Everything is on the **GRCh37** build. Licences are summarised for convenience;
check each provider's current terms before use, especially for commercial work.

| Resource | Version used | Used for | Licence / terms | How to get it |
|----------|--------------|----------|-----------------|---------------|
| Abiomix challenge bundle (`dataset.parquet`) | July 2026 sprint release | All per-case drivers | Provided to sprint participants; not redistributable | From Abiomix. Use the de-identified `dataset.parquet`, not `dataset_clinical_curated.parquet` (see `scripts/prepare_bundle.py`) |
| ClinVar `variant_summary.txt.gz` | 2026-06-28 | PS1, PM5, ClinVar class, review stars | Public domain (NCBI) | `fetch_data.sh` (NCBI FTP) |
| ClinVar GRCh37 VCF | rolling | Conflict composition, validation harness | Public domain (NCBI) | `fetch_data.sh` |
| gnomAD gene constraint | v2.1.1 | PP2 (`mis_z`), LOEUF | CC0 | `fetch_data.sh` (Google Cloud public bucket) |
| gnomAD exome allele frequencies | v2.1.1 | PM2, BS1, BA1 | CC0 | Not bulk-downloaded: `scripts/build_site_scores.py` range-reads only the needed loci (network) |
| GENCODE basic GTF | v46lift37 | NMD rule, MANE transcript (PVS1) | Open, no restrictions (EMBL-EBI terms) | `fetch_data.sh` |
| REVEL | v1.3 | PP3 / BP4 for missense | Free for non-commercial use | Manual: download the zip from the REVEL site and convert to `revel_grch37.parquet` (`fetch_data.sh` prints the steps) |
| SpliceAI precomputed scores | Illumina release | PP3 / BP7 for splicing | Free for non-commercial use; Illumina BaseSpace login | Manual (or the VEP plugin / VEP REST for a small residual set) |
| CADD | v1.7 (exon slice) | Optional supporting fallback | Free for non-commercial use | `WITH_CADD_SLICE=1 bash scripts/fetch_data.sh` |
| Human Phenotype Ontology `hp.obo` | ~2026-07 | HPO extraction | HPO licence (free use with attribution) | `fetch_data.sh`, then `python scripts/build_hpo_index.py` for `hp.index` |
| ClinGen Gene-Disease Validity + Dosage Sensitivity | downloaded 2026-07-03 | PVS1 gate, validity cap, inheritance | CC0 (ClinGen) | `fetch_data.sh` (via `acmg.clingen.download`) |
| Monarch KG | `latest` (rolling) | Phenotype → gene ranking | Aggregated; mixed upstream licences | Attached remotely by `acmg.rank` (network); pin a snapshot URL for reproducibility |
| Ensembl VEP | 105 (bioconda) / 116 (Docker) cache; REST for residuals | Consequence, gene, transcript, AF | Apache-2.0 software; open data | Optional, large (~15-26 GB cache): `WITH_VEP=1 bash scripts/fetch_data.sh` prints the steps; `scripts/vep116.sh` |
| ClassifyCNV + bedtools | commit 148757c / 2.31.1 | CNV classification (Riggs 2020) | See upstream repositories | `bash scripts/setup_classifycnv.sh` |
| OpenSpliceAI models | 0.0.7 | Optional on-device splice scores | See upstream repository | `pip install -e '.[splice]'` + `scripts/setup_openspliceai.sh` (heavy, PyTorch) |
| Europe PMC, LitVar2, MARRVEL, DECIPHER APIs | live | Literature arm | Each service's terms | Called live by `acmg.evidence` (network) |

---

## Glossary

**Classification framework**

| Term | Meaning |
|------|---------|
| ACMG/AMP | The 2015 American College of Medical Genetics and Genomics / Association for Molecular Pathology guideline for classifying variants as Pathogenic, Likely Pathogenic, Uncertain Significance (VUS), Likely Benign or Benign by combining coded evidence criteria. |
| Criterion codes | First letters give direction and default strength: PVS = Pathogenic Very Strong, PS = Pathogenic Strong, PM = Pathogenic Moderate, PP = Pathogenic Supporting, BA = Benign stand-Alone, BS = Benign Strong, BP = Benign Supporting. A suffix such as `_Supporting` or `_Strong` means the criterion was applied at a different strength than its default. |
| Tavtigian points | A Bayesian points system (Tavtigian et al. 2018, 2020) used here to combine criteria: Supporting ±1, Moderate ±2, Strong ±4, Very Strong ±8. |
| Abstain | The kernel applies no criterion when the input it needs is missing (for example no allele frequency) rather than guessing. |
| VCEP | ClinGen Variant Curation Expert Panel; publishes gene-specific adjustments to the ACMG rules (`acmg_thresholds` holds such overrides as rows). |
| ClinGen SVI | ClinGen Sequence Variant Interpretation working group; its recommendations refine PVS1, PM2, PP3/BP4 and splicing. |

**Criteria used by the kernel**

| Code | Direction, points | What it means here |
|------|-------------------|--------------------|
| PVS1 | Pathogenic, +8 (+4 as `PVS1_Strong` if the transcript escapes NMD) | Predicted loss-of-function variant (stop-gained, frameshift, canonical splice site) in a gene where loss of function is an established disease mechanism (ClinGen haploinsufficiency score 3). |
| PS1 | Pathogenic, +4 | Same amino-acid change as an established pathogenic variant (ClinVar). |
| PM1 | Pathogenic, +2 (opt-in, off by default) | Missense in a mutational hotspot or critical domain. |
| PM2 (`PM2_Supporting`) | Pathogenic, +1 | Absent or very rare in population databases (gnomAD). Applied at Supporting strength, as ClinGen SVI recommends. |
| PM4 | Pathogenic, +2 | Protein length change: in-frame deletion or insertion, or stop-loss. |
| PM5 | Pathogenic, +2 | A different pathogenic missense change at the same amino-acid position (ClinVar). |
| PP2 | Pathogenic, +1 | Missense in a gene with low tolerance to missense variation (gnomAD missense Z ≥ 3.09). |
| PP3 | Pathogenic, +1 to +4 | Computational evidence of damage: REVEL for missense, SpliceAI for splicing (CADD only as an opt-in fallback). |
| BA1 | Benign, stand-alone | Allele frequency above 5% in a population: benign on its own (with an exception list for known common pathogenic variants such as HbS). |
| BS1 | Benign, -4 | Allele frequency higher than expected for the disorder (default > 1%). |
| BP4 | Benign, -1 to -8 | Computational evidence of no impact (low REVEL; CADD fallback). |
| BP7 | Benign, -1 | Synonymous change with no predicted splicing impact. |

**Criteria the kernel does not apply** (they need literature, functional data or the family): PS2 / PM6
(de novo), PS3 / BS3 (functional studies), PS4 (more frequent in affected people), PM3 (in trans with a
pathogenic variant, recessive), PP1 / BS4 (segregation with disease in the family), PP4 (phenotype highly
specific for the gene). The literature arm can propose some of these for a human to approve.

**Genomics terms**

| Term | Meaning |
|------|---------|
| Variant / SNV / indel | A difference from the reference genome. SNV = single-nucleotide variant; indel = small insertion or deletion. |
| CNV | Copy-number variant: a deletion or duplication of a larger stretch of DNA. Classified separately (Riggs et al. 2020), not by the SNV kernel. |
| GRCh37 / GRCh38, liftover | Two versions of the human reference genome. Coordinates differ between them; liftover converts between builds. |
| Consequence | Predicted effect of a variant on a transcript (missense, frameshift, splice donor, ...), in Sequence Ontology terms. |
| LoF | Loss of function: the variant is expected to stop the gene product working. |
| NMD | Nonsense-mediated decay: cells destroy mRNAs that carry an early stop codon. A premature stop in the last exon (or < 50 nt before the last exon junction) escapes NMD, so PVS1 is downgraded. |
| MANE Select | One agreed representative transcript per gene (NCBI + Ensembl); the pipeline reports consequences on it. |
| VEP | Ensembl Variant Effect Predictor: annotates consequence, gene, transcript and population frequency for each variant. |
| gnomAD | Genome Aggregation Database: allele frequencies from large population cohorts, plus per-gene constraint scores. |
| Filtering allele frequency | The frequency used for PM2/BS1/BA1, taken as the maximum across populations so a variant common in one ancestry is not diluted by an average. |
| ClinVar | NCBI public archive of variant classifications submitted by laboratories and expert panels. |
| ClinGen | Clinical Genome Resource: curates gene-disease validity (Definitive ... Refuted) and dosage sensitivity (haploinsufficiency / triplosensitivity scores). |
| REVEL, SpliceAI, CADD | Computational predictors: REVEL scores missense pathogenicity, SpliceAI predicts splicing disruption, CADD is a genome-wide deleteriousness score. |
| HPO | Human Phenotype Ontology: standard terms for clinical features (for example HP:0001250 Seizure). |
| Monarch KG | Knowledge graph linking HPO phenotypes to genes and diseases, used for the phenotype-match score. |
| OMIM | Catalogue of human genes and genetic disorders. |
| Exomiser | Independent open-source tool that ranks variants by phenotype and genotype; used here as a cross-check. |

**Family and inheritance terms**

| Term | Meaning |
|------|---------|
| Proband | The affected person through whom the family came to attention (the "patient" of the case). |
| Singleton / duo / trio | Only the proband sequenced / proband plus one parent / proband plus both parents. |
| Segregation | Whether a variant tracks with the disease through the family (affected relatives carry it, unaffected ones do not). |
| De novo | A variant present in the child but in neither parent; strong evidence for dominant disorders. |
| AD / AR (MOI) | Mode of inheritance: autosomal dominant (one altered copy is enough) / autosomal recessive (both copies). |
| Homozygous / heterozygous / compound heterozygous | Same variant on both copies / one copy / two different variants, one on each copy, in the same gene. |
| Consanguinity | Parents are related; raises the chance of homozygous recessive variants. |
| Penetrance | The share of carriers who actually show the disorder. |

---

## Installation

Requires Python 3.10 or newer (CI runs 3.10 and 3.12).

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e '.[dev]'           # core (duckdb, pandas, requests) + test tooling (pytest, pyarrow, openpyxl)
pip install FastHPOCR pronto      # optional: HPO concept recognition (Stage 2); pronto only builds hp.index
# pip install -e '.[splice]'      # optional, heavy (PyTorch): on-device OpenSpliceAI lane
```

**Exact pinned versions.** `requirements.lock` pins the core dependencies (generated with `uv`). It pins
numpy 2.5 and pandas 3.0, which need **Python 3.12 or newer**:

```bash
pip install -r requirements.lock && pip install -e . --no-deps
```

Reference data (ClinVar, ClinGen, gnomAD constraint, GENCODE, ...) is fetched separately with
`scripts/fetch_data.sh` (see [Where to get the data](#where-to-get-the-data)); record the versions you use in
`docs/data_versions.md`.

## Running the Tests

```bash
python -m pytest -q    # offline; add -rs to see why tests were skipped
ruff check .           # lint (pip install ruff)
```

Tests that need optional resources skip with a reason instead of failing: the VCF-ingest tests need the
`duckhts` DuckDB community extension (downloaded on first use), and the HPO-mode tests need a built FastHPOCR
`hp.index`. CI (`.github/workflows/ci.yml`) runs ruff, pytest and the offline demo on every push and pull
request, and checks that `requirements.lock` installs.

---

## Troubleshooting

| Symptom | Cause and fix |
|---------|---------------|
| `ModuleNotFoundError: No module named 'acmg'` | The package is not installed in the active environment. Activate the venv and run `pip install -e .` from the repository root. |
| `ModuleNotFoundError: No module named 'scripts'` | Run drivers from the repository root with `PYTHONPATH=.`, for example `PYTHONPATH=. python scripts/run_all_probands.py --help`. |
| `No matching distribution found for numpy==2.5.0` | `requirements.lock` needs Python 3.12+. Use a newer Python, or `pip install -e '.[dev]'`, which picks versions that suit your Python. |
| `SKIPPED ... duckhts extension unavailable` | Expected offline or behind a proxy that blocks `community-extensions.duckdb.org`. The extension downloads automatically on first use with network. |
| `SKIPPED ... FastHPOCR hp.index not built` | Expected until you build the index: download `hp.obo`, `pip install FastHPOCR pronto`, `python scripts/build_hpo_index.py`. |
| `No module named 'pronto'` while building `hp.index` | FastHPOCR's indexer needs it but does not declare it: `pip install pronto`. |
| `run_case.py needs data that is not here yet` | It needs the challenge bundle and `.cache/` reference data; the message lists each missing file. For an offline run use `examples/demo.py`. |
| `pi CLI not found on PATH` | The literature arm, the LLM HPO modes and the `--agent` arm call an LLM through the `pi` CLI. Install and configure it, pass your own `prompt -> text` callable as `llm=`, or keep the drivers' default `--hpo-mode tool_only`. |
| `refusing to read ...dataset_clinical_curated.parquet` | That file carries raw PII. Use the de-identified `dataset.parquet` (`--bundle /path/to/dataset.parquet`); pass `--allow-pii-input` only if you really mean to read it. |
| `... sends the case's clinical-indication text to an external LLM provider` | An LLM HPO mode or `--agent` was requested without consent. Add `--allow-external-llm` (de-identified text only) or drop back to `--hpo-mode tool_only`. |
| `(Monarch unavailable: ...)`, phenotype score 0 | The Monarch KG is attached over HTTPS; check network access to `data.monarchinitiative.org`. |
| `vep_to_annotations: input has no column(s) ['POS']` | Your VEP table names coordinates differently. Pass `cols={'chrom': ..., 'pos': ..., 'ref': ..., 'alt': ...}`. |
| PVS1 never fires | PVS1 needs both gene curation saying LoF is a disease mechanism (`gene_curation=`, for example from `acmg.clingen`) and an NMD status (`nmd_escaping` 0 or 1). Without either it abstains by design. |
| PM2 never fires | PM2 abstains when the allele frequency is missing. Provide `filtering_af` (gnomAD) for each variant. |
| Everything comes out VUS | Usually missing evidence rather than a bug: check `criteria` for what fired and fill the inputs the abstaining rules need. |

---

## Repository Structure

```
.
├── acmg/                   Python package: SQL ACMG kernel (acmg/manifests/*.sql), ingest, family/segregation,
│                           HPO extraction, ClinGen/ClinVar/gnomAD evidence, NMD, CNV, reranking, GO decision,
│                           and the gated literature/agent arm
├── scripts/                Cohort-scale and per-proband drivers (prepare_bundle, annotate_cohort, run_bundle,
│                           run_proband, run_all_probands, validate_kernel, ...) and data/tool setup scripts
│                           (fetch_data.sh, build_hpo_index.py, setup_*.sh)
├── examples/               demo.py (offline), run_case.py and literature_arm.py (need network / local data)
├── tests/                  pytest suite with small synthetic fixtures (no network, no large data)
├── docs/                   data_versions.md (reference-data and tool versions), scale.md (cohort-scale notes)
├── data/, outputs/, notebooks/   Local working directories (large data is git-ignored)
├── pyproject.toml          Package metadata and optional extras (dev, splice)
├── requirements.txt / requirements.lock   Core dependencies (lock generated with uv)
└── README.md
```

---

## Reproducibility

The pipeline is designed so that the same inputs always generate the same structured outputs. Deterministic classification rules, logged tool versions, and temperature-zero language-model settings ensure that results can be regenerated and audited.

In practice, record: `requirements.lock` (Python packages), the reference-data versions in
`docs/data_versions.md` and `.cache/MANIFEST.txt`, the VEP release, and the git commit (which versions the SQL
rules in `acmg/manifests/`). Rolling sources (ClinVar, ClinGen, Monarch `latest`) change classifications over
time, which is the point of re-analysis, so note their download dates.

---

## Limitations

- Annotation used the GeneBe API (VEP-equivalent); offline VEP is preferred for production throughput and to un-abstain rules requiring loss-of-function context (PS1, PM5, NMD).
- gnomAD under-represents Middle-Eastern ancestry, so rarity estimates for this cohort should be interpreted with caution.
- Certain in-silico scores (for example REVEL, CADD) are subject to commercial-use restrictions.
- The reference-build mismatch between annotated (GRCh38) and raw (GRCh37) data requires liftover, and a minority of variant positions do not resolve cleanly.
- The kernel covers only the computable ACMG/AMP criteria, so it structurally under-calls compared with full expert curation, and its default frequency thresholds are general ones (for example PM2 at ≤ 0.5%), not disease-specific.
- The results are a proof-of-concept demonstration on synthetic data. Expert clinical sign-off is required before any real-world application.

---

## Roadmap

- Offline VEP annotation to strengthen loss-of-function classification rules.
- Zygosity-aware scoring tailored to consanguineous cohorts.
- Build harmonisation to improve positional segregation matching.
- A full adjudication agent operating under a propose-then-approve model with human oversight.
- Prior-analysis diffing to surface only newly reportable findings on re-analysis.

---

## Ethical and Data-Safety Statement

This project used only synthetic case data provided by Abiomix. No real patient data was used at any stage. All outputs are traceable to non-sensitive inputs. Language-model output is treated as decision support and never as clinical decision-making. Expert clinical sign-off remains essential before any real-world use.

---

## Licence and citation

See `LICENSE` for the code licence. Reference data keeps its own licence (see
[Where to get the data](#where-to-get-the-data)).

If you use or build on this work, please cite the repository and the guidelines it implements:

- Halireena (Rush), Nubla, Rayane. *Re-Analysing the Unsolved: a transparent, reproducible pipeline for
  prioritising unsolved rare-disease cases for genomic re-analysis.* BioConnect 2026 Consultancy Sprint
  (Abiomix challenge), University of Birmingham Dubai, 2026.
  https://github.com/halireena/abiomix-rare-disease-prioritisation
- Richards S, et al. Standards and guidelines for the interpretation of sequence variants. *Genet Med* 2015.
- Tavtigian SV, et al. Fitting a naturally scaled point system to the ACMG/AMP variant classification
  guidelines. *Hum Mutat* 2020.
- Abou Tayoun AN, et al. Recommendations for interpreting the loss of function PVS1 ACMG/AMP variant criterion.
  *Hum Mutat* 2018.
- Pejaver V, et al. Calibration of computational tools for missense variant pathogenicity classification and
  ClinGen recommendations for PP3/BP4 criteria. *Am J Hum Genet* 2022.
- Riggs ER, et al. Technical standards for the interpretation and reporting of constitutional copy-number
  variants. *Genet Med* 2020.
