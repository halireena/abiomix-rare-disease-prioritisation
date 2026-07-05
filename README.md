# Re-Analysing the Unsolved

### A transparent, reproducible pipeline for prioritising unsolved rare-disease cases for genomic re-analysis

**BioConnect 2026 Consultancy Sprint — Abiomix Challenge**
Variant Re-Analysis and Case Re-Classification at Scale

Team: Halireena (Rush), Nubla, Rayane
University of Birmingham Dubai · 2–5 July 2026

---

## Overview

Around half of rare-disease patients remain without a molecular diagnosis after their initial genomic analysis. Clinical genomics, however, is not static: new gene–disease associations, updated variant databases, refined phenotype ontologies, and improved computational tools mean that cases which were unsolvable a few years ago may be solvable today. Published studies show that systematic re-analysis increases diagnostic yield by 10–15%.

The practical bottleneck is scale. When a laboratory holds hundreds or thousands of historical cases, it must decide which ones justify limited expert time. This project delivers a proof-of-concept pipeline that converts historical case information into a structured phenotype × genotype representation and produces an explainable, case-level prioritisation — without bypassing clinical oversight or making unvalidated diagnostic claims.

All work was carried out on a synthetic dataset provided by Abiomix. No real patient data was used at any stage. All outputs are traceable to non-sensitive inputs, and any language-model output is treated strictly as decision support, not clinical decision-making.

---

## Key Results

- 138 de-identified cases across 119 families processed end to end.
- 7.7 million variant rows reduced to 430,284 unique variants after cleaning and deduplication (a 93% reduction).
- A transparent, rule-based ACMG classification kernel achieving 96.3% concordance with ClinVar across more than 84,000 variants.
- Independent cross-validation against Exomiser: agreement on the leading candidate gene in four of five phenotyped probands.
- Six probands examined in depth, each assigned a Reportable, Review, or Insufficient tier with a documented rationale.

Representative findings from the deep-dive probands:

| Case | Lead candidate | Tier | Basis |
|------|----------------|------|-------|
| CASE0133 | CCNO | Reportable | Homozygous proband, both parents carriers; textbook autosomal-recessive segregation confirmed against source data |
| CASE0008 | EBF3 / PACS1 | Reportable | Phenotype-confirmed; EBF3 re-prioritised from a variant of uncertain significance and corroborated by Exomiser |
| CASE0007 | CUL7 / COL27A1 / MACF1 | Reportable | Co-segregation with affected father and son in a consanguineous family; OMIM-confirmed |
| CASE0003 | TRIO / NIPBL / KMT2A | Review | Neurodevelopmental candidates matching phenotype; father-only duo limits segregation |
| CASE0067 | KMT2D | Review | Kabuki syndrome, Pathogenic, Exomiser top-ranked; singleton with weak phenotype |
| CASE0004 | None | Insufficient | Phenotype belongs to siblings (family history); correctly declined to avoid a false lead |

---

## Pipeline Architecture

The pipeline is organised as seven sequential stages, progressively narrowing millions of variants toward a small set of prioritised candidates.
1. Ingest        Parse the synthetic case dataset; normalise chromosomes, positions, and alleles.
2. HPO           Extract phenotype terms from clinical notes; apply negation and family-history filtering.
3. Annotate      Assign consequence, allele frequency, and clinical assertions to each variant.
4. ACMG          Classify variants through a transparent, rule-based kernel with ClinGen curation.
5. Rerank        Score candidates by phenotype match (IC-weighted, Monarch HPO to gene).
6. GO decision   Assign a case-level tier with an explainable rationale.
7. Literature    Gather supporting evidence for the leading candidate per case.


### Stage 1 — Ingest and quality control

The raw dataset was validated before any analysis. Checks confirmed zero null values across chromosome, position, reference, and alternate allele fields; all positions greater than zero; no records where the reference allele equalled the alternate allele; valid DNA bases only; and consistent chromosome naming (for example `chr1` normalised to `1`, `chrM` to `MT`). Deduplication collapsed 7.7 million variant rows to 430,284 unique variants.

Two dataset characteristics were treated as first-class considerations rather than footnotes. The cohort is predominantly of Middle-Eastern ancestry (Saudi, Egyptian, Iraqi, Yemeni), for which reference frequency resources such as gnomAD are less complete; some variants may therefore appear rarer than they are. Consanguinity was identified in 43 cases and factored into inheritance expectations.

### Stage 2 — Phenotype extraction

Phenotypes were extracted from clinical notes and normalised to Human Phenotype Ontology (HPO) terms. A rule-based extractor (FastHPOCR) was combined with a language-model validation layer that applies negation handling and, critically, family-history filtering.

The value of the validation layer is illustrated by CASE0004, where the clinical notes described cardiac findings that in fact belonged to the proband's siblings. The validation layer correctly identified these as family history rather than proband-intrinsic phenotype and returned zero proband terms, preventing a false lead downstream. A comparison of rule-based versus language-model-assisted extraction is included as part of the validation work.

### Stage 3 — Variant annotation

Variants were annotated for consequence, population allele frequency (gnomAD), known clinical assertions (ClinVar), and in-silico pathogenicity predictors (REVEL, SpliceAI). This provides the evidence base consumed by the classification kernel. Offline VEP annotation is identified as a production enhancement in the roadmap (see Limitations).

### Stage 4 — ACMG classification kernel

Variant classification is performed by a transparent, rule-based kernel implemented in SQL (`acmg/manifests/`, run under DuckDB). Each classification is expressed as a specification that fires the relevant ACMG/AMP criteria, converts them to Tavtigian points, and derives a final class. Every call records exactly which criteria fired, and the kernel abstains (returning uncertain significance) when the evidence is insufficient, rather than over-calling.

This design was a deliberate choice in favour of traceability: unlike an opaque third-party classifier, every decision can be audited against its inputs. Concordance against ClinVar is reported in the Validation section.

### Stage 5 — Phenotype × genotype reranking

Candidates are reranked by how well the gene matches the patient's phenotype, using information-content-weighted scoring via Monarch (HPO to gene). Specific phenotype terms are weighted more heavily than generic ones. Where two candidates carry equally damaging variants, phenotype relevance breaks the tie.

### Stage 6 — Case-level prioritisation and tiering

The prioritisation score combines the four evidence types named in the challenge brief: phenotype match, inheritance-model compatibility, variant evidence strength, and database signals. Inheritance compatibility is derived from family segregation where parental data is available; database signals are derived from OMIM phenotype-relevance.

Segregation and OMIM function as evidence layers that elevate candidates when present, not as filters that every gene must pass. Each case receives a Reportable, Review, or Insufficient tier and a short rationale explaining why it should or should not be prioritised for expert re-review.

### Stage 7 — Literature support

For the leading candidate in each prioritised case, supporting literature is gathered to accompany the evidence summary. This arm operates under a propose-then-approve model, keeping a human reviewer in the loop.

---

## Data Schema

The core representation is a phenotype × genotype view per case, combining:

- Case metadata: case identifier, family structure, consanguinity status, ancestry, sex.
- Phenotype: HPO term set, extraction source, negation and family-history flags.
- Variant records: genomic coordinates, reference and alternate alleles, gene symbol, consequence, allele frequency, in-silico scores, ClinVar assertion, ACMG classification and fired criteria.
- Segregation: proband, paternal, and maternal genotypes where available; inheritance interpretation.
- Prioritisation: composite score, tier, and rationale.

Two structural details of the source data informed the implementation. The annotated files use the GRCh38 build while the raw genotype data uses GRCh37, requiring liftover for segregation; and the raw data is variant-only, so the absence of a record at a position denotes a homozygous-reference (0/0) genotype rather than missing data.

---

## Prioritisation Rubric

Cases are tiered as follows:

- **Reportable** — a strong candidate supported by converging evidence (for example phenotype match with confirming segregation and OMIM relevance, or confirmed recessive segregation).
- **Review** — credible phenotype-matched candidates, but with segregation unresolved (for example a duo lacking one parent, or a singleton without family data).
- **Insufficient** — no reportable candidate; for example when the recorded phenotype is family history rather than proband-intrinsic, or when a singleton has only non-specific phenotype.

The rubric is deliberately conservative to avoid unvalidated diagnostic claims. A candidate gene scoring highly does not by itself make a case reportable; case-level confidence also depends on the supporting phenotype and family data.

---

## Validation and Benchmarking

The pipeline was evaluated on four independent axes.

**Exomiser cross-validation.** Candidates were compared against Exomiser 14.0.0, an independent phenotype-driven prioritisation tool. Agreement on the leading candidate gene was observed in four of five phenotyped probands, including a direct match on KMT2D (CASE0067) and top-five agreement on EBF3 and PACS1 (CASE0008) and FBN1 (CASE0007). CASE0133 has no recorded phenotype and was therefore not Exomiser-scored; it is validated by segregation instead.

**ClinVar concordance.** The ACMG kernel was benchmarked against ClinVar classifications across more than 84,000 variants, achieving 96.3% overall concordance. The residual discordances were investigated: after excluding conflicting ClinVar entries, they are almost entirely common population variants (allele frequencies up to 65%) that ClinVar labels as risk factors or protective rather than Mendelian-pathogenic. The kernel correctly applies the stand-alone benign frequency rule (BA1) to these, demonstrating that it distinguishes rare-disease-causing variants from common modifiers.

**Rule-based versus language-model HPO extraction.** Rule-based (FastHPOCR) and language-model-assisted extraction were compared for accuracy and completeness across the cohort.

**Reproducibility.** The pipeline is deterministic: the same inputs always produce the same structured outputs. Language-model steps are run at temperature zero with the model version logged.

---

## Repository Structure

```
.
├── acmg/                   Rule-based ACMG classification kernel (SQL manifests + Python: kernel,
│                           case-level fold, ranking, family segregation, HPO, ingest, evidence,
│                           agent/REPL layers)
├── scripts/                Pipeline entry points (annotate_cohort, annotate_cnv, run_proband,
│                           run_all_probands, run_bundle, validate_kernel, data-fetch/setup scripts)
├── tests/                  Unit tests (59 passing) covering the kernel, family segregation, CNV
│                           handling, ranking, and the case REPL, plus fixtures (trio VCF/PED, etc.)
├── examples/                Runnable examples, including a zero-setup synthetic demo (demo.py)
├── docs/                    Data-version and scale notes
├── data/                    (gitignored contents) local input bundles — not committed
├── notebooks/               Exploratory notebooks
├── outputs/                 Run logs (aggregate counts only — no patient-level data)
└── README.md
```
Real patient-level data, VCFs, reference genomes, and per-case results are intentionally **not** committed to this repository (see Ethical and Data-Safety Statement). The pipeline is fully open; the data it runs on stays local.

---

## Install & Run

```sh
git clone https://github.com/halireena/abiomix-rare-disease-prioritisation.git
cd abiomix-rare-disease-prioritisation

pip install -r requirements.lock   # exact pinned versions
pip install -e .

python -m pytest -q                # 59 tests, no external data required
python examples/demo.py            # end-to-end run on synthetic variants, no setup needed
```

`examples/demo.py` runs the full ACMG classification kernel on synthetic variants and prints each classification alongside the exact criteria that fired and the resulting point total — the fastest way to see the "glass box" reasoning without any real data or external services.

To run against a real cohort, supply your own tidy variant bundle (see `acmg/ingest.py`) and reference data via `scripts/fetch_data.sh` (ClinVar, GENCODE GTF, gnomAD constraint), then use `scripts/run_proband.py` or `scripts/run_all_probands.py`. None of this reference or patient data is stored in the repository.

---

## Limitations

- Annotation for this sprint used the GeneBe API (VEP-equivalent) as a fast path to results; offline VEP is preferred for production throughput, to avoid sending variants to an external service, and to un-abstain rules requiring loss-of-function context (PS1, PM5, NMD). The `acmg` package itself supports fully offline VEP annotation (`scripts/vep116.sh`, `acmg/vep_map.py`).
- gnomAD under-represents Middle-Eastern ancestry, so rarity estimates for this cohort should be interpreted with caution.
- Certain in-silico scores (for example REVEL, CADD) are subject to commercial-use restrictions.
- The reference-build mismatch between annotated (GRCh38) and raw (GRCh37) data requires liftover, and a minority of variant positions do not resolve cleanly.
- The results are a proof-of-concept demonstration on synthetic data. Expert clinical sign-off is required before any real-world application.

---

## Roadmap

- Offline VEP annotation as the default path (removing reliance on any external annotation API).
- Zygosity-aware scoring tailored to consanguineous cohorts.
- Build harmonisation to improve positional segregation matching.
- A full adjudication agent operating under a propose-then-approve model with human oversight.
- Prior-analysis diffing to surface only newly reportable findings on re-analysis.

---

## Ethical and Data-Safety Statement

This project used only synthetic case data provided by Abiomix. No real patient data was used at any stage. All outputs are traceable to non-sensitive inputs. Language-model output is treated as decision support and never as clinical decision-making. Expert clinical sign-off remains essential before any real-world use. Real patient data, if ever used with this pipeline, must never be committed to version control; the `.gitignore` in this repository excludes `.cache/`, `*.parquet`, `*.csv`, `*.vcf`-adjacent large data, and reference genomes by design.
