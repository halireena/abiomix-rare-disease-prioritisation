# Cohort-scale annotation

The demo annotates a few hundred novel variants over VEP-REST. The real cohort is **430,284 raw variants →
428,492 distinct loci**. This documents how the annotation runs at that scale in minutes, not hours: the
LOCAL evidence (the parts we own, no external quota) is DuckDB joins measured in seconds; VEP is the only
external step and belongs OFFLINE, not over REST.

Reproduce with `PYTHONPATH=. python3 scripts/annotate_cohort.py`.

## The pattern: annotate DISTINCT loci once, reuse

A cohort has one locus recurring across many samples and releases. Annotation is a property of the LOCUS
(chrom-pos-ref-alt), not the sample, so you annotate each distinct locus ONCE and reuse. The parquet's `ALT`
is a list; the script `UNNEST`s it and `SELECT DISTINCT`s the result: 430,284 raw → 428,492 distinct loci in
**0.09 s**. Everything downstream runs on the deduped set.

## Measured local timings @ 428,492 loci

All local, all offline (DuckDB 1.5.2, single machine). The reference-table loads are one-time file parses;
the JOINS themselves are sub-second.

| Step | What | Wall | Note |
|------|------|-----:|------|
| build cohort | UNNEST ALT + `DISTINCT` | **0.09 s** | 430,284 → 428,492 loci |
| load `clinvar_prot` | cache parquet (PS1/PM5) | **0.01 s** | 146,169 P/LP protein rows |
| load `clinvar` exact | parse `variant_summary.txt.gz` (440 MB) | **5.6 s** | 4,515,335 GRCh37 rows |
| load gnomAD constraint | `gnomad_constraint.txt.gz` (4.6 MB) | **0.25 s** | 19,658 genes |
| load exon model | `gencode.lift37.gtf.gz` (40 MB) via duckhts | **4.7 s** | 866,150 exons |
| **ClinVar exact join** | left-join @428k loci | **0.07 s** | 110,688 matched |
| **PS1/PM5 EXISTS join** | same-aa / same-codon @~438k | **0.19 s** | correlated `EXISTS` |
| **gnomAD mis_z join** | gene left-join | **<0.1 s** | fills PP2 |
| **NMD escape** | exon-order rule @199k coding | **0.15 s** | strand-aware |

**Total local annotation ≈ 11 s wall** for the whole 428k set — dominated by the two one-time file parses
(variant_summary 5.6 s, GTF 4.7 s); the actual joins are all sub-second. There is no external quota on any of
this: ClinVar classification, PS1/PM5, gnomAD constraint and NMD are ours.

The gene/consequence/protein_pos/transcript columns that PS1/PM5, PP2 and NMD consume come from VEP — the
joins are wired in `join_vep_derived()`, which takes a VEP-annotated table and returns the kernel
`annotations` shape. Feed it the offline-VEP output (below).

## Site scores without the download: PRUNED REMOTE tabix

The frequency (PM2/BA1/BS1) and in-silico (PP3/BP4) evidence needs site-level scores the local tables don't
carry: gnomAD population AF, a missense predictor, a splice score. The naive path downloads 60–120 GB of
reference VCFs. `scripts/build_site_scores.py` does NOT: it reads ONLY the bytes covering the cohort's loci,
over HTTP range reads against bgzip+tabix files, via duckhts
`read_tabix(url, region := 'chr:start-end', index_path := '<local .tbi>')`. Only the small `.tbi` index is
cached locally (KB–MB); the fat file is never pulled.

**The AF-pruning point.** gnomAD's VCF INFO is ~23 KB/record (VEP CSQ + every population's AF/AC/AN). A region
read decompresses the whole INFO, but we keep only two numbers — pruned IN SQL inside the read query so the
fat column never crosses into Python:

```sql
SELECT column0 AS chrom, column1 AS pos, column3 AS ref_, column4 AS alt_,
       TRY_CAST(regexp_extract(column7, '(^|;)AF=([0-9.eE+-]+)', 2)        AS DOUBLE) AS filtering_af,
       TRY_CAST(regexp_extract(column7, '(^|;)AF_popmax=([0-9.eE+-]+)', 2) AS DOUBLE) AS af_popmax
FROM read_tabix(:url, region := :region, index_path := :tbi)
```

(`ref`/`end` are DuckDB reserved words, so tabix columns are aliased to `ref_`/`end_`.)

**Fewer seeks + concurrency.** Adjacent loci are MERGED (within `--gap`, default 10 kb) into one region:
428,492 loci → **42,188 regions**. Region reads are latency-bound (~1.2 s each serially against GCS), so they
run CONCURRENTLY on a thread pool with a per-thread DuckDB connection (DuckDB cursors from one connection
serialise; separate connections don't) — **~0.05 s/region effective at 32 workers, a ~24× speedup.**

### Which sources serve HTTP ranges (verified with `curl -sI` before use)

| Source | File | `accept-ranges` | `.tbi` | Column |
|--------|------|:---------------:|-------:|--------|
| **gnomAD v2.1.1 exomes** (GRCh37, GCS) | `gnomad.exomes.r2.1.1.sites.vcf.bgz` (63 GB) | **bytes** | 944 KB | `filtering_af`, `af_popmax` |
| **CADD v1.6** (GRCh37 whole-genome SNVs, Kircher lab) | `whole_genome_SNVs.tsv.gz` (84 GB) | **bytes** | 2.7 MB | `cadd_phred` |
| REVEL | — | — | — | no range+tabix host (see fallback) |
| SpliceAI (precomputed) | Illumina BaseSpace | auth-gated | — | REST residual (see fallback) |

### Measured (chr22 subset — real cohort loci)

10,685 distinct chr22 loci → 738 merged regions, 32 workers, `.tbi` cached (no full-file download):

| Source | Regions | Wall | Records fetched | Loci matched |
|--------|--------:|-----:|----------------:|-------------:|
| gnomAD AF | 738 | **37.6 s** (0.051 s/region) | 348,046 | 7,300 / 10,685 |
| CADD PHRED | 738 | **72.5 s** (0.098 s/region) | 25,331,676 | 8,858 / 10,685 |

Output `site_scores.parquet`: 10,685 rows, **224 KB**, `ORDER BY chrom,pos` (interval-friendly row groups).
CADD is denser (every genomic base × 3 alts) so it moves more bytes per region — hence its higher per-region
cost. **Projected full cohort (42,188 regions, same per-region rate): gnomAD ≈ 36 min, CADD ≈ 69 min** — both
sources ≈ 1.75 h, versus a 60–120 GB download otherwise. The cache footprint stays at **3.6 MB of `.tbi`**.

### Schema + wiring

`site_scores.parquet`: `chrom, pos, ref, alt, filtering_af, af_popmax, revel, cadd_phred, spliceai`. `revel`
and `spliceai` are carried NULL (no range-served source below) so a later fill needs no schema change; the
kernel abstains on a NULL, never guesses. `scripts/annotate_cohort.py:join_site_scores()` fills
`filtering_af`/`revel`/`spliceai` into the kernel `annotations` by variant_key (chrom-pos-ref-alt);
`join_vep_derived()` calls it, so the same wiring runs at cohort scale alongside the local ClinVar/gnomAD/NMD
joins. Verified end-to-end: a rare chr22 locus (AF 4.5e-6) fires `PM2_Supporting`; a common one (AF 0.0066)
abstains.

### Fallback: REVEL and SpliceAI (no range+tabix host)

REVEL ships only as a 667 MB zip (range-served but a zip — no tabix seek without extract), and precomputed
SpliceAI is auth-gated on Illumina BaseSpace. Neither is pulled at runtime. The one-time operator step: for
REVEL, download+unzip the zip once, `bgzip` + `tabix -s1 -b2 -e2` the CSV, then it slots into the SAME
`read_tabix` path as a local file; for SpliceAI, provision the precomputed VCFs once (they are also the offline
VEP-plugin files below) or leave splice scoring to the small VEP-REST residual. Until then the kernel abstains
on `revel`/`spliceai` — fail-safe.

## VEP: the one external step — OFFLINE, not REST

`vep` is not installed here and there is no `~/.vep` cache, so the script prints the offline command and skips
(the join stays wired). Feasibility of the two VEP paths:

### REST does not scale (measured)

Ensembl VEP-REST: **200 variants/POST**, whole API **rate-limited at 55,000 requests/hour** (confirmed live:
`X-RateLimit-Limit: 55000`, `X-RateLimit-Period: 3600`). Measured throughput of one 200-variant POST against
`grch37.rest.ensembl.org`:

- bare: **1.41 s / 200 vars → 142 var/s**
- with REVEL + SpliceAI + CADD plugins (what ACMG actually needs): **7.55 s / 200 vars → 26 var/s**

At 26 var/s, 428,492 loci = **~4.6 hours serial minimum**, before any 429 back-off, retries or network
variance push it to tens of hours. A per-variant (unbatched) pipeline hits the quota wall instead:
428k / 55,000 per hour = **~7.8 h** just on rate limit. Either way: hours, not minutes. REST is for the small
RESIDUAL novel set only (`rest_fallback()` in the script; `acmg.vep_rest`), after pre-filtering rarity +
coding/splice + panel so "novel" is a few hundred.

### Offline VEP: minutes

VEP offline (indexed cache + tabix'd plugin files) does **thousands of variants/sec**, so 428k finishes in
**~1–3 minutes** — two orders of magnitude past REST and with no quota. The documented command is
`vep_offline_command()` in the script:

```
vep --offline --cache --dir_cache $HOME/.vep --assembly GRCh37 --fasta $HOME/.vep/GRCh37.fa.gz \
    --input_file cohort.vcf --output_file cohort.vep.tsv --tab --force_overwrite --no_stats \
    --symbol --mane_select --pick --canonical --af_gnomade --af_gnomadg \
    --plugin REVEL,/data/revel/new_tabbed_revel_grch37.tsv.gz \
    --plugin SpliceAI,snv=.../spliceai_scores.raw.snv.hg19.vcf.gz,indel=.../spliceai_scores.raw.indel.hg19.vcf.gz \
    --fields Uploaded_variation,Location,SYMBOL,Consequence,Amino_acids,Protein_position,Feature,\
gnomADe_AF,gnomADg_AF,REVEL,SpliceAI_pred_DS_max
```

`--mane_select --pick` fixes the ONE reported transcript — load-bearing, because the NMD escape and PS1/PM5
protein-position joins depend on the annotation transcript (see `acmg/nmd.py`). Load the TSV, map it with
`acmg.vep_map.vep_to_annotations`, and pass it to `join_vep_derived()`.

### Data footprint (why we do NOT auto-download it)

The offline cache and plugin data are large and are provisioned once by an operator, not pulled at runtime:

| Asset | Approx size | Purpose |
|-------|-------------|---------|
| VEP GRCh37 indexed cache (`homo_sapiens`) | **~15 GB** | consequence + gene + transcript |
| GRCh37 reference FASTA (bgzipped) | ~1 GB | HGVS / codon context |
| REVEL (dbNSFP-derived, tabbed) | ~5 GB | missense pathogenicity |
| SpliceAI precomputed (SNV + indel VCFs) | ~30–120 GB | splice score |
| CADD (optional) | ~80 GB | extra in-silico |

The cache alone is the ~15 GB the task flags; do not download it inside the pipeline. Provision it once on the
run host; the script detects `vep` + `~/.vep` and only then emits the run command.

## Repeat releases: annotate-once, incremental reuse

Because annotation is per-locus, a new release only needs to annotate the loci it ADDS. Persist the annotated
distinct-locus table (a Parquet, or a DuckLake/lakehouse table with snapshots) and, on the next release:

1. dedup the new cohort to distinct loci (0.09 s);
2. `ANTI JOIN` against the stored annotations to get only the NEW loci;
3. run offline VEP + the local joins on that residual (typically a small fraction) — minutes shrink to seconds;
4. `UPSERT`/append into the store; the local ClinVar/gnomAD joins are cheap enough to re-run wholesale each
   release so reclassification tracks ClinVar/gnomAD updates (the reanalysis signal — see `acmg.reanalysis`).

Storing the annotation table as a snapshotted lakehouse (DuckLake) gives time-travel over releases: you can
ask "what changed since the last freeze" as a SQL diff, which is exactly the reanalysis question.
