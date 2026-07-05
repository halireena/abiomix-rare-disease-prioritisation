#!/usr/bin/env bash
# Reproducible reference-data fetch for the ACMG kernel + capstone. ONE place for every download that was
# previously done ad hoc, so nothing goes unnoticed or unversioned. Idempotent: skips a file that already
# exists (delete to refresh). Records .cache/MANIFEST.txt (sha256 + size + date) = what a run actually used.
# Pairs with docs/data_versions.md (the version table). All GRCh37.
#
#   bash scripts/fetch_data.sh            # core tables (ClinVar, gnomAD, GENCODE, REVEL, HPO, ClinGen)
#   WITH_CADD_SLICE=1 bash scripts/fetch_data.sh   # + exon±500bp CADD v1.7 slice (heavy; see build_cadd_slice.sh)
#   WITH_VEP=1 bash scripts/fetch_data.sh           # + offline VEP 105 cache note (micromamba)
set -euo pipefail
cd "$(dirname "$0")/.."
CACHE="${CACHE:-.cache}"; mkdir -p "$CACHE"

# fetch URL DEST -- skip if present, curl otherwise. get "url" dest.gz "what it's for"
get() { local url="$1" dest="$CACHE/$2" note="$3"
  if [[ -s "$dest" ]]; then echo "  ok   $2  ($(du -h "$dest" | cut -f1))  [$note]"; return; fi
  echo "  GET  $2  <- $url  [$note]"; curl -fL --retry 3 -o "$dest" "$url"; }

echo "== ClinVar (NCBI FTP; rolling — record the download date) =="
# variant_summary: aggregate classifications + ReviewStatus/NumberSubmitters -> PS1/PM5, class, review stars
get "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz" variant_summary.txt.gz "PS1/PM5, class, review stars"
# clinvar.vcf.gz: INFO_MC/CLNSIG/CLNREVSTAT/CLNSIGCONF/AF_* -> CLNSIGCONF conflict-composition + BA1 exceptions + validation harness truth
get "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz" clinvar_grch37.vcf.gz "CLNSIGCONF, BA1 exceptions, harness truth"
get "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi" clinvar_grch37.vcf.gz.tbi "clinvar VCF tabix index"

echo "== gnomAD v2.1.1 (GRCh37-native) =="
# by-gene constraint: mis_z (PP2), oe_lof_upper/LOEUF (opt-in constraint-PVS1) -- NOT population AF (that's remote-range, build_site_scores.py)
get "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz" gnomad_constraint.txt.gz "PP2 mis_z + LOEUF"
# Per-site AF (PM2/BA1/BS1) is NOT bulk-downloaded: scripts/build_site_scores.py range-reads the cohort's loci
# only, over remote tabix, using this local .tbi. The 60-120GB site VCF never lands.
get "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi" gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi "remote-tabix index for build_site_scores"

echo "== GENCODE v46lift37 (basic, MANE-tagged; NMD/PVS1 transcript model) =="
get "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh37_mapping/gencode.v46lift37.basic.annotation.gtf.gz" gencode.lift37.gtf.gz "NMD / MANE transcript (PVS1)"

echo "== REVEL v1.3 (PP3/BP4 missense; single calibrated predictor, Pejaver 2022) =="
# Upstream ships a 667MB zip (no range+tabix host). We keep the derived, column-pruned parquet (chrom/pos/ref/alt/revel).
if [[ -s "$CACHE/revel_grch37.parquet" ]]; then echo "  ok   revel_grch37.parquet  ($(du -h "$CACHE/revel_grch37.parquet"|cut -f1))  [PP3/BP4 missense]"
else echo "  MISS revel_grch37.parquet — derive from upstream zip (one-off):"
     echo "       curl -fLo revel.zip https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip && unzip revel.zip"
     echo "       then load revel_with_transcript_ids CSV -> parquet(chrom,pos,ref,alt,revel) on GRCh37 (grch37_pos col)"; fi

echo "== HPO (phenotype ranking; FastHPOCR index) =="
get "https://purl.obolibrary.org/obo/hp.obo" hp.obo "HPO ontology"
if [[ -s "$CACHE/hp.index" ]]; then echo "  ok   hp.index  ($(du -h "$CACHE/hp.index"|cut -f1))  [FastHPOCR index]"
else echo "  BUILD hp.index <- python $CACHE/_build_index.py (FastHPOCR over hp.obo)"; python "$CACHE/_build_index.py" || true; fi

echo "== ClinGen (gene_curation: PVS1 gate, validity cap, MOI; rolling — record file-created date) =="
# search.clinicalgenome.org gene-validity + gene-dosage downloads, via acmg.clingen (handles CSV headers)
python -c "from acmg import clingen; import os; clingen.download(os.environ.get('CACHE','.cache'))" \
  && echo "  ok   clingen_gene_validity.csv + clingen_dosage.csv"

echo "== Monarch KG (phenotype->gene; rolling — ATTACHed remote, snapshot for a reproducible run) =="
echo "  note monarch-kg/latest attached remotely by acmg.rank (no local copy); pin the snapshot URL per run."

if [[ "${WITH_CADD_SLICE:-0}" == "1" ]]; then
  echo "== CADD v1.7 exon±500bp SLICE (optional non-missense breadth; heavy remote tabix) =="
  bash scripts/build_cadd_slice.sh
fi

if [[ "${WITH_VEP:-0}" == "1" ]]; then
  echo "== Offline VEP 105 (bioconda; cohort annotation path) =="
  echo "  micromamba create -n vep -c bioconda -c conda-forge ensembl-vep=105.0"
  echo "  vep_install -a cf -s homo_sapiens -y GRCh37 -c \$HOME/.vep --NO_HTGS  (or fetch homo_sapiens_vep_105_GRCh37)"
fi

echo "== ClassifyCNV + bedtools (CNV path; pinned separately) =="
echo "  bash scripts/setup_classifycnv.sh   (INSTALL_BEDTOOLS=1 to pin bedtools 2.31.1 via bioconda)"

echo; echo "== recording .cache/MANIFEST.txt (sha256 + size + date used) =="
{ echo "# fetched $(date -u +%FT%TZ)  (git $(git rev-parse --short HEAD 2>/dev/null || echo '?'))";
  for f in "$CACHE"/*.gz "$CACHE"/*.bgz.tbi "$CACHE"/*.parquet "$CACHE"/*.csv "$CACHE"/*.obo "$CACHE"/hp.index; do
    [[ -s "$f" ]] || continue; printf "%s  %s  %s\n" "$(sha256sum "$f"|cut -c1-16)" "$(du -h "$f"|cut -f1)" "$(basename "$f")"; done
} | tee "$CACHE/MANIFEST.txt"
echo "done. version table: docs/data_versions.md"
