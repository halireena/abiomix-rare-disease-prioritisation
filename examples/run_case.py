"""End-to-end on ONE real proband, live: ingest -> panel pre-filter -> HYBRID annotate -> classify
-> HPO -> phenotype x genotype rank -> GO decision.

Uses the real challenge bundle + the reference data in .cache/. Bounded to a gene panel so the VEP-REST
budget stays small (production pre-filters by rarity + panel, then annotates the residual novel set).

Run:  python examples/run_case.py CASE0003 [--bundle /root/bioconnect/dataset.parquet]

Needs (see README "Where to get the data"): the Abiomix challenge bundle (not public, not in this repo),
`bash scripts/fetch_data.sh` reference data in .cache/ and network (Ensembl VEP REST, Monarch KG, the duckhts
DuckDB extension).

Privacy defaults (acmg.privacy):
  - reads the DE-IDENTIFIED `dataset.parquet` (or --bundle / $BUNDLE). The raw curated
    `dataset_clinical_curated.parquet` carries PII and is refused unless --allow-pii-input is passed.
  - HPO extraction is `--hpo-mode tool_only` (local FastHPOCR, nothing leaves the machine). The LLM modes
    (augment_select, candidates_model, model_only) send the clinical-indication text to an LLM provider via the
    `pi` CLI and need --allow-external-llm.
"""
import argparse, os, sys, duckdb, pandas as pd
from acmg.clinvar import load_clinvar
from acmg.constraint import load_constraint
from acmg.nmd import load_exons
from acmg.annotate import annotate_hybrid
from acmg.kernel import classify
from acmg.hpo import extract_hpo
from acmg import rank, decision, privacy

CACHE = os.path.join(os.path.dirname(__file__), "..", ".cache")
_ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
_ap.add_argument("case", nargs="?", default="CASE0003")
_ap.add_argument("--bundle", default=os.environ.get("BUNDLE", privacy.DEID_BUNDLE),
                 help="challenge bundle parquet (default: $BUNDLE or the de-identified dataset.parquet)")
privacy.add_privacy_args(_ap)
ARGS = _ap.parse_args()
CASE = ARGS.case
BUNDLE = privacy.check_bundle(ARGS.bundle, ARGS.allow_pii_input)
privacy.check_external_llm(ARGS.hpo_mode, ARGS.allow_external_llm, allow_pii_input=ARGS.allow_pii_input)

# a neurodevelopmental / epilepsy gene panel (the clinical pre-filter for this proband)
PANEL = ["SCN1A","SCN2A","SCN8A","STXBP1","KCNQ2","KCNT1","CDKL5","MECP2","FOXG1","SYNGAP1","GRIN2B",
         "GRIN1","GABRB3","PCDH19","DEPDC5","DNM1","GNAO1","CACNA1A","TSC1","TSC2","PTEN","ARID1B","ADNP","SLC2A1"]
# minimal ClinGen-style curation so PVS1 can fire for LoF in established haploinsufficient genes
GENE_CURATION = pd.DataFrame([
    dict(gene=g, disease_id="MONDO:panel", mode_of_inheritance="AD", hi_score=3,
         lof_mechanism=True, gene_disease_validity="Definitive", source="panel")
    for g in ["SCN1A","STXBP1","SYNGAP1","CDKL5","FOXG1","KCNQ2"]
])

# Pre-flight: say what is missing and where to get it, instead of a DuckDB "No files found" traceback.
_needed = {BUNDLE: "the Abiomix challenge bundle (pass --bundle /path/to/dataset.parquet)"}
for _f in ("variant_summary.txt.gz", "gnomad_constraint.txt.gz", "gencode.lift37.gtf.gz", "hp.index"):
    _needed[os.path.join(CACHE, _f)] = "reference data: run `bash scripts/fetch_data.sh` (needs network)"
_missing = [f"  {os.path.normpath(p)}  <- {how}" for p, how in _needed.items() if not os.path.exists(p)]
if _missing:
    sys.exit("run_case.py needs data that is not here yet:\n" + "\n".join(_missing)
             + "\nFor an offline run with no downloads, use: python examples/demo.py")

con = duckdb.connect()
print("[1/6] reference tables ...")
load_clinvar(con, f"{CACHE}/variant_summary.txt.gz", cache_parquet=f"{CACHE}/clinvar_prot.parquet")
load_constraint(con, f"{CACHE}/gnomad_constraint.txt.gz")
load_exons(con, f"{CACHE}/gencode.lift37.gtf.gz")

print(f"[2/6] {CASE}: proband carried SNVs in the panel ...")
gene_list = ",".join(f"'{g}'" for g in PANEL)
cands = con.execute(f"""
    WITH span AS (SELECT gene, chrom, min(start_pos) lo, max(end_pos) hi FROM exon
                  WHERE gene IN ({gene_list}) GROUP BY gene, chrom),
    v AS (SELECT DISTINCT replace(CAST(CHROM AS VARCHAR),'chr','') chrom, CAST(POS AS BIGINT) pos,
                 CAST(REF AS VARCHAR) ref_a, CAST(ALT[1] AS VARCHAR) alt_a
          FROM read_parquet('{BUNDLE}')
          WHERE student_case_id='{CASE}' AND variant_kind='snv'
            AND FORMAT_GT NOT IN ('0/0','0|0','./.','.|.'))
    SELECT DISTINCT v.* FROM v JOIN span s ON s.chrom=v.chrom AND v.pos BETWEEN s.lo AND s.hi
    ORDER BY chrom, pos
""").df()
print(f"      {len(cands)} candidate variants in panel genes")
vcf = [f"{r.chrom} {r.pos} . {r.ref_a} {r.alt_a} . . ." for r in cands.itertuples()]

print("[3/6] hybrid annotate (VEP-REST + local ClinVar PS1/PM5 + NMD + gnomAD constraint) ...")
ann = annotate_hybrid(con, vcf)
print(f"      annotated {len(ann)}; coding: {ann['consequence'].notna().sum()}")

print("[4/6] ACMG classify ...")
cls = classify(ann, con=duckdb.connect(), gene_curation=GENE_CURATION)
cls = cls[cls["acmg_class"] != "Not evaluated (non-SNV/indel — see Riggs 2020)"]

print("[5/6] HPO from clinical text + phenotype x genotype rank (Monarch) ...")
txt = con.execute(f"SELECT DISTINCT clinical_indication_text FROM read_parquet('{BUNDLE}') WHERE student_case_id='{CASE}'").fetchone()[0]
hpo = extract_hpo(txt or "", f"{CACHE}/hp.index", mode=ARGS.hpo_mode, case_id=CASE)  # default tool_only (local); persisted
print(f"      observed HPO: {hpo.observed}  excluded: {hpo.excluded}  family: {hpo.family_scope}")
try:
    rank.attach_monarch(con); rank.monarch_gene_phenotype(con, hpo.observed)
    pheno = rank.phenotype_scores(con, hpo.observed)
except Exception as e:
    print(f"      (Monarch unavailable: {e}); phenotype = 0"); pheno = pd.DataFrame(columns=["gene","phenotype_score"])
ranked = rank.rerank(cls, pheno)

print("[6/6] GO decision\n")
go = decision.case_decision(ranked, case_id=CASE)
cols = ["gene","variant_key","acmg_class","total_points","phenotype_norm","combined_score","decision"]
print(ranked.merge(pd.DataFrame(go["variants"])[["variant_key","decision"]], on="variant_key", how="left")[cols].head(10).to_string(index=False))
print(f"\nGO: {go['go']}  | autonomous: {go['autonomous_assessment']}  | review: {go['human_review']} ({go['review_reason']})")
