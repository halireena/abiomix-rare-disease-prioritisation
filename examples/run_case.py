"""End-to-end on ONE real proband, live: ingest -> panel pre-filter -> HYBRID annotate -> classify
-> HPO -> phenotype x genotype rank -> GO decision.

Uses the real challenge bundle + the reference data in .cache/. Bounded to a gene panel so the VEP-REST
budget stays small (production pre-filters by rarity + panel, then annotates the residual novel set).

Run:  BUNDLE=/root/bioconnect/dataset_clinical_curated.parquet python examples/run_case.py CASE0003
"""
import os, sys, duckdb, pandas as pd
from acmg.clinvar import load_clinvar
from acmg.constraint import load_constraint
from acmg.nmd import load_exons
from acmg.annotate import annotate_hybrid
from acmg.kernel import classify
from acmg.hpo import extract_hpo
from acmg import rank, decision

CACHE = os.path.join(os.path.dirname(__file__), "..", ".cache")
BUNDLE = os.environ.get("BUNDLE", "/root/bioconnect/dataset_clinical_curated.parquet")
CASE = sys.argv[1] if len(sys.argv) > 1 else "CASE0003"

# a neurodevelopmental / epilepsy gene panel (the clinical pre-filter for this proband)
PANEL = ["SCN1A","SCN2A","SCN8A","STXBP1","KCNQ2","KCNT1","CDKL5","MECP2","FOXG1","SYNGAP1","GRIN2B",
         "GRIN1","GABRB3","PCDH19","DEPDC5","DNM1","GNAO1","CACNA1A","TSC1","TSC2","PTEN","ARID1B","ADNP","SLC2A1"]
# minimal ClinGen-style curation so PVS1 can fire for LoF in established haploinsufficient genes
GENE_CURATION = pd.DataFrame([
    dict(gene=g, disease_id="MONDO:panel", mode_of_inheritance="AD", hi_score=3,
         lof_mechanism=True, gene_disease_validity="Definitive", source="panel")
    for g in ["SCN1A","STXBP1","SYNGAP1","CDKL5","FOXG1","KCNQ2"]
])

con = duckdb.connect()
print(f"[1/6] reference tables ...")
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

print(f"[3/6] hybrid annotate (VEP-REST + local ClinVar PS1/PM5 + NMD + gnomAD constraint) ...")
ann = annotate_hybrid(con, vcf)
print(f"      annotated {len(ann)}; coding: {ann['consequence'].notna().sum()}")

print(f"[4/6] ACMG classify ...")
cls = classify(ann, con=duckdb.connect(), gene_curation=GENE_CURATION)
cls = cls[cls["acmg_class"] != "Not evaluated (non-SNV/indel — see Riggs 2020)"]

print(f"[5/6] HPO from clinical text + phenotype x genotype rank (Monarch) ...")
txt = con.execute(f"SELECT DISTINCT clinical_indication_text FROM read_parquet('{BUNDLE}') WHERE student_case_id='{CASE}'").fetchone()[0]
hpo = extract_hpo(txt or "", f"{CACHE}/hp.index", case_id=CASE)  # default augment_select: LLM augment + FastHPOCR ground + LLM select (spark); persisted
print(f"      observed HPO: {hpo.observed}  excluded: {hpo.excluded}  family: {hpo.family_scope}")
try:
    rank.attach_monarch(con); rank.monarch_gene_phenotype(con, hpo.observed)
    pheno = rank.phenotype_scores(con, hpo.observed)
except Exception as e:
    print(f"      (Monarch unavailable: {e}); phenotype = 0"); pheno = pd.DataFrame(columns=["gene","phenotype_score"])
ranked = rank.rerank(cls, pheno)

print(f"[6/6] GO decision\n")
go = decision.case_decision(ranked, case_id=CASE)
cols = ["gene","variant_key","acmg_class","total_points","phenotype_norm","combined_score","decision"]
print(ranked.merge(pd.DataFrame(go["variants"])[["variant_key","decision"]], on="variant_key", how="left")[cols].head(10).to_string(index=False))
print(f"\nGO: {go['go']}  | autonomous: {go['autonomous_assessment']}  | review: {go['human_review']} ({go['review_reason']})")
