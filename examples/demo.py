"""End-to-end demo that runs WITHOUT VEP installed: a synthetic VEP-shaped table -> ACMG classes.

Replace `vep_df` with a real VEP GRCh37 output (see README) and it works unchanged.
"""
import pandas as pd
from acmg import vep_to_annotations, classify

# What a VEP GRCh37 run (with REVEL + SpliceAI plugins + gnomAD v2 exome AF) gives you, per variant:
vep_df = pd.DataFrame([
    # a frameshift LoF, rare in gnomAD
    dict(CHROM="17", POS=43093464, REF="AC", ALT="A", SYMBOL="BRCA1",
         Consequence="frameshift_variant", REVEL=None, SpliceAI_pred_DS_max=None, gnomADe_AF=5e-6),
    # a missense, high REVEL, rare
    dict(CHROM="1", POS=100000, REF="A", ALT="G", SYMBOL="SCN1A",
         Consequence="missense_variant", REVEL=0.95, SpliceAI_pred_DS_max=0.01, gnomADe_AF=1e-5),
    # a common missense -> stand-alone benign by frequency
    dict(CHROM="6", POS=26090951, REF="C", ALT="G", SYMBOL="HFE",
         Consequence="missense_variant", REVEL=0.31, SpliceAI_pred_DS_max=0.0, gnomADe_AF=0.135),
    # a missense with unknown frequency -> PM2 ABSTAINS (no data != rare)
    dict(CHROM="2", POS=200000, REF="C", ALT="T", SYMBOL="XYZ",
         Consequence="missense_variant", REVEL=0.7, SpliceAI_pred_DS_max=None, gnomADe_AF=None),
])

annotations = vep_to_annotations(vep_df)

# VEP does not tell you NMD status, so `nmd_escaping` is None and PVS1 ABSTAINS (fail-safe). A real
# pipeline fills it from the NMD rule (last exon / <50nt from the last exon-exon junction). Here we set
# the frameshift to 0 (NMD-triggering) to show PVS1 firing:
annotations.loc[annotations["consequence"] == "frameshift", "nmd_escaping"] = 0

print("annotations (kernel schema):")
print(annotations.to_string(), "\n")

# For real PVS1 on the LoF, pass real ClinGen curation. Illustrative here:
gene_curation = pd.DataFrame([
    dict(gene="BRCA1", disease_id="MONDO:0007254", mode_of_inheritance="AD", hi_score=3,
         lof_mechanism=True, gene_disease_validity="Definitive", source="ClinGen"),
])

result = classify(annotations, gene_curation=gene_curation)
print("ACMG classification:")
print(result.to_string())
