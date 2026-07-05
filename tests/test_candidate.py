"""Unified SNV+CNV candidate model (acmg.candidate): both tracks in one ranked table, common 5-tier currency."""
import pandas as pd
from acmg import candidate


def test_pathogenic_cnv_outranks_vus_snv():
    snv = pd.DataFrame([{"variant_key": "1-1-A-G", "gene": "AAA", "case_class": "VUS", "variant_kind": "snv"}])
    cnv = pd.DataFrame([{"cnv_key": "X-1-2-DEL", "gene": "DMD", "cnv_class": "Pathogenic"}])
    t = candidate.candidate_table(snv=snv, cnv=cnv)
    assert list(t["kind"]) == ["cnv", "snv"]                 # Pathogenic CNV first, VUS SNV second
    assert t.iloc[0]["candidate_id"] == "X-1-2-DEL" and t.iloc[0]["class_tier"] == 5


def test_cnv_genotype_fit():
    cnv = pd.DataFrame([
        {"cnv_key": "1-a", "gene": "G1", "svtype": "DEL", "cn": 0},   # homozygous loss
        {"cnv_key": "2-a", "gene": "HIG", "svtype": "DEL", "cn": 1},  # het del of a haploinsufficient gene
        {"cnv_key": "3-a", "gene": "TSG", "svtype": "DUP", "cn": 3},  # dup of a triplosensitive gene
        {"cnv_key": "4-a", "gene": "X", "svtype": "DEL", "cn": 1},    # het del, not HI -> uncertain
    ])
    fit = candidate.cnv_genotype_fit(cnv, hi_genes={"HIG"}, ts_genes={"TSG"})
    assert list(fit) == ["biallelic_loss", "dominant_loss", "dominant_gain", "uncertain_dosage"]


def test_phenotype_and_genotype_break_ties_within_a_tier():
    # two VUS SNVs (same tier); the phenotype-matched + biallelic one ranks first
    snv = pd.DataFrame([
        {"variant_key": "1-1-A-G", "gene": "M", "case_class": "VUS", "phenotype_score": 2.0, "genotype_fit": "biallelic_homozygous"},
        {"variant_key": "2-1-C-T", "gene": "N", "case_class": "VUS", "phenotype_score": 0.0, "genotype_fit": "monoallelic_carrier"},
    ])
    t = candidate.candidate_table(snv=snv)
    assert t.iloc[0]["candidate_id"] == "1-1-A-G"


def test_empty_and_single_track():
    assert len(candidate.candidate_table()) == 0
    only_cnv = candidate.candidate_table(cnv=pd.DataFrame([{"cnv_key": "c", "gene": "G", "cnv_class": "Likely Pathogenic"}]))
    assert len(only_cnv) == 1 and only_cnv.iloc[0]["kind"] == "cnv"


def test_named_gene_boost_surfaces_targeted_gene():
    """A clinician-named gene (targeted re-test) lifts its candidate above same-tier non-named ones, without
    overriding a genuinely higher ACMG class."""
    snv = pd.DataFrame([
        {"variant_key": "1-1-A-G", "gene": "SCN5A", "case_class": "VUS", "phenotype_score": 0.0, "genotype_fit": "monoallelic_dominant"},
        {"variant_key": "2-1-C-T", "gene": "OTHER", "case_class": "VUS", "phenotype_score": 0.0, "genotype_fit": "monoallelic_dominant"},
        {"variant_key": "3-1-C-T", "gene": "BIGGER", "case_class": "Likely Pathogenic", "phenotype_score": 0.0, "genotype_fit": ""},
    ])
    t = candidate.candidate_table(snv=snv, named_genes={"SCN5A"})
    ids = list(t["candidate_id"])
    assert ids[0] == "3-1-C-T"          # LP still tops (named-gene boost < tier gap)
    assert ids[1] == "1-1-A-G"          # named-gene VUS surfaces above the non-named VUS
    # without the boost, the two VUS tie (stable order) -> SCN5A not preferred
    t0 = candidate.candidate_table(snv=snv)
    assert t0.set_index("candidate_id").loc["1-1-A-G", "combined_score"] == t0.set_index("candidate_id").loc["2-1-C-T", "combined_score"]
