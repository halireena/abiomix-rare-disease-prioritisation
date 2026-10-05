"""GO decision: per-variant triage buckets (incl. the AR-carrier guard and ClinVar-conflict escalation) and the
case-level call. Pure Python/pandas, no network. Run: `PYTHONPATH=. python3 tests/test_decision.py`."""
import pandas as pd
from acmg import decision
from acmg.decision import MIXED, NOT_ENOUGH, REPORTABLE, case_decision, triage_variant


def _bucket(**row):
    return triage_variant(row)[0]


def test_plp_with_phenotype_is_reportable():
    assert _bucket(acmg_class="Pathogenic", phenotype_norm=0.8, gene="G") == REPORTABLE
    assert _bucket(acmg_class="Likely Pathogenic", phenotype_norm=0.1, gene="G") == REPORTABLE


def test_plp_without_phenotype_is_possible_incidental():
    bucket, reason = triage_variant(dict(acmg_class="Pathogenic", phenotype_norm=0.0))
    assert bucket == MIXED and "incidental" in reason


def test_strong_phenotype_vus_is_reclassification_candidate():
    assert _bucket(acmg_class="VUS", phenotype_norm=decision.PHENO_STRONG) == MIXED
    assert _bucket(acmg_class="VUS", phenotype_norm=decision.PHENO_STRONG - 0.01) == NOT_ENOUGH


def test_clinvar_conflict_always_escalates():
    row = dict(acmg_class="Pathogenic", phenotype_norm=1.0,
               clinvar_classification="Conflicting classifications of pathogenicity")
    assert _bucket(**row) == MIXED


def test_falls_back_to_raw_phenotype_score():
    assert _bucket(acmg_class="Pathogenic", phenotype_score=3.0) == REPORTABLE
    assert _bucket(acmg_class="Pathogenic") == MIXED   # no phenotype info at all -> treated as no match


def test_ar_carrier_guard():
    base = dict(acmg_class="Pathogenic", phenotype_norm=0.9, gene="CCNO", gene_ar=True, gene_ad=False)
    bucket, reason = triage_variant({**base, "zygosity": "het"})
    assert bucket == MIXED and "MONOALLELIC" in reason          # single het in an AR-only gene = carrier
    assert _bucket(**base, zygosity="hom") == REPORTABLE         # homozygous -> diagnosis candidate
    assert _bucket(**base, zygosity="het", biallelic=True) == REPORTABLE   # compound-het pair
    assert _bucket(**{**base, "gene_ad": True}, zygosity="het") == REPORTABLE  # AR+AD gene: het can be causal


def test_case_decision_empty_shortlist():
    go = case_decision(pd.DataFrame(), case_id="C0")
    assert go["go"] == NOT_ENOUGH and go["variants"] == [] and go["case_id"] == "C0"


def test_case_decision_picks_highest_bucket_then_score():
    shortlist = pd.DataFrame([
        dict(gene="VUSGENE", variant_key="1-1-A-G", acmg_class="VUS", phenotype_norm=0.1, combined_score=0.9),
        dict(gene="LEAD", variant_key="1-2-A-G", acmg_class="Pathogenic", phenotype_norm=0.6, combined_score=0.8),
        dict(gene="OTHER", variant_key="1-3-A-G", acmg_class="Likely Pathogenic", phenotype_norm=0.2,
             combined_score=0.5),
    ])
    go = case_decision(shortlist, case_id="C1")
    assert go["go"] == REPORTABLE
    assert go["lead"] == {"gene": "LEAD", "variant_key": "1-2-A-G", "acmg_class": "Pathogenic"}
    assert go["human_review"] == "confirm"
    assert [v["decision"] for v in go["variants"]] == [NOT_ENOUGH, REPORTABLE, REPORTABLE]


def test_case_decision_nothing_reportable_gets_targeted_look():
    shortlist = pd.DataFrame([dict(gene="X", variant_key="1-1-A-G", acmg_class="VUS", phenotype_norm=0.0,
                                   combined_score=0.2)])
    go = case_decision(shortlist)
    assert go["go"] == NOT_ENOUGH
    assert go["human_review"] == "targeted"
    assert go["variants"][0]["review"] == "none"   # per-variant: no review; case-level: a targeted look


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name)
    print("all passed")
