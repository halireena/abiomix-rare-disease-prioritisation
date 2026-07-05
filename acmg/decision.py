"""Final shortlist + GO decision: an AUTONOMOUS assessment, with TARGETED human review.

The nuance (after Ma et al. 2025, doi:10.1101/2025.06.03.25328923): the framework is *fully automated* and
concords with human curators 99.3% of the time — so the system should COMMIT to an assessment, not punt
every case to a human. It errs toward UNDER-calling (in that study, almost only underestimation of
pathogenicity, one overcall in 300) — so the safe autonomous behaviour is: state a call, and escalate to a
human only where it genuinely matters. Each candidate gets:

  - autonomous_call : the system's own assessment (reportable / uncertain / unlikely);
  - review          : whether a HUMAN is needed, and why — reserved for the cases that need it, not all.

Buckets:
  - reportable          : P/LP in a phenotype-matched gene, inheritance-consistent -> autonomous call = report;
                          human review = confirm-before-report (clinical sign-off, not re-adjudication).
  - mixed               : signals disagree (P/LP but no phenotype match = possible incidental; strong phenotype
                          but only VUS = reclassification candidate; ClinVar conflict) -> human review REQUIRED.
  - not_enough_evidence : no P/LP, weak/absent phenotype -> autonomous call = unlikely; NO human review needed
                          unless nothing anywhere in the case is reportable (then a targeted look).

The literature arm (acmg.agent, gated) attaches to `mixed`: a recorded, human-approved model judgment
(PS3/PP4/PM3) can move a VUS toward reportable — but only after sign-off, never automatically.
"""
from __future__ import annotations
import pandas as pd

REPORTABLE, MIXED, NOT_ENOUGH = "reportable", "mixed", "not_enough_evidence"
# the system's own autonomous assessment per bucket (what it would say without a human)
_AUTONOMOUS_CALL = {REPORTABLE: "reportable", MIXED: "uncertain", NOT_ENOUGH: "unlikely"}
# does this bucket REQUIRE a human, and at what role
_REVIEW = {
    REPORTABLE: ("confirm", "clinical sign-off before reporting"),
    MIXED: ("adjudicate", "expert adjudication of conflicting evidence"),
    NOT_ENOUGH: ("none", "no review needed"),
}
PHENO_STRONG = 0.5  # normalised phenotype match considered "strong" (tunable)
_PLP = ("Pathogenic", "Likely Pathogenic")


def triage_variant(row: dict) -> tuple[str, str]:
    cls = row.get("acmg_class", "")
    pheno = float(row.get("phenotype_norm", row.get("phenotype_score", 0)) or 0)
    clinvar = str(row.get("clinvar_classification", "") or "").lower()
    plp = cls in _PLP
    # zygosity / inheritance for the carrier guard. A PURELY-recessive gene (AR, not also AD) with a single
    # MONOALLELIC het P/LP and no second allele is a CARRIER, not a reportable diagnosis (pi review #3).
    zyg = str(row.get("zygosity", "") or "")
    biallelic = zyg.startswith("hom") or bool(row.get("biallelic"))  # hom/hemi, or a compound-het pair
    ar_only = bool(row.get("gene_ar")) and not bool(row.get("gene_ad"))

    if "conflict" in clinvar:
        return MIXED, "ClinVar reports conflicting classifications"
    if plp and pheno > 0:
        if ar_only and not biallelic:
            return MIXED, f"{cls} in recessive gene {row.get('gene')} but MONOALLELIC (carrier — needs a second allele)"
        return REPORTABLE, f"{cls} in a phenotype-matched gene"
    if plp and pheno == 0:
        return MIXED, f"{cls} but no phenotype match (possible incidental / secondary finding)"
    if not plp and pheno >= PHENO_STRONG:
        return MIXED, "strong phenotype match but only VUS (candidate for reclassification)"
    return NOT_ENOUGH, "no P/LP and weak/absent phenotype support"


_PRIORITY = {REPORTABLE: 2, MIXED: 1, NOT_ENOUGH: 0}


def case_decision(shortlist: pd.DataFrame, *, case_id: str | None = None) -> dict:
    """Autonomous per-variant call + a case-level GO that says whether — and why — a human is needed."""
    rows = []
    for r in shortlist.to_dict("records"):
        bucket, reason = triage_variant(r)
        review_action, _ = _REVIEW[bucket]
        rows.append({**r, "decision": bucket, "autonomous_call": _AUTONOMOUS_CALL[bucket],
                     "review": review_action, "decision_reason": reason})
    if not rows:
        return {"case_id": case_id, "go": NOT_ENOUGH, "autonomous_assessment": "no candidate variants",
                "human_review": "none", "review_reason": "nothing to review", "variants": []}

    go = max(rows, key=lambda r: (_PRIORITY[r["decision"]], r.get("combined_score", 0)))
    bucket = go["decision"]
    review_action, review_why = _REVIEW[bucket]
    # even a "not_enough" case with no candidate gets a targeted look (nothing reportable found)
    if bucket == NOT_ENOUGH:
        review_action, review_why = "targeted", "no candidate reached reportable — a targeted human look"
    return {
        "case_id": case_id,
        "go": bucket,
        "autonomous_assessment": f"{_AUTONOMOUS_CALL[bucket]}: {go.get('gene')} {go.get('acmg_class')}",
        "human_review": review_action,   # confirm | adjudicate | targeted | none
        "review_reason": review_why,
        "lead": {"gene": go.get("gene"), "variant_key": go.get("variant_key"), "acmg_class": go.get("acmg_class")},
        "variants": rows,
    }
