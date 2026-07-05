"""Unified SNV+CNV candidate model — ONE per-proband ranked table, so a pathogenic CNV competes with the SNVs.

The two tracks classify on DIFFERENT scales (SNV: Tavtigian points; CNV: ClinGen/Riggs score), so the common
currency is the 5-tier CLASS, not raw points. Each candidate carries its (phenotype-free or -aware) class, a
`genotype_fit`, and a phenotype match; they rank TOGETHER by (class tier, phenotype, genotype fit). This is what
`acmg.repl.CaseSession` and the per-proband top-10 operate on — SNVs and CNVs side by side.
"""
from __future__ import annotations
import pandas as pd

CLASS_TIER = {"Pathogenic": 5, "Likely Pathogenic": 4, "VUS": 3, "Likely Benign": 2, "Benign": 1}

# genotype_fit -> ranking boost (unifies the SNV zygosity×MOI fits and the CNV copy-number×dosage fits).
GENOTYPE_BOOST = {
    # SNV (acmg.case.genotype_fit)
    "biallelic_homozygous": 3, "biallelic_candidate": 2, "monoallelic_dominant": 1, "monoallelic_carrier": -1,
    # CNV (cnv_genotype_fit)
    "biallelic_loss": 3, "dominant_loss": 3, "dominant_gain": 2, "uncertain_dosage": 0,
}


def cnv_genotype_fit(cnv: pd.DataFrame, *, hi_genes=None, ts_genes=None) -> pd.Series:
    """CNV reportability from copy number × the spanned gene's dosage sensitivity — the CNV parallel to the SNV
    zygosity×MOI fit:
      - DEL to CN=0 (homozygous loss)                        -> 'biallelic_loss';
      - DEL to CN=1 of a HAPLOINSUFFICIENT gene (ClinGen HI) -> 'dominant_loss'  (reportable);
      - DUP of a TRIPLOSENSITIVE gene                        -> 'dominant_gain'  (reportable);
      - else                                                 -> 'uncertain_dosage'.
    hi_genes / ts_genes: ClinGen haploinsufficient / triplosensitive gene sets (from gene_curation dosage)."""
    hi, ts = set(hi_genes or []), set(ts_genes or [])

    def _fit(r):
        sv = str(r.get("svtype", "")).upper()
        cn = r.get("cn")
        g = r.get("gene")
        if sv == "DEL":
            if cn == 0:
                return "biallelic_loss"
            if g in hi:
                return "dominant_loss"
        if sv == "DUP" and g in ts:
            return "dominant_gain"
        return "uncertain_dosage"

    return cnv.apply(_fit, axis=1)


# A gene the clinician NAMED in the referral (targeted re-test) is a strong case-level prior — surface its
# candidates high (above same-tier non-named ones) without overriding a genuinely higher ACMG class. This is the
# ONE lever for gene-targeted indications ("test for SCN5A"), where there is no phenotype to rank on.
NAMED_GENE_BOOST = 8


def candidate_table(snv: "pd.DataFrame | None" = None, cnv: "pd.DataFrame | None" = None, *,
                    named_genes=None) -> pd.DataFrame:
    """Unify pre-classified SNV + CNV candidates into ONE ranked per-proband table.

    snv: DataFrame with variant_key, gene, `case_class` (or acmg_class), and optionally phenotype_score / genotype_fit
         (from acmg.case.classify_case).
    cnv: DataFrame with cnv_key, gene, `cnv_class`, and optionally phenotype_score / genotype_fit (cnv_genotype_fit).
    named_genes: gene symbols the clinician named in the referral -> a rank boost (see NAMED_GENE_BOOST).
    Returns: candidate_id, kind, gene, acmg_class, class_tier, phenotype_score, genotype_fit, combined_score (ranked).
    Ranking: class tier dominates (x10), then phenotype match, then the genotype-fit boost, + named-gene prior."""
    parts = []
    if snv is not None and len(snv):
        s = snv.copy()
        base = pd.DataFrame({
            "candidate_id": s["variant_key"].values,
            "kind": (s["variant_kind"] if "variant_kind" in s else pd.Series("snv", index=s.index)).values,
            "gene": s["gene"].values,
            "acmg_class": (s["case_class"] if "case_class" in s else s["acmg_class"]).values,
            "phenotype_score": (s["phenotype_score"] if "phenotype_score" in s else pd.Series(0.0, index=s.index)).values,
            "genotype_fit": (s["genotype_fit"] if "genotype_fit" in s else pd.Series("", index=s.index)).values,
        })
        for c in ("case_points", "case_criteria", "case_total_points"):   # carry the SNV fold detail through (provenance)
            if c in s:
                base[c] = s[c].values
        parts.append(base)
    if cnv is not None and len(cnv):
        c = cnv.copy()
        parts.append(pd.DataFrame({
            "candidate_id": c["cnv_key"].values, "kind": "cnv", "gene": c["gene"].values,
            "acmg_class": c["cnv_class"].values,
            "phenotype_score": (c["phenotype_score"] if "phenotype_score" in c else pd.Series(0.0, index=c.index)).values,
            "genotype_fit": (c["genotype_fit"] if "genotype_fit" in c else pd.Series("", index=c.index)).values,
        }))
    cols = ["candidate_id", "kind", "gene", "acmg_class", "class_tier", "phenotype_score", "genotype_fit", "combined_score"]
    if not parts:
        return pd.DataFrame(columns=cols)
    df = pd.concat(parts, ignore_index=True)
    df["class_tier"] = df["acmg_class"].map(CLASS_TIER).fillna(3).astype(int)
    df["phenotype_score"] = df["phenotype_score"].fillna(0.0)   # raw IC sum kept for display/provenance
    pmax = df["phenotype_score"].max()
    pnorm = df["phenotype_score"] / pmax if pmax and pmax > 0 else 0.0   # relative match strength within the case
    boost = df["genotype_fit"].map(lambda x: GENOTYPE_BOOST.get(x, 0))
    # class tier STRICTLY dominates: normalised phenotype (<=5) + genotype boost (<=3) < the 10-pt tier gap, so the
    # ACMG class is never overturned by phenotype — phenotype+genotype only re-order WITHIN a tier (the phenotype-
    # aware refinement of the two-layer model). Raw unbounded IC would let a matched VUS leapfrog a no-match LP.
    df["combined_score"] = df["class_tier"] * 10 + pnorm * 5 + boost
    if named_genes:                                             # clinician-named gene = strong case-level prior
        df["combined_score"] += df["gene"].isin(set(named_genes)).astype(int) * NAMED_GENE_BOOST
    ordered = cols + [c for c in df.columns if c not in cols]   # keep carried-through SNV-fold detail (case_points, ...)
    return df.sort_values("combined_score", ascending=False, kind="stable").reset_index(drop=True)[ordered]
