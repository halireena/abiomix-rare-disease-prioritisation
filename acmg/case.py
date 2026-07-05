"""Phenotype-AWARE per-proband classification — the case layer ON TOP of the phenotype-free kernel.

`acmg.kernel.classify` is CASE-UNAWARE: it fires only variant-intrinsic (but gene-disease-aware) criteria —
frequency (BA1/BS1/PM2), computational (PP3/BP4), PVS1, PS1/PM5, PP2. That is the reusable, reproducible floor,
classified ONCE per variant across every proband.

This module folds the CASE-OBSERVATIONAL criteria — which ACMG/ClinGen (Tavtigian 2020 points; SVI) treats as
requiring the specific case/family/functional data — onto those points to re-derive a PER-PROBAND class:

  - INHERITANCE (deterministic, from acmg.family + the pedigree): PS2 de novo, PM3 in-trans, PP1 co-segregation,
    BS4 non-segregation. Pure SQL over the trio/multiplex structure — no judgment, recomputable on demand.
  - PHENOTYPE + LITERATURE (gated agent facts, acmg.agent): PP4 phenotype-specificity, PS3 functional, PS4
    prevalence. Each is an APPROVED (input, judgment) point that counts ONLY for its exact approved digest.

Because this fold and acmg.rank are both deterministic + recomputable, the whole phenotype-aware layer can be
re-run on demand — the substrate for an iterative REPL (refine HPO -> re-rank -> re-fold), with the LLM
querying the case variant table and proposing gated case facts, never editing the class directly.
"""
from __future__ import annotations
import pandas as pd

# Tavtigian 2020 points for the CASE-OBSERVATIONAL criteria (VCEP-tunable; same point scale as the kernel).
CASE_POINTS = {"PS2": 4, "PM6": 2, "PM3": 2, "PP1": 1, "PP1_Moderate": 2, "PP1_Strong": 4,
               "BS4": -4, "PP4": 1, "PP4_Moderate": 2, "PS3": 4, "PS3_Moderate": 2, "PS3_Supporting": 1, "PS4": 4}

_BANDS = ((10, "Pathogenic"), (6, "Likely Pathogenic"), (0, "VUS"), (-6, "Likely Benign"))

# Case structures that have >1 sequenced member -> the pedigree can support family criteria (PS2/PM3/PP1).
_FAMILY_STRUCTURES = {"duo", "trio", "quad", "multiplex"}

# Inheritance criteria that REQUIRE a supporting pedigree — deterministic from the trio/multiplex, never a
# judgment. An agent fact must not assert these on a singleton (no parents to confirm de novo, no phasing for
# in-trans); they are gated on case_structure exactly like the deterministic path.
_FAMILY_CRITERIA = {"PS2", "PM6", "PM3", "PP1", "PP1_Moderate", "PP1_Strong", "BS4"}


def band(points: float) -> str:
    """Tavtigian 2020 points -> 5-tier class (same thresholds as manifests/02_combine.sql)."""
    for thr, cls in _BANDS:
        if points >= thr:
            return cls
    return "Benign"


def _is_ar(moi: str) -> bool:
    return any(m in str(moi or "").lower() for m in ("recessive", "autosomal recessive", "ar", "xlr"))


def _is_ad(moi: str) -> bool:
    """Dominant-capable MOI (incl. semidominant). A gene curated for BOTH AD and AR (moi like 'AD,AR') is
    dominant-capable, so a lone het there is a reportable monoallelic candidate — NOT merely a carrier."""
    toks = {t.strip().lower() for t in str(moi or "").split(",")}
    return bool(toks & {"ad", "sd", "xld"}) or any("dominant" in t for t in toks)


def genotype_fit(case_variants: pd.DataFrame, zygosity, gene_moi) -> pd.DataFrame:
    """Per-variant REPORTABILITY from the PROBAND's zygosity x the gene's MOI — the signal that carries the
    SINGLETON (proband-only) recessive diagnoses, where no family criterion can apply. This is a ranking/decision
    signal, NOT ACMG points (zygosity does not change a variant's pathogenicity):
      - AR + homozygous                              -> 'biallelic_homozygous'  (directly reportable; no family);
      - AR + het WITH another het in the same gene   -> 'biallelic_candidate'   (comp-het CANDIDATE, unphased -> NOT PM3);
      - AR + single het                              -> 'monoallelic_carrier'   (not a diagnosis alone);
      - AD/XL (or unknown MOI)                       -> 'monoallelic_dominant'  (a single allele is reportable).
    `zygosity`: dict/Series variant_key -> 'hom'|'het'. `gene_moi`: dict gene -> MOI string (from gene_curation)."""
    df = case_variants.copy()
    z = (lambda k: str((zygosity.get(k, "") if hasattr(zygosity, "get") else "")).lower())
    df["_zyg"] = df["variant_key"].map(z)
    df["_ar"] = df["gene"].map(lambda g: _is_ar(gene_moi.get(g) if hasattr(gene_moi, "get") else None))
    df["_ad"] = df["gene"].map(lambda g: _is_ad(gene_moi.get(g) if hasattr(gene_moi, "get") else None))
    het_ar = df[(df["_zyg"] == "het") & df["_ar"]].groupby("gene").size()
    multi_het_ar = set(het_ar[het_ar >= 2].index)     # AR gene with >=2 het variants in this proband

    def _fit(r):
        if r["_ar"]:
            if r["_zyg"] == "hom":
                return "biallelic_homozygous"
            if r["_zyg"] == "het":
                if r["gene"] in multi_het_ar:
                    return "biallelic_candidate"
                # lone het in an AR gene that is ALSO dominant-capable -> reportable dominant, not just a carrier
                return "monoallelic_dominant" if r["_ad"] else "monoallelic_carrier"
        return "monoallelic_dominant"
    df["genotype_fit"] = df.apply(_fit, axis=1)
    return df.drop(columns=["_zyg", "_ar", "_ad"])


def classify_case(case_variants: pd.DataFrame, *, case_structure: str = "singleton", zygosity=None, gene_moi=None,
                  de_novo=None, comp_het=None, segregating=None, agent_facts=None) -> pd.DataFrame:
    """Fold CASE-LEVEL criteria onto the phenotype-FREE points -> per-proband class, CASE-STRUCTURE-AWARE.

    ~93% of this cohort is proband-only (singletons), so family criteria (PS2/PM3/PP1) rarely apply; the
    singleton recessive diagnoses ride on proband ZYGOSITY x gene MOI instead (genotype_fit).

    case_variants: DataFrame with variant_key, gene, total_points (from acmg.kernel.classify), criteria, acmg_class.
    case_structure: 'singleton' | 'duo' | 'trio' | 'quad' | 'multiplex' — family criteria fire ONLY for the
                    non-singleton structures (a singleton has no parents to confirm de novo, no trio to phase).
    zygosity / gene_moi: proband zygosity per variant_key + gene->MOI -> genotype_fit (all probands; see above).
    de_novo / comp_het / segregating: variant_key iterables from acmg.family -> deterministic PS2 / PM3 / PP1
                    (applied only when case_structure supports them).
    agent_facts: list of {variant_key, criterion, points?[, proposal, ledger]} -> GATED phenotype/literature points
                 (PP4/PS3/PS4); with a proposal+ledger, ONLY an APPROVED exact-digest proposal contributes.
    Adds: genotype_fit, case_points, case_criteria, case_total_points, case_class (BA1 stand-alone never overridden)."""
    df = genotype_fit(case_variants, zygosity, gene_moi) if (zygosity is not None and gene_moi is not None) else case_variants.copy()
    df["case_points"] = 0.0
    df["case_criteria"] = ""

    def _add(mask, crit, pts):
        df.loc[mask, "case_points"] += pts
        df.loc[mask, "case_criteria"] = df.loc[mask, "case_criteria"] + crit + ","

    # 1) DETERMINISTIC inheritance criteria (from the pedigree) — ONLY when the case structure supports them.
    #    A singleton can't confirm de novo (no parents) or PM3-in-trans (no phasing) -> those stay CANDIDATES
    #    (surfaced by genotype_fit='biallelic_candidate'), never awarded here.
    if case_structure in _FAMILY_STRUCTURES:
        for keys, crit in [(de_novo, "PS2"), (comp_het, "PM3"), (segregating, "PP1")]:
            if keys is not None:
                _add(df["variant_key"].isin(set(keys)), crit, CASE_POINTS[crit])

    # 2) GATED agent facts (phenotype / literature) — count only when approved for the exact digest.
    #    Inheritance criteria (PS2/PM3/PP1/...) are pedigree-deterministic, NOT agent judgments: a fact asserting
    #    one is dropped unless the case structure actually supports it (closes the singleton gate-bypass).
    for f in (agent_facts or []):
        if f.get("criterion") in _FAMILY_CRITERIA and case_structure not in _FAMILY_STRUCTURES:
            continue
        pts = f.get("points", CASE_POINTS.get(f.get("criterion"), 0))
        if "proposal" in f:
            from .agent import approved_points
            pts = approved_points(f["proposal"], pts, ledger=f.get("ledger"))
        if pts:
            _add(df["variant_key"] == f["variant_key"], f["criterion"], pts)

    df["case_total_points"] = df["total_points"].fillna(0) + df["case_points"]
    df["case_class"] = df["case_total_points"].map(band)
    df.loc[df["criteria"].fillna("").str.contains("BA1"), "case_class"] = "Benign"  # stand-alone benign, never overridden
    df["case_criteria"] = df["case_criteria"].str.rstrip(",")
    return df
