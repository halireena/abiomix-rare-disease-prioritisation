"""Map Ensembl VEP output onto the ACMG kernel's `annotations` schema.

The kernel (SQL) reads one flat table with these columns per variant. `vep_to_annotations`
produces the base frame from a VEP-annotated DataFrame. The core fields (gene, consequence,
frequency, REVEL, SpliceAI) map directly from VEP; the evidence columns that come from OTHER
sources start as None here and are ENRICHED downstream by dedicated modules (acmg.constraint →
gnomad_mis_z, acmg.clinvar → PS1/PM5, acmg.nmd → nmd_escaping). None is deliberate: the kernel
ABSTAINS on a missing value rather than guessing, so a partially-enriched frame is always safe.
"""
from __future__ import annotations
import pandas as pd

# The columns the SQL kernel reads (manifests/00_spec.sql, 01_acmg_rules.sql).
REQUIRED_COLS = [
    "variant_key", "gene", "consequence", "variant_kind", "filtering_af",
    "gnomad_mis_z", "revel", "spliceai", "clinvar_same_aa", "clinvar_same_codon_lp", "nmd_escaping", "pm1", "loeuf",
    "cadd_phred", "clinvar_classification",   # CADD = opt-in supporting fallback (non-REVEL/non-splice/non-ClinVar tail)
]

# Sequence Ontology consequence -> the kernel's consequence vocabulary (same map as the production system).
SO_TO_KERNEL = {
    "missense_variant": "missense", "synonymous_variant": "synonymous",
    "frameshift_variant": "frameshift", "stop_gained": "stop_gained",
    "start_lost": "start_lost", "stop_lost": "stop_lost",
    "splice_donor_variant": "splice_donor", "splice_acceptor_variant": "splice_acceptor",
    "splice_region_variant": "splice_region",
    "inframe_deletion": "inframe_deletion", "inframe_insertion": "inframe_insertion",
    "intron_variant": "intron", "5_prime_UTR_variant": "utr", "3_prime_UTR_variant": "utr",
}
# Rough severity order (most severe first) to pick ONE consequence per variant.
_SEVERITY = [
    "frameshift", "stop_gained", "splice_donor", "splice_acceptor", "start_lost", "stop_lost",
    "inframe_deletion", "inframe_insertion", "missense", "splice_region", "synonymous", "intron", "utr",
]

# Default VEP column names -> our fields. Override per your VEP flags/plugins.
# (VEP field names depend on --tab/--json, --plugin REVEL/SpliceAI, and the frequency source you enable.)
DEFAULT_VEP_COLS = {
    "chrom": "CHROM", "pos": "POS", "ref": "REF", "alt": "ALT",
    "gene": "SYMBOL", "consequence": "Consequence",
    "revel": "REVEL", "spliceai": "SpliceAI_pred_DS_max",
    "gnomad_af": "gnomADe_AF",  # gnomAD v2 exome AF (GRCh37-native) — a good default for ACMG frequency
}


def variant_kind(ref: str, alt: str) -> str:
    return "snv" if (len(ref) == 1 and len(alt) == 1 and ref not in "-" and alt not in "-") else "indel"


def most_severe_kernel_consequence(so_field: str) -> str | None:
    """Pick the most severe kernel consequence from a VEP Consequence string (comma/&-separated SO terms)."""
    if not so_field or pd.isna(so_field):
        return None
    terms = str(so_field).replace("&", ",").split(",")
    mapped = [SO_TO_KERNEL[t.strip()] for t in terms if t.strip() in SO_TO_KERNEL]
    for c in _SEVERITY:
        if c in mapped:
            return c
    return None


def _num(v):
    try:
        return float(v) if v not in (None, "", ".", "-") and not pd.isna(v) else None
    except (TypeError, ValueError):
        return None


def vep_to_annotations(vep_df: pd.DataFrame, cols: dict | None = None) -> pd.DataFrame:
    """VEP-annotated rows -> the kernel's `annotations` table (one row per variant).

    `cols` maps our field names to your VEP column names (see DEFAULT_VEP_COLS). Fields the kernel
    needs but VEP does not provide out of the box default to None and the kernel abstains on them:
      - gnomad_mis_z: gene missense-constraint Z (gnomAD gene-constraint table; join by gene) -> PP2
      - clinvar_same_aa / clinvar_same_codon_lp: ClinVar cross-reference at the same aa / codon -> PS1 / PM5
      - nmd_escaping: 0/1 from the NMD rule (last exon / <50nt from last exon-exon junction) -> PVS1 strength
    """
    c = {**DEFAULT_VEP_COLS, **(cols or {})}
    g = lambda row, key: row[c[key]] if c[key] in row else None
    out = []
    for _, row in vep_df.iterrows():
        ref, alt = str(g(row, "ref")), str(g(row, "alt"))
        vk = f"{str(g(row, 'chrom')).replace('chr', '')}-{int(g(row, 'pos'))}-{ref}-{alt}"
        out.append({
            "variant_key": vk,
            "gene": g(row, "gene"),
            "consequence": most_severe_kernel_consequence(g(row, "consequence")),
            "variant_kind": variant_kind(ref, alt),
            "filtering_af": _num(g(row, "gnomad_af")),
            "gnomad_mis_z": None,           # enriched by acmg.constraint.add_mis_z (join by gene)
            "revel": _num(g(row, "revel")),
            "spliceai": _num(g(row, "spliceai")),
            "clinvar_same_aa": None,        # enriched by acmg.clinvar.ps1_pm5 (same-aa P/LP cross-ref)
            "clinvar_same_codon_lp": None,  # enriched by acmg.clinvar.ps1_pm5 (same-codon P/LP)
            "nmd_escaping": None,           # enriched by acmg.nmd.nmd_escaping (transcript position)
        })
    return pd.DataFrame(out, columns=REQUIRED_COLS)
