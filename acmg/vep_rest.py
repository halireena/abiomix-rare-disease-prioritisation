"""Thin Ensembl VEP REST client (GRCh37) — the "hybrid" annotation path, for NOVEL variants ONLY.

There is no well-maintained pip wrapper for VEP-over-REST (pyensembl is LOCAL Ensembl data, not REST;
`ensembl_rest` is thin/unmaintained), and the endpoint is simple, so this is a small `requests` client we
control.

REST LIMITS (Ensembl): **200 variants per POST**, and the whole REST API is rate-limited at **55,000
requests/hour**. So REST is NOT the way to annotate a whole cohort — for large-scale analysis the OFFLINE
VEP tool (cache + plugins) is the right choice. Use this client only for the HANDFUL of novel variants not
found in your local bulk tables (gnomAD / dbNSFP / ClinVar). Pre-filter first (rarity + coding/splice +
panel) so "novel" is a few hundred, not tens of thousands.

Returns the kernel's `annotations` columns plus protein_pos/alt_aa1 for the ClinVar PS1/PM5 join
(acmg.clinvar). GRCh37-native — no liftover.
"""
from __future__ import annotations
import time
import requests
from .vep_map import SO_TO_KERNEL, variant_kind, REQUIRED_COLS

GRCH37 = "https://grch37.rest.ensembl.org"
MAX_POST = 200          # Ensembl: max variants per POST /vep/human/region
_MIN_INTERVAL = 3600 / 55000  # ~0.065s between requests to stay under 55,000/hour
_SEVERITY = ["frameshift", "stop_gained", "splice_donor", "splice_acceptor", "start_lost", "stop_lost",
             "inframe_deletion", "inframe_insertion", "missense", "splice_region", "synonymous", "intron", "utr"]


def vep_region(variants: list[str], *, server: str = GRCH37, revel=True, spliceai=True, cadd=True,
               timeout: int = 60, retries: int = 3) -> list[dict]:
    """POST a batch (<=200) of VCF-style strings ('17 43093464 . C A . . .') to /vep/human/region.
    Returns the raw VEP JSON records."""
    params = {"content-type": "application/json"}
    if revel: params["REVEL"] = 1
    if spliceai: params["SpliceAI"] = 1
    if cadd: params["CADD"] = 1
    url = f"{server}/vep/human/region"
    for attempt in range(retries):
        r = requests.post(url, params=params, json={"variants": variants},
                          headers={"Content-Type": "application/json", "Accept": "application/json"}, timeout=timeout)
        if r.status_code == 429:
            time.sleep(float(r.headers.get("Retry-After", 1.0)))
            continue
        r.raise_for_status()
        return r.json()
    raise RuntimeError("VEP REST: rate-limited after retries")


def vep_annotate(variants: list[str], *, server: str = GRCH37, **kw) -> "pd.DataFrame":
    """Annotate an arbitrary number of NOVEL variants: chunk into <=200 per POST and pace under 55,000/hr.
    For a whole cohort use the OFFLINE VEP tool instead — this is for the residual novel set only."""
    import pandas as pd
    frames = []
    for i in range(0, len(variants), MAX_POST):
        frames.append(to_annotations(vep_region(variants[i:i + MAX_POST], server=server, **kw)))
        if i + MAX_POST < len(variants):
            time.sleep(_MIN_INTERVAL)
    return pd.concat(frames, ignore_index=True) if frames else to_annotations([])


def _pick(rec: dict) -> dict | None:
    """Pick the transcript consequence to report (MANE Select > canonical > first)."""
    tcs = rec.get("transcript_consequences") or []
    if not tcs:
        return None
    return (next((t for t in tcs if t.get("mane_select")), None)
            or next((t for t in tcs if t.get("canonical")), None) or tcs[0])


def _consequence(rec: dict) -> str | None:
    mapped = [SO_TO_KERNEL[t] for t in (rec.get("most_severe_consequence", "").split(",")) if t in SO_TO_KERNEL]
    for c in _SEVERITY:
        if c in mapped:
            return c
    return None


def to_annotations(records: list[dict]) -> "pd.DataFrame":
    """VEP REST JSON -> the kernel's `annotations` columns + protein_pos/alt_aa1 (for ClinVar PS1/PM5).
    `filtering_af` = the max of the gnomAD exome/genome AFs VEP returns (max-population is refined via the
    kernel's variant_frequency table); REVEL/SpliceAI from the plugin fields."""
    import pandas as pd
    rows = []
    for rec in records:
        alleles = rec.get("input", "").split()
        ref, alt = (alleles[3], alleles[4]) if len(alleles) >= 5 else ("N", "N")
        tc = _pick(rec) or {}
        af = None
        for k in ("gnomADe_AF", "gnomADg_AF", "gnomad_exomes_af", "af"):
            v = (rec.get("colocated_variants") or [{}])[0].get(k) if rec.get("colocated_variants") else None
            if v is not None:
                af = float(v); break
        aa = (tc.get("amino_acids") or "/").split("/")
        rows.append({
            "variant_key": f"{rec.get('seq_region_name')}-{rec.get('start')}-{ref}-{alt}",
            "gene": tc.get("gene_symbol"),
            "consequence": _consequence(rec),
            "variant_kind": variant_kind(ref, alt),
            "filtering_af": af,
            "gnomad_mis_z": None,   # join a gnomAD gene-constraint table by gene (PP2) — see README
            "revel": tc.get("revel"),
            "spliceai": max([v for k, v in tc.items() if k.startswith("spliceai_pred_ds") and v is not None], default=None),
            "clinvar_same_aa": None,        # filled by acmg.clinvar.ps1_pm5
            "clinvar_same_codon_lp": None,  # filled by acmg.clinvar.ps1_pm5
            "nmd_escaping": None,           # filled by acmg.nmd.nmd_escaping (needs transcript_id)
            "protein_pos": tc.get("protein_start"),
            "alt_aa1": aa[1] if len(aa) > 1 and aa[1] else None,
            "transcript_id": tc.get("transcript_id"),
        })
    return pd.DataFrame(rows)
