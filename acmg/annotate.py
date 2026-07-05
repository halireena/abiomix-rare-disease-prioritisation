"""Hybrid annotation glue: VEP (REST for the novel set) + LOCAL ClinVar PS1/PM5 + NMD + gnomAD constraint,
composed into the kernel's `annotations` table.

BUILD-AGNOSTIC. Nothing here is tied to GRCh37: pass a GRCh38 VEP server (`rest.ensembl.org`), a GRCh38 GTF
(`acmg.nmd.load_exons`), and `assembly='GRCh38'` to `acmg.clinvar.load_clinvar`, and the same code annotates
on GRCh38. GRCh37 is just this bundle's build.

Prerequisites (load once): acmg.clinvar.load_clinvar (-> clinvar_prot), acmg.constraint.load_constraint
(-> gene_constraint), acmg.nmd.load_exons (-> exon).
"""
from __future__ import annotations
import os
import pandas as pd
from .vep_rest import vep_annotate, GRCH37
from .clinvar import ps1_pm5, pm1_hotspot
from .constraint import add_mis_z
from .nmd import nmd_escaping

# --- Local OpenSpliceAI splice scoring (SpliceAI-equivalent, no Illumina precompute) -------------
# GRCh37 1000G reference + bundled GENCODE-V24lift37 canonical annotation ("grch37") + the released
# OSAI_MANE PyTorch models. Fills the same `spliceai` column (max delta across acceptor/donor
# gain/loss) that vep_map/kernel consume for splice PP3 (>=0.2) / BP7 (<0.2), Walker 2023.
GRCH37_FASTA = "/root/GRCh37/human_g1k_v37.fasta"
_DEFAULT_MODEL_DIR = os.environ.get(
    "OSAI_MODEL_DIR",
    os.path.join(os.path.dirname(os.path.dirname(__file__)), ".cache", "openspliceai_models", "mane_10000nt"),
)


def _spliceai_max_delta(score_strings) -> float:
    """Max delta across the four OpenSpliceAI fields (DS_AG|DS_AL|DS_DG|DS_DL) over all genes/alleles."""
    best = 0.0
    for s in score_strings or ():
        f = s.split("|")
        for x in f[2:6]:
            if x not in (".", ""):
                v = float(x)
                if v > best:
                    best = v
    return best


def spliceai_local(variants, fasta: str = GRCH37_FASTA, *, model_dir: str | None = None,
                   annotation: str = "grch37", flanking_size: int = 10000, distance: int = 50,
                   mask: int = 0, num_models: int | None = None, batch_size: int = 128,
                   threads: int | None = None, chunk_size: int = 1500) -> pd.DataFrame:
    """Score `variants` with the local OpenSpliceAI PyTorch ensemble; return a DataFrame with
    columns [chrom, pos, ref, alt, variant_key, spliceai] where `spliceai` is the max delta score
    (SpliceAI's max across acceptor/donor gain+loss). This is a drop-in local replacement for the
    VEP-REST SpliceAI column, needing no auth-gated Illumina precompute.

    `variants`: iterable of (chrom, pos, ref, alt) tuples (pos 1-based, like the VCF/parquet).
    `flanking_size`: 80|400|2000|10000 -> ±40/±200/±1000/±5000 nt context; 10000 (±5000) matches
    the SpliceAI-10k / Walker calibration. `num_models`: use only the first N of the ensemble
    (speed knob; default = all models in `model_dir`). Variants not overlapping an annotated gene
    get spliceai=0.0 (SpliceAI only scores within genes)."""
    import glob
    try:
        import torch
        import pysam
        from openspliceai.variant.utils import Annotator, get_delta_scores_batched
    except ImportError as e:
        raise ImportError(
            "on-device splicing needs the optional [splice] extra (PyTorch + OpenSpliceAI), kept out of the core "
            "install. Run `pip install -e '.[splice]'` or `bash scripts/setup_openspliceai.sh` (also fetches models)."
        ) from e

    if threads:
        torch.set_num_threads(int(threads))
    model_dir = model_dir or _DEFAULT_MODEL_DIR

    rows = [(str(c), int(p), str(r), str(a)) for (c, p, r, a) in variants]
    cols = ["chrom", "pos", "ref", "alt"]
    if not rows:
        return pd.DataFrame(columns=cols + ["variant_key", "spliceai"])

    # Optionally restrict the ensemble to the first N models (symlinked into a temp dir).
    mdl = model_dir
    if num_models is not None:
        files = sorted(glob.glob(os.path.join(model_dir, "*.pt")))[:num_models]
        mdl = os.path.join("/tmp", "osai_models_%d_%d" % (num_models, flanking_size))
        os.makedirs(mdl, exist_ok=True)
        for f in files:
            dst = os.path.join(mdl, os.path.basename(f))
            if not os.path.exists(dst):
                os.symlink(f, dst)

    ann = Annotator(fasta, annotation, mdl, "pytorch", flanking_size)

    hdr = pysam.VariantHeader()
    contigs = set()
    with open(fasta + ".fai") as fh:
        for line in fh:
            name, length = line.split("\t")[:2]
            hdr.contigs.add(name, length=int(length))
            contigs.add(name)

    # Only chroms present in the reference can be scored; others (e.g. 'M' vs 'MT',
    # unplaced contigs) get spliceai=0.0. Keep every input row in the output.
    # Score in bounded chunks: get_delta_scores_batched buffers every window for the
    # whole call, so a large `variants` list would exhaust memory — chunk_size caps it.
    delta = [0.0] * len(rows)
    scored_idx, recs = [], []
    for i, (c, p, r, a) in enumerate(rows):
        if c not in contigs:
            continue
        rec = hdr.new_record()
        rec.chrom, rec.pos, rec.ref, rec.alts = c, p, r, (a,)
        recs.append(rec)
        scored_idx.append(i)

    for s in range(0, len(recs), chunk_size):
        chunk = recs[s:s + chunk_size]
        scored = get_delta_scores_batched(chunk, ann, distance, mask, flanking_size, 2, batch_size)
        for j, sc in enumerate(scored):
            delta[scored_idx[s + j]] = _spliceai_max_delta(sc)

    out = pd.DataFrame(rows, columns=cols)
    out["variant_key"] = out["chrom"] + "-" + out["pos"].astype(str) + "-" + out["ref"] + "-" + out["alt"]
    out["spliceai"] = delta
    return out


def annotate_hybrid(con, vcf_variants: list[str], *, server: str = GRCH37,
                    clinvar: bool = True, constraint: bool = True, nmd: bool = True) -> pd.DataFrame:
    """`vcf_variants`: 'chrom pos . ref alt . . .' strings. Returns an `annotations` DataFrame ready for
    acmg.kernel.classify, with PS1/PM5, mis_z and nmd_escaping filled from LOCAL tables (VEP only supplies
    consequence/gene/REVEL/SpliceAI/AF/protein/transcript)."""
    ann = vep_annotate(vcf_variants, server=server)
    if not len(ann):
        return ann

    if clinvar:
        q = ann.dropna(subset=["gene", "protein_pos", "alt_aa1"])[["variant_key", "gene", "protein_pos", "alt_aa1"]].copy()
        if len(q):
            q["protein_pos"] = q["protein_pos"].astype(int)
            cv = ps1_pm5(con, q).set_index("variant_key")
            ann["clinvar_same_aa"] = ann["variant_key"].map(cv["clinvar_same_aa"]).fillna(0).astype(int)
            ann["clinvar_same_codon_lp"] = ann["variant_key"].map(cv["clinvar_same_codon_lp"]).fillna(0).astype(int)
            pm1 = pm1_hotspot(con, q[["variant_key", "gene", "protein_pos"]]).set_index("variant_key")  # PM1 hotspot
            ann["pm1"] = ann["variant_key"].map(pm1["pm1"]).fillna(0).astype(int)

    if constraint:
        ann = add_mis_z(con, ann)  # gene -> gnomad_mis_z (PP2)

    if nmd:
        nq = ann.dropna(subset=["transcript_id"]).copy()
        if len(nq):
            parts = nq["variant_key"].str.split("-", n=3, expand=True)
            nq["chrom"], nq["pos"] = parts[0], parts[1].astype("int64")
            res = nmd_escaping(con, nq[["variant_key", "chrom", "pos", "transcript_id"]]).set_index("variant_key")
            ann["nmd_escaping"] = ann["variant_key"].map(res["nmd_escaping"])

    return ann
