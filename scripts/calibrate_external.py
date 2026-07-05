"""EXTERNAL calibration of the kernel's numerical output — no ClinVar truth labels needed.

The kernel emits a continuous score (Tavtigian `total_points`). An orthogonal, NOT-ClinVar-supervised
predictor (AlphaMissense — trained on population frequency + protein structure, not ClinVar labels; or the
CADD slice, which is NOT one of our rules) gives an independent axis. If our aggregate points track it
monotonically, the aggregation is meaningful independent of any label set. This is calibration, not accuracy.

HONEST caveat this reports: our missense points are largely REVEL-driven, and REVEL correlates with
AlphaMissense — so raw Spearman is partly built-in. We therefore ALSO report the REVEL-ABLATED correlation
(kernel re-run with REVEL zeroed): the residual points (PVS1/PS1/PM5/frequency/constraint/PM2) tracking the
external predictor is the genuinely-independent signal.

Source-pluggable:
  --source parquet:.cache/cadd_exome_slice.parquet:cadd_phred     (local; CADD is orthogonal to our rules)
  --source tabix:https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg19.tsv.gz:9   (am_pathogenicity col)

Run:  PYTHONPATH=. python3 scripts/calibrate_external.py --limit 3000 --source parquet:.cache/cadd_exome_slice.parquet:cadd_phred
"""
from __future__ import annotations
import argparse
import os
import duckdb
import pandas as pd
from acmg.clinvar import load_clinvar, ps1_pm5
from acmg.constraint import load_constraint, add_mis_z
from acmg.nmd import load_exons, select_transcript, nmd_escaping
from acmg.kernel import classify
from acmg import clingen
from scripts.validate_kernel import build_truth, _LOF

CACHE = os.path.join(os.path.dirname(__file__), "..", ".cache")


def _annotate(con, truth, drop_revel=False):
    """Build the kernel annotation frame from the harness truth set (REVEL, PS1/PM5, mis_z, NMD)."""
    ann = truth[["variant_key", "gene", "consequence", "variant_kind", "filtering_af"]].copy()
    rp = f"{CACHE}/revel_grch37.parquet"
    if os.path.exists(rp) and not drop_revel:
        con.register("_ak", ann[["variant_key"]])
        rev = dict(con.execute(f"SELECT r.chrom||'-'||r.pos||'-'||r.ref||'-'||r.alt, r.revel FROM read_parquet('{rp}') r JOIN _ak k ON k.variant_key = r.chrom||'-'||r.pos||'-'||r.ref||'-'||r.alt").fetchall())
        con.unregister("_ak"); ann["revel"] = ann["variant_key"].map(rev)
    ann["clinvar_same_aa"] = 0; ann["clinvar_same_codon_lp"] = 0
    q = truth.dropna(subset=["protein_pos", "alt_aa1"])[["variant_key", "gene", "protein_pos", "alt_aa1"]].copy()
    if len(q):
        q["protein_pos"] = q["protein_pos"].astype(int)
        cv = ps1_pm5(con, q).set_index("variant_key")
        ann["clinvar_same_aa"] = ann["variant_key"].map(cv["clinvar_same_aa"]).fillna(0).astype(int)
        ann["clinvar_same_codon_lp"] = ann["variant_key"].map(cv["clinvar_same_codon_lp"]).fillna(0).astype(int)
    ann["gnomad_mis_z"] = None
    ann = add_mis_z(con, ann)
    lof = ann[ann["consequence"].isin(_LOF)].copy()
    ann["nmd_escaping"] = None
    if len(lof):
        tx = {g: select_transcript(con, g) for g in lof["gene"].dropna().unique()}
        lof["transcript_id"] = lof["gene"].map(tx); lof = lof.dropna(subset=["transcript_id"])
        lof["chrom"] = lof["variant_key"].str.split("-").str[0]; lof["pos"] = lof["variant_key"].str.split("-").str[1].astype(int)
        if len(lof):
            nm = nmd_escaping(con, lof[["variant_key", "chrom", "pos", "transcript_id"]]).set_index("variant_key")
            ann["nmd_escaping"] = ann["variant_key"].map(nm["nmd_escaping"])
    return ann


def _external_scores(con, variant_keys, source):
    """source = 'parquet:PATH:COL' (local) or 'tabix:URL:COLIDX' (remote range-read). -> {variant_key: score}."""
    kind, ref, col = source.split(":", 2)
    con.register("_vk", pd.DataFrame({"variant_key": list(variant_keys)}))
    if kind == "parquet":
        return dict(con.execute(f"""
            SELECT p.chrom||'-'||p.pos||'-'||p.ref||'-'||p.alt AS vk, p.{col}
            FROM read_parquet('{ref}') p JOIN _vk k ON k.variant_key = p.chrom||'-'||p.pos||'-'||p.ref||'-'||p.alt
        """).fetchall())
    if kind == "tabix":
        con.execute("LOAD duckhts")
        loci = con.execute("SELECT DISTINCT split_part(variant_key,'-',1) c, CAST(split_part(variant_key,'-',2) AS BIGINT) p FROM _vk").fetchall()
        by = {}
        for c, p in loci: by.setdefault(c, []).append(p)
        out = {}
        idx = f"{CACHE}/{os.path.basename(ref)}.tbi"  # local .tbi must be fetched (small)
        for c, ps in by.items():
            lo, hi = min(ps), max(ps)
            try:
                rows = con.execute(f"SELECT column0,column1,column2,column3,column{col} FROM read_tabix('{ref}', region := '{c}:{lo}-{hi}', index_path := '{idx}')").fetchall()
                for r in rows:
                    out[f"{r[0].replace('chr','')}-{r[1]}-{r[2]}-{r[3]}"] = float(r[4]) if r[4] not in (None, "") else None
            except Exception as e:
                print(f"  (tabix {c} skipped: {e})")
        return out
    raise SystemExit(f"unknown source kind: {kind}")


def run(limit, source):
    con = duckdb.connect()
    print("loading ClinVar + constraint + exons + ClinGen ...", flush=True)
    load_clinvar(con, f"{CACHE}/variant_summary.txt.gz")
    load_constraint(con, f"{CACHE}/gnomad_constraint.txt.gz")
    load_exons(con, f"{CACHE}/gencode.lift37.gtf.gz")
    cg = clingen.download(CACHE); clingen.load_gene_curation(con, cg["gene_validity"], cg["dosage"])
    gene_curation = con.execute("SELECT * FROM gene_curation").df()
    truth = build_truth(con, limit)
    truth = truth[truth.consequence == "missense"].copy()  # calibrate on the missense concept (where REVEL/AM live)
    print(f"missense calibration set: {len(truth)} variants", flush=True)

    ext = _external_scores(con, truth.variant_key, source)
    truth["ext"] = truth.variant_key.map(ext)
    have = truth.dropna(subset=["ext"]).copy()
    print(f"external score covered: {len(have)}/{len(truth)}", flush=True)
    if len(have) < 20:
        raise SystemExit("too few external-score hits to calibrate (fetch the source / .tbi first)")

    for label, drop in [("full kernel (REVEL included — partly built-in)", False),
                        ("REVEL-ABLATED (independent residual signal)", True)]:
        ann = _annotate(con, have, drop_revel=drop)
        cls = classify(ann, con=duckdb.connect(), gene_curation=gene_curation, pvs1_constraint=True).set_index("variant_key")
        j = have.set_index("variant_key").join(cls[["total_points"]]).dropna(subset=["total_points", "ext"])
        rho = j["total_points"].corr(j["ext"], method="spearman")
        pear = j["total_points"].corr(j["ext"], method="pearson")
        print(f"\n[{label}]  n={len(j)}  Spearman rho={rho:+.3f}  Pearson={pear:+.3f}")
    print("\ncalibration reading: rho>0 means higher ACMG points track higher external pathogenicity. The "
          "REVEL-ablated rho is the leakage-free one — it uses NO predictor that is also in the kernel.")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--limit", type=int, default=3000)
    ap.add_argument("--source", required=True, help="parquet:PATH:COL or tabix:URL:COLIDX")
    a = ap.parse_args()
    run(a.limit, a.source)


if __name__ == "__main__":
    main()
