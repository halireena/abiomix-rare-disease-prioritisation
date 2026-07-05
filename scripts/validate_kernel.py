"""Phenotype-free VALIDATION harness for the ACMG kernel: differential-test our classifications against
ClinVar EXPERT-PANEL (3-star) variants — the automatable-criteria gold standard, no VEP needed.

Everything comes from the ClinVar GRCh37 VCF via duckhts read_bcf (INFO_MC = molecular consequence, GENEINFO,
CLNHGVS = protein change, CLNSIG = truth label, CLNREVSTAT = stars), joined to our local evidence
(ClinVar PS1/PM5, gnomAD constraint, frequency). We then run acmg.kernel.classify and compare.

HONEST caveats this reports, not hides:
  - Our kernel is CASE-UNAWARE + literature-INDEPENDENT (no PS3/PP4/PM3/PS2/PP1 functional/segregation/case
    evidence), so it structurally UNDER-calls vs full expert curation. The harness measures that gap.
  - Mild ClinVar circularity (PS1/PM5 use ClinVar); PS1 self is excluded, but the signal isn't fully
    independent. Read concordance as "reproduces expert P/LP from automatable evidence", not ground truth.

Run:  PYTHONPATH=. python3 scripts/validate_kernel.py [--limit N] [--pm1]
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

CACHE = os.path.join(os.path.dirname(__file__), "..", ".cache")
VCF = os.path.join(CACHE, "clinvar_grch37.vcf.gz")

# ClinVar molecular-consequence (SO) -> the kernel's consequence vocabulary (acmg.vep_map).
_MC = {
    "missense_variant": "missense", "synonymous_variant": "synonymous", "nonsense": "stop_gained",
    "stop_gained": "stop_gained", "frameshift_variant": "frameshift", "splice_donor_variant": "splice_donor",
    "splice_acceptor_variant": "splice_acceptor", "splice_region_variant": "splice_region",
    "initiator_codon_variant": "start_lost", "inframe_deletion": "inframe_deletion",
    "inframe_insertion": "inframe_insertion", "stop_lost": "stop_lost", "intron_variant": "intron",
    "3_prime_UTR_variant": "utr", "5_prime_UTR_variant": "utr",
}
_TRUTH_PLP = ("pathogenic", "likely pathogenic", "pathogenic/likely pathogenic")
_TRUTH_BLB = ("benign", "likely benign", "benign/likely benign")


def build_truth(con: duckdb.DuckDBPyConnection, limit: int | None) -> pd.DataFrame:
    """3-star (expert panel) germline biallelic SNVs from the ClinVar VCF -> annotations + truth label."""
    con.execute("LOAD duckhts")
    lim = f"LIMIT {limit}" if limit else ""
    df = con.execute(f"""
        WITH v AS (
          SELECT replace(CAST(CHROM AS VARCHAR),'chr','')||'-'||CAST(POS AS VARCHAR)||'-'||REF||'-'||ALT[1] AS variant_key,
                 split_part(CAST(INFO_GENEINFO AS VARCHAR), ':', 1) AS gene,
                 lower(array_to_string(INFO_CLNSIG,';')) AS clnsig,
                 array_to_string(INFO_MC, ';') AS mc,
                 array_to_string(INFO_CLNHGVS, ';') AS hgvs,
                 array_to_string(INFO_CLNREVSTAT, ';') AS rev,
                 greatest(coalesce(TRY_CAST(CAST(INFO_AF_ESP AS VARCHAR) AS DOUBLE),0),
                          coalesce(TRY_CAST(CAST(INFO_AF_EXAC AS VARCHAR) AS DOUBLE),0),
                          coalesce(TRY_CAST(CAST(INFO_AF_TGP AS VARCHAR) AS DOUBLE),0)) AS filtering_af
          FROM read_bcf('{VCF}')
          WHERE len(ALT) = 1 AND len(REF) = 1 AND len(ALT[1]) = 1               -- biallelic SNV
        )
        SELECT variant_key, gene, clnsig, mc, filtering_af,
               TRY_CAST(regexp_extract(hgvs, 'p\\.[A-Za-z]{{3}}([0-9]+)[A-Za-z]{{3}}', 1) AS INTEGER) AS protein_pos,
               regexp_extract(hgvs, 'p\\.[A-Za-z]{{3}}[0-9]+([A-Za-z]{{3}})', 1) AS alt_aa3
        FROM v WHERE rev ILIKE '%expert%panel%' AND gene <> '' {lim}
    """).df()
    # map MC -> our consequence (first recognised SO term)
    def cons(mc):
        for term in str(mc).replace("|", ";").split(";"):
            key = term.split(",")[-1].strip()
            if key in _MC:
                return _MC[key]
        return None
    df["consequence"] = df["mc"].map(cons)
    df["variant_kind"] = "snv"
    df["truth"] = df["clnsig"].map(lambda c: "P/LP" if any(t in c for t in _TRUTH_PLP)
                                   else ("B/LB" if any(t in c for t in _TRUTH_BLB) else "other"))
    from acmg.clinvar import _AA3_TO_1
    df["alt_aa1"] = df["alt_aa3"].map(_AA3_TO_1)
    return df[df["truth"].isin(["P/LP", "B/LB"]) & df["consequence"].notna()].copy()


_LOF = ("stop_gained", "frameshift", "splice_donor", "splice_acceptor")


def run(limit: int | None, pm1: bool, leakage_free: bool) -> None:
    con = duckdb.connect()
    print("loading ClinVar + constraint + exon model + ClinGen ...", flush=True)
    load_clinvar(con, f"{CACHE}/variant_summary.txt.gz")
    load_constraint(con, f"{CACHE}/gnomad_constraint.txt.gz")
    load_exons(con, f"{CACHE}/gencode.lift37.gtf.gz")                       # NMD (PVS1 strength)
    cg = clingen.download(CACHE); clingen.load_gene_curation(con, cg["gene_validity"], cg["dosage"])
    gene_curation = con.execute("SELECT * FROM gene_curation").df()          # PVS1 gate + validity cap + MOI
    truth = build_truth(con, limit)
    print(f"expert-panel (3-star) test set: {len(truth)} variants "
          f"({(truth.truth=='P/LP').sum()} P/LP, {(truth.truth=='B/LB').sum()} B/LB)", flush=True)

    ann = truth[["variant_key", "gene", "consequence", "variant_kind", "filtering_af"]].copy()
    # local REVEL (PP3/BP4 for missense) by variant_key
    rp = f"{CACHE}/revel_grch37.parquet"
    if os.path.exists(rp):
        con.register("_ak", ann[["variant_key"]])
        rev = dict(con.execute(f"SELECT r.chrom||'-'||r.pos||'-'||r.ref||'-'||r.alt, r.revel FROM read_parquet('{rp}') r JOIN _ak k ON k.variant_key = r.chrom||'-'||r.pos||'-'||r.ref||'-'||r.alt").fetchall())
        con.unregister("_ak"); ann["revel"] = ann["variant_key"].map(rev)
    # PS1/PM5 are KEPT even in leakage-free: they transfer evidence from a DIFFERENT variant (the query's own
    # ClinVar entry is self-excluded), which is the ACMG-sanctioned mechanism, not problematic circularity.
    ann["clinvar_same_aa"] = 0; ann["clinvar_same_codon_lp"] = 0
    q = truth.dropna(subset=["protein_pos", "alt_aa1"])[["variant_key", "gene", "protein_pos", "alt_aa1"]].copy()
    if len(q):
        q["protein_pos"] = q["protein_pos"].astype(int)
        cv = ps1_pm5(con, q).set_index("variant_key")
        ann["clinvar_same_aa"] = ann["variant_key"].map(cv["clinvar_same_aa"]).fillna(0).astype(int)
        ann["clinvar_same_codon_lp"] = ann["variant_key"].map(cv["clinvar_same_codon_lp"]).fillna(0).astype(int)
    ann["gnomad_mis_z"] = None  # add_mis_z fills it by gene; must exist first
    ann = add_mis_z(con, ann)
    # NMD (PVS1 strength) for LoF variants: map each gene -> its MANE transcript, then run the exon-model rule.
    lof = ann[ann["consequence"].isin(_LOF)].copy()
    ann["nmd_escaping"] = None
    if len(lof):
        tx = {g: select_transcript(con, g) for g in lof["gene"].dropna().unique()}
        lof["transcript_id"] = lof["gene"].map(tx)
        lof = lof.dropna(subset=["transcript_id"])
        lof["chrom"] = lof["variant_key"].str.split("-").str[0]
        lof["pos"] = lof["variant_key"].str.split("-").str[1].astype(int)
        if len(lof):
            nm = nmd_escaping(con, lof[["variant_key", "chrom", "pos", "transcript_id"]]).set_index("variant_key")
            ann["nmd_escaping"] = ann["variant_key"].map(nm["nmd_escaping"])
    # leakage-free also drops ClinGen gene_curation: the ClinVar EXPERT PANELS **are** ClinGen VCEPs, so
    # gene-disease validity + dosage (PVS1 gate, validity cap) are correlated with the labels. Dropping it
    # makes PVS1 abstain — leaving only ClinVar+ClinGen-INDEPENDENT evidence (frequency, REVEL, constraint).
    cls = classify(ann, con=duckdb.connect(), gene_curation=(None if leakage_free else gene_curation),
                   pm1_enabled=pm1, pvs1_constraint=True).set_index("variant_key")  # constraint PVS1 for the harness
    print(f"(mode: {'LEAKAGE-FREE (PS1/PM5 kept; no ClinGen gene_curation; PVS1 via gnomAD LOEUF constraint)' if leakage_free else 'full'})", flush=True)

    m = truth.set_index("variant_key").join(cls[["acmg_class"]])
    m["called"] = m["acmg_class"].map(lambda c: "P/LP" if c in ("Pathogenic", "Likely Pathogenic")
                                      else ("B/LB" if c in ("Benign", "Likely Benign") else "VUS"))
    print("\n=== concordance (rows = ClinVar expert truth, cols = our kernel call) ===")
    print(pd.crosstab(m["truth"], m["called"]).to_string())
    plp = m[m.truth == "P/LP"]
    print(f"\nOf expert P/LP: we call P/LP {100*(plp.called=='P/LP').mean():.1f}%, "
          f"VUS {100*(plp.called=='VUS').mean():.1f}% (the case-unaware/literature-independent GAP), "
          f"B/LB {100*(plp.called=='B/LB').mean():.1f}% (should be ~0 — any is a real bug).")
    contra = m[((m.truth == 'P/LP') & (m.called == 'B/LB')) | ((m.truth == 'B/LB') & (m.called == 'P/LP'))]
    print(f"CONTRADICTIONS (opposite call): {len(contra)}  <- these are the bugs to chase")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--limit", type=int, default=None, help="cap the test set (fast iteration)")
    ap.add_argument("--pm1", action="store_true", help="enable PM1 scoring (default off)")
    ap.add_argument("--leakage-free", action="store_true", help="zero ClinVar-derived PS1/PM5 (avoid circularity)")
    a = ap.parse_args()
    run(a.limit, a.pm1, a.leakage_free)


if __name__ == "__main__":
    main()
