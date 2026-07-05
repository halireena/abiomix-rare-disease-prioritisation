"""Per-proband command ON THE DATALAKES: load one proband's SNV + CNV candidates, join the phenotype-free
classes, derive the case metadata (zygosity, MOI, HI/TS gene sets, structure, HPO), and run the phenotype-aware,
CNV-inclusive fold -> the top-10 (the challenge deliverable) + an open REPL session.

Composition of shipped pieces (acmg.proband.run_proband). Needs: the SNV datalake
(.cache/cohort_annotations_classified.parquet from scripts/annotate_cohort.py --classify), the CNV datalake
(.cache/cnv_store.parquet from scripts/annotate_cnv.py), and the prepared bundle (prepared/{snv,cnv,pedigree_
corrected}.parquet). ClinGen gene_curation is fetched for MOI + HI/TS.

Run:  PYTHONPATH=. python3 scripts/run_proband.py CASE0007 [--out CASE0007_top10.csv]
"""
from __future__ import annotations
import argparse
import os
import re
import duckdb
import pandas as pd
from acmg.proband import run_proband
from acmg.case import _FAMILY_STRUCTURES

PREP = "/root/bioconnect/prepared"


def _cnv_gene(dosage_genes, svtype, hi_genes, ts_genes):
    """Pick a CNV's representative dosage-sensitive gene from ClassifyCNV's own gene list: the first spanned gene
    that is HI (for a DEL) or TS (for a DUP), else the first listed gene. Keys genotype_fit + phenotype rank."""
    if not dosage_genes or str(dosage_genes).lower() in ("", "nan", "none"):
        return None
    genes = [g.strip() for g in re.split(r"[,;|\s]+", str(dosage_genes)) if g.strip()]
    want = hi_genes if str(svtype).upper() == "DEL" else ts_genes
    for g in genes:
        if g in want:
            return g
    return genes[0] if genes else None


def _zygosity(gt) -> str:
    a = str(gt or "").replace("|", "/").split("/")
    return "hom" if len(a) == 2 and a[0] == a[1] and a[0] not in ("0", ".", "") else "het"


def _structure(ped: pd.DataFrame, case_id: str) -> str:
    """Case structure from the pedigree. A family structure is returned ONLY when there is genuine family
    evidence — a SEQUENCED parent of the proband, or >=2 affected members (multiplex). A 2-member 'family' with
    no parent and not multi-affected (e.g. a primary + repeat of the same patient) is singleton-like: family
    criteria can't apply, so we don't open the family gate for it."""
    if ped.empty or case_id not in set(ped.get("individual_id", [])):
        return "singleton"
    row = ped.loc[ped["individual_id"] == case_id].iloc[0]
    members = ped.loc[ped["family_id"] == row["family_id"]]
    ids = set(members["individual_id"])
    # HARD evidence only. is_multiplex is NOT trusted on its own: some families are flagged multiplex yet have <2
    # affected coded and no sequenced parent (e.g. FAM0001: 4 'child' members, status 'unclear', no parents) — no
    # family criterion can fire there, so opening the gate is wrong. A sequenced parent (father/mother_id is itself
    # a member; PLINK '0' is not) enables de novo/in-trans; >=2 affected (PLINK affected==2) enables segregation.
    father_present = row.get("father_id") in ids
    mother_present = row.get("mother_id") in ids
    n_affected = int((members["affected"] == 2).sum()) if "affected" in members else 0
    if father_present or mother_present:
        if len(members) >= 4:
            return "quad"
        return "trio" if (father_present and mother_present) else "duo"
    if n_affected >= 2:
        return "multiplex"                      # affected-sib segregation without sequenced parents
    return "singleton"                          # >=2 members but no parent & not >=2 affected


def _clinical_text(con, case_id: str) -> str:
    src = "/root/bioconnect/dataset_clinical_curated.deid.parquet"
    src = src if os.path.exists(src) else "/root/bioconnect/dataset.parquet"
    row = con.execute(f"SELECT DISTINCT clinical_indication_text FROM read_parquet('{src}') WHERE student_case_id = ?",
                      [case_id]).fetchone()
    return row[0] if row and row[0] else ""


def _hpo(con, case_id: str, cache: str) -> list:
    idx = os.path.join(cache, "hp.index")
    text = _clinical_text(con, case_id)
    if not text or not os.path.exists(idx):
        return []
    from acmg.hpo import extract_hpo
    return extract_hpo(text, idx, case_id=case_id).observed   # LLM-aided, persisted


def _referral_genes(text: str, known_genes) -> set:
    """Gene symbols the clinician NAMED in the referral (a targeted re-test, e.g. 'test for SCN5A'). Match
    HGNC-symbol-shaped tokens against the known-gene universe so we never pick up ordinary words. Legitimate
    clinical CONTEXT (the referral), not the answer key."""
    if not text or not known_genes:
        return set()
    toks = set(re.findall(r"\b[A-Z][A-Z0-9]{1,9}\b", str(text)))   # ALL-CAPS symbol-shaped tokens
    return {t for t in toks if t in known_genes}


def _pedigree_facts(con, cid: str, prepared: str, snv_class: str, structure: str, gene_moi: dict) -> dict:
    """Derive the deterministic inheritance facts (PS2 de novo / PM3 in-trans / PP1 segregation) for a FAMILY
    case, from acmg.family over the family's genotypes. Empty for singletons (no parents). gene comes from the
    VEP annotation (the prepared bundle's gene is empty), which acmg.family.compound_het needs per call."""
    if structure not in _FAMILY_STRUCTURES:
        return {}
    from acmg import family
    fam = con.execute(f"SELECT family_id FROM read_parquet('{prepared}/snv.parquet') WHERE case_id = ? LIMIT 1", [cid]).fetchone()
    if not fam:
        return {}
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE sample_call AS
        SELECT s.case_id, s.family_id, s.member_id, s.proband_status, s.relationship, s.family_size,
               s.variant_key, c.gene, s.chrom, s.pos, s.gt, s.gq, s.dp
        FROM read_parquet('{prepared}/snv.parquet') s
        LEFT JOIN read_parquet('{snv_class}') c USING (variant_key)
        WHERE s.family_id = '{fam[0]}'
    """)
    family.build_pedigree(con)
    family.build_carriage(con)
    dn = family.de_novo_candidates(con).df()
    ch = family.compound_het_candidates(con).df()
    seg = family.affected_shared_candidates(con).df()
    # GUARD (load-bearing): family criteria must not touch benign-frequency variants — a VUS in-trans with a
    # benign partner is spurious PM3, a common homozygote segregates in every affected but isn't causal. Restrict
    # all three to the non-benign candidate set (VUS or better); BA1/BS1 rarity already lands those in LB/Benign.
    cand = set(con.execute(f"""SELECT variant_key FROM read_parquet('{snv_class}')
        WHERE acmg_class IN ('Pathogenic','Likely Pathogenic','VUS')""").df()["variant_key"])
    from acmg.case import _is_ar
    dnk = [k for k in (dn.loc[dn["case_id"] == cid, "variant_key"] if len(dn) else []) if k in cand]
    # PM3 is a RECESSIVE (biallelic) criterion: keep an in-trans pair only when the GENE is AR and BOTH alleles are
    # real candidates (a VUS in-trans with a benign partner, or two hets in an AD/uncurated gene, is not PM3).
    chk = sorted({k for _, gene, a, b in ch[["gene", "variant_a", "variant_b"]].itertuples()
                  for k in ((a, b) if a in cand and b in cand and _is_ar(gene_moi.get(gene)) else ())}) if len(ch) else []
    segk = [k for k in (seg["variant_key"] if len(seg) else []) if k in cand]
    return {"de_novo": dnk, "comp_het": chk, "segregating": segk}


def load_gene_sets(con, cache: str):
    """SHARED setup, done ONCE (reused across probands by the batch runner): ClinGen gene_curation ->
    deterministic per-gene MOI aggregate + HI/TS gene sets from the FULL dosage feed (not the gene-validity-scoped
    table, which drops dosage-only genes). Returns (gene_moi, hi_genes, ts_genes)."""
    from acmg import clingen
    from acmg.clingen import _read
    cg = clingen.download(cache)
    clingen.load_gene_curation(con, cg["gene_validity"], cg["dosage"])
    gene_moi = dict(con.execute("""SELECT gene, string_agg(DISTINCT mode_of_inheritance, ',' ORDER BY mode_of_inheritance)
        FROM gene_curation WHERE mode_of_inheritance IS NOT NULL GROUP BY gene""").fetchall())
    dos = con.execute(f"""SELECT "GENE SYMBOL" AS gene, lower(trim("HAPLOINSUFFICIENCY")) AS hi,
                                 lower(trim("TRIPLOSENSITIVITY")) AS ts FROM ({_read(cg['dosage'])})""").df()
    hi_genes = set(dos.loc[dos["hi"].str.startswith("sufficient evidence", na=False), "gene"])
    ts_genes = set(dos.loc[dos["ts"].str.startswith("sufficient evidence", na=False), "gene"])
    return gene_moi, hi_genes, ts_genes


def load_known_genes(con, snv_class: str) -> set:
    """The gene-symbol universe (distinct genes in the annotation datalake) — the dictionary _referral_genes
    matches referral tokens against, so only real genes are treated as clinician-named. Loaded ONCE."""
    return set(con.execute(f"SELECT DISTINCT gene FROM read_parquet('{snv_class}') WHERE gene IS NOT NULL").df()["gene"])


def prioritize_case(con, cid: str, *, prepared: str, snv_class: str, cnv_store: str, cache: str, ped,
                    gene_moi, hi_genes, ts_genes, known_genes=None):
    """Load one proband's SNV+CNV candidates, derive case metadata + pedigree facts + HPO, and run the fold ->
    (top_df, session, meta). Shared gene sets are passed in (built once by load_gene_sets). One con is reusable
    across cases (the Monarch ATTACH is idempotent; gene_phenotype is rebuilt per case)."""
    structure = _structure(ped, cid)
    named_genes = _referral_genes(_clinical_text(con, cid), known_genes)   # clinician-named genes (targeted re-test)

    # SNV candidates: gene from the annotation (c.gene = VEP SYMBOL; the bundle's gene is empty). No class pre-filter
    # (would hide the distribution + drop a case/agent-rescuable variant); only BA1 stand-alone is provably-never a
    # candidate (classify_case pins it Benign regardless of added points), so it is the sole safe pre-drop — EXCEPT
    # in a clinician-NAMED gene, whose variants we always surface (a targeted re-test overrides frequency; e.g. the
    # GALT deletion falsely BA1'd by tandem-repeat AF miscalling).
    named_sql = (" OR c.gene IN (" + ",".join("'" + g.replace("'", "") + "'" for g in named_genes) + ")") if named_genes else ""
    snv = con.execute(f"""
        SELECT s.variant_key, c.gene, s.gt, c.total_points, c.criteria, c.acmg_class
        FROM read_parquet('{prepared}/snv.parquet') s
        JOIN read_parquet('{snv_class}') c USING (variant_key)
        WHERE s.case_id = '{cid}' AND (NOT coalesce(c.criteria, '') LIKE '%BA1%'{named_sql})
    """).df()
    zyg = {r.variant_key: _zygosity(r.gt) for r in snv.itertuples()}

    # CNV candidates -> Riggs class + ClassifyCNV's own dosage genes; representative gene = first spanned HI(DEL)/TS(DUP)
    cnv = None
    if os.path.exists(cnv_store):
        has_genes = "dosage_genes" in con.execute(f"SELECT * FROM read_parquet('{cnv_store}') LIMIT 0").df().columns
        gsel = "s.dosage_genes" if has_genes else "CAST(NULL AS VARCHAR) AS dosage_genes"
        cnv = con.execute(f"""
            WITH p AS (SELECT replace(CAST(chrom AS VARCHAR),'chr','')||'-'||CAST(start AS VARCHAR)||'-'||CAST("end" AS VARCHAR)||'-'||
                              (CASE WHEN cn<2 THEN 'DEL' WHEN cn>2 THEN 'DUP' END) AS cnv_key,
                              (CASE WHEN cn<2 THEN 'DEL' WHEN cn>2 THEN 'DUP' END) AS svtype, cn
                       FROM read_parquet('{prepared}/cnv.parquet') WHERE case_id='{cid}' AND cn<>2)
            SELECT DISTINCT p.cnv_key, p.svtype, p.cn, s.cnv_class, {gsel}
            FROM p LEFT JOIN read_parquet('{cnv_store}') s USING (cnv_key)
        """).df()
        cnv = cnv[cnv["cnv_class"].notna()].copy() if len(cnv) else None
        if cnv is not None and len(cnv):
            cnv["gene"] = [_cnv_gene(r.dosage_genes, r.svtype, hi_genes, ts_genes) for r in cnv.itertuples()]
        else:
            cnv = None

    hpo = _hpo(con, cid, cache)
    if hpo:                                    # populate gene_phenotype from Monarch for THIS case's HPO
        from acmg import rank
        rank.attach_monarch(con)               # idempotent — safe to call per case
        rank.monarch_gene_phenotype(con, hpo)
    ped_facts = _pedigree_facts(con, cid, prepared, snv_class, structure, gene_moi)

    top, session = run_proband(con, cid, snv, cnv=cnv, hpo=hpo, zygosity=zyg, gene_moi=gene_moi,
                               hi_genes=hi_genes, ts_genes=ts_genes, case_structure=structure,
                               pedigree_facts=ped_facts, named_genes=named_genes, store=None)
    meta = {"structure": structure, "n_snv": len(snv), "n_cnv": 0 if cnv is None else len(cnv),
            "n_hpo": len(hpo), "named_genes": sorted(named_genes), "family_facts": {k: len(v) for k, v in ped_facts.items()}}
    return top, session, meta


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_id")
    ap.add_argument("--prepared", default=PREP)
    ap.add_argument("--cache", default=".cache")
    ap.add_argument("--snv-class", default=".cache/cohort_annotations_classified.parquet")
    ap.add_argument("--cnv-store", default=".cache/cnv_store.parquet")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    con = duckdb.connect()
    pedp = f"{a.prepared}/pedigree_corrected.parquet"
    ped = con.execute(f"SELECT * FROM read_parquet('{pedp}')").df() if os.path.exists(pedp) else pd.DataFrame()
    gene_moi, hi_genes, ts_genes = load_gene_sets(con, a.cache)
    known = load_known_genes(con, a.snv_class)
    top, _, meta = prioritize_case(con, a.case_id, prepared=a.prepared, snv_class=a.snv_class, cnv_store=a.cnv_store,
                                   cache=a.cache, ped=ped, gene_moi=gene_moi, hi_genes=hi_genes, ts_genes=ts_genes,
                                   known_genes=known)
    out = a.out or f"{a.case_id}_top10.csv"
    top.head(10).to_csv(out, index=False)
    print(f"{a.case_id}: structure={meta['structure']}, {meta['n_snv']} SNV + {meta['n_cnv']} CNV candidates, "
          f"{meta['n_hpo']} HPO, family_facts={meta['family_facts']} -> {out}")
    print(top.head(10)[["candidate_id", "kind", "gene", "acmg_class", "combined_score"]].to_string(index=False))


if __name__ == "__main__":
    main()
