"""ClinVar evidence for PS1 / PM5 — real, from ClinVar's GRCh37 variant_summary (no liftover, no stub).

The paper (Ma et al. 2025, doi:10.1101/2025.06.03.25328923, Methods) applies these from ClinVar:
  - PS1: an established pathogenic variant causing the SAME amino-acid change.
  - PM5: a DIFFERENT pathogenic missense affecting the SAME codon (protein position).

ClinVar's `variant_summary.txt.gz` carries GeneSymbol + ClinicalSignificance + GRCh37 coords + the protein
change in the `Name` field (e.g. `...(FOXN1):c.880G>A (p.Val294Ile)`). We parse gene + protein position +
ref/alt amino acid, keep the P/LP set, and DuckDB does the join. `clinvar_classification` is an exact
chrom-pos-ref-alt lookup. GRCh37-native — same build as the bundle, no liftover.

Download once: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
"""
from __future__ import annotations
from pathlib import Path
import duckdb

_AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
    "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
    "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*", "Sec": "U",
}

# ClinVar ReviewStatus -> gold stars (0-4). ClinVar evidence is only as good as its review level: a 0-star
# "Pathogenic" (no assertion criteria) must NOT grant PS1/PM5 or rescue a common variant from BA1. rs334 is
# 4-star so it is trusted; low-review assertions are not.
_STARS = ("CASE WHEN \"ReviewStatus\" ILIKE '%practice guideline%' THEN 4 "
          "WHEN \"ReviewStatus\" ILIKE '%expert panel%' THEN 3 "
          "WHEN \"ReviewStatus\" ILIKE '%multiple submitters%' AND \"ReviewStatus\" NOT ILIKE '%conflict%' THEN 2 "
          "WHEN \"ReviewStatus\" ILIKE '%criteria provided%' THEN 1 ELSE 0 END")


def load_clinvar(con: duckdb.DuckDBPyConnection, variant_summary_gz: str, *, assembly: str = "GRCh37",
                 cache_parquet: str | None = None, min_review_stars: int = 1, clinvar_vcf: str | None = None) -> None:
    """Build two tables from ClinVar variant_summary (build-agnostic: `variant_summary` carries BOTH
    assemblies, so pass assembly='GRCh38' to run on GRCh38 with no other change):
      - clinvar (chrom,pos,ref,alt,clnsig) for exact classification lookups;
      - clinvar_prot (gene,protein_pos,alt_aa1,clnsig,is_plp) for PS1/PM5, P/LP missense only.
    Set cache_parquet to persist clinvar_prot so you parse the 440MB file only once."""
    if cache_parquet and Path(cache_parquet).exists():
        con.execute(f"CREATE OR REPLACE TABLE clinvar_prot AS SELECT * FROM read_parquet('{cache_parquet}')")
    else:
        con.execute("CREATE OR REPLACE TABLE _aa(aa3 VARCHAR, aa1 VARCHAR)")
        con.executemany("INSERT INTO _aa VALUES (?, ?)", list(_AA3_TO_1.items()))
        # P/LP eligibility for PS1/PM5. Default (no VCF): strict — exclude ALL 'conflicting'. With a ClinVar VCF,
        # use CLNSIGCONF (parsed via duckhts read_bcf) to be conflict-COMPOSITION aware: a P/LP-vs-VUS conflict
        # is still a credible pathogenic reference (no benign contradiction) and is RECOVERED; only a
        # P/LP-vs-Benign conflict stays excluded. ~25k valid refs recovered, ~1k benign-contradicted kept out.
        if clinvar_vcf:
            con.execute("LOAD duckhts")
            con.execute(f"""
                CREATE OR REPLACE TEMP TABLE _plp_ok AS
                WITH v AS (SELECT replace(CAST(CHROM AS VARCHAR),'chr','')||'-'||CAST(POS AS VARCHAR)||'-'||REF||'-'||alt_x AS variant_key,
                                  array_to_string(INFO_CLNSIG,';') sig, array_to_string(INFO_CLNSIGCONF,';') conf
                           FROM read_bcf('{clinvar_vcf}'), UNNEST(ALT) AS t(alt_x))
                SELECT DISTINCT variant_key FROM v
                WHERE (sig ILIKE '%pathogenic%' AND sig NOT ILIKE '%conflicting%' AND sig NOT ILIKE '%benign%')
                   OR (sig ILIKE '%Conflicting%' AND conf ILIKE '%pathogenic%' AND conf NOT ILIKE '%benign%')
            """)
            _pre = "clnsig ILIKE '%pathogenic%'"                                   # loose: P/LP + 'conflicting...pathogenicity'
            _elig = "cv.variant_key IN (SELECT variant_key FROM _plp_ok)"          # VCF conflict-composition gate
        else:
            _pre = "clnsig ILIKE '%pathogenic%' AND clnsig NOT ILIKE '%conflict%' AND clnsig NOT ILIKE '%benign%'"
            _elig = "TRUE"
        con.execute(f"""
            CREATE OR REPLACE TABLE clinvar_prot AS
            WITH cv AS (
              SELECT "GeneSymbol" AS gene, "ClinicalSignificance" AS clnsig,
                     replace(CAST("Chromosome" AS VARCHAR),'chr','')||'-'||TRY_CAST("PositionVCF" AS BIGINT)||'-'||
                       "ReferenceAlleleVCF"||'-'||"AlternateAlleleVCF" AS variant_key,
                     {_STARS} AS review_stars,
                     regexp_extract("Name", 'p\\.([A-Za-z]{{3}})([0-9]+)([A-Za-z]{{3}})', 1) AS ref_aa3,
                     TRY_CAST(regexp_extract("Name", 'p\\.([A-Za-z]{{3}})([0-9]+)([A-Za-z]{{3}})', 2) AS INTEGER) AS protein_pos,
                     regexp_extract("Name", 'p\\.([A-Za-z]{{3}})([0-9]+)([A-Za-z]{{3}})', 3) AS alt_aa3
              FROM read_csv('{variant_summary_gz}', delim='\t', header=true, quote='', ignore_errors=true,
                            types={{'GeneSymbol':'VARCHAR','ClinicalSignificance':'VARCHAR','Name':'VARCHAR','Assembly':'VARCHAR','ReviewStatus':'VARCHAR',
                                    'Chromosome':'VARCHAR','PositionVCF':'VARCHAR','ReferenceAlleleVCF':'VARCHAR','AlternateAlleleVCF':'VARCHAR'}})
              WHERE "Assembly" = '{assembly}' AND {_pre}
            )
            -- keep ref_aa1 + source variant_key so PS1/PM5 require a MISSENSE source + EXCLUDE self, and only
            -- trust ClinVar entries with assertion criteria (>=1 star): a 0-star "Pathogenic" must not grant PS1/PM5.
            SELECT cv.gene, cv.protein_pos, r.aa1 AS ref_aa1, a.aa1 AS alt_aa1, cv.variant_key, cv.clnsig,
                   cv.review_stars, TRUE AS is_plp
            FROM cv JOIN _aa a ON a.aa3 = cv.alt_aa3 JOIN _aa r ON r.aa3 = cv.ref_aa3
            WHERE cv.protein_pos IS NOT NULL AND cv.gene <> '' AND cv.review_stars >= {min_review_stars} AND {_elig}
        """)
        if cache_parquet:
            con.execute(f"COPY clinvar_prot TO '{cache_parquet}' (FORMAT parquet)")
        # PM1 substrate: per (gene, protein_pos) counts of P/LP and B/LB missense (>=1 star). A mutational
        # hotspot / critical region = a window with pathogenic clustering AND no benign variation — computed
        # windowed at classification time by pm1_hotspot(). Pure DuckDB compute over ClinVar we already parse.
        con.execute(f"""
            CREATE OR REPLACE TABLE clinvar_aa AS
            SELECT gene, protein_pos,
                   count(*) FILTER (WHERE is_plp) AS n_plp,
                   count(*) FILTER (WHERE is_blb) AS n_blb
            FROM (
              SELECT "GeneSymbol" AS gene,
                     TRY_CAST(regexp_extract("Name",'p\\.[A-Za-z]{{3}}([0-9]+)[A-Za-z]{{3}}',1) AS INTEGER) AS protein_pos,
                     ({_STARS}) AS stars,
                     ("ClinicalSignificance" ILIKE '%pathogenic%' AND "ClinicalSignificance" NOT ILIKE '%conflict%'
                        AND "ClinicalSignificance" NOT ILIKE '%benign%') AS is_plp,
                     ("ClinicalSignificance" ILIKE '%benign%' AND "ClinicalSignificance" NOT ILIKE '%conflict%') AS is_blb
              FROM read_csv('{variant_summary_gz}', delim='\t', header=true, quote='', ignore_errors=true,
                            types={{'GeneSymbol':'VARCHAR','ClinicalSignificance':'VARCHAR','Name':'VARCHAR','Assembly':'VARCHAR','ReviewStatus':'VARCHAR'}})
              WHERE "Assembly" = '{assembly}'
            ) WHERE protein_pos IS NOT NULL AND gene <> '' AND stars >= {min_review_stars} AND (is_plp OR is_blb)
            GROUP BY gene, protein_pos
        """)
    # exact-lookup table
    con.execute(f"""
        CREATE OR REPLACE TABLE clinvar AS
        SELECT "Chromosome" AS chrom, TRY_CAST("PositionVCF" AS BIGINT) AS pos,
               "ReferenceAlleleVCF" AS ref, "AlternateAlleleVCF" AS alt, "ClinicalSignificance" AS clnsig,
               {_STARS} AS review_stars
        FROM read_csv('{variant_summary_gz}', delim='\t', header=true, quote='', ignore_errors=true,
                      types={{'Chromosome':'VARCHAR','PositionVCF':'VARCHAR','ReferenceAlleleVCF':'VARCHAR','AlternateAlleleVCF':'VARCHAR','ClinicalSignificance':'VARCHAR','ReviewStatus':'VARCHAR','Assembly':'VARCHAR'}})
        WHERE "Assembly" = '{assembly}'
    """) if not (cache_parquet and Path(cache_parquet).exists()) else None


def ps1_pm5(con: duckdb.DuckDBPyConnection, query: "pd.DataFrame"):
    """query columns: variant_key, gene, protein_pos (int), alt_aa1 (1-letter, from VEP Amino_acids).
    Returns variant_key + clinvar_same_aa (PS1) + clinvar_same_codon_lp (PM5), both 0/1."""
    con.register("_q", query)
    # A qualifying ClinVar source must be a real MISSENSE (ref_aa <> alt_aa, and not a nonsense '*') and a
    # DIFFERENT nucleotide variant than the query — else PS1 is circular (the variant's own assertion) and PM5
    # could count a nonsense at the codon.
    out = con.execute("""
        SELECT q.variant_key,
               -- the QUERY must itself be a missense too (a stop-gain '*' at a codon must not get PS1/PM5).
               CASE WHEN q.alt_aa1 <> '*' AND EXISTS (SELECT 1 FROM clinvar_prot c
                       WHERE c.gene=q.gene AND c.protein_pos=q.protein_pos AND c.alt_aa1=q.alt_aa1
                         AND c.variant_key <> q.variant_key
                         AND c.ref_aa1 <> c.alt_aa1 AND c.alt_aa1 <> '*' AND c.ref_aa1 <> '*') THEN 1 ELSE 0 END AS clinvar_same_aa,
               CASE WHEN q.alt_aa1 <> '*' AND EXISTS (SELECT 1 FROM clinvar_prot c
                       WHERE c.gene=q.gene AND c.protein_pos=q.protein_pos AND c.alt_aa1<>q.alt_aa1
                         AND c.variant_key <> q.variant_key
                         AND c.ref_aa1 <> c.alt_aa1 AND c.alt_aa1 <> '*' AND c.ref_aa1 <> '*') THEN 1 ELSE 0 END AS clinvar_same_codon_lp
        FROM _q q
    """).df()
    con.unregister("_q")
    return out


def pm1_hotspot(con: duckdb.DuckDBPyConnection, query: "pd.DataFrame", *, window: int = 3, min_plp: int = 2):
    """PM1 (mutational hotspot / critical region). query columns: variant_key, gene, protein_pos. A residue is
    a hotspot if within +/- `window` codons there are >= `min_plp` pathogenic ClinVar missense AND ZERO benign
    (the 'without benign variation' clause is literal). Returns variant_key + pm1 (0/1). Empirical, data-driven
    from clinvar_aa — the ClinGen VCEP domain/hotspot specs (CSpec) can override/augment as more rows."""
    con.register("_qm", query)
    out = con.execute(f"""
        SELECT q.variant_key,
          CASE WHEN
            (SELECT COALESCE(sum(a.n_plp),0) FROM clinvar_aa a
               WHERE a.gene=q.gene AND abs(a.protein_pos - q.protein_pos) <= {window}) >= {min_plp}
            AND (SELECT COALESCE(sum(a.n_blb),0) FROM clinvar_aa a
               WHERE a.gene=q.gene AND abs(a.protein_pos - q.protein_pos) <= {window}) = 0
          THEN 1 ELSE 0 END AS pm1
        FROM _qm q
    """).df()
    con.unregister("_qm")
    return out


def ba1_exceptions_from_clinvar(con: duckdb.DuckDBPyConnection, clinvar_vcf: str, *, min_af: float = 0.01):
    """Systematic BA1/BS1 exception list — the Ghosh-2018 class GENERALISED beyond rs334: well-reviewed
    (>=2-star: expert panel / practice guideline / multiple concordant submitters) ClinVar P/LP variants that
    are COMMON (AF > min_af). These are established-pathogenic-despite-frequency, so BA1/BS1 must NOT auto-
    benign them. `min_af` defaults to the BS1 threshold (0.01), NOT BA1 (0.05): founder pathogenic variants
    like rs334 (HbS, global AF ~2.7%) trip BS1 not BA1, because the ClinVar VCF's ESP/EXAC/1000G AF is a
    GLOBAL frequency — rs334 is ~10-25% in African ancestry but diluted globally. So this is a global-AF
    PROXY; max-population/FAF (gnomAD grpmax) would be sharper. Data-driven from the ClinVar VCF via duckhts.
    Returns variant_key, gene, note — pass to classify(ba1_exceptions=...). ClinVar-derived: production only,
    NOT the leakage-free harness."""
    con.execute("LOAD duckhts")
    return con.execute(f"""
        WITH v AS (
          SELECT replace(CAST(CHROM AS VARCHAR),'chr','')||'-'||CAST(POS AS VARCHAR)||'-'||REF||'-'||ALT[1] AS variant_key,
                 split_part(CAST(INFO_GENEINFO AS VARCHAR),':',1) AS gene,
                 array_to_string(INFO_CLNSIG,';') AS sig, array_to_string(INFO_CLNREVSTAT,';') AS rev,
                 greatest(coalesce(TRY_CAST(CAST(INFO_AF_ESP AS VARCHAR) AS DOUBLE),0),
                          coalesce(TRY_CAST(CAST(INFO_AF_EXAC AS VARCHAR) AS DOUBLE),0),
                          coalesce(TRY_CAST(CAST(INFO_AF_TGP AS VARCHAR) AS DOUBLE),0)) AS af
          FROM read_bcf('{clinvar_vcf}') WHERE len(ALT)=1 AND len(REF)=1
        )
        SELECT variant_key, gene,
               'ClinVar >=2-star P/LP but common (AF '||round(af,3)||' > '||{min_af}||') — Ghosh-2018 BA1 exception' AS note
        FROM v
        WHERE sig ILIKE '%pathogenic%' AND sig NOT ILIKE '%conflict%' AND sig NOT ILIKE '%benign%'
          AND af > {min_af}
          AND (rev ILIKE '%expert%panel%' OR rev ILIKE '%practice%guideline%'
               OR (rev ILIKE '%multiple%submitters%' AND rev NOT ILIKE '%conflict%'))
    """).df()
