"""ClinGen gene curation -> the kernel's `gene_curation` table, from ClinGen's DOWNLOADABLE, versioned feeds.

This is the sourced, reproducible replacement for a hand-curated gene panel: expert per-gene-disease
knowledge (validity, mode of inheritance, haploinsufficiency) as DATA the ACMG kernel consumes, not a
hardcoded list. It answers "what genes carry expert curation, and what does PVS1/validity say" — the
legitimate content of a clinical panel (see the panel discussion) without reifying scope or fiefdoms.

Sources (stable URLs; record the FILE CREATED date for reproducibility -> docs/data_versions.md):
  - Gene-Disease Validity : https://search.clinicalgenome.org/kb/gene-validity/download  (gene, MONDO
        disease, MOI, Definitive/Strong/Moderate/Limited/...) -> gene_curation gene-disease rows.
  - Dosage Sensitivity    : https://search.clinicalgenome.org/kb/gene-dosage/download    (haploinsufficiency
        / triplosensitivity) -> hi_score + lof_mechanism.

The ClinGen CSpec registry (cspec.genome.network) publishes the per-gene VCEP *criterion* specifications
(BA1/BS1/PM2 thresholds, PVS1 decision trees). Its API is undocumented, so those per-gene frequency
thresholds are a follow-on (they land in the kernel's `acmg_thresholds` table, keyed by gene). This module
covers the reliably-downloadable core.
"""
from __future__ import annotations
import os
import urllib.request
import duckdb

GENE_VALIDITY_URL = "https://search.clinicalgenome.org/kb/gene-validity/download"
DOSAGE_URL = "https://search.clinicalgenome.org/kb/gene-dosage/download"


def download(dest_dir: str, *, force: bool = False) -> dict[str, str]:
    """Fetch the two ClinGen CSVs into dest_dir (cache). Returns {'gene_validity':path,'dosage':path}."""
    os.makedirs(dest_dir, exist_ok=True)
    out = {}
    for name, url in (("gene_validity", GENE_VALIDITY_URL), ("dosage", DOSAGE_URL)):
        path = os.path.join(dest_dir, f"clingen_{name}.csv")
        if force or not os.path.exists(path):
            urllib.request.urlretrieve(url, path)
        out[name] = path
    return out


# ClinGen files carry a 4-line preamble, then the column header, then a '++++' separator row before data.
def _read(csv_path: str) -> str:
    return (f"SELECT * FROM read_csv('{csv_path}', skip=4, header=true, all_varchar=true, quote='\"', "
            f"ignore_errors=true) WHERE \"GENE SYMBOL\" IS NOT NULL AND \"GENE SYMBOL\" NOT LIKE '+%'")


def load_gene_curation(con: duckdb.DuckDBPyConnection, gene_validity_csv: str, dosage_csv: str) -> None:
    """Build `gene_curation (gene, disease_id, mode_of_inheritance, hi_score, ts_score, lof_mechanism,
    gene_disease_validity, source)` from the two feeds. MOI is normalised to the leading token (AD/AR/XL/
    MT/SD); haploinsufficiency + triplosensitivity text -> the ClinGen HI/TS numeric scores (3 = sufficient);
    sufficient HI => lof_mechanism true."""
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE _dosage AS
        SELECT "GENE SYMBOL" AS gene, lower(trim("HAPLOINSUFFICIENCY")) AS hi_text,
               lower(trim("TRIPLOSENSITIVITY")) AS ts_text
        FROM ({_read(dosage_csv)})
    """)
    con.execute(f"""
        CREATE OR REPLACE TABLE gene_curation AS
        WITH v AS ({_read(gene_validity_csv)}),
        hi AS (
            SELECT gene,
                CASE
                    WHEN hi_text LIKE 'sufficient evidence%'  THEN 3
                    WHEN hi_text LIKE 'emerging evidence%' OR hi_text LIKE 'some evidence%' THEN 2
                    WHEN hi_text LIKE 'little evidence%'      THEN 1
                    WHEN hi_text LIKE 'no evidence%'          THEN 0
                    WHEN hi_text LIKE '%autosomal recessive%' THEN 30
                    WHEN hi_text LIKE 'dosage sensitivity unlikely%' THEN 40
                    ELSE NULL END AS hi_score,
                CASE                                   -- triplosensitivity, same ClinGen evidence banding (3 = TS)
                    WHEN ts_text LIKE 'sufficient evidence%'  THEN 3
                    WHEN ts_text LIKE 'emerging evidence%' OR ts_text LIKE 'some evidence%' THEN 2
                    WHEN ts_text LIKE 'little evidence%'      THEN 1
                    WHEN ts_text LIKE 'no evidence%'          THEN 0
                    WHEN ts_text LIKE 'dosage sensitivity unlikely%' THEN 40
                    ELSE NULL END AS ts_score
            FROM _dosage
        )
        SELECT v."GENE SYMBOL" AS gene,
               v."DISEASE ID (MONDO)" AS disease_id,
               regexp_extract(v."MOI", '^(AD|AR|XL|MT|SD)', 1) AS mode_of_inheritance,
               hi.hi_score AS hi_score,
               hi.ts_score AS ts_score,                     -- ClinGen triplosensitivity (3 = sufficient -> TS gene)
               (hi.hi_score IN (3, 30)) AS lof_mechanism,  -- LoF is a disease mechanism for DOMINANT haploinsuff
               -- (HI 3) AND for RECESSIVE genes (30 = 'associated with AR phenotype'); zygosity is enforced
               -- downstream (AR needs biallelic to be causal). 40 = dosage-insensitive -> not LoF.
               v."CLASSIFICATION" AS gene_disease_validity,
               'ClinGen-GDV+Dosage' AS source
        FROM v LEFT JOIN hi ON hi.gene = v."GENE SYMBOL"
    """)


def gene_curation_df(con: duckdb.DuckDBPyConnection, cache_dir: str, *, force: bool = False):
    """Convenience: download (cached) + build + return the gene_curation DataFrame for acmg.kernel.classify."""
    paths = download(cache_dir, force=force)
    load_gene_curation(con, paths["gene_validity"], paths["dosage"])
    return con.execute("SELECT * FROM gene_curation").df()
