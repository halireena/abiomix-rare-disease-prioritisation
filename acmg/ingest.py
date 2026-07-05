"""Ingest abstraction: ANY variant source -> ONE canonical `sample_call` table (+ a `pedigree` table).

Downstream code (family logic, kernel) never touches the source format. Only this module knows raw names.
Everything is SQL over DuckDB — it scales to the whole cohort without loading millions of rows into pandas.

Sources supported by `ingest_variants`:
  - tidy parquet / TSV / CSV / Excel : a CONFIGURABLE column mapping (default = the challenge bundle names).
  - multi-sample VCF / BCF           : via duckhts `read_bcf` (reads .vcf/.vcf.gz/.bcf), unpivoted to one
                                       row per sample x variant. FORMAT GT/GQ/DP/PS are parsed per sample.
Pedigree by `ingest_pedigree`:
  - a PLINK PED/FAM (6-col)          : parsed to family_id/proband_id/mother_id/father_id directly.
  - the curated relationship fields  : delegated to acmg.family.build_pedigree over `sample_call`.
"""
from __future__ import annotations
import os
import duckdb

# Canonical per-call columns the rest of the package relies on. KEEP THIS STABLE.
CANONICAL = [
    "case_id", "family_id", "member_id", "proband_status", "relationship", "family_size",
    "variant_kind", "variant_key", "chrom", "pos", "ref", "alt", "gene",
    "gt", "gq", "dp", "phase_set",
]

# Default mapping canonical-column -> SOURCE column name, for the challenge tidy bundle. Override any
# subset for a differently-named tabular source (e.g. {"chrom": "Chromosome", "gt": "genotype"}).
BUNDLE_COLUMNS = {
    "case_id": "student_case_id",
    "family_id": "student_family_id",
    "member_id": "student_member_id",
    "proband_status": "proband_status_curated",
    "relationship": "family_relationship_curated",
    "family_size": "family_size",
    "variant_kind": "variant_kind",
    "chrom": "CHROM",
    "pos": "POS",
    "ref": "REF",
    "alt": "ALT",
    "gt": "FORMAT_GT",
    "gq": "FORMAT_GQ",
    "dp": "FORMAT_DP",
    "phase_set": "FORMAT_PS",
}

# canonical fields that must resolve for a variant call; the rest are optional (NULL when absent).
_REQUIRED = ("chrom", "pos", "ref", "alt", "gt")

# canonical column -> SQL type, so an ABSENT column becomes a correctly-typed NULL (a bare NULL is INT32,
# which then breaks a later join/UPDATE against a VARCHAR id).
_TYPES = {
    "case_id": "VARCHAR", "family_id": "VARCHAR", "member_id": "VARCHAR", "proband_status": "VARCHAR",
    "relationship": "VARCHAR", "family_size": "INTEGER", "variant_kind": "VARCHAR", "chrom": "VARCHAR",
    "pos": "BIGINT", "ref": "VARCHAR", "alt": "VARCHAR", "gt": "VARCHAR",
    "gq": "INTEGER", "dp": "INTEGER", "phase_set": "INTEGER",
}


# --------------------------------------------------------------------------- format detection / readers
def _detect_fmt(source: str) -> str:
    s = source.lower()
    for suf, fmt in ((".bcf", "vcf"), (".vcf", "vcf"), (".vcf.gz", "vcf"),
                     (".parquet", "parquet"), (".pq", "parquet"),
                     (".tsv", "tsv"), (".txt", "tsv"), (".csv", "csv"),
                     (".xlsx", "excel"), (".xls", "excel")):
        if s.endswith(suf):
            return fmt
    raise ValueError(f"cannot infer format from {source!r}; pass fmt= explicitly")


def _table_sql(con: duckdb.DuckDBPyConnection, source: str, fmt: str) -> str:
    """A SQL table expression for a tabular source. Excel prefers duckdb's reader, else pandas+openpyxl."""
    if fmt == "parquet":
        return f"read_parquet('{source}')"
    if fmt in ("tsv", "csv"):
        delim = "\\t" if fmt == "tsv" else ","
        return f"read_csv('{source}', delim='{delim}', header=true, auto_detect=true)"
    if fmt == "excel":
        try:
            con.execute("INSTALL excel; LOAD excel;")
            con.execute(f"DESCRIBE SELECT * FROM read_xlsx('{source}')")
            return f"read_xlsx('{source}')"
        except Exception:
            import pandas as pd  # boundary-only: hand the sheet to DuckDB as a registered view
            df = pd.read_excel(source)
            con.register("_excel_src", df)
            return "_excel_src"
    raise ValueError(f"not a tabular format: {fmt}")


def _column_types(con: duckdb.DuckDBPyConnection, table_sql: str) -> dict[str, str]:
    rows = con.execute(f"DESCRIBE SELECT * FROM {table_sql}").fetchall()
    return {r[0]: r[1] for r in rows}


# --------------------------------------------------------------------------- canonical column expressions
def _norm_expr(canonical: str, col: str, coltype: str) -> str:
    """SQL that turns a source column into its canonical, normalized value."""
    q = f'"{col}"'
    if canonical == "chrom":
        return f"replace(CAST({q} AS VARCHAR), 'chr', '')"
    if canonical == "pos":
        return f"CAST({q} AS BIGINT)"
    if canonical in ("ref", "alt"):
        base = f"{q}[1]" if coltype.endswith("[]") else q  # VCF/bundle ALT is a list; TSV alt is scalar
        return f"CAST({base} AS VARCHAR)"
    if canonical in ("gq", "dp", "phase_set", "family_size"):
        return f"TRY_CAST({q} AS INTEGER)"
    return f"CAST({q} AS VARCHAR)"  # ids, labels, gt, variant_kind


def _canonical_select(colmap: dict[str, str], types: dict[str, str]) -> str:
    """Build the inner SELECT list: each canonical field -> normalized source expr, or NULL if absent."""
    parts = []
    for canonical in ("case_id", "family_id", "member_id", "proband_status", "relationship",
                      "family_size", "variant_kind", "chrom", "pos", "ref", "alt",
                      "gt", "gq", "dp", "phase_set"):
        col = colmap.get(canonical)
        if col and col in types:
            parts.append(f"{_norm_expr(canonical, col, types[col])} AS {canonical}")
        elif canonical in _REQUIRED:
            raise ValueError(f"required column {canonical!r} (mapped to {col!r}) not found in source")
        else:
            parts.append(f"CAST(NULL AS {_TYPES[canonical]}) AS {canonical}")
    return ",\n            ".join(parts)


# CTE that lifts an inner normalized SELECT (columns: the canonical fields except variant_key/gene) to
# the FULL canonical schema: derive variant_kind when absent, derive variant_key, add gene=NULL.
def _finalize_sql(inner: str, snv_only: bool) -> str:
    kind = ("COALESCE(variant_kind, CASE WHEN length(ref)=1 AND length(alt)=1 THEN 'snv' ELSE 'indel' END)"
            " AS variant_kind")
    key = ("chrom || '-' || CAST(pos AS VARCHAR) || '-' || ref || '-' || alt AS variant_key")
    body = f"""
        WITH _n AS (
            {inner}
        ),
        _k AS (
            SELECT case_id, family_id, member_id, proband_status, relationship, family_size,
                   {kind},
                   {key},
                   chrom, pos, ref, alt,
                   CAST(NULL AS VARCHAR) AS gene,
                   gt, gq, dp, phase_set
            FROM _n
        )
        SELECT {', '.join(CANONICAL)} FROM _k
    """
    if snv_only:
        body += " WHERE variant_kind = 'snv'"
    return body


# --------------------------------------------------------------------------- public: variants
def ingest_variants(con: duckdb.DuckDBPyConnection, source: str, fmt: str | None = None,
                    columns: dict[str, str] | None = None, *, snv_only: bool = True,
                    table: str = "sample_call") -> None:
    """Create the canonical `sample_call` table from ANY supported source.

    source   : a path/glob DuckDB can read.
    fmt      : 'parquet'|'tsv'|'csv'|'excel'|'vcf' (auto-detected from the extension when None).
    columns  : canonical->source column-name map for TABULAR sources (default BUNDLE_COLUMNS). Ignored
               for VCF (sample columns are discovered from the file). Coordinates stay on their native
               build (no liftover); chr-prefix is stripped; ALT lists are indexed to the first allele.
    """
    fmt = fmt or _detect_fmt(source)
    if fmt in ("vcf", "bcf"):
        inner = _vcf_inner(con, source)
    else:
        colmap = {**BUNDLE_COLUMNS, **(columns or {})}
        table_sql = _table_sql(con, source, fmt)
        types = _column_types(con, table_sql)
        inner = f"SELECT {_canonical_select(colmap, types)}\n            FROM {table_sql}"
    con.execute(f"CREATE OR REPLACE TABLE {table} AS " + _finalize_sql(inner, snv_only))


def _vcf_inner(con: duckdb.DuckDBPyConnection, source: str) -> str:
    """Inner normalized SELECT for a multi-sample VCF/BCF: unpivot duckhts `read_bcf`'s wide per-sample
    FORMAT_<field>_<sample> columns to one row per sample x variant. The VCF has no case/family metadata,
    so the sample id doubles as member_id AND case_id; pedigree comes from a PED (see ingest_pedigree)."""
    con.execute("INSTALL duckhts FROM community; LOAD duckhts;")
    types = _column_types(con, f"read_bcf('{source}')")
    samples = [c[len("FORMAT_GT_"):] for c in types if c.startswith("FORMAT_GT_")]
    if not samples:
        raise ValueError(f"no FORMAT/GT sample columns found in {source!r}")

    def fmt_col(field: str, sample: str, cast_int: bool) -> str:
        col = f"FORMAT_{field}_{sample}"
        if col not in types:
            return "NULL"
        return f'TRY_CAST("{col}" AS INTEGER)' if cast_int else f'CAST("{col}" AS VARCHAR)'

    blocks = []
    for s in samples:
        lit = s.replace("'", "''")
        blocks.append(f"""
            SELECT
                '{lit}' AS case_id,
                CAST(NULL AS VARCHAR) AS family_id,
                '{lit}' AS member_id,
                CAST(NULL AS VARCHAR) AS proband_status,
                CAST(NULL AS VARCHAR) AS relationship,
                CAST(NULL AS INTEGER) AS family_size,
                CAST(NULL AS VARCHAR) AS variant_kind,
                replace(CAST(CHROM AS VARCHAR), 'chr', '') AS chrom,
                CAST(POS AS BIGINT) AS pos,
                CAST(REF AS VARCHAR) AS ref,
                CAST(ALT[1] AS VARCHAR) AS alt,
                {fmt_col('GT', s, False)} AS gt,
                {fmt_col('GQ', s, True)} AS gq,
                {fmt_col('DP', s, True)} AS dp,
                {fmt_col('PS', s, True)} AS phase_set
            FROM read_bcf('{source}')""")
    return "\n            UNION ALL\n".join(blocks)


def load_bundle(con: duckdb.DuckDBPyConnection, source: str, *, snv_only: bool = True) -> None:
    """Thin wrapper: ingest the challenge tidy parquet bundle with its default column mapping."""
    ingest_variants(con, source, fmt="parquet", columns=BUNDLE_COLUMNS, snv_only=snv_only)


# --------------------------------------------------------------------------- public: pedigree
def ingest_pedigree(con: duckdb.DuckDBPyConnection, source: str | None = None) -> None:
    """Create the `pedigree` table (family_id, proband_id, mother_id, father_id[, family_size]).

    source is a PED/FAM path -> parse the 6-column PLINK pedigree. source is None -> derive from the
    curated relationship/proband_status fields already on `sample_call` (delegates to family.build_pedigree).
    """
    if source is None:
        from .family import build_pedigree
        build_pedigree(con)
        return
    _pedigree_from_ped(con, source)


def _pedigree_from_ped(con: duckdb.DuckDBPyConnection, path: str) -> None:
    """PLINK PED/FAM: family_id, individual_id, father_id, mother_id, sex, phenotype (whitespace-delim,
    no header, '0' = missing parent, phenotype '2' = affected). The proband is the affected child (has a
    parent); mother/father are read from that child's parent columns. Also backfills sample_call.family_id."""
    delim = "\\t" if path.lower().endswith(".tsv") else " "  # PLINK .ped/.fam is space/tab-delimited
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE _ped AS
        SELECT
            CAST(column0 AS VARCHAR) AS family_id,
            CAST(column1 AS VARCHAR) AS individual_id,
            NULLIF(CAST(column2 AS VARCHAR), '0') AS father_id,
            NULLIF(CAST(column3 AS VARCHAR), '0') AS mother_id,
            CAST(column5 AS VARCHAR) AS phenotype
        FROM read_csv('{path}', delim='{delim}', header=false, columns={{
            'column0':'VARCHAR','column1':'VARCHAR','column2':'VARCHAR',
            'column3':'VARCHAR','column4':'VARCHAR','column5':'VARCHAR'}})
    """)
    con.execute("""
        CREATE OR REPLACE TABLE pedigree AS
        WITH ranked AS (
            SELECT *, row_number() OVER (
                PARTITION BY family_id
                ORDER BY (father_id IS NOT NULL OR mother_id IS NOT NULL) DESC,
                         (phenotype = '2') DESC, individual_id
            ) AS rn
            FROM _ped
        )
        SELECT family_id, individual_id AS proband_id, mother_id, father_id
        FROM ranked WHERE rn = 1
    """)
    # backfill family_id onto sample_call so the cohort groups coherently (VCF ingest left it NULL).
    con.execute("""
        UPDATE sample_call SET family_id = p.family_id
        FROM _ped p WHERE p.individual_id = sample_call.member_id
    """)


# --------------------------------------------------------------------------- annotation helper (unchanged)
def unique_variants_for_annotation(con: duckdb.DuckDBPyConnection) -> duckdb.DuckDBPyRelation:
    """The DISTINCT loci to annotate ONCE (dedup across the whole cohort), so the annotator is never
    called twice for the same variant. Pre-filter (rarity + coding/splice + panel) BEFORE this to keep
    the count in the hundreds, not tens of thousands per exome."""
    return con.sql("SELECT DISTINCT chrom, pos, ref, alt, variant_key FROM sample_call")
