"""Run the ACMG kernel (the SQL manifests) from Python via DuckDB.

The rules ARE the SQL in `manifests/` — this module does not re-implement them, it runs them.
One source of truth, mirroring Abiomix's internal production system. `classify()` takes a
DataFrame in the `annotations` schema (see acmg.vep_map.REQUIRED_COLS) and returns one ACMG
classification per variant, with the criteria that fired and the Tavtigian points.
"""
from __future__ import annotations
from pathlib import Path
import duckdb
import pandas as pd
from .vep_map import REQUIRED_COLS

_MANIFEST_DIR = Path(__file__).parent / "manifests"
_STAGES = ["00_spec.sql", "01_acmg_rules.sql", "02_combine.sql"]


def _sql(name: str) -> str:
    return (_MANIFEST_DIR / name).read_text()


def _override(con: duckdb.DuckDBPyConnection, table: str, df: pd.DataFrame | None) -> None:
    """Replace an illustrative spec table (gene_curation / acmg_thresholds / variant_frequency /
    ba1_exceptions) with real curation. DELETE+INSERT keeps the table identity so views stay valid."""
    if df is None:
        return
    con.register("_override_df", df)
    con.execute(f"DELETE FROM {table}")
    con.execute(f"INSERT INTO {table} BY NAME SELECT * FROM _override_df")
    con.unregister("_override_df")


def classify(
    annotations: pd.DataFrame,
    con: duckdb.DuckDBPyConnection | None = None,
    gene_curation: pd.DataFrame | None = None,
    acmg_thresholds: pd.DataFrame | None = None,
    variant_frequency: pd.DataFrame | None = None,
    ba1_exceptions: pd.DataFrame | None = None,
    pm1_enabled: bool = False,
    pvs1_constraint: bool = False,
    cadd_supporting: bool = False,
) -> pd.DataFrame:
    """Classify variants with the ACMG kernel.

    annotations: DataFrame in the kernel schema (acmg.vep_map.REQUIRED_COLS). Missing optional
                 columns are filled with NULL and the kernel abstains on them.
    gene_curation etc.: pass real ClinGen/ACGS curation to override the illustrative defaults baked
                 into manifests/00_spec.sql. Without real gene_curation, PVS1 abstains (a null variant
                 in an uncurated gene gets no PVS1) and the gene-disease validity cap is a no-op —
                 fail-safe, not silently wrong.

    Returns columns: variant_key, gene, variant_kind, total_points, criteria, acmg_class.
    """
    con = con or duckdb.connect()
    ann = annotations.copy()
    for col in REQUIRED_COLS:
        if col not in ann.columns:
            ann[col] = None
    ann = ann[REQUIRED_COLS]

    # Cast to a fixed schema: an all-NULL column would otherwise be inferred as DOUBLE and break the
    # string comparisons in the rules (e.g. `consequence IN ('missense', ...)`).
    con.register("_ann_in", ann)
    con.execute("""CREATE OR REPLACE TABLE annotations AS SELECT
        CAST(variant_key AS VARCHAR)            AS variant_key,
        CAST(gene AS VARCHAR)                   AS gene,
        CAST(consequence AS VARCHAR)            AS consequence,
        CAST(variant_kind AS VARCHAR)           AS variant_kind,
        CAST(filtering_af AS DOUBLE)            AS filtering_af,
        CAST(gnomad_mis_z AS DOUBLE)            AS gnomad_mis_z,
        CAST(revel AS DOUBLE)                   AS revel,
        CAST(spliceai AS DOUBLE)                AS spliceai,
        CAST(clinvar_same_aa AS INTEGER)        AS clinvar_same_aa,
        CAST(clinvar_same_codon_lp AS INTEGER)  AS clinvar_same_codon_lp,
        CAST(nmd_escaping AS INTEGER)           AS nmd_escaping,
        CAST(pm1 AS INTEGER)                    AS pm1,
        CAST(loeuf AS DOUBLE)                   AS loeuf,
        CAST(cadd_phred AS DOUBLE)              AS cadd_phred,
        CAST(clinvar_classification AS VARCHAR) AS clinvar_classification
      FROM _ann_in""")
    con.unregister("_ann_in")

    con.execute(_sql("00_spec.sql"))  # spec tables + annotations_spec view (illustrative defaults)
    _override(con, "gene_curation", gene_curation)
    _override(con, "acmg_thresholds", acmg_thresholds)
    _override(con, "variant_frequency", variant_frequency)
    _override(con, "ba1_exceptions", ba1_exceptions)
    if pm1_enabled:  # PM1 is OPT-IN (contested / VCEP-specific / overlaps PP2+PP3); off by default
        con.execute("UPDATE pm1_config SET enabled = TRUE")
    if pvs1_constraint:  # OPT-IN constraint-based PVS1 (LOEUF); over-calls, off by default (ACMG-faithful)
        con.execute("UPDATE pvs1_constraint_config SET enabled = TRUE")
    if cadd_supporting:  # OPT-IN: CADD as a SUPPORTING-only computational fallback (NOT ClinGen-calibrated like
        con.execute("UPDATE cadd_config SET enabled = TRUE")  # REVEL); fires only on the no-other-evidence tail
    con.execute(_sql("01_acmg_rules.sql"))   # the 13 rules -> applied_criteria
    con.execute(_sql("02_combine.sql"))      # Tavtigian points combine -> classification

    # `pm1_hotspot` is surfaced as CONTEXT for the agent/reviewer even when PM1 is not scored (opt-in off).
    return con.execute(
        "SELECT c.variant_key, c.gene, c.variant_kind, c.total_points, c.criteria, c.acmg_class, "
        "       COALESCE(a.pm1, 0) AS pm1_hotspot "
        "FROM classification c LEFT JOIN annotations a USING (variant_key) "
        "ORDER BY c.total_points DESC, c.variant_key"
    ).df()
