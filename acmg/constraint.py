"""gnomAD gene constraint for PP2 — the missense-constraint Z score (`mis_z`) is just a per-gene table.

PP2 (the kernel) fires for a missense variant when its gene is missense-constrained: gnomAD v2 constraint
`mis_z >= 3.09`. gnomAD v2 is GRCh37-native, so this pairs with the rest of the pipeline with no liftover.

Source (download once): gnomAD v2.1.1 gene-level LoF metrics —
  https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
(bgzip is gzip-compatible; DuckDB reads it directly.)
"""
from __future__ import annotations
import duckdb


def load_constraint(con: duckdb.DuckDBPyConnection, by_gene_gz: str) -> None:
    """Build `gene_constraint(gene, mis_z, oe_lof_upper, pli)` from the gnomAD by-gene metrics."""
    con.execute(f"""
        CREATE OR REPLACE TABLE gene_constraint AS
        SELECT gene,
               TRY_CAST(mis_z AS DOUBLE)          AS mis_z,
               TRY_CAST(oe_lof_upper AS DOUBLE)   AS oe_lof_upper,
               TRY_CAST(pLI AS DOUBLE)            AS pli
        FROM read_csv('{by_gene_gz}', delim='\t', header=true, compression='gzip', ignore_errors=true,
                      types={{'gene':'VARCHAR','mis_z':'VARCHAR','oe_lof_upper':'VARCHAR','pLI':'VARCHAR'}})
        WHERE gene IS NOT NULL
        QUALIFY row_number() OVER (PARTITION BY gene ORDER BY mis_z DESC NULLS LAST) = 1
    """)


def add_mis_z(con: duckdb.DuckDBPyConnection, annotations: "pd.DataFrame") -> "pd.DataFrame":
    """Left-join `gnomad_mis_z` (PP2) AND `loeuf` (oe_lof_upper, the LEAKAGE-FREE LoF-intolerance metric for
    the PVS1 gate) onto an annotations DataFrame by gene (requires load_constraint first). A gene absent from
    the table stays NULL and the corresponding criterion abstains."""
    con.register("_ann", annotations)
    out = con.execute("""
        SELECT a.* REPLACE (COALESCE(gc.mis_z, a.gnomad_mis_z) AS gnomad_mis_z)
        FROM _ann a LEFT JOIN gene_constraint gc ON gc.gene = a.gene
    """).df()
    con.unregister("_ann")
    lo = dict(con.execute("SELECT gene, oe_lof_upper FROM gene_constraint WHERE oe_lof_upper IS NOT NULL").fetchall())
    out["loeuf"] = out["gene"].map(lo)   # by gene; pandas map so it works whether or not loeuf pre-exists
    return out
