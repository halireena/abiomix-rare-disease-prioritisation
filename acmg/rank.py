"""Phenotype x genotype reranking, Exomiser-style (transparent, not a black box).

Exomiser combines a PHENOTYPE score (patient HPO vs the gene's known phenotypes) with a GENOTYPE /
pathogenicity score, per gene. We do the same with readable parts:

  - phenotype score: over a `gene_phenotype(gene, hpo_id)` edge table (populate from Monarch, or a
    fixture). **IC-weighted** by default: a specific HPO term (few genes) counts more than a generic one
    (many genes). This is the fix for the classic failure where a huge pleiotropic gene (e.g. TTN) wins
    just by being annotated to many generic terms — raw counts over-reward gene size; IC does not.
  - genotype score: the ACMG class from our kernel, mapped to [0,1].
  - inheritance score: a small boost for de novo / biallelic (from acmg.family).

combined = wp*phenotype_norm + wg*pathogenicity + wi*inheritance  (weights explicit and tunable).
The point is auditability: every term is a number you can trace, not a learned weight.
"""
from __future__ import annotations
import math
import duckdb
import pandas as pd

PATHOGENICITY = {
    "Pathogenic": 1.0, "Likely Pathogenic": 0.9, "VUS": 0.5,
    "Likely Benign": 0.1, "Benign": 0.0,
}
DEFAULT_WEIGHTS = {"phenotype": 0.5, "genotype": 0.4, "inheritance": 0.1}


def ic_weights(con: duckdb.DuckDBPyConnection, hpo_terms: list[str], table: str = "gene_phenotype") -> dict:
    """Information content per term: -log2(genes_with_term / total_genes). Specific term -> high IC."""
    if not hpo_terms:
        return {}
    # IC denominator = the FULL gene universe of the KG, not just genes linked to this case's terms. When
    # gene_phenotype was materialised for only the case's HPO (monarch_gene_phenotype), count(DISTINCT gene) over
    # it is the matched subset -> deflated IC. Prefer the recorded universe size if present (fixtures self-contain).
    has_universe = con.execute(
        "SELECT count(*) FROM information_schema.tables WHERE table_name = 'gene_phenotype_universe'").fetchone()[0]
    total = (con.execute("SELECT total_genes FROM gene_phenotype_universe").fetchone()[0] if has_universe
             else con.execute(f"SELECT count(DISTINCT gene) FROM {table}").fetchone()[0]) or 1
    q = ",".join(f"'{t}'" for t in hpo_terms)
    rows = con.execute(
        f"SELECT hpo_id, count(DISTINCT gene) AS g FROM {table} WHERE hpo_id IN ({q}) GROUP BY hpo_id"
    ).df()
    return {r.hpo_id: -math.log2(max(r.g, 1) / total) for r in rows.itertuples()}


def phenotype_scores(con: duckdb.DuckDBPyConnection, hpo_terms: list[str], *,
                     ic_weighted: bool = True, table: str = "gene_phenotype") -> pd.DataFrame:
    """Per-gene phenotype score = sum over the patient's HPO terms linked to that gene of the term weight
    (IC if ic_weighted, else 1). Returns columns: gene, phenotype_score."""
    if not hpo_terms:
        return pd.DataFrame(columns=["gene", "phenotype_score"])
    w = ic_weights(con, hpo_terms, table) if ic_weighted else {t: 1.0 for t in hpo_terms}
    con.register("_w", pd.DataFrame({"hpo_id": list(w), "wt": list(w.values())}))
    q = ",".join(f"'{t}'" for t in hpo_terms)
    out = con.execute(f"""
        SELECT gp.gene, sum(w.wt) AS phenotype_score
        FROM {table} gp JOIN _w w ON w.hpo_id = gp.hpo_id
        WHERE gp.hpo_id IN ({q})
        GROUP BY gp.gene
    """).df()
    con.unregister("_w")
    return out


def rank_fn(con: duckdb.DuckDBPyConnection, table: str = "gene_phenotype", *, ic_weighted: bool = True):
    """A rank_fn(candidate_df, hpo) for acmg.repl.CaseSession / acmg.proband — maps each candidate's gene to its
    IC-weighted phenotype score over the `gene_phenotype` table (Monarch KG in production, a fixture in tests).
    `hpo`: the session's live
    HPO set (HP ids), so refining HPO deterministically re-ranks. Genes with no phenotype link score 0."""
    def _fn(df: pd.DataFrame, hpo) -> list:
        terms = list(hpo) if hpo else []
        ps = phenotype_scores(con, terms, ic_weighted=ic_weighted, table=table)
        m = dict(zip(ps["gene"], ps["phenotype_score"]))
        return [float(m.get(g, 0.0)) for g in df["gene"]]
    return _fn


def rerank(variants: pd.DataFrame, phenotype: pd.DataFrame, *,
           inheritance: dict | None = None, weights: dict | None = None) -> pd.DataFrame:
    """Combine per-variant ACMG class (`acmg_class`, `gene`) with the per-gene phenotype score and an
    optional per-variant_key inheritance boost in [0,1]. Returns variants + scores, sorted desc."""
    w = {**DEFAULT_WEIGHTS, **(weights or {})}
    v = variants.merge(phenotype, on="gene", how="left")
    v["phenotype_score"] = v["phenotype_score"].fillna(0.0)
    # normalise phenotype to [0,1] across this case's candidates (relative match strength)
    pmax = v["phenotype_score"].max()
    v["phenotype_norm"] = v["phenotype_score"] / pmax if pmax and pmax > 0 else 0.0
    v["pathogenicity"] = v["acmg_class"].map(PATHOGENICITY).fillna(0.5)
    inh = inheritance or {}
    v["inheritance_score"] = v["variant_key"].map(lambda k: inh.get(k, 0.0))
    v["combined_score"] = (
        w["phenotype"] * v["phenotype_norm"]
        + w["genotype"] * v["pathogenicity"]
        + w["inheritance"] * v["inheritance_score"]
    )
    return v.sort_values("combined_score", ascending=False).reset_index(drop=True)


# --- Monarch helpers (network) — populate gene_phenotype from the Monarch KG, like the sprint pipeline ---

def attach_monarch(con: duckdb.DuckDBPyConnection) -> None:
    # IDEMPOTENT: a second ATTACH of 'monarch' throws, and a per-family caller's broad except then blanks
    # phenotype for every family after the first (pi review #2). Only attach if not already attached.
    con.execute("INSTALL httpfs; LOAD httpfs;")
    if not con.execute("SELECT count(*) FROM duckdb_databases() WHERE database_name='monarch'").fetchone()[0]:
        con.execute("ATTACH 'https://data.monarchinitiative.org/monarch-kg/latest/monarch-kg.duckdb' AS monarch (READ_ONLY)")


def monarch_gene_phenotype(con: duckdb.DuckDBPyConnection, hpo_terms: list[str]) -> None:
    """Materialise `gene_phenotype(gene, hpo_id)` for the patient's HPO terms from Monarch (direct
    gene->HPO edges). Requires attach_monarch(con) first."""
    q = ",".join(f"'{t}'" for t in hpo_terms)
    con.execute(f"""
        CREATE OR REPLACE TABLE gene_phenotype AS
        SELECT DISTINCT subject_label AS gene, object AS hpo_id
        FROM monarch.denormalized_edges
        WHERE predicate = 'biolink:has_phenotype' AND subject_category = 'biolink:Gene'
          AND object IN ({q})
    """)
    # record the FULL gene universe (all phenotype-annotated genes) so ic_weights' denominator is the KG, not the
    # per-case subset — otherwise a one-term case makes the denominator == numerator and every IC collapses to 0.
    con.execute("""
        CREATE OR REPLACE TABLE gene_phenotype_universe AS
        SELECT count(DISTINCT subject_label) AS total_genes
        FROM monarch.denormalized_edges
        WHERE predicate = 'biolink:has_phenotype' AND subject_category = 'biolink:Gene'
    """)
