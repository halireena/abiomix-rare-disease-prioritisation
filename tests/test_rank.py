"""IC weighting — the denominator must be the KG's full gene universe, not the per-case matched subset."""
import math
import duckdb
import pandas as pd
from acmg import rank


def _con(edges, universe=None):
    con = duckdb.connect()
    con.register("_e", pd.DataFrame(edges, columns=["gene", "hpo_id"]))
    con.execute("CREATE TABLE gene_phenotype AS SELECT * FROM _e")
    if universe is not None:
        con.execute(f"CREATE TABLE gene_phenotype_universe AS SELECT {universe} AS total_genes")
    return con


def test_ic_uses_recorded_universe_not_subset():
    # gene_phenotype holds only the case's matched edges (2 genes), but the KG has 20000 genes.
    con = _con([("AAA", "HP:1"), ("BBB", "HP:1")], universe=20000)
    w = rank.ic_weights(con, ["HP:1"])
    # IC = -log2(genes_with_term / total). Subset would give -log2(2/2)=0 (collapsed); universe gives a real value.
    assert w["HP:1"] == math.log2(20000 / 2)
    assert w["HP:1"] > 10


def test_ic_falls_back_to_table_when_no_universe():
    con = _con([("AAA", "HP:1"), ("BBB", "HP:2")])   # fixture is self-contained -> subset IS the universe
    w = rank.ic_weights(con, ["HP:1"])
    assert w["HP:1"] == math.log2(2 / 1)             # 2 genes total, 1 has the term
