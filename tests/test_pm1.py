"""PM1 mutational-hotspot clustering: >=min_plp pathogenic within +/-window codons AND zero benign there
('without benign variation'). Run: PYTHONPATH=. python3 tests/test_pm1.py"""
import duckdb, pandas as pd
from acmg.clinvar import pm1_hotspot


def _con(rows):
    con = duckdb.connect()
    con.register("_a", pd.DataFrame(rows, columns=["gene", "protein_pos", "n_plp", "n_blb"]))
    con.execute("CREATE TABLE clinvar_aa AS SELECT * FROM _a"); con.unregister("_a")
    return con


def test_hotspot_fires_and_benign_blocks():
    # G@100: 3 P/LP within +/-3, no benign -> PM1. G@200: pathogenic cluster BUT a benign at 201 -> blocked.
    con = _con([("G", 100, 2, 0), ("G", 101, 1, 0), ("G", 200, 3, 0), ("G", 201, 0, 1)])
    q = pd.DataFrame([dict(variant_key="a", gene="G", protein_pos=100),
                      dict(variant_key="b", gene="G", protein_pos=200),
                      dict(variant_key="c", gene="G", protein_pos=500)])
    r = pm1_hotspot(con, q).set_index("variant_key")["pm1"].to_dict()
    assert r["a"] == 1, "hotspot with no benign should fire PM1"
    assert r["b"] == 0, "benign variation in the window must block PM1"
    assert r["c"] == 0, "cold position: no PM1"


if __name__ == "__main__":
    test_hotspot_fires_and_benign_blocks(); print("all passed")
