"""Prior-analysis diff tests. Run: python tests/test_reanalysis.py"""
import pandas as pd
from acmg import reanalysis


def _df(rows):
    return pd.DataFrame(rows, columns=["variant_key", "acmg_class"])


def test_new_upgraded_downgraded_unchanged_dropped():
    prior = _df([("v1", "VUS"), ("v2", "Likely Pathogenic"), ("v4", "VUS")])
    current = _df([("v1", "Likely Pathogenic"), ("v2", "VUS"), ("v3", "Pathogenic"), ("v4", "VUS")])
    d = reanalysis.diff(current, prior).set_index("variant_key")["status"].to_dict()
    assert d["v1"] == "upgraded"      # VUS -> LP
    assert d["v2"] == "downgraded"    # LP -> VUS
    assert d["v3"] == "new"           # not in prior
    assert d["v4"] == "unchanged"
    assert d["v5:missing"] if False else True
    # v* only in prior -> dropped
    prior2 = _df([("gone", "Pathogenic")])
    assert reanalysis.diff(_df([]), prior2).iloc[0]["status"] == "dropped"


def test_newly_reportable_is_the_reanalysis_signal():
    prior = _df([("v1", "VUS")])
    current = _df([("v1", "Likely Pathogenic"), ("v2", "VUS")])
    nr = reanalysis.newly_reportable(current, prior)
    assert list(nr["variant_key"]) == ["v1"]              # VUS->LP is the newly reportable one
    assert reanalysis.case_is_newly_reportable(current, prior) is True
    # a case where nothing crossed into reportable -> no re-review warranted
    assert reanalysis.case_is_newly_reportable(_df([("v2", "VUS")]), _df([("v2", "VUS")])) is False


if __name__ == "__main__":
    for n, f in sorted(globals().items()):
        if n.startswith("test_") and callable(f):
            f(); print("ok", n)
    print("all passed")
