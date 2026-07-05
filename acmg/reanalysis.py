"""Prior-analysis diffing — the "re" in re-analysis (the piece Talos has that a one-shot classifier lacks).

Re-analysis is only worth a clinician's time when something has CHANGED since the last look. Given the
current run and a prior result per case, this surfaces the DELTA: variants that are NEW, or whose class was
UPGRADED (e.g. VUS -> Likely Pathogenic because a gene-disease association or a ClinVar entry appeared), and
flags a case as "newly reportable" when any new/upgraded variant reaches P/LP. Everything a curator already
saw and dismissed stays out of the way.

`current` / `prior` are DataFrames with a variant key and an ACMG class column (as produced by
acmg.kernel.classify or acmg.rank.rerank). No I/O, no state — you supply the prior (from your datalake /
DuckLake snapshot of the last release).
"""
from __future__ import annotations
import pandas as pd

# ACMG 5-tier ordering, so an "upgrade" is a well-defined increase in pathogenicity.
RANK = {"Benign": 0, "Likely Benign": 1, "VUS": 2, "Likely Pathogenic": 3, "Pathogenic": 4}
REPORTABLE = {"Likely Pathogenic", "Pathogenic"}


def _rank(cls) -> int:
    return RANK.get(cls, -1)  # unknown / 'Not evaluated' / None sort below Benign


def diff(current: pd.DataFrame, prior: pd.DataFrame, *, key: str = "variant_key", cls: str = "acmg_class") -> pd.DataFrame:
    """Per-variant change vs the prior run. status ∈ {new, upgraded, downgraded, unchanged, dropped}.
    Columns: <key>, prior_class, current_class, status."""
    p = dict(zip(prior[key], prior[cls]))
    rows = []
    seen = set()
    for vk, cur in zip(current[key], current[cls]):
        seen.add(vk)
        pri = p.get(vk)
        if pri is None:
            status = "new"
        elif _rank(cur) > _rank(pri):
            status = "upgraded"
        elif _rank(cur) < _rank(pri):
            status = "downgraded"
        else:
            status = "unchanged"
        rows.append({key: vk, "prior_class": pri, "current_class": cur, "status": status})
    for vk, pri in p.items():
        if vk not in seen:
            rows.append({key: vk, "prior_class": pri, "current_class": None, "status": "dropped"})
    return pd.DataFrame(rows, columns=[key, "prior_class", "current_class", "status"])


def newly_reportable(current: pd.DataFrame, prior: pd.DataFrame, **kw) -> pd.DataFrame:
    """The variants that make a re-analysis worthwhile: NEW or UPGRADED and now P/LP."""
    d = diff(current, prior, **kw)
    return d[d["status"].isin(("new", "upgraded")) & d["current_class"].isin(REPORTABLE)].reset_index(drop=True)


def case_is_newly_reportable(current: pd.DataFrame, prior: pd.DataFrame, **kw) -> bool:
    """True if re-review is warranted for this case (something new/upgraded reached reportable)."""
    return len(newly_reportable(current, prior, **kw)) > 0
