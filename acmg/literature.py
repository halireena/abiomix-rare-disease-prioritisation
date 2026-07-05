"""Literature retrieval for the gated literature arm (PS3 / PP4 / PM3 / …).

Mirrors the paper's flow (Ma et al. 2025: retrieve, then an LLM summarises and proposes the criterion):
this is the retrieve step, from Europe PMC's free REST API (no key). The retrieved abstracts are handed
to acmg.agent.propose(..., llm=acmg.pi.as_llm()) so the model PROPOSES a criterion from real text; the
proposal is recorded with a digest and a human approves it before it can move the class (never automatic).

Europe PMC: https://europepmc.org/RestfulWebService
"""
from __future__ import annotations
import requests

EUROPEPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


def search(query: str, *, n: int = 5, timeout: int = 30) -> list[dict]:
    """Return up to `n` records ({pmid, title, year, abstract}) for a free-text query."""
    r = requests.get(EUROPEPMC, params={"query": query, "format": "json", "resultType": "core", "pageSize": n},
                     timeout=timeout)
    r.raise_for_status()
    out = []
    for x in r.json().get("resultList", {}).get("result", []):
        out.append({
            "pmid": x.get("pmid"),
            "title": x.get("title"),
            "year": x.get("pubYear"),
            "abstract": (x.get("abstractText") or "")[:1500],
        })
    return out


def evidence_for(gene: str, hpo_labels: list[str] | None = None, *, n: int = 5) -> list[dict]:
    """Retrieve literature for a candidate gene, phenotype-focused when HPO labels are given."""
    pheno = " OR ".join(f'"{p}"' for p in (hpo_labels or [])[:4])
    q = f'{gene} AND ({pheno})' if pheno else f'{gene} AND (variant OR mutation OR pathogenic)'
    return search(q, n=n)
