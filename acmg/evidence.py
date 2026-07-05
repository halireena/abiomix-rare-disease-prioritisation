"""Evidence sources for the shortlist / literature arm — Europe PMC alone is not enough.

The agent gathers VARIANT- and GENE-level evidence from several GRCh37-friendly sources, then PROPOSES an
ACMG criterion over it (recorded + gated via acmg.agent). Sources:

  - Ensembl Variation (GRCh37): clinical significance + phenotypes + citations per rsid.
  - LitVar2 (NCBI): the variant -> publication index that Europe PMC lacks (finds papers that mention the
    exact variant, by rsid / HGVS).
  - MARRVEL (api.marrvel.org, hg19, no key): aggregated OMIM / gnomAD / dbNSFP / model-organism (DIOPT)
    evidence for a gene/variant.
  - Europe PMC / NCBI eutils (acmg.literature): reading the actual abstracts/full text of the linked PMIDs.
  - DECIPHER: PATIENT CASE data. For an SNV shortlist candidate, other patients with variants in the same
    gene + their phenotypes are case-level PP4 / PS4 / matchmaking evidence. Consent-gated: needs a
    registered DECIPHER API token; the population/aggregate tier is downloadable for the CNV track.
  - Monarch (attached in acmg.rank): cross-species / model-organism gene-phenotype.

Everything returns plain dicts/lists so the agent can fold them into one context; nothing here decides.
"""
from __future__ import annotations
import requests
from . import literature

ENSEMBL_GRCH37 = "https://grch37.rest.ensembl.org"
MARRVEL = "http://api.marrvel.org/data"
LITVAR2 = "https://www.ncbi.nlm.nih.gov/research/litvar2-api"
DECIPHER = "https://www.deciphergenomics.org/api"


def ensembl_variation(rsid: str, *, server: str = ENSEMBL_GRCH37, timeout: int = 30) -> dict:
    """Clinical significance + phenotypes (ClinVar traits, HP terms) for a variant by rsid (GRCh37)."""
    r = requests.get(f"{server}/variation/human/{rsid}", params={"content-type": "application/json", "phenotypes": 1}, timeout=timeout)
    r.raise_for_status()
    d = r.json()
    return {"clinical_significance": d.get("clinical_significance", []),
            "phenotypes": [{"trait": p.get("trait"), "source": p.get("source"), "genes": p.get("genes")}
                           for p in d.get("phenotypes", [])]}


def litvar2(query: str, *, timeout: int = 30) -> list[dict]:
    """Variant -> literature index: match a variant (rsid / 'GENE p.Xxx') to its publication footprint."""
    r = requests.get(f"{LITVAR2}/variant/autocomplete/", params={"query": query}, timeout=timeout)
    r.raise_for_status()
    return [{"rsid": x.get("rsid"), "gene": x.get("gene"), "hgvs": x.get("hgvs"),
             "pmids_count": x.get("pmids_count"), "clinical_significance": x.get("data_clinical_significance")}
            for x in r.json()][:5]


def marrvel_dbnsfp(variant_key: str, *, timeout: int = 30) -> dict:
    """MARRVEL/dbNSFP functional prediction scores for a variant ('chr-pos-ref-alt', hg19)."""
    r = requests.get(f"{MARRVEL}/dbnsfp/variant/{variant_key}", timeout=timeout)
    r.raise_for_status()
    return r.json()


def decipher_gene(gene: str, *, token: str | None = None, timeout: int = 30) -> dict:
    """DECIPHER patient CASE data for a gene (other patients' variants + phenotypes) — SNV-relevant
    matchmaking evidence. Requires a registered DECIPHER API token (consent-gated). Returns {} without one."""
    if not token:
        return {"note": "DECIPHER case data requires a registered API token (consent-gated); skipped."}
    r = requests.get(f"{DECIPHER}/genes/{gene}/variants", headers={"authorization": f"Bearer {token}"}, timeout=timeout)
    r.raise_for_status()
    return r.json()


def gather(*, variant_key: str, gene: str, rsid: str | None = None, hpo_labels: list[str] | None = None,
           decipher_token: str | None = None) -> dict:
    """Assemble a multi-source evidence context for the agent. Best-effort: a source that errors is skipped."""
    ev: dict = {"variant_key": variant_key, "gene": gene}
    def _try(name, fn):
        try:
            ev[name] = fn()
        except Exception as e:
            ev[name] = {"error": str(e)[:120]}
    if rsid:
        _try("ensembl_variation", lambda: ensembl_variation(rsid))
        _try("litvar2", lambda: litvar2(rsid))
    _try("marrvel_dbnsfp", lambda: marrvel_dbnsfp(variant_key))
    _try("decipher_cases", lambda: decipher_gene(gene, token=decipher_token))
    _try("literature", lambda: literature.evidence_for(gene, hpo_labels, n=4))
    return ev
