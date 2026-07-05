"""Case REPL — a persistent, deterministic session for the PHENOTYPE-AWARE layer.

One proband's variant table lives in DuckDB (phenotype-free class x all annotations x genotype_fit), addressable
by SQL so the context lives OUTSIDE the prompt (no context rot; counting is a GROUP BY, not a token scan). A
session is a sequence of actions — query the table, refine HPO -> deterministic re-rank, propose a GATED case
fact -> re-fold — and EVERY action is persisted to an append-only JSONL ledger (SQL-queryable via read_json_auto).
So the session is auditable and REPLAYABLE (deterministic: the case + its action log reproduce the same state).

The LLM drives it (its tools ARE these methods); the SQL / fold (acmg.case.classify_case) / rank (acmg.rank) are
deterministic; the case facts are gated. It never edits the class directly — it proposes, and only an approved
digest counts. This is the RLM/SQL-REPL-over-context idea made concrete for one case.
"""
from __future__ import annotations
import hashlib
import json
import os
import time
import uuid
import pandas as pd
from .case import classify_case, _FAMILY_CRITERIA, _FAMILY_STRUCTURES
from .candidate import candidate_table, cnv_genotype_fit
from . import agent


def _default_rank(case_df: pd.DataFrame, hpo) -> list:
    """Deterministic placeholder phenotype rank: +2 for a variant whose gene is HPO-implicated. Real use injects
    rank_fn (acmg.rank / Monarch HPO->gene); kept trivial so the REPL is testable offline. `hpo` may be a set of
    HP terms (real) or {'genes': [...]} (test convenience)."""
    genes = set(hpo.get("genes", [])) if isinstance(hpo, dict) else set()
    return [2.0 if g in genes else 0.0 for g in case_df["gene"]]


_GENOTYPE_BOOST = {"biallelic_homozygous": 3, "biallelic_candidate": 2, "monoallelic_dominant": 1, "monoallelic_carrier": -1}


class CaseSession:
    """A persistent phenotype-aware session for one proband. See module docstring."""

    def __init__(self, con, case_id: str, phenotype_free: pd.DataFrame, *, cnv: "pd.DataFrame | None" = None,
                 hi_genes=None, ts_genes=None, store: "str | None" = ".cache/case_sessions.jsonl",
                 case_structure: str = "singleton", zygosity=None, gene_moi=None, pedigree_facts=None,
                 rank_fn=_default_rank, hpo=None, named_genes=None):
        self.con = con
        self.case_id = case_id
        self.session_id = uuid.uuid4().hex[:12]
        self.store = store
        self.free = phenotype_free.copy()
        self.cnv = cnv                                  # classified CNVs (cnv_key, gene, cnv_class[, svtype, cn])
        self.hi_genes = hi_genes
        self.ts_genes = ts_genes
        self.case_structure = case_structure
        self.zygosity = zygosity
        self.gene_moi = gene_moi
        self.pedigree_facts = pedigree_facts or {}      # {'de_novo':[...], 'comp_het':[...], 'segregating':[...]}
        self.rank_fn = rank_fn
        self.hpo = set(hpo or [])
        self.named_genes = set(named_genes or [])       # genes the clinician named in the referral (rank prior)
        self.facts: list[dict] = []                     # gated case-fact proposals (PP4/PS3/PS4)
        self._refresh()
        self._persist("open", {"case_structure": case_structure, "n_snv": int(len(self.free)),
                               "n_cnv": int(len(self.cnv)) if self.cnv is not None else 0})

    # --- deterministic recompute: SNV fold + CNV -> ONE unified candidate table (case_variants) ---
    def _refresh(self):
        approved = [{"variant_key": f["variant_key"], "criterion": f["criterion"], "points": f["points"]}
                    for f in self.facts if f["status"] == "approved"]
        snv = classify_case(self.free, case_structure=self.case_structure, zygosity=self.zygosity,
                            gene_moi=self.gene_moi, agent_facts=approved, **self.pedigree_facts)
        snv["phenotype_score"] = list(self.rank_fn(snv, self.hpo))
        cnv = None
        if self.cnv is not None and len(self.cnv):
            cnv = self.cnv.copy()
            cnv["genotype_fit"] = cnv_genotype_fit(cnv, hi_genes=self.hi_genes, ts_genes=self.ts_genes)
            cnv["phenotype_score"] = list(self.rank_fn(cnv, self.hpo))
        self.table = candidate_table(snv=snv, cnv=cnv, named_genes=self.named_genes)   # SNVs+CNVs ranked TOGETHER
        self.con.register("case_variants", self.table)

    def _state_digest(self) -> str:
        payload = {"hpo": sorted(self.hpo), "facts": sorted(f["digest"] for f in self.facts if f["status"] == "approved")}
        return "sha256:" + hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]

    def _persist(self, action: str, detail: dict):
        if not self.store:
            return
        os.makedirs(os.path.dirname(self.store) or ".", exist_ok=True)
        with open(self.store, "a") as fh:
            fh.write(json.dumps({"ts": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "session_id": self.session_id,
                                 "case_id": self.case_id, "action": action, "detail": detail,
                                 "state": self._state_digest()}, ensure_ascii=False) + "\n")

    # --- the REPL tools (LLM-callable; every one persisted) ---
    def sql(self, query: str) -> pd.DataFrame:
        """Read-only SQL over the `case_variants` table (the proband's context, addressable outside the prompt)."""
        df = self.con.execute(query).df()
        self._persist("sql", {"query": query, "rows": int(len(df))})
        return df

    def refine_hpo(self, add=(), remove=()) -> pd.DataFrame:
        """Add/remove HPO terms -> deterministic re-rank + re-fold. The whole reason it's a REPL."""
        self.hpo |= set(add)
        self.hpo -= set(remove)
        self._refresh()
        self._persist("refine_hpo", {"add": list(add), "remove": list(remove), "hpo": sorted(self.hpo)})
        return self.top()

    def propose_fact(self, variant_key: str, criterion: str, points: float, *, evidence=None) -> str:
        """Propose a GATED case fact (PP4/PS3/PS4) -> recorded pending; counts only after approve_fact (its exact
        digest). The model proposes; the human-approved digest is what moves the class."""
        if criterion in _FAMILY_CRITERIA and self.case_structure not in _FAMILY_STRUCTURES:
            raise ValueError(f"{criterion} is a pedigree criterion; a {self.case_structure} case cannot support it "
                             "(no parents/phasing). It is derived deterministically from the pedigree, not proposed.")
        payload = {"variant_key": variant_key, "criterion": criterion, "points": points, "evidence": evidence or {}}
        digest = agent._digest(payload)
        self.facts.append({**payload, "digest": digest, "status": "pending"})
        self._persist("propose_fact", {**{k: payload[k] for k in ("variant_key", "criterion", "points")}, "digest": digest})
        return digest

    def approve_fact(self, digest: str, approver: str) -> pd.DataFrame:
        for f in self.facts:
            if f["digest"] == digest:
                f["status"] = "approved"
        self._refresh()
        self._persist("approve_fact", {"digest": digest, "approver": approver})
        return self.top()

    def top(self, n: int = 10) -> pd.DataFrame:
        """Current top-N per the deterministic fold + rank (recomputed on every action)."""
        return self.table.head(n)
