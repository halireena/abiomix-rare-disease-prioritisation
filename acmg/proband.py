"""run_proband — the per-proband command: fold the phenotype-aware layer over SNV *and* CNV candidates and
return the ranked top-N (+ an open REPL session). Pure composition of shipped pieces (kernel.classify -> case.
classify_case + candidate.candidate_table + rank.rank_fn + repl.CaseSession); no new abstraction.

The caller supplies the proband's candidates already joined to the datalakes:
  snv_free : the proband's SNVs joined to the phenotype-FREE datalake — variant_key, gene, total_points, criteria,
             acmg_class (from acmg.kernel.classify over the cohort).
  cnv      : the proband's classified CNVs — cnv_key, gene, cnv_class, svtype, cn (from acmg.cnv / annotate_cnv).
plus the case metadata (zygosity, gene_moi, HI/TS gene sets, case_structure, family facts) derived from the
pedigree + gene_curation. A scripts/run_proband.py wrapper does the datalake I/O; this stays testable.
"""
from __future__ import annotations
import duckdb
from .repl import CaseSession
from .rank import rank_fn as _make_rank_fn


def _has_table(con, name: str) -> bool:
    return con.execute("SELECT count(*) FROM information_schema.tables WHERE table_name = ?", [name]).fetchone()[0] > 0


def run_proband(con: duckdb.DuckDBPyConnection, case_id: str, snv_free, *, cnv=None, hpo=(), zygosity=None,
                gene_moi=None, hi_genes=None, ts_genes=None, case_structure: str = "singleton",
                pedigree_facts=None, rank_table: str = "gene_phenotype", named_genes=None,
                store: "str | None" = ".cache/case_sessions.jsonl", top_n: int = 10):
    """Phenotype-aware, CNV-inclusive per-proband run. Returns (top_df, session) — `session` is the open
    acmg.repl.CaseSession for interactive refinement (refine HPO -> re-rank, propose gated facts -> re-fold).
    A pathogenic CNV competes with the SNVs in the same ranked table. IC-weighted phenotype ranking reads the
    `gene_phenotype` table — the caller populates it from the Monarch KG (rank.monarch_gene_phenotype); tests
    inject a fixture table. Absent (nothing populated it), phenotype_score stays 0 and rank is class-only."""
    rf = _make_rank_fn(con, rank_table) if _has_table(con, rank_table) else None
    session = CaseSession(con, case_id, snv_free, cnv=cnv, hi_genes=hi_genes, ts_genes=ts_genes, store=store,
                          case_structure=case_structure, zygosity=zygosity, gene_moi=gene_moi,
                          pedigree_facts=pedigree_facts or {}, hpo=hpo, named_genes=named_genes,
                          **({"rank_fn": rf} if rf else {}))
    return session.top(top_n), session
