"""Batch: run EVERY proband end-to-end -> a top-10 CSV per case (the challenge deliverable), plus a run summary.

Shared setup (ClinGen gene sets, pedigree, one DuckDB connection, the idempotent Monarch ATTACH) is done ONCE and
reused across all cases via scripts.run_proband.prioritize_case — no per-case re-setup. Writes to an OUTPUT dir that
is NOT the student-submission dir (never touches /root/bioconnect/*_top10.csv). No truth/judge material is read:
the pipeline sees only the student dataset + public reference DBs.

Run:  PYTHONPATH=. python3 scripts/run_all_probands.py [--out-dir .cache/top10] [--limit N]
"""
from __future__ import annotations
import argparse
import os
import time
import duckdb
import pandas as pd
from scripts.run_proband import load_gene_sets, prioritize_case, load_known_genes

PREP = "/root/bioconnect/prepared"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prepared", default=PREP)
    ap.add_argument("--cache", default=".cache")
    ap.add_argument("--snv-class", default=".cache/cohort_annotations_classified.parquet")
    ap.add_argument("--cnv-store", default=".cache/cnv_store.parquet")
    ap.add_argument("--out-dir", default=".cache/top10", help="OUR outputs — NOT the student submission dir")
    ap.add_argument("--limit", type=int, default=0, help="first N cases (0 = all)")
    ap.add_argument("--cases", default="", help="comma-separated case_ids to run (default: all)")
    ap.add_argument("--agent", action="store_true",
                    help="AGENT-AUGMENTED arm: drive the CaseSession with the LLM (refine HPO + gated PP4) per case")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    con = duckdb.connect()

    pedp = f"{a.prepared}/pedigree_corrected.parquet"
    ped = con.execute(f"SELECT * FROM read_parquet('{pedp}')").df() if os.path.exists(pedp) else pd.DataFrame()
    gene_moi, hi_genes, ts_genes = load_gene_sets(con, a.cache)
    known = load_known_genes(con, a.snv_class)

    cases = [r[0] for r in con.execute(
        f"SELECT DISTINCT case_id FROM read_parquet('{a.prepared}/snv.parquet') ORDER BY case_id").fetchall()]
    if a.cases:
        want = set(a.cases.split(","))
        cases = [c for c in cases if c in want]
    if a.limit:
        cases = cases[:a.limit]

    def _clinical_text(cid):
        src = "/root/bioconnect/dataset_clinical_curated.deid.parquet"
        src = src if os.path.exists(src) else "/root/bioconnect/dataset.parquet"
        row = con.execute(f"SELECT DISTINCT clinical_indication_text FROM read_parquet('{src}') WHERE student_case_id = ?",
                          [cid]).fetchone()
        return row[0] if row and row[0] else ""

    rows, t0 = [], time.time()
    for i, cid in enumerate(cases, 1):
        try:
            top, session, meta = prioritize_case(con, cid, prepared=a.prepared, snv_class=a.snv_class,
                                                 cnv_store=a.cnv_store, cache=a.cache, ped=ped,
                                                 gene_moi=gene_moi, hi_genes=hi_genes, ts_genes=ts_genes,
                                                 known_genes=known)
            if a.agent:                       # AGENT arm: LLM refines HPO + proposes gated PP4, deterministically re-folds
                from acmg.session_agent import agent_review
                top = agent_review(session, _clinical_text(cid), cr_index=f"{a.cache}/hp.index")
            top.head(10).to_csv(f"{a.out_dir}/{cid}_top10.csv", index=False)
            t1 = top.iloc[0] if len(top) else None
            rows.append({"case_id": cid, **meta, "ok": True,
                         "top1_kind": None if t1 is None else t1["kind"], "top1_gene": None if t1 is None else t1["gene"],
                         "top1_class": None if t1 is None else t1["acmg_class"]})
            print(f"[{i}/{len(cases)}] {cid} {meta['structure']} snv={meta['n_snv']} cnv={meta['n_cnv']} "
                  f"hpo={meta['n_hpo']} -> top1={None if t1 is None else t1['gene']} ({None if t1 is None else t1['kind']})")
        except Exception as e:                                   # one bad case must not sink the batch
            rows.append({"case_id": cid, "ok": False, "error": f"{type(e).__name__}: {e}"})
            print(f"[{i}/{len(cases)}] {cid} FAILED: {type(e).__name__}: {e}")
    pd.DataFrame(rows).to_csv(f"{a.out_dir}/_summary.csv", index=False)
    ok = sum(1 for r in rows if r.get("ok"))
    print(f"\n{ok}/{len(cases)} cases -> {a.out_dir}/  ({time.time() - t0:.0f}s); summary -> {a.out_dir}/_summary.csv")


if __name__ == "__main__":
    main()
