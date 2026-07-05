"""Produce the OpenSpliceAI splice scores for the annotation datalake — the compute behind the `spliceai` source.

On-device SpliceAI-equivalent scoring (no auth-gated Illumina precompute) for the cohort's SPLICE-RELEVANT loci
only (consequence in splice_region/donor/acceptor/intron/synonymous — where the kernel's splice PP3/BP7 applies;
missense already gets REVEL, LoF gets PVS1). Reads the loci from the annotation datalake (the VEP source gives
consequence), scores ONLY the SpliceAI increment (loci without a spliceai score yet), and writes
.cache/spliceai_openspliceai.parquet — which scripts/annotate_cohort.py then merges as the `spliceai` source.

Parallelism: N SHARDS x T torch-threads ~= nproc (torch conv saturates ~4 threads, so 10x2 or 5x4 fill 20 cores;
a single process at threads=2 would be ~50h). Needs the [splice] extra (bash scripts/setup_openspliceai.sh). GRCh37.

Run:  PYTHONPATH=. python3 scripts/run_openspliceai.py --shards 10 --threads 2
"""
from __future__ import annotations
import argparse
import os
import concurrent.futures as cf
from pathlib import Path
import duckdb
import pandas as pd

SPLICE_CONSEQUENCES = ("splice_region", "splice_donor", "splice_acceptor", "intron", "synonymous")
SPLICEAI_RELEASE = "openspliceai"   # matches SOURCE_RELEASE['spliceai'] in annotate_cohort


def splice_relevant_increment(con, store: str, out: str, *, all_loci: bool = False) -> pd.DataFrame:
    """Loci NOT yet scored in `out` -> the SpliceAI increment. Default = splice-relevant only (by VEP consequence,
    where the kernel's splice PP3/BP7 applies); `all_loci=True` scores every locus in the datalake."""
    have = f"read_parquet('{out}')" if Path(out).exists() else None
    where_new = f"AND s.variant_key NOT IN (SELECT variant_key FROM {have})" if have else ""
    scope = "" if all_loci else f"AND s.consequence IN {SPLICE_CONSEQUENCES}"
    return con.execute(f"""
        SELECT s.variant_key,
               split_part(s.variant_key,'-',1) AS chrom, CAST(split_part(s.variant_key,'-',2) AS BIGINT) AS pos,
               split_part(s.variant_key,'-',3) AS ref, split_part(s.variant_key,'-',4) AS alt
        FROM read_parquet('{store}') s
        WHERE 1=1 {scope} {where_new}
    """).df()


def _score_shard(rows, fasta, model_dir, threads, chunk_size, flanking_size):
    from acmg.annotate import spliceai_local
    variants = [(r[0], int(r[1]), r[2], r[3]) for r in rows]   # (chrom,pos,ref,alt); no Python data loop downstream
    return spliceai_local(variants, fasta=fasta, model_dir=model_dir, threads=threads, chunk_size=chunk_size,
                          flanking_size=flanking_size)   # num_models handled once in main (avoid per-shard symlink race)


def _subset_models(model_dir: str, num_models: int) -> str:
    """Build the first-N-model subdir ONCE (before sharding) so 10 shards don't race on the same symlink dir."""
    import glob
    sub = os.path.join(model_dir, f"_first{num_models}")
    os.makedirs(sub, exist_ok=True)
    for f in sorted(glob.glob(os.path.join(model_dir, "*.pt")))[:num_models]:
        dst = os.path.join(sub, os.path.basename(f))
        if not os.path.lexists(dst):
            os.symlink(os.path.abspath(f), dst)
    return sub


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--store", default=".cache/anno_store.parquet", help="datalake with the VEP consequence source")
    ap.add_argument("--out", default=".cache/spliceai_openspliceai.parquet")
    ap.add_argument("--fasta", default="/root/GRCh37/human_g1k_v37.fasta")
    ap.add_argument("--model-dir", default=None, help="default: .cache/openspliceai_models/mane_<flanking>nt (must match --flanking)")
    ap.add_argument("--all", action="store_true", help="score ALL loci in the datalake, not just splice-relevant consequences")
    ap.add_argument("--shards", type=int, default=10, help="concurrent worker processes (N x --threads ~= nproc)")
    ap.add_argument("--threads", type=int, default=2, help="torch intra-op threads PER shard (conv saturates ~4)")
    ap.add_argument("--chunk-size", type=int, default=1500, help="per-call variant chunk (bounds memory; un-chunked OOMs)")
    ap.add_argument("--flanking", type=int, default=10000,
                    help="context nt (+-flanking/2). 10000=ClinGen/Walker-calibrated & matches the 10000nt models; "
                         "smaller = faster but feeds the 10000nt-trained model less context (speed/accuracy trade — measure)")
    ap.add_argument("--num-models", type=int, default=None, help="use only the first N of the 5-model ensemble (speed knob)")
    ap.add_argument("--limit", type=int, default=None, help="cap the increment (smoke test)")
    args = ap.parse_args()

    con = duckdb.connect()
    inc = splice_relevant_increment(con, args.store, args.out, all_loci=args.all)
    if args.limit:
        inc = inc.head(args.limit)
    model_base = args.model_dir or f".cache/openspliceai_models/mane_{args.flanking}nt"
    print(f"SpliceAI increment ({'ALL loci' if args.all else 'splice-relevant'}, unscored): {len(inc)} loci; "
          f"{args.shards} shards x {args.threads} threads; flanking={args.flanking} (+-{args.flanking//2}nt)"
          f"{' '+str(args.num_models)+'-model' if args.num_models else ' 5-model'}; models={model_base}")
    if not len(inc):
        print("increment empty — every splice-relevant locus already scored (full reuse)."); return

    model_dir = _subset_models(model_base, args.num_models) if args.num_models else model_base
    rows = list(inc[["chrom", "pos", "ref", "alt"]].itertuples(index=False, name=None))
    shards = [rows[i::args.shards] for i in range(args.shards)]      # round-robin split
    results = []
    with cf.ProcessPoolExecutor(max_workers=args.shards) as ex:
        futs = [ex.submit(_score_shard, sh, args.fasta, model_dir, args.threads, args.chunk_size, args.flanking)
                for sh in shards if sh]
        for i, f in enumerate(cf.as_completed(futs)):
            df = f.result(); results.append(df)
            print(f"  shard {i+1}/{len(futs)} done: {len(df)} scored")
    scored = pd.concat(results, ignore_index=True) if results else pd.DataFrame(columns=["chrom","pos","ref","alt","variant_key","spliceai"])

    # merge with any prior scores; write the source parquet the datalake consumes
    con.register("_new", scored)
    if Path(args.out).exists():
        con.execute(f"""COPY (
            SELECT * FROM read_parquet('{args.out}') WHERE variant_key NOT IN (SELECT variant_key FROM _new)
            UNION ALL BY NAME SELECT * FROM _new) TO '{args.out}.tmp' (FORMAT parquet)""")
        Path(f"{args.out}.tmp").replace(args.out)
    else:
        con.execute(f"COPY (SELECT * FROM _new) TO '{args.out}' (FORMAT parquet)")
    n = con.execute(f"SELECT count(*) FROM read_parquet('{args.out}')").fetchone()[0]
    nz = int((scored["spliceai"] > 0).sum()); ge2 = int((scored["spliceai"] >= 0.2).sum())
    skipped = int((scored["spliceai"] == 0).sum())   # zero = out-of-gene OR ref-mismatch (~3% of true SNVs; NOT 17%
    #                                                   — that earlier figure conflated indels naively compared to a
    #                                                   single ref base. True biallelic SNVs are ~96.8% ref-concordant.)
    print(f"scored {len(scored)} ({ge2} >=0.2 splice-PP3, {nz} >0, {skipped} zero incl ref-mismatch) -> {args.out} ({n} total)")
    print("  re-run scripts/annotate_cohort.py to merge the spliceai source into the datalake + reclassify.")


if __name__ == "__main__":
    main()
