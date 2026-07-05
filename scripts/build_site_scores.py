"""Cohort site scores from PRUNED REMOTE tabix joins — no multi-GB downloads.

The cohort's distinct loci need three site-level scores the local tables don't carry:
gnomAD population AF (PM2/BA1/BS1), a missense in-silico score (PP3/BP4), and a splice
score (splice PP3/BP7). The naive path downloads 60–120 GB reference files. This does NOT:
it fetches ONLY the bytes covering the cohort's loci, over HTTP range reads against
bgzip+tabix files, via duckhts `read_tabix(url, region := ..., index_path := <local .tbi>)`.

Mechanism (proven): cache the small `.tbi` locally (KB–MB), then for each cohort region issue
`read_tabix` with a `region :=` — duckhts reads the tabix index, seeks the matching bgzip blocks,
and returns ONLY that region's records. The fat VCF INFO (~23 KB/record for gnomAD) is PRUNED in
SQL inside the read query (`regexp_extract` the AF fields) so only the needed columns cross Python.

To minimise the number of tabix seeks we MERGE adjacent loci (within `--gap` bp, default 10 kb) into
one region, and issue the region reads CONCURRENTLY (latency-bound; a thread pool with a per-thread
DuckDB connection turns ~1.2 s/region serial into ~0.05 s/region effective).

Sources (range support verified with `curl -sI` before use — see docs/scale.md):
  - gnomAD v2.1.1 exomes (GRCh37) on GCS: AF, AF_popmax          -> filtering_af, af_popmax  [REMOTE]
  - CADD v1.6 GRCh37 whole-genome SNVs (Kircher lab): PHRED       -> cadd_phred               [REMOTE]
  - REVEL: no host serves a range+tabix copy (only a 667 MB zip) -> revel stays NULL; see FALLBACK
  - SpliceAI precomputed: auth-gated (Illumina BaseSpace)        -> spliceai stays NULL; REST residual

`ref`/`end` are DuckDB reserved words, so tabix columns are aliased to non-reserved names.

Run:
    PYTHONPATH=. python3 scripts/build_site_scores.py \
        --parquet /root/bioconnect/unique_snvs_for_annotation_v3.parquet --cache .cache \
        --out .cache/site_scores.parquet [--chrom 22] [--gap 10000] [--workers 32]
"""
from __future__ import annotations
import argparse
import threading
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import duckdb
import pandas as pd

# --------------------------------------------------------------------------------------------------
# Remote sources. `verify_url` is checked with an HTTP HEAD (accept-ranges: bytes) before any read.
# `prune_sql` selects ONLY the needed columns from a `read_tabix(...)` region read — {call} is the
# read_tabix(...) call, aliased `rt`. Output columns MUST be: chrom, pos, ref_, alt_, <score cols>.
# --------------------------------------------------------------------------------------------------
GNOMAD_URL = ("https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/"
              "exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz")
CADD_URL = "https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz"

SOURCES = {
    "gnomad": {
        "data_url": GNOMAD_URL,
        "tbi_url": GNOMAD_URL + ".tbi",
        "tbi_name": "gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi",
        # VCF cols: 0 chrom, 1 pos, 3 ref, 4 alt, 7 INFO. AF/AF_popmax pruned out of the fat INFO.
        "prune_sql": """
            SELECT rt.column0 AS chrom, rt.column1 AS pos, rt.column3 AS ref_, rt.column4 AS alt_,
                   TRY_CAST(regexp_extract(rt.column7, '(^|;)AF=([0-9.eE+-]+)', 2)        AS DOUBLE) AS filtering_af,
                   TRY_CAST(regexp_extract(rt.column7, '(^|;)AF_popmax=([0-9.eE+-]+)', 2) AS DOUBLE) AS af_popmax
            FROM {call} AS rt
        """,
        "score_cols": ["filtering_af", "af_popmax"],
    },
    "cadd": {
        "data_url": CADD_URL,
        "tbi_url": CADD_URL + ".tbi",
        "tbi_name": "cadd_v1.6_GRCh37_whole_genome_SNVs.tsv.gz.tbi",
        # CADD cols: 0 chrom, 1 pos, 2 ref, 3 alt, 4 RawScore, 5 PHRED. Allele-specific (ref+alt match).
        "prune_sql": """
            SELECT rt.column0 AS chrom, rt.column1 AS pos, rt.column2 AS ref_, rt.column3 AS alt_,
                   TRY_CAST(rt.column5 AS DOUBLE) AS cadd_phred
            FROM {call} AS rt
        """,
        "score_cols": ["cadd_phred"],
    },
}

# Output schema. revel/spliceai are carried (NULL) so a later local-extract fills them without a
# schema change; the kernel abstains on a NULL (fail-safe), never guesses.
OUT_COLS = ["chrom", "pos", "ref", "alt", "filtering_af", "af_popmax", "revel", "cadd_phred", "spliceai"]


class _T:
    """Wall-clock timer that prints and accumulates."""

    def __init__(self):
        self.rows: list[tuple[str, float, str]] = []

    def step(self, label: str, fn):
        t = time.time()
        note = fn()
        dt = time.time() - t
        self.rows.append((label, dt, note or ""))
        print(f"  [{dt:7.2f}s] {label}{('  ' + note) if note else ''}")
        return note

    def total(self) -> float:
        return sum(dt for _, dt, _ in self.rows)


def http_accepts_ranges(url: str) -> bool:
    """HEAD the URL and confirm `accept-ranges: bytes` (the precondition for a tabix region read)."""
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return (r.headers.get("Accept-Ranges", "").lower() == "bytes")
    except Exception:
        return False


def cache_tbi(src: dict, cache: Path) -> str:
    """Cache the source's `.tbi` locally (KB–MB), once. Returns the local path."""
    dst = cache / src["tbi_name"]
    if not dst.exists():
        cache.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(src["tbi_url"], dst)
    return str(dst)


# --------------------------------------------------------------------------------------------------
# 1. cohort distinct loci  ->  2. merged regions
# --------------------------------------------------------------------------------------------------
def build_cohort(con: duckdb.DuckDBPyConnection, parquet: str, chrom: str | None) -> int:
    """Distinct cohort loci (chrom,pos,ref,alt). Optionally restrict to one chromosome (benchmark)."""
    where = f"WHERE replace(CAST(CHROM AS VARCHAR),'chr','') = '{chrom}'" if chrom else ""
    con.execute(f"""
        CREATE OR REPLACE TABLE cohort AS
        SELECT DISTINCT
            replace(CAST(CHROM AS VARCHAR), 'chr', '') AS chrom,
            CAST(POS AS BIGINT)                        AS pos,
            REF                                        AS ref,
            t.alt_x                                    AS alt
        FROM read_parquet('{parquet}'), UNNEST(ALT) AS t(alt_x)
        {where}
    """)
    return con.execute("SELECT count(*) FROM cohort").fetchone()[0]


def merge_regions(con: duckdb.DuckDBPyConnection, gap: int) -> pd.DataFrame:
    """Bin adjacent loci within `gap` bp into one region per chromosome. Fewer regions = fewer seeks.
    Standard gaps-and-islands: a new island starts when the gap to the previous position exceeds `gap`."""
    return con.execute(f"""
        WITH pts AS (SELECT DISTINCT chrom, pos FROM cohort),
        o AS (SELECT chrom, pos, LAG(pos) OVER (PARTITION BY chrom ORDER BY pos) AS prev FROM pts),
        f AS (SELECT chrom, pos, CASE WHEN prev IS NULL OR pos - prev > {gap} THEN 1 ELSE 0 END AS ng FROM o),
        g AS (SELECT chrom, pos, SUM(ng) OVER (PARTITION BY chrom ORDER BY pos) AS gid FROM f)
        SELECT chrom, min(pos) AS start, max(pos) AS end_
        FROM g GROUP BY chrom, gid ORDER BY chrom, start
    """).df()


# --------------------------------------------------------------------------------------------------
# 3. pruned remote region reads (concurrent)
# --------------------------------------------------------------------------------------------------
def fetch_source(src: dict, regions: pd.DataFrame, tbi: str, workers: int, retries: int = 2) -> pd.DataFrame:
    """Read every merged region from the remote bgzip+tabix file, PRUNED to the needed columns, and
    concatenate. Region reads run concurrently on a thread pool (I/O-bound); each thread owns its own
    DuckDB connection (DuckDB cursors from one connection serialise, so we open one per thread). A
    failed region is retried a few times then returned empty (error-as-value; one bad seek != abort)."""
    tl = threading.local()
    call = f"read_tabix('{src['data_url']}', region := ?, index_path := '{tbi}')"
    sql = src["prune_sql"].format(call=call)

    def read_region(reg: str) -> pd.DataFrame:
        if not hasattr(tl, "con"):
            tl.con = duckdb.connect()
            tl.con.execute("LOAD duckhts")
        for attempt in range(retries + 1):
            try:
                return tl.con.execute(sql, [reg]).df()
            except Exception:
                if attempt == retries:
                    return pd.DataFrame(columns=["chrom", "pos", "ref_", "alt_"] + src["score_cols"])
                time.sleep(0.5 * (attempt + 1))

    reg_strs = [f"{r.chrom}:{int(r.start)}-{int(r.end_)}" for r in regions.itertuples()]
    with ThreadPoolExecutor(max_workers=workers) as ex:
        parts = list(ex.map(read_region, reg_strs))
    parts = [p for p in parts if len(p)]
    if not parts:
        return pd.DataFrame(columns=["chrom", "pos", "ref_", "alt_"] + src["score_cols"])
    return pd.concat(parts, ignore_index=True)


def join_source(con: duckdb.DuckDBPyConnection, name: str, pruned: pd.DataFrame, score_cols: list[str]) -> None:
    """Left-join a source's pruned records onto the exact cohort loci (chrom,pos,ref,alt). The region
    read is a superset (it returns every record in the window); the join keeps only exact-locus matches."""
    con.register(f"_{name}_raw", pruned)
    cols = ", ".join(f"s.{c}" for c in score_cols)
    con.execute(f"""
        CREATE OR REPLACE TABLE {name}_scores AS
        SELECT c.chrom, c.pos, c.ref, c.alt, {cols}
        FROM cohort c
        LEFT JOIN (SELECT DISTINCT chrom, CAST(pos AS BIGINT) AS pos, ref_, alt_, {', '.join(score_cols)}
                   FROM _{name}_raw) s
          ON s.chrom = c.chrom AND s.pos = c.pos AND s.ref_ = c.ref AND s.alt_ = c.alt
    """)
    con.unregister(f"_{name}_raw")


# --------------------------------------------------------------------------------------------------
def build(parquet: str, cache: str, out: str, *, chrom: str | None = None,
          gap: int = 10_000, workers: int = 32, verify: bool = True, revel: str | None = None) -> dict:
    """Full build. Returns a small stats dict (also used by the benchmark)."""
    con = duckdb.connect()
    con.execute("LOAD duckhts")
    t = _T()
    cache_p = Path(cache)

    print("== cohort ==")
    n_loci = int(t.step("build cohort (distinct loci)",
                        lambda: f"loci={build_cohort(con, parquet, chrom)}").split("=")[1])
    regions = merge_regions(con, gap)
    print(f"  merged {n_loci} loci -> {len(regions)} regions (gap={gap} bp)")

    active: list[str] = []
    for name, src in SOURCES.items():
        print(f"== {name} (remote tabix) ==")
        if verify and not http_accepts_ranges(src["data_url"]):
            print(f"  SKIP {name}: no accept-ranges: bytes at {src['data_url']}")
            continue
        tbi = cache_tbi(src, cache_p)
        t0 = time.time()
        pruned = fetch_source(src, regions, tbi, workers)
        dt = time.time() - t0
        t.rows.append((f"read+prune {name}", dt, f"records={len(pruned)}"))
        print(f"  [{dt:7.2f}s] read+prune {name} @{len(regions)} regions (W={workers})  records={len(pruned)}")
        join_source(con, name, pruned, src["score_cols"])
        m = con.execute(f"SELECT count(*) FROM {name}_scores WHERE {src['score_cols'][0]} IS NOT NULL").fetchone()[0]
        print(f"  matched {m}/{n_loci} loci")
        active.append(name)

    # Local REVEL: no host serves a range+tabix REVEL copy, so it's a one-time local Parquet (77.9M loci,
    # ~172 MB, sorted). The `WHERE chrom IN (...)` forces row-group min/max pushdown so only the needed
    # chromosomes' row groups are scanned. Fills PP3/BP4, which otherwise abstains at cohort scale.
    revel_pq = revel or str(cache_p / "revel_grch37.parquet")
    have_revel = Path(revel_pq).exists()
    if have_revel:
        con.execute(f"""CREATE OR REPLACE TABLE revel_scores AS
            SELECT c.chrom, c.pos, c.ref, c.alt, r.revel
            FROM cohort c JOIN read_parquet('{revel_pq}') r
              ON r.chrom = c.chrom AND r.pos = c.pos AND r.ref = c.ref AND r.alt = c.alt
            WHERE r.chrom IN (SELECT DISTINCT chrom FROM cohort)""")
        m = con.execute("SELECT count(*) FROM revel_scores").fetchone()[0]
        print(f"== revel (local parquet, chrom-pushdown) ==  matched {m}/{n_loci}")

    # Assemble the compact site_scores table. Remote sources + local REVEL are LEFT JOINed; spliceai has no
    # source here (NULL -> the kernel abstains, fail-safe).
    print("== assemble site_scores ==")
    have = {c for name in active for c in SOURCES[name]["score_cols"]}
    joins = [f"LEFT JOIN {name}_scores USING (chrom, pos, ref, alt)" for name in active]
    selects = [col for name in active for col in SOURCES[name]["score_cols"]]
    if have_revel:
        joins.append("LEFT JOIN revel_scores USING (chrom, pos, ref, alt)")
        selects.append("revel"); have.add("revel")
    join_sql = " ".join(joins)
    score_select = ", ".join(selects)
    extra_nulls = ", ".join(f"CAST(NULL AS DOUBLE) AS {c}" for c in ["revel", "spliceai"] if c not in have)
    pieces = ["c.chrom", "c.pos", "c.ref", "c.alt"]
    if score_select:
        pieces.append(score_select)
    if extra_nulls:
        pieces.append(extra_nulls)
    con.execute(f"""
        CREATE OR REPLACE TABLE site_scores AS
        SELECT {', '.join(pieces)}
        FROM cohort c {join_sql}
    """)
    # Reorder to the canonical schema and write ORDER BY chrom,pos (interval-friendly row groups).
    have_cols = [r[0] for r in con.execute("DESCRIBE site_scores").fetchall()]
    ordered = [c for c in OUT_COLS if c in have_cols]
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    con.execute(f"""
        COPY (SELECT {', '.join(ordered)} FROM site_scores
              ORDER BY CAST(regexp_replace(chrom,'[^0-9]','','g') AS INTEGER) NULLS LAST, chrom, pos)
        TO '{out}' (FORMAT parquet)
    """)
    stats = {
        "loci": n_loci, "regions": len(regions), "gap": gap, "workers": workers,
        "sources_active": active, "wall": t.total(), "out": out,
        "af_matched": con.execute("SELECT count(*) FROM site_scores WHERE filtering_af IS NOT NULL").fetchone()[0]
                      if "gnomad" in active else 0,
        "cadd_matched": con.execute("SELECT count(*) FROM site_scores WHERE cadd_phred IS NOT NULL").fetchone()[0]
                        if "cadd" in active else 0,
    }
    print(f"\nwrote {out}  ({n_loci} loci, {t.total():.1f}s wall, sources={active})")
    print(f"  schema: {ordered}")
    return stats


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--parquet", default="/root/bioconnect/unique_snvs_for_annotation_v3.parquet")
    ap.add_argument("--cache", default=".cache")
    ap.add_argument("--out", default=".cache/site_scores.parquet")
    ap.add_argument("--chrom", default=None, help="restrict to one chromosome (benchmark subset)")
    ap.add_argument("--gap", type=int, default=10_000, help="merge loci within this many bp into one region")
    ap.add_argument("--workers", type=int, default=32, help="concurrent region reads")
    ap.add_argument("--no-verify", action="store_true", help="skip the accept-ranges HEAD check")
    ap.add_argument("--revel", default=None, help="local REVEL parquet (default: <cache>/revel_grch37.parquet)")
    args = ap.parse_args()
    build(args.parquet, args.cache, args.out, chrom=args.chrom, gap=args.gap,
          workers=args.workers, verify=not args.no_verify, revel=args.revel)


if __name__ == "__main__":
    main()
