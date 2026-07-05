"""Cohort-scale annotation: annotate DISTINCT loci once, accumulate in an INCREMENTAL datalake, reuse.

Dedup the cohort's ~7.67M sample-calls to ~477k DISTINCT loci and annotate each ONCE. Intrinsic annotations
(variant-intrinsic + expensive) accumulate in a MULTI-SOURCE release datalake keyed by variant_key, where each
source writes its own columns tagged with its own source-release:

    VEP 116  -> gene/consequence/protein_pos/alt_aa1/transcript_id   (offline docker; the one tool that RUNS)
    REVEL v1.3 -> revel     OpenSpliceAI -> spliceai     CADD v1.7 -> cadd_phred   (static-file / on-device lookups)

The INCREMENT is PER SOURCE: only loci a given source has not yet annotated (or annotated under a stale release)
are (re)annotated; everything else is reused — across every cohort and run. The cheap, current-release reference
joins (ClinVar PS1/PM5, gnomAD mis_z, NMD escape, gnomAD AF) are recomputed each run and combined into the kernel
`annotations` shape (acmg.vep_map.REQUIRED_COLS) for acmg.kernel.classify. The SAME store serves CNVs (a parallel
datalake keyed by cnv_key).

REST is NOT a scale component: with the genome-wide offline cache the datalake leaves no "novel" set for REST to
fill. `rest_fallback` is the interactive single-variant / no-setup path + the rsID evidence arm (see its docstring).

Run (dockerized VEP 116 is the default --mode; each source annotates only its increment):
    PYTHONPATH=. python3 scripts/annotate_cohort.py --parquet /root/bioconnect/prepared/snv.parquet --cache .cache --classify
See docs/scale.md for measured timings.
"""
from __future__ import annotations
import argparse
import shutil
import time
from pathlib import Path

import duckdb
import pandas as pd

from acmg.clinvar import load_clinvar, ps1_pm5
from acmg.constraint import load_constraint, add_mis_z
from acmg.nmd import load_exons, nmd_escaping
from acmg.vep_map import REQUIRED_COLS, most_severe_kernel_consequence


class _T:
    """Tiny wall-clock timer that prints and accumulates."""

    def __init__(self):
        self.rows: list[tuple[str, float, str]] = []

    def step(self, label: str, fn):
        t = time.time()
        note = fn()
        dt = time.time() - t
        self.rows.append((label, dt, note or ""))
        print(f"  [{dt:6.2f}s] {label}{('  ' + note) if note else ''}")
        return dt

    def total(self) -> float:
        return sum(dt for _, dt, _ in self.rows)


# --------------------------------------------------------------------------------------------------
# 1. cohort: distinct loci (annotate once)
# --------------------------------------------------------------------------------------------------
def build_cohort(con: duckdb.DuckDBPyConnection, parquet: str) -> int:
    """Distinct loci to annotate ONCE -> `cohort(variant_key,chrom,pos,ref,alt)`. Derive from OUR OWN
    de-identified canonical bundle (acmg.ingest schema: variant_key/chrom/pos/ref/alt), NOT an external
    'unique variants' file of unknown provenance/filtering. Over the ~7.67M sample-calls there are ~477k
    DISTINCT SNV loci; dedup is what makes annotation cheap — annotate each locus once, reuse across every
    sample and release. (A legacy CHROM/POS/REF/ALT[list] source is still accepted.)"""
    cols = [c[0] for c in con.execute(f"SELECT * FROM read_parquet('{parquet}') LIMIT 0").description]
    if "variant_key" in cols:  # our canonical bundle
        kind = "AND (variant_kind = 'snv' OR variant_kind IS NULL)" if "variant_kind" in cols else ""
        con.execute(f"""
            CREATE OR REPLACE TABLE cohort AS
            SELECT DISTINCT variant_key,
                   replace(CAST(chrom AS VARCHAR), 'chr', '') AS chrom, CAST(pos AS BIGINT) AS pos, ref, alt
            FROM read_parquet('{parquet}') WHERE ref IS NOT NULL AND alt IS NOT NULL {kind}
        """)
    else:  # legacy CHROM/POS/REF/ALT[list]
        con.execute(f"""
            CREATE OR REPLACE TABLE cohort AS
            SELECT DISTINCT replace(CAST(CHROM AS VARCHAR),'chr','') AS chrom, CAST(POS AS BIGINT) AS pos,
                   REF AS ref, t.alt_x AS alt,
                   replace(CAST(CHROM AS VARCHAR),'chr','')||'-'||CAST(POS AS BIGINT)||'-'||REF||'-'||t.alt_x AS variant_key
            FROM read_parquet('{parquet}'), UNNEST(ALT) AS t(alt_x)
        """)
    return con.execute("SELECT count(*) FROM cohort").fetchone()[0]


# --------------------------------------------------------------------------------------------------
# 2. reference tables (loaded once) + the LOCAL bulk classification
# --------------------------------------------------------------------------------------------------
def load_reference(con: duckdb.DuckDBPyConnection, cache: str, t: _T) -> None:
    """Load every LOCAL reference table from the cache. The exact `clinvar` lookup table is parsed from
    variant_summary; `clinvar_prot` is read from the cache parquet (fast). gnomAD constraint + exon model
    follow. All owned, all offline — no external quota."""
    cache = Path(cache)
    prot = cache / "clinvar_prot.parquet"
    vs = cache / "variant_summary.txt.gz"

    # clinvar_prot from cache parquet (PS1/PM5) — load_clinvar's fast path when the parquet exists.
    t.step("load clinvar_prot (cache parquet)",
           lambda: (load_clinvar(con, str(vs), cache_parquet=str(prot)),
                    f"prot_rows={con.execute('SELECT count(*) FROM clinvar_prot').fetchone()[0]}")[1])

    # exact-lookup `clinvar` table: load_clinvar skips it on the cache-parquet path, so build it here from
    # variant_summary (the one 440MB parse). Mirrors acmg.clinvar's exact-table SQL.
    def _exact():
        from acmg.clinvar import _STARS   # ReviewStatus -> gold stars 0-4 (ClinVar evidence is only as good as its review level)
        con.execute(f"""
            CREATE OR REPLACE TABLE clinvar AS
            SELECT "Chromosome" AS chrom, TRY_CAST("PositionVCF" AS BIGINT) AS pos,
                   "ReferenceAlleleVCF" AS ref, "AlternateAlleleVCF" AS alt, "ClinicalSignificance" AS clnsig,
                   {_STARS} AS review_stars, TRY_CAST("NumberSubmitters" AS INTEGER) AS n_submitters,
                   ("ClinicalSignificance" ILIKE '%conflict%') AS is_conflict
            FROM read_csv('{vs}', delim='\t', header=true, quote='', ignore_errors=true,
                          types={{'Chromosome':'VARCHAR','PositionVCF':'VARCHAR','ReferenceAlleleVCF':'VARCHAR',
                                  'AlternateAlleleVCF':'VARCHAR','ClinicalSignificance':'VARCHAR','Assembly':'VARCHAR',
                                  'ReviewStatus':'VARCHAR','NumberSubmitters':'VARCHAR'}})
            WHERE "Assembly" = 'GRCh37'
        """)
        return f"clinvar_rows={con.execute('SELECT count(*) FROM clinvar').fetchone()[0]}"
    t.step("load clinvar exact + review stars (parse variant_summary 440MB gz)", _exact)

    # CLNSIGCONF conflict COMPOSITION from the ClinVar VCF: what a 'Conflicting' call conflicts BETWEEN — the
    # subtlety that separates a P/LP-vs-VUS conflict (still informative) from a P-vs-Benign one. Needs the VCF.
    vcf = cache / "clinvar_grch37.vcf.gz"
    if vcf.exists():
        def _conf():
            con.execute("LOAD duckhts")
            con.execute(f"""
                CREATE OR REPLACE TABLE clinvar_conf AS
                SELECT replace(CAST(CHROM AS VARCHAR),'chr','')||'-'||CAST(POS AS VARCHAR)||'-'||REF||'-'||ALT[1] AS variant_key,
                       array_to_string(INFO_CLNSIGCONF,';') AS clnsigconf
                FROM read_bcf('{vcf}') WHERE len(ALT)=1 AND len(array_to_string(INFO_CLNSIGCONF,';')) > 0
            """)
            return f"conflict_comp_rows={con.execute('SELECT count(*) FROM clinvar_conf').fetchone()[0]}"
        t.step("ClinVar CLNSIGCONF conflict composition (VCF)", _conf)

    t.step("load gnomAD constraint",
           lambda: (load_constraint(con, str(cache / "gnomad_constraint.txt.gz")),
                    f"genes={con.execute('SELECT count(*) FROM gene_constraint').fetchone()[0]}")[1])

    t.step("load exon model (GTF via duckhts)",
           lambda: (load_exons(con, str(cache / "gencode.lift37.gtf.gz")),
                    f"exons={con.execute('SELECT count(*) FROM exon').fetchone()[0]}")[1])


def local_classify(con: duckdb.DuckDBPyConnection, t: _T) -> pd.DataFrame:
    """The part that needs NO VEP: join the exact ClinVar classification onto every cohort locus. Returns a
    partial `annotations` frame (variant_key + clinvar_clnsig) that seeds the kernel table; the VEP-derived
    columns are filled by join_vep_derived()."""
    def _j():
        has_conf = con.execute("SELECT count(*) FROM information_schema.tables WHERE table_name='clinvar_conf'").fetchone()[0] > 0
        conf_join = "LEFT JOIN clinvar_conf cf USING (variant_key)" if has_conf else ""
        conf_col = "cf.clnsigconf" if has_conf else "CAST(NULL AS VARCHAR)"
        # A variant can have several ClinVar records (multiple VariationIDs / the conflict join at one locus), so
        # this join FANS OUT — dedupe to ONE row per variant_key, keeping the most authoritative (highest review
        # stars, then most submitters). Without this, cohort_clinvar has duplicate variant_keys and the downstream
        # `local.set_index('variant_key')` -> Series.map raises InvalidIndexError (non-unique index).
        con.execute(f"""
            CREATE OR REPLACE TABLE cohort_clinvar AS
            SELECT variant_key, clinvar_clnsig, clinvar_review_stars, clinvar_n_submitters,
                   clinvar_is_conflict, clinvar_conflict_composition
            FROM (
                SELECT c.variant_key, cv.clnsig AS clinvar_clnsig, cv.review_stars AS clinvar_review_stars,
                       cv.n_submitters AS clinvar_n_submitters, cv.is_conflict AS clinvar_is_conflict,
                       {conf_col} AS clinvar_conflict_composition,
                       row_number() OVER (PARTITION BY c.variant_key
                           ORDER BY cv.review_stars DESC NULLS LAST, cv.n_submitters DESC NULLS LAST) AS _rn
                FROM cohort c LEFT JOIN clinvar cv
                  ON cv.chrom = c.chrom AND cv.pos = c.pos AND cv.ref = c.ref AND cv.alt = c.alt
                {conf_join}
            ) WHERE _rn = 1
        """)
        m = con.execute("SELECT count(*) FROM cohort_clinvar WHERE clinvar_clnsig IS NOT NULL").fetchone()[0]
        s2 = con.execute("SELECT count(*) FROM cohort_clinvar WHERE clinvar_review_stars >= 2").fetchone()[0]
        return f"clinvar_matched={m} (>=2-star: {s2})"
    t.step("ClinVar exact classification + stars + conflict-composition join @cohort", _j)
    return con.execute("""SELECT variant_key, clinvar_clnsig, clinvar_review_stars, clinvar_n_submitters,
                          clinvar_is_conflict, clinvar_conflict_composition FROM cohort_clinvar""").df()


# --------------------------------------------------------------------------------------------------
# 3. VEP-derived joins (wired; consume a VEP-annotated table)
# --------------------------------------------------------------------------------------------------
def join_site_scores(con: duckdb.DuckDBPyConnection, ann: pd.DataFrame,
                     site_scores: str = ".cache/site_scores.parquet") -> pd.DataFrame:
    """Fill `filtering_af`, `revel`, `spliceai` from the remote-tabix site_scores.parquet
    (scripts/build_site_scores.py — gnomAD AF over HTTP range reads, etc.), keyed by variant_key
    (chrom-pos-ref-alt). These reference site scores are authoritative for the frequency and in-silico
    concepts, so they fill/override the VEP-provided values where present; a missing score stays NULL and
    the kernel abstains. No-op if the parquet is absent (VEP values, if any, are kept)."""
    if not Path(site_scores).exists():
        return ann
    ann = ann.copy()
    ss = con.execute(f"""
        SELECT chrom || '-' || CAST(pos AS VARCHAR) || '-' || ref || '-' || alt AS variant_key,
               filtering_af, revel, spliceai
        FROM read_parquet('{site_scores}')
    """).df().drop_duplicates("variant_key").set_index("variant_key")
    for col in ("filtering_af", "revel", "spliceai"):
        vals = ann["variant_key"].map(ss[col])
        # store/VEP values are authoritative (VEP MAX_AF = max-pop); site_scores only FILLS gaps (fallback)
        ann[col] = vals if col not in ann.columns else ann[col].where(ann[col].notna(), vals)
    return ann


def join_vep_derived(con: duckdb.DuckDBPyConnection, vep: pd.DataFrame,
                     site_scores: str = ".cache/site_scores.parquet") -> pd.DataFrame:
    """Given a VEP-annotated frame with columns
        variant_key, gene, consequence, filtering_af, revel, spliceai, protein_pos, alt_aa1, transcript_id
    fill the LOCAL evidence (PS1/PM5, gnomAD mis_z, NMD escape) via DuckDB joins and return the full kernel
    `annotations` frame. This is the same wiring as acmg.annotate.annotate_hybrid, run at cohort scale.
    The site-level scores (filtering_af/revel/spliceai) are filled from the remote-tabix site_scores.parquet
    when present (join_site_scores)."""
    ann = join_site_scores(con, vep.copy(), site_scores)
    for col in REQUIRED_COLS:
        if col not in ann.columns:
            ann[col] = None

    # PS1/PM5 — same-aa / same-codon P-LP (needs gene + protein_pos + alt_aa1)
    q = ann.dropna(subset=["gene", "protein_pos", "alt_aa1"])[["variant_key", "gene", "protein_pos", "alt_aa1"]].copy()
    if len(q):
        q["protein_pos"] = q["protein_pos"].astype(int)
        cv = ps1_pm5(con, q).drop_duplicates("variant_key").set_index("variant_key")   # fail-safe: unique index for .map
        ann["clinvar_same_aa"] = ann["variant_key"].map(cv["clinvar_same_aa"]).fillna(0).astype(int)
        ann["clinvar_same_codon_lp"] = ann["variant_key"].map(cv["clinvar_same_codon_lp"]).fillna(0).astype(int)

    # PP2 — gnomAD mis_z by gene
    ann = add_mis_z(con, ann)

    # PVS1 strength — NMD escape. VEP's NMD plugin already filled the PTC (stop_gained/frameshift) variants on its
    # own picked transcript (parse_vep_tsv); acmg.nmd fills only the REST (splice LoF / any gaps) from the GTF model.
    need = ann["nmd_escaping"].isna() if "nmd_escaping" in ann.columns else pd.Series(True, index=ann.index)
    nq = ann[need].dropna(subset=["transcript_id"]).copy()
    if len(nq):
        parts = nq["variant_key"].str.split("-", n=3, expand=True)
        nq["chrom"], nq["pos"] = parts[0], parts[1].astype("int64")
        # drop_duplicates: nmd_escaping LEFT-JOINs exon `hit` on variant_key, which can fan out (>1 exon row per
        # variant on some data) -> non-unique index -> Series.map InvalidIndexError. Keep one row per variant, fail-safe.
        res = nmd_escaping(con, nq[["variant_key", "chrom", "pos", "transcript_id"]]).drop_duplicates("variant_key").set_index("variant_key")
        filled = ann["variant_key"].map(res["nmd_escaping"])
        ann["nmd_escaping"] = ann["nmd_escaping"].where(ann["nmd_escaping"].notna(), filled) if "nmd_escaping" in ann.columns else filled

    return ann[REQUIRED_COLS]


# --------------------------------------------------------------------------------------------------
# VEP offline hook + REST fallback
# --------------------------------------------------------------------------------------------------
def vep_offline_command(vcf: str, out_tsv: str, *, cache: str = "/root/.vep", assembly: str = "GRCh37",
                        fasta: str = "", vep_bin: str = "vep", fork: int = 8, mode: str = "database",
                        gtf: str = "", host: str = "grch37.ensembl.org", port: int = 3337) -> str:
    """VEP command producing the consequence/gene/protein/transcript TSV to join. We fetch ONLY the columns we
    need (`--fields`), over OUR ~477k exomic loci (the input VCF), never the whole genome. Three source modes:

      - "database" (DEFAULT): `--database --host grch37.ensembl.org --port 3337` queries Ensembl's public
        GRCh37 MySQL live — NO 15GB cache download. Best for a bounded exomic set; network-bound, so cache the
        result. (DuckDB's `mysql` extension can also ATTACH this host for pure-SQL model lookups.)
      - "gtf": `--gtf <our GENCODE v46lift37 gtf> --fasta <GRCh37 fasta>` — ONE gene model, identical to the
        NMD/PS1/PM5 joins (no cache/DB version mismatch). The "gff + supplementary" path.
      - "cache": `--offline --cache` (needs the ~15GB cache on disk).

    Frequency is NOT taken from VEP here (VEP's `gnomADe_AF` is a GLOBAL AF; the kernel needs MAX-POPULATION /
    FAF — the rs334 point). filtering_af is filled from our `variant_frequency` (per-ancestry) instead.
    `--mane_select --pick` fixes the ONE reported transcript. `--fork {fork}` parallelises (nproc≈20).
    `vep_bin` may point at a bioconda env, e.g. "micromamba run -n vep vep"."""
    if mode == "database":
        source = f"--database --host {host} --port {port} "
    elif mode == "gtf":  # branch on mode ONLY — the gtf arg has a default, so `or gtf` hijacked --mode cache
        source = f"--gtf {gtf} --fasta {fasta} "
    else:
        source = f"--offline --cache --dir_cache {cache} --assembly {assembly} " + (f"--fasta {fasta} " if fasta else "")
    forkflag = f"--fork {fork} " if fork and fork > 1 else ""
    return (
        f"{vep_bin} {source}{forkflag}"
        f"--input_file {vcf} --output_file {out_tsv} --tab --force_overwrite --no_stats "
        f"--symbol --mane_select --pick --canonical "
        f"--fields Uploaded_variation,Location,SYMBOL,Consequence,Amino_acids,Protein_position,Feature"
    )


def write_cohort_vcf(con: duckdb.DuckDBPyConnection, out_vcf: str) -> str:
    """Emit the distinct cohort loci as a minimal VCF for `vep --offline`."""
    # Sort GENOMICALLY (chrom then numeric pos). `ORDER BY 1` on the concatenated line string sorted positions
    # LEXICALLY ("1000" < "999"), which VEP rejects as unsorted — see write_loci_vcf.
    con.execute(f"""
        COPY (
          WITH v(o, c_num, c_str, p, line) AS (
            SELECT 0, NULL, NULL, NULL, '##fileformat=VCFv4.2'
            UNION ALL SELECT 1, NULL, NULL, NULL, '#CHROM'||chr(9)||'POS'||chr(9)||'ID'||chr(9)||'REF'||chr(9)||'ALT'||chr(9)||'QUAL'||chr(9)||'FILTER'||chr(9)||'INFO'
            UNION ALL SELECT 2, TRY_CAST(replace(CAST(chrom AS VARCHAR),'chr','') AS INTEGER),
                               replace(CAST(chrom AS VARCHAR),'chr',''), CAST(pos AS BIGINT),
                               chrom || chr(9) || CAST(pos AS VARCHAR) || chr(9) || '.' || chr(9) || ref || chr(9) || alt
                               || chr(9) || '.' || chr(9) || '.' || chr(9) || '.'
            FROM cohort
          ) SELECT line FROM v ORDER BY o, c_num NULLS LAST, c_str, p
        ) TO '{out_vcf}' (HEADER false, QUOTE '', DELIMITER '\n')
    """)
    return out_vcf


# --------------------------------------------------------------------------------------------------
# Dockerized VEP 116 (release matches the merged cache) + INCREMENTAL annotation store
# --------------------------------------------------------------------------------------------------
VEP_IMAGE = "ensemblorg/ensembl-vep:release_116.0"
VEP_RELEASE = "116"          # source-release tag stored with each cached VEP annotation
# MAX_AF = the MAX-POPULATION allele frequency (1000G/ESP/gnomAD) — the FAF-like value BA1/BS1 need (rs334 point),
# NOT the diluted global AF; served offline from the merged cache. gnomADe_AF/MAX_AF_POPS kept for provenance.
_VEP_FIELDS = "Uploaded_variation,Location,SYMBOL,Consequence,Amino_acids,Protein_position,Feature,MAX_AF,gnomADe_AF,MAX_AF_POPS,NMD"


def write_loci_vcf(con: duckdb.DuckDBPyConnection, loci: pd.DataFrame, out_vcf: str) -> str:
    """Minimal GRCh37 VCF (header + chrom-without-'chr' body) written by DuckDB COPY — no Python row loop."""
    con.register("_loci", loci)
    # VEP REQUIRES the input sorted by chromosome then position (the NMD plugin / merged cache reject an unsorted
    # file with exit 255). Carry chrom/pos sort keys and order the body genomically: numeric chromosomes first in
    # order (1..22), then X/Y/MT, positions ascending within each — no interleaving, positions monotone.
    con.execute(f"""COPY (
        WITH v(o, c_num, c_str, p, line) AS (
            SELECT 0, NULL, NULL, NULL, '##fileformat=VCFv4.2'
            UNION ALL SELECT 1, NULL, NULL, NULL, '#CHROM'||chr(9)||'POS'||chr(9)||'ID'||chr(9)||'REF'||chr(9)||'ALT'||chr(9)||'QUAL'||chr(9)||'FILTER'||chr(9)||'INFO'
            UNION ALL SELECT 2, TRY_CAST(replace(CAST(chrom AS VARCHAR),'chr','') AS INTEGER),
                               replace(CAST(chrom AS VARCHAR),'chr',''), CAST(pos AS BIGINT),
                               replace(CAST(chrom AS VARCHAR),'chr','')||chr(9)||CAST(pos AS VARCHAR)||chr(9)||'.'||chr(9)
                               ||ref||chr(9)||alt||chr(9)||'.'||chr(9)||'.'||chr(9)||'.' FROM _loci
        ) SELECT line FROM v ORDER BY o, c_num NULLS LAST, c_str, p
    ) TO '{out_vcf}' (HEADER false, QUOTE '')""")
    con.unregister("_loci")
    return out_vcf


def run_vep_docker(vcf: str, out_tsv: str, *, cache: str = ".cache", image: str = VEP_IMAGE,
                   assembly: str = "GRCh37", cache_version: int = 116, fork: int = 16) -> str:
    """Run Ensembl VEP OFFLINE via Docker over our loci — reproducible: same image + merged GRCh37 cache under
    <cache>/vep, --mane_select --pick, only the columns the joins need. `vcf`/`out_tsv` live under <cache>
    (mounted /data). The scripted, cohort-scale form of scripts/vep116.sh."""
    import subprocess
    cdir = str(Path(cache).resolve())
    subprocess.run([
        "docker", "run", "--rm", "-u", "0:0", "-v", f"{cdir}:/data", image,
        "vep", "--offline", "--cache", "--merged", "--dir_cache", "/data/vep",
        "--assembly", assembly, "--cache_version", str(cache_version), "--use_given_ref",
        "--mane_select", "--pick", "--symbol", "--canonical", "--fork", str(fork),
        "--af_gnomade", "--max_af",   # max-population AF (BA1/BS1) + gnomAD exomes AF, from the offline cache
        "--dir_plugins", "/plugins", "--plugin", "NMD",   # NMD escape (PVS1 strength) on VEP's own picked transcript
        "--tab", "--no_stats", "--force_overwrite", "--fields", _VEP_FIELDS,
        "--input_file", f"/data/{Path(vcf).name}", "--output_file", f"/data/{Path(out_tsv).name}",
    ], check=True)
    return out_tsv


def parse_vep_tsv(con: duckdb.DuckDBPyConnection, path: str) -> pd.DataFrame:
    """Parse a VEP --tab TSV (our --fields) into the join frame IN DuckDB (read_csv + SQL string fns, no Python
    line loop). Columns: variant_key, gene, consequence, protein_pos, alt_aa1, transcript_id. `##`/`#` header
    lines are null-padded and filtered. (Only the VEP-term -> kernel-vocab map stays a vectorised .map.)"""
    cols = {f"column{i}": "VARCHAR" for i in range(11)}   # uv, loc, sym, csq, aa, pp, feat, MAX_AF, gnomADe_AF, MAX_AF_POPS, NMD
    df = con.execute(f"""
        WITH raw AS (
            SELECT column0 AS uv, column2 AS sym, column3 AS csq, column4 AS aa, column5 AS pp, column6 AS feat,
                   column7 AS max_af, column8 AS gnomade_af, column9 AS max_af_pop, column10 AS nmd
            FROM read_csv('{path}', delim=chr(9), header=false, all_varchar=true, null_padding=true,
                          ignore_errors=true, columns={cols})
            WHERE column0 NOT LIKE '#%'
        )
        SELECT split_part(uv,'_',1)||'-'||split_part(uv,'_',2)||'-'
               ||split_part(split_part(uv,'_',3),'/',1)||'-'||split_part(split_part(uv,'_',3),'/',2) AS variant_key,
               nullif(sym,'-') AS gene, csq,
               TRY_CAST(nullif(pp,'-') AS INTEGER) AS protein_pos,
               CASE WHEN aa LIKE '%/%' THEN split_part(aa,'/',2) END AS alt_aa1,
               nullif(feat,'-') AS transcript_id,
               TRY_CAST(nullif(split_part(max_af,'&',1),'-') AS DOUBLE) AS filtering_af,   -- MAX-POPULATION AF (BA1/BS1)
               TRY_CAST(nullif(gnomade_af,'-') AS DOUBLE) AS gnomad_af,
               nullif(max_af_pop,'-') AS max_af_pop,
               -- NMD (PVS1 strength) from VEP's NMD plugin on the picked transcript. Only PTC consequences get it:
               -- escaping -> 1, PTC-but-triggering -> 0, non-PTC (splice/other) -> NULL (acmg.nmd fills as fallback).
               CASE WHEN nmd = 'NMD_escaping_variant' THEN 1
                    WHEN csq LIKE '%stop_gained%' OR csq LIKE '%frameshift%' THEN 0
                    ELSE NULL END AS nmd_escaping
        FROM raw
    """).df()
    df["consequence"] = df["csq"].map(most_severe_kernel_consequence)
    return df[["variant_key", "gene", "consequence", "protein_pos", "alt_aa1", "transcript_id",
               "filtering_af", "gnomad_af", "max_af_pop", "nmd_escaping"]]


# --------------------------------------------------------------------------------------------------
# The annotation RELEASE datalake (MULTI-SOURCE, incremental). ONE row per variant_key; EACH intrinsic
# annotation source writes its OWN columns tagged with its OWN source release:
#   VEP 116      -> gene/consequence/protein_pos/alt_aa1/transcript_id (+ vep_release)
#   REVEL v1.3   -> revel (+ revel_release)      OpenSpliceAI -> spliceai (+ spliceai_release)
#   CADD v1.7    -> cadd_phred (+ cadd_release)
# The INCREMENT is PER SOURCE: a variant VEP-annotated but not yet OpenSpliceAI-scored is in the OpenSpliceAI
# increment. Each source runs on its OWN increment and MERGES its columns in; the datalake is reused across
# cohorts/runs. The SAME functions serve CNVs — a parallel store keyed by cnv_key (pass key=/cohort=).
# Generalises the README's "fixed source release, only the increment sent to the tool, accumulate in a datalake
# versioned per release" to EVERY source. (This is what the old REST 'novel-set' fallback tried to be when the
# local cache was minimal; with a genome-wide OFFLINE cache there is NO novel set, so REST demotes to the
# interactive single-variant path + the rsID evidence arm — see rest_fallback.)
SOURCE_RELEASE = {"vep": VEP_RELEASE, "revel": "v1.3", "spliceai": "openspliceai", "cadd": "v1.7"}


def store_increment(con: duckdb.DuckDBPyConnection, store: str, release_col: str, release: str,
                    *, cohort: str = "cohort", key: str = "variant_key") -> pd.DataFrame:
    """Rows of `cohort` a SOURCE has not annotated at `release` — its increment (missing row / missing this
    source's columns / stale source release). `release_col` e.g. 'vep_release' / 'revel_release' / 'cadd_release'."""
    if not Path(store).exists():
        return con.execute(f"SELECT * FROM {cohort}").df()
    cols = list(con.execute(f"SELECT * FROM read_parquet('{store}') LIMIT 0").df().columns)
    if release_col not in cols:                       # this source has never run against the store
        return con.execute(f"SELECT * FROM {cohort}").df()
    return con.execute(f"""
        SELECT c.* FROM {cohort} c
        LEFT JOIN read_parquet('{store}') s ON s.{key} = c.{key}
        WHERE s.{key} IS NULL OR s.{release_col} IS NULL OR s.{release_col} <> '{release}'
    """).df()


def store_merge(con: duckdb.DuckDBPyConnection, store: str, new_df: pd.DataFrame, release_col: str,
                release: str, *, key: str = "variant_key") -> int:
    """Merge ONE source's columns (new_df = key + that source's columns) into the datalake by `key`, tagging
    `release_col`. New keys get a row; existing keys get these columns updated; other sources' columns preserved."""
    new_df = new_df.copy(); new_df[release_col] = release
    con.register("_new", new_df)
    if not Path(store).exists():
        con.execute(f"COPY (SELECT * FROM _new) TO '{store}' (FORMAT parquet)")
    else:
        con.execute(f"CREATE OR REPLACE VIEW _old AS SELECT * FROM read_parquet('{store}')")
        oldc = list(con.execute("SELECT * FROM _old LIMIT 0").df().columns)
        newc = list(new_df.columns)
        sel = [f'coalesce(n."{key}", o."{key}") AS "{key}"']
        for c in oldc + [c for c in newc if c not in oldc]:
            if c == key:
                continue
            if c in newc and c in oldc:
                sel.append(f'coalesce(n."{c}", o."{c}") AS "{c}"')     # the source's own columns win where present
            elif c in newc:
                sel.append(f'n."{c}" AS "{c}"')
            else:
                sel.append(f'o."{c}" AS "{c}"')                        # a different source's columns, preserved
        con.execute(f"""COPY (SELECT {', '.join(sel)}
            FROM _old o FULL OUTER JOIN _new n ON n."{key}" = o."{key}"
        ) TO '{store}.tmp' (FORMAT parquet)""")
        Path(f"{store}.tmp").replace(store)
    con.unregister("_new")
    return con.execute(f"SELECT count(*) FROM read_parquet('{store}')").fetchone()[0]


def store_read_cohort(con: duckdb.DuckDBPyConnection, store: str, *, cohort: str = "cohort",
                      key: str = "variant_key", cols: "list | None" = None) -> pd.DataFrame:
    """Read the cohort's merged intrinsic annotations from the datalake (all sources)."""
    have = list(con.execute(f"SELECT * FROM read_parquet('{store}') LIMIT 0").df().columns)
    want = [c for c in (cols or ["variant_key", "gene", "consequence", "protein_pos", "alt_aa1", "transcript_id",
                                 "filtering_af", "gnomad_af", "max_af_pop", "nmd_escaping", "revel", "spliceai", "cadd_phred"]) if c in have]
    return con.execute(f"""SELECT {', '.join('s."' + c + '"' for c in want)}
        FROM read_parquet('{store}') s JOIN {cohort} c ON c.{key} = s.{key}""").df()


def annotate_lookup(con: duckdb.DuckDBPyConnection, loci: pd.DataFrame, parquet: str, col: str) -> pd.DataFrame:
    """Annotate the increment `loci` by looking up a precomputed static file (REVEL v1.3 / OpenSpliceAI scores /
    CADD v1.7 slice), keyed chrom-pos-ref-alt -> (variant_key, <col>). A lookup IS the 'annotation' for these
    sources; they still live in the incremental datalake so the next cohort reuses them. Empty frame if absent."""
    if not Path(parquet).exists() or not len(loci):
        return pd.DataFrame(columns=["variant_key", col])
    con.register("_inc", loci[["variant_key"]])
    return con.execute(f"""SELECT p.chrom||'-'||p.pos||'-'||p.ref||'-'||p.alt AS variant_key, p."{col}"
        FROM read_parquet('{parquet}') p JOIN _inc i ON i.variant_key = p.chrom||'-'||p.pos||'-'||p.ref||'-'||p.alt""").df()


def rest_fallback(novel_variants: list[str], **kw) -> pd.DataFrame:
    """Ensembl VEP-REST — the INTERACTIVE / no-setup single-variant path, NOT a scale component. With the
    genome-wide OFFLINE cache + the incremental datalake there is no 'novel' consequence set for REST to fill;
    REST is for a quick one-off lookup without the 25GB cache + Docker, and the evidence arm (acmg.evidence,
    Variation REST by rsID) uses it for clinical-significance / phenotype / citation gathering. See docs/scale.md."""
    from acmg.vep_rest import vep_annotate
    return vep_annotate(novel_variants, **kw)


# --------------------------------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--parquet", default="/root/bioconnect/prepared/snv.parquet")
    ap.add_argument("--cache", default=".cache")
    ap.add_argument("--site-scores", default=".cache/site_scores.parquet",
                    help="remote-tabix site scores (scripts/build_site_scores.py) — fills filtering_af/revel/spliceai")
    ap.add_argument("--out", default=".cache/cohort_annotations.parquet", help="write the kernel annotations parquet here")
    ap.add_argument("--store", default=".cache/anno_store.parquet",
                    help="MULTI-SOURCE incremental annotation datalake (variant_key x VEP/REVEL/SpliceAI/CADD, each "
                         "with its own release tag); each source annotates only its increment; reused across runs")
    ap.add_argument("--revel", default=".cache/revel_grch37.parquet", help="local REVEL v1.3 parquet (PP3/BP4 source)")
    ap.add_argument("--spliceai", default=".cache/spliceai_openspliceai.parquet", help="OpenSpliceAI scores parquet (splice PP3/BP7)")
    ap.add_argument("--cadd", default=".cache/cadd_exome_slice.parquet", help="local CADD v1.7 exon slice (scripts/build_cadd_slice.sh)")
    ap.add_argument("--classify", action="store_true", help="run the SQL kernel and write classifications too")
    ap.add_argument("--fork", type=int, default=16, help="offline-VEP parallel workers (--fork); ~nproc (20 here)")
    ap.add_argument("--vep-bin", default="vep", help='native vep launcher for --mode database/gtf/cache')
    ap.add_argument("--mode", default="docker", choices=["docker", "database", "gtf", "cache"],
                    help="VEP source: docker (dockerized VEP 116 offline over <cache>/vep — default, reproducible), "
                         "database (Ensembl public MySQL), gtf (our GENCODE), cache (native offline)")
    ap.add_argument("--gtf", default=".cache/gencode.lift37.gtf.gz", help="GENCODE GTF for --mode gtf")
    args = ap.parse_args()

    con = duckdb.connect()
    t = _T()

    print("== cohort ==")
    n = t.step("build cohort (dedup to distinct loci)",
               lambda: f"distinct_loci={build_cohort(con, args.parquet)}")
    n_loci = con.execute("SELECT count(*) FROM cohort").fetchone()[0]

    print("== load LOCAL reference tables (owned, offline, no quota) ==")
    load_reference(con, args.cache, t)

    # ClinGen gene_curation (Gene-Disease Validity + Dosage) — the PVS1 ENTRY GATE (LoF must be an ESTABLISHED
    # disease mechanism, ClinGen HI3) + the gene-disease-validity cap. Without it PVS1 ABSTAINS for every LoF
    # variant (no gate opens) — so it MUST be loaded + passed to classify(), else PVS1 never fires.
    from acmg import clingen
    gene_curation = None
    try:
        cg = clingen.download(args.cache)
        clingen.load_gene_curation(con, cg["gene_validity"], cg["dosage"])
        gene_curation = con.execute("SELECT * FROM gene_curation").df()
        print(f"  ClinGen gene_curation: {len(gene_curation)} gene-disease rows (PVS1 gate + validity cap)")
    except Exception as e:
        print(f"  ClinGen gene_curation unavailable ({e}); PVS1 will abstain (fail-safe)")

    print(f"== LOCAL bulk classification @{n_loci} loci ==")
    local = local_classify(con, t)

    print(f"\nLOCAL annotation total: {t.total():.2f}s wall for {n_loci} distinct loci "
          f"({int(local['clinvar_clnsig'].notna().sum())} with a ClinVar exact classification).")

    # === INCREMENTAL MULTI-SOURCE annotation: each source annotates only ITS increment, merges into the datalake ===
    print(f"\n== incremental annotation datalake ({args.store}) — per-source increments ==")
    # 1) VEP consequence/gene/protein/transcript — the one source that RUNS a tool (offline docker)
    if args.mode == "docker":
        inc = store_increment(con, args.store, "vep_release", SOURCE_RELEASE["vep"])
        print(f"  VEP      increment {len(inc)}/{n_loci} new loci")
        if len(inc):
            vcf = write_loci_vcf(con, inc, str(Path(args.cache) / "increment.vcf"))
            tsv = str(Path(args.cache) / "increment.vep.tsv")
            t.step(f"VEP offline (docker {VEP_IMAGE}, --fork {args.fork}) x{len(inc)}",
                   lambda: (run_vep_docker(vcf, tsv, cache=args.cache, fork=args.fork), "ok")[1])
            store_merge(con, args.store, parse_vep_tsv(con, tsv), "vep_release", SOURCE_RELEASE["vep"])
    else:
        vcf = write_cohort_vcf(con, str(Path(args.cache) / "cohort.vcf"))
        print(f"  --mode {args.mode}: run this reproducible command, then re-run with --mode docker (or feed the TSV):")
        print("   ", vep_offline_command(vcf, str(Path(args.cache) / "cohort.vep.tsv"), vep_bin=args.vep_bin,
                                         fork=args.fork, mode=args.mode, gtf=args.gtf))
    # 2) precomputed / on-device score sources — each looks up ITS increment from a static parquet, merges in
    for name, parquet, col, rc, rel in [("REVEL", args.revel, "revel", "revel_release", SOURCE_RELEASE["revel"]),
                                        ("SpliceAI", args.spliceai, "spliceai", "spliceai_release", SOURCE_RELEASE["spliceai"]),
                                        ("CADD", args.cadd, "cadd_phred", "cadd_release", SOURCE_RELEASE["cadd"])]:
        inc = store_increment(con, args.store, rc, rel)
        got = annotate_lookup(con, inc, parquet, col)
        avail = Path(parquet).name if Path(parquet).exists() else "(absent)"
        print(f"  {name:8} increment {len(inc)}/{n_loci}; merged {len(got)} from {avail}")
        if len(got):
            store_merge(con, args.store, got, rc, rel)

    vep = store_read_cohort(con, args.store) if Path(args.store).exists() else pd.DataFrame(
        columns=["variant_key", "gene", "consequence", "protein_pos", "alt_aa1", "transcript_id"])

    # === cheap CURRENT-release joins (PS1/PM5, mis_z, NMD, gnomAD AF) -> full kernel annotations ===
    print("\n== local evidence joins (PS1/PM5, mis_z, NMD) + gnomAD AF (site scores) ==")
    ann = join_vep_derived(con, vep, site_scores=args.site_scores) if len(vep) else join_site_scores(con, local, args.site_scores)
    # ClinVar aggregate classification (exact chrom-pos-ref-alt match) as a FIRST-CLASS column — the known-answer
    # overlay for prioritization (NOT a formal ACMG criterion; PS1/PM5 remain the mechanical same-aa/codon transfer).
    # Current-release join, recomputed each run (ClinVar is rolling).
    if "clinvar_clnsig" in local.columns:
        li = local.set_index("variant_key")
        ann["clinvar_classification"] = ann["variant_key"].map(li["clinvar_clnsig"])
        ann["clinvar_review_stars"] = ann["variant_key"].map(li["clinvar_review_stars"])          # 0-4 gold stars
        ann["clinvar_conflict_composition"] = ann["variant_key"].map(li["clinvar_conflict_composition"])  # CLNSIGCONF (P-vs-VUS vs P-vs-B)
    print(f"  annotations: {len(ann)} loci"
          + "".join(f"; {c} on {int(ann[c].notna().sum())}" for c in ("filtering_af", "revel", "spliceai", "clinvar_classification") if c in ann))

    con.register("_ann", ann)
    con.execute(f"COPY (SELECT * FROM _ann) TO '{args.out}' (FORMAT parquet)")
    print(f"\nwrote kernel annotations -> {args.out}  ({t.total():.1f}s total)")

    if args.classify and len(vep):
        from acmg.kernel import classify
        cls = classify(ann, con=duckdb.connect(), gene_curation=gene_curation,  # gene_curation opens the PVS1 gate
                       cadd_supporting=True)  # CADD supporting fallback on the non-REVEL/non-splice/non-ClinVar tail
        if "clinvar_classification" in ann.columns:   # carry the known-answer overlay (+ stars + conflict) into the classified output
            ai = ann.set_index("variant_key")
            for c in ("clinvar_classification", "clinvar_review_stars", "clinvar_conflict_composition"):
                cls[c] = cls["variant_key"].map(ai[c])
        out_cls = args.out.replace(".parquet", "") + "_classified.parquet"
        con.register("_cls", cls); con.execute(f"COPY (SELECT * FROM _cls) TO '{out_cls}' (FORMAT parquet)")
        pvs1 = int(cls["criteria"].str.contains("PVS1", na=False).sum())
        print(f"classified {len(cls)} -> {cls['acmg_class'].value_counts().to_dict()}; PVS1 fired on {pvs1}\nwrote -> {out_cls}")


if __name__ == "__main__":
    main()
