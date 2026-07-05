"""CNV / structural-variant track — ClinGen/ACMG Riggs 2020 dosage scoring, DuckDB-composable.

The SNV/indel kernel returns "Not evaluated" for CNVs; this handles them, in the same shape as the SNV
kernel: an external scorer emits per-criterion evidence, we load it into DuckDB, and a Riggs points
*manifest* bands it. The scorer is **ClassifyCNV** (Genotek) — Python, offline, native GRCh37, and its
`Scoresheet.txt` has one column per Riggs criterion with its point value (the CNV analogue of
`applied_criteria`) plus the summed score. Academic/research license — fine for the sprint / a paper.

Setup once (not vendored — it ships a ~120 MB ClinGen dosage map):
    git clone https://github.com/Genotek/ClassifyCNV.git   # needs bedtools on PATH
Point `classifycnv_dir` at it (or set CLASSIFYCNV_DIR).
"""
from __future__ import annotations
import os
import subprocess
from pathlib import Path
import duckdb
from .ingest import _table_sql, _detect_fmt  # reuse the generic (format-agnostic) readers

_MANIFEST = Path(__file__).parent / "manifests" / "cnv_riggs.sql"
_DEFAULT_DIR = os.environ.get("CLASSIFYCNV_DIR", "/root/ClassifyCNV")

# Default map canonical -> source column for the challenge bundle's tidy CNV rows. Override for other layouts.
BUNDLE_CNV_COLUMNS = {"chrom": "CHROM", "start": "POS", "end": "INFO_END",
                      "svtype": "INFO_SVTYPE", "cn": "FORMAT_CN"}


def cnvs_from_source(con: duckdb.DuckDBPyConnection, source: str, columns: dict | None = None,
                     *, fmt: str | None = None, kind_filter: bool = True):
    """Read CNV rows from ANY source (parquet / tsv / csv / vcf / bcf) into {chrom,start,end,svtype,cn}.
    Format-agnostic like acmg.ingest; default column map = the bundle's tidy CNV names. Returns a DataFrame
    ready for classify_cnv / cnv_to_bed."""
    c = {**BUNDLE_CNV_COLUMNS, **(columns or {})}
    tbl = _table_sql(con, source, fmt or _detect_fmt(source))
    types = {r[0]: r[1] for r in con.execute(f"DESCRIBE SELECT * FROM {tbl}").fetchall()}
    sel = (f"replace(CAST(\"{c['chrom']}\" AS VARCHAR),'chr','') AS chrom, "
           f"CAST(\"{c['start']}\" AS BIGINT) AS start, CAST(\"{c['end']}\" AS BIGINT) AS \"end\"")
    sel += (f", CAST(\"{c['svtype']}\" AS VARCHAR) AS svtype" if c["svtype"] in types else ", CAST(NULL AS VARCHAR) AS svtype")
    sel += (f", TRY_CAST(\"{c['cn']}\" AS INTEGER) AS cn" if c["cn"] in types else ", CAST(NULL AS INTEGER) AS cn")
    where = "WHERE variant_kind = 'cnv'" if kind_filter and "variant_kind" in types else ""
    return con.execute(f"SELECT {sel} FROM {tbl} {where}").df()


def cnv_to_bed(cnvs, bed_path: str) -> int:
    """Write a ClassifyCNV BED (chr, start, end, DEL|DUP) from CNV rows. `cnvs`: an iterable of dicts /
    a DataFrame with chrom, start, end and either `svtype` (DEL/DUP) or `cn` (copy number; <2 -> DEL, >2 -> DUP)."""
    rows = cnvs.to_dict("records") if hasattr(cnvs, "to_dict") else list(cnvs)
    n = 0
    with open(bed_path, "w") as fh:
        for r in rows:
            sv = str(r.get("svtype") or "").upper()
            if sv not in ("DEL", "DUP"):
                cn = r.get("cn")
                sv = "DEL" if cn is not None and cn < 2 else "DUP" if cn is not None and cn > 2 else None
            if sv is None or r.get("start") is None or r.get("end") is None:
                continue
            chrom = str(r["chrom"]); chrom = chrom if chrom.startswith("chr") else "chr" + chrom
            fh.write(f"{chrom}\t{int(r['start'])}\t{int(r['end'])}\t{sv}\n")
            n += 1
    return n


def run_classifycnv(bed_path: str, outdir: str, *, build: str = "hg19", classifycnv_dir: str | None = None) -> str:
    """Run ClassifyCNV; returns the path to Scoresheet.txt."""
    d = classifycnv_dir or _DEFAULT_DIR
    subprocess.run(["python3", "ClassifyCNV.py", "--infile", os.path.abspath(bed_path),
                    "--GenomeBuild", build, "--outdir", os.path.abspath(outdir)],
                   cwd=d, check=True, capture_output=True, text=True)
    return os.path.join(outdir, "Scoresheet.txt")


def load_cnv_evidence(con: duckdb.DuckDBPyConnection, scoresheet: str) -> None:
    """Load a ClassifyCNV Scoresheet as `cnv_evidence` (one row per CNV, one column per Riggs criterion)."""
    con.execute(f"""CREATE OR REPLACE TABLE cnv_evidence AS
        SELECT * RENAME ("Total score" AS total_score)
        FROM read_csv('{scoresheet}', delim='\t', header=true, all_varchar=false)""")


def classify_cnv(cnvs, *, con: duckdb.DuckDBPyConnection | None = None, build: str = "hg19",
                 classifycnv_dir: str | None = None, workdir: str | None = None):
    """CNV rows -> Riggs classification. Returns a DataFrame: chromosome, start, end, type, total_score,
    acmg_class (5-tier, banded by the manifest). Runs ClassifyCNV then the Riggs SQL manifest."""
    import tempfile
    con = con or duckdb.connect()
    wd = workdir or tempfile.mkdtemp(prefix="cnv_")
    bed = os.path.join(wd, "cnvs.bed")
    if cnv_to_bed(cnvs, bed) == 0:
        return con.sql("SELECT NULL WHERE FALSE").df()
    sheet = run_classifycnv(bed, os.path.join(wd, "out"), build=build, classifycnv_dir=classifycnv_dir)
    load_cnv_evidence(con, sheet)
    con.execute(_MANIFEST.read_text())
    return con.execute("SELECT * FROM cnv_classification ORDER BY total_score DESC").df()
