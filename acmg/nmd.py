"""NMD determination for PVS1 strength — from the gene model (duckhts read_gtf), in SQL.

PVS1 strength (Tayoun 2018 / autoPVS1, simplified): a predicted null (stop_gained / frameshift / splice)
is NMD-TRIGGERING -> PVS1_VeryStrong, unless it ESCAPES NMD -> PVS1_Strong. The simplified escape rule:
the premature stop is in the LAST exon, or within the last 50 nt of the PENULTIMATE exon (i.e. <50 nt from
the final exon-exon junction); a single-exon transcript escapes (no junction). `nmd_escaping` = 0 triggering
/ 1 escaping, which is exactly the column the kernel's PVS1 rule reads.

TRANSCRIPT CHOICE IS LOAD-BEARING. NMD depends on the transcript's exon structure, so the answer changes
with the transcript. Drive this with the SAME transcript the annotation used — MANE Select, else MANE Plus
Clinical, else GENCODE Primary / canonical (run VEP with --mane --mane_plus_clinical --pick). The MANE
transcript also fixes the RefSeq(NM_) vs Ensembl(ENST) split: MANE is the single per-gene transcript whose
Ensembl and RefSeq versions share IDENTICAL CDS coordinates — so a MANE-driven pipeline makes VEP's protein
positions (ENST) and ClinVar's (NM_ in variant_summary) directly comparable for PS1/PM5. Off MANE, protein
numbering can differ between transcripts and the PS1/PM5 join should be treated with caution.

MANE ON GRCh37: use the **GENCODE lift37** GTF (e.g. `gencode.v46lift37.basic.annotation.gtf.gz`) — it maps
the GRCh38 annotation, INCLUDING `tag "MANE_Select"` / `"MANE_Plus_Clinical"`, onto GRCh37 via gencode-
backmap. (A plain Ensembl-87 GRCh37 GTF predates MANE and has no tags.) With lift37, `select_transcript`
returns the true clinical transcript — e.g. BRCA1 -> ENST00000357654 (MANE Select).
"""
from __future__ import annotations
import duckdb

_ATTR = lambda key: f"regexp_extract(attributes, '{key} \"([^\"]+)\"', 1)"


def load_exons(con: duckdb.DuckDBPyConnection, gtf: str) -> None:
    """Build an `exon` table (+ transcript tags) from a GTF via duckhts read_gtf."""
    con.execute("LOAD duckhts")
    con.execute(f"""
        CREATE OR REPLACE TABLE exon AS
        SELECT split_part({_ATTR('transcript_id')}, '.', 1)  AS transcript_id,  -- bare ENST (join VEP's)
               {_ATTR('gene_name')}          AS gene,
               TRY_CAST({_ATTR('exon_number')} AS INTEGER) AS exon_number,
               replace(CAST(seqname AS VARCHAR),'chr','') AS chrom,
               CAST(start AS BIGINT) AS start_pos, CAST("end" AS BIGINT) AS end_pos, strand,
               (attributes LIKE '%MANE_Select%')         AS mane_select,
               (attributes LIKE '%MANE_Plus_Clinical%')  AS mane_plus_clinical,
               (attributes LIKE '%tag "basic"%')          AS basic,
               (attributes LIKE '%appris_principal%')     AS appris_principal
        FROM read_gtf('{gtf}')
        WHERE feature = 'exon'
    """)


def select_transcript(con: duckdb.DuckDBPyConnection, gene: str) -> str | None:
    """MANE Select -> MANE Plus Clinical -> APPRIS principal / basic -> most exons. Best-effort on GRCh37."""
    row = con.execute(f"""
        SELECT transcript_id FROM exon WHERE gene = ?
        GROUP BY transcript_id, mane_select, mane_plus_clinical, appris_principal, basic
        ORDER BY bool_or(mane_select) DESC, bool_or(mane_plus_clinical) DESC,
                 bool_or(appris_principal) DESC, bool_or(basic) DESC, count(*) DESC
        LIMIT 1
    """, [gene]).fetchone()
    return row[0] if row else None


def nmd_escaping(con: duckdb.DuckDBPyConnection, variants: "pd.DataFrame") -> "pd.DataFrame":
    """variants: variant_key, chrom, pos, transcript_id (the annotation transcript — MANE-preferred).
    Returns variant_key + nmd_escaping (1 escaping / 0 triggering / NULL if the transcript/exon isn't found).
    Transcription-order aware (strand); single-exon transcripts escape."""
    con.register("_v", variants)
    out = con.execute("""
        WITH ord AS (
            SELECT e.transcript_id, e.strand, e.start_pos, e.end_pos,
                   row_number() OVER (PARTITION BY e.transcript_id
                       ORDER BY CASE WHEN e.strand = '+' THEN e.start_pos ELSE -e.start_pos END) AS tx_order,
                   count(*) OVER (PARTITION BY e.transcript_id) AS n
            FROM exon e JOIN (SELECT DISTINCT transcript_id FROM _v) t USING (transcript_id)
        ),
        hit AS (  -- the exon that contains each variant
            SELECT v.variant_key, o.tx_order, o.n, o.strand, o.start_pos, o.end_pos
            FROM _v v JOIN ord o ON o.transcript_id = v.transcript_id
                  AND CAST(v.pos AS BIGINT) BETWEEN o.start_pos AND o.end_pos
        )
        SELECT v.variant_key,
            CASE
                WHEN h.tx_order IS NULL THEN NULL                       -- variant not in a coding exon of this tx
                WHEN h.n = 1 THEN 1                                     -- single-exon transcript escapes NMD
                WHEN h.tx_order = h.n THEN 1                            -- in the last exon -> escapes
                WHEN h.tx_order = h.n - 1 AND (
                        (h.strand = '+' AND CAST(v.pos AS BIGINT) >= h.end_pos - 50) OR
                        (h.strand = '-' AND CAST(v.pos AS BIGINT) <= h.start_pos + 50)
                     ) THEN 1                                           -- <50nt from the last junction -> escapes
                ELSE 0                                                  -- NMD-triggering
            END AS nmd_escaping
        FROM _v v LEFT JOIN hit h ON h.variant_key = v.variant_key
    """).df()
    con.unregister("_v")
    return out
