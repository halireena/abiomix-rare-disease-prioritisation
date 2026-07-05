-- bioconnect-acmg: the points-based combine (Tavtigian 2020). SUM signed points per variant, band, then apply two
-- gates: (a) SNV/indel only — CNV/SV are NOT this kernel (Riggs 2020); (b) gene-disease validity CAP (PMC9035475):
-- Limited/Disputed/Refuted → ≤ VUS, Moderate → ≤ Likely Pathogenic. Starts from `annotations` (LEFT JOIN) so a
-- variant with NO criterion still classifies (VUS), not vanishes. BA1 is a stand-alone benign override.
CREATE OR REPLACE TABLE classification AS
WITH crit AS (
  -- DEDUP by (variant_key, criterion, points): each criterion counts ONCE. A duplicate spec/panel row (or any
  -- fan-out in annotations_spec) would otherwise emit the same criterion twice and SUM(points) would double it
  -- (e.g. PP3_Strong +4 -> +8) while `DISTINCT criterion` hid it. One criterion = one contribution.
  SELECT DISTINCT variant_key, criterion, points FROM applied_criteria
),
agg AS (
  SELECT a.variant_key, a.gene, a.variant_kind,
         COALESCE(SUM(c.points), 0)                    AS total_points,
         COALESCE(BOOL_OR(c.criterion = 'BA1'), FALSE) AS has_ba1,
         array_to_string(
           list_sort(array_agg(DISTINCT c.criterion) FILTER (WHERE c.criterion IS NOT NULL)), ', '
         )                                             AS criteria
  FROM annotations a
  LEFT JOIN crit c USING (variant_key)
  GROUP BY a.variant_key, a.gene, a.variant_kind
),
scored AS (
  SELECT variant_key, gene, variant_kind, total_points,
    CASE WHEN criteria IS NULL OR criteria = '' THEN '(none — abstained)' ELSE criteria END AS criteria,
    CASE
      WHEN has_ba1            THEN 'Benign'             -- BA1 stand-alone
      WHEN total_points >= 10 THEN 'Pathogenic'
      WHEN total_points >= 6  THEN 'Likely Pathogenic'
      WHEN total_points >= 0  THEN 'VUS'
      WHEN total_points >= -6 THEN 'Likely Benign'
      ELSE 'Benign'
    END AS raw_class,
    -- gene-disease validity for the disease CONTEXT. Stage 1 has no case-disease column, so it resolves the
    -- STRONGEST validity across the gene's diseases (Definitive > Strong > Moderate > Limited > Disputed >
    -- Refuted). This is fail-safe: the validity CAP downweights genes with only WEAK disease association, so a
    -- gene that is Definitive for ANY disease must NOT be capped by an unrelated Limited/Disputed row. (Picking
    -- an arbitrary disease — e.g. ORDER BY disease_id — wrongly caps TTN/MUTYH P/LP variants to VUS.) Production
    -- passes the case's disease_id and selects that exact gene-disease row. Scalar subquery = dup-safe.
    (SELECT gc.gene_disease_validity FROM gene_curation gc WHERE gc.gene = agg.gene
     ORDER BY CASE gc.gene_disease_validity
        WHEN 'Definitive' THEN 6 WHEN 'Strong' THEN 5 WHEN 'Moderate' THEN 4
        WHEN 'Limited' THEN 3 WHEN 'Disputed' THEN 2 WHEN 'Refuted' THEN 1 ELSE 0 END DESC
     LIMIT 1) AS validity
  FROM agg
)
SELECT variant_key, gene, variant_kind, total_points, criteria,
  CASE
    -- (a) SNV/indel gate: this kernel is not for CNV/SV (a NULL kind is treated as SNV for back-compat)
    WHEN variant_kind IS NOT NULL AND variant_kind NOT IN ('snv','indel')
      THEN 'Not evaluated (non-SNV/indel — see Riggs 2020)'
    -- (b) gene-disease validity cap: unknown/Definitive/Strong → uncapped
    WHEN validity IN ('Limited','Disputed','Refuted') AND raw_class IN ('Pathogenic','Likely Pathogenic')
      THEN 'VUS'
    WHEN validity = 'Moderate' AND raw_class = 'Pathogenic'
      THEN 'Likely Pathogenic'
    ELSE raw_class
  END AS acmg_class
FROM scored;
