-- CNV / structural-variant classification (ClinGen/ACMG Riggs 2020 dosage) banded from a ClassifyCNV
-- Scoresheet loaded as `cnv_evidence`. The per-criterion columns (1A-B, 2A, 2B, ...) are the applied Riggs
-- evidence (the CNV analogue of the SNV kernel's applied_criteria); `total_score` is their signed sum.
-- Banding cutoffs are the ClinGen Riggs point thresholds (P >= 0.99, LP 0.90-0.98, VUS -0.89..0.89,
-- LB -0.90..-0.98, B <= -0.99) — data here, so a VCEP/panel could tune them, same as the SNV combine.
CREATE OR REPLACE VIEW cnv_classification AS
SELECT "VariantID" AS variant_id, "Chromosome" AS chromosome,
       "Start" AS start_pos, "End" AS end_pos, "Type" AS type, total_score,
       "Known or predicted dosage-sensitive genes" AS dosage_genes,   -- the HI/TS genes ClassifyCNV scored (its own
                                                                      -- interval x gene-model overlap; no separate GTF)
  CASE
    WHEN total_score >= 0.99 THEN 'Pathogenic'
    WHEN total_score >= 0.90 THEN 'Likely Pathogenic'
    WHEN total_score >  -0.90 THEN 'VUS'
    WHEN total_score >  -0.99 THEN 'Likely Benign'
    ELSE 'Benign'
  END AS acmg_class
FROM cnv_evidence;
