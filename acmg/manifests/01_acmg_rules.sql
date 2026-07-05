-- bioconnect-acmg: the 13 literature-INDEPENDENT ACMG/AMP rules as ONE SQL statement (the rule kernel).
-- Reads `annotations_spec` (annotations + per-gene resolved thresholds from 00_spec.sql), emits
-- `applied_criteria(variant_key, criterion, direction, points, evidence)`.
-- Point values = Tavtigian 2020 (Supporting +/-1, Moderate +/-2, Strong +/-4, Very Strong +/-8).
-- Frequency thresholds are SPEC-RESOLVED per gene (ba1_maf/bs1_maf/pm2_maf) so a VCEP/panel exception is a data row.
-- A NULL frequency makes PM2 ABSTAIN — "no frequency data" is not "rare".
CREATE OR REPLACE TABLE applied_criteria AS
-- BA1: MAF above the (gene-resolved) stand-alone-benign threshold. The exception list (Ghosh 2018) is checked FIRST:
-- a known-pathogenic common variant (rs334 etc.) is NOT auto-benigned. Established pathogenicity beats frequency.
SELECT variant_key, 'BA1' AS criterion, 'B' AS direction, -8 AS points,
       'filtering AF (max-pop) ' || filtering_af || ' > ' || ba1_maf || ' (stand-alone benign)' AS evidence
FROM annotations_spec WHERE filtering_af > ba1_maf AND NOT EXISTS (SELECT 1 FROM ba1_exceptions e WHERE e.variant_key = annotations_spec.variant_key)
UNION ALL
-- BS1: between the strong-benign and stand-alone thresholds (same exception-list guard)
SELECT variant_key, 'BS1', 'B', -4, 'filtering AF (max-pop) ' || filtering_af || ' > ' || bs1_maf
FROM annotations_spec WHERE filtering_af > bs1_maf AND filtering_af <= ba1_maf AND NOT EXISTS (SELECT 1 FROM ba1_exceptions e WHERE e.variant_key = annotations_spec.variant_key)
UNION ALL
-- PM2_Supporting: at/below the (gene-resolved) rare threshold; ABSTAINS when frequency is unknown
SELECT variant_key, 'PM2_Supporting', 'P', 1, 'filtering AF (max-pop) ' || filtering_af || ' <= ' || pm2_maf
FROM annotations_spec WHERE filtering_af IS NOT NULL AND filtering_af <= pm2_maf
UNION ALL
-- PP2: missense in a missense-constrained gene (gnomAD v2 constraint Z >= 3.09)
SELECT variant_key, 'PP2', 'P', 1, 'missense constraint Z ' || gnomad_mis_z || ' >= 3.09'
FROM annotations_spec WHERE consequence = 'missense' AND gnomad_mis_z IS NOT NULL AND gnomad_mis_z >= 3.09
UNION ALL
-- PP3/BP4 for MISSENSE via a SINGLE genome-wide calibrated predictor (ClinGen/Pejaver 2022, PMC9748256): REVEL.
-- The current guideline uses ONE calibrated tool, NOT a consensus of several — stacking predictors to find the
-- strongest evidence is a documented bias. REVEL's calibrated thresholds reach Strong (pathogenic) / Very Strong
-- (benign). The indeterminate zone (0.290, 0.644) applies NO criterion. Swap `revel` for another calibrated tool's
-- score + thresholds (BayesDel, VEST4, AlphaMissense, ESM1b, VARITY) to change the predictor — one, not many.
SELECT variant_key,'BP4_VeryStrong','B',-8,'REVEL '||revel||' <= 0.003' FROM annotations_spec WHERE consequence='missense' AND revel IS NOT NULL AND revel <= 0.003
UNION ALL
SELECT variant_key,'BP4_Strong','B',-4,'REVEL '||revel||' in (0.003,0.016]' FROM annotations_spec WHERE consequence='missense' AND revel > 0.003 AND revel <= 0.016
UNION ALL
SELECT variant_key,'BP4_Moderate','B',-2,'REVEL '||revel||' in (0.016,0.183]' FROM annotations_spec WHERE consequence='missense' AND revel > 0.016 AND revel <= 0.183
UNION ALL
SELECT variant_key,'BP4','B',-1,'REVEL '||revel||' in (0.183,0.290]' FROM annotations_spec WHERE consequence='missense' AND revel > 0.183 AND revel <= 0.290
UNION ALL
SELECT variant_key,'PP3','P',1,'REVEL '||revel||' in [0.644,0.773)' FROM annotations_spec WHERE consequence='missense' AND revel >= 0.644 AND revel < 0.773
UNION ALL
SELECT variant_key,'PP3_Moderate','P',2,'REVEL '||revel||' in [0.773,0.932)' FROM annotations_spec WHERE consequence='missense' AND revel >= 0.773 AND revel < 0.932
UNION ALL
SELECT variant_key,'PP3_Strong','P',4,'REVEL '||revel||' >= 0.932' FROM annotations_spec WHERE consequence='missense' AND revel >= 0.932
UNION ALL
-- Splicing is a SEPARATE concept: its own single calibrated predictor (SpliceAI; ClinGen SVI Splicing, Walker 2023).
-- PP3 if SpliceAI >= 0.2, but ONLY for splice-relevant consequences — a missense already gets its ONE computational
-- PP3 from REVEL, so SpliceAI must not stack a second PP3 on it (concept-cap; one calibrated tool per concept).
SELECT variant_key,'PP3','P',1,'SpliceAI '||spliceai||' >= 0.2 (splice concept)'
FROM annotations_spec WHERE consequence IN ('splice_region','splice_donor','splice_acceptor','intron','synonymous') AND spliceai IS NOT NULL AND spliceai >= 0.2
UNION ALL
-- BP7: synonymous with SpliceAI < 0.2
SELECT variant_key,'BP7','B',-1,'synonymous, SpliceAI '||spliceai||' < 0.2' FROM annotations_spec WHERE consequence='synonymous' AND spliceai IS NOT NULL AND spliceai < 0.2
UNION ALL
-- PM4: in-frame protein-length change (derived from consequence, not a precomputed flag)
SELECT variant_key,'PM4','P',2,'in-frame protein-length change ('||consequence||')'
FROM annotations_spec WHERE consequence IN ('inframe_deletion','inframe_insertion','stop_lost')
UNION ALL
-- PS1: same amino-acid change as a known pathogenic variant (ClinVar)
SELECT variant_key,'PS1','P',4,'ClinVar: same aa change is pathogenic' FROM annotations_spec WHERE clinvar_same_aa = 1
UNION ALL
-- PM5: a different pathogenic missense at the same codon (ClinVar)
SELECT variant_key,'PM5','P',2,'ClinVar: pathogenic missense at same codon' FROM annotations_spec WHERE clinvar_same_codon_lp = 1
UNION ALL
-- PM1: missense in a mutational HOT SPOT / critical functional region — empirically, >=2 pathogenic ClinVar
-- missense within +/-3 codons AND ZERO benign there ('without benign variation'). Precomputed in
-- acmg.clinvar.pm1_hotspot (ClinVar clustering; ClinGen VCEP domain specs can augment). Abstains when NULL/0.
-- OPT-IN (00_spec.sql pm1_config): scores only when enabled; the pm1 flag is surfaced as context regardless.
SELECT variant_key,'PM1','P',2,'mutational hotspot / critical region (pathogenic ClinVar clustering, no benign)'
FROM annotations_spec WHERE consequence='missense' AND pm1 = 1 AND (SELECT enabled FROM pm1_config)
UNION ALL
-- PVS1: predicted null variant. Two things the kernel gets right that a naive rule misses:
--  (a) ENTRY GATE — PVS1 applies only where LoF is an ESTABLISHED disease mechanism for the gene (Tayoun 2018;
--      PMC9035475), i.e. ClinGen HI score 3. The JOIN to gene_curation enforces it; a gene that is uncurated or
--      LoF-tolerant gets NO PVS1 (abstain, never misfire). A null variant in a random gene is NOT PVS1.
--  (b) STRENGTH — derived by the kernel from the null mechanism + NMD, NOT consumed from a precomputed autoPVS1 black
--      box. Stage-1 simplification of the Tayoun 2018 tree: NMD-triggering → Very Strong; NMD-escaping → Strong. The
--      full tree also weighs %-protein-lost / critical-domain / last-exon / single-exon-gene — extend here as SQL.
-- The gate is disease-aware and dup-safe: EXISTS a gene-disease context where LoF is an established mechanism
-- (lof_mechanism AND HI3). One PVS1 row max, even if the gene has several disease curations.
SELECT s.variant_key,'PVS1_VeryStrong','P',8,'null ('||s.consequence||'), NMD-triggering, LoF is a disease mechanism (HI3)'
FROM annotations_spec s
-- start_lost excluded: Tayoun's start-loss decision path is separate, not generic NMD-PVS1
WHERE s.consequence IN ('stop_gained','frameshift','splice_donor','splice_acceptor')
  AND s.nmd_escaping = 0
  AND (EXISTS (SELECT 1 FROM gene_curation gc WHERE gc.gene = s.gene AND gc.lof_mechanism)  -- LoF a disease mechanism (ClinGen)
       OR (s.loeuf < 0.35 AND (SELECT enabled FROM pvs1_constraint_config)))  -- OPT-IN: LoF-intolerant by gnomAD LOEUF (over-calls; off by default)
UNION ALL
SELECT s.variant_key,'PVS1_Strong','P',4,'null ('||s.consequence||'), NMD-escaping (downgraded), LoF mechanism (HI3)'
FROM annotations_spec s
-- start_lost excluded: Tayoun's start-loss decision path is separate, not generic NMD-PVS1
WHERE s.consequence IN ('stop_gained','frameshift','splice_donor','splice_acceptor')
  AND s.nmd_escaping = 1
  AND (EXISTS (SELECT 1 FROM gene_curation gc WHERE gc.gene = s.gene AND gc.lof_mechanism)  -- LoF a disease mechanism (ClinGen)
       OR (s.loeuf < 0.35 AND (SELECT enabled FROM pvs1_constraint_config)))  -- OPT-IN: LoF-intolerant by gnomAD LOEUF (over-calls; off by default);
UNION ALL
-- CADD (OPT-IN, cadd_config): a SUPPORTING-only computational FALLBACK for the tail with NO other evidence. CADD
-- PHRED is not ClinGen-calibrated like REVEL, so it NEVER stacks on a calibrated tool or a known answer — it fires
-- ONLY when the variant is (a) not REVEL-scored (a missense with REVEL is covered), (b) not SpliceAI-scored (a
-- splice-concept variant), (c) not in ClinVar, and (d) not a PVS1/PM4 mechanism consequence (a null/in-frame call
-- already carries its own criterion — PP3 there would double-count). PHRED>=20 -> PP3(+1); <=10 -> BP4(-1);
-- the (10,20) zone abstains. Supporting strength only; gated behind cadd_config (default off).
SELECT variant_key,'PP3','P',1,'CADD PHRED '||cadd_phred||' >= 20 (fallback: no REVEL/SpliceAI/ClinVar)'
FROM annotations_spec
WHERE (SELECT enabled FROM cadd_config) AND clinvar_classification IS NULL
  AND cadd_phred IS NOT NULL AND cadd_phred >= 20
  AND NOT (consequence = 'missense' AND revel IS NOT NULL)
  AND NOT (consequence IN ('splice_region','splice_donor','splice_acceptor','intron','synonymous') AND spliceai IS NOT NULL)
  AND consequence NOT IN ('stop_gained','frameshift','splice_donor','splice_acceptor','inframe_deletion','inframe_insertion','stop_lost')
UNION ALL
SELECT variant_key,'BP4','B',-1,'CADD PHRED '||cadd_phred||' <= 10 (fallback: no REVEL/SpliceAI/ClinVar)'
FROM annotations_spec
WHERE (SELECT enabled FROM cadd_config) AND clinvar_classification IS NULL
  AND cadd_phred IS NOT NULL AND cadd_phred <= 10
  AND NOT (consequence = 'missense' AND revel IS NOT NULL)
  AND NOT (consequence IN ('splice_region','splice_donor','splice_acceptor','intron','synonymous') AND spliceai IS NOT NULL)
  AND consequence NOT IN ('stop_gained','frameshift','splice_donor','splice_acceptor','inframe_deletion','inframe_insertion','stop_lost')
