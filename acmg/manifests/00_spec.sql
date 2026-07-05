-- The ACMG rule SPECIFICATION + gene/variant KNOWLEDGE as DATA, so a ClinGen VCEP / national panel (ACGS) /
-- disease-specific exception is a ROW, not a code fork, and the frequency/gene sources are pluggable.

-- (1) THRESHOLDS per gene: `gene IS NULL` is the general default; a gene row OVERRIDES it. A VCEP exception is an
-- INSERT. Stage 1 parameterizes the most-overridden knobs (BA1/BS1/PM2 cutoffs); the pattern extends to REVEL/PS4.
CREATE OR REPLACE TABLE acmg_thresholds (panel TEXT, gene TEXT, ba1_maf DOUBLE, bs1_maf DOUBLE, pm2_maf DOUBLE);
INSERT INTO acmg_thresholds VALUES
  ('ACMG-general',      NULL,   0.05, 0.01,   0.005),    -- the default spec: applies to every gene
  ('ClinGen-VCEP:PAH',  'PAH',  0.05, 0.02,   0.002),    -- illustrative VCEP override (PAH is NOT in the sample)
  ('ACGS-UK:MYH7',      'MYH7', 0.05, 0.0004, 0.00004);  -- illustrative disease/national panel override

-- (2) FREQUENCY is a PROVIDER-AGNOSTIC concept, not "gnomad_af". The rule input is the ACMG filtering allele
-- frequency, resolved as the MAX over populations (ancestry-deconvoluted) with AN >= 2000 — because a variant common
-- in ONE ancestry (rs334 in African populations) must not be diluted by a cohort-wide average. The SOURCE is
-- pluggable: gnomAD per-ancestry / grpmax / FAF, or an internal cohort frequency table — just more rows here.
CREATE OR REPLACE TABLE variant_frequency (variant_key TEXT, source TEXT, population TEXT, af DOUBLE, an BIGINT);
INSERT INTO variant_frequency VALUES
  ('11-5227002-T-A', 'gnomAD-v4', 'afr', 0.06,   20000),   -- rs334 HbS: high in African ancestry (malaria selection)
  ('11-5227002-T-A', 'gnomAD-v4', 'nfe', 0.0001, 50000);   -- rare in non-Finnish European — max-pop is what counts

-- Per-variant resolved thresholds + resolved filtering AF. `filtering_af` = max-population AF from `variant_frequency`
-- (AN>=2000) if present, else the variant's own pre-summarized `filtering_af` annotation column. The kernel reads
-- this generic column, never a provider-specific one.
CREATE OR REPLACE VIEW annotations_spec AS
SELECT a.* EXCLUDE (filtering_af),
       COALESCE(
         (SELECT max(vf.af) FROM variant_frequency vf WHERE vf.variant_key = a.variant_key AND vf.af IS NOT NULL AND vf.an >= 2000),
         a.filtering_af
       ) AS filtering_af,
       COALESCE(g.ba1_maf, d.ba1_maf) AS ba1_maf,
       COALESCE(g.bs1_maf, d.bs1_maf) AS bs1_maf,
       COALESCE(g.pm2_maf, d.pm2_maf) AS pm2_maf
FROM annotations a
LEFT JOIN acmg_thresholds g ON g.gene = a.gene
CROSS JOIN (SELECT ba1_maf, bs1_maf, pm2_maf FROM acmg_thresholds WHERE gene IS NULL) d;

-- (3) GENE-LEVEL CURATION as data (ClinGen Dosage haploinsufficiency + Gene-Disease Validity; harmonized via GenCC
-- across ClinGen / Gene2Phenotype / PanelApp; disease IDs via Monarch/Mondo). PVS1 ENTRY GATE: LoF must be an
-- ESTABLISHED disease mechanism = ClinGen HI score 3 (Tayoun 2018; PMC9035475). HI 0/1/2 = insufficient, 30 =
-- recessive/biallelic, 40 = dosage-insensitive (refutes it). An uncurated gene → PVS1 ABSTAINS. Gene-disease validity
-- (Definitive/Strong/Moderate/Limited/…) CAPS the class (Moderate → <=LP, Limited → <=VUS).
-- KEYED BY (gene, disease_id): the SAME gene can have different validity, LoF mechanism, AND mode of inheritance
-- across diseases (AD in one disorder, AR in another; LoF-haploinsufficiency for one disease, GoF for another). A
-- variant is classified in a disease CONTEXT — the case's suspected disease selects the row. `lof_mechanism` and
-- `hi_score` are per gene-disease, so the PVS1 entry gate and the validity cap are disease-aware, not gene-global.
CREATE OR REPLACE TABLE gene_curation (
  gene TEXT, disease_id TEXT, mode_of_inheritance TEXT, hi_score INTEGER, ts_score INTEGER,
  lof_mechanism BOOLEAN, gene_disease_validity TEXT, source TEXT
);
-- ts_score (ClinGen triplosensitivity) is unused by the SNV kernel; it is here so a real ClinGen gene_curation
-- (which carries it, for CNV DUP dosage) INSERTs BY NAME without a column mismatch.
INSERT INTO gene_curation VALUES
  ('GENEB', 'MONDO:GENEB-disorder', 'AD', 3,  NULL, TRUE,  'Definitive', 'ClinGen-dosage:HI3'),  -- LoF mechanism → PVS1 ok
  ('HBB',   'MONDO:0011382',        'AR', 30, NULL, FALSE, 'Definitive', 'ClinGen-GDV');          -- sickle: AR, not haploinsuff

-- (4) BA1/BS1 EXCEPTION LIST (Ghosh 2018, PMC6188666): established-pathogenic variants COMMON in some population
-- (balancing selection / founder). Frequency-benign must NOT fire for these. rs334 (HBB HbS) is the textbook case.
CREATE OR REPLACE TABLE ba1_exceptions (variant_key TEXT, gene TEXT, note TEXT);
INSERT INTO ba1_exceptions VALUES
  ('11-5227002-T-A', 'HBB', 'rs334 HbS: pathogenic despite high African-population AF (Ghosh 2018 BA1 exception)');

-- (5) PM1 is OPT-IN (default OFF). PM1 (mutational hotspot / critical domain) is the most subjective,
-- VCEP-customised ACMG criterion; it overlaps PP2 (constrained gene) + PP3 (in-silico), risking correlated
-- double-counting, and its status in the emerging points-based ("v4") framework is not a settled global rule.
-- So it does NOT score by default. The `pm1` hotspot flag is STILL computed and surfaced as CONTEXT for the
-- agent/reviewer regardless. Turn on scoring per-run (classify(pm1_enabled=True)) or per-VCEP.
CREATE OR REPLACE TABLE pm1_config (enabled BOOLEAN);
INSERT INTO pm1_config VALUES (FALSE);

-- (6) Constraint-based PVS1 is OPT-IN (default OFF). gnomAD LoF-intolerance (LOEUF) is EVIDENCE that LoF is
-- deleterious, but NOT the ACMG PVS1 requirement of an ESTABLISHED gene-disease LoF mechanism — a constrained
-- essential gene may have no Mendelian disease. So by default PVS1 gates ONLY on ClinGen lof_mechanism.
-- Turn on (classify(pvs1_constraint=True)) to ALSO fire PVS1 on LOEUF<0.35 — a leakage-free proxy for
-- research/validation or coverage of constrained genes lacking ClinGen dosage. It over-calls; use knowingly.
CREATE OR REPLACE TABLE pvs1_constraint_config (enabled BOOLEAN);
INSERT INTO pvs1_constraint_config VALUES (FALSE);

-- CADD as a SUPPORTING-only computational FALLBACK (opt-in, classify(cadd_supporting=True)). CADD PHRED is NOT
-- ClinGen/Pejaver-calibrated the way REVEL is, so it never competes with a calibrated tool: it fires ONLY on the
-- tail with no other evidence — not REVEL-scored (missense w/ REVEL), not SpliceAI-scored (splice concept), not in
-- ClinVar, and not a PVS1/PM4 mechanism consequence. Supporting strength only (+/-1). Off by default (uncalibrated).
CREATE OR REPLACE TABLE cadd_config (enabled BOOLEAN);
INSERT INTO cadd_config VALUES (FALSE);
