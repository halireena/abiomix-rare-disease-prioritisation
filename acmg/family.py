"""Family abstraction: derive pedigree, carriage, and segregation from the canonical `sample_call`
table — from DATA, never hardcoded per case.

Two engineering points this handles carefully (and a common pitfall it avoids):
  1. de novo is QC-GATED. 'Absent in a parent' only counts when that parent has a CALLED genotype with
     adequate GQ/DP — a `./.` no-call or a low-coverage site is NOT evidence of absence. (The common
     bug is treating missing == homozygous-reference, which invents de novo events everywhere.)
  2. This bundle is UNPHASED (FORMAT_PS is all NULL), so compound-het cannot use phase. We infer it by
     PARENTAL ORIGIN from a trio: two hets in the same gene, one seen only in the father and one only
     in the mother, are in trans. With phased data you would read cis/trans from the phase set instead.

Coordinates are GRCh37 throughout (annotate with VEP GRCh37) so proband and parents compare directly —
no liftover.
"""
from __future__ import annotations
import duckdb

# The UNIVERSAL pedigree path is a standard PLINK PED (paternal_id / maternal_id / affected columns) — see
# acmg.ingest.ingest_pedigree; that is format- and vocabulary-independent. The vocabulary below is only for
# sources (like this challenge's curated fields) that carry free-text relationship labels instead of parent
# columns, and it is exposed as PARAMETERS, not hardcoded — pass your own for a differently-labelled source.
PROBAND_STATUS_PRIORITY = ("proband", "single_member_case", "unclear_multi_affected", "unclear")
PROBAND_LABELS = ("proband_child", "proband_son", "affected_son", "child", "affected_child")
MOTHER_LABELS = ("mother",)
FATHER_LABELS = ("father", "affected_father")
AMBIGUOUS_STATUS = ("unclear", "unclear_multi_affected")


def _inlist(xs) -> str:
    return ", ".join("'" + str(x).lower() + "'" for x in xs)


def build_pedigree(con: duckdb.DuckDBPyConnection, *,
                   proband_status_priority=PROBAND_STATUS_PRIORITY, proband_labels=PROBAND_LABELS,
                   mother_labels=MOTHER_LABELS, father_labels=FATHER_LABELS,
                   ambiguous_status=AMBIGUOUS_STATUS) -> None:
    """ONE row per family, resolving proband/mother/father from curated STATUS + RELATIONSHIP labels. The
    vocabulary is PARAMETERISED (defaults = this challenge's curated values) — pass your own for another
    source. **For standard pedigrees prefer acmg.ingest.ingest_pedigree (a PLINK PED): that is the universal,
    label-free path.** Proband priority: explicit 'proband' -> singleton -> a multiplex/unclear affected
    member (flagged `proband_ambiguous`). Robust to singletons + messy multiplex families (no dup rows,
    none dropped)."""
    proband_coalesce = ",\n                ".join(
        [f"max(case_id) FILTER (WHERE proband_status = '{s}')" for s in proband_status_priority]
        + [f"max(case_id) FILTER (WHERE lower(relationship) IN ({_inlist(proband_labels)}))"]
    )
    con.execute(f"""
        CREATE OR REPLACE TABLE pedigree AS
        WITH members AS (
            SELECT DISTINCT case_id, family_id, member_id, proband_status, relationship, family_size
            FROM sample_call
        )
        SELECT
            family_id,
            COALESCE(
                {proband_coalesce}
            ) AS proband_id,
            max(case_id) FILTER (WHERE lower(relationship) IN ({_inlist(mother_labels)})) AS mother_id,
            max(case_id) FILTER (WHERE lower(relationship) IN ({_inlist(father_labels)})) AS father_id,
            any_value(family_size) AS family_size,
            (bool_or(proband_status IN ({_inlist(ambiguous_status)}))
             AND NOT bool_or(proband_status = 'proband')) AS proband_ambiguous
        FROM members
        GROUP BY family_id
    """)


def build_carriage(con: duckdb.DuckDBPyConnection, *, min_gq: int = 20, min_dp: int = 10) -> None:
    """`carried` = confident non-reference calls. `called` = any confident call (ref or alt) — the
    denominator for 'confidently absent'. QC thresholds are explicit and tunable."""
    con.execute(f"""
        CREATE OR REPLACE VIEW called AS
        SELECT * FROM sample_call
        WHERE gt NOT IN ('./.', '.|.')
          AND (gq IS NULL OR gq >= {min_gq})
          AND (dp IS NULL OR dp >= {min_dp})
    """)
    con.execute("""
        CREATE OR REPLACE VIEW carried AS
        SELECT * FROM called WHERE gt NOT IN ('0/0', '0|0')
    """)


def de_novo_candidates(con: duckdb.DuckDBPyConnection) -> duckdb.DuckDBPyRelation:
    """Proband-carried variants where BOTH parents have a confident CALLED genotype and NEITHER carries
    it. Requiring called parents is the QC gate that avoids false de novo from no-calls / low coverage."""
    return con.sql("""
        SELECT p.proband_id AS case_id, cv.variant_key, cv.gene
        FROM pedigree p
        JOIN carried cv ON cv.case_id = p.proband_id
        WHERE p.mother_id IS NOT NULL AND p.father_id IS NOT NULL
          AND EXISTS (SELECT 1 FROM called m WHERE m.case_id = p.mother_id AND m.variant_key = cv.variant_key)
          AND EXISTS (SELECT 1 FROM called f WHERE f.case_id = p.father_id AND f.variant_key = cv.variant_key)
          AND NOT EXISTS (SELECT 1 FROM carried m WHERE m.case_id = p.mother_id AND m.variant_key = cv.variant_key)
          AND NOT EXISTS (SELECT 1 FROM carried f WHERE f.case_id = p.father_id AND f.variant_key = cv.variant_key)
    """)


def affected_shared_candidates(con: duckdb.DuckDBPyConnection, *,
                               affected_labels=("affected_son", "affected_child", "affected_father", "affected_daughter"),
                               affected_status=("proband", "unclear_multi_affected"),
                               min_affected: int = 2) -> duckdb.DuckDBPyRelation:
    """MULTIPLEX / multi-affected-sibling segregation: variants carried HOMOZYGOUS by >= `min_affected`
    AFFECTED members of a family (recessive; X-linked shows as '1/1' in hemizygous males). A genotype-driven
    reanalysis signal that needs NO phenotype — the segregation itself is the evidence. Requires build_carriage.

    IMPORTANT (ancestry-aware): a variant common in the patient's population (e.g. G6PD in Middle Eastern /
    Mediterranean populations) will segregate in every affected sib yet be an EXPECTED / incidental finding,
    not causal. Always filter these candidates by MAX-POPULATION filtering AF (annotations_spec), not global
    gnomAD — see the kernel's frequency handling. This view finds the segregating set; rarity + pathogenicity
    (the kernel) decide which one matters."""
    who = ", ".join("'" + x.lower() + "'" for x in affected_labels)
    stat = ", ".join("'" + x + "'" for x in affected_status)
    return con.sql(f"""
        WITH affected AS (
            SELECT DISTINCT case_id, family_id FROM sample_call
            WHERE lower(relationship) IN ({who}) OR proband_status IN ({stat})
        )
        SELECT c.family_id, c.variant_key, any_value(c.gene) AS gene,
               any_value(c.chrom) AS chrom, any_value(c.pos) AS pos,
               count(DISTINCT c.case_id) AS n_affected_hom
        FROM carried c JOIN affected a USING (case_id, family_id)
        WHERE c.gt IN ('1/1', '1|1')
        GROUP BY c.family_id, c.variant_key
        HAVING count(DISTINCT c.case_id) >= {min_affected}
    """)


def compound_het_candidates(con: duckdb.DuckDBPyConnection) -> duckdb.DuckDBPyRelation:
    """Trio, UNPHASED: two heterozygous variants in the SAME gene in the proband, one inherited only
    from the father and the other only from the mother -> in trans (biallelic). Requires a `gene` on
    each call (join the VEP annotation onto sample_call first). With phased data, use the phase set."""
    # SCOPED PER FAMILY: carry family_id + this proband's own father_id/mother_id, and check parental origin
    # against THOSE ids only. Joining carried->pedigree on any father_id/mother_id matched any family's parent,
    # inventing cross-family in-trans pairs (and masking true ones). Correlate to the same pedigree row.
    return con.sql("""
        WITH prob AS (
            SELECT p.family_id, p.father_id, p.mother_id, c.variant_key, c.gene
            FROM carried c JOIN pedigree p ON c.case_id = p.proband_id
            WHERE p.mother_id IS NOT NULL AND p.father_id IS NOT NULL
              AND c.gt IN ('0/1', '1/0', '0|1', '1|0') AND c.gene IS NOT NULL
        ),
        origin AS (
            SELECT prob.family_id, prob.variant_key, prob.gene,
                   EXISTS (SELECT 1 FROM carried f WHERE f.case_id = prob.father_id AND f.variant_key = prob.variant_key) AS from_father,
                   EXISTS (SELECT 1 FROM carried m WHERE m.case_id = prob.mother_id AND m.variant_key = prob.variant_key) AS from_mother
            FROM prob
        )
        SELECT a.family_id, a.gene, a.variant_key AS variant_a, b.variant_key AS variant_b
        FROM origin a JOIN origin b ON a.family_id = b.family_id AND a.gene = b.gene AND a.variant_key < b.variant_key
        WHERE a.from_father AND NOT a.from_mother
          AND b.from_mother AND NOT b.from_father
    """)
