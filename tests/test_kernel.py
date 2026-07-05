"""Kernel behavior tests. Run: `python -m pytest` (or `python tests/test_kernel.py`)."""
import pandas as pd
from acmg import classify, vep_to_annotations

GENE_CURATION = pd.DataFrame([
    dict(gene="LOFGENE", disease_id="MONDO:1", mode_of_inheritance="AD", hi_score=3,
         lof_mechanism=True, gene_disease_validity="Definitive", source="test"),
])


def _classify(rows, **kw):
    df = classify(pd.DataFrame(rows), **kw)
    return {r.variant_key: r for r in df.itertuples()}


def test_lof_in_hi3_gene_gets_pvs1():
    r = _classify([dict(variant_key="1-1-A-AC", gene="LOFGENE", consequence="frameshift",
                        variant_kind="indel", filtering_af=1e-6, nmd_escaping=0)],
                  gene_curation=GENE_CURATION)
    assert "PVS1_VeryStrong" in r["1-1-A-AC"].criteria
    assert r["1-1-A-AC"].acmg_class == "Likely Pathogenic"  # PVS1(8)+PM2(1)=9


def test_lof_in_uncurated_gene_abstains_on_pvs1():
    r = _classify([dict(variant_key="1-2-A-AC", gene="RANDOM", consequence="frameshift",
                        variant_kind="indel", filtering_af=1e-6, nmd_escaping=0)])
    assert "PVS1" not in r["1-2-A-AC"].criteria  # no PVS1 without an established LoF mechanism


def test_common_variant_is_benign_by_ba1():
    r = _classify([dict(variant_key="1-3-A-G", gene="X", consequence="missense",
                        variant_kind="snv", filtering_af=0.2, revel=0.3)])
    assert r["1-3-A-G"].acmg_class == "Benign" and "BA1" in r["1-3-A-G"].criteria


def test_pm2_abstains_on_unknown_frequency():
    r = _classify([dict(variant_key="1-4-C-T", gene="X", consequence="missense",
                        variant_kind="snv", filtering_af=None, revel=0.7)])
    assert "PM2_Supporting" not in r["1-4-C-T"].criteria  # no data != rare


def test_cnv_is_not_evaluated_by_snv_indel_kernel():
    r = _classify([dict(variant_key="1-5-DEL", gene="X", consequence=None,
                        variant_kind="cnv", filtering_af=1e-6)])
    assert r["1-5-DEL"].acmg_class.startswith("Not evaluated")


def test_vep_mapping_shape():
    ann = vep_to_annotations(pd.DataFrame([
        dict(CHROM="chr2", POS=100, REF="A", ALT="G", SYMBOL="TTN",
             Consequence="missense_variant&splice_region_variant", REVEL=0.4,
             SpliceAI_pred_DS_max=0.1, gnomADe_AF=1e-4)]))
    row = ann.iloc[0]
    assert row.variant_key == "2-100-A-G"
    assert row.consequence == "missense"   # most-severe pick among mapped SO terms
    assert row.variant_kind == "snv"
    assert abs(row.filtering_af - 1e-4) < 1e-12


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name)
    print("all passed")


def test_classify_accepts_gene_curation_with_ts_score():
    """A REAL ClinGen gene_curation carries ts_score (triplosensitivity, for CNV DUP dosage). The kernel doesn't
    use it, but its spec table must include the column so `INSERT ... BY NAME` doesn't reject the real feed —
    the annotate_cohort pipeline passes exactly this. Regression for the schema-mismatch break."""
    gc = pd.DataFrame([
        dict(gene="LOFGENE", disease_id="MONDO:1", mode_of_inheritance="AD", hi_score=3, ts_score=None,
             lof_mechanism=True, gene_disease_validity="Definitive", source="test"),
    ])
    ann = pd.DataFrame([dict(variant_key="1-1-A-G", gene="LOFGENE", consequence="missense_variant", filtering_af=0.0)])
    out = classify(ann, gene_curation=gc)   # must not raise on the extra ts_score column
    assert len(out) == 1


def test_cadd_supporting_fallback_only_on_the_uncovered_tail():
    """CADD is an OPT-IN supporting fallback: it fires PP3/BP4 (+/-1) ONLY where no calibrated tool and no ClinVar
    speaks. Off by default; never stacks on REVEL/SpliceAI/ClinVar/PVS1-PM4. All rows have revel=None, so any PP3
    here can only be the CADD PP3 (a REVEL PP3 would need a revel score)."""
    rows = [
        # (a) missense WITHOUT revel, high CADD, not in ClinVar -> CADD PP3 fires (the tail CADD is for)
        dict(variant_key="1-1-A-G", gene="G", consequence="missense", revel=None, cadd_phred=28.0),
        # (c) in ClinVar -> CADD gated out
        dict(variant_key="3-1-C-T", gene="G", consequence="missense", revel=None, cadd_phred=30.0, clinvar_classification="Pathogenic"),
        # (d) stop_gained (PVS1/mechanism consequence) -> CADD excluded (no double-count)
        dict(variant_key="4-1-C-T", gene="G", consequence="stop_gained", revel=None, cadd_phred=35.0),
        # (e) CADD in the indeterminate zone (10,20) -> abstain
        dict(variant_key="5-1-C-T", gene="G", consequence="missense", revel=None, cadd_phred=15.0),
        # (f) low CADD on the tail -> BP4 fires
        dict(variant_key="6-1-C-T", gene="G", consequence="missense", revel=None, cadd_phred=5.0),
    ]
    ann = pd.DataFrame(rows)
    off = classify(ann).set_index("variant_key")
    assert "PP3" not in (off.loc["1-1-A-G", "criteria"] or "")          # default off -> no CADD anywhere
    on = classify(ann, cadd_supporting=True).set_index("variant_key")
    assert "PP3" in on.loc["1-1-A-G", "criteria"]                        # (a) CADD PP3 fires (supporting, +1)
    assert on.loc["1-1-A-G", "total_points"] == 1
    assert "BP4" in on.loc["6-1-C-T", "criteria"]                        # (f) CADD BP4 fires (-1)
    for vk in ("3-1-C-T", "4-1-C-T", "5-1-C-T"):                         # gated / excluded / indeterminate
        crit = on.loc[vk, "criteria"] or ""
        assert "PP3" not in crit and "BP4" not in crit
