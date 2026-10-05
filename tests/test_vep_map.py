"""VEP -> kernel `annotations` mapping: consequence picking, variant kind, numeric parsing and the output schema.
Pure pandas, no network. Run: `PYTHONPATH=. python3 tests/test_vep_map.py`."""
import math
import pandas as pd
from acmg import vep_map
from acmg.vep_map import REQUIRED_COLS, most_severe_kernel_consequence, variant_kind, vep_to_annotations


def test_variant_kind():
    assert variant_kind("A", "G") == "snv"
    assert variant_kind("AC", "A") == "indel"
    assert variant_kind("A", "AT") == "indel"
    assert variant_kind("-", "A") == "indel"   # VCF-less dash notation is never an SNV
    assert variant_kind("A", "-") == "indel"


def test_most_severe_consequence_picks_by_severity_not_order():
    assert most_severe_kernel_consequence("missense_variant,stop_gained") == "stop_gained"
    assert most_severe_kernel_consequence("splice_region_variant&synonymous_variant") == "splice_region"
    assert most_severe_kernel_consequence("intron_variant, missense_variant") == "missense"  # whitespace tolerated


def test_most_severe_consequence_abstains_on_unmapped_or_missing():
    assert most_severe_kernel_consequence("intergenic_variant") is None
    assert most_severe_kernel_consequence("") is None
    assert most_severe_kernel_consequence(None) is None
    assert most_severe_kernel_consequence(float("nan")) is None


def test_num_treats_vep_placeholders_as_missing():
    assert vep_map._num("0.5") == 0.5
    assert vep_map._num(0) == 0.0
    for missing in (None, "", ".", "-", float("nan"), "not-a-number"):
        assert vep_map._num(missing) is None


def test_vep_to_annotations_schema_and_values():
    vep_df = pd.DataFrame([
        dict(CHROM="chr17", POS=43093464, REF="AC", ALT="A", SYMBOL="BRCA1",
             Consequence="frameshift_variant", REVEL=None, SpliceAI_pred_DS_max=None, gnomADe_AF="5e-06"),
        dict(CHROM="1", POS=100000, REF="A", ALT="G", SYMBOL="SCN1A",
             Consequence="missense_variant&splice_region_variant", REVEL="0.95", SpliceAI_pred_DS_max=".",
             gnomADe_AF="."),
    ])
    ann = vep_to_annotations(vep_df)
    assert list(ann.columns) == REQUIRED_COLS
    assert list(ann.variant_key) == ["17-43093464-AC-A", "1-100000-A-G"]   # 'chr' stripped
    assert list(ann.consequence) == ["frameshift", "missense"]
    assert list(ann.variant_kind) == ["indel", "snv"]
    assert ann.loc[0, "filtering_af"] == 5e-06
    assert math.isnan(ann.loc[1, "filtering_af"])   # '.' -> missing, so PM2 abstains downstream
    assert ann.loc[1, "revel"] == 0.95
    # evidence columns from other sources start empty: the kernel abstains rather than guesses
    for col in ("gnomad_mis_z", "clinvar_same_aa", "clinvar_same_codon_lp", "nmd_escaping"):
        assert ann[col].isna().all(), col


def test_vep_to_annotations_custom_column_map():
    vep_df = pd.DataFrame([dict(chr="2", start=5, ref="C", alt="T", gene_name="XYZ",
                                csq="synonymous_variant", af=0.2)])
    ann = vep_to_annotations(vep_df, cols={"chrom": "chr", "pos": "start", "ref": "ref", "alt": "alt",
                                           "gene": "gene_name", "consequence": "csq", "gnomad_af": "af"})
    row = ann.iloc[0]
    assert (row.variant_key, row.gene, row.consequence, row.filtering_af) == ("2-5-C-T", "XYZ", "synonymous", 0.2)
    assert pd.isna(row.revel) and pd.isna(row.spliceai)   # absent columns -> None, not a KeyError


def test_vep_to_annotations_missing_coordinates_says_which_column():
    # a VEP table whose position column is named differently: name it and point at `cols=`, not int(None)
    vep_df = pd.DataFrame([dict(CHROM="1", Start=5, REF="C", ALT="T", SYMBOL="XYZ")])
    try:
        vep_to_annotations(vep_df)
    except ValueError as e:
        assert "POS" in str(e) and "cols=" in str(e)
    else:
        raise AssertionError("expected a ValueError naming the missing column")
    assert vep_to_annotations(pd.DataFrame()).empty   # an empty input is still just an empty frame


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name)
    print("all passed")
