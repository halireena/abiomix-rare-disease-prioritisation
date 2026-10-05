from acmg.vep_rest import to_annotations


def _rec(colocated):
    return {"input": "1 1000 . A G", "seq_region_name": "1", "start": 1000,
            "most_severe_consequence": "missense_variant",
            "transcript_consequences": [{"gene_symbol": "GENE1", "canonical": 1}],
            "colocated_variants": colocated}


def test_gnomad_af_read_from_nested_frequencies_and_maxed():
    rec = _rec([
        {"id": "COSV123"},  # somatic entry first, no frequencies
        {"id": "rs1", "frequencies": {"G": {"gnomade": 0.001, "gnomadg": 0.003, "af": 0.5}}},
    ])
    assert to_annotations([rec])["filtering_af"].iloc[0] == 0.003


def test_gnomad_af_ignores_other_alleles():
    rec = _rec([{"id": "rs1", "frequencies": {"T": {"gnomade": 0.2}}}])
    assert to_annotations([rec])["filtering_af"].isna().iloc[0]


def test_gnomad_af_legacy_top_level_keys():
    rec = _rec([{"id": "rs1", "gnomADe_AF": 0.01}])
    assert to_annotations([rec])["filtering_af"].iloc[0] == 0.01
