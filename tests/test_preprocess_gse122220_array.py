from __future__ import annotations

from scripts.preprocess.preprocess_gse122220_array import (
    parse_gpl10558_annotation,
    parse_gse122220_series_matrix,
)


def test_parse_gse122220_series_matrix_maps_array_probes_and_recist_labels():
    annot_text = "\n".join(
        [
            "!platform_table_begin",
            "ID\tGene symbol",
            "ILMN_1\tMAP4K1",
            "ILMN_2\tTBX3",
            "ILMN_3\tTBX3",
            "ILMN_4\tAXL",
            "!platform_table_end",
        ]
    )
    matrix_text = "\n".join(
        [
            '!Sample_title\t"108"\t"19"\t"38"',
            '!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"',
            '!Sample_source_name_ch1\t"Melanoma tumor biopsy"\t"Melanoma tumor biopsy"\t"Melanoma tumor biopsy"',
            '!Sample_characteristics_ch1\t"timepoint: Pre-treatment"\t"timepoint: Pre-treatment"\t"timepoint: Pre-treatment"',
            '!Sample_characteristics_ch1\t"treatment: Ipilimumab+PD1"\t"treatment: PD1"\t"treatment: PD1"',
            '!Sample_characteristics_ch1\t"previous ipilumimab: No"\t"previous ipilumimab: Yes"\t"previous ipilumimab: Yes"',
            '!Sample_characteristics_ch1\t"response (pd = progressive disease; sd = stable disease; pr = partial response; cr = complete response): PR"\t"response (pd = progressive disease; sd = stable disease; pr = partial response; cr = complete response): PD"\t"response (pd = progressive disease; sd = stable disease; pr = partial response; cr = complete response): SD"',
            '!Sample_characteristics_ch1\t"age: 31"\t"age: 60"\t"age: 52"',
            '!Sample_characteristics_ch1\t"gender: F"\t"gender: F"\t"gender: M"',
            "!series_matrix_table_begin",
            '"ID_REF"\t"GSM1"\t"GSM2"\t"GSM3"',
            '"ILMN_1"\t1023\t255\t511',
            '"ILMN_2"\t15\t31\t63',
            '"ILMN_3"\t31\t63\t127',
            '"ILMN_4"\t255\t511\t1023',
            "!series_matrix_table_end",
        ]
    )

    annotation = parse_gpl10558_annotation(annot_text)
    expression, metadata = parse_gse122220_series_matrix(matrix_text, annotation)

    assert list(expression.index) == ["GSM1", "GSM2", "GSM3"]
    assert {"MAP4K1", "TBX3", "AXL"}.issubset(expression.columns)
    assert expression.loc["GSM1", "TBX3"] == expression.loc["GSM1", ["TBX3"]].iloc[0]
    assert metadata.loc[0, "response_raw"] == "PR"
    assert metadata.loc[1, "response_raw"] == "PD"
    assert metadata.loc[2, "response_raw"] == "SD"
    assert metadata.loc[0, "label"] == 1
    assert metadata.loc[2, "label"] == 0
    assert metadata.loc[0, "therapy"] == "anti-CTLA-4 plus anti-PD-1 combination"
    assert metadata.loc[1, "therapy"] == "anti-PD-1 monotherapy"
