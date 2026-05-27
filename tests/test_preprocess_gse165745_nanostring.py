from __future__ import annotations

from scripts.preprocess.preprocess_gse165745_nanostring import parse_gse165745_series_matrix


def test_parse_gse165745_series_matrix_maps_labels_and_gene_aliases():
    text = "\n".join(
        [
            '!Sample_title\t"R001: Responder_001"\t"NR001: Nonresponder_001"',
            '!Sample_geo_accession\t"GSM1"\t"GSM2"',
            '!Sample_source_name_ch1\t"Melanoma"\t"Melanoma"',
            '!Sample_characteristics_ch1\t"phenotype: Responder"\t"phenotype: Nonresponder"',
            '!Sample_characteristics_ch1\t"age: 82"\t"age: 55"',
            '!Sample_characteristics_ch1\t"disease site: Skin metastasis"\t"disease site: Lung metastasis"',
            '!Sample_characteristics_ch1\t"gender: F"\t"gender: M"',
            '!Sample_characteristics_ch1\t"cell type: Melanoma"\t"cell type: Melanoma"',
            '!Sample_description\t"taken prior to treatment"\t"taken prior to treatment"',
            '!Sample_data_processing\t"nSolver normalized"\t"nSolver normalized"',
            "!series_matrix_table_begin",
            '"ID_REF"\t"GSM1"\t"GSM2"',
            '"CD8a"\t3\t7',
            '"Arg1"\t1\t3',
            '"IFNG"\t15\t31',
            "!series_matrix_table_end",
        ]
    )

    expression, metadata = parse_gse165745_series_matrix(text)

    assert list(expression.index) == ["GSM1", "GSM2"]
    assert {"CD8A", "ARG1", "IFNG"}.issubset(expression.columns)
    assert metadata.loc[0, "label"] == 1
    assert metadata.loc[1, "label"] == 0
    assert metadata.loc[0, "label_status"] == "source_binary_responder_nonresponder_not_recist"
