from __future__ import annotations

import pandas as pd

from scripts.preprocess.fetch_cbioportal_melanoma_icb import CbioStudySpec, build_metadata, requested_gene_symbols, response_to_raw


def test_cbioportal_response_mapping_normalizes_recist_terms():
    assert response_to_raw("Complete Response") == "CR"
    assert response_to_raw("Partial Response") == "PR"
    assert response_to_raw("Stable Disease") == "SD"
    assert response_to_raw("Progressive Disease") == "PD"
    assert response_to_raw("Mixed Response") == "MR"


def test_cbioportal_requested_genes_include_locked_rescue_head():
    genes = set(requested_gene_symbols())

    assert {"MAP4K1", "TBX3", "AXL", "PLA2G2D", "PIK3CD"}.issubset(genes)


def test_cbioportal_metadata_filters_pre_samples_and_joins_liu_patient_attrs():
    spec = CbioStudySpec(
        study_id="mel_iatlas_liu_2019",
        molecular_profile_id="profile",
        sample_list_id="sample_list",
        output_cohort="CBIO_IATLAS_LIU_2019_PRE",
        response_attribute="RESPONSE",
        treatment_attribute="ICI_RX",
        sample_treatment_attribute="SAMPLE_TREATMENT",
        response_source_level="sample",
        treatment_source_level="patient",
    )
    expr = pd.DataFrame({"IFNG": [1.0, 2.0]}, index=["Liu_Sample1", "Liu_Sample2"])
    sample_clinical = pd.DataFrame(
        {
            "RESPONSE": ["Partial Response", "Progressive Disease"],
            "SAMPLE_TREATMENT": ["Pre", "On"],
        },
        index=["Liu_Sample1", "Liu_Sample2"],
    )
    patient_clinical = pd.DataFrame(
        {"ICI_RX": ["Pembrolizumab", "Nivolumab"]},
        index=["Patient1", "Patient2"],
    )

    metadata = build_metadata(spec, expr, sample_clinical, patient_clinical)

    assert metadata["sample_id"].tolist() == ["Liu_Sample1"]
    assert metadata["patient_id"].tolist() == ["Patient1"]
    assert metadata["response_raw"].tolist() == ["PR"]
    assert metadata["label"].tolist() == [1]
    assert metadata["treatment"].tolist() == ["Pembrolizumab"]
