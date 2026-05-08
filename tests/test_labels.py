import pandas as pd

from econiche.labels import harmonize_response, harmonize_metadata


def test_primary_recist_maps_cr_pr_to_response_and_sd_pd_to_nonresponse():
    assert harmonize_response("CR", endpoint="primary_recist") == 0
    assert harmonize_response("partial response", endpoint="primary_recist") == 0
    assert harmonize_response("SD", endpoint="primary_recist") == 1
    assert harmonize_response("progressive disease", endpoint="primary_recist") == 1


def test_strict_recist_excludes_stable_disease():
    assert harmonize_response("SD", endpoint="strict_recist") is None
    assert harmonize_response("PD", endpoint="strict_recist") == 1


def test_clinical_benefit_maps_dcb_ndb():
    assert harmonize_response("DCB", endpoint="clinical_benefit") == 0
    assert harmonize_response("NDB", endpoint="clinical_benefit") == 1


def test_harmonize_metadata_adds_label_column_and_keeps_unknowns_as_missing():
    meta = pd.DataFrame({"response_raw": ["PR", "PD", "not reported"]})

    out = harmonize_metadata(meta, endpoint="primary_recist")

    assert out["label"].tolist()[:2] == [0.0, 1.0]
    assert pd.isna(out.loc[2, "label"])
