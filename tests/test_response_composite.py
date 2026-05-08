import numpy as np
import pandas as pd

from econiche_opt.model.response_composite import (
    build_candidate_scores,
    build_signature_features,
    response_labels_from_nonresponse,
    run_nested_response_composite,
)


def test_response_labels_invert_nonresponse_labels():
    labels = pd.Series([0, 1, 1, 0], index=list("abcd"))
    response = response_labels_from_nonresponse(labels)
    assert response.tolist() == [1, 0, 0, 1]
    assert response.index.tolist() == labels.index.tolist()


def test_signature_candidates_include_immune_top5():
    X = pd.DataFrame(
        {
            "IFNG": [1.0, 2.0, 3.0],
            "CXCL9": [1.0, 2.0, 4.0],
            "CXCL10": [1.0, 2.0, 5.0],
            "STAT1": [1.0, 2.0, 6.0],
            "PDCD1LG2": [1.0, 3.0, 5.0],
            "CD8A": [1.0, 3.0, 6.0],
            "GZMA": [1.0, 4.0, 7.0],
            "GZMB": [1.0, 4.0, 8.0],
            "IDO1": [1.0, 5.0, 9.0],
            "PDCD1": [1.0, 3.0, 6.0],
            "LAG3": [1.0, 3.0, 6.0],
            "HAVCR2": [1.0, 3.0, 6.0],
            "TOX": [1.0, 3.0, 6.0],
        },
        index=["s1", "s2", "s3"],
    )
    features = build_signature_features(X)
    candidates = build_candidate_scores(features)
    assert "immune_top5" in candidates
    assert "ifn_core_pdcd1lg2_weighted" in candidates
    assert np.isfinite(candidates["immune_top5"]).all()


def test_nested_response_composite_outputs_lodo_predictions():
    X_by_cohort = {}
    y_by_cohort = {}
    metadata_by_cohort = {}
    for idx, cohort in enumerate(["c1", "c2", "c3"]):
        samples = [f"{cohort}_s{i}" for i in range(6)]
        response = np.array([0, 0, 0, 1, 1, 1])
        signal = response + idx * 0.05
        X_by_cohort[cohort] = pd.DataFrame(
            {
                "IFNG": signal,
                "CXCL9": signal,
                "CXCL10": signal,
                "STAT1": signal,
                "PDCD1LG2": signal,
                "CD8A": signal,
                "GZMA": signal,
                "GZMB": signal,
                "IDO1": signal,
                "PDCD1": signal,
                "LAG3": signal,
                "HAVCR2": signal,
                "TOX": signal,
            },
            index=samples,
        )
        y_by_cohort[cohort] = pd.Series(1 - response, index=samples)
        metadata_by_cohort[cohort] = pd.DataFrame(
            {
                "sample_id": samples,
                "patient_id": samples,
                "cohort": cohort,
                "label": 1 - response,
            },
            index=samples,
        )
    result = run_nested_response_composite(X_by_cohort, y_by_cohort, metadata_by_cohort, preferred_tolerance=1.0)
    assert set(result.metrics["cohort"]) == {"c1", "c2", "c3"}
    assert len(result.predictions) == 18
    assert {"response_probability", "true_response_label", "selected_model"}.issubset(result.predictions.columns)
