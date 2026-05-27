from __future__ import annotations

import pandas as pd

from scripts.analysis.run_no_leakage_ensemble_search import (
    ENSEMBLE_MODEL,
    compare_to_baseline,
    ensemble_specs,
    feature_family_columns,
)


def test_ensemble_specs_cover_feature_families_and_regularization():
    specs = ensemble_specs()
    families = {spec["feature_family"] for spec in specs}
    penalties = {spec["penalty"] for spec in specs}

    assert {"module", "signature", "module_signature_bio", "module_signature_edge", "all"}.issubset(families)
    assert penalties == {"l1", "l2"}


def test_feature_family_columns_selects_expected_prefixes():
    features = pd.DataFrame(
        columns=[
            "module__ifn_t_cell_inflamed",
            "signature__IFNG",
            "signature__IPRES",
            "bio_candidate__bio_prior",
            "edge_feature__edge__ifn_t_cell_inflamed__x__myeloid_suppression",
        ]
    )

    assert feature_family_columns(features, "module") == ["module__ifn_t_cell_inflamed"]
    assert "signature__IFNG" in feature_family_columns(features, "compact")
    assert "bio_candidate__bio_prior" in feature_family_columns(features, "module_signature_bio")
    assert "edge_feature__edge__ifn_t_cell_inflamed__x__myeloid_suppression" in feature_family_columns(features, "all")


def test_compare_to_baseline_reports_positive_delta():
    ensemble = pd.DataFrame(
        {
            "endpoint": ["primary_recist"] * 10,
            "stratum": ["melanoma_core_high_evidence"] * 10,
            "cohort": ["A"] * 5 + ["B"] * 5,
            "sample_id": [f"s{i}" for i in range(10)],
            "model_name": [ENSEMBLE_MODEL] * 10,
            "true_response_label": [0, 1, 1, 0, 1, 0, 1, 0, 1, 0],
            "response_probability": [0.1, 0.9, 0.8, 0.2, 0.7, 0.3, 0.75, 0.25, 0.85, 0.15],
        }
    )
    baseline = pd.DataFrame(
        {
            "endpoint": ["primary_recist"] * 10,
            "stratum": ["melanoma_core_high_evidence"] * 10,
            "cohort": ["A"] * 5 + ["B"] * 5,
            "sample_id": [f"s{i}" for i in range(10)],
            "model_name": ["EcoNiche-Opt-ModulePriorFixed"] * 10,
            "true_response_label": [0, 1, 1, 0, 1, 0, 1, 0, 1, 0],
            "response_probability": [0.3, 0.6, 0.35, 0.4, 0.5, 0.55, 0.52, 0.58, 0.48, 0.42],
        }
    )

    comparison = compare_to_baseline(ensemble, baseline)

    assert len(comparison) == 1
    assert comparison.iloc[0]["target_model"] == ENSEMBLE_MODEL
    assert comparison.iloc[0]["delta_AUROC"] > 0
