from __future__ import annotations

import numpy as np
import pandas as pd

from econiche_opt.model.ecology_optimizer import HeuristicEcologyConfig
from econiche_opt.model.endpoint_modules import (
    WORD_FULL_GRAPH_MODEL,
    build_word_ecology_features,
    build_module_features,
    default_strata,
    endpoint_response_label,
    evaluate_module_model,
    module_prior_score,
)


def test_endpoint_response_label_sensitivity_rules():
    assert endpoint_response_label("PR", "strict_recist") == 1
    assert np.isnan(endpoint_response_label("MR", "strict_recist"))
    assert endpoint_response_label("MR", "primary_recist") == 1
    assert endpoint_response_label("PD", "strict_recist") == 0
    assert np.isnan(endpoint_response_label("SD", "strict_recist"))
    assert endpoint_response_label("SD", "primary_recist") == 0
    assert endpoint_response_label("SD", "clinical_benefit") == 1
    assert endpoint_response_label("DCB", "primary_recist") == 1
    assert endpoint_response_label("NDB", "clinical_benefit") == 0


def test_module_features_are_finite_with_missing_genes():
    X = pd.DataFrame(
        {
            "IFNG": [1.0, 2.0, 3.0],
            "CXCL9": [3.0, 2.0, 1.0],
            "CD8A": [0.5, 1.0, 1.5],
            "COL1A1": [8.0, 8.0, 8.0],
        },
        index=["s1", "s2", "s3"],
    )
    features, coverage = build_module_features(X)
    assert set(["ifn_t_cell_inflamed", "cytotoxic_cd8", "stromal_exclusion"]).issubset(features.columns)
    assert np.isfinite(features.to_numpy()).all()
    assert coverage["n_genes_available"].ge(0).all()
    score = module_prior_score(features)
    assert score.index.tolist() == ["s1", "s2", "s3"]
    assert np.isfinite(score.to_numpy()).all()


def test_melanoma_recist_supported_stratum_excludes_binary_response_stress_cohorts():
    strata = default_strata(
        [
            "GSE91061",
            "GSE78220",
            "GSE168204",
            "GSE115821",
            "GSE145996",
            "PRJEB23709_PD1_PRE",
        ]
    )
    assert strata["melanoma_recist_supported_primary"]["cohorts"] == [
        "GSE91061",
        "GSE78220",
        "GSE145996",
        "PRJEB23709_PD1_PRE",
    ]
    assert strata["melanoma_binary_response_stress"]["cohorts"] == ["GSE168204", "GSE115821"]


def test_word_ecology_features_include_signed_states_and_interactions():
    X = pd.DataFrame(
        {
            "IFNG": [1.0, 2.0, 3.0, 4.0],
            "CXCL9": [1.0, 2.0, 3.0, 4.0],
            "HLA-A": [1.0, 1.5, 2.0, 2.5],
            "B2M": [1.0, 1.5, 2.0, 2.5],
            "COL1A1": [4.0, 3.0, 2.0, 1.0],
            "FN1": [4.0, 3.0, 2.0, 1.0],
            "S100A8": [4.0, 3.0, 2.0, 1.0],
            "S100A9": [4.0, 3.0, 2.0, 1.0],
        },
        index=["s1", "s2", "s3", "s4"],
    )
    features, coverage = build_word_ecology_features(X, gene_directions={"COL1A1": -1, "FN1": -1})
    assert {"tnk_effector", "caf_ecm_exclusion", "interaction__antigen_presentation_mhc__tnk_effector"}.issubset(features.columns)
    assert np.isfinite(features.to_numpy()).all()
    assert {"word_state", "word_interaction"}.issubset(set(coverage["feature_type"]))


def test_evaluate_module_model_emits_word_full_graph_ablation_rows():
    X1 = pd.DataFrame(
        {
            "IFNG": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            "CXCL9": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            "GZMB": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            "PRF1": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            "COL1A1": [6.0, 5.0, 4.0, 3.0, 2.0, 1.0],
            "FN1": [6.0, 5.0, 4.0, 3.0, 2.0, 1.0],
            "S100A8": [6.0, 5.0, 4.0, 3.0, 2.0, 1.0],
            "S100A9": [6.0, 5.0, 4.0, 3.0, 2.0, 1.0],
        },
        index=[f"a{i}" for i in range(6)],
    )
    X2 = X1.copy()
    X2.index = [f"b{i}" for i in range(6)]
    X3 = X1.copy()
    X3.index = [f"c{i}" for i in range(6)]
    X_by_cohort = {"c1": X1, "c2": X2, "c3": X3}
    y_by_cohort = {
        "c1": pd.Series([0, 0, 0, 1, 1, 1], index=X1.index),
        "c2": pd.Series([0, 0, 0, 1, 1, 1], index=X2.index),
        "c3": pd.Series([0, 0, 0, 1, 1, 1], index=X3.index),
    }
    metadata = {
        cohort: pd.DataFrame(
            {
                "sample_id": frame.index,
                "patient_id": frame.index,
                "cohort": cohort,
                "response_raw": ["PD", "PD", "PD", "PR", "PR", "PR"],
                "treatment": "anti-PD1",
            },
            index=frame.index,
        )
        for cohort, frame in X_by_cohort.items()
    }
    module_features = {cohort: build_module_features(frame)[0] for cohort, frame in X_by_cohort.items()}
    result = evaluate_module_model(
        module_features,
        y_by_cohort,
        metadata,
        endpoint="primary_recist",
        stratum="toy",
        train_pool=["c1", "c2", "c3"],
        holdouts=["c1"],
        raw_X_by_cohort=X_by_cohort,
        optimizer_config=HeuristicEcologyConfig(population_size=4, generations=2, use_gpu=False),
    )
    assert WORD_FULL_GRAPH_MODEL in set(result.metrics["model_name"])
    assert WORD_FULL_GRAPH_MODEL in set(result.predictions["model_name"])
    assert "optimized_module_gene" in set(result.feature_weights["weight_type"])
    assert "optimizer_history" in set(result.feature_weights["weight_type"])
