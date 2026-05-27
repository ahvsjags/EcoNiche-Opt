import numpy as np
import pandas as pd

from scripts.analysis.add_auroc_confidence_intervals import bootstrap_auroc_ci, primary_ci


def test_bootstrap_auroc_ci_has_bounds():
    stats = bootstrap_auroc_ci(
        pd.Series([0, 0, 0, 1, 1, 1]),
        pd.Series([0.1, 0.2, 0.3, 0.7, 0.8, 0.9]),
        n_bootstrap=50,
        seed=7,
    )
    assert stats["AUROC"] == 1.0
    assert 0.0 <= stats["AUROC_ci_low"] <= stats["AUROC_ci_high"] <= 1.0
    assert stats["AUROC_ci_bootstrap_n"] > 0


def test_primary_ci_groups_by_endpoint_stratum_model():
    predictions = pd.DataFrame(
        {
            "endpoint": ["primary_recist"] * 6,
            "stratum": ["melanoma_core_high_evidence"] * 6,
            "model_name": ["EcoNiche-Opt-HeuristicEcology"] * 6,
            "true_response_label": [0, 0, 0, 1, 1, 1],
            "response_probability": [0.1, 0.2, 0.4, 0.6, 0.7, 0.8],
        }
    )
    out = primary_ci(predictions, n_bootstrap=25)
    assert len(out) == 1
    assert np.isfinite(out.loc[0, "AUROC_ci_low"])
    assert out.loc[0, "n_responders_ci"] == 3
