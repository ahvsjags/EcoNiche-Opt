from __future__ import annotations

import json

import pandas as pd

from scripts.analysis.run_cbioportal_gpu_bioprior_external_validation import score_gpu_combo


def test_gpu_combo_scoring_requires_all_locked_genes():
    spec = {
        "rescue_combo_id": "test",
        "locked_score": {
            "weight_base": 0.8,
            "weight_component": 0.2,
            "component_training_direction_sign": 1,
            "component_gene": "PLA2G2D",
            "locked_threshold": 0.4,
        },
    }
    X = pd.DataFrame(
        {
            "MAP4K1": [3, 4, 5],
            "TBX3": [5, 4, 3],
            "AXL": [4, 4, 2],
        },
        index=["S1", "S2", "S3"],
    )

    try:
        score_gpu_combo(X, spec)
    except ValueError as exc:
        assert "PLA2G2D" in str(exc)
    else:
        raise AssertionError("Expected missing PLA2G2D to raise")


def test_gpu_combo_scoring_returns_bounded_response_probability():
    spec = {
        "rescue_combo_id": "test",
        "locked_score": {
            "weight_base": 0.8,
            "weight_component": 0.2,
            "component_training_direction_sign": 1,
            "component_gene": "PLA2G2D",
            "locked_threshold": 0.4,
        },
    }
    X = pd.DataFrame(
        {
            "MAP4K1": [1.0, 2.0, 5.0, 7.0],
            "TBX3": [5.0, 4.0, 2.0, 1.0],
            "AXL": [4.0, 4.0, 2.0, 1.0],
            "PLA2G2D": [1.0, 2.0, 3.0, 8.0],
        },
        index=["S1", "S2", "S3", "S4"],
    )

    score = score_gpu_combo(X, json.loads(json.dumps(spec)))

    assert score.index.tolist() == ["S1", "S2", "S3", "S4"]
    assert float(score.min()) >= 0.0
    assert float(score.max()) <= 1.0
    assert score.loc["S4"] > score.loc["S1"]
