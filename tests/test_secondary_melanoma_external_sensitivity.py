from __future__ import annotations

import pandas as pd

from scripts.analysis.run_secondary_melanoma_external_sensitivity import _candidate_score


def test_candidate_score_higher_for_response_like_axis():
    X = pd.DataFrame(
        {
            "MAP4K1": [10.0, 1.0, 5.0],
            "TBX3": [1.0, 10.0, 5.0],
            "AXL": [2.0, 9.0, 5.0],
        },
        index=["response_like", "resistance_like", "middle"],
    )
    spec = {
        "transform": "cohort_gene_percentile",
        "positive_genes": ["MAP4K1"],
        "negative_genes": ["TBX3", "AXL"],
        "negative_weight": 1.25,
    }

    score = _candidate_score(X, spec)

    assert score.loc["response_like"] > score.loc["middle"] > score.loc["resistance_like"]


def test_candidate_score_supports_robust_zscore_transform():
    X = pd.DataFrame(
        {
            "MAP4K1": [5.0, 1.0, 3.0, 100.0],
            "TBX3": [1.0, 5.0, 3.0, 100.0],
            "AXL": [1.0, 4.0, 3.0, 100.0],
        },
        index=["response_like", "resistance_like", "middle", "outlier"],
    )
    spec = {
        "transform": "cohort_robust_zscore",
        "positive_genes": ["MAP4K1"],
        "negative_genes": ["TBX3", "AXL"],
        "negative_weight": 1.0,
    }

    score = _candidate_score(X, spec)

    assert score.loc["response_like"] > score.loc["resistance_like"]
