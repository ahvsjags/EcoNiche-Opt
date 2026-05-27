from __future__ import annotations

import numpy as np
import pandas as pd

from scripts.analysis.run_ecological_polarity_candidate_audit import score_axis


def test_score_axis_orients_response_minus_resistance_and_scales():
    transformed = pd.DataFrame(
        {
            "MAP4K1": [0.9, 0.1, 0.7],
            "TBX3": [0.1, 0.8, 0.4],
            "AXL": [0.2, 0.9, 0.5],
        },
        index=["s1", "s2", "s3"],
    )

    score, coverage = score_axis(transformed, ["MAP4K1"], ["TBX3", "AXL"], 1.0)

    assert list(score.index) == ["s1", "s2", "s3"]
    assert np.isclose(score.max(), 1.0)
    assert np.isclose(score.min(), 0.0)
    assert score.loc["s1"] > score.loc["s3"] > score.loc["s2"]
    assert coverage["n_positive_genes_available"] == 1
    assert coverage["n_negative_genes_available"] == 2
