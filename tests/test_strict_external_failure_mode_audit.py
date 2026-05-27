from __future__ import annotations

import pandas as pd

from scripts.analysis.run_strict_external_failure_mode_audit import blend_score, build_method_scores


def test_blend_score_keeps_response_like_sample_high():
    X = pd.DataFrame(
        {
            "MAP4K1": [8.0, 1.0, 4.0],
            "TBX3": [1.0, 8.0, 4.0],
            "AXL": [1.0, 7.0, 4.0],
        },
        index=["response_like", "resistance_like", "middle"],
    )
    scores = build_method_scores({"cohort_a": X})
    spec = {
        "blend_id": "0.95*cohort_robust_zscore+0.05*cohort_zscore",
        "blend_type": "pair_grid",
        "cohort_robust_zscore": 0.95,
        "cohort_zscore": 0.05,
    }

    blended = blend_score(scores, "cohort_a", spec)

    assert blended.loc["response_like"] > blended.loc["middle"] > blended.loc["resistance_like"]
