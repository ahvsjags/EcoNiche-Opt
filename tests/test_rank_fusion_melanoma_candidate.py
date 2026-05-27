from __future__ import annotations

import pandas as pd

from scripts.analysis.run_rank_fusion_melanoma_candidate import build_rank_fusion_scores, candidate_specs, rank_percentile


def test_rank_percentile_returns_half_for_degenerate_series():
    ranked = rank_percentile(pd.Series([5.0, 5.0, 5.0], index=["a", "b", "c"]))

    assert ranked.between(0, 1).all()
    assert ranked.loc["a"] == ranked.loc["b"]


def test_rank_fusion_scores_apply_cytotoxic_penalty():
    fixed_scores = {
        "cohort_a": {
            "EcoNiche-Opt-ModulePriorFixed": pd.Series([1.0, 2.0, 3.0], index=["s1", "s2", "s3"]),
            "CYT": pd.Series([3.0, 2.0, 1.0], index=["s1", "s2", "s3"]),
        }
    }
    spec = {"candidate": "toy", "weights": {"EcoNiche-Opt-ModulePriorFixed": 1.0, "CYT": -0.5}}

    fused = build_rank_fusion_scores(fixed_scores, spec)["cohort_a"]

    assert fused.loc["s3"] > fused.loc["s1"]
    assert any("cytotoxic_penalty" in str(item["candidate"]) for item in candidate_specs())
