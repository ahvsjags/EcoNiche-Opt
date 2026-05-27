from __future__ import annotations

import pandas as pd

from scripts.analysis.run_map4k1_tbx3_transform_audit import AXES, transform_axis_score


def test_map4k1_tbx3_axl_axis_is_literature_prior_extension():
    positive, negative, rationale = AXES["MAP4K1_minus_TBX3_AXL"]

    assert positive == ["MAP4K1"]
    assert negative == ["TBX3", "AXL"]
    assert "ipres" in rationale


def test_transform_axis_score_supports_zscore_robust_zscore_and_gene_percentile():
    frame = pd.DataFrame(
        {
            "MAP4K1": [5.0, 1.0, 3.0],
            "TBX3": [1.0, 5.0, 3.0],
            "AXL": [1.0, 4.0, 3.0],
        },
        index=["responder_like", "resistant_like", "middle"],
    )

    zscore, zcoverage = transform_axis_score(frame, ["MAP4K1"], ["TBX3", "AXL"], "cohort_zscore")
    robust, rcoverage = transform_axis_score(frame, ["MAP4K1"], ["TBX3", "AXL"], "cohort_robust_zscore")
    percentile, pcoverage = transform_axis_score(frame, ["MAP4K1"], ["TBX3", "AXL"], "cohort_gene_percentile")

    assert zscore.loc["responder_like"] > zscore.loc["resistant_like"]
    assert robust.loc["responder_like"] > robust.loc["resistant_like"]
    assert percentile.loc["responder_like"] > percentile.loc["resistant_like"]
    assert zcoverage["n_negative_genes_available"] == 2
    assert rcoverage["n_negative_genes_available"] == 2
    assert pcoverage["n_positive_genes_available"] == 1
