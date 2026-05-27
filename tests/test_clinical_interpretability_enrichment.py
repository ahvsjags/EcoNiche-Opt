from __future__ import annotations

import pandas as pd

from scripts.analysis.run_clinical_interpretability_enrichment import (
    build_outputs,
    classify_therapy_context,
    high_score_enrichment,
)


def test_high_score_enrichment_uses_locked_threshold_counts():
    frame = pd.DataFrame(
        {
            "source_context": ["primary_lodo"] * 6,
            "endpoint": ["primary_recist"] * 6,
            "model_name": ["EcoNiche-Opt-HeuristicEcology"] * 6,
            "true_response_label": [1, 1, 1, 0, 0, 0],
            "response_probability": [0.9, 0.8, 0.7, 0.6, 0.2, 0.1],
            "locked_or_fold_threshold": [0.65] * 6,
        }
    )

    row = high_score_enrichment(frame, "unit", "standard_context")

    assert row["n_high_score"] == 3
    assert row["responders_high_score"] == 3
    assert row["response_rate_high_score"] == 1.0
    assert row["response_rate_low_score"] == 0.0
    assert row["odds_ratio_haldane"] > 1.0


def test_build_outputs_includes_strict_external_context_and_subgroup_axes():
    primary = pd.DataFrame(
        {
            "model_name": ["EcoNiche-Opt-HeuristicEcology"] * 8,
            "endpoint": ["primary_recist"] * 8,
            "stratum": ["melanoma_core_high_evidence"] * 8,
            "cohort": ["GSE91061"] * 4 + ["GSE78220"] * 4,
            "sample_id": [f"p{i}" for i in range(8)],
            "true_response_label": [1, 0, 1, 0, 1, 0, 1, 0],
            "response_probability": [0.8, 0.2, 0.7, 0.3, 0.75, 0.25, 0.65, 0.35],
            "threshold": [0.5] * 8,
            "pred_response_label": [1, 0, 1, 0, 1, 0, 1, 0],
            "timepoint": ["pretreatment"] * 8,
            "treatment": ["Nivolumab"] * 8,
        }
    )
    external = pd.DataFrame(
        {
            "model_name": ["EcoNiche-Opt-HeuristicEcology-LockedPanel"] * 10,
            "endpoint": ["strict_recist"] * 10,
            "analysis_type": ["locked_external_melanoma_pd1_recist"] * 5 + ["locked_external_melanoma_pd1_like"] * 5,
            "cohort": ["GSE145996"] * 5 + ["PHS000452_LIU_LIKE_PRE"] * 5,
            "sample_id": [f"e{i}" for i in range(10)],
            "true_response_label": [1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
            "response_probability": [0.8, 0.2, 0.75, 0.25, 0.7, 0.3, 0.65, 0.35, 0.6, 0.4],
            "locked_threshold": [0.5] * 10,
            "predicted_response_label": [1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
            "threshold_source": ["discovery_only"] * 10,
            "platform": [""] * 10,
            "treatment": ["anti-PD-1"] * 10,
        }
    )

    outputs = build_outputs(primary, external)

    context_ids = set(outputs["high_score_enrichment"]["context_id"])
    assert "locked_external|strict_recist|strict_melanoma_pd1_like_pooled" in context_ids
    assert "primary_lodo|primary_recist|melanoma_core_high_evidence" in context_ids
    subgroup_axes = set(outputs["subgroup_metrics"]["subgroup_axis"])
    assert {"cohort", "therapy_context", "platform_context", "sampling_context", "analysis_type"}.issubset(subgroup_axes)
    assert not outputs["threshold_operating_points"].empty
    assert not outputs["calibration_bins"].empty


def test_classify_therapy_context_separates_pd1_and_combo():
    assert classify_therapy_context("Nivolumab") == "anti_pd1_monotherapy"
    assert classify_therapy_context("Ipilimumab + Nivolumab") == "combo_therapy"
