from __future__ import annotations

import pandas as pd

from scripts.analysis.run_threshold_recalibration_audit import evaluate_primary_thresholds, select_midrange_policy


def _toy_predictions() -> pd.DataFrame:
    rows = []
    values = {
        "A": [(0, 0.10), (0, 0.20), (1, 0.45), (1, 0.80)],
        "B": [(0, 0.15), (0, 0.35), (1, 0.50), (1, 0.85)],
        "C": [(0, 0.25), (0, 0.55), (1, 0.70), (1, 0.90)],
    }
    for fold, pairs in values.items():
        for idx, (label, prob) in enumerate(pairs):
            rows.append(
                {
                    "endpoint": "primary_recist",
                    "stratum": "melanoma_core_high_evidence",
                    "model_name": "EcoNiche-Opt-HeuristicEcology",
                    "fold": fold,
                    "cohort": fold,
                    "sample_id": f"{fold}_{idx}",
                    "true_response_label": label,
                    "response_probability": prob,
                    "pred_response_label": int(prob >= 0.5),
                    "threshold": 0.5,
                }
            )
    return pd.DataFrame(rows)


def test_select_midrange_policy_uses_training_rows_only():
    frame = _toy_predictions()
    best, selection = select_midrange_policy(frame[frame["fold"].isin(["A", "B"])])

    assert best in {"fixed_0.40", "fixed_0.50", "fixed_0.60"}
    assert set(selection["candidate_policy"]) == {"fixed_0.40", "fixed_0.50", "fixed_0.60"}
    assert selection["inner_n_validations"].min() == 2


def test_evaluate_primary_thresholds_adds_nested_policy_without_external_labels():
    frame = _toy_predictions()
    summary, predictions, selection = evaluate_primary_thresholds(frame)

    policies = set(summary["threshold_policy"])
    assert "original_fold_threshold" in policies
    assert "nested_midrange_fixed_grid" in policies
    nested = summary[summary["threshold_policy"].eq("nested_midrange_fixed_grid")].iloc[0]
    assert nested["selection_boundary"] == "outer_lodo_training_cohorts_only"
    assert nested["n_samples"] == len(frame)
    assert not predictions.empty
    assert not selection.empty
    assert set(predictions["threshold_policy"]) >= {"nested_midrange_fixed_grid"}
