from __future__ import annotations

import pandas as pd

from scripts.analysis.run_training_only_melanoma_candidate_search import (
    _summarize_external_predictions,
    _summarize_lodo_predictions,
)


def test_summarize_lodo_predictions_reports_prevalence_margin():
    predictions = pd.DataFrame(
        {
            "endpoint": ["primary_recist"] * 6,
            "stratum": ["melanoma_core_high_evidence"] * 6,
            "cohort": ["A", "A", "B", "B", "C", "C"],
            "true_response_label": [0, 1, 0, 1, 0, 1],
            "response_probability": [0.1, 0.9, 0.2, 0.8, 0.3, 0.7],
            "selected_candidate": ["module_prior_composite"] * 6,
        }
    )
    summary = _summarize_lodo_predictions(predictions)
    assert len(summary) == 1
    assert summary.loc[0, "AUROC"] == 1.0
    assert summary.loc[0, "response_prevalence"] == 0.5
    assert summary.loc[0, "AUPRC_minus_prevalence"] > 0


def test_summarize_external_predictions_preserves_selection_boundary():
    predictions = pd.DataFrame(
        {
            "endpoint": ["strict_recist"] * 6,
            "cohort": ["GSE145996", "GSE145996", "PHS000452_LIU_LIKE_PRE", "PHS000452_LIU_LIKE_PRE", "PHS000452_LIU_LIKE_PRE", "PHS000452_LIU_LIKE_PRE"],
            "true_response_label": [0, 1, 0, 1, 0, 1],
            "response_probability": [0.1, 0.9, 0.2, 0.8, 0.4, 0.6],
            "selected_candidate": ["module_prior_composite"] * 6,
        }
    )
    summary = _summarize_external_predictions(predictions)
    assert len(summary) == 1
    assert summary.loc[0, "cohort_set"] == "GSE145996+PHS000452_LIU_LIKE_PRE"
    assert summary.loc[0, "selection_boundary"] == "discovery_only_inner_lodo_no_external_selection"
    assert summary.loc[0, "AUROC"] == 1.0
