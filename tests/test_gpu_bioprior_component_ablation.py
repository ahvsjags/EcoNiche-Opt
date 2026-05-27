from __future__ import annotations

import pandas as pd

from scripts.analysis.run_gpu_bioprior_component_ablation import _paired_metric_rows, _selected_candidate


def test_selected_candidate_requires_primary_lodo_boundary(tmp_path):
    selection = pd.DataFrame(
        [
            {
                "selection_id": "gpu_lipid_pi3k_primary_selected_rescue_combo",
                "candidate": "0.80*base+0.20*rz__PLA2G2D",
                "strict_external_AUROC": 0.713,
                "selection_boundary": "candidate_and_weight_selected_by_primary_lodo_only_with_biological_prior",
            }
        ]
    )
    path = tmp_path / "selection.tsv"
    selection.to_csv(path, sep="\t", index=False)

    assert _selected_candidate(path) == "0.80*base+0.20*rz__PLA2G2D"


def test_paired_metric_rows_report_directional_gain():
    predictions = pd.DataFrame(
        [
            {"cohort": "C1", "sample_id": "S1", "candidate": "target", "true_response_label": 1, "response_probability": 0.9, "threshold": 0.5},
            {"cohort": "C1", "sample_id": "S2", "candidate": "target", "true_response_label": 1, "response_probability": 0.7, "threshold": 0.5},
            {"cohort": "C1", "sample_id": "S3", "candidate": "target", "true_response_label": 0, "response_probability": 0.4, "threshold": 0.5},
            {"cohort": "C1", "sample_id": "S4", "candidate": "target", "true_response_label": 0, "response_probability": 0.2, "threshold": 0.5},
            {"cohort": "C1", "sample_id": "S1", "candidate": "base", "true_response_label": 1, "response_probability": 0.8, "threshold": 0.5},
            {"cohort": "C1", "sample_id": "S2", "candidate": "base", "true_response_label": 1, "response_probability": 0.3, "threshold": 0.5},
            {"cohort": "C1", "sample_id": "S3", "candidate": "base", "true_response_label": 0, "response_probability": 0.6, "threshold": 0.5},
            {"cohort": "C1", "sample_id": "S4", "candidate": "base", "true_response_label": 0, "response_probability": 0.1, "threshold": 0.5},
        ]
    )
    summary = pd.DataFrame(
        [
            {"candidate": "target", "AUROC": 1.0, "AUPRC": 1.0},
            {"candidate": "base", "AUROC": 0.625, "AUPRC": 0.75},
        ]
    )

    rows = _paired_metric_rows(predictions, summary, "unit", "target", "base", 200)
    by_metric = {row["metric"]: row for row in rows}

    assert by_metric["AUROC"]["delta"] > 0
    assert by_metric["AUPRC"]["delta"] > 0
    assert by_metric["balanced_accuracy"]["target_value"] > by_metric["balanced_accuracy"]["ablated_value"]
