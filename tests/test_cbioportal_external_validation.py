from __future__ import annotations

import pandas as pd

from scripts.analysis.run_cbioportal_melanoma_external_validation import TARGET_MODEL, family_comparison, summarize_metrics


def _prediction_rows(model_name: str, probabilities: list[float]) -> list[dict[str, object]]:
    labels = [0, 0, 0, 0, 1, 1, 1, 1]
    return [
        {
            "endpoint": "strict_recist",
            "group_id": "cbio_liu_dfci_only",
            "cohort": "CBIO_LIU_DFCI_2019_PRE",
            "claim_status": "independent_cbioportal_liu_external",
            "model_name": model_name,
            "sample_id": f"S{i}",
            "true_response_label": labels[i],
            "response_probability": probabilities[i],
            "threshold": 0.5,
        }
        for i in range(8)
    ]


def test_cbioportal_external_summary_and_family_gate():
    rows = []
    rows.extend(_prediction_rows(TARGET_MODEL, [0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9]))
    for baseline in ["IFNG", "CXCL9", "TIG", "TIDE_dysfunction", "TIDE_exclusion", "CYT", "APM", "IPRES"]:
        rows.extend(_prediction_rows(baseline, [0.4, 0.3, 0.2, 0.1, 0.5, 0.6, 0.7, 0.8]))
    predictions = pd.DataFrame(rows)

    metrics = summarize_metrics(predictions)
    comparison = family_comparison(predictions)

    target_metric = metrics[metrics["model_name"].eq(TARGET_MODEL)].iloc[0]
    assert target_metric["AUROC"] == 1.0
    assert not comparison.empty
    assert comparison.iloc[0]["n_signatures"] == 8
    assert comparison.iloc[0]["target_AUROC"] == 1.0
