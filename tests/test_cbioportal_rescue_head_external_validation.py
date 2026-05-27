from __future__ import annotations

import pandas as pd

from scripts.analysis.run_cbioportal_rescue_head_external_validation import build_selection


def test_cbioportal_rescue_selection_preserves_no_external_selection_boundary():
    primary = pd.DataFrame(
        [
            {
                "blend_id": "cohort_zscore",
                "blend_type": "single",
                "AUROC": 0.74,
                "AUPRC": 0.70,
                "balanced_accuracy": 0.66,
            },
            {
                "blend_id": "0.05*cohort_zscore+0.95*cohort_robust_zscore",
                "blend_type": "pair_grid",
                "AUROC": 0.73,
                "AUPRC": 0.69,
                "balanced_accuracy": 0.65,
            },
        ]
    )
    external = pd.DataFrame(
        [
            {
                "group_id": "cbio_liu_dfci_only",
                "blend_id": "cohort_zscore",
                "blend_type": "single",
                "AUROC": 0.60,
                "AUPRC": 0.61,
                "balanced_accuracy": 0.58,
                "ECE": 0.20,
            },
            {
                "group_id": "cbio_liu_dfci_only",
                "blend_id": "0.05*cohort_zscore+0.95*cohort_robust_zscore",
                "blend_type": "pair_grid",
                "AUROC": 0.65,
                "AUPRC": 0.64,
                "balanced_accuracy": 0.62,
                "ECE": 0.18,
            },
        ]
    )
    per_cohort = external.assign(cohort="CBIO_LIU_DFCI_2019_PRE")

    selection = build_selection(primary, external, per_cohort)

    boundaries = set(selection["claim_boundary"])
    assert "selected_by_primary_lodo_only_not_by_cbio_external" in boundaries
    assert "fixed_robust_transform_candidate_no_cbio_external_label_fit" in boundaries
    assert "diagnostic_current_cbio_external_stress_screen_not_selection_claim" in boundaries
