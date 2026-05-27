from __future__ import annotations

import pandas as pd

from scripts.analysis.run_aligned_biological_objective_audit import (
    BIO_MODEL,
    NO_BIO_MODEL,
    biological_alignment_terms,
    paired_comparison,
    panel_weight_candidates,
)


def test_panel_weight_candidates_include_bio_prior_and_drop_variants():
    candidates = panel_weight_candidates()
    names = set(candidates["candidate"])

    assert "bio_prior" in names
    assert "drop_ifn_t_cell_inflamed" in names
    assert {"candidate", "candidate_family", "ifn_t_cell_inflamed", "stromal_exclusion"}.issubset(candidates.columns)
    assert len(candidates) >= 10


def test_biological_alignment_terms_reward_prior_direction():
    candidates = panel_weight_candidates().set_index("candidate")
    prior_terms = biological_alignment_terms(candidates.loc["bio_prior"])
    inverted = candidates.loc["bio_prior"].copy()
    for column in ["ifn_t_cell_inflamed", "cytotoxic_cd8", "antigen_presentation", "myeloid_suppression"]:
        inverted[column] = -float(inverted[column])
    inverted_terms = biological_alignment_terms(inverted)

    assert prior_terms["bio_objective_bonus"] > inverted_terms["bio_objective_bonus"]
    assert prior_terms["bio_prior_cosine"] > 0.99


def test_paired_comparison_reports_bio_vs_no_bio_delta():
    sample_ids = [f"s{i}" for i in range(10)]
    labels = [0, 1, 1, 0, 1, 0, 1, 0, 1, 0]
    predictions = pd.DataFrame(
        {
            "endpoint": ["primary_recist"] * 20,
            "stratum": ["melanoma_core_high_evidence"] * 20,
            "cohort": (["A"] * 5 + ["B"] * 5) * 2,
            "sample_id": sample_ids * 2,
            "true_response_label": labels * 2,
            "model_name": [BIO_MODEL] * 10 + [NO_BIO_MODEL] * 10,
            "selection_mode": ["bio_objective"] * 10 + ["no_bio_objective"] * 10,
            "selected_candidate": ["bio_prior"] * 20,
            "response_probability": [0.1, 0.9, 0.8, 0.2, 0.7, 0.3, 0.75, 0.25, 0.85, 0.15, 0.3, 0.6, 0.35, 0.4, 0.5, 0.55, 0.52, 0.58, 0.48, 0.42],
        }
    )

    comparison = paired_comparison(predictions)

    assert len(comparison) == 1
    assert comparison.iloc[0]["target_model"] == BIO_MODEL
    assert comparison.iloc[0]["ablation_model"] == NO_BIO_MODEL
    assert comparison.iloc[0]["delta_AUROC"] > 0
