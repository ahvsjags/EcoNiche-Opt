from __future__ import annotations

import pandas as pd

from scripts.analysis.run_signed_rank_module_audit import estimate_gene_directions, signed_rank_module_features


def test_estimate_gene_directions_uses_training_labels_only():
    ranked = {
        "train_a": pd.DataFrame(
            {
                "IFNG": [0.1, 0.2, 1.0, 1.2, 1.4, 1.6],
                "COL1A1": [1.5, 1.3, 0.2, 0.1, 0.0, -0.1],
            },
            index=[f"S{i}" for i in range(6)],
        )
    }
    y = {"train_a": pd.Series([0, 0, 1, 1, 1, 1], index=ranked["train_a"].index)}

    directions = estimate_gene_directions(ranked, y, ["train_a"])

    assert directions["IFNG"] == 1
    assert directions["COL1A1"] == -1


def test_signed_rank_module_features_returns_expected_modules():
    ranked = pd.DataFrame(
        {
            "IFNG": [1.0, 2.0],
            "CXCL9": [2.0, 3.0],
            "COL1A1": [4.0, 5.0],
        },
        index=["S1", "S2"],
    )
    features, coverage = signed_rank_module_features(ranked, directions={"COL1A1": -1}, signed=True)

    assert "ifn_t_cell_inflamed" in features.columns
    assert "stromal_exclusion" in features.columns
    assert coverage["n_genes_available"].max() >= 1
    assert features.loc["S1", "stromal_exclusion"] < 0
