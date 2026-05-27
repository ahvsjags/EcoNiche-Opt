from __future__ import annotations

import pandas as pd

from scripts.analysis.run_tumor_immune_balance_audit import BALANCE_AXES, balance_score


def test_map4k1_tbx3_axis_is_registered():
    assert "MAP4K1_minus_TBX3" in BALANCE_AXES
    positive, negative = BALANCE_AXES["MAP4K1_minus_TBX3"]
    assert positive == ["MAP4K1"]
    assert negative == ["TBX3"]


def test_balance_score_increases_with_positive_gene_and_decreases_with_negative_gene():
    ranked = pd.DataFrame(
        {
            "MAP4K1": [0.9, 0.1],
            "TBX3": [0.1, 0.9],
        },
        index=["immune_high", "tumor_high"],
    )

    score, coverage = balance_score(ranked, ["MAP4K1"], ["TBX3"])

    assert score.loc["immune_high"] > score.loc["tumor_high"]
    assert coverage["n_positive_genes_available"] == 1
    assert coverage["n_negative_genes_available"] == 1
