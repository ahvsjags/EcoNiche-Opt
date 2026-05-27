from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from scripts.reporting.make_melanoma_rescue_head_package import build_package


def test_rescue_head_package_freezes_primary_selected_candidate(tmp_path: Path):
    selection = pd.DataFrame(
        [
            {
                "selection_id": "primary_selected_candidate",
                "claim_boundary": "selected_by_primary_lodo_only_not_by_external",
                "method": "cohort_gene_percentile",
                "axis": "MAP4K1_minus_TBX3_AXL",
                "primary_AUROC": 0.785,
                "primary_AUPRC": 0.693,
                "primary_balanced_accuracy": 0.725,
                "strict_external_AUROC": 0.661,
                "strict_external_AUPRC": 0.610,
                "strict_external_balanced_accuracy": 0.634,
            },
            {
                "selection_id": "current_external_stress_best",
                "claim_boundary": "current_external_stress_screen_not_a_locked_selection_claim",
                "method": "cohort_zscore",
                "axis": "MAP4K1_minus_TBX3_AXL",
                "primary_AUROC": 0.779,
                "primary_AUPRC": 0.670,
                "primary_balanced_accuracy": 0.722,
                "strict_external_AUROC": 0.679,
                "strict_external_AUPRC": 0.684,
                "strict_external_balanced_accuracy": 0.602,
            },
        ]
    )
    selection_path = tmp_path / "selection.tsv"
    selection.to_csv(selection_path, sep="\t", index=False)
    out = tmp_path / "rescue"

    summary = build_package(selection_path, out, "v-test")

    assert summary["locked_method"] == "cohort_gene_percentile"
    assert summary["locked_axis"] == "MAP4K1_minus_TBX3_AXL"
    assert (out / "melanoma_rescue_head_genes.tsv").exists()
    assert (out / "melanoma_rescue_head_scoring_spec.sha256").exists()
    spec = json.loads((out / "melanoma_rescue_head_scoring_spec.json").read_text(encoding="utf-8"))
    assert spec["model_status"] == "frozen_extension_for_future_locked_external_validation"
    assert spec["primary_development_evidence"]["primary_AUROC"] == 0.785
    assert spec["stress_screen_not_for_locked_selection"]["method"] == "cohort_zscore"
    assert "External labels must not enter transform selection." in spec["no_leakage_rules"]
    genes = pd.read_csv(out / "melanoma_rescue_head_genes.tsv", sep="\t")
    assert set(genes["gene_symbol"]) == {"MAP4K1", "TBX3", "AXL"}
