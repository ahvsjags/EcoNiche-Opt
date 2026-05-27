from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from scripts.reporting import make_gpu_bioprior_rescue_combo_package as package


def test_gpu_bioprior_package_freezes_candidate_and_threshold(tmp_path: Path, monkeypatch):
    selection = pd.DataFrame(
        [
            {
                "selection_id": "gpu_lipid_pi3k_primary_selected_rescue_combo",
                "candidate": "0.80*base+0.20*rz__PLA2G2D",
                "prior": "lipid_pi3k",
                "transform_policy": "robust_only",
                "device": "cuda",
                "gpu_name": "NVIDIA GeForce RTX 4070 Laptop GPU",
                "selection_boundary": "candidate_and_weight_selected_by_primary_lodo_only_with_biological_prior",
                "primary_AUROC": 0.777,
                "primary_AUPRC": 0.665,
                "primary_balanced_accuracy": 0.659,
                "strict_external_AUROC": 0.713,
                "strict_external_AUPRC": 0.726,
                "strict_external_balanced_accuracy": 0.629,
                "strict_external_ECE": 0.128,
                "family_mean_AUROC": 0.567,
                "delta_vs_family_mean": 0.146,
                "two_sided_fdr_q": 0.0,
                "claim_level": "strict_external_family_FDR_supported_numeric_target_met",
            }
        ]
    )
    predictions = pd.DataFrame(
        [
            {
                "candidate": "0.80*base+0.20*rz__PLA2G2D",
                "threshold": 0.42,
                "cohort": "GSE145996",
                "sample_id": "S1",
                "response_probability": 0.8,
                "true_response_label": 1,
            }
        ]
    )
    family_gate = pd.DataFrame(
        [
            {
                "target_model": "0.80*base+0.20*rz__PLA2G2D",
                "family_mean_AUROC": 0.567,
                "best_signature": "APM",
                "best_signature_AUROC": 0.64,
                "two_sided_fdr_q": 0.0,
            }
        ]
    )
    selection_path = tmp_path / "selection.tsv"
    predictions_path = tmp_path / "predictions.tsv"
    family_gate_path = tmp_path / "family_gate.tsv"
    selection.to_csv(selection_path, sep="\t", index=False)
    predictions.to_csv(predictions_path, sep="\t", index=False)
    family_gate.to_csv(family_gate_path, sep="\t", index=False)
    monkeypatch.setattr(package, "_component_sign", lambda processed_dir, gene, method: (1, 0.31))

    out = tmp_path / "package"
    summary = package.build_package(selection_path, predictions_path, family_gate_path, tmp_path, out, "v-test")

    assert summary["locked_candidate"] == "0.80*base+0.20*rz__PLA2G2D"
    assert summary["locked_threshold"] == 0.42
    assert summary["component_direction_sign"] == 1
    assert (out / "gpu_bioprior_rescue_combo_scoring_spec.sha256").exists()
    spec = json.loads((out / "gpu_bioprior_rescue_combo_scoring_spec.json").read_text(encoding="utf-8"))
    assert spec["model_status"] == "frozen_no_leakage_lipid_pi3k_rescue_combo_for_future_validation"
    assert spec["locked_score"]["component_gene"] == "PLA2G2D"
    assert spec["locked_score"]["weight_base"] == 0.8
    assert spec["performance_evidence"]["strict_external_AUROC"] == 0.713
    assert "Strict external labels must not enter candidate-weight selection." in spec["no_leakage_rules"]
    genes = pd.read_csv(out / "gpu_bioprior_rescue_combo_genes.tsv", sep="\t")
    assert set(genes["gene_symbol"]) == {"MAP4K1", "TBX3", "AXL", "PLA2G2D"}
