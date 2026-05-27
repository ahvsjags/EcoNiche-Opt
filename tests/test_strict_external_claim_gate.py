import pandas as pd

from scripts.validation.audit_strict_external_claim_gate import build_strict_external_claim_gate


def test_strict_external_gate_blocks_weak_primary_claim():
    metrics = pd.DataFrame(
        [
            {
                "endpoint": "strict_recist",
                "cohort": "GSE145996",
                "model_name": "EcoNiche-Opt-HeuristicEcology-LockedPanel",
                "n_samples": 13,
                "AUROC": 0.575,
            },
            {
                "endpoint": "strict_recist",
                "cohort": "PHS000452_LIU_LIKE_PRE",
                "model_name": "EcoNiche-Opt-HeuristicEcology-LockedPanel",
                "n_samples": 103,
                "AUROC": 0.573,
            },
        ]
    )
    family = pd.DataFrame(
        [
            {
                "endpoint": "strict_recist",
                "validation_family": "strict_pd1_like_external",
                "n_samples": 116,
                "target_AUROC": 0.571,
                "mean_signature_AUROC": 0.542,
                "best_signature_AUROC": 0.590,
                "mean_delta_vs_signature_family": 0.029,
                "two_sided_fdr_q": 0.309,
            }
        ]
    )
    rescue = pd.DataFrame(
        [
            {
                "endpoint": "strict_recist",
                "cohort": "GSE145996+PHS000452_LIU_LIKE_PRE",
                "model_name": "EcoNiche-Opt-PD1LikeTransferHead",
                "n_samples": 116,
                "AUROC": 0.595,
            }
        ]
    )
    report, headline = build_strict_external_claim_gate(metrics, family, rescue)
    strict = report[report["gate_id"] == "strict_family_strict_recist"].iloc[0]
    assert strict["claim_status"] == "modest_point_estimate_only"
    assert "Do not claim AUROC >=0.70" in strict["blocked_claim"]
    assert "not a high-strength clinical validation claim" in headline


def test_strict_external_gate_allows_strong_supported_external_claim():
    metrics = pd.DataFrame()
    family = pd.DataFrame(
        [
            {
                "endpoint": "strict_recist",
                "validation_family": "strict_pd1_like_external",
                "n_samples": 120,
                "target_AUROC": 0.72,
                "mean_signature_AUROC": 0.61,
                "best_signature_AUROC": 0.66,
                "mean_delta_vs_signature_family": 0.11,
                "two_sided_fdr_q": 0.01,
            }
        ]
    )
    rescue = pd.DataFrame()
    report, _ = build_strict_external_claim_gate(metrics, family, rescue)
    strict = report[report["gate_id"] == "strict_family_strict_recist"].iloc[0]
    assert strict["claim_status"] == "primary_external_claim_supported"
    assert strict["blocked_claim"] == ""
