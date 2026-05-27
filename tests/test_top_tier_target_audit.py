from __future__ import annotations

import pandas as pd

from scripts.validation.audit_top_tier_targets import build_audit


def test_top_tier_audit_flags_primary_and_external_gaps():
    primary = pd.DataFrame(
        [
            {
                "endpoint": "primary_recist",
                "stratum": "melanoma_core_high_evidence",
                "model_name": "EcoNiche-Opt-HeuristicEcology",
                "n_samples": 100,
                "n_responders": 40,
                "pooled_AUROC": 0.705,
                "pooled_AUPRC": 0.56,
                "pooled_balanced_accuracy": 0.63,
                "pooled_AUROC_ci_low": 0.60,
                "pooled_AUROC_ci_high": 0.80,
            }
        ]
    )
    family = pd.DataFrame(
        [
            {
                "endpoint": "primary_recist",
                "stratum": "melanoma_core_high_evidence",
                "target_model": "EcoNiche-Opt-HeuristicEcology",
                "baseline_family": "eight_strong_signatures",
                "n_signatures": 8,
                "mean_delta_vs_signature_family": 0.07,
                "two_sided_fdr_q": 0.002,
            }
        ]
    )
    strict_gate = pd.DataFrame(
        [
            {
                "gate_id": "strict_family_strict_recist",
                "target_AUROC": 0.571,
                "family_mean_AUROC": 0.542,
                "two_sided_fdr_q": 0.309,
                "claim_status": "modest_point_estimate_only",
            }
        ]
    )
    ablation = pd.DataFrame(
        [
            {"ablation_model": "EcoNiche-Opt-NoResponseModules", "claim_level": "FDR_supported_component_gain"},
            {"ablation_model": "EcoNiche-Opt-NoResistanceModules", "claim_level": "FDR_supported_component_gain"},
            {"ablation_model": "EcoNiche-Opt-UnsignedStateDirection", "claim_level": "FDR_supported_component_gain"},
            {"ablation_model": "EcoNiche-Opt-AlignedPanelCalibrated", "claim_level": "calibration_improves_ECE_with_discrimination_tradeoff"},
        ]
    )
    decision_curve = pd.DataFrame([{"model_name": "EcoNiche-Opt-HeuristicEcology", "threshold": 0.5}])
    training_search = pd.DataFrame(
        [{"selected_candidate": "module_prior_composite", "AUROC": 0.572, "AUPRC": 0.60, "ECE": 0.26}]
    )
    interaction_edge = pd.DataFrame(
        [
            {
                "stratum": "melanoma_core_high_evidence",
                "delta_AUROC": -0.05,
                "delta_ECE": -0.10,
                "claim_level": "interaction_edge_calibration_and_stability_tradeoff",
            }
        ]
    )
    data_candidates = pd.DataFrame(
        [{"accession": "EGAS00001001552", "normalized_access_status": "controlled"}]
    )
    biological_objective = pd.DataFrame(
        [
            {
                "endpoint": "primary_recist",
                "stratum": "melanoma_core_high_evidence",
                "target_model": "EcoNiche-Opt-BioObjectivePanelSearch",
                "ablation_model": "EcoNiche-Opt-NoBioObjectivePanelSearch",
                "delta_AUROC": 0.012,
                "delta_balanced_accuracy": 0.007,
                "delta_ECE": 0.002,
                "fdr_q": 0.417,
            },
            {
                "endpoint": "strict_recist",
                "stratum": "strict_melanoma_pd1_like_external",
                "target_model": "EcoNiche-Opt-BioObjectivePanelSearch",
                "ablation_model": "EcoNiche-Opt-NoBioObjectivePanelSearch",
                "delta_AUROC": -0.007,
                "delta_balanced_accuracy": 0.024,
                "delta_ECE": -0.050,
                "fdr_q": 0.734,
            },
        ]
    )
    public_external_leads = pd.DataFrame(
        [
            {"lead_id": "LIU_MGSP_PHS000452", "eligibility_status": "eligible_processed_duplicate_sensitive"},
            {"lead_id": "GIDE_PRJEB23709", "eligibility_status": "eligible_processed"},
            {"lead_id": "RIAZ_GSE91061", "eligibility_status": "eligible_processed"},
            {"lead_id": "HUGO_GSE78220", "eligibility_status": "eligible_processed"},
            {"lead_id": "LEE_RIZOS_EGAS00001001552", "eligibility_status": "controlled_access_required"},
            {"lead_id": "ABRIL_RODRIGUEZ_PHS001919", "eligibility_status": "controlled_access_required"},
            {"lead_id": "MGH_HACOHEN_PHS002683", "eligibility_status": "controlled_access_required"},
            {"lead_id": "GSE123728", "eligibility_status": "panel_transfer_metadata_pending"},
            {"lead_id": "GSE165745", "eligibility_status": "panel_transfer_public"},
            {"lead_id": "GSE122220", "eligibility_status": "low_n_array_public_processed"},
            {"lead_id": "GSE93157", "eligibility_status": "panel_transfer_processed"},
            {"lead_id": "IMVIGOR210", "eligibility_status": "ineligible_for_melanoma_primary"},
            {"lead_id": "ICBATLAS_TIGER_AGGREGATES", "eligibility_status": "aggregate_duplicate_screen"},
        ]
    )
    gse165745_panel_qc = pd.DataFrame(
        [
            {
                "n_samples": 24,
                "n_responders": 12,
                "locked_panel_overlap": 8,
                "locked_panel_genes": 62,
                "locked_panel_overlap_fraction": 8 / 62,
            }
        ]
    )
    gse165745_panel_metrics = pd.DataFrame(
        [
            {
                "endpoint": "primary_recist",
                "status": "completed",
                "AUROC": 0.542,
                "balanced_accuracy": 0.458,
                "ECE": 0.163,
            }
        ]
    )
    rank_fusion_primary = pd.DataFrame(
        [
            {
                "endpoint": "primary_recist",
                "stratum": "melanoma_core_high_evidence",
                "AUROC": 0.681,
                "balanced_accuracy": 0.638,
                "selected_candidates": "rank_module_cytotoxic_penalty_0.625",
            }
        ]
    )
    rank_fusion_external = pd.DataFrame(
        [
            {
                "endpoint": "strict_recist",
                "stratum": "strict_melanoma_pd1_like_external",
                "AUROC": 0.541,
                "selected_candidates": "rank_module_dysfunction_0.25_cytpen_0.75",
            }
        ]
    )
    tumor_immune_primary = pd.DataFrame(
        [
            {
                "endpoint": "primary_recist",
                "stratum": "melanoma_core_high_evidence",
                "model_name": "EcoNiche-Opt-TumorImmuneBalancePair",
                "AUROC": 0.721,
                "balanced_accuracy": 0.68,
            }
        ]
    )
    tumor_immune_external = pd.DataFrame(
        [
            {
                "endpoint": "strict_recist",
                "stratum": "strict_melanoma_pd1_like_external",
                "model_name": "EcoNiche-Opt-TumorImmuneBalancePair",
                "AUROC": 0.653,
                "balanced_accuracy": 0.59,
            }
        ]
    )
    map4k1_transform_selection = pd.DataFrame(
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

    audit = build_audit(
        primary,
        family,
        strict_gate,
        ablation,
        decision_curve,
        training_search,
        interaction_edge,
        data_candidates,
        biological_objective,
        public_external_leads=public_external_leads,
        gse165745_panel_qc=gse165745_panel_qc,
        gse165745_panel_metrics=gse165745_panel_metrics,
        rank_fusion_primary=rank_fusion_primary,
        rank_fusion_external=rank_fusion_external,
        tumor_immune_primary=tumor_immune_primary,
        tumor_immune_external=tumor_immune_external,
        map4k1_transform_selection=map4k1_transform_selection,
    )
    by_requirement = audit.set_index("requirement")["status"].to_dict()

    assert by_requirement["Primary melanoma LODO AUROC >=0.72"] == "PASS"
    assert by_requirement["Primary melanoma balanced accuracy >=0.65"] == "PARTIAL"
    assert by_requirement["High-value restricted melanoma external candidates are registered"] == "PASS"
    assert by_requirement["Public external melanoma ICB lead triage is complete"] == "PASS"
    assert by_requirement["Controlled external access package is ready for strict melanoma validation"] in {"PASS", "GAP"}
    assert by_requirement["GSE165745 NanoString panel-transfer sensitivity scoring is registered"] == "PASS"
    assert by_requirement["No-leakage rank-fusion ecological negative audit is registered"] == "PASS"
    assert by_requirement["Tumor-immune balance rescue head improves primary AUROC and strict external point estimate"] == "PASS"
    assert by_requirement["MAP4K1-TBX3 transform audit improves primary-selected score and current external stress test"] == "PASS"
    assert by_requirement["Eight-signature family comparison FDR q<0.05"] == "PASS"
    assert by_requirement["Strict melanoma external AUROC >=0.70 with FDR support"] == "PARTIAL"
    assert by_requirement["Interaction edges improve prediction or stability"] == "PASS"
    assert by_requirement["Biological objective improves prediction, calibration, or stability"] == "PASS"
    assert by_requirement["Training-only candidate search negative audit is registered"] == "PASS"


def test_top_tier_audit_accepts_nested_threshold_ba_only_for_operating_point():
    primary = pd.DataFrame(
        [
            {
                "endpoint": "primary_recist",
                "stratum": "melanoma_core_high_evidence",
                "model_name": "EcoNiche-Opt-HeuristicEcology",
                "n_samples": 100,
                "n_responders": 40,
                "pooled_AUROC": 0.705,
                "pooled_AUPRC": 0.56,
                "pooled_balanced_accuracy": 0.63,
            }
        ]
    )
    threshold_recalibration = pd.DataFrame(
        [
            {
                "endpoint": "primary_recist",
                "stratum": "melanoma_core_high_evidence",
                "threshold_policy": "nested_midrange_fixed_grid",
                "AUROC": 0.705,
                "AUPRC": 0.56,
                "balanced_accuracy": 0.651,
                "selected_policies": "fixed_0.40",
            }
        ]
    )
    audit = build_audit(
        primary,
        pd.DataFrame(),
        pd.DataFrame(),
        pd.DataFrame(),
        pd.DataFrame(),
        pd.DataFrame(),
        pd.DataFrame(),
        pd.DataFrame(),
        threshold_recalibration=threshold_recalibration,
    )
    by_requirement = audit.set_index("requirement")["status"].to_dict()

    assert by_requirement["Primary melanoma balanced accuracy >=0.65"] == "PASS"
    assert by_requirement["Primary melanoma LODO AUROC >=0.72"] == "PARTIAL"
    assert by_requirement["No-leakage threshold recalibration reaches primary balanced-accuracy target"] == "PASS"
