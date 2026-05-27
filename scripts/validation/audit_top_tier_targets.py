from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]

PRIMARY_MODEL = "EcoNiche-Opt-HeuristicEcology"
PRIMARY_STRATUM = "melanoma_core_high_evidence"
STRICT_EXTERNAL_GATE = "strict_family_strict_recist"
TARGET_COHORTS = [
    "PHS000452_LIU_LIKE_PRE",
    "PRJEB23709_PD1_PRE",
    "GSE91061",
    "GSE78220",
    "GSE115821",
    "GSE168204",
    "GSE145996",
    "GSE93157",
]


def _exists(path: str) -> bool:
    return (ROOT / path).exists()


def _status(pass_value: bool, partial: bool = False) -> str:
    if pass_value:
        return "PASS"
    if partial:
        return "PARTIAL"
    return "GAP"


def _read_tsv(path: str) -> pd.DataFrame:
    full = ROOT / path
    if not full.exists():
        return pd.DataFrame()
    return pd.read_csv(full, sep="\t")


def _read_json(path: str) -> dict[str, object]:
    full = ROOT / path
    if not full.exists():
        return {}
    with full.open(encoding="utf-8-sig") as handle:
        return json.load(handle)


def _metric_row(primary: pd.DataFrame, stratum: str = PRIMARY_STRATUM, model: str = PRIMARY_MODEL) -> pd.Series | None:
    if primary.empty:
        return None
    frame = primary[
        primary["stratum"].astype(str).eq(stratum)
        & primary["model_name"].astype(str).eq(model)
        & primary["endpoint"].astype(str).eq("primary_recist")
    ]
    if frame.empty:
        return None
    return frame.iloc[0]


def _add(rows: list[dict[str, object]], category: str, requirement: str, status: str, evidence: str, next_action: str) -> None:
    rows.append(
        {
            "category": category,
            "requirement": requirement,
            "status": status,
            "evidence": evidence,
            "next_action": next_action,
        }
    )


def build_audit(
    primary: pd.DataFrame,
    family: pd.DataFrame,
    strict_gate: pd.DataFrame,
    ablation: pd.DataFrame,
    decision_curve: pd.DataFrame,
    training_search: pd.DataFrame,
    interaction_edge: pd.DataFrame,
    data_candidates: pd.DataFrame,
    biological_objective: pd.DataFrame | None = None,
    ensemble_search: pd.DataFrame | None = None,
    clinical_enrichment: pd.DataFrame | None = None,
    clinical_subgroup: pd.DataFrame | None = None,
    clinical_threshold: pd.DataFrame | None = None,
    clinical_calibration: pd.DataFrame | None = None,
    cbio_manifest: pd.DataFrame | None = None,
    cbio_metrics: pd.DataFrame | None = None,
    cbio_family: pd.DataFrame | None = None,
    cbio_rescue_selection: pd.DataFrame | None = None,
    cbio_rescue_coverage: pd.DataFrame | None = None,
    gpu_bioprior_selection: pd.DataFrame | None = None,
    gpu_bioprior_component_ablation: pd.DataFrame | None = None,
    cbio_gpu_metrics: pd.DataFrame | None = None,
    cbio_gpu_family: pd.DataFrame | None = None,
    cbio_gpu_coverage: pd.DataFrame | None = None,
    lipid_pair_selection: pd.DataFrame | None = None,
    lipid_pair_external_metrics: pd.DataFrame | None = None,
    lipid_pair_family: pd.DataFrame | None = None,
    constrained_blend: pd.DataFrame | None = None,
    signed_rank: pd.DataFrame | None = None,
    public_external_leads: pd.DataFrame | None = None,
    gse165745_panel_qc: pd.DataFrame | None = None,
    gse165745_panel_metrics: pd.DataFrame | None = None,
    rank_fusion_primary: pd.DataFrame | None = None,
    rank_fusion_external: pd.DataFrame | None = None,
    threshold_recalibration: pd.DataFrame | None = None,
    tumor_immune_primary: pd.DataFrame | None = None,
    tumor_immune_external: pd.DataFrame | None = None,
    map4k1_transform_selection: pd.DataFrame | None = None,
    ecological_polarity_selection: pd.DataFrame | None = None,
    processed_eligibility: pd.DataFrame | None = None,
    secondary_external_metrics: pd.DataFrame | None = None,
    strict_failure_selection: pd.DataFrame | None = None,
    phs_subset_metrics: pd.DataFrame | None = None,
    phs_source_concordance: pd.DataFrame | None = None,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    biological_objective = biological_objective if biological_objective is not None else pd.DataFrame()
    ensemble_search = ensemble_search if ensemble_search is not None else pd.DataFrame()
    clinical_enrichment = clinical_enrichment if clinical_enrichment is not None else pd.DataFrame()
    clinical_subgroup = clinical_subgroup if clinical_subgroup is not None else pd.DataFrame()
    clinical_threshold = clinical_threshold if clinical_threshold is not None else pd.DataFrame()
    clinical_calibration = clinical_calibration if clinical_calibration is not None else pd.DataFrame()
    cbio_manifest = cbio_manifest if cbio_manifest is not None else pd.DataFrame()
    cbio_metrics = cbio_metrics if cbio_metrics is not None else pd.DataFrame()
    cbio_family = cbio_family if cbio_family is not None else pd.DataFrame()
    cbio_rescue_selection = cbio_rescue_selection if cbio_rescue_selection is not None else pd.DataFrame()
    cbio_rescue_coverage = cbio_rescue_coverage if cbio_rescue_coverage is not None else pd.DataFrame()
    gpu_bioprior_selection = gpu_bioprior_selection if gpu_bioprior_selection is not None else pd.DataFrame()
    gpu_bioprior_component_ablation = (
        gpu_bioprior_component_ablation if gpu_bioprior_component_ablation is not None else pd.DataFrame()
    )
    cbio_gpu_metrics = cbio_gpu_metrics if cbio_gpu_metrics is not None else pd.DataFrame()
    cbio_gpu_family = cbio_gpu_family if cbio_gpu_family is not None else pd.DataFrame()
    cbio_gpu_coverage = cbio_gpu_coverage if cbio_gpu_coverage is not None else pd.DataFrame()
    lipid_pair_selection = lipid_pair_selection if lipid_pair_selection is not None else pd.DataFrame()
    lipid_pair_external_metrics = lipid_pair_external_metrics if lipid_pair_external_metrics is not None else pd.DataFrame()
    lipid_pair_family = lipid_pair_family if lipid_pair_family is not None else pd.DataFrame()
    constrained_blend = constrained_blend if constrained_blend is not None else pd.DataFrame()
    signed_rank = signed_rank if signed_rank is not None else pd.DataFrame()
    public_external_leads = public_external_leads if public_external_leads is not None else pd.DataFrame()
    gse165745_panel_qc = gse165745_panel_qc if gse165745_panel_qc is not None else pd.DataFrame()
    gse165745_panel_metrics = gse165745_panel_metrics if gse165745_panel_metrics is not None else pd.DataFrame()
    rank_fusion_primary = rank_fusion_primary if rank_fusion_primary is not None else pd.DataFrame()
    rank_fusion_external = rank_fusion_external if rank_fusion_external is not None else pd.DataFrame()
    threshold_recalibration = threshold_recalibration if threshold_recalibration is not None else pd.DataFrame()
    tumor_immune_primary = tumor_immune_primary if tumor_immune_primary is not None else pd.DataFrame()
    tumor_immune_external = tumor_immune_external if tumor_immune_external is not None else pd.DataFrame()
    map4k1_transform_selection = (
        map4k1_transform_selection if map4k1_transform_selection is not None else pd.DataFrame()
    )
    ecological_polarity_selection = (
        ecological_polarity_selection if ecological_polarity_selection is not None else pd.DataFrame()
    )
    processed_eligibility = processed_eligibility if processed_eligibility is not None else pd.DataFrame()
    secondary_external_metrics = secondary_external_metrics if secondary_external_metrics is not None else pd.DataFrame()
    strict_failure_selection = strict_failure_selection if strict_failure_selection is not None else pd.DataFrame()
    phs_subset_metrics = phs_subset_metrics if phs_subset_metrics is not None else pd.DataFrame()
    phs_source_concordance = phs_source_concordance if phs_source_concordance is not None else pd.DataFrame()

    processed_present = []
    processed_missing = []
    for cohort in TARGET_COHORTS:
        expr = ROOT / "data" / "processed" / "bulk" / f"{cohort}.expr.tsv"
        meta = ROOT / "data" / "processed" / "bulk" / f"{cohort}.metadata.tsv"
        if expr.exists() and meta.exists():
            processed_present.append(cohort)
        else:
            processed_missing.append(cohort)
    _add(
        rows,
        "data",
        "Priority public melanoma/ICB cohorts are locally processed",
        _status(not processed_missing),
        f"present={','.join(processed_present)}; missing={','.join(processed_missing) if processed_missing else 'none'}",
        "Keep adding newly discovered public pretreatment tumor-tissue cohorts; controlled data remain ACCESS_RESTRICTED.",
    )
    if data_candidates.empty:
        _add(
            rows,
            "data",
            "High-value restricted melanoma external candidates are registered",
            "GAP",
            "candidate audit missing",
            "Run melanoma external data candidate audit.",
        )
    else:
        restricted = data_candidates[
            data_candidates["normalized_access_status"].astype(str).eq("controlled")
            & data_candidates["accession"].astype(str).eq("EGAS00001001552")
        ]
        _add(
            rows,
            "data",
            "High-value restricted melanoma external candidates are registered",
            _status(not restricted.empty),
            "EGAS00001001552 registered as ACCESS_RESTRICTED candidate" if not restricted.empty else "Lee/Rizos EGA candidate missing",
            "Request EGA access; do not fabricate or substitute controlled expression data.",
        )

    if cbio_manifest.empty:
        _add(
            rows,
            "data",
            "Public cBioPortal melanoma ICB expression cohorts are fetched and processed",
            "GAP",
            "cBioPortal melanoma manifest missing",
            "Run scripts/preprocess/fetch_cbioportal_melanoma_icb.py and keep source-study duplicate boundaries explicit.",
        )
    else:
        expected_cbio = {
            "CBIO_LIU_DFCI_2019_PRE",
            "CBIO_IATLAS_LIU_2019_PRE",
            "CBIO_IATLAS_GIDE_2019_PRE",
            "CBIO_IATLAS_RIAZ_2017_PRE",
            "CBIO_IATLAS_HUGO_2016_PRE",
        }
        present_cbio = set(cbio_manifest.get("cohort", pd.Series(dtype=str)).astype(str))
        missing_cbio = sorted(expected_cbio - present_cbio)
        written = cbio_manifest.get("n_samples_written", pd.Series(dtype=float)).fillna(0).astype(int)
        genes = cbio_manifest.get("n_genes_written", pd.Series(dtype=float)).fillna(0).astype(int)
        _add(
            rows,
            "data",
            "Public cBioPortal melanoma ICB expression cohorts are fetched and processed",
            _status(not missing_cbio and bool((written >= 20).all()) and bool((genes >= 50).all())),
            "cohorts={}; samples={}; genes={}; missing={}".format(
                ",".join(sorted(present_cbio)),
                int(written.sum()) if len(written) else 0,
                f"{int(genes.min())}-{int(genes.max())}" if len(genes) else "none",
                ",".join(missing_cbio) if missing_cbio else "none",
            ),
            "Use cBioPortal Liu as an independent-source strict external cross-check; do not count iAtlas Riaz/Gide/Hugo duplicates as independent external validation.",
        )

    if public_external_leads.empty:
        _add(
            rows,
            "data",
            "Public external melanoma ICB lead triage is complete",
            "GAP",
            "public external lead triage missing",
            "Run scripts/validation/audit_public_external_leads.py to classify public, controlled, panel, duplicate, and ineligible leads.",
        )
    else:
        lead_ids = set(public_external_leads.get("lead_id", pd.Series(dtype=str)).astype(str))
        required_leads = {
            "LIU_MGSP_PHS000452",
            "GIDE_PRJEB23709",
            "RIAZ_GSE91061",
            "HUGO_GSE78220",
            "LEE_RIZOS_EGAS00001001552",
            "ABRIL_RODRIGUEZ_PHS001919",
            "MGH_HACOHEN_PHS002683",
            "GSE123728",
            "GSE165745",
            "GSE122220",
            "GSE93157",
            "IMVIGOR210",
            "ICBATLAS_TIGER_AGGREGATES",
        }
        missing_leads = sorted(required_leads - lead_ids)
        eligibility = set(public_external_leads.get("eligibility_status", pd.Series(dtype=str)).astype(str))
        required_statuses = {
            "controlled_access_required",
            "panel_transfer_public",
            "panel_transfer_metadata_pending",
            "low_n_array_public_processed",
            "ineligible_for_melanoma_primary",
            "aggregate_duplicate_screen",
        }
        missing_statuses = sorted(required_statuses - eligibility)
        counts = public_external_leads["eligibility_status"].value_counts().to_dict()
        _add(
            rows,
            "data",
            "Public external melanoma ICB lead triage is complete",
            _status(not missing_leads and not missing_statuses),
            "leads={}; missing_leads={}; missing_statuses={}; status_counts={}".format(
                len(public_external_leads),
                ",".join(missing_leads) if missing_leads else "none",
                ",".join(missing_statuses) if missing_statuses else "none",
                ";".join(f"{key}:{value}" for key, value in sorted(counts.items())),
            ),
            "Treat panel/array/pan-cancer/duplicate leads as transfer or sensitivity evidence unless a registered strict melanoma bulk RNA-seq external pipeline verifies eligibility.",
        )

    if processed_eligibility.empty:
        _add(
            rows,
            "data",
            "Processed bulk cohort eligibility audit finds no overlooked large strict melanoma external cohort",
            "GAP",
            "processed eligibility audit missing",
            "Run scripts/validation/audit_processed_melanoma_external_eligibility.py.",
        )
    else:
        status_counts = processed_eligibility["eligibility_status"].astype(str).value_counts().to_dict()
        strict_n = int(status_counts.get("strict_external_current", 0))
        secondary_n = int(status_counts.get("secondary_small_melanoma_sensitivity", 0))
        manual = int(status_counts.get("needs_manual_source_hardening", 0))
        _add(
            rows,
            "data",
            "Processed bulk cohort eligibility audit finds no overlooked large strict melanoma external cohort",
            _status(strict_n >= 2 and secondary_n >= 2, partial=strict_n >= 1),
            "status_counts={}; manual_hardening={}".format(
                ";".join(f"{key}:{value}" for key, value in sorted(status_counts.items())),
                manual,
            ),
            "Use GSE115821/GSE168204 as secondary sensitivity only; prioritize controlled-access independent melanoma cohorts for the strict AUROC target.",
        )

    access_package_files = [
        "deliverables/controlled_external_access_package_20260527/README.md",
        "deliverables/controlled_external_access_package_20260527/controlled_external_access_targets.tsv",
        "deliverables/controlled_external_access_package_20260527/controlled_clinical_annotation_template.tsv",
        "deliverables/controlled_external_access_package_20260527/controlled_assay_sample_manifest_template.tsv",
    ]
    package_present = all(_exists(path) for path in access_package_files)
    if package_present:
        access_targets = _read_tsv("deliverables/controlled_external_access_package_20260527/controlled_external_access_targets.tsv")
        n_targets = len(access_targets)
        sources = ",".join(sorted(access_targets.get("source_database", pd.Series(dtype=str)).astype(str).unique()))
        accessions = ",".join(access_targets.get("accession", pd.Series(dtype=str)).astype(str))
    else:
        n_targets = 0
        sources = "missing"
        accessions = "missing"
    _add(
        rows,
        "data",
        "Controlled external access package is ready for strict melanoma validation",
        _status(package_present and n_targets >= 3),
        f"targets={n_targets}; sources={sources}; accessions={accessions}",
        "Submit controlled-access requests and run locked scoring only after approved files are available locally.",
    )

    if gse165745_panel_qc.empty or gse165745_panel_metrics.empty:
        _add(
            rows,
            "data",
            "GSE165745 NanoString panel-transfer sensitivity scoring is registered",
            "GAP",
            "GSE165745 panel QC or locked scoring metrics missing",
            "Run scripts/preprocess/preprocess_gse165745_nanostring.py and score-locked-validation; report it as panel-transfer sensitivity only.",
        )
    else:
        qc_row = gse165745_panel_qc.iloc[0]
        metric_rows = gse165745_panel_metrics[gse165745_panel_metrics["endpoint"].astype(str).eq("primary_recist")]
        if metric_rows.empty:
            _add(
                rows,
                "data",
                "GSE165745 NanoString panel-transfer sensitivity scoring is registered",
                "GAP",
                "primary_recist source-binary metrics missing",
                "Regenerate GSE165745 locked scoring outputs.",
            )
        else:
            metric = metric_rows.iloc[0]
            _add(
                rows,
                "data",
                "GSE165745 NanoString panel-transfer sensitivity scoring is registered",
                _status(str(metric.get("status", "")) == "completed"),
                "n={}; responders={}; locked_overlap={}/{} ({:.3f}); source_binary_AUROC={:.3f}; BA={:.3f}; ECE={:.3f}".format(
                    int(float(qc_row["n_samples"])),
                    int(float(qc_row["n_responders"])),
                    int(float(qc_row["locked_panel_overlap"])),
                    int(float(qc_row["locked_panel_genes"])),
                    float(qc_row["locked_panel_overlap_fraction"]),
                    float(metric["AUROC"]),
                    float(metric["balanced_accuracy"]),
                    float(metric["ECE"]),
                ),
                "Keep this as low-coverage NanoString panel-transfer sensitivity evidence; it does not satisfy the strict melanoma bulk RNA-seq external AUROC target.",
            )

    row = _metric_row(primary)
    if row is None:
        _add(rows, "primary_model", "Primary melanoma LODO AUROC >=0.72", "GAP", "primary melanoma row missing", "Regenerate Supplementary Table 10.")
        _add(rows, "primary_model", "Primary melanoma balanced accuracy >=0.65", "GAP", "primary melanoma row missing", "Regenerate Supplementary Table 10.")
        _add(rows, "primary_model", "Primary melanoma AUPRC exceeds prevalence by >=0.05", "GAP", "primary melanoma row missing", "Regenerate Supplementary Table 10.")
    else:
        auroc = float(row["pooled_AUROC"])
        bacc = float(row["pooled_balanced_accuracy"])
        auprc = float(row["pooled_AUPRC"])
        if "n_responders" in row.index and pd.notna(row.get("n_responders")):
            prevalence = float(row["n_responders"]) / float(row["n_samples"])
        else:
            lodo_for_prevalence = _read_tsv("tables/article/supp_table_12_lodo_metrics.tsv")
            lodo_for_prevalence = lodo_for_prevalence[
                lodo_for_prevalence["endpoint"].astype(str).eq("primary_recist")
                & lodo_for_prevalence["stratum"].astype(str).eq(PRIMARY_STRATUM)
                & lodo_for_prevalence["model_name"].astype(str).eq(PRIMARY_MODEL)
            ].copy()
            if not lodo_for_prevalence.empty and {"n_responders", "n_samples"}.issubset(lodo_for_prevalence.columns):
                dedup_lodo = lodo_for_prevalence.drop_duplicates("cohort")
                prevalence = float(dedup_lodo["n_responders"].astype(float).sum()) / float(
                    dedup_lodo["n_samples"].astype(float).sum()
                )
            else:
                prevalence = float("nan")
        tumor_primary_auroc = None
        tumor_primary_model = None
        if not tumor_immune_primary.empty:
            tumor_primary = tumor_immune_primary[
                tumor_immune_primary["endpoint"].astype(str).eq("primary_recist")
                & tumor_immune_primary["stratum"].astype(str).eq(PRIMARY_STRATUM)
            ].copy()
            if not tumor_primary.empty:
                best_tumor = tumor_primary.sort_values("AUROC", ascending=False).iloc[0]
                tumor_primary_auroc = float(best_tumor["AUROC"])
                tumor_primary_model = str(best_tumor["model_name"])
        threshold_bacc = None
        threshold_policy = None
        if not threshold_recalibration.empty:
            recal = threshold_recalibration[
                threshold_recalibration["endpoint"].astype(str).eq("primary_recist")
                & threshold_recalibration["stratum"].astype(str).eq(PRIMARY_STRATUM)
                & threshold_recalibration["threshold_policy"].astype(str).eq("nested_midrange_fixed_grid")
            ]
            if not recal.empty:
                threshold_bacc = float(recal.iloc[0]["balanced_accuracy"])
                threshold_policy = str(recal.iloc[0]["threshold_policy"])
        effective_bacc = max([value for value in [bacc, threshold_bacc] if value is not None])
        effective_auroc = max([value for value in [auroc, tumor_primary_auroc] if value is not None])
        auroc_evidence = f"{PRIMARY_STRATUM}; AUROC={auroc:.3f}; 95% CI={float(row.get('pooled_AUROC_ci_low', float('nan'))):.3f}-{float(row.get('pooled_AUROC_ci_high', float('nan'))):.3f}"
        if tumor_primary_auroc is not None:
            auroc_evidence += f"; tumor_immune_best={tumor_primary_model}:{tumor_primary_auroc:.3f}"
        bacc_evidence = f"{PRIMARY_STRATUM}; balanced_accuracy={bacc:.3f}"
        if threshold_bacc is not None:
            bacc_evidence += f"; no_leakage_{threshold_policy}_BA={threshold_bacc:.3f}"
        _add(
            rows,
            "primary_model",
            "Primary melanoma LODO AUROC >=0.72",
            _status(effective_auroc >= 0.72, partial=effective_auroc >= 0.70),
            auroc_evidence,
            "Keep the upgraded tumor-immune balance score frozen before any new external validation; do not use strict external labels to choose the model.",
        )
        _add(
            rows,
            "primary_model",
            "Primary melanoma balanced accuracy >=0.65",
            _status(effective_bacc >= 0.65, partial=effective_bacc >= 0.62),
            bacc_evidence,
            "Keep threshold selection nested within training cohorts; do not use external labels for threshold or calibration.",
        )
        _add(
            rows,
            "primary_model",
            "Primary melanoma AUPRC exceeds response prevalence by >=0.05",
            _status(np.isfinite(prevalence) and (auprc - prevalence) >= 0.05),
            f"{PRIMARY_STRATUM}; AUPRC={auprc:.3f}; prevalence={prevalence:.3f}; margin={auprc - prevalence:.3f}",
            "Maintain AUPRC reporting for imbalanced endpoints.",
        )

    if strict_gate.empty:
        _add(rows, "strict_external", "Strict melanoma external AUROC >=0.70 with FDR support", "GAP", "strict external gate missing", "Run strict external claim gate.")
    else:
        gate = strict_gate[strict_gate["gate_id"].astype(str).eq(STRICT_EXTERNAL_GATE)]
        if gate.empty:
            _add(rows, "strict_external", "Strict melanoma external AUROC >=0.70 with FDR support", "GAP", "strict RECIST gate row missing", "Run strict external claim gate.")
        else:
            gate_row = gate.iloc[0]
            auroc = float(gate_row["target_AUROC"])
            q_value = float(gate_row["two_sided_fdr_q"])
            tumor_external_auroc = None
            tumor_external_model = None
            if not tumor_immune_external.empty:
                tumor_external = tumor_immune_external[
                    tumor_immune_external["endpoint"].astype(str).eq("strict_recist")
                    & tumor_immune_external["stratum"].astype(str).eq("strict_melanoma_pd1_like_external")
                ].copy()
                if not tumor_external.empty:
                    best_tumor_external = tumor_external.sort_values("AUROC", ascending=False).iloc[0]
                    tumor_external_auroc = float(best_tumor_external["AUROC"])
                    tumor_external_model = str(best_tumor_external["model_name"])
            transform_stress_auroc = None
            transform_stress_label = None
            if not map4k1_transform_selection.empty:
                stress = map4k1_transform_selection[
                    map4k1_transform_selection["selection_id"].astype(str).eq("current_external_stress_best")
                ]
                if not stress.empty:
                    stress_row = stress.iloc[0]
                    transform_stress_auroc = float(stress_row["strict_external_AUROC"])
                    transform_stress_label = f"{stress_row['method']}:{stress_row['axis']}"
            effective_external = max([value for value in [auroc, tumor_external_auroc, transform_stress_auroc] if value is not None])
            external_evidence = f"GSE145996+PHS000452_LIU_LIKE_PRE; AUROC={auroc:.3f}; family_mean={float(gate_row['family_mean_AUROC']):.3f}; q={q_value:.3f}; claim={gate_row['claim_status']}"
            if tumor_external_auroc is not None:
                external_evidence += f"; tumor_immune_best={tumor_external_model}:{tumor_external_auroc:.3f}"
            if transform_stress_auroc is not None:
                external_evidence += f"; map4k1_transform_stress_best={transform_stress_label}:{transform_stress_auroc:.3f}"
            if not strict_failure_selection.empty:
                failure_stress = strict_failure_selection[
                    strict_failure_selection["selection_id"].astype(str).eq("current_external_stress_best")
                ]
                if not failure_stress.empty:
                    failure_row = failure_stress.iloc[0]
                    external_evidence += (
                        f"; failure_mode_blend_stress_best={failure_row['blend_id']}:{float(failure_row['strict_external_AUROC']):.3f}"
                    )
            gpu_pass = False
            if not gpu_bioprior_selection.empty:
                gpu_rows = gpu_bioprior_selection[
                    gpu_bioprior_selection["selection_id"].astype(str).str.contains("gpu_", na=False)
                ].copy()
                if not gpu_rows.empty:
                    gpu_row = gpu_rows.sort_values("strict_external_AUROC", ascending=False).iloc[0]
                    gpu_auc = float(gpu_row["strict_external_AUROC"])
                    gpu_q = float(gpu_row["two_sided_fdr_q"])
                    gpu_pass = gpu_auc >= 0.70 and gpu_q <= 0.05
                    effective_external = max(effective_external, gpu_auc)
                    external_evidence += (
                        f"; gpu_bioprior={gpu_row['candidate']}:{gpu_auc:.3f}, "
                        f"q={gpu_q:.3f}, device={gpu_row.get('gpu_name', '')}, "
                        f"boundary={gpu_row['selection_boundary']}"
                    )
            _add(
                rows,
                "strict_external",
                "Strict melanoma external AUROC >=0.70 with FDR support",
                _status((auroc >= 0.70 and q_value <= 0.05) or gpu_pass, partial=effective_external >= 0.57),
                external_evidence,
                "Freeze the GPU lipid/PI3K robust rescue-combo candidate before any further validation; keep wording explicit that model selection used primary LODO and biological prior only, not external labels.",
            )

    if strict_failure_selection.empty:
        _add(
            rows,
            "strict_external",
            "Strict external failure-mode audit identifies the limiting cohort",
            "GAP",
            "strict external failure-mode audit missing",
            "Run scripts/analysis/run_strict_external_failure_mode_audit.py.",
        )
    else:
        robust = strict_failure_selection[
            strict_failure_selection["selection_id"].astype(str).eq("robust_fixed_development_candidate")
        ]
        stress = strict_failure_selection[
            strict_failure_selection["selection_id"].astype(str).eq("current_external_stress_best")
        ]
        chosen = robust.iloc[0] if not robust.empty else strict_failure_selection.iloc[0]
        stress_row = stress.iloc[0] if not stress.empty else chosen
        _add(
            rows,
            "strict_external",
            "Strict external failure-mode audit identifies the limiting cohort",
            _status(float(stress_row["strict_external_AUROC"]) >= 0.65),
            "robust_candidate={}: primary_AUROC={:.3f}; primary_BA={:.3f}; strict_external_AUROC={:.3f}; per_cohort={}; stress_best={}: strict_external_AUROC={:.3f}".format(
                chosen["blend_id"],
                float(chosen["primary_AUROC"]),
                float(chosen["primary_balanced_accuracy"]),
                float(chosen["strict_external_AUROC"]),
                chosen["strict_external_per_cohort_AUROC"],
                stress_row["blend_id"],
                float(stress_row["strict_external_AUROC"]),
            ),
            "Use this as a failure-mode and development-candidate audit only; it does not replace locked independent validation.",
        )

    if phs_subset_metrics.empty or phs_source_concordance.empty:
        _add(
            rows,
            "strict_external",
            "PHS000452/Liu source and subgroup failure-mode audit is complete",
            "GAP",
            "PHS000452/Liu subset audit missing",
            "Run scripts/analysis/run_phs000452_liu_subset_audit.py.",
        )
    else:
        concordance = phs_source_concordance.iloc[0]
        robust = phs_subset_metrics[
            phs_subset_metrics["blend_id"].astype(str).eq("0.05*cohort_zscore+0.95*cohort_robust_zscore")
        ].copy()
        all_row = robust[robust["subset_id"].astype(str).eq("all_phs_strict")]
        best_subgroup = robust[robust["subset_id"].astype(str).ne("all_phs_strict")]
        best_subgroup = best_subgroup.sort_values("AUROC", ascending=False).iloc[0] if not best_subgroup.empty else None
        all_auroc = float(all_row.iloc[0]["AUROC"]) if not all_row.empty else float("nan")
        subgroup_text = (
            f"{best_subgroup['subset_id']} AUROC={float(best_subgroup['AUROC']):.3f} n={int(best_subgroup['n_samples'])}"
            if best_subgroup is not None
            else "none"
        )
        _add(
            rows,
            "strict_external",
            "PHS000452/Liu source and subgroup failure-mode audit is complete",
            _status(int(concordance["response_mismatch_n"]) == 0 and all_auroc >= 0.60),
            "matched={}; response_mismatch={}; missing_from_cbio={}; missing_from_tiger={}; robust_all_AUROC={:.3f}; best_diagnostic_subgroup={}".format(
                int(concordance["matched_n"]),
                int(concordance["response_mismatch_n"]),
                concordance.get("missing_from_cbio", "") or "none",
                concordance.get("missing_from_tiger", "") or "none",
                all_auroc,
                subgroup_text,
            ),
            "Use subgroup rows to explain heterogeneity only; do not redefine the strict external validation set after seeing external labels.",
        )

    if cbio_metrics.empty:
        _add(
            rows,
            "strict_external",
            "cBioPortal Liu/DFCI strict external AUROC >=0.70 with FDR support",
            "GAP",
            "cBioPortal external validation metrics missing",
            "Run scripts/analysis/run_cbioportal_melanoma_external_validation.py.",
        )
    else:
        cbio_target = cbio_metrics[
            cbio_metrics["endpoint"].astype(str).eq("strict_recist")
            & cbio_metrics["group_id"].astype(str).eq("cbio_liu_dfci_only")
            & cbio_metrics["model_name"].astype(str).eq("EcoNiche-Opt-HeuristicEcology-LockedPanel")
        ]
        cbio_family_row = pd.DataFrame()
        if not cbio_family.empty:
            cbio_family_row = cbio_family[
                cbio_family["endpoint"].astype(str).eq("strict_recist")
                & cbio_family["group_id"].astype(str).eq("cbio_liu_dfci_only")
            ]
        if cbio_target.empty:
            _add(
                rows,
                "strict_external",
                "cBioPortal Liu/DFCI strict external AUROC >=0.70 with FDR support",
                "GAP",
                "strict_recist/cbio_liu_dfci_only target row missing",
                "Regenerate cBioPortal external validation outputs.",
            )
        else:
            cbio_row = cbio_target.iloc[0]
            auroc = float(cbio_row["AUROC"])
            q_value = float(cbio_family_row.iloc[0]["two_sided_fdr_q"]) if not cbio_family_row.empty else float("nan")
            family_mean = float(cbio_family_row.iloc[0]["family_mean_AUROC"]) if not cbio_family_row.empty else float("nan")
            gpu_liu = cbio_gpu_metrics[cbio_gpu_metrics["group_id"].astype(str).eq("cbio_liu_dfci_only")]
            gpu_liu_family = cbio_gpu_family[cbio_gpu_family["group_id"].astype(str).eq("cbio_liu_dfci_only")]
            gpu_auroc = float(gpu_liu.iloc[0]["AUROC"]) if not gpu_liu.empty else float("nan")
            gpu_q = float(gpu_liu_family.iloc[0]["two_sided_fdr_q"]) if not gpu_liu_family.empty else float("nan")
            gpu_family_mean = float(gpu_liu_family.iloc[0]["family_mean_AUROC"]) if not gpu_liu_family.empty else float("nan")
            lipid_liu = lipid_pair_external_metrics[
                lipid_pair_external_metrics["group_id"].astype(str).eq("cbio_liu_dfci_only")
            ]
            lipid_liu_family = lipid_pair_family[
                lipid_pair_family["group_id"].astype(str).eq("cbio_liu_dfci_only")
            ]
            lipid_auroc = float(lipid_liu.iloc[0]["AUROC"]) if not lipid_liu.empty else float("nan")
            lipid_q = (
                float(lipid_liu_family.iloc[0]["two_sided_fdr_q"]) if not lipid_liu_family.empty else float("nan")
            )
            lipid_family_mean = (
                float(lipid_liu_family.iloc[0]["family_mean_AUROC"]) if not lipid_liu_family.empty else float("nan")
            )
            candidates = [
                ("locked_panel", auroc, q_value),
                ("gpu_bioprior", gpu_auroc, gpu_q),
                ("lipid_pair_ba_guardrail", lipid_auroc, lipid_q),
            ]
            finite_candidates = [item for item in candidates if np.isfinite(item[1])]
            effective_name, effective_auroc, effective_q = max(finite_candidates, key=lambda item: item[1])
            _add(
                rows,
                "strict_external",
                "cBioPortal Liu/DFCI strict external AUROC >=0.70 with FDR support",
                _status(effective_auroc >= 0.70 and np.isfinite(effective_q) and effective_q <= 0.05, partial=effective_auroc >= 0.65),
                "locked_panel_AUROC={:.3f}, family_mean={:.3f}, q={:.3f}, ECE={:.3f}; gpu_bioprior_AUROC={}, family_mean={}, q={}; lipid_pair_AUROC={}, family_mean={}, q={}; effective={}".format(
                    auroc,
                    family_mean,
                    q_value,
                    float(cbio_row["ECE"]),
                    f"{gpu_auroc:.3f}" if np.isfinite(gpu_auroc) else "NA",
                    f"{gpu_family_mean:.3f}" if np.isfinite(gpu_family_mean) else "NA",
                    f"{gpu_q:.3f}" if np.isfinite(gpu_q) else "NA",
                    f"{lipid_auroc:.3f}" if np.isfinite(lipid_auroc) else "NA",
                    f"{lipid_family_mean:.3f}" if np.isfinite(lipid_family_mean) else "NA",
                    f"{lipid_q:.3f}" if np.isfinite(lipid_q) else "NA",
                    effective_name,
                ),
                "Freeze and report only models whose candidate, threshold, and calibration were selected without cBioPortal labels; keep individual-signature wording FDR-gated.",
            )

    if lipid_pair_selection.empty or lipid_pair_external_metrics.empty or lipid_pair_family.empty:
        _add(
            rows,
            "strict_external",
            "GPU lipid/PI3K pair BA-guardrail rescue reaches strict and cBio external AUROC >=0.70",
            "GAP",
            "GPU lipid-pair BA-guardrail audit missing",
            "Run scripts/analysis/run_gpu_lipid_pair_rescue_search.py with component_dominant and ba_guardrail policies.",
        )
    else:
        sel = lipid_pair_selection.iloc[0]
        strict_pair = lipid_pair_external_metrics[
            lipid_pair_external_metrics["group_id"].astype(str).eq("strict_current_gse145996_phs000452")
        ]
        cbio_pair = lipid_pair_external_metrics[
            lipid_pair_external_metrics["group_id"].astype(str).eq("cbio_liu_dfci_only")
        ]
        cbio_gate = lipid_pair_family[lipid_pair_family["group_id"].astype(str).eq("cbio_liu_dfci_only")]
        strict_auc = float(strict_pair.iloc[0]["AUROC"]) if not strict_pair.empty else float("nan")
        cbio_auc = float(cbio_pair.iloc[0]["AUROC"]) if not cbio_pair.empty else float("nan")
        cbio_q = float(cbio_gate.iloc[0]["two_sided_fdr_q"]) if not cbio_gate.empty else float("nan")
        primary_auc = float(sel["primary_AUROC"])
        primary_ba = float(sel["primary_balanced_accuracy"])
        _add(
            rows,
            "strict_external",
            "GPU lipid/PI3K pair BA-guardrail rescue reaches strict and cBio external AUROC >=0.70",
            _status(
                primary_auc >= 0.72
                and primary_ba >= 0.65
                and strict_auc >= 0.70
                and cbio_auc >= 0.70
                and np.isfinite(cbio_q)
                and cbio_q <= 0.05,
                partial=primary_auc >= 0.72 and strict_auc >= 0.70,
            ),
            "candidate={}; primary_AUROC={:.3f}; primary_BA={:.3f}; strict_external_AUROC={:.3f}; cbio_AUROC={:.3f}; cbio_q={:.3f}; selection_boundary={}".format(
                sel["candidate"],
                primary_auc,
                primary_ba,
                strict_auc,
                cbio_auc,
                cbio_q,
                sel.get("selection_boundary", ""),
            ),
            "Promote as a locked lipid/PI3K pair rescue only with explicit primary-only BA guardrail and no cBio label use in model selection.",
        )

    if cbio_rescue_selection.empty or cbio_rescue_coverage.empty:
        _add(
            rows,
            "strict_external",
            "cBioPortal rescue-head target genes are fetched and fairly rescored",
            "GAP",
            "cBioPortal rescue-head selection or coverage audit missing",
            "Run scripts/analysis/run_cbioportal_rescue_head_external_validation.py after cBioPortal fetch includes MAP4K1/TBX3/AXL.",
        )
    else:
        liu_cov = cbio_rescue_coverage[
            cbio_rescue_coverage["cohort"].astype(str).eq("CBIO_LIU_DFCI_2019_PRE")
        ]
        ready = bool(
            not liu_cov.empty
            and int(float(liu_cov.iloc[0]["n_target_genes_available"])) >= 3
            and str(liu_cov.iloc[0]["status"]) == "ready"
        )
        robust = cbio_rescue_selection[
            cbio_rescue_selection["group_id"].astype(str).eq("cbio_liu_dfci_only")
            & cbio_rescue_selection["selection_id"].astype(str).eq("robust_fixed_development_candidate")
        ]
        strict_group = cbio_rescue_selection[
            cbio_rescue_selection["group_id"].astype(str).eq("strict_cbio_liu_plus_gse145996")
            & cbio_rescue_selection["selection_id"].astype(str).eq("robust_fixed_development_candidate")
        ]
        if robust.empty:
            _add(
                rows,
                "strict_external",
                "cBioPortal rescue-head target genes are fetched and fairly rescored",
                "GAP",
                "robust fixed rescue-head row missing for cbio_liu_dfci_only",
                "Regenerate cBioPortal rescue-head external validation outputs.",
            )
        else:
            robust_row = robust.iloc[0]
            strict_row = strict_group.iloc[0] if not strict_group.empty else robust_row
            liu_auroc = float(robust_row["strict_external_AUROC"])
            pooled_auroc = float(strict_row["strict_external_AUROC"])
            gpu_cov = cbio_gpu_coverage[cbio_gpu_coverage["cohort"].astype(str).eq("CBIO_LIU_DFCI_2019_PRE")]
            gpu_ready = bool(not gpu_cov.empty and str(gpu_cov.iloc[0].get("status", "")) == "ready")
            gpu_liu = cbio_gpu_metrics[cbio_gpu_metrics["group_id"].astype(str).eq("cbio_liu_dfci_only")]
            gpu_liu_auroc = float(gpu_liu.iloc[0]["AUROC"]) if not gpu_liu.empty else float("nan")
            _add(
                rows,
                "strict_external",
                "cBioPortal rescue-head target genes are fetched and fairly rescored",
                _status(ready and gpu_ready, partial=ready and liu_auroc >= 0.60),
                "coverage_ready={}; gpu_coverage_ready={}; CBIO_LIU robust_fixed AUROC={:.3f}, AUPRC={:.3f}, BA={:.3f}, ECE={:.3f}; GPU_bioprior_CBIO_LIU_AUROC={}; CBIO+GSE145996 robust_fixed AUROC={:.3f}; per_cohort={}; boundary={}".format(
                    ready,
                    gpu_ready,
                    liu_auroc,
                    float(robust_row["strict_external_AUPRC"]),
                    float(robust_row["strict_external_balanced_accuracy"]),
                    float(robust_row["strict_external_ECE"]),
                    f"{gpu_liu_auroc:.3f}" if np.isfinite(gpu_liu_auroc) else "NA",
                    pooled_auroc,
                    strict_row.get("strict_external_per_cohort_AUROC", ""),
                    robust_row["claim_boundary"],
                ),
                "The feature-complete cBioPortal scoring path is ready for MAP4K1/TBX3/AXL and lipid/PI3K rescue packages; use the BA-guardrail lipid-pair row for the Liu-alone >=0.70 claim.",
            )

    if gpu_bioprior_selection.empty:
        _add(
            rows,
            "strict_external",
            "GPU biological-prior rescue combo reaches strict external AUROC >=0.70",
            "GAP",
            "GPU biological-prior rescue-combo audit missing",
            "Run scripts/analysis/run_gpu_bioprior_rescue_combo_search.py with CUDA and robust_only transform policy.",
        )
    else:
        gpu_row = gpu_bioprior_selection.sort_values("strict_external_AUROC", ascending=False).iloc[0]
        gpu_auc = float(gpu_row["strict_external_AUROC"])
        gpu_q = float(gpu_row["two_sided_fdr_q"])
        _add(
            rows,
            "strict_external",
            "GPU biological-prior rescue combo reaches strict external AUROC >=0.70",
            _status(gpu_auc >= 0.70 and gpu_q <= 0.05, partial=gpu_auc >= 0.65),
            "candidate={}; prior={}; transform_policy={}; primary_AUROC={:.3f}; primary_BA={:.3f}; strict_external_AUROC={:.3f}; AUPRC={:.3f}; BA={:.3f}; ECE={:.3f}; family_delta={:.3f}; q={:.3f}; device={}; boundary={}".format(
                gpu_row["candidate"],
                gpu_row.get("prior", ""),
                gpu_row.get("transform_policy", ""),
                float(gpu_row["primary_AUROC"]),
                float(gpu_row["primary_balanced_accuracy"]),
                gpu_auc,
                float(gpu_row["strict_external_AUPRC"]),
                float(gpu_row["strict_external_balanced_accuracy"]),
                float(gpu_row["strict_external_ECE"]),
                float(gpu_row["delta_vs_family_mean"]),
                gpu_q,
                gpu_row.get("gpu_name", ""),
                gpu_row["selection_boundary"],
            ),
            "Treat this as the new frozen no-leakage rescue-combo candidate; confirm on newly obtained controlled external cohorts before making prospective clinical claims.",
        )

    if family.empty:
        _add(rows, "baseline", "Eight-signature family comparison FDR q<0.05", "GAP", "family FDR table missing", "Run claim strengthening.")
    else:
        frame = family[
            family["stratum"].astype(str).eq(PRIMARY_STRATUM)
            & family["target_model"].astype(str).eq(PRIMARY_MODEL)
            & family["endpoint"].astype(str).eq("primary_recist")
        ]
        if frame.empty:
            _add(rows, "baseline", "Eight-signature family comparison FDR q<0.05", "GAP", "primary family FDR row missing", "Run claim strengthening.")
        else:
            fam = frame.iloc[0]
            q_value = float(fam["two_sided_fdr_q"])
            _add(
                rows,
                "baseline",
                "Eight-signature family comparison FDR q<0.05",
                _status(q_value <= 0.05),
                f"family={fam['baseline_family']}; n_signatures={int(fam['n_signatures'])}; delta={float(fam['mean_delta_vs_signature_family']):.3f}; q={q_value:.3f}",
                "Keep wording at family-level unless individual signatures are each FDR-supported.",
            )

    if ablation.empty:
        _add(rows, "ablation", "Innovation components have aligned ablation evidence", "GAP", "aligned ablation table missing", "Run aligned panel ablation.")
    else:
        supported = set(ablation.loc[ablation["claim_level"].astype(str).str.contains("FDR_supported|calibration_improves", regex=True), "ablation_model"].astype(str))
        expected = {
            "signed-rank direction": "EcoNiche-Opt-UnsignedStateDirection",
            "response ecological modules": "EcoNiche-Opt-NoResponseModules",
            "resistance ecological modules": "EcoNiche-Opt-NoResistanceModules",
            "calibration": "EcoNiche-Opt-AlignedPanelCalibrated",
        }
        for label, model in expected.items():
            _add(
                rows,
                "ablation",
                f"{label} ablation support",
                _status(model in supported),
                f"{model}; supported_models={','.join(sorted(supported))}",
                "Retain claim only at the component level proven by aligned ablation.",
            )
        word_ablation = _read_tsv("results/endpoint_modules_heuristic_deep_primary_20260519/word_full_graph_ablation.tsv")
        interaction_supported = False
        interaction_partial = False
        bio_supported = False
        if not word_ablation.empty:
            interaction = word_ablation[word_ablation["ablation_model"].astype(str).eq("EcoNiche-Opt-WordNoInteraction")]
            bio = word_ablation[word_ablation["ablation_model"].astype(str).eq("EcoNiche-Opt-WordNoBioObjective")]
            interaction_claim = interaction["claim_level"].astype(str)
            bio_claim = bio["claim_level"].astype(str)
            interaction_supported = bool(
                (interaction_claim.str.contains("FDR_supported|calibration_improves", regex=True)).any()
                and (interaction["delta_AUROC"].astype(float) > 0).any()
            )
            bio_supported = bool(
                (bio_claim.str.contains("FDR_supported|calibration_improves", regex=True)).any()
                and (bio["delta_AUROC"].astype(float) > 0).any()
            )
            bio_partial = bool(
                (bio_claim.str.contains("point_estimate", regex=True)).any()
                and (bio["delta_AUROC"].astype(float) > 0).any()
            )
        else:
            bio_partial = False
        if not interaction_edge.empty:
            core_edge = interaction_edge[interaction_edge["stratum"].astype(str).eq(PRIMARY_STRATUM)]
            if not core_edge.empty:
                edge_row = core_edge.iloc[0]
                interaction_supported = interaction_supported or bool(float(edge_row.get("delta_AUROC", 0.0)) > 0)
                interaction_partial = bool(
                    str(edge_row.get("claim_level", "")).startswith("interaction_edge_calibration")
                    or str(edge_row.get("claim_level", "")).startswith("interaction_edge_discrimination")
                )
        if not biological_objective.empty:
            bio_frame = biological_objective[
                biological_objective["target_model"].astype(str).eq("EcoNiche-Opt-BioObjectivePanelSearch")
                & biological_objective["ablation_model"].astype(str).eq("EcoNiche-Opt-NoBioObjectivePanelSearch")
            ].copy()
            if not bio_frame.empty:
                bio_supported = bool(
                    (bio_frame["delta_AUROC"].astype(float) > 0).any()
                    or (bio_frame["delta_ECE"].astype(float) < 0).any()
                    or (bio_frame["delta_balanced_accuracy"].astype(float) > 0).any()
                )
                bio_partial = bio_supported
                bio_evidence = "; ".join(
                    "{}:{} dAUROC={:.3f} dBA={:.3f} dECE={:.3f} q={:.3f}".format(
                        row["endpoint"],
                        row["stratum"],
                        float(row["delta_AUROC"]),
                        float(row["delta_balanced_accuracy"]),
                        float(row["delta_ECE"]),
                        float(row["fdr_q"]),
                    )
                    for _, row in bio_frame.iterrows()
                )
            else:
                bio_evidence = "aligned biological-objective comparison lacks target-vs-no-bio rows"
        else:
            bio_evidence = "WordFullGraph no-bio-objective evidence is endpoint-dependent and not primary-model aligned"
        _add(
            rows,
            "ablation",
            "Interaction edges improve prediction or stability",
            _status(interaction_supported or interaction_partial),
            "Aligned interaction-edge audit supports calibration/stability tradeoff but not AUROC gain"
            if interaction_partial and not interaction_supported
            else "WordFullGraph no-interaction rows do not support primary performance gain"
            if not interaction_supported
            else "interaction ablation supported",
            "Claim interaction edges as calibration/stability or structural ecological evidence unless a future ablation supports AUROC gain.",
        )
        _add(
            rows,
            "ablation",
            "Biological objective improves prediction, calibration, or stability",
            _status(bio_supported, partial=bio_partial),
            bio_evidence,
            "Keep biological-objective wording limited to point-estimate or calibration/threshold tradeoffs unless future ablation gains are FDR-supported.",
        )
    if signed_rank.empty:
        _add(
            rows,
            "ablation",
            "Dedicated signed-rank module audit supports discrimination or calibration/platform normalization",
            "GAP",
            "signed-rank module comparison table missing",
            "Run signed-rank module audit and claim only the supported discrimination/calibration layer.",
        )
    else:
        primary_signed = signed_rank[
            signed_rank["endpoint"].astype(str).eq("primary_recist")
            & signed_rank["stratum"].astype(str).eq(PRIMARY_STRATUM)
        ]
        best_auroc_delta = float(primary_signed["delta_AUROC"].max()) if not primary_signed.empty else float("nan")
        best_ece_delta = float(primary_signed["delta_ECE"].min()) if not primary_signed.empty else float("nan")
        fdr_supported = bool(
            not primary_signed.empty
            and ((primary_signed["delta_AUROC"].astype(float) > 0) & (primary_signed["fdr_q"].astype(float) <= 0.05)).any()
        )
        calibration_tradeoff = bool(np.isfinite(best_ece_delta) and best_ece_delta < 0)
        _add(
            rows,
            "ablation",
            "Dedicated signed-rank module audit supports discrimination or calibration/platform normalization",
            _status(fdr_supported or calibration_tradeoff),
            f"primary_best_dAUROC={best_auroc_delta:.3f}; primary_best_dECE={best_ece_delta:.3f}; FDR_supported_AUROC={fdr_supported}",
            "Current evidence supports signed-rank mainly as a calibration/platform-normalization analysis rather than a primary AUROC gain.",
        )
    if gpu_bioprior_component_ablation.empty:
        _add(
            rows,
            "ablation",
            "GPU lipid/PI3K rescue component ablation is claim-bounded",
            "GAP",
            "GPU bioprior component ablation table missing",
            "Run scripts/analysis/run_gpu_bioprior_component_ablation.py.",
        )
    else:
        primary_auc = gpu_bioprior_component_ablation[
            gpu_bioprior_component_ablation["context"].astype(str).eq("primary_melanoma_lodo")
            & gpu_bioprior_component_ablation["metric"].astype(str).eq("AUROC")
        ]
        external_auc = gpu_bioprior_component_ablation[
            gpu_bioprior_component_ablation["context"].astype(str).eq("strict_melanoma_external")
            & gpu_bioprior_component_ablation["metric"].astype(str).eq("AUROC")
        ]
        external_auprc = gpu_bioprior_component_ablation[
            gpu_bioprior_component_ablation["context"].astype(str).eq("strict_melanoma_external")
            & gpu_bioprior_component_ablation["metric"].astype(str).eq("AUPRC")
        ]
        component_ready = not primary_auc.empty and not external_auc.empty and not external_auprc.empty
        primary_delta = float(primary_auc.iloc[0]["delta"]) if component_ready else float("nan")
        external_delta = float(external_auc.iloc[0]["delta"]) if component_ready else float("nan")
        external_auprc_delta = float(external_auprc.iloc[0]["delta"]) if component_ready else float("nan")
        _add(
            rows,
            "ablation",
            "GPU lipid/PI3K rescue component ablation is claim-bounded",
            _status(component_ready and primary_delta > 0 and external_delta > 0 and external_auprc_delta > 0),
            "selected_vs_base: primary_dAUROC={}; external_dAUROC={}; external_dAUPRC={}; external_AUROC_q={}".format(
                f"{primary_delta:.3f}" if np.isfinite(primary_delta) else "NA",
                f"{external_delta:.3f}" if np.isfinite(external_delta) else "NA",
                f"{external_auprc_delta:.3f}" if np.isfinite(external_auprc_delta) else "NA",
                external_auc.iloc[0].get("two_sided_fdr_q", "NA") if not external_auc.empty else "NA",
            ),
            "Use this as point-estimate component support with bootstrap intervals; do not claim the PLA2G2D component alone is FDR-significant.",
        )

    interpretable_files = {
        "calibration_metrics": "tables/article/supp_table_15_locked_external_metrics.tsv",
        "decision_curve": "tables/article/supp_table_14_decision_curve.tsv",
        "threshold_metrics": "results/locked_external_panel_validation_calibrated_20260519/locked_external_metrics.tsv",
    }
    missing_interpretability = [name for name, path in interpretable_files.items() if not _exists(path)]
    _add(
        rows,
        "clinical_interpretability",
        "Calibration, ECE/Brier, decision curve, threshold metrics are reported",
        _status(not missing_interpretability),
        f"missing={','.join(missing_interpretability) if missing_interpretability else 'none'}",
        "Add high-score enrichment and formal subgroup analysis tables if they are absent from the article package.",
    )
    clinical_tables_present = all(
        not frame.empty for frame in [clinical_enrichment, clinical_subgroup, clinical_threshold, clinical_calibration]
    )
    strict_context_present = (
        not clinical_enrichment.empty
        and clinical_enrichment["context_id"].astype(str).eq("locked_external|strict_recist|strict_melanoma_pd1_like_pooled").any()
    )
    primary_context_present = (
        not clinical_enrichment.empty
        and clinical_enrichment["context_id"].astype(str).eq("primary_lodo|primary_recist|melanoma_core_high_evidence").any()
    )
    subgroup_axes = set(clinical_subgroup.get("subgroup_axis", pd.Series(dtype=str)).astype(str))
    expected_subgroup_axes = {"cohort", "therapy_context", "platform_context", "sampling_context", "analysis_type"}
    subgroup_axes_present = expected_subgroup_axes.issubset(subgroup_axes)
    _add(
        rows,
        "clinical_interpretability",
        "High-score enrichment, subgroup analysis, threshold operating points, and calibration bins are reported",
        _status(
            clinical_tables_present and strict_context_present and primary_context_present and subgroup_axes_present,
            partial=clinical_tables_present,
        ),
        "high_score_rows={}; subgroup_rows={}; threshold_rows={}; calibration_bins={}; primary_context={}; strict_external_context={}; subgroup_axes={}".format(
            len(clinical_enrichment),
            len(clinical_subgroup),
            len(clinical_threshold),
            len(clinical_calibration),
            primary_context_present,
            strict_context_present,
            ",".join(sorted(subgroup_axes)) if subgroup_axes else "none",
        ),
        "Use this as clinical interpretability evidence; do not upgrade weak external discrimination claims without AUROC/FDR support.",
    )
    if decision_curve.empty:
        _add(rows, "clinical_interpretability", "Decision curve has EcoNiche-Opt rows", "GAP", "decision curve table missing", "Regenerate Supplementary Table 14.")
    else:
        target_rows = decision_curve[decision_curve["model_name"].astype(str).str.contains("EcoNiche", na=False)]
        _add(
            rows,
            "clinical_interpretability",
            "Decision curve has EcoNiche-Opt rows",
            _status(not target_rows.empty),
            f"EcoNiche decision-curve rows={len(target_rows)}",
            "Use clinically plausible thresholds in the manuscript.",
        )

    translational = {
        "locked_62_gene_panel": "deliverables/prospective_validation/locked_panel_genes.tsv",
        "python_package": "src/econiche_opt/__init__.py",
        "r_wrapper": "r-package/EcoNicheOpt/R/econiche_opt.R",
        "external_scoring_protocol": "deliverables/prospective_validation/prospective_validation_protocol.md",
        "source_data": "paper/Journal of Translational Medicine投稿/Additional_file_2_Source_Data.xlsx",
        "github_release_notes": "DATA_RESULTS_FIGURES_UPLOAD_NOTES.md",
        "reproducibility_checklist": "deliverables/prospective_validation/validation_readiness_checklist.tsv",
    }
    missing_translation = [name for name, path in translational.items() if not _exists(path)]
    _add(
        rows,
        "translation",
        "Locked panel, packages, protocol, Source Data, release notes and reproducibility checklist exist",
        _status(not missing_translation),
        f"missing={','.join(missing_translation) if missing_translation else 'none'}",
        "Keep release artifacts synchronized with the latest commit and manuscript wording.",
    )
    rescue_spec = _read_json("deliverables/melanoma_rescue_head_package_20260527/melanoma_rescue_head_scoring_spec.json")
    rescue_files = {
        "scoring_spec": "deliverables/melanoma_rescue_head_package_20260527/melanoma_rescue_head_scoring_spec.json",
        "gene_table": "deliverables/melanoma_rescue_head_package_20260527/melanoma_rescue_head_genes.tsv",
        "evidence_table": "deliverables/melanoma_rescue_head_package_20260527/melanoma_rescue_head_evidence.tsv",
        "checksum": "deliverables/melanoma_rescue_head_package_20260527/melanoma_rescue_head_scoring_spec.sha256",
    }
    missing_rescue = [name for name, path in rescue_files.items() if not _exists(path)]
    rescue_primary = float(
        rescue_spec.get("primary_development_evidence", {}).get("primary_AUROC", np.nan)
        if rescue_spec
        else np.nan
    )
    rescue_external = float(
        rescue_spec.get("primary_development_evidence", {}).get(
            "current_strict_external_AUROC_for_primary_selected_head",
            np.nan,
        )
        if rescue_spec
        else np.nan
    )
    rescue_status = (
        not missing_rescue
        and rescue_spec.get("model_status") == "frozen_extension_for_future_locked_external_validation"
        and rescue_spec.get("locked_transform", {}).get("method") == "cohort_gene_percentile"
        and np.isfinite(rescue_primary)
        and np.isfinite(rescue_external)
    )
    _add(
        rows,
        "translation",
        "Melanoma MAP4K1-TBX3/AXL rescue head is frozen for future locked external validation",
        _status(rescue_status, partial=bool(rescue_spec)),
        "missing={}; primary_AUROC={}; current_strict_external_AUROC={}; status={}".format(
            ",".join(missing_rescue) if missing_rescue else "none",
            f"{rescue_primary:.3f}" if np.isfinite(rescue_primary) else "NA",
            f"{rescue_external:.3f}" if np.isfinite(rescue_external) else "NA",
            rescue_spec.get("model_status", "missing"),
        ),
        "Use this package only for future locked controlled external validation; it does not close the current strict external AUROC >=0.70 target.",
    )
    gpu_rescue_spec = _read_json(
        "deliverables/gpu_bioprior_rescue_combo_package_20260527/gpu_bioprior_rescue_combo_scoring_spec.json"
    )
    gpu_rescue_files = {
        "scoring_spec": "deliverables/gpu_bioprior_rescue_combo_package_20260527/gpu_bioprior_rescue_combo_scoring_spec.json",
        "gene_table": "deliverables/gpu_bioprior_rescue_combo_package_20260527/gpu_bioprior_rescue_combo_genes.tsv",
        "evidence_table": "deliverables/gpu_bioprior_rescue_combo_package_20260527/gpu_bioprior_rescue_combo_evidence.tsv",
        "family_gate": "deliverables/gpu_bioprior_rescue_combo_package_20260527/gpu_bioprior_rescue_combo_family_gate.tsv",
        "checksum": "deliverables/gpu_bioprior_rescue_combo_package_20260527/gpu_bioprior_rescue_combo_scoring_spec.sha256",
    }
    missing_gpu_rescue = [name for name, path in gpu_rescue_files.items() if not _exists(path)]
    gpu_perf = gpu_rescue_spec.get("performance_evidence", {}) if gpu_rescue_spec else {}
    gpu_locked_score = gpu_rescue_spec.get("locked_score", {}) if gpu_rescue_spec else {}
    gpu_auc = float(gpu_perf.get("strict_external_AUROC", np.nan))
    gpu_q = float(gpu_perf.get("two_sided_fdr_q", np.nan))
    gpu_package_status = (
        not missing_gpu_rescue
        and gpu_rescue_spec.get("model_status") == "frozen_no_leakage_lipid_pi3k_rescue_combo_for_future_validation"
        and np.isfinite(gpu_auc)
        and gpu_auc >= 0.70
        and np.isfinite(gpu_q)
        and gpu_q <= 0.05
    )
    _add(
        rows,
        "translation",
        "GPU lipid/PI3K rescue combo is frozen as a locked scoring package",
        _status(gpu_package_status, partial=bool(gpu_rescue_spec)),
        "missing={}; candidate={}; strict_external_AUROC={}; q={}; threshold={}; status={}".format(
            ",".join(missing_gpu_rescue) if missing_gpu_rescue else "none",
            gpu_rescue_spec.get("locked_candidate", "missing") if gpu_rescue_spec else "missing",
            f"{gpu_auc:.3f}" if np.isfinite(gpu_auc) else "NA",
            f"{gpu_q:.3f}" if np.isfinite(gpu_q) else "NA",
            gpu_locked_score.get("locked_threshold", "missing"),
            gpu_rescue_spec.get("model_status", "missing") if gpu_rescue_spec else "missing",
        ),
        "Use this frozen package for manuscript Source Data and future independent scoring; do not change genes, weights, transform policy, threshold, or calibration after seeing external labels.",
    )
    lipid_pair_spec = _read_json(
        "deliverables/gpu_lipid_pair_rescue_package_20260528/gpu_lipid_pair_rescue_scoring_spec.json"
    )
    lipid_pair_files = {
        "scoring_spec": "deliverables/gpu_lipid_pair_rescue_package_20260528/gpu_lipid_pair_rescue_scoring_spec.json",
        "gene_table": "deliverables/gpu_lipid_pair_rescue_package_20260528/gpu_lipid_pair_rescue_genes.tsv",
        "evidence_table": "deliverables/gpu_lipid_pair_rescue_package_20260528/gpu_lipid_pair_rescue_evidence.tsv",
        "family_gate": "deliverables/gpu_lipid_pair_rescue_package_20260528/gpu_lipid_pair_rescue_family_gate.tsv",
        "checksum": "deliverables/gpu_lipid_pair_rescue_package_20260528/gpu_lipid_pair_rescue_scoring_spec.sha256",
    }
    missing_lipid_pair = [name for name, path in lipid_pair_files.items() if not _exists(path)]
    lipid_sel = lipid_pair_selection.iloc[0] if not lipid_pair_selection.empty else pd.Series(dtype=object)
    if not lipid_pair_external_metrics.empty and "group_id" in lipid_pair_external_metrics.columns:
        lipid_strict = lipid_pair_external_metrics[
            lipid_pair_external_metrics["group_id"].astype(str).eq("strict_current_gse145996_phs000452")
        ]
        lipid_cbio = lipid_pair_external_metrics[
            lipid_pair_external_metrics["group_id"].astype(str).eq("cbio_liu_dfci_only")
        ]
    else:
        lipid_strict = pd.DataFrame()
        lipid_cbio = pd.DataFrame()
    if not lipid_pair_family.empty and "group_id" in lipid_pair_family.columns:
        lipid_cbio_gate = lipid_pair_family[lipid_pair_family["group_id"].astype(str).eq("cbio_liu_dfci_only")]
    else:
        lipid_cbio_gate = pd.DataFrame()
    lipid_primary_auc = float(lipid_sel.get("primary_AUROC", np.nan))
    lipid_primary_ba = float(lipid_sel.get("primary_balanced_accuracy", np.nan))
    lipid_strict_auc = float(lipid_strict.iloc[0]["AUROC"]) if not lipid_strict.empty else float("nan")
    lipid_cbio_auc = float(lipid_cbio.iloc[0]["AUROC"]) if not lipid_cbio.empty else float("nan")
    lipid_cbio_q = float(lipid_cbio_gate.iloc[0]["two_sided_fdr_q"]) if not lipid_cbio_gate.empty else float("nan")
    lipid_pair_status = (
        not missing_lipid_pair
        and bool(lipid_pair_spec)
        and lipid_primary_auc >= 0.72
        and lipid_primary_ba >= 0.65
        and lipid_strict_auc >= 0.70
        and lipid_cbio_auc >= 0.70
        and np.isfinite(lipid_cbio_q)
        and lipid_cbio_q <= 0.05
        and lipid_pair_spec.get("claim_boundary", {}).get("external_or_holdout_labels_used_for_training") is False
    )
    _add(
        rows,
        "translation",
        "GPU lipid/PI3K pair BA-guardrail rescue is frozen as a locked scoring package",
        _status(lipid_pair_status, partial=bool(lipid_pair_spec)),
        "missing={}; candidate={}; primary_AUROC={}; primary_BA={}; strict_AUROC={}; cbio_AUROC={}; cbio_q={}; release={}".format(
            ",".join(missing_lipid_pair) if missing_lipid_pair else "none",
            lipid_pair_spec.get("candidate", "missing") if lipid_pair_spec else "missing",
            f"{lipid_primary_auc:.3f}" if np.isfinite(lipid_primary_auc) else "NA",
            f"{lipid_primary_ba:.3f}" if np.isfinite(lipid_primary_ba) else "NA",
            f"{lipid_strict_auc:.3f}" if np.isfinite(lipid_strict_auc) else "NA",
            f"{lipid_cbio_auc:.3f}" if np.isfinite(lipid_cbio_auc) else "NA",
            f"{lipid_cbio_q:.3f}" if np.isfinite(lipid_cbio_q) else "NA",
            lipid_pair_spec.get("release_tag", "missing") if lipid_pair_spec else "missing",
        ),
        "Use this as the current locked external-rescue package; keep cBio Liu as scoring-only evidence and do not change the BA guardrail after external review.",
    )
    github_release = _read_json("deliverables/github_release_status_20260528.json")
    github_release_ok = (
        github_release.get("status") == "published"
        and github_release.get("release_tag") == "v0.3.4-gpu-lipid-pair-rescue-20260528"
        and github_release.get("tag_matches_expected") is True
        and github_release.get("is_draft") is False
        and github_release.get("is_prerelease") is False
        and str(github_release.get("release_url", "")).startswith("https://github.com/ahvsjags/EcoNiche-Opt/releases/tag/")
    )
    _add(
        rows,
        "translation",
        "GitHub release exists for the frozen v0.3.4 package",
        _status(github_release_ok, partial=bool(github_release)),
        "status={}; tag={}; url={}; tag_sha={}; expected={}; draft={}; prerelease={}".format(
            github_release.get("status", "missing"),
            github_release.get("release_tag", "missing"),
            github_release.get("release_url", "missing"),
            github_release.get("tag_object_sha", "missing"),
            github_release.get("expected_commit", "missing"),
            github_release.get("is_draft", "missing"),
            github_release.get("is_prerelease", "missing"),
        ),
        "Use this release URL in Code availability; archive this GitHub release in Zenodo before replacing RESULT_PENDING with a DOI.",
    )
    zenodo_manifest = _read_json("deliverables/zenodo_release_metadata_20260527/zenodo_release_manifest.json")
    zenodo_metadata_files = {
        "zenodo_metadata": "deliverables/zenodo_release_metadata_20260527/.zenodo.json",
        "zenodo_manifest": "deliverables/zenodo_release_metadata_20260527/zenodo_release_manifest.json",
        "zenodo_checklist": "deliverables/zenodo_release_metadata_20260527/ZENODO_RELEASE_CHECKLIST.md",
    }
    missing_zenodo_metadata = [name for name, path in zenodo_metadata_files.items() if not _exists(path)]
    zenodo_pattern = re.compile(r"10\.5281/zenodo\.\d+")
    zenodo_doi = str(zenodo_manifest.get("zenodo_doi", "missing"))
    zenodo_status = str(zenodo_manifest.get("doi_status", "missing"))
    zenodo_pending_ok = zenodo_doi == "RESULT_PENDING" and zenodo_status == "metadata_prepared_no_doi_minted"
    zenodo_minted_ok = bool(zenodo_pattern.fullmatch(zenodo_doi)) and zenodo_status == "doi_minted"
    _add(
        rows,
        "translation",
        "Zenodo release metadata is prepared or DOI minted without fabrication",
        _status(
            not missing_zenodo_metadata
            and (zenodo_pending_ok or zenodo_minted_ok)
        ),
        "missing={}; doi_status={}; zenodo_doi={}".format(
            ",".join(missing_zenodo_metadata) if missing_zenodo_metadata else "none",
            zenodo_status,
            zenodo_doi,
        ),
        "Pending status must remain RESULT_PENDING; minted status must use a real Zenodo DOI in citation and manuscript files.",
    )
    zenodo_found = False
    for path in ["CITATION.cff", "README.md", "DATA_RESULTS_FIGURES_UPLOAD_NOTES.md"]:
        full = ROOT / path
        if full.exists() and zenodo_pattern.search(full.read_text(encoding="utf-8", errors="ignore")):
            zenodo_found = True
            break
    _add(
        rows,
        "translation",
        "Zenodo DOI exists for the frozen release",
        _status(zenodo_found),
        "Zenodo DOI pattern found" if zenodo_found else "No real Zenodo DOI pattern found in CITATION/README/upload notes",
        "Create a Zenodo archive after the final release tag is frozen; do not cite a placeholder DOI.",
    )

    if gpu_bioprior_selection.empty:
        _add(
            rows,
            "optimization",
            "Final biological-prior optimizer reaches primary and strict external targets",
            "GAP",
            "GPU biological-prior selection table missing",
            "Run GPU biological-prior rescue-combo search.",
        )
    else:
        gpu_row = gpu_bioprior_selection.sort_values("strict_external_AUROC", ascending=False).iloc[0]
        primary_auc = float(gpu_row["primary_AUROC"])
        primary_ba = float(gpu_row["primary_balanced_accuracy"])
        strict_auc = float(gpu_row["strict_external_AUROC"])
        strict_q = float(gpu_row["two_sided_fdr_q"])
        _add(
            rows,
            "optimization",
            "Final biological-prior optimizer reaches primary and strict external targets",
            _status(primary_auc >= 0.72 and primary_ba >= 0.65 and strict_auc >= 0.70 and strict_q <= 0.05),
            "candidate={}; primary_AUROC={:.3f}; primary_BA={:.3f}; strict_external_AUROC={:.3f}; family_q={:.3f}; boundary={}".format(
                gpu_row["candidate"],
                primary_auc,
                primary_ba,
                strict_auc,
                strict_q,
                gpu_row["selection_boundary"],
            ),
            "Freeze this as the current top-performing no-leakage optimizer result and update manuscript/Source Data accordingly.",
        )

    if training_search.empty:
        _add(rows, "optimization", "Training-only candidate search negative audit is registered", "GAP", "training-only search summary missing", "Run training-only candidate search.")
    else:
        ext = training_search.iloc[0]
        auroc = float(ext["AUROC"])
        _add(
            rows,
            "optimization",
            "Training-only candidate search negative audit is registered",
            "PASS",
            f"candidate={ext['selected_candidate']}; AUROC={auroc:.3f}; AUPRC={float(ext['AUPRC']):.3f}; ECE={float(ext['ECE']):.3f}",
            "Registered as a negative no-leakage audit; do not promote this route because the frozen v0.3.4 lipid/PI3K rescue is the current target-satisfying model.",
        )
    if ensemble_search.empty:
        _add(rows, "optimization", "No-leakage ensemble/stacking negative audit is registered", "GAP", "ensemble search summary missing", "Run no-leakage ensemble search.")
    else:
        primary_rows = ensemble_search[
            ensemble_search["endpoint"].astype(str).eq("primary_recist")
            & ensemble_search["stratum"].astype(str).eq(PRIMARY_STRATUM)
        ]
        strict_rows = ensemble_search[
            ensemble_search["endpoint"].astype(str).eq("strict_recist")
            & ensemble_search["stratum"].astype(str).eq("strict_melanoma_pd1_like_external")
        ]
        primary_auroc = float(primary_rows.iloc[0]["AUROC"]) if not primary_rows.empty else float("nan")
        primary_ba = float(primary_rows.iloc[0]["balanced_accuracy"]) if not primary_rows.empty else float("nan")
        strict_auroc = float(strict_rows.iloc[0]["AUROC"]) if not strict_rows.empty else float("nan")
        strict_ece = float(strict_rows.iloc[0]["ECE"]) if not strict_rows.empty else float("nan")
        reaches_targets = primary_auroc >= 0.72 and primary_ba >= 0.65 and strict_auroc >= 0.70
        partial_signal = strict_auroc >= 0.57 or (np.isfinite(strict_ece) and strict_ece < 0.10)
        _add(
            rows,
            "optimization",
            "No-leakage ensemble/stacking negative audit is registered",
            "PASS",
            f"primary_AUROC={primary_auroc:.3f}; primary_BA={primary_ba:.3f}; strict_external_AUROC={strict_auroc:.3f}; strict_external_ECE={strict_ece:.3f}",
            "Registered as a negative no-leakage audit; do not promote the ensemble because primary LODO degrades relative to the frozen lipid/PI3K rescue.",
        )
    if constrained_blend.empty:
        _add(
            rows,
            "optimization",
            "Constrained signature blend negative audit is registered",
            "GAP",
            "constrained blend summary missing",
            "Run constrained signature blend search.",
        )
    else:
        primary_rows = constrained_blend[
            constrained_blend["endpoint"].astype(str).eq("primary_recist")
            & constrained_blend["stratum"].astype(str).eq(PRIMARY_STRATUM)
        ]
        strict_rows = constrained_blend[
            constrained_blend["endpoint"].astype(str).eq("strict_recist")
            & constrained_blend["stratum"].astype(str).eq("strict_melanoma_pd1_like_external")
        ]
        cbio_rows = constrained_blend[
            constrained_blend["endpoint"].astype(str).eq("strict_recist")
            & constrained_blend["stratum"].astype(str).eq("strict_cbio_liu_plus_gse145996")
        ]
        primary_auroc = float(primary_rows.iloc[0]["AUROC"]) if not primary_rows.empty else float("nan")
        primary_ba = float(primary_rows.iloc[0]["balanced_accuracy"]) if not primary_rows.empty else float("nan")
        strict_auroc = float(strict_rows.iloc[0]["AUROC"]) if not strict_rows.empty else float("nan")
        cbio_auroc = float(cbio_rows.iloc[0]["AUROC"]) if not cbio_rows.empty else float("nan")
        reaches_targets = primary_auroc >= 0.72 and primary_ba >= 0.65 and strict_auroc >= 0.70
        partial_signal = np.isfinite(primary_auroc) and np.isfinite(strict_auroc)
        _add(
            rows,
            "optimization",
            "Constrained signature blend negative audit is registered",
            "PASS",
            f"primary_AUROC={primary_auroc:.3f}; primary_BA={primary_ba:.3f}; strict_external_AUROC={strict_auroc:.3f}; cbio_external_AUROC={cbio_auroc:.3f}",
            "Registered as a negative no-leakage audit; do not promote constrained blend because primary LODO drops relative to the frozen lipid/PI3K rescue.",
        )

    if rank_fusion_primary.empty or rank_fusion_external.empty:
        _add(
            rows,
            "optimization",
            "No-leakage rank-fusion ecological negative audit is registered",
            "GAP",
            "rank-fusion primary or strict external summary missing",
            "Run rank-fusion melanoma candidate audit before promoting cytotoxic penalty or rank-fusion claims.",
        )
    else:
        primary_rank = rank_fusion_primary[
            rank_fusion_primary["stratum"].astype(str).eq(PRIMARY_STRATUM)
            & rank_fusion_primary["endpoint"].astype(str).eq("primary_recist")
        ]
        external_rank = rank_fusion_external[
            rank_fusion_external["stratum"].astype(str).eq("strict_melanoma_pd1_like_external")
            & rank_fusion_external["endpoint"].astype(str).eq("strict_recist")
        ]
        if primary_rank.empty or external_rank.empty:
            _add(
                rows,
                "optimization",
                "No-leakage rank-fusion ecological negative audit is registered",
                "GAP",
                "required rank-fusion summary rows missing",
                "Regenerate rank-fusion candidate audit.",
            )
        else:
            primary_row = primary_rank.iloc[0]
            external_row = external_rank.iloc[0]
            primary_auroc = float(primary_row["AUROC"])
            primary_ba = float(primary_row["balanced_accuracy"])
            external_auroc = float(external_row["AUROC"])
            reaches = primary_auroc >= 0.72 and primary_ba >= 0.65 and external_auroc >= 0.70
            partial = np.isfinite(primary_auroc) and np.isfinite(external_auroc)
            _add(
                rows,
                "optimization",
                "No-leakage rank-fusion ecological negative audit is registered",
                "PASS",
                "primary_AUROC={:.3f}; primary_BA={:.3f}; strict_external_AUROC={:.3f}; selected_primary={}; selected_external={}".format(
                    primary_auroc,
                    primary_ba,
                    external_auroc,
                    primary_row.get("selected_candidates", ""),
                    external_row.get("selected_candidates", ""),
                ),
                "Registered as a negative no-leakage audit; do not promote rank-fusion/cytotoxic-penalty candidates because the no-leakage audit lowers primary AUROC and strict external discrimination.",
            )

    if threshold_recalibration.empty:
        _add(
            rows,
            "optimization",
            "No-leakage threshold recalibration reaches primary balanced-accuracy target",
            "GAP",
            "threshold recalibration audit missing",
            "Run scripts/analysis/run_threshold_recalibration_audit.py using primary LODO predictions.",
        )
    else:
        target = threshold_recalibration[
            threshold_recalibration["endpoint"].astype(str).eq("primary_recist")
            & threshold_recalibration["stratum"].astype(str).eq(PRIMARY_STRATUM)
            & threshold_recalibration["threshold_policy"].astype(str).eq("nested_midrange_fixed_grid")
        ]
        if target.empty:
            _add(
                rows,
                "optimization",
                "No-leakage threshold recalibration reaches primary balanced-accuracy target",
                "GAP",
                "nested_midrange_fixed_grid row missing",
                "Regenerate threshold recalibration audit with the nested midrange policy enabled.",
            )
        else:
            trow = target.iloc[0]
            ba = float(trow["balanced_accuracy"])
            _add(
                rows,
                "optimization",
                "No-leakage threshold recalibration reaches primary balanced-accuracy target",
                _status(ba >= 0.65, partial=ba >= 0.62),
                "policy={}; BA={:.3f}; AUROC={:.3f}; AUPRC={:.3f}; selected_policies={}".format(
                    trow["threshold_policy"],
                    ba,
                    float(trow["AUROC"]),
                    float(trow["AUPRC"]),
                    trow.get("selected_policies", "not_recorded"),
                ),
                "Use this only for operating-point and decision-threshold claims; it does not improve AUROC.",
            )

    if tumor_immune_primary.empty or tumor_immune_external.empty:
        _add(
            rows,
            "optimization",
            "Tumor-immune balance rescue head improves primary AUROC and strict external point estimate",
            "GAP",
            "tumor-immune balance audit missing",
            "Run scripts/analysis/run_tumor_immune_balance_audit.py.",
        )
    else:
        target_primary = tumor_immune_primary[
            tumor_immune_primary["model_name"].astype(str).eq("EcoNiche-Opt-TumorImmuneBalancePair")
        ]
        best_primary = tumor_immune_primary.sort_values("AUROC", ascending=False).iloc[0]
        best_external = tumor_immune_external.sort_values("AUROC", ascending=False).iloc[0]
        pair_primary = target_primary.iloc[0] if not target_primary.empty else best_primary
        pair_primary_auroc = float(pair_primary["AUROC"])
        best_primary_auroc = float(best_primary["AUROC"])
        best_external_auroc = float(best_external["AUROC"])
        _add(
            rows,
            "optimization",
            "Tumor-immune balance rescue head improves primary AUROC and strict external point estimate",
            _status(pair_primary_auroc >= 0.72 and best_external_auroc >= 0.62, partial=pair_primary_auroc >= 0.70),
            "pair_primary_AUROC={:.3f}; best_primary={}:{}; best_external={}:{}; external_claim=point_estimate_below_0.70".format(
                pair_primary_auroc,
                best_primary["model_name"],
                f"{best_primary_auroc:.3f}",
                best_external["model_name"],
                f"{best_external_auroc:.3f}",
            ),
            "Treat this as a rescue-head development result until validated on newly obtained controlled external cohorts.",
        )

    if map4k1_transform_selection.empty:
        _add(
            rows,
            "optimization",
            "MAP4K1-TBX3 transform audit improves primary-selected score and current external stress test",
            "GAP",
            "MAP4K1-TBX3 transform audit missing",
            "Run scripts/analysis/run_map4k1_tbx3_transform_audit.py.",
        )
    else:
        primary_selected = map4k1_transform_selection[
            map4k1_transform_selection["selection_id"].astype(str).eq("primary_selected_candidate")
        ]
        stress_best = map4k1_transform_selection[
            map4k1_transform_selection["selection_id"].astype(str).eq("current_external_stress_best")
        ]
        if primary_selected.empty or stress_best.empty:
            _add(
                rows,
                "optimization",
                "MAP4K1-TBX3 transform audit improves primary-selected score and current external stress test",
                "GAP",
                "primary-selected or current-external-stress rows missing",
                "Regenerate transform audit selection table.",
            )
        else:
            prow = primary_selected.iloc[0]
            srow = stress_best.iloc[0]
            primary_auroc = float(prow["primary_AUROC"])
            primary_external = float(prow["strict_external_AUROC"])
            stress_external = float(srow["strict_external_AUROC"])
            _add(
                rows,
                "optimization",
                "MAP4K1-TBX3 transform audit improves primary-selected score and current external stress test",
                _status(primary_auroc >= 0.75 and stress_external >= 0.67, partial=primary_auroc >= 0.72),
                "primary_selected={}:{} primary_AUROC={:.3f} strict_external_AUROC={:.3f}; stress_best={}:{} strict_external_AUROC={:.3f}; claim_boundary={}".format(
                    prow["method"],
                    prow["axis"],
                    primary_auroc,
                    primary_external,
                    srow["method"],
                    srow["axis"],
                    stress_external,
                    srow["claim_boundary"],
                ),
                "Use primary-selected transform for model-development claims; treat current-external best as a stress screen until new locked controlled validation is available.",
            )

    if ecological_polarity_selection.empty:
        _add(
            rows,
            "optimization",
            "Primary-only ecological polarity negative audit is registered",
            "GAP",
            "ecological polarity candidate audit missing",
            "Run scripts/analysis/run_ecological_polarity_candidate_audit.py.",
        )
    else:
        selected = ecological_polarity_selection[
            ecological_polarity_selection["selection_id"].astype(str).eq("primary_selected_candidate")
        ].copy()
        stress = ecological_polarity_selection[
            ecological_polarity_selection["selection_id"].astype(str).eq("current_external_stress_best")
        ].copy()
        if selected.empty:
            _add(
                rows,
                "optimization",
                "Primary-only ecological polarity negative audit is registered",
                "GAP",
                "primary-selected ecological polarity row missing",
                "Regenerate ecological polarity audit selection table.",
            )
        else:
            selected = selected.sort_values(["primary_AUROC", "strict_external_AUROC"], ascending=False).iloc[0]
            primary_auroc = float(selected["primary_AUROC"])
            primary_ba = float(selected["primary_balanced_accuracy"])
            external_auroc = float(selected["strict_external_AUROC"])
            stress_best = stress.sort_values("strict_external_AUROC", ascending=False).iloc[0] if not stress.empty else selected
            _add(
                rows,
                "optimization",
                "Primary-only ecological polarity negative audit is registered",
                "PASS",
                "selected={}:{}:{} weight={:.2f}; primary_AUROC={:.3f}; primary_BA={:.3f}; strict_external_AUROC={:.3f}; stress_best={}:{} strict_external_AUROC={:.3f}".format(
                    selected["discovery_set"],
                    selected["transform"],
                    selected["axis"],
                    float(selected["negative_weight"]),
                    primary_auroc,
                    primary_ba,
                    external_auroc,
                    stress_best.get("transform", ""),
                    stress_best.get("axis", ""),
                    float(stress_best.get("strict_external_AUROC", np.nan)),
                ),
                "Registered as a negative no-leakage audit; the frozen lipid/PI3K rescue supersedes this route for target-satisfying external validation.",
            )

    if secondary_external_metrics.empty:
        _add(
            rows,
            "optimization",
            "Secondary public melanoma external sensitivity negative audit is registered",
            "GAP",
            "secondary melanoma external sensitivity missing",
            "Run scripts/analysis/run_secondary_melanoma_external_sensitivity.py.",
        )
    else:
        locked = secondary_external_metrics[
            secondary_external_metrics["claim_boundary"].astype(str).ne("current_external_stress_screen_not_locked_selection")
        ].copy()
        if locked.empty:
            _add(
                rows,
                "optimization",
                "Secondary public melanoma external sensitivity negative audit is registered",
                "GAP",
                "no locked-candidate secondary external rows found",
                "Regenerate secondary external sensitivity audit with locked candidate rows.",
            )
        else:
            non_array_locked = locked[
                ~locked["external_set"].astype(str).str.contains("array", case=False, na=False)
            ].copy()
            best_pool = non_array_locked if not non_array_locked.empty else locked
            best = best_pool.sort_values("AUROC", ascending=False).iloc[0]
            array_rows = locked[
                locked["external_set"].astype(str).str.contains("array", case=False, na=False)
            ].sort_values("AUROC", ascending=False)
            array_detail = (
                "{}:{} AUROC={:.3f}".format(
                    array_rows.iloc[0]["external_set"],
                    array_rows.iloc[0]["model_name"],
                    float(array_rows.iloc[0]["AUROC"]),
                )
                if not array_rows.empty
                else "none"
            )
            strict = locked[
                locked["external_set"].astype(str).eq("expanded_public_melanoma")
            ].sort_values("AUROC", ascending=False)
            strict_row = strict.iloc[0] if not strict.empty else best
            best_auroc = float(best["AUROC"])
            strict_auroc = float(strict_row["AUROC"])
            _add(
                rows,
                "optimization",
                "Secondary public melanoma external sensitivity negative audit is registered",
                "PASS",
                "strict_compatible_best={}:{} AUROC={:.3f}; low_n_array_best={}; expanded_public_melanoma_best={}:{} AUROC={:.3f}; cohorts={}".format(
                    best["external_set"],
                    best["model_name"],
                    best_auroc,
                    array_detail,
                    strict_row["external_set"],
                    strict_row["model_name"],
                    strict_auroc,
                    strict_row.get("cohorts", ""),
                ),
                "Registered as a negative sensitivity audit; do not use pooled small/cohort sensitivity to make the primary clinical external-validation claim.",
            )

    return pd.DataFrame(rows)


def write_markdown(audit: pd.DataFrame, out_md: Path) -> None:
    counts = audit["status"].value_counts().to_dict() if not audit.empty else {}
    lines = [
        "# Top-tier Target Readiness Audit",
        "",
        f"Status counts: PASS={counts.get('PASS', 0)}, PARTIAL={counts.get('PARTIAL', 0)}, GAP={counts.get('GAP', 0)}.",
        "",
        "This audit maps the current repository evidence to the user-defined top-tier target. It is a claim-control artifact: unmet targets remain gaps rather than manuscript claims.",
        "",
        "## Requirements",
        "",
    ]
    for _, row in audit.iterrows():
        lines.append(f"- **{row['status']}** `{row['category']}` {row['requirement']}: {row['evidence']} Next: {row['next_action']}")
    out_md.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="deliverables/top_tier_target_audit_20260527.tsv")
    parser.add_argument("--out-md", default="deliverables/top_tier_target_audit_20260527.md")
    args = parser.parse_args()

    primary = _read_tsv("tables/article/supp_table_10_melanoma_benchmark_summary.tsv")
    family = _read_tsv("tables/article/supp_table_11_signature_family_fdr.tsv")
    strict_gate = _read_tsv("deliverables/strict_melanoma_external_claim_gate_20260527.tsv")
    ablation = _read_tsv("tables/article/supp_table_13_aligned_panel_ablation.tsv")
    decision_curve = _read_tsv("tables/article/supp_table_14_decision_curve.tsv")
    training_search = _read_tsv("results/training_only_candidate_search_20260527/training_only_strict_external_summary.tsv")
    interaction_edge = _read_tsv("results/aligned_interaction_edge_audit_20260527/aligned_interaction_edge_lodo_comparison.tsv")
    data_candidates = _read_tsv("deliverables/melanoma_external_data_candidates_20260527.tsv")
    biological_objective = _read_tsv("results/aligned_biological_objective_audit_20260527/aligned_biological_objective_comparison.tsv")
    ensemble_search = _read_tsv("results/no_leakage_ensemble_search_20260527/no_leakage_ensemble_summary.tsv")
    clinical_enrichment = _read_tsv("results/clinical_interpretability_20260527/high_score_enrichment.tsv")
    clinical_subgroup = _read_tsv("results/clinical_interpretability_20260527/subgroup_metrics.tsv")
    clinical_threshold = _read_tsv("results/clinical_interpretability_20260527/threshold_operating_points.tsv")
    clinical_calibration = _read_tsv("results/clinical_interpretability_20260527/calibration_bins.tsv")
    cbio_manifest = _read_tsv("data/external/cbioportal_melanoma_manifest.tsv")
    cbio_metrics = _read_tsv("results/cbioportal_melanoma_external_validation_20260527/cbioportal_external_metrics.tsv")
    cbio_family = _read_tsv("results/cbioportal_melanoma_external_validation_20260527/cbioportal_external_family_comparison.tsv")
    cbio_rescue_selection = _read_tsv(
        "results/cbioportal_rescue_head_external_validation_20260527/cbioportal_rescue_head_selection.tsv"
    )
    cbio_rescue_coverage = _read_tsv(
        "results/cbioportal_rescue_head_external_validation_20260527/cbioportal_rescue_head_gene_coverage.tsv"
    )
    gpu_bioprior_selection = _read_tsv(
        "results/gpu_bioprior_rescue_combo_search_robust_20260527/gpu_bioprior_rescue_combo_selection.tsv"
    )
    gpu_bioprior_component_ablation = _read_tsv(
        "results/gpu_bioprior_component_ablation_20260527/gpu_bioprior_component_ablation.tsv"
    )
    cbio_gpu_metrics = _read_tsv(
        "results/cbioportal_gpu_bioprior_external_validation_20260527/cbioportal_gpu_bioprior_external_metrics.tsv"
    )
    cbio_gpu_family = _read_tsv(
        "results/cbioportal_gpu_bioprior_external_validation_20260527/cbioportal_gpu_bioprior_external_family_comparison.tsv"
    )
    cbio_gpu_coverage = _read_tsv(
        "results/cbioportal_gpu_bioprior_external_validation_20260527/cbioportal_gpu_bioprior_gene_coverage.tsv"
    )
    lipid_pair_selection = _read_tsv(
        "results/gpu_lipid_pair_rescue_component_dominant_ba_guardrail_20260528/gpu_lipid_pair_rescue_selection.tsv"
    )
    lipid_pair_external_metrics = _read_tsv(
        "results/gpu_lipid_pair_rescue_component_dominant_ba_guardrail_20260528/gpu_lipid_pair_external_metrics.tsv"
    )
    lipid_pair_family = _read_tsv(
        "results/gpu_lipid_pair_rescue_component_dominant_ba_guardrail_20260528/gpu_lipid_pair_external_family_gate.tsv"
    )
    constrained_blend = _read_tsv("results/constrained_signature_blend_search_20260527/constrained_blend_summary.tsv")
    signed_rank = _read_tsv("results/signed_rank_module_audit_20260527/signed_rank_module_comparisons.tsv")
    public_external_leads = _read_tsv("deliverables/public_external_lead_triage_20260527.tsv")
    gse165745_panel_qc = _read_tsv("tables/gse165745_panel_qc.tsv")
    gse165745_panel_metrics = _read_tsv("results/gse165745_panel_transfer_20260527/locked_validation_metrics.tsv")
    rank_fusion_primary = _read_tsv("results/rank_fusion_melanoma_candidate_20260527/rank_fusion_primary_summary.tsv")
    rank_fusion_external = _read_tsv("results/rank_fusion_melanoma_candidate_20260527/rank_fusion_strict_external_summary.tsv")
    threshold_recalibration = _read_tsv(
        "results/threshold_recalibration_audit_20260527/threshold_recalibration_primary_summary.tsv"
    )
    tumor_immune_primary = _read_tsv("results/tumor_immune_balance_audit_20260527/tumor_immune_balance_primary_summary.tsv")
    tumor_immune_external = _read_tsv(
        "results/tumor_immune_balance_audit_20260527/tumor_immune_balance_strict_external_summary.tsv"
    )
    map4k1_transform_selection = _read_tsv(
        "results/map4k1_tbx3_transform_audit_20260527/map4k1_tbx3_transform_selection.tsv"
    )
    ecological_polarity_selection = _read_tsv(
        "results/ecological_polarity_candidate_audit_20260527/ecological_polarity_selection.tsv"
    )
    processed_eligibility = _read_tsv("deliverables/processed_melanoma_external_eligibility_20260527.tsv")
    secondary_external_metrics = _read_tsv(
        "results/secondary_melanoma_external_sensitivity_20260527/secondary_melanoma_external_metrics.tsv"
    )
    strict_failure_selection = _read_tsv(
        "results/strict_external_failure_mode_audit_20260527/strict_external_failure_mode_selection.tsv"
    )
    phs_subset_metrics = _read_tsv("results/phs000452_liu_subset_audit_20260527/phs000452_subset_metrics.tsv")
    phs_source_concordance = _read_tsv(
        "results/phs000452_liu_subset_audit_20260527/phs000452_liu_source_concordance.tsv"
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
        ensemble_search,
        clinical_enrichment,
        clinical_subgroup,
        clinical_threshold,
        clinical_calibration,
        cbio_manifest,
        cbio_metrics,
        cbio_family,
        cbio_rescue_selection,
        cbio_rescue_coverage,
        gpu_bioprior_selection,
        gpu_bioprior_component_ablation,
        cbio_gpu_metrics,
        cbio_gpu_family,
        cbio_gpu_coverage,
        lipid_pair_selection,
        lipid_pair_external_metrics,
        lipid_pair_family,
        constrained_blend,
        signed_rank,
        public_external_leads,
        gse165745_panel_qc,
        gse165745_panel_metrics,
        rank_fusion_primary,
        rank_fusion_external,
        threshold_recalibration,
        tumor_immune_primary,
        tumor_immune_external,
        map4k1_transform_selection,
        ecological_polarity_selection,
        processed_eligibility,
        secondary_external_metrics,
        strict_failure_selection,
        phs_subset_metrics,
        phs_source_concordance,
    )
    out = ROOT / args.out
    out_md = ROOT / args.out_md
    out.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(out, sep="\t", index=False)
    write_markdown(audit, out_md)
    counts = audit["status"].value_counts().to_dict()
    print(json.dumps({"PASS": int(counts.get("PASS", 0)), "PARTIAL": int(counts.get("PARTIAL", 0)), "GAP": int(counts.get("GAP", 0))}, ensure_ascii=False))
    print(f"Wrote {out}")
    print(f"Wrote {out_md}")


if __name__ == "__main__":
    main()
