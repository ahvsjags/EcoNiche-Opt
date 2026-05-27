from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn import metrics as sk_metrics
from sklearn.linear_model import LogisticRegression

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.io import load_processed_bulk
from econiche.metrics import compute_binary_metrics
from econiche.statistics import benjamini_hochberg, paired_bootstrap_delta
from econiche_opt.model.endpoint_modules import (
    MODULE_GENE_SETS,
    OPTIMIZED_ADAPTIVE_MODEL,
    build_fixed_scores_by_cohort,
    build_module_features_by_cohort,
    endpoint_label_series,
    select_threshold,
    sigmoid,
)


DISCOVERY_COHORTS = ["GSE91061", "GSE78220", "PRJEB23709_PD1_PRE"]
LOCKED_EXTERNAL_COHORTS = [
    "GSE145996",
    "PHS000452_LIU_LIKE_PRE",
    "PRJEB23709_COMBO_PRE",
    "GSE93157",
    "GSE140901",
]
EIGHT_SIGNATURES = ["IFNG", "CXCL9", "TIG", "TIDE_dysfunction", "APM", "CYT", "IPRES", "TIDE_exclusion"]
TARGET_PANEL_MODEL = f"{OPTIMIZED_ADAPTIVE_MODEL}-LockedPanel"
NANOSTRING_COHORTS = {"GSE93157", "GSE140901"}


def _analysis_type(cohort: str) -> str:
    if cohort == "GSE145996":
        return "locked_external_melanoma_pd1_recist"
    if cohort == "PHS000452_LIU_LIKE_PRE":
        return "locked_external_melanoma_pd1_like"
    if cohort == "PRJEB23709_COMBO_PRE":
        return "locked_external_melanoma_combination_transfer"
    if cohort == "GSE93157":
        return "nanostring_melanoma_clinical_assay_transfer"
    if cohort == "GSE140901":
        return "nanostring_hcc_clinical_assay_transfer"
    return "locked_external_other"


def _cohort_meta_value(meta: pd.DataFrame, column: str) -> str:
    if column not in meta.columns:
        return ""
    values = meta[column].dropna().astype(str).unique()
    return ";".join(sorted(values[:6]))


def _endpoint_labels(metadata_by_cohort: dict[str, pd.DataFrame], cohorts: list[str], endpoint: str) -> dict[str, pd.Series]:
    labels: dict[str, pd.Series] = {}
    for cohort in cohorts:
        if cohort not in metadata_by_cohort:
            continue
        meta = metadata_by_cohort[cohort]
        if "response_raw" not in meta.columns:
            continue
        y = endpoint_label_series(meta["response_raw"], endpoint)
        y = y.dropna().astype(int)
        if len(y) >= 4 and y.nunique() == 2:
            labels[cohort] = y
    return labels


def _available_scores(
    fixed_scores: dict[str, dict[str, pd.Series]],
    module_features_by_cohort: dict[str, pd.DataFrame],
    model_names: list[str],
) -> dict[str, dict[str, pd.Series]]:
    out: dict[str, dict[str, pd.Series]] = {}
    for cohort, module_features in module_features_by_cohort.items():
        cohort_scores: dict[str, pd.Series] = {
            TARGET_PANEL_MODEL: pd.Series(module_features.mul(
                pd.Series(
                    {
                        "ifn_t_cell_inflamed": 1.0,
                        "cytotoxic_cd8": 0.5,
                        "exhaustion_checkpoint": 0.25,
                        "antigen_presentation": 0.5,
                        "myeloid_suppression": -0.5,
                        "stromal_exclusion": -0.5,
                        "trm_tls": 0.25,
                    }
                ),
                axis=1,
            ).sum(axis=1), index=module_features.index).fillna(0.0)
        }
        for model in model_names:
            if cohort in fixed_scores and model in fixed_scores[cohort]:
                cohort_scores[model] = fixed_scores[cohort][model].fillna(0.0)
        out[cohort] = cohort_scores
    return out


def _fit_thresholds(
    labels_by_cohort: dict[str, pd.Series],
    scores_by_cohort: dict[str, dict[str, pd.Series]],
    discovery_cohorts: list[str],
    model_names: list[str],
) -> dict[str, dict[str, object]]:
    thresholds: dict[str, dict[str, object]] = {}
    for model in model_names:
        train_scores = []
        y_parts = []
        p_parts = []
        train_cohorts = []
        for cohort in discovery_cohorts:
            if cohort not in labels_by_cohort or cohort not in scores_by_cohort or model not in scores_by_cohort[cohort]:
                continue
            y = labels_by_cohort[cohort]
            score = scores_by_cohort[cohort][model].reindex(y.index).dropna()
            common = y.index.intersection(score.index)
            if len(common) < 4:
                continue
            y_parts.append(y.loc[common])
            train_scores.append(score.loc[common])
            p_parts.append(pd.Series(sigmoid(score.loc[common]), index=common))
            train_cohorts.append(cohort)
        if not y_parts:
            continue
        y_train = pd.concat(y_parts)
        raw_train_score = pd.concat(train_scores).reindex(y_train.index)
        calibrator = _fit_monotone_calibrator(raw_train_score, y_train)
        p_train = _predict_with_calibrator(raw_train_score, calibrator)
        if y_train.nunique() < 2:
            continue
        threshold = select_threshold(y_train.to_numpy(dtype=int), p_train.to_numpy(dtype=float))
        train_metrics = compute_binary_metrics(y_train, p_train, threshold=threshold)
        calibration_coef = float(calibrator.coef_[0, 0]) if calibrator is not None else 1.0
        calibration_intercept = float(calibrator.intercept_[0]) if calibrator is not None else 0.0
        thresholds[model] = {
            "threshold": float(threshold),
            "calibrator": calibrator,
            "calibration": "discovery_only_platt",
            "calibration_coef": calibration_coef,
            "calibration_intercept": calibration_intercept,
            "training_n": int(len(y_train)),
            "training_responders": int(y_train.sum()),
            "training_nonresponders": int((y_train == 0).sum()),
            "training_AUROC": float(train_metrics.get("AUROC", np.nan)),
            "training_balanced_accuracy": float(train_metrics.get("balanced_accuracy", np.nan)),
            "training_cohorts": ",".join(train_cohorts),
        }
    return thresholds


def _fit_monotone_calibrator(score: pd.Series, y: pd.Series) -> LogisticRegression | None:
    """Fit a one-dimensional Platt calibrator using discovery labels only."""
    common = score.index.intersection(y.index)
    if len(common) < 8 or y.loc[common].nunique() < 2:
        return None
    model = LogisticRegression(C=1.0, class_weight="balanced", solver="lbfgs", max_iter=5000, random_state=20260519)
    model.fit(score.loc[common].astype(float).to_numpy().reshape(-1, 1), y.loc[common].astype(int).to_numpy())
    if float(model.coef_[0, 0]) <= 0:
        return None
    return model


def _predict_training_only_calibrated(score: pd.Series, y_train: pd.Series) -> pd.Series:
    calibrator = _fit_monotone_calibrator(score, y_train)
    return _predict_with_calibrator(score, calibrator)


def _predict_with_calibrator(score: pd.Series, calibrator: LogisticRegression | None) -> pd.Series:
    if calibrator is None:
        return pd.Series(sigmoid(score), index=score.index)
    prob = calibrator.predict_proba(score.astype(float).to_numpy().reshape(-1, 1))[:, 1]
    return pd.Series(prob, index=score.index)


def _metric_rows(
    labels_by_cohort: dict[str, pd.Series],
    metadata_by_cohort: dict[str, pd.DataFrame],
    scores_by_cohort: dict[str, dict[str, pd.Series]],
    thresholds: dict[str, dict[str, object]],
    cohorts: list[str],
    model_names: list[str],
    endpoint: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    prediction_rows: list[dict[str, object]] = []
    metric_rows: list[dict[str, object]] = []
    for cohort in cohorts:
        if cohort not in labels_by_cohort:
            continue
        y = labels_by_cohort[cohort]
        meta = metadata_by_cohort[cohort].reindex(y.index)
        for model in model_names:
            if model not in thresholds or cohort not in scores_by_cohort or model not in scores_by_cohort[cohort]:
                continue
            score = scores_by_cohort[cohort][model].reindex(y.index).dropna()
            common = y.index.intersection(score.index)
            if len(common) < 4 or y.loc[common].nunique() < 2:
                continue
            prob = _predict_with_calibrator(score.loc[common], thresholds[model].get("calibrator"))
            threshold = float(thresholds[model]["threshold"])
            metrics = compute_binary_metrics(y.loc[common], prob, threshold=threshold)
            row = {
                "endpoint": endpoint,
                "cohort": cohort,
                "analysis_type": _analysis_type(cohort),
                "model_name": model,
                "n_samples": int(len(common)),
                "n_responders": int(y.loc[common].sum()),
                "n_nonresponders": int((y.loc[common] == 0).sum()),
                "threshold_source": "discovery_only",
                "threshold_training_cohorts": thresholds[model]["training_cohorts"],
                "threshold": threshold,
                "calibration": thresholds[model].get("calibration", "raw_sigmoid"),
                "training_AUROC": thresholds[model]["training_AUROC"],
                "training_balanced_accuracy": thresholds[model]["training_balanced_accuracy"],
                "platform": _cohort_meta_value(meta.loc[common], "platform"),
                "treatment": _cohort_meta_value(meta.loc[common], "treatment"),
                **metrics,
            }
            metric_rows.append(row)
            for sample_id in common:
                prediction_rows.append(
                    {
                        "endpoint": endpoint,
                        "cohort": cohort,
                        "analysis_type": _analysis_type(cohort),
                        "model_name": model,
                        "sample_id": sample_id,
                        "true_response_label": int(y.loc[sample_id]),
                        "response_raw": meta.loc[sample_id, "response_raw"] if "response_raw" in meta.columns else "",
                        "response_probability": float(prob.loc[sample_id]),
                        "locked_threshold": threshold,
                        "predicted_response_label": int(prob.loc[sample_id] >= threshold),
                        "threshold_source": "discovery_only",
                        "threshold_training_cohorts": thresholds[model]["training_cohorts"],
                        "platform": meta.loc[sample_id, "platform"] if "platform" in meta.columns else "",
                        "treatment": meta.loc[sample_id, "treatment"] if "treatment" in meta.columns else "",
                    }
                )
    return pd.DataFrame(metric_rows), pd.DataFrame(prediction_rows)


def _pairwise_rows(predictions: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (endpoint, analysis_type, cohort), frame in predictions.groupby(["endpoint", "analysis_type", "cohort"]):
        target = frame[frame["model_name"] == TARGET_PANEL_MODEL]
        if target.empty:
            continue
        for baseline in EIGHT_SIGNATURES:
            base = frame[frame["model_name"] == baseline]
            if base.empty:
                continue
            merged = target.merge(
                base[["sample_id", "true_response_label", "response_probability"]],
                on=["sample_id", "true_response_label"],
                suffixes=("_target", "_baseline"),
            )
            if len(merged) < 8 or merged["true_response_label"].nunique() < 2:
                continue
            y = merged["true_response_label"].astype(int)
            target_prob = merged["response_probability_target"].astype(float)
            baseline_prob = merged["response_probability_baseline"].astype(float)
            delta = sk_metrics.roc_auc_score(y, target_prob) - sk_metrics.roc_auc_score(y, baseline_prob)
            stats = paired_bootstrap_delta(y, target_prob, baseline_prob, n_bootstrap=1000, random_state=20260507)
            rows.append(
                {
                    "endpoint": endpoint,
                    "analysis_type": analysis_type,
                    "cohort": cohort,
                    "target_model": TARGET_PANEL_MODEL,
                    "baseline_model": baseline,
                    "n_samples": int(len(merged)),
                    "target_AUROC": float(sk_metrics.roc_auc_score(y, target_prob)),
                    "baseline_AUROC": float(sk_metrics.roc_auc_score(y, baseline_prob)),
                    "delta_AUROC": float(delta),
                    **stats,
                }
            )
    pairwise = pd.DataFrame(rows)
    if pairwise.empty:
        return pairwise
    pairwise["fdr_q"] = 1.0
    for _, idx in pairwise.groupby(["endpoint", "analysis_type", "cohort"]).groups.items():
        pairwise.loc[idx, "fdr_q"] = benjamini_hochberg(pairwise.loc[idx, "p_value"].fillna(1.0))
    pairwise["claim_level"] = np.where(
        (pairwise["delta_AUROC"] > 0) & (pairwise["fdr_q"] <= 0.05),
        "FDR_supported",
        np.where(pairwise["delta_AUROC"] > 0, "point_estimate_only", "not_superior"),
    )
    return pairwise


def _family_bootstrap(
    y_true: pd.Series,
    target_prob: pd.Series,
    baseline_probs: pd.DataFrame,
    n_bootstrap: int = 1000,
) -> dict[str, float]:
    y = y_true.to_numpy(dtype=int)
    target = target_prob.to_numpy(dtype=float)
    baselines = baseline_probs.to_numpy(dtype=float)
    rng = np.random.default_rng(20260507)
    deltas: list[float] = []
    for _ in range(n_bootstrap):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        target_auc = sk_metrics.roc_auc_score(y[idx], target[idx])
        baseline_aucs = [sk_metrics.roc_auc_score(y[idx], baselines[idx, col]) for col in range(baselines.shape[1])]
        deltas.append(float(target_auc - np.mean(baseline_aucs)))
    if not deltas:
        return {}
    arr = np.asarray(deltas)
    observed_target_auc = float(sk_metrics.roc_auc_score(y, target))
    observed_baselines = [sk_metrics.roc_auc_score(y, baseline_probs[col]) for col in baseline_probs.columns]
    one_sided = float((arr <= 0).mean())
    return {
        "target_AUROC": observed_target_auc,
        "mean_signature_AUROC": float(np.mean(observed_baselines)),
        "best_signature_AUROC": float(np.max(observed_baselines)),
        "mean_delta_vs_signature_family": float(observed_target_auc - np.mean(observed_baselines)),
        "ci_low": float(np.quantile(arr, 0.025)),
        "ci_high": float(np.quantile(arr, 0.975)),
        "one_sided_p": one_sided,
        "two_sided_p": float(min(1.0, 2 * min(one_sided, float((arr >= 0).mean())))),
        "n_signatures": int(baseline_probs.shape[1]),
    }


def _family_omnibus_rows(predictions: pd.DataFrame) -> pd.DataFrame:
    if predictions.empty:
        return pd.DataFrame()
    family_masks = {
        "all_locked_external_and_panel": lambda frame: pd.Series(True, index=frame.index),
        "melanoma_external_and_panel": lambda frame: frame["cohort"].isin(
            ["GSE145996", "PHS000452_LIU_LIKE_PRE", "PRJEB23709_COMBO_PRE", "GSE93157"]
        ),
        "strict_pd1_like_external": lambda frame: frame["cohort"].isin(["GSE145996", "PHS000452_LIU_LIKE_PRE"]),
        "nanostring_panel_transfer": lambda frame: frame["cohort"].isin(["GSE93157", "GSE140901"]),
    }
    rows: list[dict[str, object]] = []
    for endpoint, endpoint_frame in predictions.groupby("endpoint"):
        for family_name, mask_fn in family_masks.items():
            frame = endpoint_frame[mask_fn(endpoint_frame)]
            target = frame[frame["model_name"] == TARGET_PANEL_MODEL]
            if target.empty:
                continue
            y_ref = None
            target_ref = None
            baseline_series: dict[str, pd.Series] = {}
            for baseline in EIGHT_SIGNATURES:
                base = frame[frame["model_name"] == baseline]
                merged = target.merge(
                    base[["cohort", "sample_id", "true_response_label", "response_probability"]],
                    on=["cohort", "sample_id", "true_response_label"],
                    suffixes=("_target", "_baseline"),
                )
                if len(merged) < 8 or merged["true_response_label"].nunique() < 2:
                    continue
                key = merged["cohort"].astype(str) + "::" + merged["sample_id"].astype(str)
                baseline_series[baseline] = pd.Series(merged["response_probability_baseline"].to_numpy(dtype=float), index=key)
                if y_ref is None:
                    y_ref = pd.Series(merged["true_response_label"].to_numpy(dtype=int), index=key)
                    target_ref = pd.Series(merged["response_probability_target"].to_numpy(dtype=float), index=key)
            if len(baseline_series) < 4 or y_ref is None or target_ref is None:
                continue
            baseline_frame = pd.DataFrame(baseline_series).dropna(axis=0)
            common = baseline_frame.index.intersection(y_ref.index).intersection(target_ref.index)
            if len(common) < 12 or y_ref.loc[common].nunique() < 2:
                continue
            stats = _family_bootstrap(y_ref.loc[common], target_ref.loc[common], baseline_frame.loc[common])
            if not stats:
                continue
            rows.append(
                {
                    "endpoint": endpoint,
                    "validation_family": family_name,
                    "target_model": TARGET_PANEL_MODEL,
                    "baseline_family": "eight_strong_signatures",
                    "n_samples": int(len(common)),
                    **stats,
                }
            )
    omnibus = pd.DataFrame(rows)
    if omnibus.empty:
        return omnibus
    omnibus["two_sided_fdr_q"] = 1.0
    omnibus["one_sided_fdr_q"] = 1.0
    for _, idx in omnibus.groupby("validation_family").groups.items():
        omnibus.loc[idx, "two_sided_fdr_q"] = benjamini_hochberg(omnibus.loc[idx, "two_sided_p"].fillna(1.0))
        omnibus.loc[idx, "one_sided_fdr_q"] = benjamini_hochberg(omnibus.loc[idx, "one_sided_p"].fillna(1.0))
    omnibus["claim_level"] = np.where(
        (omnibus["mean_delta_vs_signature_family"] > 0) & (omnibus["two_sided_fdr_q"] <= 0.05),
        "family_two_sided_FDR_supported",
        np.where(
            (omnibus["mean_delta_vs_signature_family"] > 0) & (omnibus["one_sided_fdr_q"] <= 0.05),
            "family_pre_directional_FDR_supported",
            np.where(omnibus["mean_delta_vs_signature_family"] > 0, "family_point_estimate_only", "family_not_superior"),
        ),
    )
    return omnibus


def _panel_coverage_rows(
    X_by_cohort: dict[str, pd.DataFrame],
    coverage: pd.DataFrame,
    cohorts: list[str],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for cohort in cohorts:
        if cohort not in X_by_cohort:
            continue
        columns = set(X_by_cohort[cohort].columns)
        for module, genes in MODULE_GENE_SETS.items():
            available = [gene for gene in genes if gene in columns]
            rows.append(
                {
                    "cohort": cohort,
                    "analysis_type": _analysis_type(cohort),
                    "module": module,
                    "n_module_genes": len(genes),
                    "n_available_genes": len(available),
                    "coverage_fraction": float(len(available) / len(genes)) if genes else np.nan,
                    "available_genes": ",".join(available),
                    "is_nanostring_panel_transfer": cohort in NANOSTRING_COHORTS,
                }
            )
    if not coverage.empty:
        rows.append(
            {
                "cohort": "__all__",
                "analysis_type": "module_feature_builder_summary",
                "module": "__builder_rows__",
                "n_module_genes": int(coverage["n_genes"].sum()) if "n_genes" in coverage.columns else np.nan,
                "n_available_genes": int(coverage["n_available"].sum()) if "n_available" in coverage.columns else np.nan,
                "coverage_fraction": np.nan,
                "available_genes": "",
                "is_nanostring_panel_transfer": False,
            }
        )
    return pd.DataFrame(rows)


def _write_audit(
    out_dir: Path,
    metrics: pd.DataFrame,
    pairwise: pd.DataFrame,
    family_omnibus: pd.DataFrame,
    panel_coverage: pd.DataFrame,
    endpoints: list[str],
    discovery_cohorts: list[str],
    external_cohorts: list[str],
) -> None:
    lines = [
        "# Locked External and Clinical-Assay Panel Validation Audit",
        "",
        "This audit locks threshold selection to the discovery melanoma cohorts and evaluates untouched external cohorts and public NanoString panel-transfer cohorts. It is not a prospective wet-lab validation; true prospective validation remains a future clinical study requirement.",
        "",
        f"- Discovery-only threshold cohorts: {', '.join(discovery_cohorts)}",
        f"- Locked external/panel cohorts: {', '.join(external_cohorts)}",
        f"- Endpoints: {', '.join(endpoints)}",
        "",
        "## Target Model External Metrics",
        "",
    ]
    target_metrics = metrics[metrics["model_name"] == TARGET_PANEL_MODEL] if not metrics.empty else pd.DataFrame()
    if target_metrics.empty:
        lines.append("No target metrics were produced.")
    else:
        for _, row in target_metrics.iterrows():
            lines.append(
                f"- {row['endpoint']} / {row['cohort']} ({row['analysis_type']}): n={int(row['n_samples'])}, "
                f"AUROC={row['AUROC']:.3f}, balanced_accuracy={row['balanced_accuracy']:.3f}, "
                f"ECE={row['ECE']:.3f}, threshold={row['threshold']:.3f}."
            )
    if not pairwise.empty:
        supported = pairwise[pairwise["claim_level"] == "FDR_supported"]
        positive = pairwise[pairwise["claim_level"] == "point_estimate_only"]
        lines.extend(["", "## Baseline Comparison Boundary", ""])
        lines.append(f"- FDR-supported per-cohort comparisons: {len(supported)}")
        lines.append(f"- Positive point-estimate comparisons without FDR support: {len(positive)}")
        if not supported.empty:
            for _, row in supported.head(12).iterrows():
                lines.append(
                    f"- Supported: {row['endpoint']} / {row['cohort']} vs {row['baseline_model']}: "
                    f"delta AUROC={row['delta_AUROC']:.3f}, FDR q={row['fdr_q']:.3f}."
                )
    if not family_omnibus.empty:
        lines.extend(["", "## External Signature-Family Omnibus", ""])
        for _, row in family_omnibus.iterrows():
            lines.append(
                f"- {row['endpoint']} / {row['validation_family']}: target AUROC={row['target_AUROC']:.3f}, "
                f"mean signature AUROC={row['mean_signature_AUROC']:.3f}, best signature AUROC={row['best_signature_AUROC']:.3f}, "
                f"delta vs family mean={row['mean_delta_vs_signature_family']:.3f}, two-sided FDR q={row['two_sided_fdr_q']:.3f} "
                f"({row['claim_level']})."
            )
    panel = panel_coverage[panel_coverage["is_nanostring_panel_transfer"] == True] if not panel_coverage.empty else pd.DataFrame()
    if not panel.empty:
        lines.extend(["", "## NanoString Panel Gene Coverage", ""])
        summary = panel.groupby("cohort")["coverage_fraction"].mean().reset_index()
        for _, row in summary.iterrows():
            lines.append(f"- {row['cohort']}: mean module gene coverage={row['coverage_fraction']:.3f}.")
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "Allowed: 'locked external/panel-transfer validation was run without threshold leakage' and cohort-specific performance statements from the tables. Avoid: 'prospective clinical validation completed' unless a new prospective assay cohort is generated outside this retrospective public-data pipeline.",
            "",
        ]
    )
    (out_dir / "LOCKED_EXTERNAL_PANEL_VALIDATION_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")


def run(
    processed_dir: Path,
    out_dir: Path,
    endpoints: list[str],
    discovery_cohorts: list[str],
    external_cohorts: list[str],
) -> None:
    X_by_cohort, _, metadata_by_cohort = load_processed_bulk(processed_dir)
    selected_cohorts = [cohort for cohort in dict.fromkeys(discovery_cohorts + external_cohorts) if cohort in X_by_cohort]
    module_features_by_cohort, coverage = build_module_features_by_cohort({cohort: X_by_cohort[cohort] for cohort in selected_cohorts})
    fixed_scores = build_fixed_scores_by_cohort(
        {cohort: X_by_cohort[cohort] for cohort in selected_cohorts},
        module_features_by_cohort,
        baselines=EIGHT_SIGNATURES,
    )
    all_scores = _available_scores(fixed_scores, module_features_by_cohort, EIGHT_SIGNATURES)
    model_names = [TARGET_PANEL_MODEL] + EIGHT_SIGNATURES
    all_metrics: list[pd.DataFrame] = []
    all_predictions: list[pd.DataFrame] = []
    threshold_rows: list[dict[str, object]] = []
    for endpoint in endpoints:
        labels_by_cohort = _endpoint_labels(metadata_by_cohort, selected_cohorts, endpoint)
        thresholds = _fit_thresholds(labels_by_cohort, all_scores, discovery_cohorts, model_names)
        for model, threshold_info in thresholds.items():
            threshold_rows.append(
                {
                    "endpoint": endpoint,
                    "model_name": model,
                    **{key: value for key, value in threshold_info.items() if key != "calibrator"},
                }
            )
        metrics, predictions = _metric_rows(
            labels_by_cohort,
            metadata_by_cohort,
            all_scores,
            thresholds,
            external_cohorts,
            model_names,
            endpoint,
        )
        all_metrics.append(metrics)
        all_predictions.append(predictions)
    metrics = pd.concat([frame for frame in all_metrics if not frame.empty], axis=0, ignore_index=True) if all_metrics else pd.DataFrame()
    predictions = (
        pd.concat([frame for frame in all_predictions if not frame.empty], axis=0, ignore_index=True) if all_predictions else pd.DataFrame()
    )
    pairwise = _pairwise_rows(predictions) if not predictions.empty else pd.DataFrame()
    family_omnibus = _family_omnibus_rows(predictions) if not predictions.empty else pd.DataFrame()
    panel_coverage = _panel_coverage_rows(X_by_cohort, coverage, external_cohorts)
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(threshold_rows).to_csv(out_dir / "locked_thresholds.tsv", sep="\t", index=False)
    metrics.to_csv(out_dir / "locked_external_metrics.tsv", sep="\t", index=False)
    predictions.to_csv(out_dir / "locked_external_predictions.tsv", sep="\t", index=False)
    pairwise.to_csv(out_dir / "locked_external_pairwise.tsv", sep="\t", index=False)
    family_omnibus.to_csv(out_dir / "locked_external_signature_family_omnibus.tsv", sep="\t", index=False)
    panel_coverage.to_csv(out_dir / "clinical_assay_panel_transfer.tsv", sep="\t", index=False)
    _write_audit(out_dir, metrics, pairwise, family_omnibus, panel_coverage, endpoints, discovery_cohorts, external_cohorts)
    print(f"Wrote locked external/panel validation outputs to {out_dir}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--out", default="results/locked_external_panel_validation")
    parser.add_argument("--endpoint", action="append", dest="endpoints", default=None)
    parser.add_argument("--discovery-cohort", action="append", dest="discovery_cohorts", default=None)
    parser.add_argument("--external-cohort", action="append", dest="external_cohorts", default=None)
    args = parser.parse_args()
    endpoints = args.endpoints or ["primary_recist", "strict_recist", "clinical_benefit"]
    discovery_cohorts = args.discovery_cohorts or DISCOVERY_COHORTS
    external_cohorts = args.external_cohorts or LOCKED_EXTERNAL_COHORTS
    run(ROOT / args.processed_dir, ROOT / args.out, endpoints, discovery_cohorts, external_cohorts)


if __name__ == "__main__":
    main()
