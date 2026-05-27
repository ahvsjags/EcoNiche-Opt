from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.io import load_processed_bulk
from econiche.metrics import compute_binary_metrics
from econiche.statistics import benjamini_hochberg, paired_bootstrap_delta
from econiche_opt.model.endpoint_modules import (
    MODULE_PRIOR_WEIGHTS,
    build_module_features_by_cohort,
    default_strata,
    endpoint_label_series,
    prepare_endpoint_data,
    select_threshold,
    sigmoid,
)


TARGET_MODEL = "EcoNiche-Opt-AlignedPanelNoCalibration"
CALIBRATED_MODEL = "EcoNiche-Opt-AlignedPanelCalibrated"
ENDPOINT = "primary_recist"
STRATA = ["melanoma_core_high_evidence", "melanoma_recist_supported_primary"]


def _active_real_cohorts(X_by_cohort: dict[str, pd.DataFrame], metadata_by_cohort: dict[str, pd.DataFrame]) -> list[str]:
    cohorts = []
    for cohort in sorted(X_by_cohort):
        if cohort.startswith("demo_cohort_"):
            continue
        meta = metadata_by_cohort.get(cohort)
        if meta is not None and "response_raw" in meta.columns and meta["response_raw"].notna().any():
            cohorts.append(cohort)
    return cohorts


def _weights_for_variant(variant: str) -> dict[str, float]:
    weights = dict(MODULE_PRIOR_WEIGHTS)
    response_modules = {"ifn_t_cell_inflamed", "cytotoxic_cd8", "exhaustion_checkpoint", "antigen_presentation", "trm_tls"}
    resistance_modules = {"myeloid_suppression", "stromal_exclusion"}
    if variant == "full":
        return weights
    if variant == "no_resistance_modules":
        return {module: (0.0 if module in resistance_modules else weight) for module, weight in weights.items()}
    if variant == "no_response_modules":
        return {module: (0.0 if module in response_modules else weight) for module, weight in weights.items()}
    if variant == "unsigned_state_direction":
        return {module: abs(weight) for module, weight in weights.items()}
    if variant == "ifn_only":
        return {module: (1.0 if module == "ifn_t_cell_inflamed" else 0.0) for module in weights}
    if variant == "response_modules_only_equal":
        return {module: (1.0 / len(response_modules) if module in response_modules else 0.0) for module in weights}
    raise ValueError(f"Unknown variant: {variant}")


def _score_modules(features: pd.DataFrame, weights: dict[str, float]) -> pd.Series:
    score = pd.Series(0.0, index=features.index)
    for module, weight in weights.items():
        if module in features.columns:
            score = score + float(weight) * pd.to_numeric(features[module], errors="coerce").fillna(0.0)
    return score.fillna(0.0)


def _fit_monotone_platt(score: pd.Series, y: pd.Series):
    from sklearn.linear_model import LogisticRegression

    common = score.index.intersection(y.index)
    if len(common) < 8 or y.loc[common].nunique() < 2:
        return None
    model = LogisticRegression(C=1.0, class_weight="balanced", solver="lbfgs", max_iter=5000, random_state=20260527)
    model.fit(score.loc[common].astype(float).to_numpy().reshape(-1, 1), y.loc[common].astype(int).to_numpy())
    return model if float(model.coef_[0, 0]) > 0 else None


def _predict_prob(score: pd.Series, calibrator) -> pd.Series:
    if calibrator is None:
        return pd.Series(sigmoid(score), index=score.index)
    prob = calibrator.predict_proba(score.astype(float).to_numpy().reshape(-1, 1))[:, 1]
    return pd.Series(prob, index=score.index)


def _evaluate_stratum(
    endpoint: str,
    stratum: str,
    train_pool: list[str],
    holdouts: list[str],
    module_features_by_cohort: dict[str, pd.DataFrame],
    y_by_cohort: dict[str, pd.Series],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    prediction_rows: list[dict[str, object]] = []
    metric_rows: list[dict[str, object]] = []
    variants = [
        ("full", TARGET_MODEL, False),
        ("full", CALIBRATED_MODEL, True),
        ("no_resistance_modules", "EcoNiche-Opt-NoResistanceModules", True),
        ("no_response_modules", "EcoNiche-Opt-NoResponseModules", True),
        ("unsigned_state_direction", "EcoNiche-Opt-UnsignedStateDirection", True),
        ("ifn_only", "IFNG_ModuleOnly", True),
        ("response_modules_only_equal", "ResponseModulesEqualWeight", True),
    ]
    for holdout in holdouts:
        if holdout not in y_by_cohort or holdout not in module_features_by_cohort:
            continue
        train_cohorts = [cohort for cohort in train_pool if cohort != holdout and cohort in y_by_cohort]
        if not train_cohorts:
            continue
        y_train = pd.concat([y_by_cohort[cohort] for cohort in train_cohorts]).astype(int)
        y_test = y_by_cohort[holdout].astype(int)
        if y_train.nunique() < 2 or y_test.nunique() < 2:
            continue
        for variant, model_name, use_calibration in variants:
            weights = _weights_for_variant(variant)
            train_score = pd.concat(
                [_score_modules(module_features_by_cohort[cohort].reindex(y_by_cohort[cohort].index), weights) for cohort in train_cohorts]
            ).reindex(y_train.index)
            test_score = _score_modules(module_features_by_cohort[holdout].reindex(y_test.index), weights)
            calibrator = _fit_monotone_platt(train_score, y_train) if use_calibration else None
            train_prob = _predict_prob(train_score, calibrator)
            test_prob = _predict_prob(test_score, calibrator)
            threshold = select_threshold(y_train.to_numpy(dtype=int), train_prob.to_numpy(dtype=float))
            metrics = compute_binary_metrics(y_test, test_prob, threshold=threshold)
            metric_rows.append(
                {
                    **metrics,
                    "endpoint": endpoint,
                    "stratum": stratum,
                    "cohort": holdout,
                    "model_name": model_name,
                    "variant": variant,
                    "calibration": "training_only_platt" if use_calibration and calibrator is not None else "raw_sigmoid",
                    "threshold": float(threshold),
                    "train_cohorts": ",".join(train_cohorts),
                    "n_samples": int(len(y_test)),
                    "n_responders": int(y_test.sum()),
                    "n_nonresponders": int((y_test == 0).sum()),
                }
            )
            for sample_id in y_test.index:
                prediction_rows.append(
                    {
                        "endpoint": endpoint,
                        "stratum": stratum,
                        "cohort": holdout,
                        "model_name": model_name,
                        "variant": variant,
                        "sample_id": sample_id,
                        "true_response_label": int(y_test.loc[sample_id]),
                        "response_probability": float(test_prob.loc[sample_id]),
                        "threshold": float(threshold),
                        "pred_response_label": int(test_prob.loc[sample_id] >= threshold),
                    }
                )
    return pd.DataFrame(metric_rows), pd.DataFrame(prediction_rows)


def _summarize_predictions(predictions: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (endpoint, stratum, model_name, variant), frame in predictions.groupby(["endpoint", "stratum", "model_name", "variant"]):
        y = frame["true_response_label"].astype(int)
        p = frame["response_probability"].astype(float)
        if len(frame) < 5 or y.nunique() < 2:
            continue
        metrics = compute_binary_metrics(y, p)
        rows.append(
            {
                "endpoint": endpoint,
                "stratum": stratum,
                "model_name": model_name,
                "variant": variant,
                "n_samples": int(len(frame)),
                "n_cohorts": int(frame["cohort"].nunique()),
                "pooled_AUROC": metrics["AUROC"],
                "pooled_AUPRC": metrics["AUPRC"],
                "pooled_balanced_accuracy": metrics["balanced_accuracy"],
                "pooled_Brier": metrics["Brier"],
                "pooled_ECE": metrics["ECE"],
            }
        )
    return pd.DataFrame(rows).sort_values(["endpoint", "stratum", "pooled_AUROC"], ascending=[True, True, False])


def _paired_comparisons(predictions: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (endpoint, stratum), frame in predictions.groupby(["endpoint", "stratum"]):
        target = frame[frame["model_name"] == TARGET_MODEL]
        if target.empty:
            continue
        for model_name, baseline in frame.groupby("model_name"):
            if model_name == TARGET_MODEL:
                continue
            merged = target.merge(
                baseline[["cohort", "sample_id", "true_response_label", "response_probability"]],
                on=["cohort", "sample_id", "true_response_label"],
                suffixes=("_target", "_baseline"),
            )
            if len(merged) < 8 or merged["true_response_label"].nunique() < 2:
                continue
            y = merged["true_response_label"].astype(int)
            target_prob = merged["response_probability_target"].astype(float)
            baseline_prob = merged["response_probability_baseline"].astype(float)
            target_metrics = compute_binary_metrics(y, target_prob)
            baseline_metrics = compute_binary_metrics(y, baseline_prob)
            stats = paired_bootstrap_delta(y, target_prob, baseline_prob, n_bootstrap=1000, random_state=20260527)
            rows.append(
                {
                    "endpoint": endpoint,
                    "stratum": stratum,
                    "target_model": TARGET_MODEL,
                    "ablation_model": model_name,
                    "n_samples": int(len(merged)),
                    "target_AUROC": target_metrics["AUROC"],
                    "ablation_AUROC": baseline_metrics["AUROC"],
                    "delta_AUROC": float(target_metrics["AUROC"] - baseline_metrics["AUROC"]),
                    "target_ECE": target_metrics["ECE"],
                    "ablation_ECE": baseline_metrics["ECE"],
                    "delta_ECE": float(target_metrics["ECE"] - baseline_metrics["ECE"]),
                    **stats,
                }
            )
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result["fdr_q"] = 1.0
    for _, idx in result.groupby(["endpoint", "stratum"]).groups.items():
        result.loc[idx, "fdr_q"] = benjamini_hochberg(result.loc[idx, "p_value"].fillna(1.0))
    result["claim_level"] = np.where(
        result["ablation_model"] == CALIBRATED_MODEL,
        np.where(result["delta_ECE"] > 0, "calibration_improves_ECE_with_discrimination_tradeoff", "calibration_not_supported"),
        np.where(
            (result["delta_AUROC"] > 0) & (result["fdr_q"] <= 0.05),
            "FDR_supported_component_gain",
            np.where(result["delta_AUROC"] > 0, "point_estimate_component_gain", "component_not_performance_supported"),
        ),
    )
    return result.sort_values(["endpoint", "stratum", "delta_AUROC"], ascending=[True, True, False])


def _write_audit(out_dir: Path, summary: pd.DataFrame, pairwise: pd.DataFrame) -> None:
    lines = [
        "# Aligned Locked-Panel Ablation Audit",
        "",
        "This audit tests ablation variants on the same locked module-panel scoring family that supports the primary EcoNiche-Opt melanoma result. It replaces the older WordFullGraph-only ablation as the primary component-evidence table because WordFullGraph is not the strongest discriminative model.",
        "",
        "## Pooled Model Summary",
        "",
    ]
    if summary.empty:
        lines.append("No aligned panel ablation summary rows were produced.")
    else:
        for _, row in summary.iterrows():
            lines.append(
                f"- {row['endpoint']} / {row['stratum']} / {row['model_name']}: n={int(row['n_samples'])}, "
                f"AUROC={float(row['pooled_AUROC']):.3f}, AUPRC={float(row['pooled_AUPRC']):.3f}, "
                f"BA={float(row['pooled_balanced_accuracy']):.3f}, ECE={float(row['pooled_ECE']):.3f}."
            )
    lines.extend(["", "## Component Evidence Boundary", ""])
    if pairwise.empty:
        lines.append("No paired ablation comparisons were produced.")
    else:
        for _, row in pairwise.iterrows():
            lines.append(
                f"- {row['endpoint']} / {row['stratum']} vs {row['ablation_model']}: "
                f"delta AUROC={float(row['delta_AUROC']):.3f}, delta ECE={float(row['delta_ECE']):.3f}, "
                f"95% CI [{float(row['ci_low']):.3f}, {float(row['ci_high']):.3f}], "
                f"FDR q={float(row['fdr_q']):.3f} ({row['claim_level']})."
            )
    lines.extend(
        [
            "",
            "## Claim Rule",
            "",
            "Use performance-gain language only for variants with positive delta AUROC and FDR support. Components without performance support can still be described as biological representation, interpretation, calibration, or robustness components if their corresponding metric supports that narrower claim.",
            "",
        ]
    )
    (out_dir / "ALIGNED_PANEL_ABLATION_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")


def run(processed_dir: Path, out_dir: Path) -> None:
    X_by_cohort, _, metadata_by_cohort = load_processed_bulk(processed_dir)
    active = _active_real_cohorts(X_by_cohort, metadata_by_cohort)
    endpoint_data = prepare_endpoint_data(X_by_cohort, metadata_by_cohort, active, ENDPOINT)
    module_features_by_cohort, coverage = build_module_features_by_cohort(endpoint_data.X_by_cohort)
    strata = default_strata(endpoint_data.X_by_cohort)
    metric_frames = []
    prediction_frames = []
    for stratum in STRATA:
        spec = strata.get(stratum, {})
        train_pool = [cohort for cohort in spec.get("train_pool", []) if cohort in endpoint_data.y_response_by_cohort]
        holdouts = [cohort for cohort in spec.get("holdouts", []) if cohort in endpoint_data.y_response_by_cohort]
        metrics, predictions = _evaluate_stratum(
            ENDPOINT,
            stratum,
            train_pool,
            holdouts,
            module_features_by_cohort,
            endpoint_data.y_response_by_cohort,
        )
        metric_frames.append(metrics)
        prediction_frames.append(predictions)
    metrics = pd.concat([frame for frame in metric_frames if not frame.empty], ignore_index=True) if metric_frames else pd.DataFrame()
    predictions = (
        pd.concat([frame for frame in prediction_frames if not frame.empty], ignore_index=True) if prediction_frames else pd.DataFrame()
    )
    summary = _summarize_predictions(predictions) if not predictions.empty else pd.DataFrame()
    pairwise = _paired_comparisons(predictions) if not predictions.empty else pd.DataFrame()
    out_dir.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(out_dir / "aligned_panel_ablation_lodo_metrics.tsv", sep="\t", index=False)
    predictions.to_csv(out_dir / "aligned_panel_ablation_predictions.tsv", sep="\t", index=False)
    summary.to_csv(out_dir / "aligned_panel_ablation_summary.tsv", sep="\t", index=False)
    pairwise.to_csv(out_dir / "aligned_panel_ablation_pairwise.tsv", sep="\t", index=False)
    coverage.to_csv(out_dir / "aligned_panel_ablation_gene_coverage.tsv", sep="\t", index=False)
    _write_audit(out_dir, summary, pairwise)
    print(f"Wrote aligned panel ablation outputs to {out_dir}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--out", default="results/aligned_panel_ablation_20260527")
    args = parser.parse_args()
    run(ROOT / args.processed_dir, ROOT / args.out)


if __name__ == "__main__":
    main()
