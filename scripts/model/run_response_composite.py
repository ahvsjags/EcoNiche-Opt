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
from econiche_opt.model.response_composite import (
    NONRESPONSE_HIGH_SIGNATURES,
    RESPONSE_HIGH_SIGNATURES,
    PREFERRED_CANDIDATE,
    run_nested_response_composite,
)


def _write_global_metrics(predictions: pd.DataFrame, metrics: pd.DataFrame, out: Path) -> None:
    y = predictions["true_response_label"].astype(int)
    p = predictions["response_probability"].astype(float)
    global_metrics = compute_binary_metrics(y, p)
    rows = [
        {"scope": "global_pooled", "metric": metric, "estimate": value, "n": len(predictions)}
        for metric, value in global_metrics.items()
    ]
    for metric in ["AUROC", "AUPRC", "balanced_accuracy", "MCC", "F1", "Brier", "ECE"]:
        values = pd.to_numeric(metrics.get(metric), errors="coerce")
        rows.append(
            {
                "scope": "mean_lodo_cohort",
                "metric": metric,
                "estimate": float(values.mean()) if values.notna().any() else np.nan,
                "n": int(values.notna().sum()),
            }
        )
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)


def _baseline_response_probability(frame: pd.DataFrame, model_name: str) -> pd.Series:
    if model_name in RESPONSE_HIGH_SIGNATURES:
        return frame["pred_prob"].astype(float)
    if model_name in NONRESPONSE_HIGH_SIGNATURES:
        return 1.0 - frame["pred_prob"].astype(float)
    raise KeyError(model_name)


def _compare_with_baselines(predictions: pd.DataFrame, baseline_path: Path) -> pd.DataFrame:
    if not baseline_path.exists():
        return pd.DataFrame(
            [
                {
                    "comparison": "EcoNiche-Opt-ImmuneComposite_vs_baseline",
                    "status": "RESULT_PENDING",
                    "reason": "missing baseline predictions",
                }
            ]
        )
    baseline = pd.read_csv(baseline_path, sep="\t")
    comparable = sorted(set(RESPONSE_HIGH_SIGNATURES + NONRESPONSE_HIGH_SIGNATURES) & set(baseline["model_name"].dropna()))
    rows = []
    for model_name in comparable:
        base = baseline[baseline["model_name"] == model_name].dropna(subset=["pred_prob", "true_label"]).copy()
        if base.empty:
            continue
        base["baseline_response_probability"] = _baseline_response_probability(base, model_name)
        merged = predictions.merge(
            base[["sample_id", "cohort", "baseline_response_probability"]],
            on=["sample_id", "cohort"],
            how="inner",
        ).dropna(subset=["response_probability", "baseline_response_probability", "true_response_label"])
        if len(merged) < 5 or merged["true_response_label"].nunique() < 2:
            continue
        stats = paired_bootstrap_delta(
            merged["true_response_label"],
            merged["response_probability"],
            merged["baseline_response_probability"],
            n_bootstrap=1000,
        )
        model_metrics = compute_binary_metrics(merged["true_response_label"], merged["response_probability"])
        baseline_metrics = compute_binary_metrics(merged["true_response_label"], merged["baseline_response_probability"])
        rows.append(
            {
                **stats,
                "comparison": f"EcoNiche-Opt-ImmuneComposite_vs_{model_name}",
                "metric": "response_AUROC",
                "n_samples": len(merged),
                "optimized_AUROC": model_metrics["AUROC"],
                "baseline_AUROC": baseline_metrics["AUROC"],
                "baseline_direction": "response_high" if model_name in RESPONSE_HIGH_SIGNATURES else "nonresponse_high_inverted",
                "status": "computed_from_pipeline",
            }
        )
    result = pd.DataFrame(rows)
    if result.empty:
        return pd.DataFrame(
            [
                {
                    "comparison": "EcoNiche-Opt-ImmuneComposite_vs_baseline",
                    "status": "RESULT_PENDING",
                    "reason": "no comparable baseline rows",
                }
            ]
        )
    result["fdr_q"] = benjamini_hochberg(result["p_value"].fillna(1.0))
    return result.sort_values(["baseline_AUROC", "comparison"], ascending=[False, True])


def _compare_with_ml_baselines(predictions: pd.DataFrame, ml_predictions_path: Path) -> pd.DataFrame:
    if not ml_predictions_path.exists():
        return pd.DataFrame(
            [
                {
                    "comparison": "EcoNiche-Opt-ImmuneComposite_vs_ML_baseline",
                    "status": "RESULT_PENDING",
                    "reason": "missing ML baseline predictions",
                }
            ]
        )
    baseline = pd.read_csv(ml_predictions_path, sep="\t")
    rows = []
    for model_name, frame in baseline.groupby("model_name"):
        merged = predictions.merge(
            frame[["sample_id", "cohort", "response_probability"]],
            on=["sample_id", "cohort"],
            suffixes=("_optimized", "_baseline"),
            how="inner",
        ).dropna(subset=["response_probability_optimized", "response_probability_baseline", "true_response_label"])
        if len(merged) < 5 or merged["true_response_label"].nunique() < 2:
            continue
        stats = paired_bootstrap_delta(
            merged["true_response_label"],
            merged["response_probability_optimized"],
            merged["response_probability_baseline"],
            n_bootstrap=1000,
        )
        optimized_metrics = compute_binary_metrics(merged["true_response_label"], merged["response_probability_optimized"])
        baseline_metrics = compute_binary_metrics(merged["true_response_label"], merged["response_probability_baseline"])
        rows.append(
            {
                **stats,
                "comparison": f"EcoNiche-Opt-ImmuneComposite_vs_{model_name}",
                "metric": "response_AUROC",
                "n_samples": len(merged),
                "optimized_AUROC": optimized_metrics["AUROC"],
                "baseline_AUROC": baseline_metrics["AUROC"],
                "status": "computed_from_pipeline",
            }
        )
    result = pd.DataFrame(rows)
    if result.empty:
        return pd.DataFrame(
            [
                {
                    "comparison": "EcoNiche-Opt-ImmuneComposite_vs_ML_baseline",
                    "status": "RESULT_PENDING",
                    "reason": "no comparable ML baseline rows",
                }
            ]
        )
    result["fdr_q"] = benjamini_hochberg(result["p_value"].fillna(1.0))
    return result.sort_values(["baseline_AUROC", "comparison"], ascending=[False, True])


def _cohortwise_metric_summary(
    predictions: pd.DataFrame,
    probability_column: str,
    label_column: str,
    group_column: str = "cohort",
) -> tuple[float, pd.DataFrame]:
    rows = []
    for cohort, frame in predictions.groupby(group_column):
        frame = frame.dropna(subset=[probability_column, label_column])
        if len(frame) < 2 or frame[label_column].nunique() < 2:
            continue
        metrics = compute_binary_metrics(frame[label_column], frame[probability_column])
        rows.append({"cohort": cohort, **metrics})
    metrics_by_cohort = pd.DataFrame(rows)
    mean_auroc = float(metrics_by_cohort["AUROC"].mean()) if not metrics_by_cohort.empty else float("nan")
    return mean_auroc, metrics_by_cohort


def _compare_with_current_model(predictions: pd.DataFrame, current_predictions_path: Path) -> pd.DataFrame:
    if not current_predictions_path.exists():
        return pd.DataFrame(
            [
                {
                    "comparison": "EcoNiche-Opt-ImmuneComposite_vs_current_full_model",
                    "status": "RESULT_PENDING",
                    "reason": "missing current model predictions",
                }
            ]
        )
    current = pd.read_csv(current_predictions_path, sep="\t").dropna(subset=["pred_prob", "true_label"]).copy()
    current["true_response_label"] = 1 - current["true_label"].astype(int)
    current["current_response_probability"] = 1.0 - current["pred_prob"].astype(float)
    merged = predictions.merge(
        current[["sample_id", "cohort", "true_response_label", "current_response_probability"]],
        on=["sample_id", "cohort", "true_response_label"],
        how="inner",
    )
    if len(merged) < 5 or merged["true_response_label"].nunique() < 2:
        return pd.DataFrame(
            [
                {
                    "comparison": "EcoNiche-Opt-ImmuneComposite_vs_current_full_model",
                    "status": "RESULT_PENDING",
                    "reason": "no comparable current model rows",
                }
            ]
        )
    stats = paired_bootstrap_delta(
        merged["true_response_label"],
        merged["response_probability"],
        merged["current_response_probability"],
        n_bootstrap=1000,
    )
    optimized_mean, _ = _cohortwise_metric_summary(merged, "response_probability", "true_response_label")
    current_mean, _ = _cohortwise_metric_summary(merged, "current_response_probability", "true_response_label")
    optimized_global = compute_binary_metrics(merged["true_response_label"], merged["response_probability"])
    current_global = compute_binary_metrics(merged["true_response_label"], merged["current_response_probability"])
    return pd.DataFrame(
        [
            {
                **stats,
                "comparison": "EcoNiche-Opt-ImmuneComposite_vs_current_full_model",
                "metric": "response_AUROC",
                "n_samples": len(merged),
                "optimized_global_AUROC": optimized_global["AUROC"],
                "current_global_AUROC": current_global["AUROC"],
                "optimized_mean_lodo_AUROC": optimized_mean,
                "current_mean_lodo_AUROC": current_mean,
                "mean_lodo_delta": optimized_mean - current_mean,
                "status": "computed_from_pipeline",
            }
        ]
    )


def _write_ml_lodo_comparison(predictions: pd.DataFrame, ml_predictions_path: Path, out: Path) -> None:
    if not ml_predictions_path.exists():
        pd.DataFrame(
            [
                {
                    "comparison": "EcoNiche-Opt-ImmuneComposite_vs_ML_baseline",
                    "status": "RESULT_PENDING",
                    "reason": "missing ML baseline predictions",
                }
            ]
        ).to_csv(out, sep="\t", index=False)
        return
    optimized_mean, _ = _cohortwise_metric_summary(predictions, "response_probability", "true_response_label")
    ml = pd.read_csv(ml_predictions_path, sep="\t")
    rows = []
    for model_name, frame in ml.groupby("model_name"):
        baseline_mean, baseline_by_cohort = _cohortwise_metric_summary(frame, "response_probability", "true_response_label")
        rows.append(
            {
                "comparison": f"EcoNiche-Opt-ImmuneComposite_vs_{model_name}",
                "metric": "mean_lodo_response_AUROC",
                "optimized_mean_lodo_AUROC": optimized_mean,
                "baseline_mean_lodo_AUROC": baseline_mean,
                "mean_lodo_delta": optimized_mean - baseline_mean,
                "n_cohorts": int(len(baseline_by_cohort)),
                "status": "computed_from_pipeline",
            }
        )
    pd.DataFrame(rows).sort_values("baseline_mean_lodo_AUROC", ascending=False).to_csv(out, sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--baseline-predictions", default="results/real/baseline_predictions.tsv")
    parser.add_argument("--ml-baseline-predictions", default="results/real/ml_baseline_predictions.tsv")
    parser.add_argument("--current-predictions", default="results/real/lodo_predictions.tsv")
    parser.add_argument("--out", default="results/real_optimized")
    parser.add_argument("--include-demo", action="store_true")
    parser.add_argument("--preferred-candidate", default=PREFERRED_CANDIDATE)
    parser.add_argument("--preferred-tolerance", type=float, default=0.02)
    args = parser.parse_args()

    X_by_cohort, y_by_cohort, metadata_by_cohort = load_processed_bulk(ROOT / args.processed_dir)
    cohorts = sorted(X_by_cohort)
    if not args.include_demo:
        cohorts = [cohort for cohort in cohorts if not cohort.startswith("demo_cohort_")]
    cohorts = [cohort for cohort in cohorts if y_by_cohort[cohort].nunique() >= 2]
    if len(cohorts) < 2:
        raise SystemExit("Need at least two labeled cohorts for nested response-composite LODO.")

    result = run_nested_response_composite(
        X_by_cohort,
        y_by_cohort,
        metadata_by_cohort,
        cohorts=cohorts,
        preferred_candidate=args.preferred_candidate,
        preferred_tolerance=args.preferred_tolerance,
    )

    out = ROOT / args.out
    out.mkdir(parents=True, exist_ok=True)
    result.metrics.to_csv(out / "optimized_lodo_response_metrics.tsv", sep="\t", index=False)
    result.predictions.to_csv(out / "optimized_lodo_response_predictions.tsv", sep="\t", index=False)
    result.inner_selection.to_csv(out / "optimized_inner_selection.tsv", sep="\t", index=False)
    _write_global_metrics(result.predictions, result.metrics, out / "optimized_global_response_metrics.tsv")
    comparison = _compare_with_baselines(result.predictions, ROOT / args.baseline_predictions)
    comparison.to_csv(out / "optimized_vs_baselines_response.tsv", sep="\t", index=False)
    ml_comparison = _compare_with_ml_baselines(result.predictions, ROOT / args.ml_baseline_predictions)
    ml_comparison.to_csv(out / "optimized_vs_ml_baselines_response.tsv", sep="\t", index=False)
    _write_ml_lodo_comparison(result.predictions, ROOT / args.ml_baseline_predictions, out / "optimized_vs_ml_baselines_lodo.tsv")
    current_comparison = _compare_with_current_model(result.predictions, ROOT / args.current_predictions)
    current_comparison.to_csv(out / "optimized_vs_current_model_response.tsv", sep="\t", index=False)
    print(result.metrics[["cohort", "AUROC", "AUPRC", "balanced_accuracy", "selected_model", "selected_by"]].to_string(index=False))
    print(f"Wrote optimized response-composite outputs to {out}")


if __name__ == "__main__":
    main()
