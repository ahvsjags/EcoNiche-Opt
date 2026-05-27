from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.metrics import compute_binary_metrics


PRIMARY_MODEL = "EcoNiche-Opt-HeuristicEcology"
PRIMARY_ENDPOINT = "primary_recist"
PRIMARY_STRATA = ["melanoma_core_high_evidence", "melanoma_recist_supported_primary"]
MIDRANGE_FIXED_POLICIES = ["fixed_0.40", "fixed_0.50", "fixed_0.60"]


def _read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t")


def _parse_fixed_policy(policy: str) -> float:
    if not policy.startswith("fixed_"):
        raise ValueError(f"Not a fixed-threshold policy: {policy}")
    return float(policy.removeprefix("fixed_"))


def _best_training_threshold(frame: pd.DataFrame) -> float:
    y = frame["true_response_label"].astype(int).to_numpy()
    p = frame["response_probability"].astype(float).to_numpy()
    if len(np.unique(y)) < 2:
        return 0.5
    best_threshold = 0.5
    best_score = -math.inf
    for threshold in np.unique(np.clip(p, 0.0, 1.0)):
        score = compute_binary_metrics(y, p, threshold=float(threshold))["balanced_accuracy"]
        if score > best_score:
            best_score = float(score)
            best_threshold = float(threshold)
    return best_threshold


def _prevalence_quantile_threshold(frame: pd.DataFrame) -> float:
    y = frame["true_response_label"].astype(int).to_numpy()
    p = frame["response_probability"].astype(float).to_numpy()
    if len(p) == 0:
        return 0.5
    prevalence = float(np.mean(y))
    return float(np.quantile(p, np.clip(1.0 - prevalence, 0.0, 1.0)))


def threshold_for_policy(policy: str, train: pd.DataFrame) -> float:
    if policy.startswith("fixed_"):
        return _parse_fixed_policy(policy)
    if policy == "training_youden":
        return _best_training_threshold(train)
    if policy == "training_prevalence_quantile":
        return _prevalence_quantile_threshold(train)
    raise ValueError(f"Unsupported threshold policy: {policy}")


def _metrics_row(frame: pd.DataFrame, threshold_values: pd.Series | float) -> dict[str, object]:
    y = frame["true_response_label"].astype(int)
    p = frame["response_probability"].astype(float)
    if isinstance(threshold_values, pd.Series):
        thresholds = threshold_values.astype(float)
        pred = (p >= thresholds).astype(int)
        threshold_min = float(thresholds.min())
        threshold_median = float(thresholds.median())
        threshold_max = float(thresholds.max())
    else:
        pred = (p >= float(threshold_values)).astype(int)
        threshold_min = threshold_median = threshold_max = float(threshold_values)
    metrics = compute_binary_metrics(y, p, threshold=0.5)
    metrics["balanced_accuracy"] = compute_binary_metrics(y, pred, threshold=0.5)["balanced_accuracy"]
    metrics["accuracy"] = compute_binary_metrics(y, pred, threshold=0.5)["accuracy"]
    metrics["MCC"] = compute_binary_metrics(y, pred, threshold=0.5)["MCC"]
    metrics["F1"] = compute_binary_metrics(y, pred, threshold=0.5)["F1"]
    tn_fp_fn_tp = _confusion_counts(y, pred)
    metrics.update(tn_fp_fn_tp)
    metrics["sensitivity"] = _ratio(tn_fp_fn_tp["tp"], tn_fp_fn_tp["tp"] + tn_fp_fn_tp["fn"])
    metrics["specificity"] = _ratio(tn_fp_fn_tp["tn"], tn_fp_fn_tp["tn"] + tn_fp_fn_tp["fp"])
    metrics["PPV"] = _ratio(tn_fp_fn_tp["tp"], tn_fp_fn_tp["tp"] + tn_fp_fn_tp["fp"])
    metrics["NPV"] = _ratio(tn_fp_fn_tp["tn"], tn_fp_fn_tp["tn"] + tn_fp_fn_tp["fn"])
    metrics.update(
        {
            "n_samples": int(len(frame)),
            "n_responders": int(y.sum()),
            "n_nonresponders": int((y == 0).sum()),
            "response_prevalence": float(y.mean()),
            "threshold_min": threshold_min,
            "threshold_median": threshold_median,
            "threshold_max": threshold_max,
            "n_predicted_high": int(pred.sum()),
            "AUPRC_minus_prevalence": float(metrics["AUPRC"] - y.mean()),
        }
    )
    return metrics


def _ratio(num: int, denom: int) -> float:
    return float(num / denom) if denom else math.nan


def _confusion_counts(y: pd.Series, pred: pd.Series) -> dict[str, int]:
    y_int = y.astype(int)
    pred_int = pred.astype(int)
    return {
        "tn": int(((y_int == 0) & (pred_int == 0)).sum()),
        "fp": int(((y_int == 0) & (pred_int == 1)).sum()),
        "fn": int(((y_int == 1) & (pred_int == 0)).sum()),
        "tp": int(((y_int == 1) & (pred_int == 1)).sum()),
    }


def _inner_policy_score(train: pd.DataFrame, policy: str) -> tuple[float, float, int]:
    scores: list[float] = []
    for validation_fold in sorted(train["fold"].astype(str).unique()):
        inner_train = train[~train["fold"].astype(str).eq(validation_fold)].copy()
        validation = train[train["fold"].astype(str).eq(validation_fold)].copy()
        if inner_train.empty or validation.empty or validation["true_response_label"].nunique() < 2:
            continue
        threshold = threshold_for_policy(policy, inner_train)
        pred = (validation["response_probability"].astype(float) >= threshold).astype(int)
        score = compute_binary_metrics(validation["true_response_label"].astype(int), pred, threshold=0.5)["balanced_accuracy"]
        scores.append(float(score))
    if not scores:
        return -math.inf, math.inf, 0
    return float(np.mean(scores)), float(np.std(scores, ddof=0)), len(scores)


def select_midrange_policy(train: pd.DataFrame) -> tuple[str, pd.DataFrame]:
    rows = []
    for policy in MIDRANGE_FIXED_POLICIES:
        mean_ba, sd_ba, n_inner = _inner_policy_score(train, policy)
        selection_score = mean_ba - 0.05 * sd_ba
        rows.append(
            {
                "candidate_policy": policy,
                "inner_mean_balanced_accuracy": mean_ba,
                "inner_sd_balanced_accuracy": sd_ba,
                "inner_n_validations": n_inner,
                "selection_score": selection_score,
            }
        )
    frame = pd.DataFrame(rows)
    best = frame.sort_values(["selection_score", "candidate_policy"], ascending=[False, True]).iloc[0]
    return str(best["candidate_policy"]), frame


def evaluate_primary_thresholds(predictions: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    frame = predictions[
        predictions["endpoint"].astype(str).eq(PRIMARY_ENDPOINT)
        & predictions["model_name"].astype(str).eq(PRIMARY_MODEL)
        & predictions["stratum"].astype(str).isin(PRIMARY_STRATA)
    ].copy()
    summary_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    selection_rows: list[pd.DataFrame] = []
    for stratum, stratum_frame in frame.groupby("stratum", sort=True):
        stratum_frame = stratum_frame.copy()
        original_thresholds = stratum_frame["threshold"].astype(float)
        summary_rows.append(
            {
                "endpoint": PRIMARY_ENDPOINT,
                "stratum": stratum,
                "threshold_policy": "original_fold_threshold",
                "selection_boundary": "outer_lodo_training_threshold_from_source_pipeline",
                **_metrics_row(stratum_frame, original_thresholds),
            }
        )
        for policy in [*MIDRANGE_FIXED_POLICIES, "training_youden", "training_prevalence_quantile", "nested_midrange_fixed_grid"]:
            parts = []
            selected_policies = []
            for holdout in sorted(stratum_frame["fold"].astype(str).unique()):
                train = stratum_frame[~stratum_frame["fold"].astype(str).eq(holdout)].copy()
                test = stratum_frame[stratum_frame["fold"].astype(str).eq(holdout)].copy()
                if policy == "nested_midrange_fixed_grid":
                    selected_policy, selection = select_midrange_policy(train)
                    selection.insert(0, "endpoint", PRIMARY_ENDPOINT)
                    selection.insert(1, "stratum", stratum)
                    selection.insert(2, "holdout", holdout)
                    selection["outer_policy"] = policy
                    selection["selected_policy"] = selected_policy
                    selection_rows.append(selection)
                    threshold = threshold_for_policy(selected_policy, train)
                    selected_policies.append(selected_policy)
                else:
                    threshold = threshold_for_policy(policy, train)
                    selected_policies.append(policy)
                test = test.copy()
                test["recalibrated_threshold"] = float(threshold)
                test["recalibrated_pred_response_label"] = (
                    test["response_probability"].astype(float) >= float(threshold)
                ).astype(int)
                test["threshold_policy"] = policy
                test["selected_threshold_policy"] = selected_policies[-1]
                parts.append(test)
            scored = pd.concat(parts, ignore_index=True)
            summary_rows.append(
                {
                    "endpoint": PRIMARY_ENDPOINT,
                    "stratum": stratum,
                    "threshold_policy": policy,
                    "selection_boundary": "outer_lodo_training_cohorts_only",
                    "selected_policies": ",".join(sorted(set(selected_policies))),
                    **_metrics_row(scored, scored["recalibrated_threshold"]),
                }
            )
            prediction_rows.extend(scored.to_dict("records"))
    selection_df = pd.concat(selection_rows, ignore_index=True) if selection_rows else pd.DataFrame()
    return pd.DataFrame(summary_rows), pd.DataFrame(prediction_rows), selection_df


def write_markdown(summary: pd.DataFrame, out_md: Path) -> None:
    lines = [
        "# Threshold Recalibration Audit",
        "",
        "This audit evaluates threshold policies without changing the rank ordering of EcoNiche-Opt scores.",
        "For LODO rows, each holdout threshold is selected from training cohorts only.",
        "The nested midrange policy selects among fixed 0.40, 0.50 and 0.60 using inner validation on the outer-training cohorts.",
        "",
    ]
    for _, row in summary.iterrows():
        lines.append(
            "- `{}` `{}`: BA={:.3f}, AUROC={:.3f}, AUPRC={:.3f}, threshold_median={:.3f}, boundary={}".format(
                row["stratum"],
                row["threshold_policy"],
                float(row["balanced_accuracy"]),
                float(row["AUROC"]),
                float(row["AUPRC"]),
                float(row["threshold_median"]),
                row["selection_boundary"],
            )
        )
    out_md.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--predictions",
        default="results/endpoint_modules_heuristic_deep_primary_20260519/endpoint_module_predictions.tsv",
    )
    parser.add_argument("--out", default="results/threshold_recalibration_audit_20260527")
    args = parser.parse_args()

    out = ROOT / args.out
    out.mkdir(parents=True, exist_ok=True)
    predictions = _read_tsv(ROOT / args.predictions)
    summary, recalibrated_predictions, selection = evaluate_primary_thresholds(predictions)
    summary.to_csv(out / "threshold_recalibration_primary_summary.tsv", sep="\t", index=False)
    recalibrated_predictions.to_csv(out / "threshold_recalibration_primary_predictions.tsv", sep="\t", index=False)
    selection.to_csv(out / "threshold_recalibration_policy_selection.tsv", sep="\t", index=False)
    write_markdown(summary, out / "THRESHOLD_RECALIBRATION_AUDIT.md")
    best = summary.sort_values(["stratum", "balanced_accuracy"], ascending=[True, False]).groupby("stratum").head(1)
    print(json.dumps(best[["stratum", "threshold_policy", "AUROC", "AUPRC", "balanced_accuracy"]].to_dict("records")))
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
