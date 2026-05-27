from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from sklearn import metrics as sk_metrics

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.metrics import calibration_bins, calibration_slope_intercept, expected_calibration_error
from econiche.statistics import benjamini_hochberg


PRIMARY_MODEL = "EcoNiche-Opt-HeuristicEcology"
LOCKED_EXTERNAL_MODEL = "EcoNiche-Opt-HeuristicEcology-LockedPanel"
STRICT_MELANOMA_EXTERNAL_COHORTS = {"GSE145996", "PHS000452_LIU_LIKE_PRE"}


def _read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t")


def _safe_float(value: object, default: float = math.nan) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def _nonnull_text(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    return str(value).strip()


def classify_therapy_context(treatment: object, analysis_type: object = "") -> str:
    text = f"{_nonnull_text(treatment)} {_nonnull_text(analysis_type)}".lower()
    if any(token in text for token in ["ipilimumab", "ctla4", "ctla-4", "combo", "combination", "+"]):
        return "combo_therapy"
    if any(token in text for token in ["pembro", "pem", "nivo", "nivolumab", "anti-pd-1", "pd1", "pd-1"]):
        return "anti_pd1_monotherapy"
    if text.strip():
        return "other_or_mixed_therapy"
    return "not_recorded"


def _sampling_context(row: pd.Series) -> str:
    timepoint = _nonnull_text(row.get("timepoint", "")).lower()
    if "pre" in timepoint or "baseline" in timepoint:
        return "baseline_only"
    analysis_type = _nonnull_text(row.get("analysis_type", "")).lower()
    cohort = _nonnull_text(row.get("cohort", "")).lower()
    if "_pre" in cohort or "external" in analysis_type or "transfer" in analysis_type:
        return "pretreatment_design_or_panel_transfer"
    return "not_recorded"


def normalize_primary_predictions(primary: pd.DataFrame) -> pd.DataFrame:
    frame = primary[
        primary["model_name"].astype(str).eq(PRIMARY_MODEL)
        & primary["endpoint"].astype(str).eq("primary_recist")
    ].copy()
    frame["source_context"] = "primary_lodo"
    frame["analysis_type"] = frame.get("stratum", "").astype(str)
    frame["evaluation_group"] = frame.get("stratum", "").astype(str)
    frame["locked_or_fold_threshold"] = frame["threshold"].astype(float)
    frame["predicted_response_label"] = frame["pred_response_label"].astype(int)
    frame["threshold_source"] = "training_fold_lodo"
    if "platform" not in frame.columns:
        frame["platform"] = "bulk_RNAseq_or_processed_expression"
    platform_context = frame["platform"].fillna("").astype(str)
    frame["platform_context"] = platform_context.mask(platform_context.eq(""), "bulk_RNAseq_or_processed_expression")
    frame["therapy_context"] = [
        classify_therapy_context(treatment, analysis_type) for treatment, analysis_type in zip(frame.get("treatment", ""), frame["analysis_type"])
    ]
    frame["sampling_context"] = frame.apply(_sampling_context, axis=1)
    return frame


def normalize_external_predictions(external: pd.DataFrame) -> pd.DataFrame:
    frame = external[external["model_name"].astype(str).eq(LOCKED_EXTERNAL_MODEL)].copy()
    frame["source_context"] = "locked_external"
    frame["stratum"] = frame.get("analysis_type", "").astype(str)
    frame["evaluation_group"] = frame.get("analysis_type", "").astype(str)
    frame["locked_or_fold_threshold"] = frame["locked_threshold"].astype(float)
    frame["predicted_response_label"] = frame["predicted_response_label"].astype(int)
    frame["threshold_source"] = frame.get("threshold_source", "discovery_only").astype(str)
    if "platform" not in frame.columns:
        frame["platform"] = ""
    platform_context = frame["platform"].fillna("").astype(str)
    frame["platform_context"] = platform_context.mask(platform_context.eq(""), "bulk_RNAseq_or_processed_expression")
    frame["therapy_context"] = [
        classify_therapy_context(treatment, analysis_type) for treatment, analysis_type in zip(frame.get("treatment", ""), frame["analysis_type"])
    ]
    frame["sampling_context"] = frame.apply(_sampling_context, axis=1)
    return frame


def context_frames(predictions: pd.DataFrame) -> list[tuple[str, str, pd.DataFrame]]:
    contexts: list[tuple[str, str, pd.DataFrame]] = []
    for (source, endpoint, group), frame in predictions.groupby(["source_context", "endpoint", "evaluation_group"], dropna=False):
        label = f"{source}|{endpoint}|{group}"
        contexts.append((label, "standard_context", frame.copy()))
    external_strict = predictions[
        predictions["source_context"].eq("locked_external")
        & predictions["endpoint"].eq("strict_recist")
        & predictions["cohort"].astype(str).isin(STRICT_MELANOMA_EXTERNAL_COHORTS)
    ].copy()
    if not external_strict.empty:
        contexts.append(
            (
                "locked_external|strict_recist|strict_melanoma_pd1_like_pooled",
                "predeclared_strict_external_pool",
                external_strict,
            )
        )
    for (source, endpoint), frame in predictions.groupby(["source_context", "endpoint"], dropna=False):
        contexts.append((f"{source}|{endpoint}|all_target_rows", "endpoint_pooled", frame.copy()))
    return contexts


def metrics_with_thresholds(y_true: pd.Series, prob: pd.Series, thresholds: pd.Series | float) -> dict[str, float]:
    y = pd.to_numeric(y_true, errors="coerce").to_numpy(dtype=float)
    p = pd.to_numeric(prob, errors="coerce").to_numpy(dtype=float)
    if isinstance(thresholds, pd.Series):
        t = pd.to_numeric(thresholds, errors="coerce").to_numpy(dtype=float)
    else:
        t = np.full(len(y), float(thresholds))
    mask = np.isfinite(y) & np.isfinite(p) & np.isfinite(t)
    y = y[mask].astype(int)
    p = np.clip(p[mask], 1e-6, 1 - 1e-6)
    t = t[mask]
    if len(y) == 0:
        return {name: math.nan for name in ["AUROC", "AUPRC", "balanced_accuracy", "sensitivity", "specificity", "PPV", "NPV", "Brier", "ECE"]}
    pred = (p >= t).astype(int)
    tn, fp, fn, tp = sk_metrics.confusion_matrix(y, pred, labels=[0, 1]).ravel()
    slope, intercept = calibration_slope_intercept(y, p)
    two_class = len(np.unique(y)) == 2
    return {
        "n_samples": int(len(y)),
        "n_responders": int(y.sum()),
        "n_nonresponders": int((y == 0).sum()),
        "prevalence": float(y.mean()),
        "n_predicted_high": int(pred.sum()),
        "threshold_min": float(np.min(t)),
        "threshold_median": float(np.median(t)),
        "threshold_max": float(np.max(t)),
        "AUROC": float(sk_metrics.roc_auc_score(y, p)) if two_class else math.nan,
        "AUPRC": float(sk_metrics.average_precision_score(y, p)) if two_class else math.nan,
        "balanced_accuracy": float(sk_metrics.balanced_accuracy_score(y, pred)) if two_class else math.nan,
        "accuracy": float(sk_metrics.accuracy_score(y, pred)),
        "MCC": float(sk_metrics.matthews_corrcoef(y, pred)) if two_class else math.nan,
        "F1": float(sk_metrics.f1_score(y, pred, zero_division=0)),
        "sensitivity": float(tp / (tp + fn)) if (tp + fn) else math.nan,
        "specificity": float(tn / (tn + fp)) if (tn + fp) else math.nan,
        "PPV": float(tp / (tp + fp)) if (tp + fp) else math.nan,
        "NPV": float(tn / (tn + fn)) if (tn + fn) else math.nan,
        "Brier": float(sk_metrics.brier_score_loss(y, p)),
        "ECE": expected_calibration_error(y, p),
        "calibration_slope": slope,
        "calibration_intercept": intercept,
    }


def high_score_enrichment(frame: pd.DataFrame, context_id: str, context_type: str) -> dict[str, object]:
    y = frame["true_response_label"].astype(int)
    high = frame["response_probability"].astype(float) >= frame["locked_or_fold_threshold"].astype(float)
    high_y = y[high]
    low_y = y[~high]
    n_high = int(high.sum())
    n_low = int((~high).sum())
    responders_high = int(high_y.sum())
    responders_low = int(low_y.sum())
    nonresponders_high = int(n_high - responders_high)
    nonresponders_low = int(n_low - responders_low)
    rate_high = float(responders_high / n_high) if n_high else math.nan
    rate_low = float(responders_low / n_low) if n_low else math.nan
    if n_high and n_low:
        _, p_two_sided = fisher_exact(
            [[responders_high, nonresponders_high], [responders_low, nonresponders_low]],
            alternative="two-sided",
        )
        _, p_greater = fisher_exact(
            [[responders_high, nonresponders_high], [responders_low, nonresponders_low]],
            alternative="greater",
        )
    else:
        p_two_sided = math.nan
        p_greater = math.nan
    odds_ratio_haldane = ((responders_high + 0.5) * (nonresponders_low + 0.5)) / (
        (nonresponders_high + 0.5) * (responders_low + 0.5)
    )
    return {
        "context_id": context_id,
        "context_type": context_type,
        "source_context": frame["source_context"].iloc[0],
        "endpoint": frame["endpoint"].iloc[0],
        "model_name": frame["model_name"].iloc[0],
        "n_samples": int(len(frame)),
        "n_high_score": n_high,
        "n_low_score": n_low,
        "responders_high_score": responders_high,
        "responders_low_score": responders_low,
        "response_rate_high_score": rate_high,
        "response_rate_low_score": rate_low,
        "response_rate_difference": rate_high - rate_low if math.isfinite(rate_high) and math.isfinite(rate_low) else math.nan,
        "response_rate_ratio": rate_high / rate_low if math.isfinite(rate_high) and math.isfinite(rate_low) and rate_low > 0 else math.nan,
        "odds_ratio_haldane": float(odds_ratio_haldane),
        "fisher_two_sided_p": float(p_two_sided),
        "fisher_greater_p": float(p_greater),
        "threshold_policy": "locked_or_training_fold_threshold",
    }


def build_threshold_rows(frame: pd.DataFrame, context_id: str, context_type: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    policies: list[tuple[str, pd.Series | float]] = [
        ("locked_or_training_fold_threshold", frame["locked_or_fold_threshold"].astype(float)),
        ("fixed_0.25", 0.25),
        ("fixed_0.50", 0.50),
        ("fixed_0.75", 0.75),
    ]
    for policy, threshold in policies:
        metrics = metrics_with_thresholds(frame["true_response_label"], frame["response_probability"], threshold)
        rows.append(
            {
                "context_id": context_id,
                "context_type": context_type,
                "source_context": frame["source_context"].iloc[0],
                "endpoint": frame["endpoint"].iloc[0],
                "model_name": frame["model_name"].iloc[0],
                "threshold_policy": policy,
                **metrics,
            }
        )
    return rows


def build_subgroup_rows(predictions: pd.DataFrame) -> pd.DataFrame:
    axes = ["cohort", "therapy_context", "platform_context", "sampling_context", "analysis_type"]
    rows: list[dict[str, object]] = []
    for (source, endpoint, group), frame in predictions.groupby(["source_context", "endpoint", "evaluation_group"], dropna=False):
        for axis in axes:
            if axis not in frame.columns:
                continue
            for subgroup_value, sub in frame.groupby(axis, dropna=False):
                if len(sub) < 4:
                    continue
                metrics = metrics_with_thresholds(
                    sub["true_response_label"],
                    sub["response_probability"],
                    sub["locked_or_fold_threshold"].astype(float),
                )
                rows.append(
                    {
                        "source_context": source,
                        "endpoint": endpoint,
                        "evaluation_group": group,
                        "subgroup_axis": axis,
                        "subgroup_value": _nonnull_text(subgroup_value) or "not_recorded",
                        "model_name": sub["model_name"].iloc[0],
                        "threshold_policy": "locked_or_training_fold_threshold",
                        **metrics,
                    }
                )
    return pd.DataFrame(rows)


def build_calibration_rows(contexts: list[tuple[str, str, pd.DataFrame]]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for context_id, context_type, frame in contexts:
        bins = calibration_bins(frame["true_response_label"], frame["response_probability"], n_bins=10)
        for _, row in bins.iterrows():
            rows.append(
                {
                    "context_id": context_id,
                    "context_type": context_type,
                    "source_context": frame["source_context"].iloc[0],
                    "endpoint": frame["endpoint"].iloc[0],
                    "model_name": frame["model_name"].iloc[0],
                    **row.to_dict(),
                }
            )
    return pd.DataFrame(rows)


def build_outputs(primary_predictions: pd.DataFrame, external_predictions: pd.DataFrame) -> dict[str, pd.DataFrame]:
    primary = normalize_primary_predictions(primary_predictions)
    external = normalize_external_predictions(external_predictions)
    predictions = pd.concat([primary, external], ignore_index=True, sort=False)
    contexts = context_frames(predictions)

    enrichment_rows = [high_score_enrichment(frame, context_id, context_type) for context_id, context_type, frame in contexts]
    enrichment = pd.DataFrame(enrichment_rows)
    if not enrichment.empty:
        enrichment["fisher_two_sided_q"] = benjamini_hochberg(enrichment["fisher_two_sided_p"].fillna(1.0))
        enrichment["fisher_greater_q"] = benjamini_hochberg(enrichment["fisher_greater_p"].fillna(1.0))

    threshold_rows: list[dict[str, object]] = []
    for context_id, context_type, frame in contexts:
        threshold_rows.extend(build_threshold_rows(frame, context_id, context_type))

    return {
        "normalized_predictions": predictions,
        "high_score_enrichment": enrichment,
        "threshold_operating_points": pd.DataFrame(threshold_rows),
        "subgroup_metrics": build_subgroup_rows(predictions),
        "calibration_bins": build_calibration_rows(contexts),
    }


def write_audit(outputs: dict[str, pd.DataFrame], out_dir: Path) -> None:
    enrichment = outputs["high_score_enrichment"]
    subgroup = outputs["subgroup_metrics"]
    threshold = outputs["threshold_operating_points"]
    calibration = outputs["calibration_bins"]
    strict = enrichment[enrichment["context_id"].eq("locked_external|strict_recist|strict_melanoma_pd1_like_pooled")]
    primary = enrichment[enrichment["context_id"].eq("primary_lodo|primary_recist|melanoma_core_high_evidence")]
    lines = [
        "# Clinical Interpretability Audit",
        "",
        "This audit is derived from registered prediction tables only. It does not retrain the model, choose features, alter thresholds, or use locked external labels for model selection.",
        "",
        f"- High-score enrichment contexts: {len(enrichment)}.",
        f"- Subgroup metric rows: {len(subgroup)}.",
        f"- Threshold operating-point rows: {len(threshold)}.",
        f"- Calibration-bin rows: {len(calibration)}.",
    ]
    if not primary.empty:
        row = primary.iloc[0]
        lines.append(
            f"- Primary high-evidence melanoma high-score response rate: {row['response_rate_high_score']:.3f} versus {row['response_rate_low_score']:.3f} low-score; Fisher two-sided q={row['fisher_two_sided_q']:.3g}."
        )
    if not strict.empty:
        row = strict.iloc[0]
        lines.append(
            f"- Strict melanoma PD1-like external high-score response rate: {row['response_rate_high_score']:.3f} versus {row['response_rate_low_score']:.3f} low-score; Fisher two-sided q={row['fisher_two_sided_q']:.3g}."
        )
    lines.extend(
        [
            "",
            "Files:",
            "- high_score_enrichment.tsv",
            "- subgroup_metrics.tsv",
            "- threshold_operating_points.tsv",
            "- calibration_bins.tsv",
        ]
    )
    (out_dir / "CLINICAL_INTERPRETABILITY_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--primary-predictions",
        default="results/endpoint_modules_heuristic_deep_primary_20260519/endpoint_module_predictions.tsv",
    )
    parser.add_argument(
        "--external-predictions",
        default="results/locked_external_panel_validation_calibrated_20260519/locked_external_predictions.tsv",
    )
    parser.add_argument("--out", default="results/clinical_interpretability_20260527")
    args = parser.parse_args()

    primary = _read_tsv(ROOT / args.primary_predictions)
    external = _read_tsv(ROOT / args.external_predictions)
    outputs = build_outputs(primary, external)

    out_dir = ROOT / args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, frame in outputs.items():
        frame.to_csv(out_dir / f"{name}.tsv", sep="\t", index=False)
    write_audit(outputs, out_dir)
    print(f"Wrote clinical interpretability audit to {out_dir}")


if __name__ == "__main__":
    main()
