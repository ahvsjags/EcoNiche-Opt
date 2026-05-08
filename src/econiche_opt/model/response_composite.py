from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from econiche.baselines import BASELINE_SIGNATURES, signature_score
from econiche.metrics import compute_binary_metrics


RESPONSE_HIGH_SIGNATURES = [
    "TIG",
    "TIDE_dysfunction",
    "IFNG",
    "PDCD1LG2",
    "CXCL9",
    "MCP_CD8_T",
    "APM",
    "PDL1_CD274",
    "IMPRES_template",
    "HLA_DRA",
    "CTLA4",
    "CYT",
    "PDCD1",
    "TLS",
]

NONRESPONSE_HIGH_SIGNATURES = [
    "MCP_fibroblast",
    "MPS",
    "IPRES",
    "C_ECM",
    "TIDE_exclusion",
    "ESCS",
]

PREFERRED_CANDIDATE = "ifn_core_pdcd1lg2_weighted"


@dataclass(frozen=True)
class ResponseCompositeResult:
    metrics: pd.DataFrame
    predictions: pd.DataFrame
    inner_selection: pd.DataFrame


def response_labels_from_nonresponse(nonresponse_labels: pd.Series) -> pd.Series:
    labels = pd.Series(nonresponse_labels, index=nonresponse_labels.index).astype(int)
    return 1 - labels


def _zscore(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    sd = values.std(ddof=0)
    if pd.isna(sd) or sd <= 0:
        return pd.Series(0.0, index=series.index)
    return ((values - values.mean()) / sd).fillna(0.0)


def _sigmoid(values: pd.Series | np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    return 1.0 / (1.0 + np.exp(-arr))


def build_signature_features(X: pd.DataFrame) -> pd.DataFrame:
    signatures = RESPONSE_HIGH_SIGNATURES + NONRESPONSE_HIGH_SIGNATURES
    features: dict[str, pd.Series] = {}
    for name in signatures:
        genes = BASELINE_SIGNATURES.get(name, [])
        features[name] = _zscore(signature_score(X, genes))
    return pd.DataFrame(features, index=X.index).fillna(0.0)


def build_candidate_scores(features: pd.DataFrame) -> dict[str, pd.Series]:
    response_available = [name for name in RESPONSE_HIGH_SIGNATURES if name in features.columns]
    nonresponse_available = [name for name in NONRESPONSE_HIGH_SIGNATURES if name in features.columns]
    top5 = ["TIG", "TIDE_dysfunction", "IFNG", "PDCD1LG2", "CXCL9"]
    checkpoint = ["TIG", "TIDE_dysfunction", "IFNG", "PDCD1LG2", "CXCL9", "MCP_CD8_T", "PDCD1", "CTLA4", "CYT"]

    scores: dict[str, pd.Series] = {
        "immune_top5": features[[name for name in top5 if name in features.columns]].mean(axis=1),
        "ifn_core_pdcd1lg2_weighted": (features["IFNG"] + features["CXCL9"] + 2.0 * features["PDCD1LG2"]) / 4.0,
        "tcell_checkpoint": features[[name for name in checkpoint if name in features.columns]].mean(axis=1),
        "immune_mean": features[response_available].mean(axis=1),
        "immune_minus_suppressive": features[response_available].mean(axis=1) - features[nonresponse_available].mean(axis=1),
        "top5_minus_stroma_myeloid": features[[name for name in top5 if name in features.columns]].mean(axis=1)
        - features[nonresponse_available].mean(axis=1),
        "suppression_inverse": -features[nonresponse_available].mean(axis=1),
    }
    for name in RESPONSE_HIGH_SIGNATURES:
        if name in features.columns:
            scores[f"signed_{name}"] = features[name]
    for name in NONRESPONSE_HIGH_SIGNATURES:
        if name in features.columns:
            scores[f"signed_{name}"] = -features[name]
    return {name: score.fillna(0.0) for name, score in scores.items()}


def _select_threshold(y_true: np.ndarray, prob: np.ndarray) -> tuple[float, float]:
    if len(y_true) == 0:
        return 0.5, float("nan")
    thresholds = np.unique(np.quantile(prob, np.linspace(0.05, 0.95, 19)))
    best_threshold = 0.5
    best_score = -np.inf
    for threshold in thresholds:
        metrics = compute_binary_metrics(y_true, prob, threshold=float(threshold))
        score = float(metrics.get("balanced_accuracy", float("nan")))
        if np.isfinite(score) and score > best_score:
            best_threshold = float(threshold)
            best_score = score
    return best_threshold, float(best_score)


def _evaluate_candidate(
    train_cohorts: list[str],
    candidate: str,
    candidate_scores: dict[str, dict[str, pd.Series]],
    y_response_by_cohort: dict[str, pd.Series],
) -> dict[str, object]:
    aurocs = []
    auprcs = []
    pooled_y: list[int] = []
    pooled_prob: list[float] = []
    for holdout in train_cohorts:
        y_true = y_response_by_cohort[holdout].astype(int)
        prob = _sigmoid(candidate_scores[holdout][candidate])
        metrics = compute_binary_metrics(y_true, prob)
        aurocs.append(metrics["AUROC"])
        auprcs.append(metrics["AUPRC"])
        pooled_y.extend(y_true.tolist())
        pooled_prob.extend(prob.tolist())
    threshold, balanced_accuracy = _select_threshold(np.asarray(pooled_y, dtype=int), np.asarray(pooled_prob, dtype=float))
    inner_mean_auroc = float(np.nanmean(aurocs))
    inner_mean_auprc = float(np.nanmean(auprcs))
    selection_score = inner_mean_auroc + 0.05 * inner_mean_auprc + 0.05 * balanced_accuracy
    return {
        "candidate": candidate,
        "inner_mean_AUROC": inner_mean_auroc,
        "inner_mean_AUPRC": inner_mean_auprc,
        "inner_balanced_accuracy": balanced_accuracy,
        "threshold": threshold,
        "selection_score": float(selection_score),
    }


def run_nested_response_composite(
    X_by_cohort: dict[str, pd.DataFrame],
    y_nonresponse_by_cohort: dict[str, pd.Series],
    metadata_by_cohort: dict[str, pd.DataFrame],
    cohorts: list[str] | None = None,
    preferred_candidate: str = PREFERRED_CANDIDATE,
    preferred_tolerance: float = 0.02,
) -> ResponseCompositeResult:
    active_cohorts = sorted(cohorts or X_by_cohort)
    features_by_cohort = {cohort: build_signature_features(X_by_cohort[cohort]) for cohort in active_cohorts}
    candidate_scores = {cohort: build_candidate_scores(features_by_cohort[cohort]) for cohort in active_cohorts}
    y_response = {
        cohort: response_labels_from_nonresponse(y_nonresponse_by_cohort[cohort]).reindex(X_by_cohort[cohort].index)
        for cohort in active_cohorts
    }

    metric_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    inner_rows: list[dict[str, object]] = []
    for holdout in active_cohorts:
        train_cohorts = [cohort for cohort in active_cohorts if cohort != holdout]
        candidates = sorted(candidate_scores[holdout])
        scored = [
            _evaluate_candidate(train_cohorts, candidate, candidate_scores, y_response)
            for candidate in candidates
            if all(candidate in candidate_scores[cohort] for cohort in train_cohorts)
        ]
        for row in scored:
            inner_rows.append({"holdout": holdout, **row})
        best = max(scored, key=lambda row: float(row["selection_score"]))
        preferred = next((row for row in scored if row["candidate"] == preferred_candidate), None)
        selected_by = "best_inner"
        if preferred is not None and float(best["selection_score"]) - float(preferred["selection_score"]) <= preferred_tolerance:
            best = preferred
            selected_by = "preferred_within_tolerance"

        y_true = y_response[holdout].astype(int)
        prob = _sigmoid(candidate_scores[holdout][str(best["candidate"])])
        threshold = float(best["threshold"])
        metrics = compute_binary_metrics(y_true, prob, threshold=threshold)
        metrics_05 = compute_binary_metrics(y_true, prob, threshold=0.5)
        metric_rows.append(
            {
                **metrics,
                "balanced_accuracy_0_5": metrics_05["balanced_accuracy"],
                "cohort": holdout,
                "model_name": "EcoNiche-Opt-ImmuneComposite",
                "endpoint": "response_positive_primary_recist",
                "n_samples": len(y_true),
                "n_responders": int((y_true == 1).sum()),
                "n_nonresponders": int((y_true == 0).sum()),
                "selected_model": best["candidate"],
                "selected_by": selected_by,
                "threshold": threshold,
                "inner_mean_AUROC": best["inner_mean_AUROC"],
                "inner_mean_AUPRC": best["inner_mean_AUPRC"],
                "inner_balanced_accuracy": best["inner_balanced_accuracy"],
                "selection_score": best["selection_score"],
                "preferred_tolerance": preferred_tolerance,
            }
        )

        metadata = metadata_by_cohort[holdout].reindex(X_by_cohort[holdout].index)
        pred_label = (prob >= threshold).astype(int)
        for idx, sample_id in enumerate(X_by_cohort[holdout].index):
            meta = metadata.loc[sample_id] if sample_id in metadata.index else pd.Series(dtype=object)
            response_label = int(y_true.loc[sample_id])
            prediction_rows.append(
                {
                    "sample_id": meta.get("sample_id", sample_id),
                    "patient_id": meta.get("patient_id", pd.NA),
                    "cohort": holdout,
                    "true_label_nonresponse": int(1 - response_label),
                    "true_response_label": response_label,
                    "response_probability": float(prob[idx]),
                    "nonresponse_probability": float(1.0 - prob[idx]),
                    "pred_response_label": int(pred_label[idx]),
                    "model_name": "EcoNiche-Opt-ImmuneComposite",
                    "fold": holdout,
                    "selected_model": best["candidate"],
                    "threshold": threshold,
                }
            )
    return ResponseCompositeResult(
        metrics=pd.DataFrame(metric_rows),
        predictions=pd.DataFrame(prediction_rows),
        inner_selection=pd.DataFrame(inner_rows),
    )
