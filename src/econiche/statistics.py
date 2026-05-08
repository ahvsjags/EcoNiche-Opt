from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn import metrics as sk_metrics


def make_lodo_folds(metadata: pd.DataFrame, cohort_column: str = "cohort") -> dict[str, dict[str, pd.DataFrame]]:
    if cohort_column not in metadata.columns:
        raise ValueError(f"Missing cohort column: {cohort_column}")
    if "label" not in metadata.columns:
        raise ValueError("Missing label column")

    folds: dict[str, dict[str, pd.DataFrame]] = {}
    labeled = metadata[metadata["label"].notna()].copy()
    for holdout in sorted(labeled[cohort_column].dropna().unique()):
        test = labeled[labeled[cohort_column] == holdout].copy()
        holdout_patients = set(test.get("patient_id", pd.Series(dtype=str)).astype(str))
        train = labeled[labeled[cohort_column] != holdout].copy()
        if "patient_id" in train.columns:
            train = train[~train["patient_id"].astype(str).isin(holdout_patients)].copy()
        if train.empty or test.empty:
            continue
        if train["label"].nunique() < 2 or test["label"].nunique() < 2:
            continue
        folds[str(holdout)] = {"train": train.reset_index(drop=True), "test": test.reset_index(drop=True)}
    return folds


def benjamini_hochberg(p_values) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    ranked = np.empty_like(p)
    n = len(p)
    cumulative = 1.0
    for rank, idx in reversed(list(enumerate(order, start=1))):
        cumulative = min(cumulative, p[idx] * n / rank)
        ranked[idx] = cumulative
    return np.clip(ranked, 0, 1)


def paired_bootstrap_delta(
    y_true,
    prob_a,
    prob_b,
    metric: str = "AUROC",
    n_bootstrap: int = 1000,
    random_state: int = 42,
) -> dict[str, float]:
    y = np.asarray(y_true)
    a = np.asarray(prob_a)
    b = np.asarray(prob_b)
    rng = np.random.default_rng(random_state)
    deltas = []
    for _ in range(n_bootstrap):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        if metric.upper() == "AUPRC":
            da = sk_metrics.average_precision_score(y[idx], a[idx])
            db = sk_metrics.average_precision_score(y[idx], b[idx])
        else:
            da = sk_metrics.roc_auc_score(y[idx], a[idx])
            db = sk_metrics.roc_auc_score(y[idx], b[idx])
        deltas.append(da - db)
    if not deltas:
        return {"mean_delta": float("nan"), "ci_low": float("nan"), "ci_high": float("nan"), "p_value": float("nan")}
    arr = np.asarray(deltas)
    p_value = 2 * min(float((arr <= 0).mean()), float((arr >= 0).mean()))
    return {
        "mean_delta": float(arr.mean()),
        "ci_low": float(np.quantile(arr, 0.025)),
        "ci_high": float(np.quantile(arr, 0.975)),
        "p_value": min(p_value, 1.0),
    }
