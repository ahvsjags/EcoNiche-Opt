from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn import metrics as sk_metrics
from sklearn.linear_model import LogisticRegression


def _as_arrays(y_true, prob) -> tuple[np.ndarray, np.ndarray]:
    y = np.asarray(y_true, dtype=float)
    p = np.asarray(prob, dtype=float)
    mask = np.isfinite(y) & np.isfinite(p)
    return y[mask].astype(int), np.clip(p[mask], 1e-6, 1 - 1e-6)


def expected_calibration_error(y_true, prob, n_bins: int = 10) -> float:
    y, p = _as_arrays(y_true, prob)
    if len(y) == 0:
        return float("nan")
    bins = np.linspace(0, 1, n_bins + 1)
    ece = 0.0
    for left, right in zip(bins[:-1], bins[1:]):
        if right == 1:
            mask = (p >= left) & (p <= right)
        else:
            mask = (p >= left) & (p < right)
        if mask.any():
            ece += mask.mean() * abs(float(y[mask].mean()) - float(p[mask].mean()))
    return float(ece)


def calibration_bins(y_true, prob, n_bins: int = 10) -> pd.DataFrame:
    y, p = _as_arrays(y_true, prob)
    rows = []
    bins = np.linspace(0, 1, n_bins + 1)
    for idx, (left, right) in enumerate(zip(bins[:-1], bins[1:])):
        mask = (p >= left) & (p <= right) if right == 1 else (p >= left) & (p < right)
        rows.append(
            {
                "bin": idx,
                "prob_low": left,
                "prob_high": right,
                "n": int(mask.sum()),
                "prob_mean": float(p[mask].mean()) if mask.any() else float("nan"),
                "event_rate": float(y[mask].mean()) if mask.any() else float("nan"),
            }
        )
    return pd.DataFrame(rows)


def calibration_slope_intercept(y_true, prob) -> tuple[float, float]:
    y, p = _as_arrays(y_true, prob)
    if len(np.unique(y)) < 2 or len(y) < 4:
        return float("nan"), float("nan")
    logit = np.log(p / (1 - p)).reshape(-1, 1)
    try:
        model = LogisticRegression(C=1e6, solver="lbfgs").fit(logit, y)
    except ValueError:
        return float("nan"), float("nan")
    return float(model.coef_[0, 0]), float(model.intercept_[0])


def compute_binary_metrics(y_true, prob, threshold: float = 0.5) -> dict[str, float]:
    y, p = _as_arrays(y_true, prob)
    if len(y) == 0:
        return {name: float("nan") for name in ["AUROC", "AUPRC", "balanced_accuracy", "MCC", "F1"]}
    pred = (p >= threshold).astype(int)
    tn, fp, fn, tp = sk_metrics.confusion_matrix(y, pred, labels=[0, 1]).ravel()
    slope, intercept = calibration_slope_intercept(y, p)
    return {
        "AUROC": float(sk_metrics.roc_auc_score(y, p)) if len(np.unique(y)) == 2 else float("nan"),
        "AUPRC": float(sk_metrics.average_precision_score(y, p)) if len(np.unique(y)) == 2 else float("nan"),
        "balanced_accuracy": float(sk_metrics.balanced_accuracy_score(y, pred)),
        "accuracy": float(sk_metrics.accuracy_score(y, pred)),
        "MCC": float(sk_metrics.matthews_corrcoef(y, pred)),
        "F1": float(sk_metrics.f1_score(y, pred, zero_division=0)),
        "sensitivity": float(tp / (tp + fn)) if (tp + fn) else float("nan"),
        "specificity": float(tn / (tn + fp)) if (tn + fp) else float("nan"),
        "PPV": float(tp / (tp + fp)) if (tp + fp) else float("nan"),
        "NPV": float(tn / (tn + fn)) if (tn + fn) else float("nan"),
        "Brier": float(sk_metrics.brier_score_loss(y, p)),
        "ECE": expected_calibration_error(y, p),
        "calibration_slope": slope,
        "calibration_intercept": intercept,
    }


def decision_curve(y_true, prob, thresholds: list[float] | None = None) -> pd.DataFrame:
    y, p = _as_arrays(y_true, prob)
    thresholds = thresholds or [round(x, 2) for x in np.linspace(0.05, 0.95, 19)]
    rows = []
    n = len(y)
    for threshold in thresholds:
        pred = p >= threshold
        tp = int(((y == 1) & pred).sum())
        fp = int(((y == 0) & pred).sum())
        net_benefit = (tp / n) - (fp / n) * (threshold / (1 - threshold)) if n else float("nan")
        rows.append({"threshold": threshold, "net_benefit": float(net_benefit)})
    return pd.DataFrame(rows)
