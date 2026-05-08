from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.metrics import balanced_accuracy_score


@dataclass(frozen=True)
class ThresholdResult:
    threshold: float
    metric: str
    training_only: bool


def select_threshold_training_only(y_train, prob_train, metric: str = "balanced_accuracy") -> ThresholdResult:
    y = np.asarray(y_train).astype(int)
    p = np.asarray(prob_train).astype(float)
    if len(y) != len(p):
        raise ValueError("y_train and prob_train must have the same length")
    thresholds = np.unique(np.clip(p, 0.0, 1.0))
    if len(thresholds) == 0:
        return ThresholdResult(0.5, metric, True)
    best_threshold = 0.5
    best_score = -np.inf
    for threshold in thresholds:
        pred = (p >= threshold).astype(int)
        if metric == "balanced_accuracy":
            score = balanced_accuracy_score(y, pred)
        else:
            raise ValueError(f"Unsupported threshold metric: {metric}")
        if score > best_score:
            best_score = score
            best_threshold = float(threshold)
    return ThresholdResult(best_threshold, metric, True)
