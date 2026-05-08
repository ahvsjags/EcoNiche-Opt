from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.isotonic import IsotonicRegression


@dataclass
class TrainingOnlyCalibrator:
    method: str
    model: IsotonicRegression
    training_only: bool = True

    def predict(self, prob):
        values = np.asarray(prob, dtype=float)
        return np.clip(self.model.predict(values), 0.0, 1.0)


def fit_training_only_calibrator(y_train, prob_train, method: str = "isotonic") -> TrainingOnlyCalibrator:
    if method != "isotonic":
        raise ValueError(f"Unsupported calibration method: {method}")
    y = np.asarray(y_train).astype(int)
    p = np.asarray(prob_train).astype(float)
    if len(y) != len(p):
        raise ValueError("y_train and prob_train must have the same length")
    model = IsotonicRegression(out_of_bounds="clip")
    model.fit(p, y)
    return TrainingOnlyCalibrator(method=method, model=model)
