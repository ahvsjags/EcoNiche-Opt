import numpy as np

from econiche.metrics import calibration_bins, compute_binary_metrics, expected_calibration_error


def test_compute_binary_metrics_returns_required_fields_for_good_predictions():
    y_true = np.array([0, 0, 1, 1])
    prob = np.array([0.05, 0.2, 0.8, 0.95])

    metrics = compute_binary_metrics(y_true, prob, threshold=0.5)

    for name in [
        "AUROC",
        "AUPRC",
        "balanced_accuracy",
        "MCC",
        "F1",
        "sensitivity",
        "specificity",
        "ECE",
        "Brier",
    ]:
        assert name in metrics
    assert metrics["AUROC"] == 1.0
    assert metrics["balanced_accuracy"] == 1.0


def test_expected_calibration_error_and_bins_are_finite():
    y_true = np.array([0, 1, 0, 1, 1])
    prob = np.array([0.1, 0.6, 0.4, 0.8, 0.9])

    assert expected_calibration_error(y_true, prob, n_bins=3) >= 0
    bins = calibration_bins(y_true, prob, n_bins=3)
    assert {"bin", "n", "prob_mean", "event_rate"}.issubset(bins.columns)
