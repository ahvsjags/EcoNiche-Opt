from pathlib import Path

import numpy as np
import pandas as pd

from econiche_opt.analysis.precomputed_scores import import_precomputed_scores
from econiche_opt.data.controlled_access import controlled_access_decision, is_controlled_access
from econiche_opt.model.calibration import fit_training_only_calibrator
from econiche_opt.model.io import load_model, save_model
from econiche_opt.model.thresholds import select_threshold_training_only
from econiche_opt.preprocess.gene_coverage import coverage_fraction, gene_coverage_report


def test_model_serialization_roundtrip(tmp_path: Path):
    model = {"coef": np.array([1.0, 2.0]), "metadata": {"training_only": True}}
    path = save_model(model, tmp_path / "model.joblib")
    loaded = load_model(path)
    assert np.allclose(loaded["coef"], model["coef"])
    assert loaded["metadata"]["training_only"]


def test_threshold_and_calibration_mark_training_only():
    y = np.array([0, 0, 1, 1])
    p = np.array([0.1, 0.4, 0.6, 0.9])
    threshold = select_threshold_training_only(y, p)
    calibrator = fit_training_only_calibrator(y, p)
    assert threshold.training_only
    assert calibrator.training_only
    assert calibrator.predict([0.2, 0.8]).shape == (2,)


def test_precomputed_scores_import_validates_columns(tmp_path: Path):
    path = tmp_path / "scores.tsv"
    path.write_text("sample_id\tmodel\tscore\nS1\tTIDE\t0.2\n", encoding="utf-8")
    frame = import_precomputed_scores(path)
    assert frame.loc[0, "source"] == "scores.tsv"


def test_gene_coverage_and_controlled_access():
    report = gene_coverage_report(["A", "B"], ["A", "C"])
    assert report["available"].tolist() == [True, False]
    assert coverage_fraction(["A", "B"], ["A", "C"]) == 0.5
    assert is_controlled_access("dbGaP_or_controlled_verify")
    decision = controlled_access_decision("Liu_DFCI", "dbGaP_or_controlled_verify")
    assert decision.status == "ACCESS_RESTRICTED"
