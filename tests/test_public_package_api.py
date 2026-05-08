from __future__ import annotations

import numpy as np
import pandas as pd

from econiche_opt import EcoNicheOptClassifier, load_demo_multicohort
from econiche_opt.cli import main
from econiche_opt.model.ecology_optimizer import HeuristicEcologyConfig


def test_public_classifier_fixed_word_graph_roundtrip(tmp_path):
    demo = load_demo_multicohort(n_cohorts=3, n_per_cohort=12, random_state=11)
    model = EcoNicheOptClassifier(mode="word_full_graph", calibration=None, random_state=11)
    model.fit_multicohort(demo["X_by_cohort"], demo["y_response_by_cohort"], demo["metadata_by_cohort"])

    cohort = sorted(demo["X_by_cohort"])[0]
    scores = model.score_samples(demo["X_by_cohort"][cohort])
    assert {"sample_id", "response_probability", "predicted_response_label", "threshold"}.issubset(scores.columns)
    assert np.isfinite(scores["response_probability"]).all()
    assert not model.module_table().empty
    assert not model.edge_table().empty

    saved = model.save(tmp_path / "econiche_model.joblib")
    loaded = EcoNicheOptClassifier.load(saved)
    assert np.allclose(
        model.predict_response_proba(demo["X_by_cohort"][cohort]),
        loaded.predict_response_proba(demo["X_by_cohort"][cohort]),
    )


def test_public_classifier_heuristic_optimizer_small_demo():
    demo = load_demo_multicohort(n_cohorts=2, n_per_cohort=10, random_state=13)
    cfg = HeuristicEcologyConfig(population_size=4, generations=2, use_gpu=False, random_state=13)
    model = EcoNicheOptClassifier(mode="heuristic_ecology", optimizer_config=cfg, random_state=13)
    model.fit_multicohort(demo["X_by_cohort"], demo["y_response_by_cohort"], demo["metadata_by_cohort"])

    cohort = sorted(demo["X_by_cohort"])[0]
    prob = model.predict_response_proba(demo["X_by_cohort"][cohort])
    assert prob.shape == (10,)
    assert model.package_metadata()["training_summary"]["n_cohorts"] == 2
    assert "optimized" in set(model.feature_coverage(demo["X_by_cohort"][cohort])["feature_type"].astype(str).str.split("_").str[0])


def test_package_cli_fit_and_score(tmp_path):
    demo = load_demo_multicohort(n_cohorts=2, n_per_cohort=10, random_state=17)
    expression = []
    labels = []
    for cohort, X in demo["X_by_cohort"].items():
        expression.append(X)
        y = demo["y_response_by_cohort"][cohort]
        labels.append(
            demo["metadata_by_cohort"][cohort][["sample_id", "cohort"]].assign(response_label=y.reindex(X.index).values)
        )
    expression_path = tmp_path / "expression.tsv"
    labels_path = tmp_path / "labels.tsv"
    model_path = tmp_path / "model.joblib"
    score_path = tmp_path / "scores.tsv"
    coverage_path = tmp_path / "coverage.tsv"
    out_dir = tmp_path / "artifacts"
    expression_frame = pd.concat(expression, axis=0)
    label_frame = pd.concat(labels, axis=0)
    expression_frame.to_csv(expression_path, sep="\t")
    label_frame.to_csv(labels_path, sep="\t", index=False)

    assert (
        main(
            [
                "fit-package-model",
                "--expression",
                str(expression_path),
                "--labels",
                str(labels_path),
                "--model-out",
                str(model_path),
                "--out-dir",
                str(out_dir),
            ]
        )
        == 0
    )
    assert model_path.exists()
    assert (out_dir / "econiche_module_table.tsv").exists()

    assert (
        main(
            [
                "score-package-model",
                "--model",
                str(model_path),
                "--expression",
                str(expression_path),
                "--out",
                str(score_path),
                "--coverage-out",
                str(coverage_path),
            ]
        )
        == 0
    )
    assert score_path.exists()
    assert coverage_path.exists()
