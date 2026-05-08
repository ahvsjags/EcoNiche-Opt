from pathlib import Path

import pandas as pd

from econiche_opt.preprocess.deconvolution import (
    evaluate_abundance_baselines,
    score_marker_abundance,
    score_processed_cohorts,
    summarize_abundance,
)


def test_marker_abundance_and_metrics(tmp_path: Path):
    expression = pd.DataFrame(
        {
            "CD3D": [1.0, 2.0, 8.0, 9.0],
            "CD3E": [1.0, 2.0, 7.0, 8.0],
            "CD8A": [1.0, 3.0, 8.0, 10.0],
            "CD68": [8.0, 7.0, 2.0, 1.0],
        },
        index=["S1", "S2", "S3", "S4"],
    )
    scores = score_marker_abundance(expression, cohort="C1", min_markers=2)
    assert {"sample_id", "cell_type", "abundance_score", "status"}.issubset(scores.columns)
    assert (scores.loc[scores["cell_type"] == "cd8_t_effector", "status"] == "scored_marker_z_baseline").all()

    summary = summarize_abundance(scores)
    assert "mean_abundance_score" in summary.columns

    metadata = pd.DataFrame({"sample_id": ["S1", "S2", "S3", "S4"], "label": [0, 0, 1, 1]})
    metrics = evaluate_abundance_baselines(scores, {"C1": metadata})
    cd8_all = metrics[(metrics["cohort"] == "all_public_processed") & (metrics["cell_type"] == "cd8_t_effector")]
    assert not cd8_all.empty
    assert cd8_all.iloc[0]["status"] == "PASS"


def test_processed_cohort_discovery_scores_real_not_demo(tmp_path: Path):
    input_dir = tmp_path / "bulk"
    input_dir.mkdir()
    pd.DataFrame({"CD3D": [1, 3], "CD3E": [2, 4]}, index=["S1", "S2"]).to_csv(
        input_dir / "GSE1.expr.tsv", sep="\t"
    )
    pd.DataFrame({"sample_id": ["S1", "S2"], "label": [0, 1]}).to_csv(
        input_dir / "GSE1.metadata.tsv", sep="\t", index=False
    )
    pd.DataFrame({"CD3D": [1, 3], "CD3E": [2, 4]}, index=["D1", "D2"]).to_csv(
        input_dir / "demo_cohort_1.expr.tsv", sep="\t"
    )
    scores, metadata = score_processed_cohorts(input_dir)
    assert set(scores["cohort"]) == {"GSE1"}
    assert "GSE1" in metadata
