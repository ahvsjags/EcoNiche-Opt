import pandas as pd

from econiche.qc import check_patient_leakage, split_by_timepoint_priority
from econiche.statistics import make_lodo_folds


def test_check_patient_leakage_detects_cross_split_patient():
    folds = {
        "fold_a": {
            "train": pd.DataFrame({"patient_id": ["p1", "p2"]}),
            "test": pd.DataFrame({"patient_id": ["p2", "p3"]}),
        }
    }

    report = check_patient_leakage(folds)

    assert not report.loc[0, "ok"]
    assert report.loc[0, "overlap_count"] == 1


def test_lodo_folds_have_no_holdout_cohort_in_train_and_no_patient_overlap():
    meta = pd.DataFrame(
        {
            "cohort": ["c1", "c1", "c2", "c2", "c3", "c3"],
            "patient_id": ["p1", "p2", "p3", "p4", "p5", "p6"],
            "label": [0, 1, 0, 1, 0, 1],
        }
    )

    folds = make_lodo_folds(meta)
    report = check_patient_leakage(folds)

    assert len(folds) == 3
    assert report["ok"].all()
    for holdout, fold in folds.items():
        assert holdout not in set(fold["train"]["cohort"])


def test_split_by_timepoint_priority_keeps_primary_and_separates_secondary():
    meta = pd.DataFrame(
        {
            "sample_id": ["a", "b", "c"],
            "patient_id": ["p1", "p1", "p2"],
            "timepoint": ["pretreatment", "on_treatment", "progression"],
        }
    )

    primary, secondary, progression = split_by_timepoint_priority(meta)

    assert primary["sample_id"].tolist() == ["a"]
    assert secondary["sample_id"].tolist() == ["b"]
    assert progression["sample_id"].tolist() == ["c"]
