import pandas as pd

from econiche.statistics import make_lodo_folds


def test_make_lodo_folds_require_labeled_train_and_holdout_samples():
    meta = pd.DataFrame(
        {
            "cohort": ["a", "a", "b", "b"],
            "patient_id": ["a1", "a2", "b1", "b2"],
            "label": [0, 1, 0, 1],
        }
    )

    folds = make_lodo_folds(meta)

    assert sorted(folds) == ["a", "b"]
    for fold in folds.values():
        assert fold["train"]["label"].notna().all()
        assert fold["test"]["label"].notna().all()
        assert set(fold["train"]["cohort"]).isdisjoint(set(fold["test"]["cohort"]))
