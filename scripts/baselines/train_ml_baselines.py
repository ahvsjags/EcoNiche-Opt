from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, balanced_accuracy_score, brier_score_loss, f1_score, matthews_corrcoef, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.io import load_processed_bulk
from econiche.normalize import intersect_gene_space, rank_gaussian_normalize


def _response_labels(nonresponse_labels: pd.Series) -> pd.Series:
    return 1 - pd.Series(nonresponse_labels, index=nonresponse_labels.index).astype(int)


def _build_models(random_state: int = 42) -> dict[str, object]:
    models: dict[str, object] = {
        "ML_Logistic_L2": Pipeline(
            [
                ("scale", StandardScaler()),
                (
                    "model",
                    LogisticRegression(
                        C=0.03,
                        penalty="l2",
                        solver="liblinear",
                        class_weight="balanced",
                        max_iter=5000,
                        random_state=random_state,
                    ),
                ),
            ]
        ),
        "ML_Logistic_L1": Pipeline(
            [
                ("scale", StandardScaler()),
                (
                    "model",
                    LogisticRegression(
                        C=0.03,
                        penalty="l1",
                        solver="liblinear",
                        class_weight="balanced",
                        max_iter=5000,
                        random_state=random_state,
                    ),
                ),
            ]
        ),
        "ML_RandomForest": RandomForestClassifier(
            n_estimators=300,
            max_features="sqrt",
            min_samples_leaf=3,
            class_weight="balanced",
            random_state=random_state,
            n_jobs=-1,
        ),
        "ML_ExtraTrees": ExtraTreesClassifier(
            n_estimators=300,
            max_features="sqrt",
            min_samples_leaf=3,
            class_weight="balanced",
            random_state=random_state,
            n_jobs=-1,
        ),
    }
    try:
        from xgboost import XGBClassifier

        models["ML_XGBoost"] = XGBClassifier(
            n_estimators=120,
            max_depth=2,
            learning_rate=0.03,
            subsample=0.8,
            colsample_bytree=0.8,
            reg_lambda=5.0,
            objective="binary:logistic",
            eval_metric="logloss",
            random_state=random_state,
            n_jobs=1,
        )
    except Exception:
        pass
    try:
        from lightgbm import LGBMClassifier

        models["ML_LightGBM"] = LGBMClassifier(
            n_estimators=120,
            max_depth=2,
            learning_rate=0.03,
            subsample=0.8,
            colsample_bytree=0.8,
            reg_lambda=5.0,
            class_weight="balanced",
            random_state=random_state,
            n_jobs=1,
            verbose=-1,
        )
    except Exception:
        pass
    return models


def _metrics(y_true: pd.Series, prob: np.ndarray) -> dict[str, float]:
    pred = prob >= 0.5
    return {
        "AUROC": float(roc_auc_score(y_true, prob)) if y_true.nunique() == 2 else np.nan,
        "AUPRC": float(average_precision_score(y_true, prob)) if y_true.nunique() == 2 else np.nan,
        "balanced_accuracy": float(balanced_accuracy_score(y_true, pred)),
        "MCC": float(matthews_corrcoef(y_true, pred)),
        "F1": float(f1_score(y_true, pred, zero_division=0)),
        "Brier": float(brier_score_loss(y_true, prob)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--out-dir", default="results/real")
    parser.add_argument("--include-demo", action="store_true")
    parser.add_argument("--random-state", type=int, default=42)
    args = parser.parse_args()

    X_by_cohort, y_by_cohort, metadata_by_cohort = load_processed_bulk(ROOT / args.processed_dir)
    cohorts = sorted(X_by_cohort)
    if not args.include_demo:
        cohorts = [cohort for cohort in cohorts if not cohort.startswith("demo_cohort_")]
    X_by_cohort = {cohort: X_by_cohort[cohort] for cohort in cohorts}
    y_by_cohort = {cohort: y_by_cohort[cohort] for cohort in cohorts}
    metadata_by_cohort = {cohort: metadata_by_cohort[cohort] for cohort in cohorts}
    if len(X_by_cohort) < 2:
        raise SystemExit("Need at least two cohorts for ML baseline LODO.")

    common = intersect_gene_space(X_by_cohort)
    X_norm = {cohort: rank_gaussian_normalize(frame) for cohort, frame in common.items()}
    models = _build_models(args.random_state)
    prediction_rows = []
    metric_rows = []
    for holdout in sorted(X_norm):
        train_cohorts = [cohort for cohort in sorted(X_norm) if cohort != holdout]
        X_train = pd.concat([X_norm[cohort] for cohort in train_cohorts], axis=0)
        y_train = pd.concat([_response_labels(y_by_cohort[cohort]).reindex(X_norm[cohort].index) for cohort in train_cohorts])
        X_test = X_norm[holdout]
        y_test = _response_labels(y_by_cohort[holdout]).reindex(X_test.index)
        if y_train.nunique() < 2 or y_test.nunique() < 2:
            continue
        for model_name, model in models.items():
            fitted = model.fit(X_train, y_train)
            prob = fitted.predict_proba(X_test)[:, 1]
            metrics = _metrics(y_test, prob)
            metric_rows.append(
                {
                    **metrics,
                    "cohort": holdout,
                    "model_name": model_name,
                    "endpoint": "response_positive_primary_recist",
                    "n_samples": len(y_test),
                    "n_features": X_train.shape[1],
                }
            )
            metadata = metadata_by_cohort[holdout].reindex(X_test.index)
            for sample_id, probability, response_label in zip(X_test.index, prob, y_test):
                meta = metadata.loc[sample_id] if sample_id in metadata.index else pd.Series(dtype=object)
                prediction_rows.append(
                    {
                        "sample_id": meta.get("sample_id", sample_id),
                        "patient_id": meta.get("patient_id", pd.NA),
                        "cohort": holdout,
                        "true_label_nonresponse": int(1 - response_label),
                        "true_response_label": int(response_label),
                        "response_probability": float(probability),
                        "pred_response_label": int(probability >= 0.5),
                        "model_name": model_name,
                        "fold": holdout,
                    }
                )
    out_dir = ROOT / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(prediction_rows).to_csv(out_dir / "ml_baseline_predictions.tsv", sep="\t", index=False)
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(out_dir / "ml_baseline_metrics.tsv", sep="\t", index=False)
    summary = metrics.groupby("model_name", as_index=False)[["AUROC", "AUPRC", "balanced_accuracy", "MCC", "F1", "Brier"]].mean()
    summary.to_csv(out_dir / "ml_baseline_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))
    print(f"Wrote ML baseline outputs to {out_dir}")


if __name__ == "__main__":
    main()
