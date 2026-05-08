from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from econiche_opt.api import EcoNicheOptClassifier
from econiche_opt.model.ecology_optimizer import HeuristicEcologyConfig


def read_expression_table(path: str | Path, transpose: bool = False) -> pd.DataFrame:
    """Read a TSV/CSV expression matrix with samples as rows and genes as columns."""

    path = Path(path)
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    frame = pd.read_csv(path, sep=sep, index_col=0)
    if transpose:
        frame = frame.T
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    return frame.apply(pd.to_numeric, errors="coerce")


def read_label_table(
    path: str | Path,
    sample_column: str = "sample_id",
    label_column: str = "response_label",
    cohort_column: str = "cohort",
) -> pd.DataFrame:
    """Read sample labels for package training.

    Labels must be encoded as 1=response/responder and 0=non-response.
    """

    path = Path(path)
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    labels = pd.read_csv(path, sep=sep)
    if sample_column not in labels.columns:
        first = labels.columns[0]
        labels = labels.rename(columns={first: sample_column})
    if label_column not in labels.columns:
        fallback = next((name for name in ["label", "y", "response", "response_label"] if name in labels.columns), None)
        if fallback is None:
            raise ValueError(f"Could not find label column '{label_column}' or a fallback label/y/response column")
        labels = labels.rename(columns={fallback: label_column})
    labels[sample_column] = labels[sample_column].astype(str)
    if cohort_column not in labels.columns:
        labels[cohort_column] = "training"
    labels[cohort_column] = labels[cohort_column].astype(str)
    labels[label_column] = pd.to_numeric(labels[label_column], errors="coerce")
    if labels[label_column].isna().any():
        raise ValueError("Label table contains missing or non-numeric labels")
    labels = labels.set_index(sample_column, drop=False)
    return labels


def split_training_tables(
    expression: pd.DataFrame,
    labels: pd.DataFrame,
    label_column: str = "response_label",
    cohort_column: str = "cohort",
) -> tuple[dict[str, pd.DataFrame], dict[str, pd.Series], dict[str, pd.DataFrame]]:
    common = expression.index.intersection(labels.index)
    if len(common) < 4:
        raise ValueError("Need at least four labeled samples shared by expression and labels")
    expression = expression.loc[common].copy()
    labels = labels.loc[common].copy()
    X_by_cohort: dict[str, pd.DataFrame] = {}
    y_by_cohort: dict[str, pd.Series] = {}
    metadata_by_cohort: dict[str, pd.DataFrame] = {}
    for cohort, cohort_labels in labels.groupby(cohort_column, sort=True):
        samples = cohort_labels.index.astype(str)
        X = expression.loc[samples].copy()
        y = cohort_labels[label_column].astype(int)
        if y.nunique() < 2:
            continue
        X_by_cohort[str(cohort)] = X
        y_by_cohort[str(cohort)] = pd.Series(y.values, index=X.index, name="response_label")
        metadata_by_cohort[str(cohort)] = cohort_labels.reindex(X.index).copy()
    if not X_by_cohort:
        raise ValueError("No cohort contains both response classes after label alignment")
    return X_by_cohort, y_by_cohort, metadata_by_cohort


def write_training_artifacts(model: EcoNicheOptClassifier, out_dir: str | Path) -> dict[str, str]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {
        "module_table": out / "econiche_module_table.tsv",
        "edge_table": out / "econiche_edge_table.tsv",
        "metadata": out / "econiche_model_metadata.json",
    }
    model.module_table().to_csv(paths["module_table"], sep="\t", index=False)
    model.edge_table().to_csv(paths["edge_table"], sep="\t", index=False)
    paths["metadata"].write_text(json.dumps(model.package_metadata(), indent=2, sort_keys=True), encoding="utf-8")
    return {key: str(path) for key, path in paths.items()}


def fit_package_model(
    expression_path: str | Path,
    labels_path: str | Path,
    model_path: str | Path,
    out_dir: str | Path,
    mode: str = "word_full_graph",
    transpose: bool = False,
    label_column: str = "response_label",
    cohort_column: str = "cohort",
    calibration: str | None = None,
    random_state: int = 42,
    optimizer_kwargs: dict[str, Any] | None = None,
) -> EcoNicheOptClassifier:
    expression = read_expression_table(expression_path, transpose=transpose)
    labels = read_label_table(labels_path, label_column=label_column, cohort_column=cohort_column)
    X_by_cohort, y_by_cohort, metadata_by_cohort = split_training_tables(
        expression,
        labels,
        label_column=label_column,
        cohort_column=cohort_column,
    )
    cfg = HeuristicEcologyConfig(**(optimizer_kwargs or {})) if mode == "heuristic_ecology" else None
    model = EcoNicheOptClassifier(
        mode=mode,  # type: ignore[arg-type]
        calibration=calibration,  # type: ignore[arg-type]
        random_state=random_state,
        optimizer_config=cfg,
    )
    model.fit_multicohort(X_by_cohort, y_by_cohort, metadata_by_cohort)
    model.save(model_path)
    write_training_artifacts(model, out_dir)
    return model


def score_package_model(
    model_path: str | Path,
    expression_path: str | Path,
    out_path: str | Path,
    transpose: bool = False,
    coverage_path: str | Path | None = None,
) -> pd.DataFrame:
    model = EcoNicheOptClassifier.load(model_path)
    expression = read_expression_table(expression_path, transpose=transpose)
    scores = model.score_samples(expression)
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    scores.to_csv(out, sep="\t", index=False)
    if coverage_path is not None:
        coverage = model.feature_coverage(expression)
        cov = Path(coverage_path)
        cov.parent.mkdir(parents=True, exist_ok=True)
        coverage.to_csv(cov, sep="\t", index=False)
    return scores
