from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from econiche.normalize import rank_gaussian_normalize
from econiche_opt.model.calibration import TrainingOnlyCalibrator, fit_training_only_calibrator
from econiche_opt.model.ecology_optimizer import (
    HeuristicEcologyConfig,
    HeuristicEcologyResult,
    build_ecology_features_from_module,
    optimize_ecology_module,
)
from econiche_opt.model.endpoint_modules import (
    WORD_INTERACTION_EDGES,
    WORD_STATE_GENE_SETS,
    _estimate_word_gene_directions,
    _rank_expression_by_cohort,
    build_word_ecology_features,
)
from econiche_opt.model.io import load_model, save_model
from econiche_opt.model.thresholds import ThresholdResult, select_threshold_training_only


ModelMode = Literal["word_full_graph", "heuristic_ecology"]


@dataclass(frozen=True)
class EcoNichePackageMetadata:
    """Audit metadata stored with fitted public package models."""

    mode: str
    response_positive_label: int = 1
    label_semantics: str = "1=response/responder, 0=non-response/non-responder"
    training_only_selection: bool = True
    include_interactions: bool = True
    signed: bool = True
    calibration: str | None = None
    threshold_metric: str = "balanced_accuracy"


class EcoNicheOptClassifier(BaseEstimator, ClassifierMixin):
    """Sklearn-style public interface for the EcoNiche-Opt model framework.

    Parameters
    ----------
    mode:
        ``"word_full_graph"`` uses the locked six-state Word-spec model:
        signed rank module scores plus ecological interaction edges.
        ``"heuristic_ecology"`` runs the training-only heuristic ecology
        optimizer to select state genes and interaction edges before fitting
        the final classifier.
    include_interactions:
        Include ecological interaction edge features.
    signed:
        Estimate gene directions on the training data and use signed rank
        scores. Set to ``False`` for ablation-style unsigned features.
    calibration:
        ``None`` for raw logistic probabilities or ``"isotonic"`` for an
        in-sample training-only calibrator. External validation should fit any
        calibration only inside the training set.
    optimizer_config:
        Optional ``HeuristicEcologyConfig`` used only when
        ``mode="heuristic_ecology"``.

    Notes
    -----
    This class never downloads data and never uses holdout labels during
    feature selection, thresholding, calibration, or model fitting. Users must
    pass expression matrices with samples as rows and genes as columns.
    """

    def __init__(
        self,
        mode: ModelMode = "word_full_graph",
        include_interactions: bool = True,
        signed: bool = True,
        calibration: Literal["isotonic"] | None = None,
        threshold_metric: str = "balanced_accuracy",
        class_weight: str | dict[int, float] | None = "balanced",
        random_state: int = 42,
        optimizer_config: HeuristicEcologyConfig | None = None,
    ) -> None:
        self.mode = mode
        self.include_interactions = include_interactions
        self.signed = signed
        self.calibration = calibration
        self.threshold_metric = threshold_metric
        self.class_weight = class_weight
        self.random_state = random_state
        self.optimizer_config = optimizer_config

    @staticmethod
    def _as_expression_frame(X: pd.DataFrame | np.ndarray, index: list[str] | None = None) -> pd.DataFrame:
        if isinstance(X, pd.DataFrame):
            frame = X.copy()
        else:
            frame = pd.DataFrame(X, index=index)
        if frame.index.has_duplicates:
            raise ValueError("Expression matrix sample index contains duplicates")
        frame.columns = [str(column) for column in frame.columns]
        return frame.apply(pd.to_numeric, errors="coerce")

    @staticmethod
    def _as_label_series(y: pd.Series | np.ndarray | list[int], index: pd.Index) -> pd.Series:
        labels = pd.Series(y, index=index) if not isinstance(y, pd.Series) else y.reindex(index)
        labels = pd.to_numeric(labels, errors="coerce")
        mask = labels.notna()
        if int(mask.sum()) != len(labels):
            raise ValueError("Labels must be complete for all training samples")
        labels = labels.astype(int)
        unique = set(labels.unique().tolist())
        if not unique.issubset({0, 1}):
            raise ValueError("EcoNicheOptClassifier expects binary labels encoded as 0/1")
        if labels.nunique() < 2:
            raise ValueError("Training labels must contain both classes")
        return labels

    @staticmethod
    def _default_metadata(X_by_cohort: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
        metadata: dict[str, pd.DataFrame] = {}
        for cohort, X in X_by_cohort.items():
            metadata[cohort] = pd.DataFrame(
                {
                    "sample_id": X.index.astype(str),
                    "patient_id": X.index.astype(str),
                    "cohort": cohort,
                },
                index=X.index,
            )
        return metadata

    def fit(self, X: pd.DataFrame | np.ndarray, y: pd.Series | np.ndarray | list[int]) -> "EcoNicheOptClassifier":
        """Fit on one training cohort.

        ``y`` must encode response as ``1`` and non-response as ``0``.
        For multicohort leakage-safe training, use ``fit_multicohort``.
        """

        frame = self._as_expression_frame(X)
        labels = self._as_label_series(y, frame.index)
        return self.fit_multicohort({"training": frame}, {"training": labels})

    def fit_multicohort(
        self,
        X_by_cohort: dict[str, pd.DataFrame],
        y_by_cohort: dict[str, pd.Series],
        metadata_by_cohort: dict[str, pd.DataFrame] | None = None,
        train_cohorts: list[str] | None = None,
    ) -> "EcoNicheOptClassifier":
        """Fit using multiple training cohorts.

        All feature direction estimation, module/edge selection, calibration,
        and threshold selection are restricted to ``train_cohorts``.
        """

        if self.mode not in {"word_full_graph", "heuristic_ecology"}:
            raise ValueError(f"Unsupported mode: {self.mode}")
        cohorts = list(train_cohorts or X_by_cohort.keys())
        if not cohorts:
            raise ValueError("At least one training cohort is required")

        X_clean = {cohort: self._as_expression_frame(X_by_cohort[cohort]) for cohort in cohorts}
        y_clean = {cohort: self._as_label_series(y_by_cohort[cohort], X_clean[cohort].index) for cohort in cohorts}
        metadata = metadata_by_cohort or self._default_metadata(X_clean)
        metadata = {cohort: metadata.get(cohort, pd.DataFrame(index=X_clean[cohort].index)).reindex(X_clean[cohort].index) for cohort in cohorts}

        if self.mode == "word_full_graph":
            ranked = _rank_expression_by_cohort(X_clean)
            self.gene_directions_ = (
                _estimate_word_gene_directions(ranked, y_clean, cohorts) if self.signed else {gene: 1 for genes in WORD_STATE_GENE_SETS.values() for gene in genes}
            )
            feature_by_cohort = {
                cohort: build_word_ecology_features(
                    X_clean[cohort],
                    gene_directions=self.gene_directions_,
                    signed=self.signed,
                    include_interactions=self.include_interactions,
                )[0]
                for cohort in cohorts
            }
            self.optimizer_result_ = None
        else:
            ranked = {cohort: rank_gaussian_normalize(X_clean[cohort].astype(float)) for cohort in cohorts}
            cfg = self.optimizer_config or HeuristicEcologyConfig(random_state=self.random_state)
            self.optimizer_result_ = optimize_ecology_module(
                ranked,
                y_clean,
                cohorts,
                WORD_STATE_GENE_SETS,
                WORD_INTERACTION_EDGES,
                metadata_by_cohort=metadata,
                signed=self.signed,
                include_interactions=self.include_interactions,
                cfg=cfg,
            )
            self.gene_directions_ = dict(self.optimizer_result_.gene_directions)
            feature_by_cohort = {
                cohort: build_ecology_features_from_module(
                    ranked[cohort].reindex(y_clean[cohort].index),
                    self.optimizer_result_.genes_by_state,
                    self.optimizer_result_.edges,
                    self.gene_directions_,
                    signed=self.signed,
                    include_interactions=self.include_interactions,
                )
                for cohort in cohorts
            }

        X_features = pd.concat([feature_by_cohort[cohort] for cohort in cohorts], axis=0)
        y_train = pd.concat([y_clean[cohort].reindex(feature_by_cohort[cohort].index) for cohort in cohorts], axis=0).astype(int)
        X_features = X_features.replace([np.inf, -np.inf], np.nan).fillna(0.0)
        self.feature_columns_ = list(X_features.columns)
        self.model_ = Pipeline(
            [
                ("scale", StandardScaler()),
                (
                    "logistic",
                    LogisticRegression(
                        C=0.75,
                        class_weight=self.class_weight,
                        solver="lbfgs",
                        max_iter=3000,
                        random_state=self.random_state,
                    ),
                ),
            ]
        )
        self.model_.fit(X_features.astype(float), y_train)
        train_prob = self.model_.predict_proba(X_features.astype(float))[:, 1]
        self.threshold_ = select_threshold_training_only(y_train, train_prob, metric=self.threshold_metric)
        self.calibrator_: TrainingOnlyCalibrator | None = None
        if self.calibration is not None:
            self.calibrator_ = fit_training_only_calibrator(y_train, train_prob, method=self.calibration)
            train_prob = self.calibrator_.predict(train_prob)

        self.training_metadata_ = EcoNichePackageMetadata(
            mode=self.mode,
            include_interactions=self.include_interactions,
            signed=self.signed,
            calibration=self.calibration,
            threshold_metric=self.threshold_metric,
        )
        self.training_summary_ = {
            "n_cohorts": len(cohorts),
            "cohorts": list(cohorts),
            "n_samples": int(len(y_train)),
            "n_response": int((y_train == 1).sum()),
            "n_nonresponse": int((y_train == 0).sum()),
            "threshold": float(self.threshold_.threshold),
            "feature_count": int(len(self.feature_columns_)),
        }
        return self

    def _check_is_fit(self) -> None:
        if not hasattr(self, "model_") or not hasattr(self, "feature_columns_"):
            raise RuntimeError("EcoNicheOptClassifier must be fit before prediction")

    def transform(self, X: pd.DataFrame | np.ndarray) -> pd.DataFrame:
        """Transform expression into EcoNiche-Opt model features."""

        self._check_is_fit()
        frame = self._as_expression_frame(X)
        if self.mode == "word_full_graph":
            features = build_word_ecology_features(
                frame,
                gene_directions=self.gene_directions_,
                signed=self.signed,
                include_interactions=self.include_interactions,
            )[0]
        else:
            ranked = rank_gaussian_normalize(frame.astype(float))
            result: HeuristicEcologyResult = self.optimizer_result_
            features = build_ecology_features_from_module(
                ranked,
                result.genes_by_state,
                result.edges,
                self.gene_directions_,
                signed=self.signed,
                include_interactions=self.include_interactions,
            )
        return features.reindex(columns=self.feature_columns_, fill_value=0.0).replace([np.inf, -np.inf], np.nan).fillna(0.0)

    def predict_response_proba(self, X: pd.DataFrame | np.ndarray) -> np.ndarray:
        """Return P(response) for each sample."""

        features = self.transform(X)
        prob = self.model_.predict_proba(features.astype(float))[:, 1]
        if self.calibrator_ is not None:
            prob = self.calibrator_.predict(prob)
        return np.clip(prob, 0.0, 1.0)

    def predict_proba(self, X: pd.DataFrame | np.ndarray) -> np.ndarray:
        """Return sklearn-style ``[P(non-response), P(response)]`` columns."""

        response = self.predict_response_proba(X)
        return np.column_stack([1.0 - response, response])

    def predict(self, X: pd.DataFrame | np.ndarray) -> np.ndarray:
        response = self.predict_response_proba(X)
        return (response >= float(self.threshold_.threshold)).astype(int)

    def score_samples(self, X: pd.DataFrame | np.ndarray) -> pd.DataFrame:
        """Return a tidy per-sample prediction table."""

        self._check_is_fit()
        frame = self._as_expression_frame(X)
        features = self.transform(frame)
        response = self.predict_response_proba(frame)
        decision = self.model_.decision_function(features.astype(float))
        return pd.DataFrame(
            {
                "sample_id": frame.index.astype(str),
                "response_probability": response,
                "nonresponse_probability": 1.0 - response,
                "predicted_response_label": (response >= float(self.threshold_.threshold)).astype(int),
                "decision_score": decision,
                "threshold": float(self.threshold_.threshold),
                "mode": self.mode,
            },
            index=frame.index,
        )

    def feature_coverage(self, X: pd.DataFrame | np.ndarray) -> pd.DataFrame:
        """Return feature-level gene/edge coverage for an expression matrix."""

        self._check_is_fit()
        frame = self._as_expression_frame(X)
        if self.mode == "word_full_graph":
            return build_word_ecology_features(
                frame,
                gene_directions=self.gene_directions_,
                signed=self.signed,
                include_interactions=self.include_interactions,
            )[1]
        result: HeuristicEcologyResult = self.optimizer_result_
        rows: list[dict[str, Any]] = []
        for state, genes in result.genes_by_state.items():
            available = [gene for gene in genes if gene in frame.columns]
            rows.append(
                {
                    "feature": state,
                    "feature_type": "optimized_state",
                    "state": state,
                    "n_genes_defined": len(genes),
                    "n_genes_available": len(available),
                    "genes_available": ",".join(available),
                }
            )
        for source, target, gene_a, gene_b, edge_class in result.edges:
            rows.append(
                {
                    "feature": f"interaction__{source}__{target}",
                    "feature_type": "optimized_interaction",
                    "source_state": source,
                    "target_state": target,
                    "edge_class": edge_class,
                    "gene_a": gene_a,
                    "gene_b": gene_b,
                    "edge_available": bool(gene_a in frame.columns and gene_b in frame.columns),
                }
            )
        return pd.DataFrame(rows)

    def module_table(self) -> pd.DataFrame:
        """Return genes, states, directions, and selection evidence."""

        self._check_is_fit()
        if self.mode == "heuristic_ecology" and self.optimizer_result_ is not None:
            return self.optimizer_result_.gene_table.copy()
        rows: list[dict[str, Any]] = []
        for state, genes in WORD_STATE_GENE_SETS.items():
            for gene in genes:
                rows.append(
                    {
                        "state": state,
                        "gene": gene,
                        "direction": int(self.gene_directions_.get(gene, 1)),
                        "selection_frequency": 1.0,
                        "source": "locked_word_full_graph",
                    }
                )
        return pd.DataFrame(rows)

    def edge_table(self) -> pd.DataFrame:
        """Return ecological interaction edges used by the model."""

        self._check_is_fit()
        edges = self.optimizer_result_.edges if self.mode == "heuristic_ecology" and self.optimizer_result_ is not None else WORD_INTERACTION_EDGES
        return pd.DataFrame(
            [
                {
                    "source_state": source,
                    "target_state": target,
                    "gene_a": gene_a,
                    "gene_b": gene_b,
                    "edge_class": edge_class,
                    "edge_id": f"{source}|{target}|{gene_a}|{gene_b}|{edge_class}",
                }
                for source, target, gene_a, gene_b, edge_class in edges
            ]
        )

    def package_metadata(self) -> dict[str, Any]:
        """Return reproducibility metadata for the fitted package model."""

        self._check_is_fit()
        return {
            "training_metadata": asdict(self.training_metadata_),
            "training_summary": dict(self.training_summary_),
            "threshold": asdict(self.threshold_) if isinstance(self.threshold_, ThresholdResult) else self.threshold_,
        }

    def save(self, path: str | Path) -> Path:
        """Serialize a fitted EcoNiche-Opt classifier with joblib."""

        self._check_is_fit()
        return save_model(self, path)

    @classmethod
    def load(cls, path: str | Path) -> "EcoNicheOptClassifier":
        """Load a classifier saved by ``save``."""

        loaded = load_model(path)
        if not isinstance(loaded, cls):
            raise TypeError(f"Expected {cls.__name__}, found {type(loaded).__name__}")
        return loaded
