from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

from econiche.metrics import compute_binary_metrics
from econiche.module import EcoNicheConfig, EcoNicheModule, EcoNicheResult
from econiche.normalize import intersect_gene_space, rank_gaussian_normalize
from econiche.optim import make_optimizer_history, objective_for_metrics, select_module_from_priors
from econiche.priors import make_default_cell_state_priors
from econiche.scoring import compute_state_scores, estimate_gene_directions


class EcoNicheOpt:
    def __init__(
        self,
        config: EcoNicheConfig | None = None,
        priors: pd.DataFrame | None = None,
        pathways=None,
        network_edges=None,
        lr_edges=None,
    ):
        self.config = config or EcoNicheConfig()
        self.priors = priors
        self.pathways = pathways
        self.network_edges = network_edges
        self.lr_edges = lr_edges
        self.best_module_: EcoNicheModule | None = None
        self.gene_directions_: dict[str, int] = {}
        self.model_: LogisticRegression | None = None
        self.common_genes_: list[str] = []

    def _prepare(self, X_by_cohort: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
        common = intersect_gene_space(X_by_cohort)
        self.common_genes_ = list(next(iter(common.values())).columns) if common else []
        return {cohort: rank_gaussian_normalize(X) for cohort, X in common.items()}

    def _fit_one(
        self,
        X_by_cohort: dict[str, pd.DataFrame],
        y_by_cohort: dict[str, pd.Series],
    ) -> tuple[EcoNicheModule, dict[str, int], LogisticRegression, dict[str, float]]:
        X_train = pd.concat(X_by_cohort.values(), axis=0)
        y_train = pd.concat([pd.Series(y, index=X_by_cohort[cohort].index) for cohort, y in y_by_cohort.items()], axis=0)
        gene_directions = estimate_gene_directions(X_train, y_train)
        priors = self.priors if self.priors is not None else make_default_cell_state_priors(list(X_train.columns))
        module, selection_frequency = select_module_from_priors(X_train, y_train, priors, self.config)
        state_scores = compute_state_scores(X_train, module, gene_directions, normalize=False)
        clf = LogisticRegression(C=1e6, solver="lbfgs", max_iter=1000, random_state=self.config.random_state)
        clf.fit(state_scores, y_train.astype(int))
        return module, gene_directions, clf, selection_frequency

    def _predict_with(
        self,
        X: pd.DataFrame,
        metadata: pd.DataFrame,
        module: EcoNicheModule,
        gene_directions: dict[str, int],
        clf: LogisticRegression,
        y_true: pd.Series | None = None,
    ) -> pd.DataFrame:
        X_norm = rank_gaussian_normalize(X.loc[:, [gene for gene in self.common_genes_ if gene in X.columns]])
        state_scores = compute_state_scores(X_norm, module, gene_directions, normalize=False)
        prob = clf.predict_proba(state_scores)[:, 1]
        score = clf.decision_function(state_scores)
        meta = metadata.reindex(X.index).copy() if metadata is not None else pd.DataFrame(index=X.index)
        out = pd.DataFrame(
            {
                "sample_id": meta.get("sample_id", pd.Series(X.index, index=X.index)).values,
                "patient_id": meta.get("patient_id", pd.Series(pd.NA, index=X.index)).values,
                "cohort": meta.get("cohort", pd.Series(pd.NA, index=X.index)).values,
                "true_label": y_true.reindex(X.index).values if y_true is not None else pd.NA,
                "pred_prob": prob,
                "pred_label": (prob >= 0.5).astype(int),
                "EcoNicheScore": score,
            },
            index=X.index,
        )
        for state in state_scores.columns:
            out[f"{state}_subscore"] = state_scores[state].values
        return out

    def fit(
        self,
        X_by_cohort: dict[str, pd.DataFrame],
        y_by_cohort: dict[str, pd.Series],
        metadata_by_cohort: dict[str, pd.DataFrame],
    ) -> EcoNicheResult:
        X_prepared = self._prepare(X_by_cohort)
        priors = self.priors if self.priors is not None else make_default_cell_state_priors(self.common_genes_)
        self.priors = priors

        lodo_predictions = []
        metric_rows = []
        for holdout in sorted(X_prepared):
            train_cohorts = [cohort for cohort in X_prepared if cohort != holdout]
            if len(train_cohorts) < 1:
                continue
            X_train = {cohort: X_prepared[cohort] for cohort in train_cohorts}
            y_train = {cohort: y_by_cohort[cohort] for cohort in train_cohorts}
            module, gene_directions, clf, _ = self._fit_one(X_train, y_train)
            pred = self._predict_with(
                X_prepared[holdout],
                metadata_by_cohort[holdout],
                module,
                gene_directions,
                clf,
                pd.Series(y_by_cohort[holdout], index=X_prepared[holdout].index),
            )
            pred["fold"] = holdout
            lodo_predictions.append(pred)
            metrics = compute_binary_metrics(pred["true_label"], pred["pred_prob"])
            metrics.update(
                {
                    "cohort": holdout,
                    "model_name": "EcoNiche-Opt",
                    "n_samples": len(pred),
                    "n_responders": int((pred["true_label"] == 0).sum()),
                    "n_nonresponders": int((pred["true_label"] == 1).sum()),
                }
            )
            metric_rows.append(metrics)

        self.best_module_, self.gene_directions_, self.model_, selection_frequency = self._fit_one(X_prepared, y_by_cohort)
        all_predictions = pd.concat(lodo_predictions, axis=0).reset_index(drop=True) if lodo_predictions else pd.DataFrame()
        lodo_metrics = pd.DataFrame(metric_rows)
        objective_terms = objective_for_metrics(self.config, lodo_metrics, self.best_module_, priors)
        history = make_optimizer_history(self.config, objective_terms, self.best_module_)
        coefficients = pd.DataFrame(
            {
                "feature": list(self.model_.feature_names_in_),
                "coefficient": self.model_.coef_[0],
            }
        )
        coefficients = pd.concat(
            [pd.DataFrame([{"feature": "intercept", "coefficient": float(self.model_.intercept_[0])}]), coefficients],
            ignore_index=True,
        )
        return EcoNicheResult(
            best_module=self.best_module_,
            objective_terms=objective_terms,
            lodo_metrics=lodo_metrics,
            predictions=all_predictions,
            history=history,
            coefficients=coefficients,
            gene_directions=self.gene_directions_,
            selection_frequency=selection_frequency,
            priors=priors,
        )

    def score_samples(self, X: pd.DataFrame, metadata: pd.DataFrame | None = None) -> pd.DataFrame:
        if self.best_module_ is None or self.model_ is None:
            raise RuntimeError("EcoNicheOpt must be fit before scoring samples")
        meta = metadata if metadata is not None else pd.DataFrame(index=X.index)
        return self._predict_with(X, meta, self.best_module_, self.gene_directions_, self.model_)

    def predict_proba(self, X: pd.DataFrame, metadata: pd.DataFrame | None = None) -> np.ndarray:
        return self.score_samples(X, metadata)["pred_prob"].to_numpy()
