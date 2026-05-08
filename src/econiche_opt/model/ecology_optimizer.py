from __future__ import annotations

from dataclasses import dataclass, field
import os
import random
from typing import Iterable

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from econiche.metrics import compute_binary_metrics


@dataclass(frozen=True)
class HeuristicEcologyConfig:
    min_genes_per_state: int = 3
    max_genes_per_state: int = 9
    candidate_pool_per_state: int = 28
    global_candidate_genes: int = 360
    population_size: int = 14
    generations: int = 8
    elite_fraction: float = 0.25
    mutation_rate: float = 0.35
    crossover_rate: float = 0.60
    robust_rho: float = 0.45
    max_auto_edges_per_state_pair: int = 2
    min_edge_abs_corr: float = 0.35
    n_jobs: int = 1
    random_state: int = 42
    use_gpu: bool = True
    use_bio_objective: bool = True


@dataclass
class HeuristicEcologyResult:
    genes_by_state: dict[str, list[str]]
    edges: list[tuple[str, str, str, str, str]]
    gene_directions: dict[str, int]
    feature_coefficients: dict[str, float]
    selection_frequency: dict[str, float]
    objective_terms: dict[str, float]
    history: pd.DataFrame
    gene_table: pd.DataFrame
    edge_table: pd.DataFrame
    backend: str
    backend_details: str = ""


def safe_zscore(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    sd = numeric.std(ddof=0)
    if pd.isna(sd) or sd <= 0:
        return pd.Series(0.0, index=values.index)
    return ((numeric - numeric.mean()) / sd).fillna(0.0)


def sigmoid(values: pd.Series | np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    return 1.0 / (1.0 + np.exp(-np.clip(arr, -30, 30)))


def _concat_ranked(
    ranked_expression_by_cohort: dict[str, pd.DataFrame],
    y_response: dict[str, pd.Series],
    train_cohorts: list[str],
) -> tuple[pd.DataFrame, pd.Series, pd.Series]:
    frames: list[pd.DataFrame] = []
    labels: list[pd.Series] = []
    cohort_labels: list[pd.Series] = []
    for cohort in train_cohorts:
        if cohort not in ranked_expression_by_cohort or cohort not in y_response:
            continue
        y = y_response[cohort].astype(float)
        X = ranked_expression_by_cohort[cohort].reindex(y.index)
        frames.append(X)
        labels.append(y)
        cohort_labels.append(pd.Series(cohort, index=y.index))
    if not frames:
        return pd.DataFrame(), pd.Series(dtype=float), pd.Series(dtype=str)
    X_train = pd.concat(frames, axis=0, join="outer")
    y_train = pd.concat(labels, axis=0).reindex(X_train.index)
    cohorts = pd.concat(cohort_labels, axis=0).reindex(X_train.index)
    return X_train, y_train, cohorts


def _gene_correlations(X_train: pd.DataFrame, y_train: pd.Series) -> pd.Series:
    if X_train.empty or y_train.nunique(dropna=True) < 2:
        return pd.Series(0.0, index=X_train.columns)
    corr = X_train.corrwith(y_train.astype(float), axis=0).replace([np.inf, -np.inf], np.nan)
    return corr.fillna(0.0)


def _direction_stability_by_gene(
    ranked_expression_by_cohort: dict[str, pd.DataFrame],
    y_response: dict[str, pd.Series],
    train_cohorts: list[str],
    genes: Iterable[str],
) -> dict[str, float]:
    result: dict[str, float] = {}
    for gene in genes:
        signs: list[int] = []
        for cohort in train_cohorts:
            X = ranked_expression_by_cohort.get(cohort)
            y = y_response.get(cohort)
            if X is None or y is None or gene not in X.columns:
                continue
            values = pd.to_numeric(X[gene].reindex(y.index), errors="coerce")
            valid = values.notna() & y.notna()
            if int(valid.sum()) < 4 or values.loc[valid].nunique() < 2 or y.loc[valid].nunique() < 2:
                continue
            corr = values.loc[valid].corr(y.loc[valid].astype(float))
            if pd.notna(corr):
                signs.append(1 if corr >= 0 else -1)
        result[gene] = float(max(signs.count(1), signs.count(-1)) / len(signs)) if signs else 0.5
    return result


def _candidate_universe(
    X_train: pd.DataFrame,
    correlations: pd.Series,
    state_gene_sets: dict[str, list[str]],
    seed_edges: list[tuple[str, str, str, str, str]],
    max_genes: int,
) -> list[str]:
    available = set(map(str, X_train.columns))
    seed_genes: set[str] = set()
    for genes in state_gene_sets.values():
        seed_genes.update(str(gene) for gene in genes)
    for _, _, gene_a, gene_b, _ in seed_edges:
        seed_genes.add(str(gene_a))
        seed_genes.add(str(gene_b))
    top = correlations.abs().sort_values(ascending=False).head(max_genes).index.astype(str).tolist()
    universe = [gene for gene in sorted(seed_genes) if gene in available]
    for gene in top:
        if gene in available and gene not in universe:
            universe.append(gene)
    return universe


def _build_candidate_pools(
    X_train: pd.DataFrame,
    correlations: pd.Series,
    stability: dict[str, float],
    state_gene_sets: dict[str, list[str]],
    seed_edges: list[tuple[str, str, str, str, str]],
    cfg: HeuristicEcologyConfig,
) -> tuple[dict[str, list[str]], dict[str, dict[str, float]], dict[str, float]]:
    universe = _candidate_universe(X_train, correlations, state_gene_sets, seed_edges, cfg.global_candidate_genes)
    seed_any = set().union(*(set(map(str, genes)) for genes in state_gene_sets.values())) if state_gene_sets else set()
    edge_genes = {str(gene) for _, _, gene_a, gene_b, _ in seed_edges for gene in (gene_a, gene_b)}
    abs_corr = correlations.abs().reindex(universe).fillna(0.0)
    max_corr = float(abs_corr.max()) if len(abs_corr) else 0.0
    corr_scaled = abs_corr / max_corr if max_corr > 0 else abs_corr
    pools: dict[str, list[str]] = {}
    state_scores: dict[str, dict[str, float]] = {}
    for state, seeds in state_gene_sets.items():
        seed_set = set(map(str, seeds))
        scores: dict[str, float] = {}
        for gene in universe:
            prior = 1.0 if gene in seed_set else 0.18 if gene in seed_any else 0.0
            if gene in edge_genes:
                prior += 0.10
            score = 1.20 * prior + 0.75 * float(corr_scaled.get(gene, 0.0)) + 0.25 * float(stability.get(gene, 0.5))
            scores[gene] = float(score)
        ranked = sorted(scores, key=lambda gene: scores[gene], reverse=True)
        pool = [gene for gene in ranked[: cfg.candidate_pool_per_state]]
        for gene in sorted(seed_set):
            if gene in X_train.columns and gene not in pool:
                pool.append(gene)
        if len(pool) < cfg.min_genes_per_state:
            for gene in ranked:
                if gene not in pool:
                    pool.append(gene)
                if len(pool) >= cfg.min_genes_per_state:
                    break
        pools[state] = pool
        state_scores[state] = scores
    global_scores = {gene: float(corr_scaled.get(gene, 0.0)) + 0.25 * float(stability.get(gene, 0.5)) for gene in universe}
    return pools, state_scores, global_scores


def _weighted_sample_without_replacement(
    rng: random.Random,
    genes: list[str],
    scores: dict[str, float],
    n: int,
) -> list[str]:
    selected: list[str] = []
    remaining = list(dict.fromkeys(genes))
    for _ in range(min(n, len(remaining))):
        weights = np.asarray([max(1e-6, float(scores.get(gene, 0.0))) for gene in remaining], dtype=float)
        weights = weights / weights.sum() if float(weights.sum()) > 0 else np.full(len(remaining), 1.0 / len(remaining))
        pick = rng.choices(remaining, weights=weights.tolist(), k=1)[0]
        selected.append(pick)
        remaining.remove(pick)
    return selected


def _make_individual(
    rng: random.Random,
    pools: dict[str, list[str]],
    state_scores: dict[str, dict[str, float]],
    cfg: HeuristicEcologyConfig,
) -> tuple[tuple[str, tuple[str, ...]], ...]:
    rows: list[tuple[str, tuple[str, ...]]] = []
    for state, pool in pools.items():
        if not pool:
            rows.append((state, tuple()))
            continue
        upper = min(cfg.max_genes_per_state, len(pool))
        lower = min(cfg.min_genes_per_state, upper)
        size = rng.randint(lower, upper) if upper >= lower else upper
        top = sorted(pool, key=lambda gene: state_scores[state].get(gene, 0.0), reverse=True)[: max(1, lower // 2)]
        remaining = [gene for gene in pool if gene not in top]
        sampled = top + _weighted_sample_without_replacement(rng, remaining, state_scores[state], max(0, size - len(top)))
        rows.append((state, tuple(sorted(dict.fromkeys(sampled)))))
    return tuple(rows)


def _individual_to_dict(individual: tuple[tuple[str, tuple[str, ...]], ...]) -> dict[str, list[str]]:
    return {state: list(genes) for state, genes in individual}


def _build_gene_neighbors(X_train: pd.DataFrame, genes: list[str], max_neighbors: int = 8) -> dict[str, list[str]]:
    available = [gene for gene in genes if gene in X_train.columns]
    if len(available) < 2:
        return {gene: [] for gene in genes}
    corr = X_train[available].corr().abs().replace([np.inf, -np.inf], np.nan).fillna(0.0)
    neighbors: dict[str, list[str]] = {}
    for gene in available:
        ranked = corr[gene].drop(labels=[gene], errors="ignore").sort_values(ascending=False)
        neighbors[gene] = ranked.head(max_neighbors).index.astype(str).tolist()
    return neighbors


def _mutate_individual(
    individual: tuple[tuple[str, tuple[str, ...]], ...],
    rng: random.Random,
    pools: dict[str, list[str]],
    state_scores: dict[str, dict[str, float]],
    neighbors: dict[str, list[str]],
    cfg: HeuristicEcologyConfig,
) -> tuple[tuple[str, tuple[str, ...]], ...]:
    state_genes = _individual_to_dict(individual)
    for state in list(state_genes):
        if rng.random() > cfg.mutation_rate:
            continue
        genes = list(state_genes[state])
        pool = [gene for gene in pools.get(state, []) if gene not in genes]
        action = rng.choice(["replace", "add", "drop", "neighbor"])
        if action == "drop" and len(genes) > cfg.min_genes_per_state:
            genes.remove(rng.choice(genes))
        elif action == "add" and pool and len(genes) < cfg.max_genes_per_state:
            genes.extend(_weighted_sample_without_replacement(rng, pool, state_scores[state], 1))
        elif action == "neighbor" and genes:
            source = rng.choice(genes)
            options = [gene for gene in neighbors.get(source, []) if gene in pools.get(state, []) and gene not in genes]
            if options:
                if len(genes) >= cfg.max_genes_per_state:
                    genes.remove(rng.choice(genes))
                genes.append(rng.choice(options))
            elif pool:
                genes[rng.randrange(len(genes))] = _weighted_sample_without_replacement(rng, pool, state_scores[state], 1)[0]
        elif pool and genes:
            genes[rng.randrange(len(genes))] = _weighted_sample_without_replacement(rng, pool, state_scores[state], 1)[0]
        state_genes[state] = sorted(dict.fromkeys(genes))
    return tuple((state, tuple(state_genes[state])) for state in sorted(state_genes))


def _crossover_individuals(
    left: tuple[tuple[str, tuple[str, ...]], ...],
    right: tuple[tuple[str, tuple[str, ...]], ...],
    rng: random.Random,
    state_scores: dict[str, dict[str, float]],
    cfg: HeuristicEcologyConfig,
) -> tuple[tuple[str, tuple[str, ...]], ...]:
    left_map = _individual_to_dict(left)
    right_map = _individual_to_dict(right)
    rows = []
    for state in sorted(left_map):
        merged = sorted(set(left_map[state]) | set(right_map.get(state, [])))
        upper = min(cfg.max_genes_per_state, len(merged))
        lower = min(cfg.min_genes_per_state, upper)
        target_size = rng.randint(lower, upper) if upper >= lower else upper
        selected = _weighted_sample_without_replacement(rng, merged, state_scores[state], target_size)
        rows.append((state, tuple(sorted(selected))))
    return tuple(rows)


def _interaction_edges_for_module(
    genes_by_state: dict[str, list[str]],
    seed_edges: list[tuple[str, str, str, str, str]],
    X_train: pd.DataFrame,
    cfg: HeuristicEcologyConfig,
    include_interactions: bool,
) -> list[tuple[str, str, str, str, str]]:
    if not include_interactions:
        return []
    edges: list[tuple[str, str, str, str, str]] = []
    state_lookup = {gene: state for state, genes in genes_by_state.items() for gene in genes}
    seen: set[tuple[str, str, str, str, str]] = set()
    for source, target, gene_a, gene_b, edge_class in seed_edges:
        gene_a = str(gene_a)
        gene_b = str(gene_b)
        if gene_a in state_lookup and gene_b in state_lookup:
            row = (source, target, gene_a, gene_b, edge_class)
            edges.append(row)
            seen.add(row)
    state_names = sorted(genes_by_state)
    for i, source in enumerate(state_names):
        for target in state_names[i:]:
            pairs: list[tuple[float, str, str]] = []
            for gene_a in genes_by_state[source]:
                for gene_b in genes_by_state[target]:
                    if source == target and gene_a >= gene_b:
                        continue
                    if gene_a not in X_train.columns or gene_b not in X_train.columns:
                        continue
                    corr = X_train[gene_a].corr(X_train[gene_b])
                    if pd.notna(corr) and abs(float(corr)) >= cfg.min_edge_abs_corr:
                        pairs.append((abs(float(corr)), gene_a, gene_b))
            pairs.sort(reverse=True)
            for _, gene_a, gene_b in pairs[: cfg.max_auto_edges_per_state_pair]:
                row = (source, target, gene_a, gene_b, "coexpression_network")
                if row not in seen:
                    edges.append(row)
                    seen.add(row)
    return edges


def build_ecology_features_from_module(
    ranked_X: pd.DataFrame,
    genes_by_state: dict[str, list[str]],
    edges: list[tuple[str, str, str, str, str]],
    gene_directions: dict[str, int],
    signed: bool = True,
    include_interactions: bool = True,
) -> pd.DataFrame:
    state_scores: dict[str, pd.Series] = {}
    for state, genes in genes_by_state.items():
        available = [gene for gene in genes if gene in ranked_X.columns]
        if not available:
            state_scores[state] = pd.Series(0.0, index=ranked_X.index)
            continue
        values = ranked_X[available].astype(float).copy()
        if signed:
            for gene in available:
                values[gene] = values[gene] * int(gene_directions.get(gene, 1))
        state_scores[state] = values.sum(axis=1) / np.sqrt(len(available))
    feature_frame = pd.DataFrame(state_scores, index=ranked_X.index).fillna(0.0)
    if not include_interactions or not edges:
        return feature_frame
    abundance = feature_frame.apply(safe_zscore, axis=0).apply(sigmoid, axis=0)
    grouped: dict[str, list[pd.Series]] = {}
    for source, target, gene_a, gene_b, _ in edges:
        feature = f"interaction__{source}__{target}"
        if gene_a not in ranked_X.columns or gene_b not in ranked_X.columns:
            continue
        source_abundance = abundance[source] if source in abundance.columns else pd.Series(1.0, index=ranked_X.index)
        target_abundance = abundance[target] if target in abundance.columns else pd.Series(1.0, index=ranked_X.index)
        gene_a_values = ranked_X[gene_a].astype(float) * (int(gene_directions.get(gene_a, 1)) if signed else 1)
        gene_b_values = ranked_X[gene_b].astype(float) * (int(gene_directions.get(gene_b, 1)) if signed else 1)
        grouped.setdefault(feature, []).append(gene_a_values * gene_b_values * source_abundance * target_abundance)
    for feature, values in grouped.items():
        feature_frame[feature] = safe_zscore(pd.concat(values, axis=1).mean(axis=1))
    return feature_frame.fillna(0.0)


def _fit_logistic_predict(X_train: pd.DataFrame, y_train: pd.Series, X_test: pd.DataFrame) -> np.ndarray:
    if y_train.nunique() < 2:
        return np.full(len(X_test), float(y_train.mean()) if len(y_train) else 0.5)
    X_train = X_train.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    X_test = X_test.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    model = Pipeline(
        [
            ("scale", StandardScaler()),
            (
                "logistic",
                LogisticRegression(
                    C=0.75,
                    class_weight="balanced",
                    solver="lbfgs",
                    max_iter=3000,
                    random_state=42,
                ),
            ),
        ]
    )
    model.fit(X_train.astype(float), y_train.astype(int))
    return model.predict_proba(X_test.astype(float))[:, 1]


def _try_xgboost_cuda_predict(X_train: pd.DataFrame, y_train: pd.Series, X_test: pd.DataFrame) -> tuple[np.ndarray | None, str]:
    try:
        from xgboost import XGBClassifier

        X_train = X_train.replace([np.inf, -np.inf], np.nan).fillna(0.0)
        X_test = X_test.replace([np.inf, -np.inf], np.nan).fillna(0.0)
        model = XGBClassifier(
            n_estimators=12,
            max_depth=2,
            learning_rate=0.08,
            subsample=0.9,
            colsample_bytree=0.9,
            reg_lambda=2.0,
            min_child_weight=2.0,
            tree_method="hist",
            device="cuda",
            eval_metric="logloss",
            verbosity=0,
            random_state=42,
        )
        model.fit(X_train.astype(float), y_train.astype(int))
        return model.predict_proba(X_test.astype(float))[:, 1], "xgboost_cuda"
    except Exception as exc:  # pragma: no cover - hardware dependent
        return None, f"xgboost_cuda_unavailable:{type(exc).__name__}"


def _fit_predict_backend(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    X_test: pd.DataFrame,
    use_gpu: bool,
) -> tuple[np.ndarray, str]:
    if use_gpu and len(X_train) >= 20 and X_train.shape[1] >= 4:
        pred, backend = _try_xgboost_cuda_predict(X_train, y_train, X_test)
        if pred is not None:
            return pred, backend
    return _fit_logistic_predict(X_train, y_train, X_test), "sklearn_logistic_cpu"


def _batch_dependence(feature_frame: pd.DataFrame, cohort_labels: pd.Series) -> float:
    labels = cohort_labels.reindex(feature_frame.index)
    ratios = []
    for column in feature_frame.columns:
        values = pd.to_numeric(feature_frame[column], errors="coerce")
        valid = values.notna() & labels.notna()
        if int(valid.sum()) < 4 or values.loc[valid].nunique() < 2 or labels.loc[valid].nunique() < 2:
            continue
        overall = float(values.loc[valid].var(ddof=0))
        if overall <= 0:
            continue
        ratios.append(float(values.loc[valid].groupby(labels.loc[valid]).mean().var(ddof=0) / (overall + 1e-9)))
    return float(np.clip(np.mean(ratios), 0.0, 1.0)) if ratios else 0.0


def _therapy_confounding(feature_frame: pd.DataFrame, metadata: pd.DataFrame | None) -> float:
    if metadata is None or metadata.empty:
        return 0.0
    treatment_columns = ["treatment", "therapy", "io_therapy", "antibody", "treatment_group"]
    available = [column for column in treatment_columns if column in metadata.columns]
    if not available:
        return 0.0
    treatment = metadata[available].astype("string").bfill(axis=1).iloc[:, 0].reindex(feature_frame.index)
    if treatment.nunique(dropna=True) < 2:
        return 0.0
    return _batch_dependence(feature_frame, treatment)


def _state_redundancy(feature_frame: pd.DataFrame, state_names: list[str]) -> float:
    columns = [column for column in state_names if column in feature_frame.columns]
    if len(columns) < 2:
        return 0.0
    corr = feature_frame[columns].corr().abs().to_numpy(dtype=float)
    upper = corr[np.triu_indices_from(corr, k=1)]
    upper = upper[np.isfinite(upper)]
    return float(np.mean(upper)) if len(upper) else 0.0


def _evaluate_individual(
    individual: tuple[tuple[str, tuple[str, ...]], ...],
    ranked_expression_by_cohort: dict[str, pd.DataFrame],
    y_response: dict[str, pd.Series],
    train_cohorts: list[str],
    X_train_ranked: pd.DataFrame,
    y_train: pd.Series,
    cohort_labels: pd.Series,
    metadata: pd.DataFrame | None,
    seed_edges: list[tuple[str, str, str, str, str]],
    state_gene_sets: dict[str, list[str]],
    gene_directions: dict[str, int],
    stability: dict[str, float],
    cfg: HeuristicEcologyConfig,
    signed: bool,
    include_interactions: bool,
) -> dict[str, object]:
    genes_by_state = _individual_to_dict(individual)
    edges = _interaction_edges_for_module(genes_by_state, seed_edges, X_train_ranked, cfg, include_interactions)
    feature_by_cohort = {
        cohort: build_ecology_features_from_module(
            ranked_expression_by_cohort[cohort].reindex(y_response[cohort].index),
            genes_by_state,
            edges,
            gene_directions,
            signed=signed,
            include_interactions=include_interactions,
        )
        for cohort in train_cohorts
    }
    metric_rows = []
    backends: list[str] = []
    if len(train_cohorts) >= 2:
        for inner_holdout in train_cohorts:
            inner_train = [cohort for cohort in train_cohorts if cohort != inner_holdout]
            X_inner_train = pd.concat([feature_by_cohort[cohort] for cohort in inner_train], axis=0)
            y_inner_train = pd.concat([y_response[cohort].reindex(feature_by_cohort[cohort].index) for cohort in inner_train], axis=0).astype(int)
            X_inner_test = feature_by_cohort[inner_holdout]
            y_inner_test = y_response[inner_holdout].reindex(X_inner_test.index).astype(int)
            if y_inner_train.nunique() < 2 or y_inner_test.nunique() < 2:
                continue
            columns = sorted(set(X_inner_train.columns) | set(X_inner_test.columns))
            X_inner_train = X_inner_train.reindex(columns=columns, fill_value=0.0)
            X_inner_test = X_inner_test.reindex(columns=columns, fill_value=0.0)
            prob, backend = _fit_predict_backend(X_inner_train, y_inner_train, X_inner_test, cfg.use_gpu)
            backends.append(backend)
            metric_rows.append(compute_binary_metrics(y_inner_test, prob))
    if not metric_rows:
        feature_frame = pd.concat(feature_by_cohort.values(), axis=0)
        pseudo_score = safe_zscore(feature_frame.mean(axis=1)).corr(y_train.reindex(feature_frame.index))
        metric_rows = [
            {
                "AUROC": 0.5 + 0.1 * (float(pseudo_score) if pd.notna(pseudo_score) else 0.0),
                "AUPRC": float(y_train.mean()) if len(y_train) else 0.5,
                "balanced_accuracy": 0.5,
                "ECE": 0.25,
            }
        ]
        backends.append("surrogate_single_train_pool")
    metrics = pd.DataFrame(metric_rows)
    feature_frame = pd.concat(feature_by_cohort.values(), axis=0)
    seed_by_state = {state: set(map(str, genes)) for state, genes in state_gene_sets.items()}
    selected_prior_hits = [
        1.0 if gene in seed_by_state.get(state, set()) else 0.0
        for state, genes in genes_by_state.items()
        for gene in genes
    ]
    cell_specificity = float(np.mean(selected_prior_hits)) if selected_prior_hits else 0.0
    edge_classes = [edge[-1] for edge in edges]
    lr = float(np.mean([edge_class == "ligand_receptor" for edge_class in edge_classes])) if edge_classes else 0.0
    pathway = float(np.mean([edge_class == "pathway" for edge_class in edge_classes])) if edge_classes else 0.0
    network = float(
        np.mean([edge_class in {"network", "regulatory", "checkpoint", "coexpression_network"} for edge_class in edge_classes])
    ) if edge_classes else 0.0
    selected_genes = [gene for genes in genes_by_state.values() for gene in genes]
    direction_stability = float(np.mean([stability.get(gene, 0.5) for gene in selected_genes])) if selected_genes else 0.5
    size = float(np.clip((len(selected_genes) - 18) / (72 - 18), 0.0, 1.0))
    batch = _batch_dependence(feature_frame, cohort_labels)
    redundancy = _state_redundancy(feature_frame, list(state_gene_sets))
    therapy = _therapy_confounding(feature_frame, metadata)
    bio_bonus = 0.08 * cell_specificity + 0.04 * pathway + 0.04 * network + 0.04 * lr + 0.08 * direction_stability
    bio_penalty = 0.04 * size + 0.08 * batch + 0.05 * redundancy + 0.05 * therapy
    performance = float(
        metrics["AUROC"].mean()
        - cfg.robust_rho * metrics["AUROC"].std(ddof=0)
        + 0.10 * metrics["AUPRC"].mean()
        + 0.05 * metrics["balanced_accuracy"].mean()
        - 0.15 * metrics["ECE"].mean()
    )
    score = performance + (bio_bonus - bio_penalty if cfg.use_bio_objective else 0.0)
    return {
        "individual": individual,
        "genes_by_state": genes_by_state,
        "edges": edges,
        "score": float(score),
        "performance_score": performance,
        "inner_mean_AUROC": float(metrics["AUROC"].mean()),
        "inner_sd_AUROC": float(metrics["AUROC"].std(ddof=0)),
        "inner_mean_AUPRC": float(metrics["AUPRC"].mean()),
        "inner_mean_balanced_accuracy": float(metrics["balanced_accuracy"].mean()),
        "inner_mean_ECE": float(metrics["ECE"].mean()),
        "bio_cell_specificity": cell_specificity,
        "bio_pathway": pathway,
        "bio_network": network,
        "bio_lr": lr,
        "bio_direction_stability": direction_stability,
        "penalty_size": size,
        "penalty_batch": batch,
        "penalty_redundancy": redundancy,
        "penalty_therapy_confounding": therapy,
        "bio_bonus": float(bio_bonus),
        "bio_penalty": float(bio_penalty),
        "bio_objective_delta": float(bio_bonus - bio_penalty),
        "n_genes": int(len(selected_genes)),
        "n_edges": int(len(edges)),
        "backend": ",".join(sorted(set(backends))),
    }


def optimize_ecology_module(
    ranked_expression_by_cohort: dict[str, pd.DataFrame],
    y_response: dict[str, pd.Series],
    train_cohorts: list[str],
    state_gene_sets: dict[str, list[str]],
    seed_edges: list[tuple[str, str, str, str, str]],
    metadata_by_cohort: dict[str, pd.DataFrame] | None = None,
    signed: bool = True,
    include_interactions: bool = True,
    cfg: HeuristicEcologyConfig | None = None,
) -> HeuristicEcologyResult:
    cfg = cfg or HeuristicEcologyConfig()
    rng = random.Random(cfg.random_state + abs(hash(tuple(train_cohorts))) % 100000)
    X_train, y_train, cohort_labels = _concat_ranked(ranked_expression_by_cohort, y_response, train_cohorts)
    if X_train.empty:
        raise ValueError("No training expression matrix available for ecological optimization")
    correlations = _gene_correlations(X_train, y_train)
    universe = _candidate_universe(X_train, correlations, state_gene_sets, seed_edges, cfg.global_candidate_genes)
    stability = _direction_stability_by_gene(ranked_expression_by_cohort, y_response, train_cohorts, universe)
    gene_directions = {gene: (1 if float(correlations.get(gene, 0.0)) >= 0 else -1) for gene in universe} if signed else {
        gene: 1 for gene in universe
    }
    pools, state_scores, global_scores = _build_candidate_pools(X_train, correlations, stability, state_gene_sets, seed_edges, cfg)
    neighbors = _build_gene_neighbors(X_train, universe)
    metadata = None
    if metadata_by_cohort:
        metadata = pd.concat(
            [metadata_by_cohort[cohort].reindex(y_response[cohort].index) for cohort in train_cohorts if cohort in metadata_by_cohort],
            axis=0,
        )
    population = [_make_individual(rng, pools, state_scores, cfg) for _ in range(max(2, cfg.population_size))]
    history_rows: list[dict[str, object]] = []
    selection_counts: dict[str, float] = {}
    best_record: dict[str, object] | None = None
    n_jobs = max(1, int(cfg.n_jobs))
    if n_jobs < 0:
        n_jobs = max(1, (os.cpu_count() or 2) - 1)
    for generation in range(1, max(1, cfg.generations) + 1):
        unique_population = list(dict.fromkeys(population))
        records = Parallel(n_jobs=n_jobs, prefer="threads")(
            delayed(_evaluate_individual)(
                individual,
                ranked_expression_by_cohort,
                y_response,
                train_cohorts,
                X_train,
                y_train,
                cohort_labels,
                metadata,
                seed_edges,
                state_gene_sets,
                gene_directions,
                stability,
                cfg,
                signed,
                include_interactions,
            )
            for individual in unique_population
        )
        records = sorted(records, key=lambda row: float(row["score"]), reverse=True)
        if best_record is None or float(records[0]["score"]) > float(best_record["score"]):
            best_record = records[0]
        elite_n = max(2, int(np.ceil(len(records) * cfg.elite_fraction)))
        elites = records[:elite_n]
        for record in elites:
            for genes in record["genes_by_state"].values():
                for gene in genes:
                    selection_counts[gene] = selection_counts.get(gene, 0.0) + 1.0
        history_rows.append(
            {
                "generation": generation,
                "best_score": float(records[0]["score"]),
                "mean_score": float(np.mean([float(row["score"]) for row in records])),
                "best_AUROC": float(records[0]["inner_mean_AUROC"]),
                "best_AUPRC": float(records[0]["inner_mean_AUPRC"]),
                "best_ECE": float(records[0]["inner_mean_ECE"]),
                "best_n_genes": int(records[0]["n_genes"]),
                "best_n_edges": int(records[0]["n_edges"]),
                "backend": str(records[0]["backend"]),
            }
        )
        next_population = [record["individual"] for record in elites]
        while len(next_population) < cfg.population_size:
            if rng.random() < cfg.crossover_rate and len(elites) >= 2:
                p1 = rng.choice(elites)["individual"]
                p2 = rng.choice(elites)["individual"]
                child = _crossover_individuals(p1, p2, rng, state_scores, cfg)
            else:
                child = rng.choice(elites)["individual"]
            child = _mutate_individual(child, rng, pools, state_scores, neighbors, cfg)
            next_population.append(child)
        population = next_population
    assert best_record is not None
    genes_by_state = {state: list(genes) for state, genes in best_record["genes_by_state"].items()}
    edges = list(best_record["edges"])
    feature_by_cohort = {
        cohort: build_ecology_features_from_module(
            ranked_expression_by_cohort[cohort].reindex(y_response[cohort].index),
            genes_by_state,
            edges,
            gene_directions,
            signed=signed,
            include_interactions=include_interactions,
        )
        for cohort in train_cohorts
    }
    X_feature = pd.concat(feature_by_cohort.values(), axis=0)
    X_feature = X_feature.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    y_feature = pd.concat([y_response[cohort].reindex(feature_by_cohort[cohort].index) for cohort in train_cohorts], axis=0).astype(int)
    coefficients: dict[str, float] = {}
    if y_feature.nunique() >= 2 and not X_feature.empty:
        model = Pipeline(
            [
                ("scale", StandardScaler()),
                (
                    "logistic",
                    LogisticRegression(
                        C=0.75,
                        class_weight="balanced",
                        solver="lbfgs",
                        max_iter=3000,
                        random_state=cfg.random_state,
                    ),
                ),
            ]
        )
        model.fit(X_feature.astype(float), y_feature)
        logistic = model.named_steps["logistic"]
        coefficients = {feature: float(value) for feature, value in zip(X_feature.columns, logistic.coef_[0])}
        coefficients["intercept"] = float(logistic.intercept_[0])
    max_count = max(selection_counts.values()) if selection_counts else 1.0
    selection_frequency = {gene: float(count / max_count) for gene, count in selection_counts.items()}
    gene_rows = []
    for state, genes in genes_by_state.items():
        for gene in genes:
            gene_rows.append(
                {
                    "state": state,
                    "gene": gene,
                    "direction": int(gene_directions.get(gene, 1)),
                    "selection_frequency": float(selection_frequency.get(gene, 0.0)),
                    "training_abs_correlation": float(abs(correlations.get(gene, 0.0))),
                    "direction_stability": float(stability.get(gene, 0.5)),
                    "state_prior_seed": int(gene in set(map(str, state_gene_sets.get(state, [])))),
                    "state_candidate_score": float(state_scores.get(state, {}).get(gene, global_scores.get(gene, 0.0))),
                }
            )
    edge_rows = [
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
    objective_terms = {key: value for key, value in best_record.items() if isinstance(value, (int, float, np.integer, np.floating))}
    backend = str(best_record.get("backend", "sklearn_logistic_cpu"))
    return HeuristicEcologyResult(
        genes_by_state=genes_by_state,
        edges=edges,
        gene_directions=gene_directions,
        feature_coefficients=coefficients,
        selection_frequency=selection_frequency,
        objective_terms=objective_terms,
        history=pd.DataFrame(history_rows),
        gene_table=pd.DataFrame(gene_rows),
        edge_table=pd.DataFrame(edge_rows),
        backend=backend,
    )
