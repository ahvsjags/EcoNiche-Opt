from __future__ import annotations

import numpy as np
import pandas as pd

from econiche.module import EcoNicheModule
from econiche.normalize import rank_gaussian_normalize


def estimate_gene_directions(X: pd.DataFrame, y: pd.Series | np.ndarray) -> dict[str, int]:
    labels = pd.Series(y, index=X.index, dtype=float)
    directions: dict[str, int] = {}
    for gene in X.columns:
        values = X[gene].astype(float)
        if values.nunique(dropna=True) < 2 or labels.nunique(dropna=True) < 2:
            directions[gene] = 1
            continue
        corr = values.corr(labels)
        directions[gene] = -1 if pd.notna(corr) and corr < 0 else 1
    return directions


def compute_state_scores(
    X: pd.DataFrame,
    module: EcoNicheModule,
    gene_directions: dict[str, int] | None = None,
    normalize: bool = True,
) -> pd.DataFrame:
    gene_directions = gene_directions or {}
    features = rank_gaussian_normalize(X) if normalize else X.astype(float)
    scores = pd.DataFrame(index=X.index)
    for state, genes in module.genes_by_state.items():
        available = [gene for gene in sorted(genes) if gene in features.columns]
        if not available:
            scores[state] = 0.0
            continue
        signed = features.loc[:, available].copy()
        for gene in available:
            signed[gene] = signed[gene] * int(gene_directions.get(gene, 1))
        scores[state] = signed.sum(axis=1) / np.sqrt(len(available))
    return scores


def gene_label_correlations(X: pd.DataFrame, y: pd.Series | np.ndarray) -> pd.Series:
    labels = pd.Series(y, index=X.index, dtype=float)
    values = {}
    for gene in X.columns:
        corr = X[gene].corr(labels)
        values[gene] = 0.0 if pd.isna(corr) else float(corr)
    return pd.Series(values).sort_values(key=lambda s: s.abs(), ascending=False)
