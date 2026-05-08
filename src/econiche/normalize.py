from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def rank_gaussian_normalize(X: pd.DataFrame) -> pd.DataFrame:
    ranks = X.rank(axis=1, method="average", pct=True)
    ranks = ranks.clip(lower=1e-4, upper=1 - 1e-4)
    values = stats.norm.ppf(ranks.to_numpy(dtype=float))
    return pd.DataFrame(values, index=X.index, columns=X.columns)


def log2_if_needed(X: pd.DataFrame) -> pd.DataFrame:
    numeric = X.astype(float)
    max_value = np.nanmax(numeric.to_numpy()) if numeric.size else 0
    min_value = np.nanmin(numeric.to_numpy()) if numeric.size else 0
    if max_value > 50 or min_value < 0:
        return np.log2(numeric.clip(lower=0) + 1)
    return numeric


def intersect_gene_space(matrices: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    if not matrices:
        return {}
    common = set(next(iter(matrices.values())).columns)
    for matrix in matrices.values():
        common &= set(matrix.columns)
    genes = sorted(common)
    return {cohort: matrix.loc[:, genes].copy() for cohort, matrix in matrices.items()}
