from __future__ import annotations

from pathlib import Path

import pandas as pd


def ensure_dir(path: str | Path) -> Path:
    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out


def load_processed_bulk(processed_dir: str | Path) -> tuple[dict[str, pd.DataFrame], dict[str, pd.Series], dict[str, pd.DataFrame]]:
    processed = Path(processed_dir)
    X_by_cohort: dict[str, pd.DataFrame] = {}
    y_by_cohort: dict[str, pd.Series] = {}
    metadata_by_cohort: dict[str, pd.DataFrame] = {}
    for expr_path in sorted(processed.glob("*.expr.tsv")):
        cohort = expr_path.name.replace(".expr.tsv", "")
        meta_path = processed / f"{cohort}.metadata.tsv"
        if not meta_path.exists():
            continue
        X = pd.read_csv(expr_path, sep="\t", index_col=0)
        if X.empty or X.shape[0] == 0 or X.shape[1] == 0:
            continue
        meta = pd.read_csv(meta_path, sep="\t")
        if "sample_id" not in meta.columns or "label" not in meta.columns:
            continue
        meta = meta.set_index("sample_id", drop=False).reindex(X.index)
        meta = meta[meta["label"].notna()]
        X = X.loc[meta.index]
        if X.empty or meta["label"].nunique() < 2:
            continue
        X_by_cohort[cohort] = X
        y_by_cohort[cohort] = meta["label"].astype(int)
        metadata_by_cohort[cohort] = meta
    return X_by_cohort, y_by_cohort, metadata_by_cohort
