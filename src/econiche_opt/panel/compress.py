from __future__ import annotations

import pandas as pd


def compress_module_panel(module_table: pd.DataFrame, max_genes: int = 24) -> pd.DataFrame:
    if module_table.empty:
        return module_table.copy()
    ranked = module_table.copy()
    if "selection_frequency" in ranked.columns:
        ranked["_rank"] = pd.to_numeric(ranked["selection_frequency"], errors="coerce").fillna(0.0)
    else:
        ranked["_rank"] = 1.0
    return (
        ranked.sort_values(["_rank", "state", "gene"], ascending=[False, True, True])
        .drop(columns=["_rank"])
        .head(max_genes)
        .reset_index(drop=True)
    )
