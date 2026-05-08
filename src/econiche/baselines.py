from __future__ import annotations

import numpy as np
import pandas as pd


BASELINE_SIGNATURES = {
    "PDL1_CD274": ["CD274"],
    "PDCD1": ["PDCD1"],
    "PDCD1LG2": ["PDCD1LG2"],
    "CXCL9": ["CXCL9"],
    "HLA_DRA": ["HLA-DRA"],
    "CTLA4": ["CTLA4"],
    "CYT": ["GZMA", "PRF1"],
    "IFNG": ["IFNG", "CXCL9", "CXCL10", "STAT1"],
    "TIG": ["CD8A", "GZMA", "GZMB", "CXCL9", "CXCL10", "IDO1", "IFNG", "STAT1"],
    "TLS": ["CXCL13", "MS4A1", "CD79A", "CCR7"],
    "APM": ["HLA-A", "HLA-B", "B2M", "TAP1", "TAP2"],
    "IPRES": ["AXL", "VIM", "ZEB1", "COL1A1", "FN1"],
    "TIDE_dysfunction": ["PDCD1", "LAG3", "HAVCR2", "TOX"],
    "TIDE_exclusion": ["TGFB1", "COL1A1", "CXCL12", "FAP"],
    "MPS": ["S100A8", "S100A9", "IL10", "CD163"],
    "C_ECM": ["COL1A1", "COL3A1", "FN1", "TGFB1"],
    "ESCS": ["AXL", "NGFR", "ITGA3"],
    "IMPRES_template": ["PDCD1", "CTLA4", "LAG3", "HAVCR2"],
    "MCP_CD8_T": ["CD8A", "CD8B", "GZMB", "NKG7"],
    "MCP_fibroblast": ["COL1A1", "COL3A1", "FAP"],
    "ElasticNet_template": [],
    "RandomForest_template": [],
}


def signature_score(X: pd.DataFrame, genes: list[str]) -> pd.Series:
    available = [gene for gene in genes if gene in X.columns]
    if not available:
        return pd.Series(np.nan, index=X.index)
    return X[available].astype(float).mean(axis=1)


def score_baselines(X: pd.DataFrame, metadata: pd.DataFrame, signatures: dict[str, list[str]] | None = None) -> pd.DataFrame:
    signatures = signatures or BASELINE_SIGNATURES
    rows = []
    for name, genes in signatures.items():
        score = signature_score(X, genes)
        available = [gene for gene in genes if gene in X.columns]
        if score.notna().any():
            centered = score.fillna(score.mean()) - score.mean()
            scale = score.std() if score.std() and not np.isnan(score.std()) else 1.0
            prob = 1 / (1 + np.exp(-(centered / scale)))
            status = "available"
            reason = ""
        else:
            prob = pd.Series(np.nan, index=X.index)
            status = "unavailable_with_reason"
            reason = "required_genes_absent"
        for sample_id in X.index:
            meta = metadata.loc[sample_id] if sample_id in metadata.index else {}
            rows.append(
                {
                    "sample_id": sample_id,
                    "patient_id": meta.get("patient_id", pd.NA) if hasattr(meta, "get") else pd.NA,
                    "cohort": meta.get("cohort", pd.NA) if hasattr(meta, "get") else pd.NA,
                    "model_name": name,
                    "score": score.loc[sample_id],
                    "pred_prob": prob.loc[sample_id],
                    "pred_label": int(prob.loc[sample_id] >= 0.5) if pd.notna(prob.loc[sample_id]) else pd.NA,
                    "true_label": meta.get("label", pd.NA) if hasattr(meta, "get") else pd.NA,
                    "analysis_endpoint": "primary_recist",
                    "status": status,
                    "unavailable_reason": reason,
                    "n_genes_available": len(available),
                }
            )
    return pd.DataFrame(rows)
