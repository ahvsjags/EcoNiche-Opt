from __future__ import annotations

import numpy as np

from econiche.module import EcoNicheConfig, EcoNicheModule


def module_size_penalty(module: EcoNicheModule, cfg: EcoNicheConfig) -> float:
    penalties = []
    for genes in module.genes_by_state.values():
        size = len(genes)
        below = max(0, cfg.min_genes_per_state - size) / max(cfg.min_genes_per_state, 1)
        above = max(0, size - cfg.max_genes_per_state) / max(cfg.max_genes_per_state, 1)
        penalties.append(below + above)
    return float(np.mean(penalties)) if penalties else 0.0


def robust_objective(
    cfg: EcoNicheConfig,
    auc_values,
    auprc_values,
    ece_mean: float = 0.0,
    cindex_mean: float = 0.0,
    biological_terms: dict[str, float] | None = None,
    penalties: dict[str, float] | None = None,
) -> dict[str, float]:
    biological_terms = biological_terms or {}
    penalties = penalties or {}
    auc = np.asarray(list(auc_values), dtype=float)
    auprc = np.asarray(list(auprc_values), dtype=float)
    auc_mean = float(np.nanmean(auc)) if len(auc) else float("nan")
    auc_sd = float(np.nanstd(auc)) if len(auc) else 0.0
    auprc_mean = float(np.nanmean(auprc)) if len(auprc) else 0.0
    terms = {
        "auc_mean": auc_mean,
        "auc_sd": auc_sd,
        "auprc_mean": auprc_mean,
        "cindex_mean": float(cindex_mean or 0.0),
        "ece_mean": float(ece_mean or 0.0),
        "cell_specificity": float(biological_terms.get("cell_specificity", 0.0)),
        "pathway": float(biological_terms.get("pathway", 0.0)),
        "network": float(biological_terms.get("network", 0.0)),
        "lr": float(biological_terms.get("lr", 0.0)),
        "stability": float(biological_terms.get("stability", 0.0)),
        "size": float(penalties.get("size", 0.0)),
        "batch": float(penalties.get("batch", 0.0)),
        "leakage": float(penalties.get("leakage", 0.0)),
        "redundancy": float(penalties.get("redundancy", 0.0)),
        "therapy_confounding": float(penalties.get("therapy_confounding", 0.0)),
    }
    score = (
        cfg.w_auc * (terms["auc_mean"] - cfg.robust_rho * terms["auc_sd"])
        + cfg.w_auprc * terms["auprc_mean"]
        + cfg.w_cindex * terms["cindex_mean"]
        - cfg.w_ece * terms["ece_mean"]
        + cfg.w_cell_specificity * terms["cell_specificity"]
        + cfg.w_pathway * terms["pathway"]
        + cfg.w_network * terms["network"]
        + cfg.w_lr * terms["lr"]
        + cfg.w_stability * terms["stability"]
        - cfg.w_size * terms["size"]
        - cfg.w_batch * terms["batch"]
        - cfg.w_leakage * terms["leakage"]
        - cfg.w_redundancy * terms["redundancy"]
        - cfg.w_therapy_confounding * terms["therapy_confounding"]
    )
    terms["score"] = float(score)
    return terms
