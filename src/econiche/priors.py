from __future__ import annotations

import pandas as pd

from econiche.module import DEFAULT_STATES


DEFAULT_MARKERS = {
    "tumor_dedifferentiation": ["AXL", "NGFR", "ITGA3", "VIM", "ZEB1"],
    "antigen_presentation_mhc": ["HLA-A", "HLA-B", "B2M", "TAP1", "TAP2"],
    "tnk_effector": ["GZMB", "PRF1", "NKG7", "GNLY", "IFNG"],
    "tcell_dysfunction": ["PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX"],
    "caf_ecm_exclusion": ["COL1A1", "FN1", "TGFB1", "COL3A1", "CXCL12"],
    "myeloid_suppression": ["S100A8", "S100A9", "LILRB2", "IL10", "CD163"],
}


def make_default_cell_state_priors(genes: list[str] | set[str] | None = None) -> pd.DataFrame:
    gene_set = set(genes or [])
    if not gene_set:
        for markers in DEFAULT_MARKERS.values():
            gene_set.update(markers)
    priors = pd.DataFrame(0.0, index=sorted(gene_set), columns=list(DEFAULT_STATES))
    for state, markers in DEFAULT_MARKERS.items():
        for gene in markers:
            if gene in priors.index:
                priors.loc[gene, state] = 0.9
    return priors


def load_priors(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0)
