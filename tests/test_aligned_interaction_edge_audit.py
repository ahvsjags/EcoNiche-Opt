from __future__ import annotations

import pandas as pd

from scripts.analysis.run_aligned_interaction_edge_audit import build_edge_features


def test_build_edge_features_adds_response_resistance_edges():
    features = pd.DataFrame(
        {
            "ifn_t_cell_inflamed": [1.0, 2.0],
            "cytotoxic_cd8": [0.5, 1.0],
            "antigen_presentation": [3.0, 4.0],
            "exhaustion_checkpoint": [0.2, 0.3],
            "trm_tls": [0.1, 0.2],
            "myeloid_suppression": [2.0, 3.0],
            "stromal_exclusion": [4.0, 5.0],
        },
        index=["s1", "s2"],
    )

    out = build_edge_features(features)

    assert "edge__ifn_t_cell_inflamed__x__myeloid_suppression" in out.columns
    assert out.loc["s1", "edge__ifn_t_cell_inflamed__x__myeloid_suppression"] == 2.0
    assert "contrast__cytotoxic_cd8__minus__stromal_exclusion" in out.columns
    assert out.loc["s2", "contrast__cytotoxic_cd8__minus__stromal_exclusion"] == -4.0
    assert "synergy__ifn_t_cell_inflamed__x__antigen_presentation" in out.columns
