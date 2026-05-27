from __future__ import annotations

import pandas as pd

from scripts.analysis.run_discovery_only_rescue_combo_search import combo_specs, primary_gene_universe


def test_combo_specs_are_primary_gene_screen_derived_and_include_base():
    top_genes = pd.DataFrame(
        [
            {"method": "pct", "gene": "TPI1P3", "primary_AUROC": 0.7},
            {"method": "pct", "gene": "TPI1P3", "primary_AUROC": 0.69},
            {"method": "z", "gene": "PLA2G2D", "primary_AUROC": 0.68},
        ]
    )

    specs = combo_specs(top_genes, top_k=3)
    candidates = {spec["candidate"] for spec in specs}

    assert "base_rescue_robust" in candidates
    assert "0.80*base+0.20*pct__TPI1P3" in candidates
    assert "0.80*base+0.20*z__PLA2G2D" in candidates


def test_primary_gene_universe_does_not_require_external_cohorts():
    X_by_cohort = {
        "GSE91061": pd.DataFrame(columns=["MAP4K1", "TBX3", "AXL", "TRAINONLY"]),
        "GSE78220": pd.DataFrame(columns=["MAP4K1", "TBX3", "AXL", "TRAINONLY"]),
        "PRJEB23709_PD1_PRE": pd.DataFrame(columns=["MAP4K1", "TBX3", "AXL", "TRAINONLY"]),
    }

    genes = primary_gene_universe(X_by_cohort)

    assert {"MAP4K1", "TBX3", "AXL", "TRAINONLY"}.issubset(set(genes))
