from __future__ import annotations

import numpy as np
import pandas as pd

from econiche.priors import DEFAULT_MARKERS


PROTECTIVE_STATES = {"antigen_presentation_mhc", "tnk_effector"}


def make_synthetic_data(
    n_cohorts: int = 4,
    n_per_cohort: int = 48,
    n_noise_genes: int = 70,
    random_state: int = 42,
) -> dict:
    rng = np.random.default_rng(random_state)
    planted_genes = [gene for genes in DEFAULT_MARKERS.values() for gene in genes[:3]]
    noise_genes = [f"NOISE{i:03d}" for i in range(n_noise_genes)]
    genes = planted_genes + noise_genes
    X_by_cohort: dict[str, pd.DataFrame] = {}
    y_by_cohort: dict[str, pd.Series] = {}
    metadata_by_cohort: dict[str, pd.DataFrame] = {}

    for cohort_idx in range(n_cohorts):
        cohort = f"demo_cohort_{cohort_idx + 1}"
        y = np.array([0, 1] * (n_per_cohort // 2) + ([0] if n_per_cohort % 2 else []))
        rng.shuffle(y)
        X = rng.normal(0, 1, size=(n_per_cohort, len(genes)))
        gene_to_idx = {gene: idx for idx, gene in enumerate(genes)}
        signal = 1.25 + cohort_idx * 0.05
        for state, markers in DEFAULT_MARKERS.items():
            for gene in markers[:3]:
                idx = gene_to_idx[gene]
                if state in PROTECTIVE_STATES:
                    X[:, idx] += np.where(y == 1, -signal, signal)
                else:
                    X[:, idx] += np.where(y == 1, signal, -signal)
        sample_ids = [f"{cohort}_S{i:03d}" for i in range(n_per_cohort)]
        X_df = pd.DataFrame(X, index=sample_ids, columns=genes)
        y_series = pd.Series(y, index=sample_ids, name="label")
        meta = pd.DataFrame(
            {
                "sample_id": sample_ids,
                "patient_id": [f"{cohort}_P{i:03d}" for i in range(n_per_cohort)],
                "cohort": cohort,
                "accession": cohort,
                "cancer_type": "melanoma",
                "therapy": "anti-PD1",
                "platform": "synthetic_RNAseq",
                "timepoint": "pretreatment",
                "response_raw": np.where(y == 1, "PD", "PR"),
                "label": y,
            }
        ).set_index("sample_id", drop=False)
        X_by_cohort[cohort] = X_df
        y_by_cohort[cohort] = y_series
        metadata_by_cohort[cohort] = meta
    return {
        "X_by_cohort": X_by_cohort,
        "y_by_cohort": y_by_cohort,
        "metadata_by_cohort": metadata_by_cohort,
        "genes": genes,
        "planted_genes": planted_genes,
    }
