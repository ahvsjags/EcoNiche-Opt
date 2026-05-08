from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import pandas as pd


DEFAULT_STATES = [
    "tumor_dedifferentiation",
    "antigen_presentation_mhc",
    "tnk_effector",
    "tcell_dysfunction",
    "caf_ecm_exclusion",
    "myeloid_suppression",
]


@dataclass
class EcoNicheConfig:
    states: list[str] = field(default_factory=lambda: list(DEFAULT_STATES))
    min_genes_per_state: int = 3
    max_genes_per_state: int = 25
    population_size: int = 120
    generations: int = 120
    elite_fraction: float = 0.10
    mutation_rate: float = 0.20
    crossover_rate: float = 0.50
    robust_rho: float = 0.50
    random_state: int = 42
    w_auc: float = 1.0
    w_auprc: float = 0.20
    w_cindex: float = 0.10
    w_ece: float = 0.15
    w_cell_specificity: float = 0.25
    w_pathway: float = 0.12
    w_network: float = 0.12
    w_lr: float = 0.12
    w_stability: float = 0.25
    w_size: float = 0.10
    w_batch: float = 0.20
    w_leakage: float = 1.00
    w_redundancy: float = 0.08
    w_therapy_confounding: float = 0.10


@dataclass
class EcoNicheModule:
    genes_by_state: dict[str, set[str]]
    edges_by_state_pair: dict[tuple[str, str], set[tuple[str, str]]] = field(default_factory=dict)

    def all_genes(self) -> set[str]:
        genes: set[str] = set()
        for state_genes in self.genes_by_state.values():
            genes.update(state_genes)
        return genes

    def size(self) -> int:
        return sum(len(genes) for genes in self.genes_by_state.values())

    def module_table(
        self,
        gene_directions: dict[str, int] | None = None,
        selection_frequency: dict[str, float] | None = None,
        coefficients: dict[str, float] | None = None,
        priors: pd.DataFrame | None = None,
    ) -> pd.DataFrame:
        rows: list[dict[str, Any]] = []
        gene_directions = gene_directions or {}
        selection_frequency = selection_frequency or {}
        coefficients = coefficients or {}
        for state in sorted(self.genes_by_state):
            for gene in sorted(self.genes_by_state[state]):
                prior_value = None
                if priors is not None and gene in priors.index and state in priors.columns:
                    prior_value = float(priors.loc[gene, state])
                rows.append(
                    {
                        "state": state,
                        "gene": gene,
                        "direction": int(gene_directions.get(gene, 1)),
                        "selection_frequency": float(selection_frequency.get(gene, 1.0)),
                        "cell_state_prior": prior_value,
                        "coefficient": coefficients.get(state),
                    }
                )
        return pd.DataFrame(rows)


@dataclass
class EcoNicheResult:
    best_module: EcoNicheModule
    objective_terms: dict[str, Any]
    lodo_metrics: pd.DataFrame
    predictions: pd.DataFrame
    history: pd.DataFrame
    coefficients: pd.DataFrame
    gene_directions: dict[str, int] = field(default_factory=dict)
    selection_frequency: dict[str, float] = field(default_factory=dict)
    priors: pd.DataFrame | None = None

    def module_table(self) -> pd.DataFrame:
        coefficient_map = {}
        if not self.coefficients.empty and {"feature", "coefficient"}.issubset(self.coefficients.columns):
            coefficient_map = dict(zip(self.coefficients["feature"], self.coefficients["coefficient"]))
        return self.best_module.module_table(
            gene_directions=self.gene_directions,
            selection_frequency=self.selection_frequency,
            coefficients=coefficient_map,
            priors=self.priors,
        )
