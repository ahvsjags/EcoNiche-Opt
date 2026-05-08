from __future__ import annotations

import pandas as pd

from econiche.module import EcoNicheModule


def pathway_coherence(module: EcoNicheModule, pathways: dict[str, set[str]] | None) -> float:
    if not pathways:
        return 0.0
    scores = []
    for genes in module.genes_by_state.values():
        if not genes:
            continue
        best = max((len(set(genes) & set(pathway_genes)) / len(genes) for pathway_genes in pathways.values()), default=0.0)
        scores.append(best)
    return float(sum(scores) / len(scores)) if scores else 0.0


def network_coherence(module: EcoNicheModule, edges: pd.DataFrame | set[tuple[str, str]] | None) -> float:
    if edges is None:
        return 0.0
    if isinstance(edges, pd.DataFrame):
        edge_set = {tuple(sorted((str(row[0]), str(row[1])))) for row in edges.iloc[:, :2].itertuples(index=False)}
    else:
        edge_set = {tuple(sorted(edge)) for edge in edges}
    scores = []
    for genes in module.genes_by_state.values():
        genes = sorted(genes)
        if len(genes) < 2:
            continue
        possible = len(genes) * (len(genes) - 1) / 2
        observed = sum(1 for i, g1 in enumerate(genes) for g2 in genes[i + 1 :] if tuple(sorted((g1, g2))) in edge_set)
        scores.append(observed / possible)
    return float(sum(scores) / len(scores)) if scores else 0.0


def lr_coherence(
    module: EcoNicheModule,
    lr_edges: pd.DataFrame | set[tuple[str, str]] | None,
    state_pairs: list[tuple[str, str]] | None = None,
) -> float:
    if lr_edges is None:
        return 0.0
    if isinstance(lr_edges, pd.DataFrame):
        edge_set = {(str(row[0]), str(row[1])) for row in lr_edges.iloc[:, :2].itertuples(index=False)}
    else:
        edge_set = set(lr_edges)
    pairs = state_pairs or list(module.edges_by_state_pair)
    if not pairs:
        pairs = [(a, b) for idx, a in enumerate(module.genes_by_state) for b in list(module.genes_by_state)[idx + 1 :]]
    observed = 0
    possible = 0
    for left_state, right_state in pairs:
        for left_gene in module.genes_by_state.get(left_state, set()):
            for right_gene in module.genes_by_state.get(right_state, set()):
                possible += 1
                if (left_gene, right_gene) in edge_set or (right_gene, left_gene) in edge_set:
                    observed += 1
    return float(observed / possible) if possible else 0.0
