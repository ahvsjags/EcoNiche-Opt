from __future__ import annotations

from collections.abc import Iterable

import pandas as pd


def gene_coverage_report(module_genes: Iterable[str], available_genes: Iterable[str]) -> pd.DataFrame:
    module = sorted(set(module_genes))
    available = set(available_genes)
    rows = []
    for gene in module:
        rows.append({"gene": gene, "available": gene in available})
    report = pd.DataFrame(rows)
    if report.empty:
        report["gene"] = []
        report["available"] = []
    return report


def coverage_fraction(module_genes: Iterable[str], available_genes: Iterable[str]) -> float:
    report = gene_coverage_report(module_genes, available_genes)
    if report.empty:
        return 0.0
    return float(report["available"].mean())
