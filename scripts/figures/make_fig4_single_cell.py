from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.plotting import write_pending_figure

path = ROOT / "results/scrna/cell_type_enrichment.tsv"
out = ROOT / "figures/fig4_single_cell.pdf"
if path.exists():
    data = pd.read_csv(path, sep="\t")
    if {"cell_type", "state", "mean"}.issubset(data.columns) and not data.empty:
        pivot = data.pivot_table(index="cell_type", columns="state", values="mean", aggfunc="mean").fillna(0)
        fig, ax = plt.subplots(figsize=(9, max(4, 0.35 * len(pivot))))
        im = ax.imshow(pivot.to_numpy(), aspect="auto", cmap="viridis")
        ax.set_xticks(range(len(pivot.columns)), pivot.columns, rotation=45, ha="right")
        ax.set_yticks(range(len(pivot.index)), pivot.index)
        ax.set_title("Single-cell EcoNiche module scores")
        fig.colorbar(im, ax=ax, label="Mean module score")
        fig.tight_layout()
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out)
        plt.close(fig)
    else:
        write_pending_figure(out, "Single-cell mechanism mapping", "RESULT_PENDING: run scRNA pipeline")
else:
    write_pending_figure(out, "Single-cell mechanism mapping", "RESULT_PENDING: run scRNA pipeline")
print("Wrote figures/fig4_single_cell.pdf")
