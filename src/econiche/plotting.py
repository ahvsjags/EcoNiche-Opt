from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import pandas as pd


def write_pending_figure(path: str | Path, title: str, message: str = "RESULT_PENDING") -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.axis("off")
    ax.text(0.5, 0.62, title, ha="center", va="center", fontsize=14, weight="bold")
    ax.text(0.5, 0.42, message, ha="center", va="center", fontsize=11)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def plot_metric_bar(metrics: pd.DataFrame, out: str | Path) -> None:
    if metrics.empty or "AUROC" not in metrics.columns:
        write_pending_figure(out, "Benchmark performance")
        return
    fig, ax = plt.subplots(figsize=(8, 4))
    metrics.plot.bar(x="cohort", y="AUROC", ax=ax, legend=False)
    ax.set_ylim(0, 1)
    ax.set_ylabel("AUROC")
    fig.tight_layout()
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
