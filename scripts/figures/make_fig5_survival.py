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

path = ROOT / "results/real/survival_km.tsv"
out = ROOT / "figures/fig5_survival.pdf"
if path.exists():
    data = pd.read_csv(path, sep="\t")
    if {"risk_group", "median_rfs_months"}.issubset(data.columns):
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.bar(data["risk_group"], data["median_rfs_months"], color=["#4c78a8", "#f58518"][: len(data)])
        ax.set_ylabel("Median RFS (months)")
        ax.set_title("GSE183924 survival association")
        fig.tight_layout()
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out)
        plt.close(fig)
    else:
        write_pending_figure(out, "Survival and clinical association", "RESULT_PENDING: survival columns unavailable")
else:
    write_pending_figure(out, "Survival and clinical association", "RESULT_PENDING: survival columns unavailable")
print("Wrote figures/fig5_survival.pdf")
