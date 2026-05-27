from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.plotting import plot_metric_bar

metrics_path = ROOT / "results/real/lodo_metrics.tsv"
if not metrics_path.exists():
    metrics_path = ROOT / "results/demo/lodo_metrics.tsv"
metrics = pd.read_csv(metrics_path, sep="\t") if metrics_path.exists() else pd.DataFrame()
plot_metric_bar(metrics, ROOT / "figures/fig2_benchmark.pdf")
print("Wrote figures/fig2_benchmark.pdf")
