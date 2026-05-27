from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.plotting import plot_metric_bar


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--out", default="figures/figure3_benchmark.svg")
    args = parser.parse_args()
    metrics_path = ROOT / ("results/demo/lodo_metrics.tsv" if args.demo else "results/real/lodo_metrics.tsv")
    metrics = pd.read_csv(metrics_path, sep="\t") if metrics_path.exists() else pd.DataFrame()
    plot_metric_bar(metrics, ROOT / args.out)
    print(f"Wrote {ROOT / args.out}")


if __name__ == "__main__":
    main()
