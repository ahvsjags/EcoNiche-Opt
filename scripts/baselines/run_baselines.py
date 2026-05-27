from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.baselines import score_baselines
from econiche.io import load_processed_bulk


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/baselines.yml")
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--out", default="results/real/baseline_predictions.tsv")
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--include-demo", action="store_true")
    args = parser.parse_args()
    X_by_cohort, _, metadata_by_cohort = load_processed_bulk(ROOT / args.processed_dir)
    if args.demo:
        X_by_cohort = {k: v for k, v in X_by_cohort.items() if k.startswith("demo_cohort_")}
        metadata_by_cohort = {k: v for k, v in metadata_by_cohort.items() if k in X_by_cohort}
    elif not args.include_demo:
        X_by_cohort = {k: v for k, v in X_by_cohort.items() if not k.startswith("demo_cohort_")}
        metadata_by_cohort = {k: v for k, v in metadata_by_cohort.items() if k in X_by_cohort}
    frames = []
    for cohort, X in X_by_cohort.items():
        frames.append(score_baselines(X, metadata_by_cohort[cohort]))
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    predictions = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    predictions.to_csv(out, sep="\t", index=False)
    print(f"Wrote baseline predictions to {out}")


if __name__ == "__main__":
    main()
