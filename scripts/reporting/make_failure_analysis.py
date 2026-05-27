from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    mode = "demo" if args.demo else "real"
    pred_path = ROOT / f"results/{mode}/lodo_predictions.tsv"
    out = ROOT / (args.out or f"results/{mode}_failure_analysis.md")
    out.parent.mkdir(parents=True, exist_ok=True)
    if pred_path.exists():
        pred = pd.read_csv(pred_path, sep="\t")
        if {"true_label", "pred_label"}.issubset(pred.columns):
            pred["is_error"] = pred["true_label"].astype(str) != pred["pred_label"].astype(str)
            summary = pred.groupby("cohort")["is_error"].agg(["sum", "count"]).reset_index()
            lines = ["# Failure Analysis", "", "```text", summary.to_string(index=False), "```"]
        else:
            lines = ["# Failure Analysis", "", "RESULT_PENDING: prediction labels are not available."]
    else:
        lines = ["# Failure Analysis", "", "RESULT_PENDING: predictions not found."]
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
