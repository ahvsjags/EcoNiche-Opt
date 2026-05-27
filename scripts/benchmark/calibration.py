from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.metrics import calibration_bins


def main() -> None:
    predictions = ROOT / "results/real/lodo_predictions.tsv"
    out = ROOT / "results/real/calibration_bins.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    if not predictions.exists():
        pd.DataFrame([{"status": "RESULT_PENDING", "reason": "missing lodo_predictions.tsv"}]).to_csv(out, sep="\t", index=False)
        return
    pred = pd.read_csv(predictions, sep="\t")
    frames = []
    for cohort, frame in pred.groupby("cohort"):
        bins = calibration_bins(frame["true_label"], frame["pred_prob"])
        bins["cohort"] = cohort
        frames.append(bins)
    pd.concat(frames, ignore_index=True).to_csv(out, sep="\t", index=False)
    print(f"Wrote calibration bins to {out}")


if __name__ == "__main__":
    main()
