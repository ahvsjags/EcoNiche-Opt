from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche_opt.reporting.claim_gate import sanitize_claim


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    mode = "demo" if args.demo else "real"
    metrics_path = ROOT / f"results/{mode}/lodo_metrics.tsv"
    metrics_note = "RESULT_PENDING"
    if metrics_path.exists():
        metrics = pd.read_csv(metrics_path, sep="\t")
        metrics_note = f"Pipeline generated {len(metrics)} leave-one-dataset-out metric rows."
    text = "\n".join(
        [
            "# EcoNiche-Opt Manuscript Skeleton",
            "",
            "## Abstract",
            sanitize_claim("EcoNiche-Opt was evaluated in a reproducible multicohort benchmark without claiming superiority."),
            "",
            "## Methods",
            "All splits, thresholding, calibration, and model selection are constrained to training cohorts.",
            "",
            "## Results",
            metrics_note,
            "",
            "## Data Availability",
            "Controlled or license-restricted resources are marked ACCESS_RESTRICTED and are not replaced by fabricated data.",
        ]
    )
    out = ROOT / (args.out or f"paper/{mode}_manuscript.md")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
