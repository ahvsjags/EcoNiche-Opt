from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche_opt.preprocess.deconvolution import (  # noqa: E402
    evaluate_abundance_baselines,
    score_processed_cohorts,
    summarize_abundance,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run marker-based cell abundance baselines on processed bulk cohorts.")
    parser.add_argument("--input-dir", default="data/processed/bulk")
    parser.add_argument("--out-dir", default="results/real")
    parser.add_argument("--include-demo", action="store_true")
    parser.add_argument("--min-markers", type=int, default=2)
    args = parser.parse_args()

    out_dir = ROOT / args.out_dir if not Path(args.out_dir).is_absolute() else Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    scores, metadata_by_cohort = score_processed_cohorts(
        ROOT / args.input_dir if not Path(args.input_dir).is_absolute() else Path(args.input_dir),
        include_demo=args.include_demo,
        min_markers=args.min_markers,
    )
    summary = summarize_abundance(scores)
    metrics = evaluate_abundance_baselines(scores, metadata_by_cohort)

    scores_out = out_dir / "deconvolution_scores.tsv"
    summary_out = out_dir / "deconvolution_cohort_summary.tsv"
    metrics_out = out_dir / "deconvolution_baseline_metrics.tsv"
    scores.to_csv(scores_out, sep="\t", index=False)
    summary.to_csv(summary_out, sep="\t", index=False)
    metrics.to_csv(metrics_out, sep="\t", index=False)

    if scores.empty:
        print("No processed expression cohorts were found.")
    else:
        print(f"Scored {scores['cohort'].nunique()} cohorts and {scores['sample_id'].nunique()} samples.")
    print(f"Wrote {scores_out}")
    print(f"Wrote {summary_out}")
    print(f"Wrote {metrics_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
