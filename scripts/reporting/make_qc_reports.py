from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    processed = ROOT / args.processed_dir
    rows = []
    for expr_path in sorted(processed.glob("*.expr.tsv")):
        if args.demo and not expr_path.name.startswith("demo_cohort_"):
            continue
        frame = pd.read_csv(expr_path, sep="\t", index_col=0)
        rows.append(
            {
                "cohort": expr_path.name.replace(".expr.tsv", ""),
                "n_samples": frame.shape[0],
                "n_genes": frame.shape[1],
                "missing_fraction": float(frame.isna().mean().mean()),
            }
        )
    out_dir = ROOT / (args.out or ("results/demo_qc" if args.demo else "results/real_qc"))
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_dir / "dataset_qc.tsv", sep="\t", index=False)
    print(f"Wrote {out_dir / 'dataset_qc.tsv'}")


if __name__ == "__main__":
    main()
