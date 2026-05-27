from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--report", default="tables/nanostring_qc_report.tsv")
    args = parser.parse_args()
    rows = []
    for expr_path in (ROOT / args.processed_dir).glob("*.expr.tsv"):
        cohort = expr_path.name.replace(".expr.tsv", "")
        X = pd.read_csv(expr_path, sep="\t", index_col=0)
        rows.append({"cohort": cohort, "n_samples": X.shape[0], "n_panel_genes": X.shape[1], "status": "processed"})
    out = ROOT / args.report
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=["cohort", "n_samples", "n_panel_genes", "status"]).to_csv(out, sep="\t", index=False)
    print(f"Wrote NanoString QC report to {out}")


if __name__ == "__main__":
    main()
