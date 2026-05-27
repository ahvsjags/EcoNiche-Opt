from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.priors import DEFAULT_MARKERS


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--out", default="data/priors/gene_universe.txt")
    args = parser.parse_args()
    processed = ROOT / args.processed_dir
    common = None
    for expr_path in processed.glob("*.expr.tsv"):
        X = pd.read_csv(expr_path, sep="\t", index_col=0, nrows=1)
        genes = set(X.columns)
        common = genes if common is None else common & genes
    universe = set(common or [])
    for markers in DEFAULT_MARKERS.values():
        universe.update(markers)
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(sorted(universe)) + ("\n" if universe else ""), encoding="utf-8")
    print(f"Wrote {len(universe)} genes to {out}")


if __name__ == "__main__":
    main()
