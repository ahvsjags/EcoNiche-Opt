from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.priors import make_default_cell_state_priors


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genes", default="data/priors/gene_universe.txt")
    parser.add_argument("--out", default="data/priors/cell_state_priors.tsv")
    args = parser.parse_args()
    gene_path = ROOT / args.genes
    genes = [line.strip() for line in gene_path.read_text(encoding="utf-8").splitlines() if line.strip()] if gene_path.exists() else []
    priors = make_default_cell_state_priors(genes)
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    priors.to_csv(out, sep="\t")
    print(f"Wrote cell-state priors to {out}")


if __name__ == "__main__":
    main()
