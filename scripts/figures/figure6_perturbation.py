from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.plotting import write_pending_figure


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--out", default="figures/figure6_perturbation.svg")
    args = parser.parse_args()
    write_pending_figure(ROOT / args.out, "Perturbation prioritization", "candidate hypotheses only; no treatment recommendation")
    print(f"Wrote {ROOT / args.out}")


if __name__ == "__main__":
    main()
