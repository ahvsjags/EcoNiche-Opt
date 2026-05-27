from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche_opt.reporting.reproducibility import generate_reproducibility_report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    mode = "demo" if args.demo else "real"
    out = ROOT / (args.out or f"paper/{mode}_reproducibility_report.md")
    generate_reproducibility_report(out, mode=mode)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
