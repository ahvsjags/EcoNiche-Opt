from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche_opt.panel.compress import compress_module_panel


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--max-genes", type=int, default=24)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    results_dir = ROOT / ("results/demo" if args.demo else "results/real")
    module_path = results_dir / "econiche_module.tsv"
    out_dir = ROOT / (args.out or ("results/demo_panel" if args.demo else "results/real_panel"))
    out_dir.mkdir(parents=True, exist_ok=True)
    module = pd.read_csv(module_path, sep="\t") if module_path.exists() else pd.DataFrame()
    panel = compress_module_panel(module, max_genes=args.max_genes)
    out = out_dir / "compressed_panel.tsv"
    panel.to_csv(out, sep="\t", index=False)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
