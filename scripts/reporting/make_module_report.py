from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    mode = "demo" if args.demo else "real"
    module_path = ROOT / f"results/{mode}/econiche_module.tsv"
    out = ROOT / (args.out or f"results/{mode}_module_report.md")
    out.parent.mkdir(parents=True, exist_ok=True)
    if module_path.exists():
        module = pd.read_csv(module_path, sep="\t")
        counts = module.groupby("state")["gene"].nunique().reset_index(name="n_genes")
        lines = ["# Module Report", "", "```text", counts.to_string(index=False), "```"]
    else:
        lines = ["# Module Report", "", "RESULT_PENDING: module table not found."]
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
