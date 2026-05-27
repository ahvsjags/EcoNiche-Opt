from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def _copy_if_exists(src: Path, dest: Path) -> str:
    if not src.exists():
        return "RESULT_PENDING"
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dest)
    return "copied_from_pipeline"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    source = ROOT / ("results/demo" if args.demo else "results/real")
    out_dir = ROOT / (args.out or ("tables/demo" if args.demo else "tables/real"))
    rows = []
    for name in ["lodo_metrics.tsv", "lodo_predictions.tsv", "econiche_module.tsv", "objective_history.tsv"]:
        status = _copy_if_exists(source / name, out_dir / name)
        rows.append({"table": name, "status": status})
    manifest = pd.DataFrame(rows)
    manifest.to_csv(out_dir / "table_manifest.tsv", sep="\t", index=False)
    print(f"Wrote {out_dir / 'table_manifest.tsv'}")


if __name__ == "__main__":
    main()
