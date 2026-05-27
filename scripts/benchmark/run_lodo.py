from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohorts", default=None)
    parser.add_argument("--endpoint", default="primary_recist")
    parser.add_argument("--out", default="results/real")
    parser.add_argument("--demo", action="store_true")
    args = parser.parse_args()
    cmd = [
        sys.executable,
        str(ROOT / "scripts/model/run_econiche.py"),
        "--config",
        "config/model_config.yml",
        "--out",
        args.out,
    ]
    if args.demo:
        cmd.append("--demo")
    subprocess.run(cmd, cwd=ROOT, check=True)
    out_dir = ROOT / args.out
    for src_name, dest_name in [
        ("lodo_predictions.tsv", "lodo_predictions.tsv"),
        ("lodo_metrics.tsv", "lodo_metrics_by_cohort.tsv"),
        ("objective_history.tsv", "lodo_objective_history.tsv"),
    ]:
        src = out_dir / src_name
        dest = out_dir / dest_name
        if src.exists() and src != dest:
            shutil.copyfile(src, dest)
    print(f"Wrote LODO benchmark outputs to {out_dir}")


if __name__ == "__main__":
    main()
