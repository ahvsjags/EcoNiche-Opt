from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.analysis.run_locked_external_panel_validation import run


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--out", default="results/locked_external_panel_validation")
    args = parser.parse_args()

    processed_dir = ROOT / args.processed_dir
    out_dir = ROOT / args.out
    status_path = ROOT / "results/real/locked_external_status.tsv"
    status_path.parent.mkdir(parents=True, exist_ok=True)
    if not processed_dir.exists():
        status_path.write_text(
            "status\treason\toutput\nRESULT_PENDING\tProcessed public cohorts are not available yet\t\n",
            encoding="utf-8",
        )
        print(f"Wrote {status_path}")
        return

    run(
        processed_dir=processed_dir,
        out_dir=out_dir,
        endpoints=["primary_recist", "strict_recist", "clinical_benefit"],
        discovery_cohorts=["GSE91061", "GSE78220", "PRJEB23709_PD1_PRE"],
        external_cohorts=["GSE145996", "PHS000452_LIU_LIKE_PRE", "PRJEB23709_COMBO_PRE", "GSE93157", "GSE140901"],
    )
    status_path.write_text(
        f"status\treason\toutput\ncompleted\tLocked external and clinical-assay panel validation completed\t{out_dir}\n",
        encoding="utf-8",
    )
    print(f"Wrote {status_path}")


if __name__ == "__main__":
    main()
