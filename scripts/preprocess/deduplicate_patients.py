from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.qc import split_by_timepoint_priority


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/metadata/metadata_harmonized.tsv")
    parser.add_argument("--primary", default="data/metadata/metadata_dedup_primary.tsv")
    parser.add_argument("--secondary", default="data/metadata/metadata_secondary_on_treatment.tsv")
    parser.add_argument("--progression", default="data/metadata/metadata_progression.tsv")
    parser.add_argument("--report", default="tables/deduplication_report.tsv")
    args = parser.parse_args()
    input_path = ROOT / args.input
    metadata = pd.read_csv(input_path, sep="\t") if input_path.exists() else pd.DataFrame()
    primary, secondary, progression = split_by_timepoint_priority(metadata)
    for path, frame in [(args.primary, primary), (args.secondary, secondary), (args.progression, progression)]:
        out = ROOT / path
        out.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(out, sep="\t", index=False)
    report = pd.DataFrame(
        [
            {"set": "primary", "n_samples": len(primary), "n_patients": primary.get("patient_id", pd.Series(dtype=str)).nunique()},
            {"set": "secondary_on_treatment", "n_samples": len(secondary), "n_patients": secondary.get("patient_id", pd.Series(dtype=str)).nunique()},
            {"set": "progression", "n_samples": len(progression), "n_patients": progression.get("patient_id", pd.Series(dtype=str)).nunique()},
        ]
    )
    report_path = ROOT / args.report
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(report_path, sep="\t", index=False)
    print(f"Wrote deduplication report to {report_path}")


if __name__ == "__main__":
    main()
