from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.labels import harmonize_metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/metadata/all_geo_samples_raw.tsv")
    parser.add_argument("--endpoint", default="primary_recist")
    parser.add_argument("--out", default="data/metadata/metadata_harmonized.tsv")
    parser.add_argument("--needs-manual", default="data/metadata/needs_manual_curation.tsv")
    args = parser.parse_args()
    input_path = ROOT / args.input
    out_path = ROOT / args.out
    manual_path = ROOT / args.needs_manual
    out_path.parent.mkdir(parents=True, exist_ok=True)
    manual_path.parent.mkdir(parents=True, exist_ok=True)
    if input_path.exists():
        metadata = pd.read_csv(input_path, sep="\t")
    else:
        metadata = pd.DataFrame(columns=["sample_id", "patient_id_raw", "cohort", "accession", "response_raw"])
    harmonized = harmonize_metadata(metadata, endpoint=args.endpoint)
    harmonized.to_csv(out_path, sep="\t", index=False)
    harmonized[harmonized["label"].isna()].to_csv(manual_path, sep="\t", index=False)
    print(f"Wrote harmonized metadata to {out_path}")


if __name__ == "__main__":
    main()
