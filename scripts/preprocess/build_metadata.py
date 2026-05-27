from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.geo import harmonized_geo_metadata, parse_series_matrix_metadata
from econiche.labels import harmonize_metadata
from econiche.registry import load_registry


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", default="data/raw")
    parser.add_argument("--registry", default="config/data_registry.yml")
    parser.add_argument("--out", default="data/metadata/all_geo_samples_raw.tsv")
    parser.add_argument("--needs-manual", default="data/metadata/needs_manual_curation.tsv")
    args = parser.parse_args()
    raw_dir = ROOT / args.raw_dir
    registry = load_registry(ROOT / args.registry)
    platform_by_accession = {cohort["accession"]: cohort.get("platform") for cohort in registry.get("cohorts", [])}
    frames = []
    for matrix_path in raw_dir.glob("GSE*/matrix/*series_matrix*.txt.gz"):
        accession = matrix_path.parts[-3]
        raw = parse_series_matrix_metadata(matrix_path)
        if raw.empty:
            continue
        frames.append(harmonized_geo_metadata(raw, accession, platform=platform_by_accession.get(accession)))
    for path in raw_dir.glob("*/metadata.tsv"):
        frame = pd.read_csv(path, sep="\t")
        frames.append(frame)
    columns = [
        "sample_id",
        "patient_id_raw",
        "cohort",
        "accession",
        "platform",
        "title",
        "source_name",
        "characteristics_ch1",
        "therapy",
        "timepoint",
        "response_raw",
        "response_harmonized",
        "pfs",
        "os",
        "pfs_event",
        "os_event",
    ]
    metadata = pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame(columns=columns)
    if not metadata.empty:
        metadata = metadata.drop_duplicates(["accession", "sample_id"], keep="first")
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(out, sep="\t", index=False)
    harmonized = harmonize_metadata(metadata)
    harmonized.to_csv(ROOT / "data/metadata/metadata_harmonized.tsv", sep="\t", index=False)
    needs = harmonized[harmonized.get("label", pd.Series(dtype=object)).isna()] if not harmonized.empty else harmonized
    manual = ROOT / args.needs_manual
    manual.parent.mkdir(parents=True, exist_ok=True)
    needs.to_csv(manual, sep="\t", index=False)
    print(f"Wrote raw metadata table to {out}")


if __name__ == "__main__":
    main()
