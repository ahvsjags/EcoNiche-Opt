from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.registry import load_registry, write_registry_report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", default="config/data_registry.yml")
    parser.add_argument("--out", default="tables/dataset_access_audit.tsv")
    args = parser.parse_args()
    registry = load_registry(ROOT / args.registry)
    audit = write_registry_report(registry, ROOT / args.out)
    unknown = int((audit["access_status"] == "unknown").sum())
    controlled = int((audit["access_status"] == "controlled").sum())
    print(f"Wrote registry audit with {len(audit)} cohorts, {unknown} unknown, {controlled} controlled.")


if __name__ == "__main__":
    main()
