from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche_opt.data.registry import audit_accession, load_registry


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", default="config/data_registry.yml")
    parser.add_argument("--out", default="results/audit/data_access_audit.tsv")
    args = parser.parse_args()

    report = audit_accession(load_registry(ROOT / args.registry))
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(out, sep="\t", index=False)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
