from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche_opt.perturbation.reversal import (  # noqa: E402
    REVERSAL_COLUMNS,
    aggregate_reversal_candidates,
    query_lincs_reversal,
)


def _pending(reason: str) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "query_direction": "RESULT_PENDING",
                "library": "RESULT_PENDING",
                "rank": "",
                "perturbation_id": "RESULT_PENDING",
                "perturbation_name": "RESULT_PENDING",
                "cell_line": "",
                "timepoint": "",
                "dose": "",
                "pvalue": "",
                "zscore": "",
                "combined_score": "",
                "adjusted_pvalue": "",
                "overlapping_genes": "",
                "signature_genes_submitted": "",
                "reversal_score": "",
                "status": "RESULT_PENDING",
                "interpretation": reason,
            }
        ],
        columns=REVERSAL_COLUMNS,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Query Enrichr LINCS L1000 libraries for module signature reversal.")
    parser.add_argument("--module", default="results/real/econiche_module.tsv")
    parser.add_argument("--out", default="results/perturbation/lincs_reversal.tsv")
    parser.add_argument("--summary-out", default="results/perturbation/prioritized_reversal.tsv")
    parser.add_argument("--top-genes", type=int, default=150)
    parser.add_argument("--max-results", type=int, default=100)
    parser.add_argument("--timeout", type=int, default=45)
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()

    module_path = ROOT / args.module if not Path(args.module).is_absolute() else Path(args.module)
    out = ROOT / args.out if not Path(args.out).is_absolute() else Path(args.out)
    summary_out = ROOT / args.summary_out if not Path(args.summary_out).is_absolute() else Path(args.summary_out)
    out.parent.mkdir(parents=True, exist_ok=True)
    summary_out.parent.mkdir(parents=True, exist_ok=True)

    if not module_path.exists():
        results = _pending(f"module table missing: {module_path}")
    elif args.offline:
        results = _pending("offline mode requested; Enrichr LINCS query not executed")
    else:
        module = pd.read_csv(module_path, sep="\t")
        try:
            results = query_lincs_reversal(
                module,
                top_genes=args.top_genes,
                max_results=args.max_results,
                timeout=args.timeout,
            )
            if results.empty:
                results = _pending("no Enrichr LINCS hits returned for module signature")
        except (requests.RequestException, ValueError, KeyError) as exc:
            results = _pending(f"Enrichr LINCS query failed: {exc}")

    results.to_csv(out, sep="\t", index=False)
    summary = aggregate_reversal_candidates(results[results["status"] != "RESULT_PENDING"])
    if summary.empty:
        summary = pd.DataFrame(
            [
                {
                    "perturbation_id": "RESULT_PENDING",
                    "perturbation_name": "RESULT_PENDING",
                    "target_gene": "RESULT_PENDING",
                    "mechanism": "RESULT_PENDING",
                    "reversal_score": "",
                    "depmap_score": "",
                    "dgidb_evidence": "",
                    "target_state": "RESULT_PENDING",
                    "priority_score": "",
                    "status": "RESULT_PENDING",
                    "interpretation": str(results.iloc[0]["interpretation"]) if not results.empty else "no results",
                }
            ]
        )
    summary.to_csv(summary_out, sep="\t", index=False)
    print(f"Wrote {out}")
    print(f"Wrote {summary_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
