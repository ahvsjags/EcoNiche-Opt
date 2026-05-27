from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche_opt.perturbation.reversal import aggregate_reversal_candidates, combine_with_dgidb  # noqa: E402


PRIORITY_COLUMNS = [
    "perturbation_id",
    "perturbation_name",
    "target_gene",
    "mechanism",
    "reversal_score",
    "depmap_score",
    "dgidb_evidence",
    "target_state",
    "priority_score",
    "status",
    "interpretation",
]


def _run_lincs(module_path: Path, lincs_path: Path) -> None:
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/perturbation/lincs_reversal.py"),
            "--module",
            str(module_path),
            "--out",
            str(lincs_path),
            "--summary-out",
            str(ROOT / "results/perturbation/prioritized_reversal.tsv"),
        ],
        cwd=ROOT,
        check=True,
    )


def _pending(reason: str) -> pd.DataFrame:
    return pd.DataFrame(
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
                "interpretation": reason,
            }
        ],
        columns=PRIORITY_COLUMNS,
    )


def _merge_depmap(priority: pd.DataFrame, depmap: pd.DataFrame) -> pd.DataFrame:
    if priority.empty or depmap.empty or "target_gene" not in depmap.columns or "depmap_score" not in depmap.columns:
        return priority
    depmap_scores = {
        str(row["target_gene"]).strip().upper(): pd.to_numeric(pd.Series([row["depmap_score"]]), errors="coerce").iloc[0]
        for _, row in depmap.iterrows()
    }
    updated = priority.copy()
    scores = []
    for genes in updated["target_gene"].astype(str):
        values = []
        for gene in genes.split(";"):
            value = depmap_scores.get(gene.strip().upper())
            if pd.notna(value):
                values.append(float(value))
        scores.append(sum(values) / len(values) if values else pd.NA)
    updated["depmap_score"] = scores
    numeric_priority = pd.to_numeric(updated["priority_score"], errors="coerce").fillna(0.0)
    numeric_depmap = pd.to_numeric(updated["depmap_score"], errors="coerce").fillna(0.0).clip(lower=0.0)
    updated["priority_score"] = numeric_priority + numeric_depmap * 10.0
    return updated.sort_values("priority_score", ascending=False)


def main() -> int:
    parser = argparse.ArgumentParser(description="Merge perturbation evidence into a ranked hypothesis table.")
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--out", default="results/perturbation/prioritized_perturbations.tsv")
    parser.add_argument("--lincs", default="results/perturbation/lincs_reversal.tsv")
    parser.add_argument("--dgidb", default="results/perturbation/dgidb_hits.tsv")
    parser.add_argument("--depmap", default="results/perturbation/depmap_targets.tsv")
    parser.add_argument("--refresh-lincs", action="store_true")
    args = parser.parse_args()

    module_path = ROOT / ("results/demo/econiche_module.tsv" if args.demo else "results/real/econiche_module.tsv")
    lincs_path = ROOT / args.lincs if not Path(args.lincs).is_absolute() else Path(args.lincs)
    out = ROOT / args.out if not Path(args.out).is_absolute() else Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    if args.refresh_lincs or not lincs_path.exists():
        _run_lincs(module_path, lincs_path)

    if not lincs_path.exists():
        result = _pending(f"LINCS reversal file missing: {lincs_path}")
    else:
        lincs = pd.read_csv(lincs_path, sep="\t")
        lincs = lincs[lincs.get("status", pd.Series(dtype=str)).astype(str) != "RESULT_PENDING"]
        result = aggregate_reversal_candidates(lincs)
        if result.empty:
            result = _pending("no scored LINCS reversal candidates available")

    dgidb_path = ROOT / args.dgidb if not Path(args.dgidb).is_absolute() else Path(args.dgidb)
    if dgidb_path.exists() and not result.empty and result.iloc[0]["status"] != "RESULT_PENDING":
        dgidb = pd.read_csv(dgidb_path, sep="\t")
        result = combine_with_dgidb(result, dgidb)

    depmap_path = ROOT / args.depmap if not Path(args.depmap).is_absolute() else Path(args.depmap)
    if depmap_path.exists() and not result.empty and result.iloc[0]["status"] != "RESULT_PENDING":
        depmap = pd.read_csv(depmap_path, sep="\t")
        depmap = depmap[depmap.get("status", pd.Series(dtype=str)).astype(str) != "RESULT_PENDING"]
        result = _merge_depmap(result, depmap)

    if args.demo and not result.empty:
        result["status"] = result["status"].replace({"hypothesis_only": "demo_hypothesis_only"})
    result.to_csv(out, sep="\t", index=False)
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
