from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.registry import normalize_access_status
from econiche_opt.data.registry import load_registry


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--demo", action="store_true", help="Accepted for CLI symmetry; smoke tests use real/public registry entries.")
    parser.add_argument("--registry", default="config/data_registry.yml")
    parser.add_argument("--results-dir", default="results/real")
    parser.add_argument("--raw-dir", default="data/raw")
    parser.add_argument("--out", default="results/audit/public_cohort_smoke_tests.tsv")
    args = parser.parse_args()

    registry = load_registry(ROOT / args.registry)
    results_dir = ROOT / args.results_dir
    raw_dir = ROOT / args.raw_dir
    metrics_path = results_dir / "lodo_metrics.tsv"
    prediction_path = results_dir / "lodo_predictions.tsv"
    metrics = pd.read_csv(metrics_path, sep="\t") if metrics_path.exists() else pd.DataFrame()
    prediction_cohorts = set()
    if prediction_path.exists():
        prediction_cohorts = set(pd.read_csv(prediction_path, sep="\t", usecols=["cohort"])["cohort"].astype(str))
    metric_cohorts = set(metrics.get("cohort", pd.Series(dtype=str)).astype(str))
    artifact_map = {
        "TCGA_SKCM_Xena": [ROOT / "results/real/xena_tcga_manifest.tsv", raw_dir / "TCGA_SKCM_Xena"],
        "GDC_TCGA_SKCM": [ROOT / "results/real/xena_tcga_manifest.tsv", raw_dir / "TCGA_SKCM_Xena/GDC_TCGA_SKCM_STAR_counts_manifest.tsv"],
        "LINCS_L1000": [ROOT / "results/perturbation/lincs_reversal.tsv"],
        "CMap": [ROOT / "results/perturbation/lincs_reversal.tsv"],
        "DepMap": [ROOT / "results/perturbation/depmap_targets.tsv", raw_dir / "DepMap"],
        "DGIdb": [ROOT / "results/perturbation/dgidb_hits.tsv"],
    }

    rows = []
    for cohort in registry.get("cohorts", []):
        accession = str(cohort.get("accession", "UNKNOWN"))
        access = normalize_access_status(cohort.get("access"))
        raw_exists = (raw_dir / accession).exists()
        artifact_exists = any(path.exists() for path in artifact_map.get(accession, []))
        has_metrics = accession in metric_cohorts
        has_predictions = accession in prediction_cohorts
        if access == "controlled":
            status = "ACCESS_RESTRICTED"
            reason = "controlled_or_licensed_access"
        elif has_metrics and has_predictions:
            status = "PASS"
            reason = "real_result_available"
        elif artifact_exists:
            status = "PASS"
            reason = "registered_artifact_available"
        elif raw_exists:
            status = "RESULT_PENDING"
            reason = "raw_or_metadata_available_but_not_in_current_result_table"
        elif access == "public":
            status = "RESULT_PENDING"
            reason = "public_entry_registered_but_not_downloaded_or_curated"
        else:
            status = "VERIFY_ACCESS"
            reason = "registry_access_status_requires_manual_review"
        rows.append(
            {
                "accession": accession,
                "access_status": access,
                "raw_exists": raw_exists,
                "artifact_exists": artifact_exists,
                "has_metrics": has_metrics,
                "has_predictions": has_predictions,
                "status": status,
                "reason": reason,
            }
        )

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
