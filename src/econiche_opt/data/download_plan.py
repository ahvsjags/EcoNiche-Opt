from __future__ import annotations

from pathlib import Path

import pandas as pd

from econiche_opt.data.registry import load_registry
from econiche.registry import normalize_access_status


def build_download_plan(registry_path: str | Path) -> pd.DataFrame:
    registry = load_registry(registry_path)
    rows = []
    for cohort in registry.get("cohorts", []):
        accession = cohort.get("accession", "UNKNOWN")
        access = normalize_access_status(cohort.get("access"))
        if access == "public":
            status = "READY"
            action = cohort.get("download_script", "UNSPECIFIED")
        elif access == "controlled":
            status = "ACCESS_RESTRICTED"
            action = "record_instructions_only"
        else:
            status = "VERIFY_ACCESS"
            action = "manual_access_review"
        rows.append(
            {
                "accession": accession,
                "name": cohort.get("name", accession),
                "access_status": access,
                "planned_status": status,
                "download_script": cohort.get("download_script", ""),
                "preprocessing_script": cohort.get("preprocessing_script", ""),
                "action": action,
            }
        )
    return pd.DataFrame(rows)


def write_download_plan(registry_path: str | Path, out: str | Path) -> pd.DataFrame:
    plan = build_download_plan(registry_path)
    out_path = Path(out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plan.to_csv(out_path, sep="\t", index=False)
    return plan
