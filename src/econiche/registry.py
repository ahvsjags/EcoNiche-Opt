from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import yaml


REQUIRED_FIELDS = [
    "accession",
    "layer",
    "cancer_type",
    "therapy",
    "platform",
    "timepoints",
    "endpoint",
    "role",
    "access",
]


def load_registry(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if isinstance(data, list):
        return {"cohorts": data}
    if "cohorts" not in data:
        data["cohorts"] = []
    return data


def validate_registry(registry: dict[str, Any]) -> pd.DataFrame:
    rows = []
    for cohort in registry.get("cohorts", []):
        missing = [field for field in REQUIRED_FIELDS if field not in cohort or cohort.get(field) in (None, "")]
        rows.append(
            {
                "accession": cohort.get("accession", "UNKNOWN"),
                "is_valid": len(missing) == 0,
                "missing_fields": ",".join(missing),
                "warning": "missing_required_fields" if missing else "",
            }
        )
    return pd.DataFrame(rows)


def normalize_access_status(access: str | None) -> str:
    text = (access or "").lower()
    if text == "public" or text.startswith("public_"):
        return "public"
    if any(token in text for token in ["controlled", "dbgap", "ega", "restricted"]):
        return "controlled"
    return "unknown"


def audit_accession(registry: dict[str, Any]) -> pd.DataFrame:
    rows = []
    validation = validate_registry(registry).set_index("accession") if registry.get("cohorts") else pd.DataFrame()
    for cohort in registry.get("cohorts", []):
        accession = cohort.get("accession", "UNKNOWN")
        status = normalize_access_status(cohort.get("access"))
        rows.append(
            {
                "accession": accession,
                "name": cohort.get("name", accession),
                "layer": cohort.get("layer"),
                "cancer_type": cohort.get("cancer_type"),
                "therapy": cohort.get("therapy"),
                "platform": cohort.get("platform"),
                "role": cohort.get("role"),
                "access": cohort.get("access"),
                "access_status": status,
                "priority": cohort.get("priority", "medium"),
                "warning": "ACCESS_RESTRICTED"
                if status == "controlled"
                else ("ACCESS_UNKNOWN_VERIFY_BEFORE_USE" if status == "unknown" else ""),
                "is_valid": bool(validation.loc[accession, "is_valid"]) if accession in validation.index else False,
            }
        )
    return pd.DataFrame(rows)


def write_registry_report(registry: dict[str, Any], out: str | Path) -> pd.DataFrame:
    out_path = Path(out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    audit = audit_accession(registry)
    audit.to_csv(out_path, sep="\t", index=False)
    roles = audit[["accession", "name", "layer", "cancer_type", "therapy", "role", "access_status", "priority"]]
    roles.to_csv(out_path.with_name("dataset_roles.tsv"), sep="\t", index=False)
    validation = validate_registry(registry)
    validation.to_csv(out_path.with_name("dataset_registry_validation.tsv"), sep="\t", index=False)
    return audit
