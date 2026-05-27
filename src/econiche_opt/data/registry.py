from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from econiche.registry import audit_accession as _audit_accession
from econiche.registry import load_registry as _load_registry
from econiche.registry import normalize_access_status
from econiche.registry import validate_registry as _validate_registry
from econiche.registry import write_registry_report

MINIMUM_ACCESSIONS = {
    "GSE91061",
    "GSE78220",
    "GSE145996",
    "GSE168204",
    "GSE115821",
    "GSE93157",
    "GSE244982",
    "GSE123728",
    "GSE165745",
    "GSE122220",
    "Liu_DFCI_s41591_2019",
    "EGAS00001001552",
    "Gide_PRJEB23709",
    "GSE136961",
    "GSE176307",
    "IMvigor210_EGAS00001002556",
    "GSE67501",
    "GSE121810",
    "GSE140901",
    "GSE165252",
    "GSE183924",
    "GSE115978",
    "GSE123139",
    "TCGA_SKCM_Xena",
    "GDC_TCGA_SKCM",
    "LINCS_L1000",
    "CMap",
    "DepMap",
    "DGIdb",
    "DrugBank_optional",
}

ACCESSION_ALIASES = {
    "PRJEB23709": "Gide_PRJEB23709",
    "Gide_PRJEB23709": "PRJEB23709",
    "Liu_DFCI_melanoma": "Liu_DFCI_s41591_2019",
    "Liu_DFCI_s41591_2019": "Liu_DFCI_melanoma",
    "IMvigor210": "IMvigor210_EGAS00001002556",
    "IMvigor210_EGAS00001002556": "IMvigor210",
}

RECOMMENDED_FIELDS = {
    "accession",
    "cancer_type",
    "therapy",
    "platform",
    "timepoints",
    "endpoint",
    "role",
    "access",
    "download_script",
    "preprocessing_script",
    "uses",
}


def load_registry(path: str | Path) -> dict[str, Any]:
    return _load_registry(path)


def _present_accessions(registry: dict[str, Any]) -> set[str]:
    present = {str(row.get("accession", "")) for row in registry.get("cohorts", [])}
    expanded = set(present)
    for accession in present:
        if accession in ACCESSION_ALIASES:
            expanded.add(ACCESSION_ALIASES[accession])
    return {item for item in expanded if item}


def validate_registry(
    registry: dict[str, Any],
    require_minimum: bool = True,
) -> pd.DataFrame:
    base = _validate_registry(registry).copy()
    rows: list[dict[str, Any]] = []
    for cohort in registry.get("cohorts", []):
        accession = str(cohort.get("accession", "UNKNOWN"))
        missing_recommended = sorted(
            field for field in RECOMMENDED_FIELDS if field not in cohort or cohort.get(field) in (None, "", [])
        )
        rows.append(
            {
                "accession": accession,
                "is_valid": not missing_recommended,
                "missing_fields": ",".join(missing_recommended),
                "warning": "missing_recommended_fields" if missing_recommended else "",
            }
        )
    detailed = pd.DataFrame(rows)
    if not base.empty and not detailed.empty:
        detailed = detailed.drop(columns=["is_valid"], errors="ignore").merge(
            base[["accession", "is_valid"]],
            on="accession",
            how="left",
        )
        detailed["is_valid"] = detailed["is_valid"].fillna(False) & (detailed["missing_fields"] == "")

    if require_minimum:
        present = _present_accessions(registry)
        missing = sorted(MINIMUM_ACCESSIONS - present)
        if missing:
            missing_rows = pd.DataFrame(
                {
                    "accession": missing,
                    "is_valid": False,
                    "missing_fields": "registry_entry",
                    "warning": "missing_minimum_accession",
                }
            )
            detailed = pd.concat([detailed, missing_rows], ignore_index=True)
    return detailed


def audit_accession(registry: dict[str, Any]) -> pd.DataFrame:
    audit = _audit_accession(registry).copy()
    if audit.empty:
        return audit
    audit["normalized_access_status"] = audit["access"].map(normalize_access_status)
    audit["download_action"] = audit["normalized_access_status"].map(
        {
            "public": "download_or_metadata_extract",
            "controlled": "emit_access_instructions",
            "unknown": "verify_before_use",
        }
    )
    return audit


def write_audit(registry_path: str | Path, out: str | Path) -> pd.DataFrame:
    registry = load_registry(registry_path)
    return write_registry_report(registry, out)
