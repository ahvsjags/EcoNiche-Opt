from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

REQUIRED_FILES = [
    "README.md",
    "locked_scoring_spec.json",
    "locked_scoring_spec.sha256",
    "locked_panel_genes.tsv",
    "assay_sample_manifest_template.tsv",
    "clinical_annotation_template.tsv",
    "prospective_validation_protocol.md",
    "statistical_analysis_plan.md",
    "clinical_partner_intake_checklist.md",
    "tumor_tissue_assay_sop.md",
    "validation_readiness_checklist.tsv",
]

REQUIRED_SAMPLE_COLUMNS = {
    "site_id",
    "subject_id",
    "sample_id",
    "sample_source",
    "baseline_status",
    "therapy",
    "assay_platform",
    "tumor_content_percent",
    "panel_gene_coverage_percent",
    "qc_pass",
    "locked_validation_use_flag",
    "exclusion_reason",
    "data_freeze_id",
}

REQUIRED_CLINICAL_COLUMNS = {
    "subject_id",
    "sample_id",
    "response_raw",
    "recist_version",
    "primary_recist_label",
    "strict_recist_label",
    "clinical_benefit_label",
    "label_source_document",
    "source_page_or_record_id",
}

EXPECTED_ENDPOINTS = {"primary_recist", "strict_recist", "clinical_benefit"}
EXPECTED_COMPARATORS = {"IFNG", "CXCL9", "TIG", "TIDE_dysfunction", "APM", "CYT", "IPRES", "TIDE_exclusion"}


def _row(check: str, ok: bool, detail: str = "") -> dict[str, object]:
    return {"check": check, "is_valid": bool(ok), "detail": detail}


def validate_package(package_dir: str | Path) -> pd.DataFrame:
    root = Path(package_dir)
    rows: list[dict[str, object]] = []

    for name in REQUIRED_FILES:
        rows.append(_row(f"exists:{name}", (root / name).exists()))

    panel_path = root / "locked_panel_genes.tsv"
    if panel_path.exists():
        panel = pd.read_csv(panel_path, sep="\t")
        rows.append(_row("panel_unique_gene_count", int(panel["gene_symbol"].nunique()) == 62, str(panel["gene_symbol"].nunique())))
        rows.append(_row("panel_gene_module_rows", len(panel) >= 62, str(len(panel))))
        rows.append(_row("panel_required_columns", {"module", "gene_symbol", "module_weight", "score_direction", "locked_in_model"}.issubset(panel.columns)))

    sample_path = root / "assay_sample_manifest_template.tsv"
    if sample_path.exists():
        sample = pd.read_csv(sample_path, sep="\t")
        missing = sorted(REQUIRED_SAMPLE_COLUMNS - set(sample.columns))
        rows.append(_row("sample_manifest_required_columns", not missing, ",".join(missing)))

    clinical_path = root / "clinical_annotation_template.tsv"
    if clinical_path.exists():
        clinical = pd.read_csv(clinical_path, sep="\t")
        missing = sorted(REQUIRED_CLINICAL_COLUMNS - set(clinical.columns))
        rows.append(_row("clinical_annotation_required_columns", not missing, ",".join(missing)))

    spec_path = root / "locked_scoring_spec.json"
    if spec_path.exists():
        spec = json.loads(spec_path.read_text(encoding="utf-8"))
        endpoints = {str(item.get("endpoint")) for item in spec.get("endpoint_thresholds", [])}
        comparators = set(spec.get("predeclared_comparators", []))
        rows.append(_row("spec_model_status", spec.get("model_status") == "frozen_for_external_validation", str(spec.get("model_status"))))
        rows.append(_row("spec_release_commit_present", bool(spec.get("release_commit")), str(spec.get("release_commit", ""))[:12]))
        rows.append(_row("spec_panel_unique_genes", spec.get("panel_unique_genes") == 62, str(spec.get("panel_unique_genes"))))
        rows.append(_row("spec_endpoints", endpoints == EXPECTED_ENDPOINTS, ",".join(sorted(endpoints))))
        rows.append(_row("spec_comparators", comparators == EXPECTED_COMPARATORS, ",".join(sorted(comparators))))
        missing_calibration = [
            str(item.get("endpoint"))
            for item in spec.get("endpoint_thresholds", [])
            if not {"calibration_coef", "calibration_intercept"}.issubset(item)
        ]
        rows.append(_row("spec_calibration_parameters", not missing_calibration, ",".join(missing_calibration)))

    hash_path = root / "locked_scoring_spec.sha256"
    if spec_path.exists() and hash_path.exists():
        expected = hash_path.read_text(encoding="utf-8").split()[0]
        observed = hashlib.sha256(spec_path.read_bytes()).hexdigest()
        rows.append(_row("spec_sha256", observed == expected, observed))

    readiness_path = root / "validation_readiness_checklist.tsv"
    if readiness_path.exists():
        readiness = pd.read_csv(readiness_path, sep="\t")
        rows.append(_row("readiness_required_rows", len(readiness) >= 7, str(len(readiness))))
        rows.append(_row("readiness_required_columns", {"section", "item", "required", "evidence_field", "status"}.issubset(readiness.columns)))

    text_blob = ""
    for name in ["README.md", "prospective_validation_protocol.md", "tumor_tissue_assay_sop.md"]:
        path = root / name
        if path.exists():
            text_blob += path.read_text(encoding="utf-8").lower() + "\n"
    rows.append(_row("tumor_tissue_context_documented", "pretreatment melanoma tumor tissue" in text_blob or "pretreatment melanoma tumor-tissue" in text_blob))
    rows.append(_row("blood_boundary_documented", "blood" in text_blob and "separately trained" in text_blob))
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package-dir", default="deliverables/prospective_validation")
    args = parser.parse_args()
    report = validate_package(args.package_dir)
    print(report.to_string(index=False))
    if not bool(report["is_valid"].all()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
