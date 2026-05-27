from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.registry import normalize_access_status
from econiche_opt.data.registry import load_registry


STRICT_EXTERNAL_ROLE_TOKENS = ("external_candidate", "melanoma_primary", "high_value_restricted_external_candidate")


def _as_list(value: object) -> list[str]:
    if isinstance(value, list):
        return [str(item) for item in value]
    if value in (None, ""):
        return []
    return [str(value)]


def _controlled_targets(registry_path: Path, lead_triage_path: Path | None = None) -> pd.DataFrame:
    registry = load_registry(registry_path)
    lead_triage = pd.DataFrame()
    if lead_triage_path is not None and lead_triage_path.exists():
        lead_triage = pd.read_csv(lead_triage_path, sep="\t")
    triage_by_accession: dict[str, pd.Series] = {}
    if not lead_triage.empty:
        for _, row in lead_triage.iterrows():
            for token in str(row.get("candidate_accession", "")).replace(";", ",").split(","):
                token = token.strip()
                if token:
                    triage_by_accession[token] = row
                    triage_by_accession[token.split(".")[0]] = row

    rows: list[dict[str, object]] = []
    for cohort in registry.get("cohorts", []):
        access_status = normalize_access_status(cohort.get("access"))
        role = str(cohort.get("role", ""))
        uses = _as_list(cohort.get("uses"))
        role_text = " ".join([role, *uses])
        if access_status != "controlled":
            continue
        if not any(token in role_text for token in STRICT_EXTERNAL_ROLE_TOKENS):
            continue
        accession = str(cohort.get("accession", "UNKNOWN"))
        triage = triage_by_accession.get(accession)
        if triage is None:
            triage = triage_by_accession.get(accession.split(".")[0])
        rows.append(
            {
                "accession": accession,
                "name": cohort.get("name", accession),
                "source_database": _source_database(accession, cohort.get("download_script", "")),
                "request_url": cohort.get("download_script", ""),
                "access_status": cohort.get("access", ""),
                "cancer_type": cohort.get("cancer_type", ""),
                "therapy": cohort.get("therapy", ""),
                "platform": cohort.get("platform", ""),
                "timepoints": ",".join(_as_list(cohort.get("timepoints"))),
                "endpoint": ",".join(_as_list(cohort.get("endpoint"))),
                "intended_role": role,
                "uses": ",".join(uses),
                "lead_id": "" if triage is None else triage.get("lead_id", ""),
                "lead_title": "" if triage is None else triage.get("title", ""),
                "local_status": "ACCESS_RESTRICTED_NOT_VALIDATION_EVIDENCE",
                "required_after_access": "expression_matrix,sample_manifest,clinical_annotation,response_provenance,baseline_timepoint_evidence",
                "locked_scoring_boundary": "No cohort labels may enter training, feature selection, thresholding, calibration, or model selection.",
                "import_expression_path": f"data/raw/controlled/{accession}/expression.tsv",
                "import_clinical_path": f"data/raw/controlled/{accession}/clinical_annotation.tsv",
                "import_manifest_path": f"data/raw/controlled/{accession}/assay_sample_manifest.tsv",
                "notes": cohort.get("notes", ""),
            }
        )
    return pd.DataFrame(rows).sort_values(["source_database", "accession"]).reset_index(drop=True)


def _source_database(accession: str, url: object) -> str:
    text = f"{accession} {url}".lower()
    if "dbgap" in text or accession.lower().startswith("phs"):
        return "dbGaP"
    if "ega" in text or accession.lower().startswith("egas"):
        return "EGA"
    return "controlled_repository"


def _write_templates(out: Path, targets: pd.DataFrame) -> None:
    template_cols = [
        "sample_id",
        "patient_id",
        "cohort",
        "assay_sample_id",
        "timepoint",
        "treatment",
        "response_raw",
        "response_binary",
        "recist_version",
        "response_assessment_date",
        "biopsy_date",
        "source_document",
        "curation_note",
    ]
    pd.DataFrame(columns=template_cols).to_csv(out / "controlled_clinical_annotation_template.tsv", sep="\t", index=False)
    manifest_cols = [
        "sample_id",
        "patient_id",
        "cohort",
        "raw_file",
        "library_strategy",
        "platform",
        "tumor_or_blood",
        "baseline_status",
        "normalization_status",
        "md5_or_checksum",
    ]
    pd.DataFrame(columns=manifest_cols).to_csv(out / "controlled_assay_sample_manifest_template.tsv", sep="\t", index=False)
    expr_note = (
        "Expression import format: rows are sample_id, columns are HGNC gene symbols, values are normalized expression "
        "or raw counts with a documented normalization path. Do not include labels in the expression matrix.\n"
    )
    (out / "controlled_expression_matrix_format.txt").write_text(expr_note, encoding="utf-8")

    for source, frame in targets.groupby("source_database"):
        filename = f"{source.lower()}_request_checklist.md"
        title = f"{source} Access Request Checklist"
        lines = [
            f"# {title}",
            "",
            "Purpose: obtain controlled melanoma checkpoint-blockade tumor RNA-seq and phenotype files for locked external validation of EcoNiche-Opt.",
            "",
            "Required safeguards:",
            "",
            "- Use approved institutional access only.",
            "- Store controlled files outside Git and public archives.",
            "- Import only through the documented schema in this package.",
            "- Keep the cohort fully locked: no training, feature selection, thresholding, calibration, or model selection.",
            "",
            "Targets:",
            "",
        ]
        for _, row in frame.iterrows():
            lines.append(f"- `{row['accession']}`: {row['name']} ({row['request_url']})")
        lines.extend(
            [
                "",
                "After access is approved:",
                "",
                "1. Place raw controlled files under `data/raw/controlled/<accession>/` locally only.",
                "2. Fill `controlled_assay_sample_manifest_template.tsv` and `controlled_clinical_annotation_template.tsv`.",
                "3. Run the registered preprocessing/import path and then the locked external scoring command.",
                "4. Record missing labels or missing baseline evidence as `RESULT_PENDING`, not as negative or imputed validation evidence.",
            ]
        )
        (out / filename).write_text("\n".join(lines), encoding="utf-8")


def _write_readme(out: Path, targets: pd.DataFrame) -> None:
    lines = [
        "# Controlled External Validation Access Package",
        "",
        "This package supports the strict melanoma external-validation target for EcoNiche-Opt.",
        "It does not contain controlled data and does not count any controlled cohort as validation evidence.",
        "",
        "## Why This Package Exists",
        "",
        "The current public strict melanoma external AUROC is below the top-tier target. The realistic path to a stronger locked external test is to obtain independent controlled melanoma checkpoint-blockade tumor RNA-seq cohorts, then score them without refitting.",
        "",
        "## Locked Boundary",
        "",
        "Controlled external cohort labels must never enter model training, feature selection, threshold selection, calibration, or candidate selection. They may only be used after the model, genes, weights, thresholds, and scoring script are frozen.",
        "",
        "## Targets",
        "",
    ]
    for _, row in targets.iterrows():
        lines.append(f"- `{row['accession']}` ({row['source_database']}): {row['name']}; request: {row['request_url']}")
    lines.extend(
        [
            "",
            "## Expected Local Import Files",
            "",
            "- `expression.tsv`: sample-by-gene expression matrix.",
            "- `assay_sample_manifest.tsv`: file/sample/patient/platform/timing manifest.",
            "- `clinical_annotation.tsv`: patient-level response, treatment, timing, and provenance table.",
            "",
            "## Locked Scoring Command",
            "",
            "```bash",
            "python -m econiche_opt.cli score-locked-validation \\",
            "  --package-dir deliverables/prospective_validation \\",
            "  --expression data/raw/controlled/<accession>/expression.tsv \\",
            "  --sample-manifest data/raw/controlled/<accession>/assay_sample_manifest.tsv \\",
            "  --clinical-annotation data/raw/controlled/<accession>/clinical_annotation.tsv \\",
            "  --out-dir results/controlled_external_validation/<accession>",
            "```",
            "",
            "## Claim Rule",
            "",
            "If expression, baseline timing, or response provenance is incomplete, report `RESULT_PENDING`. Do not impute controlled external validation.",
            "",
        ]
    )
    (out / "README.md").write_text("\n".join(lines), encoding="utf-8")


def build_package(registry: Path, lead_triage: Path, out: Path) -> pd.DataFrame:
    targets = _controlled_targets(registry, lead_triage)
    out.mkdir(parents=True, exist_ok=True)
    targets.to_csv(out / "controlled_external_access_targets.tsv", sep="\t", index=False)
    _write_templates(out, targets)
    _write_readme(out, targets)
    return targets


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", default="config/data_registry.yml")
    parser.add_argument("--lead-triage", default="deliverables/public_external_lead_triage_20260527.tsv")
    parser.add_argument("--out", default="deliverables/controlled_external_access_package_20260527")
    args = parser.parse_args()

    targets = build_package(ROOT / args.registry, ROOT / args.lead_triage, ROOT / args.out)
    print(json.dumps({"n_targets": int(len(targets)), "sources": sorted(targets["source_database"].unique())}, ensure_ascii=False))
    print(f"Wrote {ROOT / args.out}")


if __name__ == "__main__":
    main()
