from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche_opt.model.endpoint_modules import MODULE_GENE_SETS, MODULE_PRIOR_WEIGHTS

RELEASE_TAG = "v0.3.1-jtm-20260527"
RELEASE_COMMIT = "4371e3ac3d70133f4635bf83a5165dcc2fe4357c"
MODEL_NAME = "EcoNiche-Opt-HeuristicEcology-LockedPanel"
PRIMARY_CONTEXT = "independent pretreatment melanoma tumor-tissue anti-PD-1/anti-PD-1-based cohort"
COMPARATORS = ["IFNG", "CXCL9", "TIG", "TIDE_dysfunction", "APM", "CYT", "IPRES", "TIDE_exclusion"]


def _panel_gene_table() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for module, genes in MODULE_GENE_SETS.items():
        for gene in genes:
            rows.append(
                {
                    "module": module,
                    "gene_symbol": gene,
                    "module_weight": MODULE_PRIOR_WEIGHTS.get(module, 0.0),
                    "score_direction": "response_high" if MODULE_PRIOR_WEIGHTS.get(module, 0.0) >= 0 else "nonresponse_high",
                    "locked_in_model": True,
                }
            )
    return pd.DataFrame(rows)


def _threshold_spec(threshold_path: Path) -> list[dict[str, object]]:
    if not threshold_path.exists():
        return []
    thresholds = pd.read_csv(threshold_path, sep="\t")
    target = thresholds[thresholds["model_name"] == MODEL_NAME].copy()
    keep = [
        "endpoint",
        "threshold",
        "calibration",
        "calibration_coef",
        "calibration_intercept",
        "training_n",
        "training_responders",
        "training_nonresponders",
        "training_AUROC",
        "training_balanced_accuracy",
        "training_cohorts",
    ]
    return target[[col for col in keep if col in target.columns]].to_dict(orient="records")


def _sample_template() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "site_id": "SITE001",
                "subject_id": "SUBJECT001",
                "sample_id": "SUBJECT001_BL_RNA",
                "patient_id_at_site": "LOCAL001",
                "specimen_id": "PATH001-BLK1",
                "collection_timepoint": "baseline_pre_treatment",
                "collection_date": "YYYY-MM-DD",
                "days_before_icb_start": "",
                "cancer_type": "melanoma",
                "disease_stage": "unresectable_stage_III_or_IV",
                "sample_source": "tumor_tissue",
                "specimen_type": "FFPE_or_fresh_frozen",
                "anatomic_site": "skin/lymph_node/visceral/other",
                "baseline_status": "pretreatment_before_first_ICB_dose",
                "therapy": "anti-PD-1_or_anti-PD-1-based",
                "line_of_therapy": "",
                "icb_start_date": "YYYY-MM-DD",
                "assay_platform": "qPCR_or_NanoString_or_RNAseq",
                "rna_input_ng": "",
                "rna_dv200_or_rin": "",
                "tumor_content_percent": "",
                "necrosis_percent": "",
                "macrodissection_performed": "",
                "housekeeping_qc_pass": "",
                "panel_gene_coverage_percent": "",
                "qc_pass": "",
                "locked_validation_use_flag": "",
                "exclusion_reason": "",
                "data_freeze_id": "FREEZE001",
            }
        ]
    )


def _clinical_template() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "subject_id": "SUBJECT001",
                "sample_id": "SUBJECT001_BL_RNA",
                "icb_regimen": "anti-PD-1_or_anti-PD-1-based",
                "icb_start_date": "YYYY-MM-DD",
                "baseline_scan_date": "YYYY-MM-DD",
                "response_raw": "CR/PR/SD/PD",
                "best_overall_response_date": "YYYY-MM-DD",
                "progression_date": "YYYY-MM-DD_or_blank",
                "last_follow_up_date": "YYYY-MM-DD",
                "dcbr_6mo": "DCB/NDB_or_blank",
                "recist_version": "RECIST 1.1",
                "primary_recist_label": "1_for_CR_PR__0_for_SD_PD",
                "strict_recist_label": "1_for_CR_PR__0_for_PD__blank_for_SD",
                "clinical_benefit_label": "1_for_CR_PR_SD_or_DCB__0_for_PD_or_NDB",
                "response_assessor": "radiologist_or_clinician_id",
                "label_source_document": "source_file_or_eCRF_id",
                "source_page_or_record_id": "",
                "curation_notes": "",
            }
        ]
    )


def _scoring_spec(panel: pd.DataFrame, threshold_path: Path) -> dict[str, object]:
    return {
        "model_name": MODEL_NAME,
        "model_status": "frozen_for_external_validation",
        "release_tag": RELEASE_TAG,
        "release_commit": RELEASE_COMMIT,
        "intended_validation_context": PRIMARY_CONTEXT,
        "unsupported_context_without_new_training": "blood, plasma, PBMC, serum, on-treatment biopsy, and non-ICB treatment-only cohorts",
        "training_threshold_source": "GSE91061,GSE78220,PRJEB23709_PD1_PRE only",
        "input_matrix_format": "samples-by-genes matrix with sample_id as row identifier and HGNC gene symbols as columns",
        "input_transform": "assay values must be normalized within the registered platform workflow; EcoNiche-Opt then computes signed rank-gene and module scores using the locked gene list and response/resistance directions",
        "score_formula": "locked logistic response probability from signed immune-ecology module activities plus ecological interaction features; endpoint thresholds are discovery-only calibrated and frozen before validation scoring",
        "panel_unique_genes": int(panel["gene_symbol"].nunique()),
        "panel_gene_module_rows": int(len(panel)),
        "module_weights": MODULE_PRIOR_WEIGHTS,
        "endpoint_thresholds": _threshold_spec(threshold_path),
        "predeclared_comparators": COMPARATORS,
        "required_primary_report_metrics": [
            "AUROC",
            "AUPRC",
            "balanced_accuracy",
            "sensitivity",
            "specificity",
            "PPV",
            "NPV",
            "Brier",
            "ECE",
            "calibration_intercept",
            "calibration_slope",
            "decision_curve_net_benefit",
        ],
        "label_encoding": {"responder": 1, "nonresponder": 0},
        "allowed_primary_claim": "locked independent validation performance only after newly assembled samples are scored with this frozen specification",
        "forbidden_claim": "do not use validation labels for feature selection, thresholding, calibration, refitting, endpoint relabeling, or claim wording before the statistical analysis plan is frozen",
    }


PROTOCOL = """# Prospective Locked Validation Protocol

## Objective

Validate the frozen EcoNiche-Opt-HeuristicEcology-LockedPanel score in an independent pretreatment melanoma tumor-tissue cohort treated with anti-PD-1 or anti-PD-1-based immune-checkpoint blockade. The preferred assay is the locked 62-unique-gene qPCR/NanoString-compatible panel; RNA-seq is acceptable when all locked genes can be quantified and the registered normalization workflow is documented.

## Locked Model

The model, module gene list, gene directions, module weights, endpoint definitions, calibration method, and thresholds are frozen before new clinical samples are scored. No feature selection, coefficient refitting, threshold tuning, calibration fitting, endpoint relabeling, assay-gene substitution, or cohort-specific model selection is allowed on the independent validation cohort.

## Primary Endpoint

Primary RECIST analysis: CR/PR are responders and SD/PD are nonresponders. Sensitivity analyses use strict RECIST (CR/PR vs PD, SD excluded) and clinical benefit (CR/PR/SD vs PD or DCB vs NDB where prospectively specified).

## Inclusion Criteria

- Histologically confirmed melanoma.
- Pretreatment tumor tissue collected before the first anti-PD-1/anti-PD-1-based dose.
- Tumor RNA available from FFPE or fresh-frozen tissue with documented pathology review.
- Patient-level RECIST 1.1 best overall response or prospectively defined DCB/NDB endpoint.
- Sample, patient, therapy, scan, response, and expression identifiers are traceable through the manifest and clinical annotation templates.
- Adequate assay QC, tumor content, and locked-panel coverage.

## Exclusion Criteria

- Blood, plasma, serum, PBMC, or other non-tumor specimens, unless a separate blood-specific EcoNiche-Opt model is trained and validated.
- On-treatment biopsies when the primary validation question is pretreatment prediction.
- Missing response evidence, ambiguous subject/sample matching, failed RNA/assay QC, or duplicate samples without a predeclared patient-level deduplication rule.
- Non-ICB treatment-only cohorts.

## Statistical Analysis

Report AUROC, AUPRC, balanced accuracy at the locked threshold, sensitivity, specificity, PPV, NPV, Brier score, ECE, calibration slope/intercept, and decision-curve net benefit. Compare against IFNG, CXCL9, TIG, TIDE_dysfunction, APM, CYT, IPRES, and TIDE_exclusion using paired bootstrap or DeLong where available with Benjamini-Hochberg FDR correction.

The primary analysis is patient-level. When multiple pretreatment specimens exist for the same patient, the protocol must choose one representative specimen before response labels are viewed or use a predeclared aggregation rule. A target validation size of 50-100 patients is recommended for a credible first independent tumor-tissue estimate; smaller cohorts remain feasibility or pilot analyses.

## Leakage Guard

The validation cohort must not be used for module selection, hyperparameter tuning, threshold selection, calibration, or manuscript claim wording before the analysis is locked. Any excluded sample must retain an auditable exclusion reason.
"""


SAP = """# Statistical Analysis Plan

1. Freeze the `locked_scoring_spec.json` file and hash it before sample scoring.
2. Register the data freeze, site list, inclusion/exclusion decisions, and patient-level deduplication rule before outcome analysis.
3. Score every QC-passing pretreatment tumor-tissue sample exactly once.
4. Apply the endpoint-specific locked threshold from the discovery cohorts.
5. Run the primary RECIST analysis first, then strict RECIST and clinical benefit sensitivity analyses.
6. Use paired bootstrap or DeLong comparisons with FDR correction for superiority claims.
7. Use family-level omnibus claims only when the predeclared eight-signature family test is FDR-supported.
8. Report all failed or missing assay genes through `locked_panel_genes.tsv` coverage fields; do not impute outcome labels.
9. Report calibration and clinical-threshold behavior even when AUROC is favorable: Brier score, ECE, calibration intercept/slope, and decision-curve net benefit.
10. Mark any analysis that uses fewer than 50 independent patients as pilot validation, not definitive clinical validation.
11. Store all output tables, scoring logs, and exclusion reasons in the validation archive before manuscript claim drafting.
"""


INTAKE = """# Clinical Partner Intake Checklist

Use this checklist before accepting a new independent cohort for EcoNiche-Opt validation.

## Required Cohort Definition

- Disease: melanoma.
- Specimen: pretreatment tumor tissue, not blood, plasma, serum, or PBMC.
- Treatment: anti-PD-1 or anti-PD-1-based immune-checkpoint blockade.
- Endpoint evidence: RECIST 1.1 CR/PR/SD/PD with scan dates, or a prospectively defined DCB/NDB endpoint.
- Preferred sample size: 50-100 independent patients for a credible first validation; smaller sets are pilot/feasibility only.

## Files to Request From the Clinical Team

- De-identified sample manifest using `assay_sample_manifest_template.tsv`.
- De-identified clinical annotation using `clinical_annotation_template.tsv`.
- Expression matrix with samples as rows and HGNC gene symbols as columns.
- Pathology QC summary: tumor content, necrosis, macrodissection, FFPE/fresh-frozen status, RNA QC.
- Response evidence: source eCRF, tumor board export, RECIST worksheet, or annotated clinical Excel.
- Data dictionary explaining every non-standard field.

## Pre-Scoring Gate

Do not score the cohort until subject IDs, sample IDs, baseline status, therapy dates, assay platform, response labels, and exclusion decisions are complete. Validation labels must not be used for feature selection, calibration, threshold tuning, or model refitting.
"""


SOP = """# Tumor-Tissue Assay SOP for Independent EcoNiche-Opt Validation

## Specimen Requirement

EcoNiche-Opt is currently locked for pretreatment tumor-tissue transcriptomes. Blood-derived samples cannot directly validate this model because the score was trained and calibrated on tumor tissue and includes tumor-microenvironment states such as antigen presentation, myeloid suppression, stromal exclusion, TRM/TLS, and T/NK effector activity.

## Tissue and Pathology QC

1. Confirm melanoma diagnosis and sampling date relative to first ICB dose.
2. Prefer pretreatment FFPE or fresh-frozen tumor tissue with pathology review.
3. Record tumor content, necrosis percentage, macrodissection status, specimen block/slide ID, and anatomic site.
4. Extract RNA under the local certified workflow and record RIN or DV200, input amount, and assay batch.
5. Quantify the locked panel genes in `locked_panel_genes.tsv`; failed genes must be reported, not silently replaced.

## Expression Matrix

Rows must be sample IDs and columns must be HGNC gene symbols. Values should be normalized within the assay workflow and comparable across samples within the validation freeze. If the platform exports genes-by-samples, transpose before scoring and retain the raw export in the local audit archive.

## Locked Scoring

Use `locked_scoring_spec.json`, verify the SHA256 hash in `locked_scoring_spec.sha256`, score all QC-passing samples once, and apply the endpoint-specific locked threshold. Do not retrain or recalibrate on the independent validation labels.
"""


README = """# Prospective and Clinical-Assay Validation Package

This folder freezes the EcoNiche-Opt panel and analysis plan for independent qPCR/NanoString or RNA-seq validation in pretreatment melanoma tumor tissue. It complements, but does not replace, the retrospective locked external and NanoString transfer analyses in `results/locked_external_panel_validation_calibrated_20260519`.

Generated artifacts:

- `locked_scoring_spec.json`: frozen score formula, endpoint thresholds, and claim boundary.
- `locked_scoring_spec.sha256`: file hash to confirm the scoring spec has not changed.
- `locked_panel_genes.tsv`: 62-unique-gene qPCR/NanoString panel with gene-module rows, directions, and module weights.
- `assay_sample_manifest_template.tsv`: assay/sample metadata template.
- `clinical_annotation_template.tsv`: response-label curation template.
- `prospective_validation_protocol.md`: validation protocol.
- `statistical_analysis_plan.md`: predeclared analysis plan.
- `clinical_partner_intake_checklist.md`: what to request from a hospital or clinical collaborator.
- `tumor_tissue_assay_sop.md`: tumor-tissue processing and locked scoring SOP.
- `validation_readiness_checklist.tsv`: pre-scoring readiness gate.

The primary clinical scenario is pretreatment melanoma tumor tissue before anti-PD-1/anti-PD-1-based therapy. Blood, plasma, serum, PBMC, or on-treatment samples require a separately trained and validated model before clinical claims can be made.

## Locked Scoring Command

After the sample manifest, clinical annotation, and expression matrix pass the pre-scoring gate, score the independent cohort without retraining:

```bash
python -m econiche_opt.cli score-locked-validation \
  --package-dir deliverables/prospective_validation \
  --expression independent_expression.tsv \
  --sample-manifest assay_sample_manifest.tsv \
  --clinical-annotation clinical_annotation.tsv \
  --out-dir results/independent_locked_validation
```

The command writes sample-level endpoint probabilities, locked threshold calls, module scores, panel coverage, manifest audit rows, optional AUROC/calibration metrics when labels are supplied, and decision-curve outputs. It never performs feature selection, threshold tuning, calibration fitting, or model refitting on the independent cohort.
"""


def _readiness_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"section": "cohort", "item": "melanoma anti-PD-1/anti-PD-1-based cohort", "required": True, "evidence_field": "cancer_type;therapy", "status": "to_complete"},
            {"section": "specimen", "item": "pretreatment tumor tissue before first ICB dose", "required": True, "evidence_field": "sample_source;baseline_status;days_before_icb_start", "status": "to_complete"},
            {"section": "pathology", "item": "tumor content and RNA QC recorded", "required": True, "evidence_field": "tumor_content_percent;rna_dv200_or_rin;qc_pass", "status": "to_complete"},
            {"section": "labels", "item": "RECIST raw response and source evidence recorded", "required": True, "evidence_field": "response_raw;recist_version;label_source_document", "status": "to_complete"},
            {"section": "labels", "item": "primary/strict/clinical-benefit binary labels derived before scoring", "required": True, "evidence_field": "primary_recist_label;strict_recist_label;clinical_benefit_label", "status": "to_complete"},
            {"section": "assay", "item": "locked panel coverage reported", "required": True, "evidence_field": "panel_gene_coverage_percent", "status": "to_complete"},
            {"section": "leakage", "item": "no validation labels used for tuning or calibration", "required": True, "evidence_field": "data_freeze_id", "status": "to_complete"},
        ]
    )


def write_package(out_dir: Path, threshold_path: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    panel = _panel_gene_table()
    panel.to_csv(out_dir / "locked_panel_genes.tsv", sep="\t", index=False)
    _sample_template().to_csv(out_dir / "assay_sample_manifest_template.tsv", sep="\t", index=False)
    _clinical_template().to_csv(out_dir / "clinical_annotation_template.tsv", sep="\t", index=False)

    spec_path = out_dir / "locked_scoring_spec.json"
    spec_path.write_text(json.dumps(_scoring_spec(panel, threshold_path), indent=2), encoding="utf-8")
    (out_dir / "locked_scoring_spec.sha256").write_text(
        hashlib.sha256(spec_path.read_bytes()).hexdigest() + "  locked_scoring_spec.json\n",
        encoding="utf-8",
    )

    (out_dir / "prospective_validation_protocol.md").write_text(PROTOCOL, encoding="utf-8")
    (out_dir / "statistical_analysis_plan.md").write_text(SAP, encoding="utf-8")
    (out_dir / "clinical_partner_intake_checklist.md").write_text(INTAKE, encoding="utf-8")
    (out_dir / "tumor_tissue_assay_sop.md").write_text(SOP, encoding="utf-8")
    _readiness_table().to_csv(out_dir / "validation_readiness_checklist.tsv", sep="\t", index=False)
    (out_dir / "README.md").write_text(README, encoding="utf-8")
    print(f"Wrote prospective validation package to {out_dir}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--thresholds", default="results/locked_external_panel_validation_calibrated_20260519/locked_thresholds.tsv")
    parser.add_argument("--out", default="deliverables/prospective_validation")
    args = parser.parse_args()
    write_package(ROOT / args.out, ROOT / args.thresholds)


if __name__ == "__main__":
    main()
