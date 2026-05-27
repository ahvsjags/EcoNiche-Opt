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

from econiche_opt.model.endpoint_modules import MODULE_GENE_SETS, MODULE_PRIOR_WEIGHTS


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
    target = thresholds[thresholds["model_name"] == "EcoNiche-Opt-HeuristicEcology-LockedPanel"].copy()
    keep = [
        "endpoint",
        "threshold",
        "training_n",
        "training_responders",
        "training_nonresponders",
        "training_AUROC",
        "training_balanced_accuracy",
        "training_cohorts",
    ]
    return target[[col for col in keep if col in target.columns]].to_dict(orient="records")


def write_package(out_dir: Path, threshold_path: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    panel = _panel_gene_table()
    panel.to_csv(out_dir / "locked_panel_genes.tsv", sep="\t", index=False)

    sample_template = pd.DataFrame(
        [
            {
                "site_id": "SITE001",
                "subject_id": "SUBJECT001",
                "sample_id": "SUBJECT001_BL_RNA",
                "collection_timepoint": "baseline_pre_treatment",
                "collection_date": "YYYY-MM-DD",
                "cancer_type": "melanoma",
                "therapy": "anti-PD-1 monotherapy",
                "line_of_therapy": "",
                "assay_platform": "qPCR_or_NanoString",
                "rna_input_ng": "",
                "tumor_content_percent": "",
                "qc_pass": "",
                "locked_validation_use_flag": "",
            }
        ]
    )
    sample_template.to_csv(out_dir / "assay_sample_manifest_template.tsv", sep="\t", index=False)

    clinical_template = pd.DataFrame(
        [
            {
                "subject_id": "SUBJECT001",
                "sample_id": "SUBJECT001_BL_RNA",
                "response_raw": "CR/PR/SD/PD",
                "best_overall_response_date": "YYYY-MM-DD",
                "progression_date": "YYYY-MM-DD_or_blank",
                "last_follow_up_date": "YYYY-MM-DD",
                "dcbr_6mo": "DCB/NDB_or_blank",
                "recist_version": "RECIST 1.1",
                "label_source_document": "source_file_or_eCRF_id",
                "curation_notes": "",
            }
        ]
    )
    clinical_template.to_csv(out_dir / "clinical_annotation_template.tsv", sep="\t", index=False)

    scoring_spec = {
        "model_name": "EcoNiche-Opt-HeuristicEcology-LockedPanel",
        "model_status": "frozen_for_external_validation",
        "training_threshold_source": "GSE91061,GSE78220,PRJEB23709_PD1_PRE only",
        "input_transform": "within-sample expression matrix, gene-level values normalized by the registered pipeline; module score is the mean available gene score per module",
        "score_formula": "sigmoid(sum_q module_weight_q * module_activity_q)",
        "module_weights": MODULE_PRIOR_WEIGHTS,
        "endpoint_thresholds": _threshold_spec(threshold_path),
        "label_encoding": {"responder": 1, "nonresponder": 0},
        "allowed_primary_claim": "locked prospective validation performance only after independently collected samples are scored with this frozen specification",
        "forbidden_claim": "do not claim prospective validation from retrospective public cohorts",
    }
    (out_dir / "locked_scoring_spec.json").write_text(json.dumps(scoring_spec, indent=2), encoding="utf-8")

    protocol = """# Prospective Locked Validation Protocol

## Objective

Validate the frozen EcoNiche-Opt-HeuristicEcology-LockedPanel score in an independent baseline pretreatment melanoma anti-PD-1 cohort using a qPCR or NanoString-compatible panel.

## Locked Model

The model, module gene list, module weights, endpoint definitions, and thresholds are frozen before new clinical samples are scored. No feature selection, coefficient refitting, threshold tuning, calibration fitting, or endpoint relabeling is allowed on the prospective validation cohort.

## Primary Endpoint

Primary RECIST analysis: CR/PR are responders and SD/PD are nonresponders. Sensitivity analyses use strict RECIST (CR/PR vs PD, SD excluded) and clinical benefit (CR/PR/SD vs PD or DCB vs NDB where prospectively specified).

## Inclusion Criteria

Pretreatment tumor RNA sample, documented ICB regimen, patient-level response annotation, and sufficient assay QC. Baseline samples must be collected before first dose or before the relevant ICB cycle specified in the protocol.

## Statistical Analysis

Report AUROC, AUPRC, balanced accuracy at the locked threshold, sensitivity, specificity, PPV, NPV, Brier score, ECE, calibration slope/intercept, and decision-curve net benefit. Compare against IFNG, CXCL9, TIG, TIDE_dysfunction, APM, CYT, IPRES, and TIDE_exclusion using paired bootstrap or DeLong where available with Benjamini-Hochberg FDR correction.

## Leakage Guard

The validation cohort must not be used for module selection, hyperparameter tuning, threshold selection, calibration, or manuscript claim wording before the analysis is locked. Any excluded sample must retain an auditable exclusion reason.
"""
    (out_dir / "prospective_validation_protocol.md").write_text(protocol, encoding="utf-8")

    sap = """# Statistical Analysis Plan

1. Freeze the `locked_scoring_spec.json` file and hash it before sample scoring.
2. Score every QC-passing baseline sample exactly once.
3. Apply the endpoint-specific locked threshold from the discovery cohorts.
4. Run the primary RECIST analysis first, then strict RECIST and clinical benefit sensitivity analyses.
5. Use paired bootstrap or DeLong comparisons with FDR correction for superiority claims.
6. Use family-level omnibus claims only when the predeclared eight-signature family test is FDR-supported.
7. Report all failed or missing assay genes through `locked_panel_genes.tsv` coverage fields; do not impute outcome labels.
"""
    (out_dir / "statistical_analysis_plan.md").write_text(sap, encoding="utf-8")

    readme = """# Prospective and Clinical-Assay Validation Package

This folder freezes the EcoNiche-Opt panel and analysis plan for future independent qPCR/NanoString validation. It complements, but does not replace, the retrospective locked external and NanoString transfer analyses in `results/locked_external_panel_validation`.

Generated artifacts:

- `locked_scoring_spec.json`: frozen score formula, endpoint thresholds, and claim boundary.
- `locked_panel_genes.tsv`: qPCR/NanoString panel gene list with module weights.
- `assay_sample_manifest_template.tsv`: assay/sample metadata template.
- `clinical_annotation_template.tsv`: response-label curation template.
- `prospective_validation_protocol.md`: validation protocol.
- `statistical_analysis_plan.md`: predeclared analysis plan.

No prospective performance claim is allowed until new independently collected clinical samples are scored with this frozen package.
"""
    (out_dir / "README.md").write_text(readme, encoding="utf-8")
    print(f"Wrote prospective validation package to {out_dir}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--thresholds", default="results/locked_external_panel_validation/locked_thresholds.tsv")
    parser.add_argument("--out", default="deliverables/prospective_validation")
    args = parser.parse_args()
    write_package(ROOT / args.out, ROOT / args.thresholds)


if __name__ == "__main__":
    main()
