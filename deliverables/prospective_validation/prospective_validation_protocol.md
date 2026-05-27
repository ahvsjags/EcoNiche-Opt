# Prospective Locked Validation Protocol

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
