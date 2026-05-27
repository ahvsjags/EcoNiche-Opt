# Prospective Locked Validation Protocol

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
