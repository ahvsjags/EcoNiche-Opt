# Statistical Analysis Plan

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
