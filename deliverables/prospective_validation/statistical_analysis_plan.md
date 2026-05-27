# Statistical Analysis Plan

1. Freeze the `locked_scoring_spec.json` file and hash it before sample scoring.
2. Score every QC-passing baseline sample exactly once.
3. Apply the endpoint-specific locked threshold from the discovery cohorts.
4. Run the primary RECIST analysis first, then strict RECIST and clinical benefit sensitivity analyses.
5. Use paired bootstrap or DeLong comparisons with FDR correction for superiority claims.
6. Use family-level omnibus claims only when the predeclared eight-signature family test is FDR-supported.
7. Report all failed or missing assay genes through `locked_panel_genes.tsv` coverage fields; do not impute outcome labels.
