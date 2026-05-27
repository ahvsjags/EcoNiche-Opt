# cBioPortal Rescue-Head External Validation Audit

This audit checks whether the frozen MAP4K1-TBX3/AXL rescue head can be fairly recomputed on cBioPortal Liu/DFCI after explicitly requesting all three target genes.
Primary-model selection and thresholding use only discovery cohorts. Rows marked as stress screens are diagnostic and cannot be used as locked external model selection.

## Target-Gene Coverage

- CBIO_IATLAS_LIU_2019_PRE: 3/3 target genes (AXL,MAP4K1,TBX3); status=ready.
- CBIO_LIU_DFCI_2019_PRE: 3/3 target genes (AXL,MAP4K1,TBX3); status=ready.
- GSE145996: 3/3 target genes (AXL,MAP4K1,TBX3); status=ready.
- GSE78220: 3/3 target genes (AXL,MAP4K1,TBX3); status=ready.
- GSE91061: 3/3 target genes (AXL,MAP4K1,TBX3); status=ready.
- PRJEB23709_PD1_PRE: 3/3 target genes (AXL,MAP4K1,TBX3); status=ready.

## Selection Summary

- `primary_auc_selected_blend` / `cbio_iatlas_liu_duplicate_crosscheck`: blend=0.50*cohort_gene_percentile+0.50*cohort_zscore; primary AUROC=0.789; strict external AUROC=0.627; BA=0.551; per-cohort CBIO_IATLAS_LIU_2019_PRE:0.627; boundary=selected_by_primary_lodo_only_not_by_cbio_external
- `robust_fixed_development_candidate` / `cbio_iatlas_liu_duplicate_crosscheck`: blend=0.05*cohort_zscore+0.95*cohort_robust_zscore; primary AUROC=0.772; strict external AUROC=0.628; BA=0.542; per-cohort CBIO_IATLAS_LIU_2019_PRE:0.628; boundary=fixed_robust_transform_candidate_no_cbio_external_label_fit
- `current_cbio_external_stress_best` / `cbio_iatlas_liu_duplicate_crosscheck`: blend=cohort_zscore; primary AUROC=0.779; strict external AUROC=0.631; BA=0.548; per-cohort CBIO_IATLAS_LIU_2019_PRE:0.631; boundary=diagnostic_current_cbio_external_stress_screen_not_selection_claim
- `primary_auc_selected_blend` / `cbio_liu_dfci_only`: blend=0.50*cohort_gene_percentile+0.50*cohort_zscore; primary AUROC=0.789; strict external AUROC=0.609; BA=0.590; per-cohort CBIO_LIU_DFCI_2019_PRE:0.609; boundary=selected_by_primary_lodo_only_not_by_cbio_external
- `robust_fixed_development_candidate` / `cbio_liu_dfci_only`: blend=0.05*cohort_zscore+0.95*cohort_robust_zscore; primary AUROC=0.772; strict external AUROC=0.638; BA=0.552; per-cohort CBIO_LIU_DFCI_2019_PRE:0.638; boundary=fixed_robust_transform_candidate_no_cbio_external_label_fit
- `current_cbio_external_stress_best` / `cbio_liu_dfci_only`: blend=0.05*cohort_zscore+0.95*cohort_robust_zscore; primary AUROC=0.772; strict external AUROC=0.638; BA=0.552; per-cohort CBIO_LIU_DFCI_2019_PRE:0.638; boundary=diagnostic_current_cbio_external_stress_screen_not_selection_claim
- `primary_auc_selected_blend` / `strict_cbio_liu_plus_gse145996`: blend=0.50*cohort_gene_percentile+0.50*cohort_zscore; primary AUROC=0.789; strict external AUROC=0.646; BA=0.605; per-cohort CBIO_LIU_DFCI_2019_PRE:0.609;GSE145996:0.825; boundary=selected_by_primary_lodo_only_not_by_cbio_external
- `robust_fixed_development_candidate` / `strict_cbio_liu_plus_gse145996`: blend=0.05*cohort_zscore+0.95*cohort_robust_zscore; primary AUROC=0.772; strict external AUROC=0.661; BA=0.568; per-cohort CBIO_LIU_DFCI_2019_PRE:0.638;GSE145996:0.875; boundary=fixed_robust_transform_candidate_no_cbio_external_label_fit
- `current_cbio_external_stress_best` / `strict_cbio_liu_plus_gse145996`: blend=0.50*cohort_zscore+0.50*cohort_robust_zscore; primary AUROC=0.777; strict external AUROC=0.666; BA=0.567; per-cohort CBIO_LIU_DFCI_2019_PRE:0.634;GSE145996:0.875; boundary=diagnostic_current_cbio_external_stress_screen_not_selection_claim