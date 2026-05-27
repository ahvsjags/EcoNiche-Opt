# Strict External Failure-Mode Audit

This audit evaluates MAP4K1-TBX3/AXL transform blends and reports whether the strict external gap is driven by one external cohort or by both cohorts.
External labels are not used for the primary-selected or fixed robust-development candidates. Rows marked as stress screens are diagnostic only.

- `primary_auc_selected_blend`: blend=0.50*cohort_gene_percentile+0.50*cohort_zscore; primary AUROC=0.789, BA=0.744; strict external AUROC=0.677, BA=0.626; per-cohort GSE145996:0.825;PHS000452_LIU_LIKE_PRE:0.647; boundary=selected_by_primary_lodo_only_not_by_external
- `robust_fixed_development_candidate`: blend=0.05*cohort_zscore+0.95*cohort_robust_zscore; primary AUROC=0.772, BA=0.656; strict external AUROC=0.686, BA=0.624; per-cohort GSE145996:0.875;PHS000452_LIU_LIKE_PRE:0.652; boundary=fixed_robust_transform_candidate_no_external_label_fit
- `primary_pass_external_stress_best`: blend=0.05*cohort_zscore+0.95*cohort_robust_zscore; primary AUROC=0.772, BA=0.656; strict external AUROC=0.686, BA=0.624; per-cohort GSE145996:0.875;PHS000452_LIU_LIKE_PRE:0.652; boundary=current_external_stress_screen_not_a_locked_selection_claim
- `current_external_stress_best`: blend=cohort_robust_zscore; primary AUROC=0.771, BA=0.644; strict external AUROC=0.686, BA=0.624; per-cohort GSE145996:0.875;PHS000452_LIU_LIKE_PRE:0.652; boundary=current_external_stress_screen_not_a_locked_selection_claim

## Interpretation

The strongest development candidates show high AUROC in GSE145996 but substantially lower AUROC in PHS000452_LIU_LIKE_PRE, making the Liu/MGSP-like cohort the main limiter of the strict external AUROC target.
