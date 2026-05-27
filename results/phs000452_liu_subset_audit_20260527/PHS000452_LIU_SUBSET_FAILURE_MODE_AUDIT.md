# PHS000452/Liu Subset Failure-Mode Audit

This audit checks whether the strict external gap reflects response-label discordance, source processing differences, or identifiable metadata subgroups in the Liu/MGSP-like PHS000452 cohort.

## Source Concordance

- TIGER unique patients=121, cBioPortal unique patients=121, matched=120, response mismatches=0; missing_from_cbio=Patient192; missing_from_tiger=Patient20.

## Cohort Contrast

- `current_external_stress_best` / `GSE145996`: AUROC=0.875, AUPRC=0.921, BA=0.738, n=13.
- `current_external_stress_best` / `PHS000452_LIU_LIKE_PRE`: AUROC=0.652, AUPRC=0.638, BA=0.609, n=103.
- `primary_auc_selected_blend` / `GSE145996`: AUROC=0.825, AUPRC=0.864, BA=0.738, n=13.
- `primary_auc_selected_blend` / `PHS000452_LIU_LIKE_PRE`: AUROC=0.647, AUPRC=0.604, BA=0.611, n=103.
- `robust_fixed_development_candidate` / `GSE145996`: AUROC=0.875, AUPRC=0.921, BA=0.738, n=13.
- `robust_fixed_development_candidate` / `PHS000452_LIU_LIKE_PRE`: AUROC=0.652, AUPRC=0.638, BA=0.609, n=103.

## PHS000452 Subgroups

- `0.05*cohort_zscore+0.95*cohort_robust_zscore` / `sex::Male`: AUROC=0.703, AUPRC=0.649, BA=0.696, n=61, responses=CR:11;PD:33;PR:17; boundary=demographic_subgroup_diagnostic.
- `0.05*cohort_zscore+0.95*cohort_robust_zscore` / `patient_id_suffix::T_M`: AUROC=0.695, AUPRC=0.607, BA=0.659, n=31, responses=CR:2;PD:20;PR:9; boundary=metadata_suffix_diagnostic.
- `0.05*cohort_zscore+0.95*cohort_robust_zscore` / `vital_status::Dead`: AUROC=0.695, AUPRC=0.354, BA=0.683, n=53, responses=CR:1;PD:47;PR:5; boundary=post_outcome_diagnostic_not_predictive_claim.
- `0.05*cohort_zscore+0.95*cohort_robust_zscore` / `m_stage::M1c`: AUROC=0.683, AUPRC=0.702, BA=0.634, n=79, responses=CR:13;PD:40;PR:26; boundary=clinical_stage_diagnostic.
- `0.05*cohort_zscore+0.95*cohort_robust_zscore` / `patient_id_suffix::T_P`: AUROC=0.672, AUPRC=0.705, BA=0.625, n=28, responses=CR:4;PD:16;PR:8; boundary=metadata_suffix_diagnostic.
- `0.05*cohort_zscore+0.95*cohort_robust_zscore` / `all_phs_strict`: AUROC=0.652, AUPRC=0.638, BA=0.609, n=103, responses=CR:16;PD:55;PR:32; boundary=primary_strict_external_reference.
- `0.05*cohort_zscore+0.95*cohort_robust_zscore` / `m_stage::M1b`: AUROC=0.639, AUPRC=0.665, BA=0.583, n=12, responses=CR:1;PD:6;PR:5; boundary=clinical_stage_diagnostic.
- `0.05*cohort_zscore+0.95*cohort_robust_zscore` / `patient_id_suffix::unsuffixed`: AUROC=0.615, AUPRC=0.651, BA=0.583, n=44, responses=CR:10;PD:19;PR:15; boundary=metadata_suffix_diagnostic.
- `0.05*cohort_zscore+0.95*cohort_robust_zscore` / `sex::Female`: AUROC=0.575, AUPRC=0.632, BA=0.482, n=42, responses=CR:5;PD:22;PR:15; boundary=demographic_subgroup_diagnostic.
- `0.05*cohort_zscore+0.95*cohort_robust_zscore` / `vital_status::Alive`: AUROC=0.500, AUPRC=0.856, BA=0.533, n=50, responses=CR:15;PD:8;PR:27; boundary=post_outcome_diagnostic_not_predictive_claim.
- `0.50*cohort_gene_percentile+0.50*cohort_zscore` / `sex::Male`: AUROC=0.705, AUPRC=0.632, BA=0.701, n=61, responses=CR:11;PD:33;PR:17; boundary=demographic_subgroup_diagnostic.
- `0.50*cohort_gene_percentile+0.50*cohort_zscore` / `vital_status::Dead`: AUROC=0.699, AUPRC=0.290, BA=0.661, n=53, responses=CR:1;PD:47;PR:5; boundary=post_outcome_diagnostic_not_predictive_claim.
- `0.50*cohort_gene_percentile+0.50*cohort_zscore` / `patient_id_suffix::T_M`: AUROC=0.677, AUPRC=0.528, BA=0.659, n=31, responses=CR:2;PD:20;PR:9; boundary=metadata_suffix_diagnostic.
- `0.50*cohort_gene_percentile+0.50*cohort_zscore` / `m_stage::M1c`: AUROC=0.674, AUPRC=0.670, BA=0.622, n=79, responses=CR:13;PD:40;PR:26; boundary=clinical_stage_diagnostic.
- `0.50*cohort_gene_percentile+0.50*cohort_zscore` / `patient_id_suffix::T_P`: AUROC=0.661, AUPRC=0.647, BA=0.594, n=28, responses=CR:4;PD:16;PR:8; boundary=metadata_suffix_diagnostic.
- `0.50*cohort_gene_percentile+0.50*cohort_zscore` / `all_phs_strict`: AUROC=0.647, AUPRC=0.604, BA=0.611, n=103, responses=CR:16;PD:55;PR:32; boundary=primary_strict_external_reference.
- `0.50*cohort_gene_percentile+0.50*cohort_zscore` / `patient_id_suffix::unsuffixed`: AUROC=0.629, AUPRC=0.674, BA=0.597, n=44, responses=CR:10;PD:19;PR:15; boundary=metadata_suffix_diagnostic.
- `0.50*cohort_gene_percentile+0.50*cohort_zscore` / `m_stage::M1b`: AUROC=0.611, AUPRC=0.569, BA=0.667, n=12, responses=CR:1;PD:6;PR:5; boundary=clinical_stage_diagnostic.
- `0.50*cohort_gene_percentile+0.50*cohort_zscore` / `sex::Female`: AUROC=0.550, AUPRC=0.557, BA=0.482, n=42, responses=CR:5;PD:22;PR:15; boundary=demographic_subgroup_diagnostic.
- `0.50*cohort_gene_percentile+0.50*cohort_zscore` / `vital_status::Alive`: AUROC=0.470, AUPRC=0.820, BA=0.557, n=50, responses=CR:15;PD:8;PR:27; boundary=post_outcome_diagnostic_not_predictive_claim.
- `cohort_robust_zscore` / `sex::Male`: AUROC=0.703, AUPRC=0.649, BA=0.696, n=61, responses=CR:11;PD:33;PR:17; boundary=demographic_subgroup_diagnostic.
- `cohort_robust_zscore` / `patient_id_suffix::T_M`: AUROC=0.695, AUPRC=0.607, BA=0.659, n=31, responses=CR:2;PD:20;PR:9; boundary=metadata_suffix_diagnostic.
- `cohort_robust_zscore` / `vital_status::Dead`: AUROC=0.695, AUPRC=0.354, BA=0.683, n=53, responses=CR:1;PD:47;PR:5; boundary=post_outcome_diagnostic_not_predictive_claim.
- `cohort_robust_zscore` / `m_stage::M1c`: AUROC=0.683, AUPRC=0.703, BA=0.634, n=79, responses=CR:13;PD:40;PR:26; boundary=clinical_stage_diagnostic.
- `cohort_robust_zscore` / `patient_id_suffix::T_P`: AUROC=0.677, AUPRC=0.711, BA=0.625, n=28, responses=CR:4;PD:16;PR:8; boundary=metadata_suffix_diagnostic.
- `cohort_robust_zscore` / `all_phs_strict`: AUROC=0.652, AUPRC=0.638, BA=0.609, n=103, responses=CR:16;PD:55;PR:32; boundary=primary_strict_external_reference.
- `cohort_robust_zscore` / `m_stage::M1b`: AUROC=0.639, AUPRC=0.665, BA=0.583, n=12, responses=CR:1;PD:6;PR:5; boundary=clinical_stage_diagnostic.
- `cohort_robust_zscore` / `patient_id_suffix::unsuffixed`: AUROC=0.615, AUPRC=0.651, BA=0.583, n=44, responses=CR:10;PD:19;PR:15; boundary=metadata_suffix_diagnostic.
- `cohort_robust_zscore` / `sex::Female`: AUROC=0.573, AUPRC=0.632, BA=0.482, n=42, responses=CR:5;PD:22;PR:15; boundary=demographic_subgroup_diagnostic.
- `cohort_robust_zscore` / `vital_status::Alive`: AUROC=0.494, AUPRC=0.855, BA=0.533, n=50, responses=CR:15;PD:8;PR:27; boundary=post_outcome_diagnostic_not_predictive_claim.

## Interpretation

TIGER and cBioPortal labels are concordant for matched Liu/DFCI patients. The robust MAP4K1-TBX3/AXL candidate performs well in GSE145996 but is weaker in PHS000452 overall, with metadata subgroup heterogeneity. These subgroup rows are diagnostic and should not be used to redefine the locked strict external claim.
