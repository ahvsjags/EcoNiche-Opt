# Threshold Recalibration Audit

This audit evaluates threshold policies without changing the rank ordering of EcoNiche-Opt scores.
For LODO rows, each holdout threshold is selected from training cohorts only.
The nested midrange policy selects among fixed 0.40, 0.50 and 0.60 using inner validation on the outer-training cohorts.

- `melanoma_core_high_evidence` `original_fold_threshold`: BA=0.613, AUROC=0.705, AUPRC=0.558, threshold_median=0.436, boundary=outer_lodo_training_threshold_from_source_pipeline
- `melanoma_core_high_evidence` `fixed_0.40`: BA=0.688, AUROC=0.705, AUPRC=0.558, threshold_median=0.400, boundary=outer_lodo_training_cohorts_only
- `melanoma_core_high_evidence` `fixed_0.50`: BA=0.632, AUROC=0.705, AUPRC=0.558, threshold_median=0.500, boundary=outer_lodo_training_cohorts_only
- `melanoma_core_high_evidence` `fixed_0.60`: BA=0.654, AUROC=0.705, AUPRC=0.558, threshold_median=0.600, boundary=outer_lodo_training_cohorts_only
- `melanoma_core_high_evidence` `training_youden`: BA=0.613, AUROC=0.705, AUPRC=0.558, threshold_median=0.424, boundary=outer_lodo_training_cohorts_only
- `melanoma_core_high_evidence` `training_prevalence_quantile`: BA=0.600, AUROC=0.705, AUPRC=0.558, threshold_median=0.722, boundary=outer_lodo_training_cohorts_only
- `melanoma_core_high_evidence` `nested_midrange_fixed_grid`: BA=0.688, AUROC=0.705, AUPRC=0.558, threshold_median=0.400, boundary=outer_lodo_training_cohorts_only
- `melanoma_recist_supported_primary` `original_fold_threshold`: BA=0.594, AUROC=0.685, AUPRC=0.580, threshold_median=0.421, boundary=outer_lodo_training_threshold_from_source_pipeline
- `melanoma_recist_supported_primary` `fixed_0.40`: BA=0.647, AUROC=0.685, AUPRC=0.580, threshold_median=0.400, boundary=outer_lodo_training_cohorts_only
- `melanoma_recist_supported_primary` `fixed_0.50`: BA=0.613, AUROC=0.685, AUPRC=0.580, threshold_median=0.500, boundary=outer_lodo_training_cohorts_only
- `melanoma_recist_supported_primary` `fixed_0.60`: BA=0.631, AUROC=0.685, AUPRC=0.580, threshold_median=0.600, boundary=outer_lodo_training_cohorts_only
- `melanoma_recist_supported_primary` `training_youden`: BA=0.548, AUROC=0.685, AUPRC=0.580, threshold_median=0.414, boundary=outer_lodo_training_cohorts_only
- `melanoma_recist_supported_primary` `training_prevalence_quantile`: BA=0.590, AUROC=0.685, AUPRC=0.580, threshold_median=0.621, boundary=outer_lodo_training_cohorts_only
- `melanoma_recist_supported_primary` `nested_midrange_fixed_grid`: BA=0.622, AUROC=0.685, AUPRC=0.580, threshold_median=0.600, boundary=outer_lodo_training_cohorts_only