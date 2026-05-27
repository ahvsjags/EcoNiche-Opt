# Signed-rank Module Audit

This audit tests whether sample-wise rank Gaussian module scoring and training-only gene-direction signing improve the current raw module-prior score. External cohorts are used only for evaluation.

## Summary

- primary_recist / melanoma_core_high_evidence / EcoNiche-Opt-RankModulePrior: AUROC=0.693, AUPRC=0.557, balanced accuracy=0.538, ECE=0.237.
- primary_recist / melanoma_core_high_evidence / EcoNiche-Opt-RankModulePriorCalibrated: AUROC=0.640, AUPRC=0.485, balanced accuracy=0.622, ECE=0.187.
- primary_recist / melanoma_core_high_evidence / EcoNiche-Opt-RawModulePrior: AUROC=0.705, AUPRC=0.558, balanced accuracy=0.632, ECE=0.235.
- primary_recist / melanoma_core_high_evidence / EcoNiche-Opt-SignedRankResponseModules: AUROC=0.580, AUPRC=0.438, balanced accuracy=0.455, ECE=0.277.
- primary_recist / melanoma_core_high_evidence / EcoNiche-Opt-SignedRankResponseModulesCalibrated: AUROC=0.550, AUPRC=0.419, balanced accuracy=0.518, ECE=0.190.
- primary_recist / melanoma_recist_supported_primary / EcoNiche-Opt-RankModulePrior: AUROC=0.683, AUPRC=0.578, balanced accuracy=0.550, ECE=0.248.
- primary_recist / melanoma_recist_supported_primary / EcoNiche-Opt-RankModulePriorCalibrated: AUROC=0.650, AUPRC=0.529, balanced accuracy=0.609, ECE=0.157.
- primary_recist / melanoma_recist_supported_primary / EcoNiche-Opt-RawModulePrior: AUROC=0.685, AUPRC=0.580, balanced accuracy=0.613, ECE=0.239.
- primary_recist / melanoma_recist_supported_primary / EcoNiche-Opt-SignedRankResponseModules: AUROC=0.584, AUPRC=0.477, balanced accuracy=0.497, ECE=0.260.
- primary_recist / melanoma_recist_supported_primary / EcoNiche-Opt-SignedRankResponseModulesCalibrated: AUROC=0.575, AUPRC=0.481, balanced accuracy=0.555, ECE=0.142.
- strict_recist / strict_cbio_liu_plus_gse145996 / EcoNiche-Opt-RankModulePrior: AUROC=0.547, AUPRC=0.530, balanced accuracy=0.527, ECE=0.394.
- strict_recist / strict_cbio_liu_plus_gse145996 / EcoNiche-Opt-RankModulePriorCalibrated: AUROC=0.547, AUPRC=0.530, balanced accuracy=0.509, ECE=0.131.
- strict_recist / strict_cbio_liu_plus_gse145996 / EcoNiche-Opt-RawModulePrior: AUROC=0.556, AUPRC=0.592, balanced accuracy=0.552, ECE=0.209.
- strict_recist / strict_cbio_liu_plus_gse145996 / EcoNiche-Opt-SignedRankResponseModules: AUROC=0.559, AUPRC=0.536, balanced accuracy=0.503, ECE=0.267.
- strict_recist / strict_cbio_liu_plus_gse145996 / EcoNiche-Opt-SignedRankResponseModulesCalibrated: AUROC=0.559, AUPRC=0.536, balanced accuracy=0.501, ECE=0.144.
- strict_recist / strict_melanoma_pd1_like_external / EcoNiche-Opt-RankModulePrior: AUROC=0.546, AUPRC=0.539, balanced accuracy=0.504, ECE=0.371.
- strict_recist / strict_melanoma_pd1_like_external / EcoNiche-Opt-RankModulePriorCalibrated: AUROC=0.546, AUPRC=0.539, balanced accuracy=0.490, ECE=0.151.
- strict_recist / strict_melanoma_pd1_like_external / EcoNiche-Opt-RawModulePrior: AUROC=0.571, AUPRC=0.600, balanced accuracy=0.531, ECE=0.261.
- strict_recist / strict_melanoma_pd1_like_external / EcoNiche-Opt-SignedRankResponseModules: AUROC=0.566, AUPRC=0.546, balanced accuracy=0.532, ECE=0.223.
- strict_recist / strict_melanoma_pd1_like_external / EcoNiche-Opt-SignedRankResponseModulesCalibrated: AUROC=0.566, AUPRC=0.546, balanced accuracy=0.564, ECE=0.115.

## Baseline Comparisons

- primary_recist / melanoma_core_high_evidence / EcoNiche-Opt-RankModulePrior: dAUROC=-0.012, dBA=-0.094, dECE=0.002, q=0.596 (not_supported).
- primary_recist / melanoma_core_high_evidence / EcoNiche-Opt-RankModulePriorCalibrated: dAUROC=-0.064, dBA=-0.010, dECE=-0.047, q=0.016 (point_estimate_or_calibration_tradeoff).
- primary_recist / melanoma_core_high_evidence / EcoNiche-Opt-SignedRankResponseModules: dAUROC=-0.125, dBA=-0.177, dECE=0.042, q=0.004 (not_supported).
- primary_recist / melanoma_core_high_evidence / EcoNiche-Opt-SignedRankResponseModulesCalibrated: dAUROC=-0.155, dBA=-0.113, dECE=-0.044, q=0.000 (point_estimate_or_calibration_tradeoff).
- primary_recist / melanoma_recist_supported_primary / EcoNiche-Opt-RankModulePrior: dAUROC=-0.001, dBA=-0.063, dECE=0.009, q=0.946 (not_supported).
- primary_recist / melanoma_recist_supported_primary / EcoNiche-Opt-RankModulePriorCalibrated: dAUROC=-0.034, dBA=-0.004, dECE=-0.082, q=0.211 (point_estimate_or_calibration_tradeoff).
- primary_recist / melanoma_recist_supported_primary / EcoNiche-Opt-SignedRankResponseModules: dAUROC=-0.101, dBA=-0.116, dECE=0.021, q=0.000 (not_supported).
- primary_recist / melanoma_recist_supported_primary / EcoNiche-Opt-SignedRankResponseModulesCalibrated: dAUROC=-0.110, dBA=-0.058, dECE=-0.096, q=0.000 (point_estimate_or_calibration_tradeoff).
- strict_recist / strict_cbio_liu_plus_gse145996 / EcoNiche-Opt-RankModulePrior: dAUROC=-0.009, dBA=-0.024, dECE=0.185, q=0.942 (not_supported).
- strict_recist / strict_cbio_liu_plus_gse145996 / EcoNiche-Opt-RankModulePriorCalibrated: dAUROC=-0.009, dBA=-0.043, dECE=-0.078, q=0.942 (point_estimate_or_calibration_tradeoff).
- strict_recist / strict_cbio_liu_plus_gse145996 / EcoNiche-Opt-SignedRankResponseModules: dAUROC=0.002, dBA=-0.049, dECE=0.057, q=0.942 (point_estimate_or_calibration_tradeoff).
- strict_recist / strict_cbio_liu_plus_gse145996 / EcoNiche-Opt-SignedRankResponseModulesCalibrated: dAUROC=0.002, dBA=-0.051, dECE=-0.066, q=0.942 (point_estimate_or_calibration_tradeoff).
- strict_recist / strict_melanoma_pd1_like_external / EcoNiche-Opt-RankModulePrior: dAUROC=-0.025, dBA=-0.027, dECE=0.110, q=0.644 (not_supported).
- strict_recist / strict_melanoma_pd1_like_external / EcoNiche-Opt-RankModulePriorCalibrated: dAUROC=-0.025, dBA=-0.040, dECE=-0.110, q=0.644 (point_estimate_or_calibration_tradeoff).
- strict_recist / strict_melanoma_pd1_like_external / EcoNiche-Opt-SignedRankResponseModules: dAUROC=-0.005, dBA=0.001, dECE=-0.039, q=0.812 (point_estimate_or_calibration_tradeoff).
- strict_recist / strict_melanoma_pd1_like_external / EcoNiche-Opt-SignedRankResponseModulesCalibrated: dAUROC=-0.005, dBA=0.033, dECE=-0.146, q=0.812 (point_estimate_or_calibration_tradeoff).