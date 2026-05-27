# Aligned Locked-Panel Ablation Audit

This audit tests ablation variants on the same locked module-panel scoring family that supports the primary EcoNiche-Opt melanoma result. It replaces the older WordFullGraph-only ablation as the primary component-evidence table because WordFullGraph is not the strongest discriminative model.

## Pooled Model Summary

- primary_recist / melanoma_core_high_evidence / EcoNiche-Opt-AlignedPanelNoCalibration: n=117, AUROC=0.705, AUPRC=0.558, BA=0.632, ECE=0.235.
- primary_recist / melanoma_core_high_evidence / EcoNiche-Opt-AlignedPanelCalibrated: n=117, AUROC=0.655, AUPRC=0.488, BA=0.638, ECE=0.149.
- primary_recist / melanoma_core_high_evidence / ResponseModulesEqualWeight: n=117, AUROC=0.615, AUPRC=0.461, BA=0.627, ECE=0.151.
- primary_recist / melanoma_core_high_evidence / EcoNiche-Opt-NoResistanceModules: n=117, AUROC=0.612, AUPRC=0.458, BA=0.602, ECE=0.133.
- primary_recist / melanoma_core_high_evidence / IFNG_ModuleOnly: n=117, AUROC=0.591, AUPRC=0.439, BA=0.606, ECE=0.151.
- primary_recist / melanoma_core_high_evidence / EcoNiche-Opt-UnsignedStateDirection: n=117, AUROC=0.563, AUPRC=0.425, BA=0.611, ECE=0.166.
- primary_recist / melanoma_core_high_evidence / EcoNiche-Opt-NoResponseModules: n=117, AUROC=0.554, AUPRC=0.481, BA=0.565, ECE=0.149.
- primary_recist / melanoma_recist_supported_primary / EcoNiche-Opt-AlignedPanelNoCalibration: n=131, AUROC=0.685, AUPRC=0.580, BA=0.613, ECE=0.239.
- primary_recist / melanoma_recist_supported_primary / EcoNiche-Opt-AlignedPanelCalibrated: n=131, AUROC=0.630, AUPRC=0.517, BA=0.625, ECE=0.121.
- primary_recist / melanoma_recist_supported_primary / ResponseModulesEqualWeight: n=131, AUROC=0.605, AUPRC=0.487, BA=0.603, ECE=0.114.
- primary_recist / melanoma_recist_supported_primary / EcoNiche-Opt-NoResistanceModules: n=131, AUROC=0.599, AUPRC=0.481, BA=0.597, ECE=0.132.
- primary_recist / melanoma_recist_supported_primary / IFNG_ModuleOnly: n=131, AUROC=0.562, AUPRC=0.465, BA=0.591, ECE=0.112.
- primary_recist / melanoma_recist_supported_primary / EcoNiche-Opt-UnsignedStateDirection: n=131, AUROC=0.557, AUPRC=0.450, BA=0.594, ECE=0.107.
- primary_recist / melanoma_recist_supported_primary / EcoNiche-Opt-NoResponseModules: n=131, AUROC=0.539, AUPRC=0.400, BA=0.584, ECE=0.122.

## Component Evidence Boundary

- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-NoResponseModules: delta AUROC=0.151, delta ECE=0.086, 95% CI [0.009, 0.296], FDR q=0.028 (FDR_supported_component_gain).
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-UnsignedStateDirection: delta AUROC=0.141, delta ECE=0.069, 95% CI [0.069, 0.211], FDR q=0.000 (FDR_supported_component_gain).
- primary_recist / melanoma_core_high_evidence vs IFNG_ModuleOnly: delta AUROC=0.114, delta ECE=0.084, 95% CI [0.059, 0.175], FDR q=0.000 (FDR_supported_component_gain).
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-NoResistanceModules: delta AUROC=0.093, delta ECE=0.102, 95% CI [0.041, 0.146], FDR q=0.000 (FDR_supported_component_gain).
- primary_recist / melanoma_core_high_evidence vs ResponseModulesEqualWeight: delta AUROC=0.090, delta ECE=0.084, 95% CI [0.037, 0.143], FDR q=0.000 (FDR_supported_component_gain).
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-AlignedPanelCalibrated: delta AUROC=0.050, delta ECE=0.086, 95% CI [0.020, 0.080], FDR q=0.002 (calibration_improves_ECE_with_discrimination_tradeoff).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-NoResponseModules: delta AUROC=0.146, delta ECE=0.117, 95% CI [0.004, 0.278], FDR q=0.042 (FDR_supported_component_gain).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-UnsignedStateDirection: delta AUROC=0.127, delta ECE=0.131, 95% CI [0.062, 0.191], FDR q=0.000 (FDR_supported_component_gain).
- primary_recist / melanoma_recist_supported_primary vs IFNG_ModuleOnly: delta AUROC=0.122, delta ECE=0.126, 95% CI [0.067, 0.172], FDR q=0.000 (FDR_supported_component_gain).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-NoResistanceModules: delta AUROC=0.086, delta ECE=0.106, 95% CI [0.043, 0.127], FDR q=0.000 (FDR_supported_component_gain).
- primary_recist / melanoma_recist_supported_primary vs ResponseModulesEqualWeight: delta AUROC=0.080, delta ECE=0.125, 95% CI [0.035, 0.124], FDR q=0.002 (FDR_supported_component_gain).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-AlignedPanelCalibrated: delta AUROC=0.054, delta ECE=0.118, 95% CI [0.029, 0.084], FDR q=0.000 (calibration_improves_ECE_with_discrimination_tradeoff).

## Claim Rule

Use performance-gain language only for variants with positive delta AUROC and FDR support. Components without performance support can still be described as biological representation, interpretation, calibration, or robustness components if their corresponding metric supports that narrower claim.
