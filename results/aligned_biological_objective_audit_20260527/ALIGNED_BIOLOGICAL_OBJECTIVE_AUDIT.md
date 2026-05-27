# Aligned Biological Objective Audit

This audit tests whether adding a biological-prior term to training-only panel-weight selection changes primary melanoma and strict external performance. Candidate selection, thresholding, and calibration use training or discovery cohorts only; locked external labels are used only after the selected candidate and threshold policy are fixed.

## Model Summary

- primary_recist / melanoma_core_high_evidence / EcoNiche-Opt-BioObjectivePanelSearch: n=117, AUROC=0.637, AUPRC=0.469, BA=0.632, ECE=0.153; selected=apm_cyto_resistance,bio_prior.
- primary_recist / melanoma_core_high_evidence / EcoNiche-Opt-NoBioObjectivePanelSearch: n=117, AUROC=0.625, AUPRC=0.460, BA=0.625, ECE=0.152; selected=apm_cyto_resistance,drop_stromal_exclusion,response_equal_resistance_negative.
- primary_recist / melanoma_recist_supported_primary / EcoNiche-Opt-BioObjectivePanelSearch: n=131, AUROC=0.607, AUPRC=0.499, BA=0.620, ECE=0.129; selected=apm_cyto_resistance,bio_prior,bio_prior_strong_suppression,drop_stromal_exclusion.
- primary_recist / melanoma_recist_supported_primary / EcoNiche-Opt-NoBioObjectivePanelSearch: n=131, AUROC=0.618, AUPRC=0.506, BA=0.600, ECE=0.136; selected=apm_cyto_resistance,drop_ifn_t_cell_inflamed,drop_stromal_exclusion.
- strict_recist / strict_melanoma_pd1_like_external / EcoNiche-Opt-BioObjectivePanelSearch: n=116, AUROC=0.572, AUPRC=0.600, BA=0.531, ECE=0.102; selected=bio_prior.
- strict_recist / strict_melanoma_pd1_like_external / EcoNiche-Opt-NoBioObjectivePanelSearch: n=116, AUROC=0.579, AUPRC=0.594, BA=0.507, ECE=0.152; selected=drop_ifn_t_cell_inflamed.

## Paired Bio Objective Comparison

- primary_recist / melanoma_core_high_evidence: delta AUROC=0.013, delta BA=0.007, delta ECE=0.002, 95% bootstrap CI [-0.007, 0.032], FDR q=0.417; point_estimate_biological_objective_gain.
- primary_recist / melanoma_recist_supported_primary: delta AUROC=-0.011, delta BA=0.020, delta ECE=-0.007, 95% bootstrap CI [-0.031, 0.008], FDR q=0.417; biological_objective_calibration_or_threshold_tradeoff.
- strict_recist / strict_melanoma_pd1_like_external: delta AUROC=-0.007, delta BA=0.024, delta ECE=-0.050, 95% bootstrap CI [-0.051, 0.034], FDR q=0.734; biological_objective_calibration_or_threshold_tradeoff.

## Claim Boundary

Use biological-objective performance language only for contexts with positive delta AUROC and FDR support. Otherwise, restrict claims to point-estimate, calibration, threshold-operation, or stability tradeoffs reflected in the table.