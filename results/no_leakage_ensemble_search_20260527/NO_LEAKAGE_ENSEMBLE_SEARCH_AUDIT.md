# No-leakage Ensemble Search Audit

This registered search tests whether a training-only stacked ensemble can improve primary melanoma and strict external evidence. Feature family, regularization, threshold, and calibration are selected only within training folds or discovery cohorts. Strict external labels are used only for final scoring.

## Summary

- primary_recist / melanoma_core_high_evidence: n=117, AUROC=0.513, AUPRC=0.394, BA=0.550, ECE=0.155; selected=all_l1_C0.3,module_l1_C0.03,module_l2_C0.03.
- primary_recist / melanoma_recist_supported_primary: n=131, AUROC=0.593, AUPRC=0.497, BA=0.524, ECE=0.126; selected=module_l1_C0.3,module_l2_C0.03,module_l2_C0.3.
- strict_recist / strict_melanoma_pd1_like_external: n=116, AUROC=0.579, AUPRC=0.609, BA=0.546, ECE=0.044; selected=module_l2_C0.03.

## Comparison To ModulePriorFixed

- primary_recist / melanoma_core_high_evidence: delta AUROC=-0.192, delta AUPRC=-0.165, delta ECE=-0.080, 95% CI [-0.265, -0.117], FDR q=0.000; ensemble_not_supported.
- primary_recist / melanoma_recist_supported_primary: delta AUROC=-0.091, delta AUPRC=-0.084, delta ECE=-0.112, 95% CI [-0.134, -0.050], FDR q=0.000; ensemble_not_supported.

## Claim Boundary

Only promote the ensemble as the primary model if it improves primary LODO and does not degrade strict external validation. Otherwise, retain it as a negative optimization audit.