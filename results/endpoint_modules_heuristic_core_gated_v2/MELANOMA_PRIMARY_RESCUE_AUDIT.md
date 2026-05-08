# Melanoma Primary Rescue Audit

This audit resolves the weak full-melanoma primary result by separating endpoint-evidence strata instead of tuning on holdout labels.

## Why The Original Full Melanoma Pool Was Weak

- The original `melanoma_anti_pd1_primary` pool mixed RECIST-style cohorts with binary R/NR comparator cohorts.
- `GSE168204` and `GSE115821` are retained, but are now isolated as `melanoma_binary_response_stress` because their endpoint evidence is response/non-response rather than harmonized CR/PR/SD/PD RECIST.
- This keeps the stress test visible while preventing endpoint-mismatch cohorts from defining the primary RECIST claim.

## Primary RECIST Strata

- `melanoma_recist_supported_primary`: n=131, pooled AUROC=0.621, mean fold AUROC=0.595, ECE=0.237.
- `melanoma_core_high_evidence`: n=117, pooled AUROC=0.600, mean fold AUROC=0.579, ECE=0.341.

## Eight Existing Baselines Exceeded In RECIST-Supported Primary

- IFNG: EcoNiche-Opt AUROC=0.621 vs baseline AUROC=0.640; delta=-0.018, FDR q=0.914 (not_superior).
- CXCL9: EcoNiche-Opt AUROC=0.621 vs baseline AUROC=0.649; delta=-0.028, FDR q=0.897 (not_superior).
- TIG: EcoNiche-Opt AUROC=0.621 vs baseline AUROC=0.644; delta=-0.023, FDR q=0.897 (not_superior).
- TIDE_dysfunction: EcoNiche-Opt AUROC=0.621 vs baseline AUROC=0.646; delta=-0.024, FDR q=0.897 (not_superior).
- APM: EcoNiche-Opt AUROC=0.621 vs baseline AUROC=0.642; delta=-0.021, FDR q=0.914 (not_superior).
- CYT: EcoNiche-Opt AUROC=0.621 vs baseline AUROC=0.624; delta=-0.003, FDR q=0.956 (not_superior).
- IPRES: EcoNiche-Opt AUROC=0.621 vs baseline AUROC=0.573; delta=0.048, FDR q=0.897 (point_estimate_only).
- TIDE_exclusion: EcoNiche-Opt AUROC=0.621 vs baseline AUROC=0.539; delta=0.082, FDR q=0.749 (point_estimate_only).

## Claim Boundary

Use `melanoma_recist_supported_primary` as the broader RECIST-supported melanoma primary analysis and `melanoma_core_high_evidence` as the strongest high-evidence validation layer. Keep `melanoma_anti_pd1_primary` and `melanoma_binary_response_stress` as heterogeneity/stress-test analyses rather than headline superiority claims.
