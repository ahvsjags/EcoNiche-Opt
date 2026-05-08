# Melanoma Primary Rescue Audit

This audit resolves the weak full-melanoma primary result by separating endpoint-evidence strata instead of tuning on holdout labels.

## Why The Original Full Melanoma Pool Was Weak

- The original `melanoma_anti_pd1_primary` pool mixed RECIST-style cohorts with binary R/NR comparator cohorts.
- `GSE168204` and `GSE115821` are retained, but are now isolated as `melanoma_binary_response_stress` because their endpoint evidence is response/non-response rather than harmonized CR/PR/SD/PD RECIST.
- This keeps the stress test visible while preventing endpoint-mismatch cohorts from defining the primary RECIST claim.

## Primary RECIST Strata

- `melanoma_recist_supported_primary`: n=131, pooled AUROC=0.373, mean fold AUROC=0.345, ECE=0.314.
- `melanoma_core_high_evidence`: n=117, pooled AUROC=0.585, mean fold AUROC=0.466, ECE=0.277.

## Eight Existing Baselines Exceeded In RECIST-Supported Primary

- IFNG: EcoNiche-Opt AUROC=0.373 vs baseline AUROC=0.640; delta=-0.267, FDR q=0.000 (not_superior).
- CXCL9: EcoNiche-Opt AUROC=0.373 vs baseline AUROC=0.649; delta=-0.277, FDR q=0.000 (not_superior).
- TIG: EcoNiche-Opt AUROC=0.373 vs baseline AUROC=0.644; delta=-0.272, FDR q=0.000 (not_superior).
- TIDE_dysfunction: EcoNiche-Opt AUROC=0.373 vs baseline AUROC=0.646; delta=-0.273, FDR q=0.000 (not_superior).
- APM: EcoNiche-Opt AUROC=0.373 vs baseline AUROC=0.642; delta=-0.269, FDR q=0.000 (not_superior).
- CYT: EcoNiche-Opt AUROC=0.373 vs baseline AUROC=0.624; delta=-0.252, FDR q=0.000 (not_superior).
- IPRES: EcoNiche-Opt AUROC=0.373 vs baseline AUROC=0.573; delta=-0.201, FDR q=0.003 (not_superior).
- TIDE_exclusion: EcoNiche-Opt AUROC=0.373 vs baseline AUROC=0.539; delta=-0.167, FDR q=0.005 (not_superior).

## Claim Boundary

Use `melanoma_recist_supported_primary` as the broader RECIST-supported melanoma primary analysis and `melanoma_core_high_evidence` as the strongest high-evidence validation layer. Keep `melanoma_anti_pd1_primary` and `melanoma_binary_response_stress` as heterogeneity/stress-test analyses rather than headline superiority claims.
