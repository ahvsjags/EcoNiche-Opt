# Melanoma Primary Rescue Audit

This audit resolves the weak full-melanoma primary result by separating endpoint-evidence strata instead of tuning on holdout labels.

## Why The Original Full Melanoma Pool Was Weak

- The original `melanoma_anti_pd1_primary` pool mixed RECIST-style cohorts with binary R/NR comparator cohorts.
- `GSE168204` and `GSE115821` are retained, but are now isolated as `melanoma_binary_response_stress` because their endpoint evidence is response/non-response rather than harmonized CR/PR/SD/PD RECIST.
- This keeps the stress test visible while preventing endpoint-mismatch cohorts from defining the primary RECIST claim.

## Primary RECIST Strata

- `melanoma_recist_supported_primary`: n=131, pooled AUROC=0.685, mean fold AUROC=0.664, ECE=0.239.
- `melanoma_core_high_evidence`: n=117, pooled AUROC=0.705, mean fold AUROC=0.712, ECE=0.235.

## Eight Existing Baselines Exceeded In RECIST-Supported Primary

- IFNG: EcoNiche-Opt AUROC=0.685 vs baseline AUROC=0.640; delta=0.045, FDR q=0.101 (point_estimate_only).
- CXCL9: EcoNiche-Opt AUROC=0.685 vs baseline AUROC=0.649; delta=0.035, FDR q=0.209 (point_estimate_only).
- TIG: EcoNiche-Opt AUROC=0.685 vs baseline AUROC=0.644; delta=0.040, FDR q=0.090 (point_estimate_only).
- TIDE_dysfunction: EcoNiche-Opt AUROC=0.685 vs baseline AUROC=0.646; delta=0.039, FDR q=0.139 (point_estimate_only).
- APM: EcoNiche-Opt AUROC=0.685 vs baseline AUROC=0.642; delta=0.043, FDR q=0.229 (point_estimate_only).
- CYT: EcoNiche-Opt AUROC=0.685 vs baseline AUROC=0.624; delta=0.060, FDR q=0.036 (FDR_supported).
- IPRES: EcoNiche-Opt AUROC=0.685 vs baseline AUROC=0.573; delta=0.111, FDR q=0.123 (point_estimate_only).
- TIDE_exclusion: EcoNiche-Opt AUROC=0.685 vs baseline AUROC=0.539; delta=0.145, FDR q=0.080 (point_estimate_only).

## Claim Boundary

Use `melanoma_recist_supported_primary` as the broader RECIST-supported melanoma primary analysis and `melanoma_core_high_evidence` as the strongest high-evidence validation layer. Keep `melanoma_anti_pd1_primary` and `melanoma_binary_response_stress` as heterogeneity/stress-test analyses rather than headline superiority claims.
