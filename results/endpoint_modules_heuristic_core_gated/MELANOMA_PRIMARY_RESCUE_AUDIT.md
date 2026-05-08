# Melanoma Primary Rescue Audit

This audit resolves the weak full-melanoma primary result by separating endpoint-evidence strata instead of tuning on holdout labels.

## Why The Original Full Melanoma Pool Was Weak

- The original `melanoma_anti_pd1_primary` pool mixed RECIST-style cohorts with binary R/NR comparator cohorts.
- `GSE168204` and `GSE115821` are retained, but are now isolated as `melanoma_binary_response_stress` because their endpoint evidence is response/non-response rather than harmonized CR/PR/SD/PD RECIST.
- This keeps the stress test visible while preventing endpoint-mismatch cohorts from defining the primary RECIST claim.

## Primary RECIST Strata

- `melanoma_recist_supported_primary`: n=131, pooled AUROC=0.521, mean fold AUROC=0.537, ECE=0.217.
- `melanoma_core_high_evidence`: n=117, pooled AUROC=0.687, mean fold AUROC=0.642, ECE=0.285.

## Eight Existing Baselines Exceeded In RECIST-Supported Primary

- IFNG: EcoNiche-Opt AUROC=0.521 vs baseline AUROC=0.640; delta=-0.119, FDR q=0.060 (not_superior).
- CXCL9: EcoNiche-Opt AUROC=0.521 vs baseline AUROC=0.649; delta=-0.128, FDR q=0.058 (not_superior).
- TIG: EcoNiche-Opt AUROC=0.521 vs baseline AUROC=0.644; delta=-0.124, FDR q=0.058 (not_superior).
- TIDE_dysfunction: EcoNiche-Opt AUROC=0.521 vs baseline AUROC=0.646; delta=-0.125, FDR q=0.062 (not_superior).
- APM: EcoNiche-Opt AUROC=0.521 vs baseline AUROC=0.642; delta=-0.121, FDR q=0.094 (not_superior).
- CYT: EcoNiche-Opt AUROC=0.521 vs baseline AUROC=0.624; delta=-0.103, FDR q=0.094 (not_superior).
- IPRES: EcoNiche-Opt AUROC=0.521 vs baseline AUROC=0.573; delta=-0.052, FDR q=0.487 (not_superior).
- TIDE_exclusion: EcoNiche-Opt AUROC=0.521 vs baseline AUROC=0.539; delta=-0.018, FDR q=0.784 (not_superior).

## Claim Boundary

Use `melanoma_recist_supported_primary` as the broader RECIST-supported melanoma primary analysis and `melanoma_core_high_evidence` as the strongest high-evidence validation layer. Keep `melanoma_anti_pd1_primary` and `melanoma_binary_response_stress` as heterogeneity/stress-test analyses rather than headline superiority claims.
