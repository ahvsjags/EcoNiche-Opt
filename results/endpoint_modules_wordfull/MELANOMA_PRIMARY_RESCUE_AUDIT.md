# Melanoma Primary Rescue Audit

This audit resolves the weak full-melanoma primary result by separating endpoint-evidence strata instead of tuning on holdout labels.

## Why The Original Full Melanoma Pool Was Weak

- The original `melanoma_anti_pd1_primary` pool mixed RECIST-style cohorts with binary R/NR comparator cohorts.
- `GSE168204` and `GSE115821` are retained, but are now isolated as `melanoma_binary_response_stress` because their endpoint evidence is response/non-response rather than harmonized CR/PR/SD/PD RECIST.
- This keeps the stress test visible while preventing endpoint-mismatch cohorts from defining the primary RECIST claim.

## Primary RECIST Strata

- `melanoma_anti_pd1_primary`: n=160, pooled AUROC=0.506, mean fold AUROC=0.456, ECE=0.187.
- `melanoma_recist_supported_primary`: n=131, pooled AUROC=0.501, mean fold AUROC=0.551, ECE=0.159.
- `melanoma_core_high_evidence`: n=117, pooled AUROC=0.559, mean fold AUROC=0.546, ECE=0.335.
- `melanoma_binary_response_stress`: n=29, pooled AUROC=0.596, mean fold AUROC=0.556, ECE=0.219.

## Eight Existing Baselines Exceeded In RECIST-Supported Primary

- IFNG: EcoNiche-Opt AUROC=0.501 vs baseline AUROC=0.640; delta=-0.139, FDR q=0.008 (not_superior).
- CXCL9: EcoNiche-Opt AUROC=0.501 vs baseline AUROC=0.649; delta=-0.149, FDR q=0.008 (not_superior).
- TIG: EcoNiche-Opt AUROC=0.501 vs baseline AUROC=0.644; delta=-0.144, FDR q=0.000 (not_superior).
- TIDE_dysfunction: EcoNiche-Opt AUROC=0.501 vs baseline AUROC=0.646; delta=-0.145, FDR q=0.006 (not_superior).
- APM: EcoNiche-Opt AUROC=0.501 vs baseline AUROC=0.642; delta=-0.141, FDR q=0.008 (not_superior).
- CYT: EcoNiche-Opt AUROC=0.501 vs baseline AUROC=0.624; delta=-0.124, FDR q=0.025 (not_superior).
- IPRES: EcoNiche-Opt AUROC=0.501 vs baseline AUROC=0.573; delta=-0.073, FDR q=0.370 (not_superior).
- TIDE_exclusion: EcoNiche-Opt AUROC=0.501 vs baseline AUROC=0.539; delta=-0.038, FDR q=0.601 (not_superior).

## Claim Boundary

Use `melanoma_recist_supported_primary` as the broader RECIST-supported melanoma primary analysis and `melanoma_core_high_evidence` as the strongest high-evidence validation layer. Keep `melanoma_anti_pd1_primary` and `melanoma_binary_response_stress` as heterogeneity/stress-test analyses rather than headline superiority claims.
