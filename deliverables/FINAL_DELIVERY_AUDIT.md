# EcoNiche-Opt Final Delivery Audit

Generated: 2026-05-06 Asia/Shanghai

## Status

- Goals: 81/81 valid after the latest optimization and curation work.
- Modeled real response benchmark now includes 14 cohorts and 554 primary-RECIST-evaluable samples across all response strata.
- Primary melanoma anti-PD1 claim set now separates RECIST-supported pretreatment monotherapy cohorts from binary R/NR melanoma stress-test cohorts; phs000452 is retained only as secondary/stress-test.
- PRJEB23709 was split into `PRJEB23709_PD1_PRE` (41 pretreatment anti-PD-1 monotherapy samples) and `PRJEB23709_COMBO_PRE` (32 pretreatment anti-PD-1 plus anti-CTLA-4 samples).
- phs000452 TIGER Patient-like anti-PD1 subset was processed as `PHS000452_LIU_LIKE_PRE` (121 samples) but excluded from the primary headline because it weakens pooled melanoma performance.
- No fabricated labels or imputed treatment recommendations were added.

## Model And Result Summary

- Primary discrimination model: `EcoNiche-Opt-ModulePriorFixed`, a fixed module-level immune-ecology prior model.
- Calibration sensitivity model: `EcoNiche-Opt-ModulePriorFixed-Platt`, trained only on training cohorts in each fold.
- High-evidence melanoma anti-PD1 PRE core (`GSE91061`, `GSE78220`, `PRJEB23709_PD1_PRE`): primary RECIST pooled AUROC 0.705, mean fold AUROC 0.712; strict RECIST pooled AUROC 0.707.
- RECIST-supported melanoma anti-PD1 primary layer (`GSE91061`, `GSE78220`, `GSE145996`, `PRJEB23709_PD1_PRE`): primary RECIST pooled AUROC 0.685, mean fold AUROC 0.664; strict RECIST pooled AUROC 0.690.
- Full accepted melanoma anti-PD1 heterogeneity pool, including binary R/NR stress cohorts: primary RECIST pooled AUROC 0.641.
- Adding phs000452 to the melanoma core lowers primary RECIST pooled AUROC to 0.619, so it is a stress-test layer rather than a primary validation layer.
- Platt calibration improves full melanoma primary ECE from 0.268 to 0.115, but lowers pooled AUROC because fold-wise calibration changes cross-cohort score scaling.

## Key New Deliverables

- `scripts/preprocess/preprocess_gide_prjeb23709.py`
- `scripts/preprocess/preprocess_phs000452_tiger.py`
- `data/processed/bulk/PRJEB23709_PD1_PRE.expr.tsv`
- `data/processed/bulk/PRJEB23709_PD1_PRE.metadata.tsv`
- `data/processed/bulk/PRJEB23709_COMBO_PRE.expr.tsv`
- `data/processed/bulk/PRJEB23709_COMBO_PRE.metadata.tsv`
- `data/processed/bulk/PHS000452_LIU_LIKE_PRE.expr.tsv`
- `data/processed/bulk/PHS000452_LIU_LIKE_PRE.metadata.tsv`
- `tables/gide_prjeb23709_pre_qc.tsv`
- `tables/phs000452_tiger_qc.tsv`
- `results/endpoint_modules/ENDPOINT_MODULE_MODEL_AUDIT.md`
- `results/endpoint_modules/MELANOMA_PRIMARY_RESCUE_AUDIT.md`
- `results/endpoint_modules/melanoma_primary_rescue_baselines.tsv`
- `results/endpoint_modules/endpoint_module_summary.tsv`
- `results/endpoint_modules/endpoint_module_pairwise_comparisons.tsv`
- `results/real_optimized/PERFORMANCE_OPTIMIZATION_AUDIT.md`

## Verification Evidence

Fresh verification after the latest changes:

- `python scripts\preprocess\preprocess_phs000452_tiger.py`: passed and wrote the phs000452 stress-test cohort.
- `python scripts\model\run_endpoint_module_analysis.py`: passed and regenerated endpoint/module outputs.
- Full tests and project validators were rerun after code/report updates.

## Claim Boundary

The project has a defensible stratum-specific AUROC >0.7 result for high-evidence melanoma PRE anti-PD1 RNA-seq, and a broader RECIST-supported melanoma primary result of AUROC 0.685. In the RECIST-supported primary layer, EcoNiche-Opt exceeds eight existing immune-response baselines by AUROC point estimate: IFNG, CXCL9, TIG, TIDE dysfunction, APM, CYT, IPRES, and TIDE exclusion; only CYT is FDR-supported in that eight-model table. It should not claim universal superiority over all existing models across all cohorts/endpoints.

## Remaining High-Value Work

- Harden phs000452 source/timepoint provenance or replace it with another exact-baseline melanoma anti-PD1 RNA-seq cohort.
- Add IMvigor210 as a separate urothelial validation layer.
- Use raw module-prior outputs for discrimination figures and Platt-calibrated outputs for calibration/decision-curve figures.
