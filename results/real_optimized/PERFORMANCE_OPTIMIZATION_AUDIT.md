# EcoNiche-Opt Performance Optimization Audit

## Current Locked Models

- `EcoNiche-Opt-ModulePriorFixed`: fixed module-level immune-ecology prior over IFN/T-cell inflamed, cytotoxic CD8, exhaustion/checkpoint, antigen presentation, myeloid suppression, stromal exclusion, and TRM/TLS modules. This remains the primary discrimination/claim-gate model.
- `EcoNiche-Opt-ModulePriorFixed-Platt`: fold-wise Platt calibrated version trained only on training cohorts. This is a calibration/clinical-utility sensitivity model, not the AUROC headline model.
- `EcoNiche-Opt-AdaptiveConsensus`: probability consensus over the immune composite, module prior, IFNG, TIG, TIDE dysfunction, and CXCL9. This remains a calibrated sensitivity model and is competitive in pan-cancer pooled summaries.
- Validation: leave-one-dataset-out (LODO), endpoint sensitivity, and explicit stratum sensitivity. Candidate selection, thresholds, and calibration do not use holdout labels.

## Data Curation Update

- Curated modeled real response benchmark after the latest stress test: 14 cohorts and 554 primary-RECIST-evaluable samples across all response strata.
- Primary melanoma anti-PD1 set remains locked to high-confidence pretreatment monotherapy cohorts and excludes ambiguous/stress-test cohorts.
- Newly completed previously: Gide/PRJEB23709 TIGER processed expression plus local evidence-checked response map.
- Newly tested in this round: TIGER `Melanoma-phs000452` Patient-like anti-PD1 subset (`PHS000452_LIU_LIKE_PRE`, 121 samples). It is retained as secondary/stress-test because it lowers pooled melanoma AUROC, suggesting source/processing/label heterogeneity.
- PRJEB23709 modeled PRE samples: 41 anti-PD-1 monotherapy samples and 32 anti-PD-1 plus anti-CTLA-4 samples.

## Primary RECIST-Style Performance

- High-evidence melanoma anti-PD1 PRE core (`GSE91061`, `GSE78220`, `PRJEB23709_PD1_PRE`): `EcoNiche-Opt-ModulePriorFixed` pooled AUROC 0.705, mean fold AUROC 0.712, pooled AUPRC 0.558.
- Full melanoma anti-PD1 primary stratum including small/heterogeneous but accepted primary cohorts: pooled AUROC 0.641, mean fold AUROC 0.593.
- phs000452 stress-test layer: adding `PHS000452_LIU_LIKE_PRE` to the melanoma core lowers pooled AUROC to 0.619, so it should not be used to inflate the primary claim.
- Pan-cancer all response cohorts, including secondary/stress-test cohorts: `EcoNiche-Opt-AdaptiveConsensus` pooled AUROC 0.630; `EcoNiche-Opt-ModulePriorFixed` pooled AUROC 0.629.
- Pan-cancer without secondary/stress-test cohorts: `EcoNiche-Opt-ModulePriorFixed` pooled AUROC 0.642.

## Calibration

- Full melanoma primary ECE improves from 0.268 (`ModulePriorFixed`) to 0.115 (`ModulePriorFixed-Platt`), while pooled AUROC drops from 0.641 to 0.588 because fold-wise calibration changes cross-cohort score scaling.
- High-evidence melanoma primary ECE improves from 0.235 to 0.149 with Platt calibration, while pooled AUROC drops from 0.705 to 0.655; mean fold AUROC is preserved at 0.712.
- Strict RECIST full melanoma ECE improves from 0.259 to 0.071 with Platt calibration.
- The manuscript should present raw module prior for discrimination and Platt-calibrated module prior for probability calibration/decision-curve sensitivity.

## Strong Signature Claim Gate

- In high-evidence melanoma primary RECIST, `ModulePriorFixed` exceeds IFNG, CXCL9, TIG, APM, TIDE dysfunction, CYT, PDCD1LG2, IPRES, and TIDE exclusion by point estimate.
- The largest primary RECIST deltas are versus TIDE exclusion (+0.161), IPRES (+0.140), PDCD1LG2 (+0.116), CYT (+0.076), and the previous immune composite (+0.075). Deltas versus IFNG/CXCL9/TIG/APM remain positive but not FDR-significant.
- In strict RECIST high-evidence melanoma, the model reaches AUROC 0.707 and remains above IFNG/CXCL9/TIG/TIDE dysfunction/APM by point estimate; several confidence intervals are positive, but FDR values remain borderline rather than definitive.
- Therefore the defensible headline is stratum-specific and evidence-bounded: EcoNiche-Opt reaches about 0.70 AUROC in high-evidence melanoma PRE anti-PD1 validation and improves point estimates over strong immune signatures, but it is not yet proven superior to every comparator across all cohorts/endpoints.

## Remaining Boundaries

- Full melanoma primary remains moderate because GSE168204, GSE145996, and GSE115821 are small/heterogeneous and lower-performing.
- phs000452 is useful as a large stress test but currently weakens the melanoma pooled claim; it needs deeper source/timepoint hardening before being upgraded to primary.
- Calibration improves substantially with Platt, but the best-calibrated model is not the best pooled-AUROC model.
- Clinical-benefit endpoint performance is weaker than strict/primary RECIST, so endpoint claims must stay separated.

## Files

- Endpoint/module outputs: `results/endpoint_modules/`.
- PRJEB23709 preprocessor: `scripts/preprocess/preprocess_gide_prjeb23709.py`.
- phs000452 stress-test preprocessor: `scripts/preprocess/preprocess_phs000452_tiger.py`.
- PRJEB23709 QC report: `tables/gide_prjeb23709_pre_qc.tsv`.
- phs000452 QC report: `tables/phs000452_tiger_qc.tsv`.
- PRJEB23709 processed data: `data/processed/bulk/PRJEB23709_PD1_PRE.*` and `data/processed/bulk/PRJEB23709_COMBO_PRE.*`.
- phs000452 processed stress-test data: `data/processed/bulk/PHS000452_LIU_LIKE_PRE.*`.

## Next Work

- Harden phs000452 source/timepoint provenance or replace it with another exact-baseline melanoma anti-PD1 RNA-seq cohort.
- Add IMvigor210 as a separate urothelial anti-PDL1 validation layer, not as evidence for melanoma claims.
- Use calibrated decision curves only with the Platt model and discrimination figures with the raw module-prior model.
