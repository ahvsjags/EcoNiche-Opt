# External Cohort Gap Audit

- Audit table: `results\curation\external_cohort_gap_audit.tsv`
- Goal: prioritize cohorts that can improve melanoma primary and pan-cancer transfer without weakening traceability.

| Accession | Status | Evidence | Remaining gap | Next action |
|---|---|---|---|---|
| GSE145996 | COMPLETED_INCLUDED | strong_sample_level | Small RNA-seq subset only; SD handling must remain endpoint-stratified. | Use in melanoma/endpoint-stratified analysis; cite exact supplement table evidence. |
| PRJEB23709 | COMPLETED_INCLUDED_SPLIT_PRE | strong_expression_plus_clinical_map | TIGER processed expression was used rather than local FASTQ quantification; provenance must be described as TIGER-processed, not author-generated. | Included as `PRJEB23709_PD1_PRE` for melanoma anti-PD1 primary/core analysis and `PRJEB23709_COMBO_PRE` as secondary/confounded combo therapy. |
| phs000452 | COMPLETED_STRESS_TEST_SECONDARY | TIGER_processed_patient_like_subset | Patient-like anti-PD1 rows are large but lower-performing and have inferred baseline/source provenance; adding them to melanoma core lowers pooled AUROC. | Retain as `PHS000452_LIU_LIKE_PRE` secondary/stress-test; do not use for primary melanoma headline until source/timepoint evidence is hardened. |
| IMvigor210 | EXPRESSION_AND_CLINICAL_PENDING_ENVIRONMENT_BLOCKED | known_public_package_not_local | Need expression matrix plus response/therapy/sample annotations from Bioconductor/ExperimentHub or a verified mirror. | Install/enable Rscript or fetch package RData/RDS assets via a direct verified source; only then add as urothelial anti-PDL1 validation. |
| GSE123728 | RAW_AND_RESPONSE_PENDING | not_local | Need verify whether tumor bulk expression and binary ICB response exist at sample level. | Download GEO supplementary files and sample metadata; include only if expression columns can be matched to response evidence. |
| GSE165745 | RAW_AND_RESPONSE_PENDING | not_local | Need verify expression modality, baseline status, therapy, and binary response. | Download GEO supplementary files and clinical/sample annotations; reject if only single-cell/mechanism without patient response. |
| GSE122220 | RAW_AND_RESPONSE_PENDING | not_local | Need verify sample-level expression and response labels. | Download GEO supplementary files and audit sample-to-patient/response mapping before use. |
| GSE93157 | INCLUDED_ORDER_MATCH_NEEDS_EXTRA_SUPPLEMENT_CONFIRMATION | usable_but_order_match | A supplementary sample map would further harden the expression-column to metadata match. | Search article supplement for a NanoString sample sheet; retain current cohort but flag order-match evidence in Methods. |
| GSE165252 | INCLUDED_SECONDARY_CONFOUNDED_ORDER_MATCH | usable_secondary_but_order_match | Treatment combines chemoradiotherapy and atezolizumab; expression columns are GEO-order matched. | Keep out of primary melanoma/pan-cancer headline; use as secondary/confounded transfer validation. |
