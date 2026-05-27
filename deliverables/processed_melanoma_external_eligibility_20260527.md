# Processed Melanoma External Eligibility Audit

This audit checks all processed bulk metadata files for possible strict melanoma pretreatment anti-PD-1 tumor-tissue external-validation use.

## Status Counts

- low_n_array_platform_sensitivity: 1
- needs_manual_source_hardening: 2
- not_melanoma_primary: 5
- not_usable_empty_metadata: 3
- panel_transfer_not_bulk_strict: 1
- primary_discovery_or_lodo: 3
- secondary_small_melanoma_sensitivity: 2
- strict_external_current: 2

## Cohort Decisions

- `GSE122220`: low_n_array_platform_sensitivity (platform=Illumina_HumanHT12_V4_expression_beadchip; n=10); n=10.
- `GSE165745`: needs_manual_source_hardening (role=candidate_external_requires_metadata_curation; notes=Registry entry required by implementation contract; no benchmark use until patient labels pass QC.); n=24.
- `PRJEB23709_COMBO_PRE`: needs_manual_source_hardening (role=locked_external_if_accessible; notes=Processed expression obtained from TIGER immunotherapy download interface; PRE samples are split into anti-PD1 monothera); n=32.
- `GSE136961`: not_melanoma_primary (registry cancer_type=NSCLC); n=21.
- `GSE140901`: not_melanoma_primary (registry cancer_type=hepatocellular_carcinoma); n=24.
- `GSE165252`: not_melanoma_primary (registry cancer_type=esophageal_cancer); n=32.
- `GSE176307`: not_melanoma_primary (registry cancer_type=urothelial_cancer); n=88.
- `GSE67501`: not_melanoma_primary (registry cancer_type=renal_cell_carcinoma); n=11.
- `GSE121810`: not_usable_empty_metadata (metadata has zero rows); n=0.
- `GSE183924`: not_usable_empty_metadata (metadata has zero rows); n=0.
- `GSE244982`: not_usable_empty_metadata (metadata has zero rows); n=0.
- `GSE93157`: panel_transfer_not_bulk_strict (platform=NanoString_nCounter_immune_panel); n=65.
- `GSE78220`: primary_discovery_or_lodo (used in high-evidence primary melanoma discovery/LODO boundary); n=27.
- `GSE91061`: primary_discovery_or_lodo (used in high-evidence primary melanoma discovery/LODO boundary); n=49.
- `PRJEB23709_PD1_PRE`: primary_discovery_or_lodo (used in high-evidence primary melanoma discovery/LODO boundary); n=41.
- `GSE115821`: secondary_small_melanoma_sensitivity (eligible melanoma pretreatment ICB-like cohort, but small and response/timing evidence is less robust than primary discovery cohorts); n=8.
- `GSE168204`: secondary_small_melanoma_sensitivity (eligible melanoma pretreatment ICB-like cohort, but small and response/timing evidence is less robust than primary discovery cohorts); n=21.
- `GSE145996`: strict_external_current (used in current strict melanoma external stress test); n=14.
- `PHS000452_LIU_LIKE_PRE`: strict_external_current (used in current strict melanoma external stress test); n=121.

## Conclusion

No overlooked large public bulk RNA-seq pretreatment melanoma anti-PD-1 external cohort was found among the currently processed metadata files. The remaining path to a top-tier strict external claim is controlled-access acquisition or a newly discovered independent public cohort with hardened sample-level response evidence.
