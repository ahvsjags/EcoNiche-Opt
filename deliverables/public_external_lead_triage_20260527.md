# Public External Melanoma ICB Lead Triage

This file records public, aggregate, and controlled-access leads for the melanoma anti-PD-1 tumor-tissue validation target. It is a no-fabrication control: a row is not validation evidence until expression, response labels, sample timing, and independence are verified by a registered pipeline.

## Status Counts

- aggregate_duplicate_screen: 1
- controlled_access_required: 3
- eligible_processed: 3
- eligible_processed_duplicate_sensitive: 1
- eligible_processed_evidence_sensitive: 1
- eligible_processed_small: 2
- ineligible_for_melanoma_primary: 1
- low_n_array_public_processed: 1
- panel_transfer_metadata_pending: 1
- panel_transfer_processed: 1
- panel_transfer_public: 1

## Leads

- **LIU_MGSP_PHS000452**: Liu/MGSP pretreatment melanoma anti-PD-1 tumor RNA-seq; status=eligible_processed_duplicate_sensitive; processed=processed; suitability=usable_with_predeclared_train_external_boundary; next=Use only one Liu-derived source per locked external analysis; keep duplicate-source boundaries explicit.
- **GIDE_PRJEB23709**: Gide melanoma anti-PD-1 and anti-PD-1 plus anti-CTLA-4 pretreatment RNA-seq; status=eligible_processed; processed=processed; suitability=usable_with_predeclared_train_external_boundary; next=Keep monotherapy as a primary melanoma component or freeze before any external use; do not mix duplicate iAtlas data as independent evidence.
- **RIAZ_GSE91061**: Riaz nivolumab metastatic melanoma RNA-seq; status=eligible_processed; processed=processed; suitability=usable_with_predeclared_train_external_boundary; next=Maintain pretreatment-only filtering and patient-level fold boundaries.
- **HUGO_GSE78220**: Hugo pretreatment melanoma anti-PD-1 RNA-seq; status=eligible_processed; processed=processed; suitability=usable_with_predeclared_train_external_boundary; next=Keep as part of high-evidence primary melanoma LODO or freeze before external scoring.
- **MGH_GSE115821**: MGH melanoma anti-PD-1 pre/on-treatment RNA-seq; status=eligible_processed_small; processed=processed; suitability=usable_with_predeclared_train_external_boundary; next=Use in endpoint/timing sensitivity analyses with baseline-only filtering.
- **MGH_GSE168204**: New MGH melanoma anti-PD-1 pre/on-treatment RNA-seq; status=eligible_processed_small; processed=processed; suitability=usable_with_predeclared_train_external_boundary; next=Use as a small melanoma sensitivity layer, not as the only top-tier external validation.
- **GSE145996**: Melanoma tumor RNA-seq cohort with limited response-label evidence; status=eligible_processed_evidence_sensitive; processed=processed; suitability=usable_with_predeclared_train_external_boundary; next=Keep as strict stress-test support; harden sample-level response and baseline evidence before headline use.
- **LEE_RIZOS_EGAS00001001552**: Lee/Rizos melanoma immune checkpoint blockade RNA-seq; status=controlled_access_required; processed=not_applicable; suitability=potentially_high_value_after_access; next=Request EGA access and run registered preprocessing; do not substitute or fabricate expression matrices.
- **ABRIL_RODRIGUEZ_PHS001919**: Abril-Rodriguez melanoma anti-PD-1 RNA-seq immune-exclusion cohort; status=controlled_access_required; processed=not_applicable; suitability=potentially_high_value_after_access; next=Request dbGaP access, verify pretreatment tumor timing and RECIST response variables, then run the registered locked external scorer.
- **MGH_HACOHEN_PHS002683**: Combined tumor and immune signals in melanoma checkpoint blockade; status=controlled_access_required; processed=not_applicable; suitability=potentially_high_value_after_access; next=Request dbGaP access, isolate pretreatment anti-PD-1 or PD-1-like tumor RNA-seq samples, and keep this cohort locked outside model development.
- **GSE123728**: Neoadjuvant single-dose PD-1 blockade melanoma NanoString panel; status=panel_transfer_metadata_pending; processed=not_processed; suitability=panel_transfer_not_bulk_primary; next=Curate pre-only samples and response source if used; do not count as strict bulk RNA-seq external validation.
- **GSE165745**: Metastatic melanoma prior to anti-PD-1 NanoString Wnt pathway panel; status=panel_transfer_public; processed=processed; suitability=panel_transfer_not_bulk_primary; next=Implement a targeted-panel scorer only if enough locked panel genes overlap; keep separate from strict bulk RNA-seq validation.
- **GSE122220**: Melanoma tumor samples before checkpoint inhibitors; status=low_n_array_public_processed; processed=processed; suitability=low_n_platform_sensitivity_only; next=Run the registered array preprocessing and score only as low-n platform sensitivity, not as strict bulk RNA-seq validation.
- **GSE93157**: Mixed tumor anti-PD-1 NanoString/pan-cancer ICB panel; status=panel_transfer_processed; processed=processed; suitability=panel_transfer_not_bulk_primary; next=Use for panel compatibility and transfer analyses; separate from pure melanoma anti-PD-1 tumor RNA-seq claims.
- **IMVIGOR210**: Urothelial atezolizumab trial transcriptomic cohort; status=ineligible_for_melanoma_primary; processed=not_processed; suitability=not_melanoma_primary_validation; next=Keep for pan-cancer transfer if processed; do not use as melanoma primary validation.
- **ICBATLAS_TIGER_AGGREGATES**: Aggregated immunotherapy transcriptomic resources; status=aggregate_duplicate_screen; processed=not_applicable; suitability=curation_or_duplicate_screen_only; next=Use as a locator and provenance cross-check; trace each cohort back to original accession before claim wording.

## Claim Boundary

NanoString, microarray, pan-cancer, aggregate-resource, and duplicate-source leads can support panel transfer, provenance checks, or sensitivity analyses, but they do not by themselves satisfy the strict independent melanoma bulk RNA-seq external AUROC target.
