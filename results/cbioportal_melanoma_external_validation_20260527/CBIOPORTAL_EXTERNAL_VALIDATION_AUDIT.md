# cBioPortal Melanoma External Validation Audit

This audit adds public cBioPortal melanoma ICB expression profiles as a real-data source. Discovery-only cohorts are fixed to GSE91061, GSE78220 and PRJEB23709_PD1_PRE for calibration and thresholding. cBioPortal external labels are used only for evaluation.

## Group Results

- strict_recist / cbio_liu_dfci_only: n=103, AUROC=0.561, AUPRC=0.571, balanced accuracy=0.586, ECE=0.071 (independent_cbioportal_liu_external).
- strict_recist / strict_cbio_liu_plus_gse145996: n=116, AUROC=0.556, AUPRC=0.592, balanced accuracy=0.565, ECE=0.100 (strict_melanoma_pd1_like_external_cbioportal).
- strict_recist / current_tiger_phs_plus_gse145996: n=116, AUROC=0.571, AUPRC=0.600, balanced accuracy=0.542, ECE=0.108 (current_tiger_strict_external_reference).
- strict_recist / cbio_iatlas_liu_duplicate_crosscheck: n=104, AUROC=0.539, AUPRC=0.540, balanced accuracy=0.570, ECE=0.119 (duplicate_source_crosscheck_not_independent).
- strict_recist / cbio_discovery_overlap_crosscheck: n=121, AUROC=0.716, AUPRC=0.722, balanced accuracy=0.629, ECE=0.055 (discovery_overlap_crosscheck_not_external).
- primary_recist / cbio_liu_dfci_only: n=121, AUROC=0.528, AUPRC=0.478, balanced accuracy=0.567, ECE=0.102 (independent_cbioportal_liu_external).
- primary_recist / strict_cbio_liu_plus_gse145996: n=135, AUROC=0.520, AUPRC=0.501, balanced accuracy=0.545, ECE=0.120 (strict_melanoma_pd1_like_external_cbioportal).
- primary_recist / current_tiger_phs_plus_gse145996: n=135, AUROC=0.531, AUPRC=0.517, balanced accuracy=0.514, ECE=0.124 (current_tiger_strict_external_reference).
- primary_recist / cbio_iatlas_liu_duplicate_crosscheck: n=122, AUROC=0.516, AUPRC=0.463, balanced accuracy=0.525, ECE=0.165 (duplicate_source_crosscheck_not_independent).
- primary_recist / cbio_discovery_overlap_crosscheck: n=148, AUROC=0.707, AUPRC=0.635, balanced accuracy=0.633, ECE=0.064 (discovery_overlap_crosscheck_not_external).

## Family Claim Gate

- strict_recist / cbio_liu_dfci_only: target AUROC=0.561, family mean AUROC=0.542, delta=0.020, q=0.525, best signature=APM (family_point_estimate_only).
- strict_recist / strict_cbio_liu_plus_gse145996: target AUROC=0.556, family mean AUROC=0.542, delta=0.015, q=0.600, best signature=TIDE_dysfunction (family_point_estimate_only).
- strict_recist / current_tiger_phs_plus_gse145996: target AUROC=0.571, family mean AUROC=0.542, delta=0.030, q=0.193, best signature=CXCL9 (family_point_estimate_only).
- strict_recist / cbio_iatlas_liu_duplicate_crosscheck: target AUROC=0.539, family mean AUROC=0.533, delta=0.007, q=0.759, best signature=CXCL9 (family_point_estimate_only).
- strict_recist / cbio_discovery_overlap_crosscheck: target AUROC=0.716, family mean AUROC=0.642, delta=0.075, q=0.000, best signature=APM (family_two_sided_FDR_supported).
- primary_recist / cbio_liu_dfci_only: target AUROC=0.528, family mean AUROC=0.516, delta=0.012, q=0.666, best signature=CXCL9 (family_point_estimate_only).
- primary_recist / strict_cbio_liu_plus_gse145996: target AUROC=0.520, family mean AUROC=0.514, delta=0.007, q=0.799, best signature=TIDE_dysfunction (family_point_estimate_only).
- primary_recist / current_tiger_phs_plus_gse145996: target AUROC=0.531, family mean AUROC=0.519, delta=0.014, q=0.507, best signature=CXCL9 (family_point_estimate_only).
- primary_recist / cbio_iatlas_liu_duplicate_crosscheck: target AUROC=0.516, family mean AUROC=0.520, delta=-0.004, q=0.859, best signature=CXCL9 (family_not_superior).
- primary_recist / cbio_discovery_overlap_crosscheck: target AUROC=0.707, family mean AUROC=0.642, delta=0.066, q=0.000, best signature=APM (family_two_sided_FDR_supported).

## Coverage

- CBIO_IATLAS_GIDE_2019_PRE: total available module-gene hits=63.
- CBIO_IATLAS_HUGO_2016_PRE: total available module-gene hits=63.
- CBIO_IATLAS_LIU_2019_PRE: total available module-gene hits=63.
- CBIO_IATLAS_RIAZ_2017_PRE: total available module-gene hits=63.
- CBIO_LIU_DFCI_2019_PRE: total available module-gene hits=61.
- GSE115821: total available module-gene hits=63.
- GSE122220: total available module-gene hits=54.
- GSE136961: total available module-gene hits=47.
- GSE140901: total available module-gene hits=47.
- GSE145996: total available module-gene hits=56.
- GSE165252: total available module-gene hits=62.
- GSE165745: total available module-gene hits=8.
- GSE168204: total available module-gene hits=63.
- GSE176307: total available module-gene hits=62.
- GSE67501: total available module-gene hits=63.
- GSE78220: total available module-gene hits=63.
- GSE91061: total available module-gene hits=63.
- GSE93157: total available module-gene hits=47.
- PHS000452_LIU_LIKE_PRE: total available module-gene hits=63.
- PRJEB23709_COMBO_PRE: total available module-gene hits=63.
- PRJEB23709_PD1_PRE: total available module-gene hits=63.
- demo_cohort_1: total available module-gene hits=14.
- demo_cohort_2: total available module-gene hits=14.
- demo_cohort_3: total available module-gene hits=14.
- demo_cohort_4: total available module-gene hits=14.