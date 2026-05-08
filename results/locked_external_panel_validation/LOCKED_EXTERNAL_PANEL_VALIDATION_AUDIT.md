# Locked External and Clinical-Assay Panel Validation Audit

This audit locks threshold selection to the discovery melanoma cohorts and evaluates untouched external cohorts and public NanoString panel-transfer cohorts. It is not a prospective wet-lab validation; true prospective validation remains a future clinical study requirement.

- Discovery-only threshold cohorts: GSE91061, GSE78220, PRJEB23709_PD1_PRE
- Locked external/panel cohorts: GSE145996, PHS000452_LIU_LIKE_PRE, PRJEB23709_COMBO_PRE, GSE93157, GSE140901
- Endpoints: primary_recist, strict_recist, clinical_benefit

## Target Model External Metrics

- primary_recist / GSE145996 (locked_external_melanoma_pd1_recist): n=14, AUROC=0.521, balanced_accuracy=0.354, ECE=0.402, threshold=0.369.
- primary_recist / PHS000452_LIU_LIKE_PRE (locked_external_melanoma_pd1_like): n=121, AUROC=0.535, balanced_accuracy=0.534, ECE=0.277, threshold=0.369.
- primary_recist / PRJEB23709_COMBO_PRE (locked_external_melanoma_combination_transfer): n=32, AUROC=0.740, balanced_accuracy=0.652, ECE=0.205, threshold=0.369.
- primary_recist / GSE93157 (nanostring_melanoma_clinical_assay_transfer): n=65, AUROC=0.572, balanced_accuracy=0.500, ECE=0.330, threshold=0.369.
- primary_recist / GSE140901 (nanostring_hcc_clinical_assay_transfer): n=24, AUROC=0.759, balanced_accuracy=0.722, ECE=0.282, threshold=0.369.
- strict_recist / GSE145996 (locked_external_melanoma_pd1_recist): n=13, AUROC=0.575, balanced_accuracy=0.388, ECE=0.370, threshold=0.391.
- strict_recist / PHS000452_LIU_LIKE_PRE (locked_external_melanoma_pd1_like): n=103, AUROC=0.573, balanced_accuracy=0.562, ECE=0.266, threshold=0.391.
- strict_recist / PRJEB23709_COMBO_PRE (locked_external_melanoma_combination_transfer): n=27, AUROC=0.762, balanced_accuracy=0.667, ECE=0.247, threshold=0.391.
- strict_recist / GSE93157 (nanostring_melanoma_clinical_assay_transfer): n=49, AUROC=0.586, balanced_accuracy=0.524, ECE=0.304, threshold=0.391.
- strict_recist / GSE140901 (nanostring_hcc_clinical_assay_transfer): n=14, AUROC=0.833, balanced_accuracy=0.750, ECE=0.149, threshold=0.391.
- clinical_benefit / GSE145996 (locked_external_melanoma_pd1_recist): n=14, AUROC=0.622, balanced_accuracy=0.522, ECE=0.356, threshold=0.424.
- clinical_benefit / PHS000452_LIU_LIKE_PRE (locked_external_melanoma_pd1_like): n=121, AUROC=0.583, balanced_accuracy=0.573, ECE=0.266, threshold=0.424.
- clinical_benefit / PRJEB23709_COMBO_PRE (locked_external_melanoma_combination_transfer): n=32, AUROC=0.724, balanced_accuracy=0.622, ECE=0.308, threshold=0.424.
- clinical_benefit / GSE93157 (nanostring_melanoma_clinical_assay_transfer): n=65, AUROC=0.577, balanced_accuracy=0.530, ECE=0.249, threshold=0.424.
- clinical_benefit / GSE140901 (nanostring_hcc_clinical_assay_transfer): n=24, AUROC=0.711, balanced_accuracy=0.594, ECE=0.227, threshold=0.424.

## Baseline Comparison Boundary

- FDR-supported per-cohort comparisons: 9
- Positive point-estimate comparisons without FDR support: 58
- Supported: clinical_benefit / PRJEB23709_COMBO_PRE vs TIDE_exclusion: delta AUROC=0.404, FDR q=0.011.
- Supported: primary_recist / PRJEB23709_COMBO_PRE vs IPRES: delta AUROC=0.377, FDR q=0.000.
- Supported: primary_recist / PRJEB23709_COMBO_PRE vs TIDE_exclusion: delta AUROC=0.420, FDR q=0.000.
- Supported: primary_recist / GSE140901 vs IPRES: delta AUROC=0.491, FDR q=0.003.
- Supported: primary_recist / GSE140901 vs TIDE_exclusion: delta AUROC=0.583, FDR q=0.000.
- Supported: strict_recist / PRJEB23709_COMBO_PRE vs IPRES: delta AUROC=0.365, FDR q=0.027.
- Supported: strict_recist / PRJEB23709_COMBO_PRE vs TIDE_exclusion: delta AUROC=0.460, FDR q=0.005.
- Supported: strict_recist / GSE140901 vs IPRES: delta AUROC=0.542, FDR q=0.021.
- Supported: strict_recist / GSE140901 vs TIDE_exclusion: delta AUROC=0.562, FDR q=0.011.

## External Signature-Family Omnibus

- clinical_benefit / all_locked_external_and_panel: target AUROC=0.604, mean signature AUROC=0.566, best signature AUROC=0.628, delta vs family mean=0.038, two-sided FDR q=0.031 (family_two_sided_FDR_supported).
- clinical_benefit / melanoma_external_and_panel: target AUROC=0.594, mean signature AUROC=0.559, best signature AUROC=0.616, delta vs family mean=0.035, two-sided FDR q=0.058 (family_pre_directional_FDR_supported).
- clinical_benefit / strict_pd1_like_external: target AUROC=0.586, mean signature AUROC=0.553, best signature AUROC=0.605, delta vs family mean=0.032, two-sided FDR q=0.307 (family_point_estimate_only).
- clinical_benefit / nanostring_panel_transfer: target AUROC=0.606, mean signature AUROC=0.567, best signature AUROC=0.637, delta vs family mean=0.039, two-sided FDR q=0.172 (family_point_estimate_only).
- primary_recist / all_locked_external_and_panel: target AUROC=0.577, mean signature AUROC=0.549, best signature AUROC=0.610, delta vs family mean=0.028, two-sided FDR q=0.071 (family_pre_directional_FDR_supported).
- primary_recist / melanoma_external_and_panel: target AUROC=0.567, mean signature AUROC=0.543, best signature AUROC=0.599, delta vs family mean=0.025, two-sided FDR q=0.125 (family_point_estimate_only).
- primary_recist / strict_pd1_like_external: target AUROC=0.531, mean signature AUROC=0.519, best signature AUROC=0.556, delta vs family mean=0.013, two-sided FDR q=0.547 (family_point_estimate_only).
- primary_recist / nanostring_panel_transfer: target AUROC=0.609, mean signature AUROC=0.557, best signature AUROC=0.650, delta vs family mean=0.052, two-sided FDR q=0.172 (family_point_estimate_only).
- strict_recist / all_locked_external_and_panel: target AUROC=0.610, mean signature AUROC=0.570, best signature AUROC=0.645, delta vs family mean=0.040, two-sided FDR q=0.031 (family_two_sided_FDR_supported).
- strict_recist / melanoma_external_and_panel: target AUROC=0.596, mean signature AUROC=0.560, best signature AUROC=0.629, delta vs family mean=0.036, two-sided FDR q=0.058 (family_pre_directional_FDR_supported).
- strict_recist / strict_pd1_like_external: target AUROC=0.571, mean signature AUROC=0.542, best signature AUROC=0.590, delta vs family mean=0.030, two-sided FDR q=0.307 (family_point_estimate_only).
- strict_recist / nanostring_panel_transfer: target AUROC=0.633, mean signature AUROC=0.580, best signature AUROC=0.687, delta vs family mean=0.054, two-sided FDR q=0.172 (family_point_estimate_only).

## NanoString Panel Gene Coverage

- GSE140901: mean module gene coverage=0.746.
- GSE93157: mean module gene coverage=0.746.

## Claim Boundary

Allowed: 'locked external/panel-transfer validation was run without threshold leakage' and cohort-specific performance statements from the tables. Avoid: 'prospective clinical validation completed' unless a new prospective assay cohort is generated outside this retrospective public-data pipeline.
