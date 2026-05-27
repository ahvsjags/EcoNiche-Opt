# Prospective and Clinical-Assay Validation Package

This folder freezes the EcoNiche-Opt panel and analysis plan for independent qPCR/NanoString or RNA-seq validation in pretreatment melanoma tumor tissue. It complements, but does not replace, the retrospective locked external and NanoString transfer analyses in `results/locked_external_panel_validation_calibrated_20260519`.

Generated artifacts:

- `locked_scoring_spec.json`: frozen score formula, endpoint thresholds, and claim boundary.
- `locked_scoring_spec.sha256`: file hash to confirm the scoring spec has not changed.
- `locked_panel_genes.tsv`: 62-unique-gene qPCR/NanoString panel with gene-module rows, directions, and module weights.
- `assay_sample_manifest_template.tsv`: assay/sample metadata template.
- `clinical_annotation_template.tsv`: response-label curation template.
- `prospective_validation_protocol.md`: validation protocol.
- `statistical_analysis_plan.md`: predeclared analysis plan.
- `clinical_partner_intake_checklist.md`: what to request from a hospital or clinical collaborator.
- `tumor_tissue_assay_sop.md`: tumor-tissue processing and locked scoring SOP.
- `validation_readiness_checklist.tsv`: pre-scoring readiness gate.

The primary clinical scenario is pretreatment melanoma tumor tissue before anti-PD-1/anti-PD-1-based therapy. Blood, plasma, serum, PBMC, or on-treatment samples require a separately trained and validated model before clinical claims can be made.
