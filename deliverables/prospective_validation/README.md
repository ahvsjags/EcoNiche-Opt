# Prospective and Clinical-Assay Validation Package

This folder freezes the EcoNiche-Opt panel and analysis plan for future independent qPCR/NanoString validation. It complements, but does not replace, the retrospective locked external and NanoString transfer analyses in `results/locked_external_panel_validation`.

Generated artifacts:

- `locked_scoring_spec.json`: frozen score formula, endpoint thresholds, and claim boundary.
- `locked_panel_genes.tsv`: qPCR/NanoString panel gene list with module weights.
- `assay_sample_manifest_template.tsv`: assay/sample metadata template.
- `clinical_annotation_template.tsv`: response-label curation template.
- `prospective_validation_protocol.md`: validation protocol.
- `statistical_analysis_plan.md`: predeclared analysis plan.

No prospective performance claim is allowed until new independently collected clinical samples are scored with this frozen package.
