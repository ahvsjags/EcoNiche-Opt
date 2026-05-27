# Public data/results/figures upload notes

This repository includes public reproducibility artifacts from the EcoNiche-Opt
workspace, but intentionally does not include raw public downloads, controlled
or access-restricted data, nested third-party Git checkouts, or oversized TIFF
exports.

Included:

- `data/metadata/`: sample metadata, manual curation tables, and file manifests.
- `data/priors/`: gene universe, cell-state priors, and lightweight network resources.
- `data/interim/`: audit and pipeline manifests.
- `data/processed/`: metadata files only, not expression matrices.
- `results/`: generated benchmark, validation, claim-gate, perturbation, and audit tables, excluding duplicated processed expression inputs.
- `figures/article/`: PNG/PDF/SVG manuscript figure exports.

Authoritative current result directories for the submitted manuscript include
`results/endpoint_modules_heuristic_core_locked_gpu/`,
`results/aligned_panel_ablation_20260527/`, and
`results/locked_external_panel_validation_calibrated_20260519/`,
`results/gpu_bioprior_rescue_combo_search_robust_20260527/`,
`results/cbioportal_gpu_bioprior_external_validation_20260527/`, and
`results/gpu_lipid_pair_rescue_component_dominant_ba_guardrail_20260528/`.
The frozen package-level evidence for the current promoted external-rescue
model is in `deliverables/gpu_lipid_pair_rescue_package_20260528/`.
Older un-calibrated validation folders and no-leakage negative optimization
audits are retained as historical pipeline outputs and are not the source for
the current manuscript performance claims.

Excluded:

- `data/raw/`
- `data/external/` except README placeholders
- `*.expr.tsv` expression matrices in `data/processed/` and result `processed_inputs/`
- `figures/**/*.tiff`
- nested `.git` directories and local caches

These exclusions keep the GitHub repository usable and avoid redistributing
controlled, licensed, or oversized data. Regeneration and download paths remain
documented in the package and pipeline documentation.

Release status: the current local manuscript/package version is
`v0.3.4-gpu-lipid-pair-rescue-20260528`. Zenodo metadata are prepared with
`zenodo_doi=RESULT_PENDING`; no DOI should be cited until Zenodo or an
institutional archive mints a real identifier.
