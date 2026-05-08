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

Excluded:

- `data/raw/`
- `data/external/` except README placeholders
- `*.expr.tsv` expression matrices in `data/processed/` and result `processed_inputs/`
- `figures/**/*.tiff`
- nested `.git` directories and local caches

These exclusions keep the GitHub repository usable and avoid redistributing
controlled, licensed, or oversized data. Regeneration and download paths remain
documented in the package and pipeline documentation.
