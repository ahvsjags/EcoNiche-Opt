# Interim Data Layer

This directory records preprocessing-stage manifests and QC snapshots.

The canonical intermediate data products are stored in `data/metadata/`, `data/priors/`, and `data/processed/` because the current pipeline writes directly from raw inputs to harmonized metadata, priors, and processed matrices. The manifest files here point to those intermediate artifacts so the handoff is auditable without duplicating large matrices.
