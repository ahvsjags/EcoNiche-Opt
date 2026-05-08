# Data Directory

Raw public data, processed expression matrices, metadata, biological priors, and audit manifests are written here by the pipeline.

- `raw/`: downloaded public source files and access-status placeholders for restricted resources.
- `external/`: manifest layer for external public resources; large files remain in `raw/`.
- `interim/`: preprocessing-stage manifests and QC snapshots that point to metadata, priors, and processed matrices.
- `metadata/`: harmonized sample metadata and manual-curation reports.
- `priors/`: gene universe, cell-state priors, network, and ligand-receptor resources.
- `processed/`: model-ready expression matrices and sample metadata.

Do not commit controlled-access data. Mark controlled or unclear datasets in `config/data_registry.yml` and keep local access credentials outside this repository.
