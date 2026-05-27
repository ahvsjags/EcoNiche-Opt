# Running Public Data

1. Audit data access:

```bash
python -m econiche_opt.cli audit-dataset --registry config/data_registry.yml
```

2. Dry-run downloads:

```bash
python -m econiche_opt.cli download-data --dry-run
```

3. Run the real-data orchestration entrypoint without overwriting existing results:

```bash
python -m econiche_opt.cli run-real-pipeline
```

This writes `results/real_pipeline/pipeline_run_manifest.tsv` and records which stages passed, were skipped, or remain pending because data access or manual curation is required.

4. Download public GEO metadata or files when you intentionally want to refresh raw inputs:

```bash
python scripts/download/download_geo.py --registry config/data_registry.yml --metadata-only
python scripts/download/download_geo.py --registry config/data_registry.yml --download-matrix --download-supplementary
```

5. Preprocess, train, benchmark, and report after raw/metadata QC:

```bash
python -m econiche_opt.cli run-real-pipeline --execute-preprocess --execute-training
```

6. Download or audit public TCGA SKCM reference data through Xena/GDC:

```bash
python -m econiche_opt.cli download-xena
python -m econiche_opt.cli download-xena --download --strict --max-download-mb 200
```

The downloader writes `results/real/xena_tcga_manifest.tsv`, downloads TCGA SKCM Xena expression/clinical files when requested, and records a GDC API manifest for STAR-count expression files.

7. Run marker-based deconvolution/cell-abundance baselines:

```bash
python -m econiche_opt.cli run-deconvolution --out-dir results/real
```

This writes `results/real/deconvolution_scores.tsv`, `results/real/deconvolution_cohort_summary.tsv`, and `results/real/deconvolution_baseline_metrics.tsv`.

8. Run perturbation reversal prioritization:

```bash
python -m econiche_opt.cli run-perturbation
```

This queries Enrichr LINCS L1000 up/down perturbation libraries, annotates module genes through DGIdb when available, and writes hypothesis-only rankings under `results/perturbation/`.

9. Run real-data ablation and sensitivity reruns only after non-demo processed cohorts are curated:

```bash
python scripts/analysis/run_ablation.py --execute-real --out results/real_ablation
python scripts/analysis/run_sensitivity.py --execute-real --out results/real_sensitivity
```

Without `--execute-real`, these scripts record `RESULT_PENDING` instead of accidentally using demo processed inputs.

Controlled or licensed cohorts must remain `ACCESS_RESTRICTED` until approved files are available locally.
