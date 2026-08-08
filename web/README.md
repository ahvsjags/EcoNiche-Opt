# EcoNiche-Opt research portal

This is a static, GitHub Pages-compatible research interface for the EcoNiche-Opt repository. It provides:

- a generated view of the locked 62-gene validation specification;
- browser-side expression-matrix format and panel-coverage auditing;
- a label-free locked-score preview using the registered module-prior and frozen endpoint calibration;
- sample-level score and audit exports;
- module, endpoint and reproducibility views for model users;
- direct downloads of the locked specification and panel mapping.

The browser preview is not a replacement for the formal independent-validation command. For production analysis use:

```bash
econiche-opt score-locked-validation \
  --package-dir deliverables/prospective_validation \
  --expression independent_expression.tsv \
  --sample-manifest assay_sample_manifest.tsv \
  --clinical-annotation clinical_annotation.tsv \
  --out-dir results/independent_locked_validation
```

## Build the data bundle

From the repository root:

```bash
python scripts/web/build_portal_data.py
```

The generated `web/data/portal_manifest.json` is derived from registered repository files. It must not be hand-edited to change scientific results.

## Run locally

```bash
python -m http.server 4173 --directory web
```

Open <http://127.0.0.1:4173>.

## GitHub Pages

The workflow at `.github/workflows/deploy-pages.yml` regenerates the portal data and publishes the `web/` directory on pushes to `main`. In the repository settings, set **Pages → Build and deployment → Source** to **GitHub Actions**.

The portal uses Lucide icons from the public CDN when online. All scoring and audit logic remains local to the browser; no expression matrix is uploaded by this static page.
