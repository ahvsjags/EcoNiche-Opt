# EcoNiche-Opt

EcoNiche-Opt is a reproducible Python package and analysis repository for multicohort transcriptomic immune-checkpoint blockade (ICB) response benchmarking and six-state ecological module modeling. It exposes the model as a direct-call API for downstream users while keeping the full manuscript benchmark, figures, tables, and reproducibility checks in registered pipelines.

This repository intentionally starts with a synthetic/demo pipeline. Real public cohorts are routed through `config/data_registry.yml` and downloader/preprocessing entrypoints. Controlled or unclear access datasets are marked as restricted or unknown and are not silently treated as available.

## No Fabrication Policy

Real-data claims must come from registered download, preprocessing, model, benchmark, and reporting commands in this repository. Missing real outputs are written as `RESULT_PENDING`; controlled or licensed datasets are written as `ACCESS_RESTRICTED` with instructions. The manuscript and claim gate must not state superiority, best-model status, or clinical actionability unless paired bootstrap or DeLong evidence with FDR support is present.

## Quick Start

```bash
python -m pip install -e .
python - <<'PY'
from econiche_opt import EcoNicheOptClassifier, load_demo_multicohort

demo = load_demo_multicohort(n_cohorts=3, n_per_cohort=24, random_state=42)
model = EcoNicheOptClassifier(mode="word_full_graph", random_state=42)
model.fit_multicohort(
    demo["X_by_cohort"],
    demo["y_response_by_cohort"],
    demo["metadata_by_cohort"],
)
cohort = sorted(demo["X_by_cohort"])[0]
print(model.score_samples(demo["X_by_cohort"][cohort]).head())
print(model.module_table().head())
PY
```

For the full demo pipeline:

```bash
python scripts/make_demo_data.py
python scripts/model/run_econiche.py --config config/model_config.yml --demo
python -m econiche_opt.cli validate-project --mode demo
python -m pytest -q
```

## Public Package API

EcoNiche-Opt now has a stable public package interface:

- `EcoNicheOptClassifier(mode="word_full_graph")`: locked six-state signed-rank module model with ecological interaction edges.
- `EcoNicheOptClassifier(mode="heuristic_ecology")`: training-only heuristic optimizer that searches module genes and interaction edges under biological constraints.
- `fit_multicohort(X_by_cohort, y_by_cohort, metadata_by_cohort)`: multicohort fitting with labels encoded as `1=response/responder` and `0=non-response/non-responder`.
- `score_samples(X)`: tidy per-sample response probabilities, predicted labels, decision scores, and locked threshold.
- `module_table()`, `edge_table()`, and `feature_coverage(X)`: transparent model internals for assay translation and audit.
- `save(path)` and `EcoNicheOptClassifier.load(path)`: model serialization for deployment.

The public API is intentionally label-safe: feature directions, module selection, thresholds, and calibration are fit only on the supplied training data. External or holdout labels must not be passed into fitting.

## Command Line Use

Users who do not want to write Python can train and score from TSV/CSV files:

```bash
econiche-opt fit-package-model \
  --expression expression.tsv \
  --labels labels.tsv \
  --model-out econiche_model.joblib \
  --out-dir econiche_model_artifacts \
  --mode word_full_graph

econiche-opt score-package-model \
  --model econiche_model.joblib \
  --expression new_expression.tsv \
  --out econiche_scores.tsv \
  --coverage-out econiche_feature_coverage.tsv
```

`expression.tsv` should have samples as rows and gene symbols as columns. Use
`--transpose` if genes are rows and samples are columns. `labels.tsv` should
contain `sample_id`, `response_label`, and optionally `cohort`.

## R Interface

An R wrapper is provided in `r-package/EcoNicheOpt` through `reticulate`:

```r
install.packages("reticulate")
# devtools::install("r-package/EcoNicheOpt")

library(EcoNicheOpt)
demo <- econiche_load_demo(n_cohorts = 3, n_per_cohort = 24)
model <- econiche_classifier(mode = "word_full_graph")
econiche_fit_multicohort(model, demo$X_by_cohort, demo$y_response_by_cohort, demo$metadata_by_cohort)
scores <- econiche_score(model, demo$X_by_cohort[[1]])
```

The R package is a thin wrapper around the Python implementation, so the Python package must be installed in the active `reticulate` environment.

On systems with GNU Make:

```bash
make demo
make test
make benchmark-demo
make report-demo
make registry-audit
```

## Core Outputs

- `results/demo/econiche_module.tsv`
- `results/demo/lodo_metrics.tsv`
- `results/demo/lodo_predictions.tsv`
- `results/demo/objective_terms.json`
- `results/demo/objective_history.tsv`
- `tables/dataset_access_audit.tsv`
- `paper/manuscript.md`

## Real Data Workflow

```bash
python scripts/preprocess/audit_registry.py --registry config/data_registry.yml --out tables/dataset_access_audit.tsv
python scripts/download/download_geo.py --registry config/data_registry.yml --metadata-only
python scripts/preprocess/build_metadata.py
python scripts/preprocess/harmonize_labels.py
python scripts/preprocess/deduplicate_patients.py
python scripts/preprocess/preprocess_bulk.py
python scripts/preprocess/build_gene_universe.py
python scripts/preprocess/build_priors.py
python scripts/model/run_econiche.py --config config/model_config.yml
python scripts/baselines/run_baselines.py
python scripts/benchmark/run_lodo.py
python scripts/benchmark/bootstrap_compare.py
python scripts/figures/make_fig1_overview.py
python scripts/paper/generate_methods_text.py
```

## Analysis Guardrails

- Labels are harmonized as `0=responder/sensitive` and `1=non-responder/resistant`.
- All LODO folds estimate gene directions, modules, coefficients, thresholds, and calibration using training cohorts only.
- Patient identifiers are checked for train/test leakage.
- Missing real results are emitted as `RESULT_PENDING`; unavailable baselines are emitted as `unavailable_with_reason`.
- Perturbation outputs are candidate hypotheses only, not treatment recommendations.
