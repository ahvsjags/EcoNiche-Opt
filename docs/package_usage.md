# EcoNiche-Opt package usage

This document describes the public API for using EcoNiche-Opt as a package.
Expression matrices must have samples as rows and gene symbols as columns.
Labels are encoded as `1=response/responder` and
`0=non-response/non-responder`.

## Python

```python
from econiche_opt import EcoNicheOptClassifier, load_demo_multicohort

demo = load_demo_multicohort(n_cohorts=3, n_per_cohort=24, random_state=42)

model = EcoNicheOptClassifier(mode="word_full_graph", random_state=42)
model.fit_multicohort(
    demo["X_by_cohort"],
    demo["y_response_by_cohort"],
    demo["metadata_by_cohort"],
)

cohort = sorted(demo["X_by_cohort"])[0]
scores = model.score_samples(demo["X_by_cohort"][cohort])
module = model.module_table()
edges = model.edge_table()
coverage = model.feature_coverage(demo["X_by_cohort"][cohort])

model.save("econiche_opt_model.joblib")
loaded = EcoNicheOptClassifier.load("econiche_opt_model.joblib")
```

## Heuristic ecology optimizer

```python
from econiche_opt import EcoNicheOptClassifier, load_demo_multicohort
from econiche_opt.model.ecology_optimizer import HeuristicEcologyConfig

demo = load_demo_multicohort(n_cohorts=3, n_per_cohort=24, random_state=42)
cfg = HeuristicEcologyConfig(population_size=20, generations=8, n_jobs=-1)

model = EcoNicheOptClassifier(
    mode="heuristic_ecology",
    optimizer_config=cfg,
    random_state=42,
)
model.fit_multicohort(
    demo["X_by_cohort"],
    demo["y_response_by_cohort"],
    demo["metadata_by_cohort"],
)
```

`n_jobs=-1` uses available CPU threads. If `xgboost` with CUDA support is
installed, the optimizer can use its optional GPU backend in inner scoring;
otherwise it falls back to sklearn logistic regression.

## R

```r
install.packages("reticulate")
# devtools::install("r-package/EcoNicheOpt")

library(EcoNicheOpt)

demo <- econiche_load_demo(n_cohorts = 3, n_per_cohort = 24)
model <- econiche_classifier(mode = "word_full_graph")
econiche_fit_multicohort(
  model,
  demo$X_by_cohort,
  demo$y_response_by_cohort,
  demo$metadata_by_cohort
)

scores <- econiche_score(model, demo$X_by_cohort[[1]])
module <- econiche_module_table(model)
```

## Guardrails

- Do not fit or calibrate on locked external cohorts.
- Do not pass holdout labels into `fit` or `fit_multicohort`.
- Do not interpret perturbation rankings as treatment recommendations.
- Report real-data performance only from registered benchmark pipelines.
