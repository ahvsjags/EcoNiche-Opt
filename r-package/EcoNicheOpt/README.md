# EcoNicheOpt R wrapper

`EcoNicheOpt` is a lightweight R wrapper around the Python `econiche-opt`
package. It is intended for analysts who keep expression matrices in R but
want to call the locked EcoNiche-Opt scoring model or the training-only
heuristic ecology optimizer.

```r
install.packages("reticulate")

# After the Python package is published:
# EcoNicheOpt::install_econiche_python("https://github.com/<OWNER>/EcoNiche-Opt")

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
modules <- econiche_module_table(model)
```

Labels are encoded as `1=response/responder` and
`0=non-response/non-responder`. The wrapper does not download or bundle
restricted cohort data.
