library(EcoNicheOpt)

demo <- econiche_load_demo(n_cohorts = 3, n_per_cohort = 24, random_state = 42)

model <- econiche_classifier(mode = "word_full_graph", random_state = 42)
econiche_fit_multicohort(
  model,
  demo$X_by_cohort,
  demo$y_response_by_cohort,
  demo$metadata_by_cohort
)

scores <- econiche_score(model, demo$X_by_cohort[[1]])
coverage <- econiche_feature_coverage(model, demo$X_by_cohort[[1]])
module_table <- econiche_module_table(model)

print(head(scores))
print(head(coverage))
print(head(module_table))
