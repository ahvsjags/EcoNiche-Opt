.as_expression_frame <- function(expression) {
  frame <- as.data.frame(expression, check.names = FALSE)
  if (is.null(rownames(frame))) {
    stop("expression must have sample IDs as row names", call. = FALSE)
  }
  frame
}

econiche_python <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("The reticulate package is required.", call. = FALSE)
  }
  reticulate::import("econiche_opt", delay_load = FALSE)
}

install_econiche_python <- function(repo = NULL, envname = "r-econiche-opt") {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Install reticulate before installing the Python package.", call. = FALSE)
  }
  package <- if (is.null(repo)) "econiche-opt" else paste0("git+", repo)
  reticulate::py_install(package, envname = envname, pip = TRUE)
  invisible(TRUE)
}

econiche_classifier <- function(mode = "word_full_graph",
                                include_interactions = TRUE,
                                signed = TRUE,
                                calibration = NULL,
                                random_state = 42L) {
  py <- econiche_python()
  py$EcoNicheOptClassifier(
    mode = mode,
    include_interactions = include_interactions,
    signed = signed,
    calibration = calibration,
    random_state = as.integer(random_state)
  )
}

econiche_fit <- function(model, expression, labels) {
  expression <- .as_expression_frame(expression)
  model$fit(expression, as.integer(labels))
  invisible(model)
}

econiche_fit_multicohort <- function(model,
                                     expression_by_cohort,
                                     labels_by_cohort,
                                     metadata_by_cohort = NULL) {
  model$fit_multicohort(expression_by_cohort, labels_by_cohort, metadata_by_cohort)
  invisible(model)
}

econiche_score <- function(model, expression) {
  expression <- .as_expression_frame(expression)
  as.data.frame(model$score_samples(expression))
}

econiche_feature_coverage <- function(model, expression) {
  expression <- .as_expression_frame(expression)
  as.data.frame(model$feature_coverage(expression))
}

econiche_module_table <- function(model) {
  as.data.frame(model$module_table())
}

econiche_save <- function(model, path) {
  model$save(path)
  invisible(normalizePath(path, mustWork = FALSE))
}

econiche_load <- function(path) {
  py <- econiche_python()
  py$EcoNicheOptClassifier$load(path)
}

econiche_load_demo <- function(n_cohorts = 3L, n_per_cohort = 24L, random_state = 42L) {
  py <- econiche_python()
  py$load_demo_multicohort(
    n_cohorts = as.integer(n_cohorts),
    n_per_cohort = as.integer(n_per_cohort),
    random_state = as.integer(random_state)
  )
}
