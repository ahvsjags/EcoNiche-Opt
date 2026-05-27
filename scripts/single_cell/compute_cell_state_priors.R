dir.create("data/priors", recursive = TRUE, showWarnings = FALSE)
write.table(
  data.frame(status = "RESULT_PENDING", reason = "Use curated scRNA markers or public annotations"),
  file = "data/priors/scrna_prior_status.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
message("Wrote data/priors/scrna_prior_status.tsv")
