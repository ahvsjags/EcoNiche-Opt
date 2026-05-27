dir.create("results/scrna", recursive = TRUE, showWarnings = FALSE)
write.table(
  data.frame(status = "RESULT_PENDING", reason = "Configure scRNA count matrix and annotations"),
  file = "results/scrna/preprocess_status.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
message("Wrote results/scrna/preprocess_status.tsv")
