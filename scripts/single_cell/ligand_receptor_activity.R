dir.create("results/scrna", recursive = TRUE, showWarnings = FALSE)
write.table(
  data.frame(ligand = "RESULT_PENDING", receptor = "RESULT_PENDING", activity = NA, status = "RESULT_PENDING"),
  file = "results/scrna/ligand_receptor_activity.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
message("Wrote results/scrna/ligand_receptor_activity.tsv")
