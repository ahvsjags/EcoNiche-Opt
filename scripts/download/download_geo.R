args <- commandArgs(trailingOnly = TRUE)
registry <- if (length(args) >= 1) args[[1]] else "config/data_registry.yml"
message("GEO R downloader placeholder. Prefer scripts/download/download_geo.py for metadata extraction.")
message("Registry: ", registry)
dir.create("data/metadata", recursive = TRUE, showWarnings = FALSE)
write.table(
  data.frame(status = "RESULT_PENDING", reason = "Run python scripts/download/download_geo.py --metadata-only"),
  file = "data/metadata/geo_download_R_placeholder.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
