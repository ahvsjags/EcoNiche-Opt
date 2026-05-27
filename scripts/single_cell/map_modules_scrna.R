status <- system2("python", c("scripts/single_cell/map_modules_scrna.py"))
if (status != 0) {
  stop("Python scRNA module mapping failed")
}
