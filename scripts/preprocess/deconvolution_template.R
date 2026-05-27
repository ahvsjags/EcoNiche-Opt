args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- normalizePath(sub(file_arg, "", args[grepl(file_arg, args)]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
runner <- file.path(root, "scripts", "preprocess", "run_deconvolution.py")
status <- system2("python", c(runner), stdout = TRUE, stderr = TRUE)
cat(paste(status, collapse = "\n"), "\n")
