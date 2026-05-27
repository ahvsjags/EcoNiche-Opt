args <- commandArgs(trailingOnly = TRUE)
script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- normalizePath(sub(file_arg, "", script_args[grepl(file_arg, script_args)]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
runner <- file.path(root, "scripts", "download", "download_xena.py")
status <- system2("python", c(runner, args), stdout = TRUE, stderr = TRUE)
cat(paste(status, collapse = "\n"), "\n")
