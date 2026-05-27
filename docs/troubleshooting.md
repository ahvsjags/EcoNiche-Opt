# Troubleshooting

- `Rscript` missing: install R for Windows and refresh PATH. In the verified local environment, `Rscript.exe` was found under `C:\Program Files\R\R-4.6.0\bin`.
- `make` missing on Windows: install GnuWin32 Make or run the equivalent `python -m econiche_opt.cli ...` commands. In the verified local environment, `make.exe` was found under `C:\Program Files (x86)\GnuWin32\bin`.
- `snakemake` missing: install `snakemake>=9.20` with pip, then run `snakemake -n -s workflow/Snakefile --configfile workflow/config/demo.yml`.
- Registry validation fails: add a cohort entry with access status, downloader, preprocessing script, role, endpoint, and notes.
- Real-data result missing: write `RESULT_PENDING` and document the reason rather than inventing values.
