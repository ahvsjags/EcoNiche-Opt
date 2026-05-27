# License, Citation, and Data-Use Notes

This repository is MIT licensed for code. Public datasets remain governed by their original sources and citations.

Controlled or licensed datasets, including optional DrugBank-style resources and dbGaP/EGA-like sources, require independent approval. The repository stores no credentials and does not redistribute controlled raw data.

All source records are tracked in `config/source_registry.yml`. Run:

```bash
python -m econiche_opt.cli validate-sources --source-registry config/source_registry.yml
```
