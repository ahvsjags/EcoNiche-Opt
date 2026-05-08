# GitHub release checklist

This repository contains public package code and local analysis outputs.
Before publishing to GitHub, use a clean repository that includes source code,
tests, package metadata, documentation, and small demo interfaces, but excludes
raw data, processed public-data matrices, controlled-access files, large
result artifacts, and manuscript working binaries.

## Safe publish scope

Include:

- `src/`
- `tests/`
- `config/`
- `schemas/`
- `workflow/`
- `.github/workflows/test.yml`
- `r-package/EcoNicheOpt/`
- `docs/`
- `scripts/` that implement registered pipelines
- `README.md`, `LICENSE`, `CITATION.cff`, `pyproject.toml`, `requirements.txt`, `environment.yml`, `MANIFEST.in`, `AGENTS.md`

Exclude:

- `data/raw/`
- `data/processed/`
- `data/interim/` except README placeholders
- `data/external/` except README placeholders
- `results/`
- `deliverables/*.zip`
- `tmp/`
- `.snakemake/`
- rendered Word/PDF manuscript binaries unless intentionally released

## Commands when GitHub CLI is available

```bash
git init
git add README.md LICENSE CITATION.cff pyproject.toml requirements.txt environment.yml MANIFEST.in AGENTS.md
git add src tests config schemas workflow docs scripts r-package .github
git commit -m "package EcoNiche-Opt public API"
gh repo create EcoNiche-Opt --public --source . --remote origin --push
```

If the repository already exists:

```bash
git remote add origin https://github.com/<OWNER>/EcoNiche-Opt.git
git push -u origin main
```
