# EcoNiche-Opt Agent Instructions

## Project Goal

Build a reproducible EcoNiche-Opt repository for multicohort ICB response benchmarking, six-state ecological module optimization, validation, figures, tables, manuscript skeletons, and reproducibility reporting.

## Non-Negotiable Policies

- No fabrication: real-data outputs must come from registered pipelines. Missing outputs are `RESULT_PENDING`.
- Controlled or licensed data are `ACCESS_RESTRICTED`; provide instructions and interfaces, not substitute data.
- All training, feature selection, thresholding, calibration, and model selection must use patient-level training data only.
- Locked external or holdout cohorts must not leak into training, priors, thresholds, calibration, or claim wording.
- Superiority claims require paired bootstrap or DeLong evidence with FDR support through `econiche_opt.reporting.claim_gate`.

## Common Commands

```bash
python -m pip install -e .
python -m econiche_opt.cli make-demo
python -m pytest -q
python -m econiche_opt.cli validate-goals --goal-file docs/goal_status.yml
python -m econiche_opt.cli validate-project --mode demo
```

## Goal Tracking

Update `docs/goal_status.yml` whenever a GOAL changes. Each GOAL entry records status, files changed, tests run, notes, and blocking issues. Use `interface_completed` for real-data tasks whose code path exists but cannot fully run because access, R/Bioconductor, or an external service is unavailable.

## Deliverables

Generate figures, tables, manuscript, and reproducibility reports from scripts or CLI commands. Do not place hand-edited results in `results/` as evidence.
