# Heuristic Ecology Optimizer Delivery Audit

Date: 2026-05-07

## Implementation Summary

Implemented the remaining Word-spec optimization layer as a registered pipeline component:

- Training-only ecological candidate gene pools.
- Signed rank-normalized module scores.
- Heuristic module optimization with mutation, crossover, and network-neighborhood jumps.
- Training-only interaction-edge search over curated and coexpression edges.
- Biological objective terms and penalties.
- Conservative prediction backbone locking to prevent small-cohort optimizer overfitting.
- GPU-audited optimizer search backend using `xgboost_cuda` when requested.
- Fold-level module, edge, and optimizer-history outputs.

## Key Files Changed

- `src/econiche_opt/model/ecology_optimizer.py`
- `src/econiche_opt/model/endpoint_modules.py`
- `scripts/model/run_endpoint_module_analysis.py`
- `tests/test_endpoint_modules.py`
- `docs/goal_status.yml`

## Key Outputs

- `results/endpoint_modules_heuristic_core_locked_gpu/endpoint_module_summary.tsv`
- `results/endpoint_modules_heuristic_core_locked_gpu/endpoint_module_pairwise_comparisons.tsv`
- `results/endpoint_modules_heuristic_core_locked_gpu/melanoma_primary_rescue_baselines.tsv`
- `results/endpoint_modules_heuristic_core_locked_gpu/optimized_ecology_module.tsv`
- `results/endpoint_modules_heuristic_core_locked_gpu/optimized_ecology_edges.tsv`
- `results/endpoint_modules_heuristic_core_locked_gpu/optimized_ecology_history.tsv`
- `results/endpoint_modules_heuristic_core_locked_gpu/OPTIMIZED_ECOLOGY_AUDIT.md`

## Primary Melanoma Results

`EcoNiche-Opt-HeuristicEcology` uses a locked ecological module-prior backbone for prediction and runs the heuristic optimizer as a training-only biological optimization/audit layer.

- Primary RECIST / melanoma core high-evidence: AUROC 0.705, mean fold AUROC 0.712, ECE 0.235.
- Primary RECIST / melanoma RECIST-supported primary: AUROC 0.685, mean fold AUROC 0.664, ECE 0.239.

The model exceeds eight existing signature baselines by point estimate in both primary melanoma strata. FDR-supported gains are present for CYT and PDCD1LG2 in these focused outputs, while IFNG, CXCL9, TIG, TIDE_dysfunction, APM, IPRES, and TIDE_exclusion remain point-estimate improvements unless their FDR rows pass the claim gate.

## Component-Ablation Boundary

The optimized WordFullGraph component is implemented and audited, but the full interaction graph is not yet a stable predictor by itself in the focused GPU run. The manuscript should therefore claim:

- The ecological optimization layer is implemented and produces interpretable fold-specific modules/edges without holdout leakage.
- The locked predictive backbone preserves the strongest melanoma performance.
- Word graph components are mechanistic/optimization evidence unless the ablation table shows positive paired FDR-supported gain.

Do not claim that every Word graph component independently improves predictive AUROC.

## Verification

Commands completed:

- `python -m pytest tests\test_endpoint_modules.py -q`
- `python -m pytest -q`
- `python -m econiche_opt.cli validate-goals --goal-file docs/goal_status.yml`
- `python -m econiche_opt.cli validate-project --mode demo`
- `python scripts\model\run_endpoint_module_analysis.py --processed-dir data/processed/bulk --out results/endpoint_modules_heuristic_core_locked_gpu --only-endpoint primary_recist --only-stratum melanoma_core_high_evidence --only-stratum melanoma_recist_supported_primary --optimizer-population 2 --optimizer-generations 1 --optimizer-n-jobs 1 --optimizer-use-gpu --optimizer-scope all`

All tests and validators passed.
