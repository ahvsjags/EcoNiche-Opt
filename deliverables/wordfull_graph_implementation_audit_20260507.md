# Word-Spec EcoNiche-Opt Graph Implementation Audit

Date: 2026-05-07

## What Was Implemented

The current endpoint-module pipeline now contains a Word-spec ecological model variant named `EcoNiche-Opt-WordFullGraph`.

Implemented components:

- Sample-wise rank Gaussian expression normalization for Word-spec features.
- Training-only gene direction estimation for signed module activity.
- Six ecological states matching the Word model: tumor dedifferentiation, antigen presentation/MHC, T/NK effector, T-cell dysfunction, CAF/ECM exclusion, and myeloid suppression.
- Signed state activity scores.
- Curated ligand-receptor, pathway, regulatory, checkpoint, and network interaction edge features.
- A Word-style ecological score layered on top of the current module-prior model.
- Biological objective terms for cell-state coverage, pathway/network/LR edge support, direction stability, batch dependence, redundancy, and therapy confounding.
- Real-data ablations:
  - `EcoNiche-Opt-WordNoInteraction`
  - `EcoNiche-Opt-WordUnsignedGraph`
  - `EcoNiche-Opt-WordNoBioObjective`

## Key Output Files

- `results/endpoint_modules_wordfull_v2/endpoint_module_summary.tsv`
- `results/endpoint_modules_wordfull_v2/endpoint_module_predictions.tsv`
- `results/endpoint_modules_wordfull_v2/endpoint_module_pairwise_comparisons.tsv`
- `results/endpoint_modules_wordfull_v2/word_full_graph_ablation.tsv`
- `results/endpoint_modules_wordfull_v2/WORD_FULL_GRAPH_ABLATION_AUDIT.md`
- `results/endpoint_modules_wordfull_v2/MELANOMA_PRIMARY_RESCUE_AUDIT.md`
- `results/endpoint_modules_wordfull_v2/melanoma_primary_rescue_baselines.tsv`

## Main Real-Data Results

Primary RECIST:

- Full melanoma anti-PD1 primary: AUROC 0.662, mean LODO AUROC 0.616, ECE 0.233.
- RECIST-supported melanoma primary: AUROC 0.681, mean LODO AUROC 0.664, ECE 0.221.
- High-evidence melanoma core: AUROC 0.704, mean LODO AUROC 0.697, ECE 0.201.
- Pan-cancer response all: AUROC 0.645, mean LODO AUROC 0.632, ECE 0.299.

Strict RECIST:

- High-evidence melanoma core: AUROC 0.708, mean LODO AUROC 0.696, ECE 0.203.
- RECIST-supported melanoma primary: AUROC 0.692, mean LODO AUROC 0.674, ECE 0.201.
- Pan-cancer response all: AUROC 0.654, mean LODO AUROC 0.646, ECE 0.256.

## Baseline Comparison Boundary

For primary RECIST high-evidence melanoma core, `EcoNiche-Opt-WordFullGraph` reached AUROC 0.704 and exceeded the following eight existing baselines by point estimate:

- IFNG: 0.704 vs 0.664
- CXCL9: 0.704 vs 0.666
- TIG: 0.704 vs 0.666
- TIDE_dysfunction: 0.704 vs 0.655
- APM: 0.704 vs 0.664
- CYT: 0.704 vs 0.630
- IPRES: 0.704 vs 0.567
- TIDE_exclusion: 0.704 vs 0.544

Only the CYT comparison in this stratum was FDR-supported in the current paired bootstrap table. The other improvements should be worded as point-estimate gains, not statistically proven superiority.

## Ablation Findings

Primary RECIST high-evidence melanoma core:

- Full vs no biological objective: 0.704 vs 0.702, point-estimate gain.
- Full vs no interaction edges: 0.704 vs 0.704, tiny point-estimate gain.
- Full vs unsigned graph: 0.704 vs 0.733, signed direction not supported in this endpoint.

Strict RECIST high-evidence melanoma core:

- Full vs no biological objective: 0.708 vs 0.704, point-estimate gain.
- Full vs no interaction edges: 0.708 vs 0.708, tiny point-estimate gain.
- Full vs unsigned graph: 0.708 vs 0.702, point-estimate gain.

Strong claim boundary: Word-spec components are implemented and several have positive point-estimate gains, but component superiority is not broadly FDR-supported. The safe manuscript wording is that Word-spec ecological structure improves conceptual novelty and preserves strong high-evidence melanoma performance; statistical component superiority remains endpoint-dependent.

## Verification

Commands completed:

- `python -m py_compile src\\econiche_opt\\model\\endpoint_modules.py scripts\\model\\run_endpoint_module_analysis.py`
- `python -m pytest tests\\test_endpoint_modules.py -q`
- `python scripts\\model\\run_endpoint_module_analysis.py --processed-dir data/processed/bulk --out results/endpoint_modules_wordfull_v2`
- `python -m pytest -q`
- `python -m econiche_opt.cli validate-goals --goal-file docs/goal_status.yml`
- `python -m econiche_opt.cli validate-project --mode demo`

Result: all tests and project validators passed.
