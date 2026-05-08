# Optimized Ecological Module Audit

The WordFullGraph model now optimizes state genes and interaction edges inside each training fold using only training cohorts. The optimizer uses signed rank module scores, curated priors, mutation/crossover, network-neighborhood jumps, biological objective terms, and component ablations.

## Output Files

- `optimized_ecology_module.tsv`: fold-specific selected genes, directions, training correlations, and selection frequencies.
- `optimized_ecology_edges.tsv`: fold-specific curated and coexpression interaction edges.
- `optimized_ecology_history.tsv`: generation-level optimizer diagnostics and compute backend.

## Target Snapshot

- primary_recist / melanoma_core_high_evidence: n=117, pooled AUROC=0.600, mean fold AUROC=0.579, ECE=0.341.
- primary_recist / melanoma_recist_supported_primary: n=131, pooled AUROC=0.621, mean fold AUROC=0.595, ECE=0.237.

## Compute Backend

- Backends observed: sklearn_logistic_cpu

## Claim Boundary

Optimizer gains are claimable only through the paired baseline and ablation outputs generated in this same run. Do not describe fold-specific selected genes as locked clinical markers until a final training-only lockbox/panel analysis is run.
