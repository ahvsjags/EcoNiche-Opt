# Word-Full EcoNiche Graph Ablation

This report tests whether the Word-spec components add value over ablations: removing interaction edges, removing signed gene directions, or removing the biological objective from inner model selection.

## Full Model Snapshot

- primary_recist / melanoma_core_high_evidence: n=117, pooled AUROC=0.600, mean fold AUROC=0.579, pooled ECE=0.341.
- primary_recist / melanoma_recist_supported_primary: n=131, pooled AUROC=0.621, mean fold AUROC=0.595, pooled ECE=0.237.

## Component Ablations

- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.600, ablation AUROC=0.475, delta=0.125, FDR q=0.369 (point_estimate_component_gain).
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.600, ablation AUROC=0.556, delta=0.044, FDR q=0.468 (point_estimate_component_gain).
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.600, ablation AUROC=0.355, delta=0.245, FDR q=0.036 (FDR_supported_component_gain).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.621, ablation AUROC=0.438, delta=0.183, FDR q=0.000 (FDR_supported_component_gain).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.621, ablation AUROC=0.557, delta=0.064, FDR q=0.749 (point_estimate_component_gain).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.621, ablation AUROC=0.445, delta=0.176, FDR q=0.018 (FDR_supported_component_gain).

## Claim Boundary

Treat component gains as supported only when paired bootstrap deltas are positive and FDR-supported in the pre-specified stratum. Point-estimate gains can motivate mechanistic interpretation but should not be worded as statistically proven superiority.
