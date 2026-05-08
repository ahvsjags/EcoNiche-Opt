# Word-Full EcoNiche Graph Ablation

This report tests whether the Word-spec components add value over ablations: removing interaction edges, removing signed gene directions, or removing the biological objective from inner model selection.

## Full Model Snapshot

- primary_recist / melanoma_core_high_evidence: n=117, pooled AUROC=0.687, mean fold AUROC=0.642, pooled ECE=0.285.
- primary_recist / melanoma_recist_supported_primary: n=131, pooled AUROC=0.521, mean fold AUROC=0.537, pooled ECE=0.217.

## Component Ablations

- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.687, ablation AUROC=0.355, delta=0.332, FDR q=0.018 (FDR_supported_component_gain).
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.687, ablation AUROC=0.701, delta=-0.014, FDR q=0.550 (component_not_supported).
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.687, ablation AUROC=0.471, delta=0.215, FDR q=0.000 (FDR_supported_component_gain).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.521, ablation AUROC=0.575, delta=-0.054, FDR q=0.414 (component_not_supported).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.521, ablation AUROC=0.506, delta=0.015, FDR q=0.784 (point_estimate_component_gain).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.521, ablation AUROC=0.453, delta=0.067, FDR q=0.288 (point_estimate_component_gain).

## Claim Boundary

Treat component gains as supported only when paired bootstrap deltas are positive and FDR-supported in the pre-specified stratum. Point-estimate gains can motivate mechanistic interpretation but should not be worded as statistically proven superiority.
