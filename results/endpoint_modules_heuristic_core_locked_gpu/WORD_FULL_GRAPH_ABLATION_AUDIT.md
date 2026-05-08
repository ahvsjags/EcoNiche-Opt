# Word-Full EcoNiche Graph Ablation

This report tests whether the Word-spec components add value over ablations: removing interaction edges, removing signed gene directions, or removing the biological objective from inner model selection.

## Full Model Snapshot

- primary_recist / melanoma_core_high_evidence: n=117, pooled AUROC=0.705, mean fold AUROC=0.712, pooled ECE=0.235.
- primary_recist / melanoma_recist_supported_primary: n=131, pooled AUROC=0.685, mean fold AUROC=0.664, pooled ECE=0.239.

## Component Ablations

- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.519, ablation AUROC=0.467, delta=0.053, FDR q=0.478 (point_estimate_component_gain).
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.519, ablation AUROC=0.503, delta=0.017, FDR q=0.762 (point_estimate_component_gain).
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.519, ablation AUROC=0.488, delta=0.032, FDR q=0.608 (point_estimate_component_gain).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.478, ablation AUROC=0.429, delta=0.049, FDR q=0.355 (point_estimate_component_gain).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.478, ablation AUROC=0.585, delta=-0.106, FDR q=0.030 (component_not_supported).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.478, ablation AUROC=0.500, delta=-0.022, FDR q=0.712 (component_not_supported).

## Claim Boundary

Treat component gains as supported only when paired bootstrap deltas are positive and FDR-supported in the pre-specified stratum. Point-estimate gains can motivate mechanistic interpretation but should not be worded as statistically proven superiority.
