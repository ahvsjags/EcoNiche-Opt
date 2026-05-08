# Word-Full EcoNiche Graph Ablation

This report tests whether the Word-spec components add value over ablations: removing interaction edges, removing signed gene directions, or removing the biological objective from inner model selection.

## Full Model Snapshot

- primary_recist / melanoma_core_high_evidence: n=117, pooled AUROC=0.705, mean fold AUROC=0.712, pooled ECE=0.235.
- primary_recist / melanoma_recist_supported_primary: n=131, pooled AUROC=0.685, mean fold AUROC=0.664, pooled ECE=0.239.

## Component Ablations

- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.705, ablation AUROC=0.609, delta=0.096, FDR q=0.137 (point_estimate_component_gain).
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.705, ablation AUROC=0.545, delta=0.160, FDR q=0.036 (FDR_supported_component_gain).
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.705, ablation AUROC=0.549, delta=0.156, FDR q=0.041 (FDR_supported_component_gain).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.685, ablation AUROC=0.548, delta=0.136, FDR q=0.088 (point_estimate_component_gain).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.685, ablation AUROC=0.550, delta=0.134, FDR q=0.012 (FDR_supported_component_gain).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.685, ablation AUROC=0.536, delta=0.149, FDR q=0.000 (FDR_supported_component_gain).

## Claim Boundary

Treat component gains as supported only when paired bootstrap deltas are positive and FDR-supported in the pre-specified stratum. Point-estimate gains can motivate mechanistic interpretation but should not be worded as statistically proven superiority.
