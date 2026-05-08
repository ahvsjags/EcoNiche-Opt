# Word-Full EcoNiche Graph Ablation

This report tests whether the Word-spec components add value over ablations: removing interaction edges, removing signed gene directions, or removing the biological objective from inner model selection.

## Full Model Snapshot

- primary_recist / melanoma_core_high_evidence: n=117, pooled AUROC=0.585, mean fold AUROC=0.466, pooled ECE=0.277.
- primary_recist / melanoma_recist_supported_primary: n=131, pooled AUROC=0.373, mean fold AUROC=0.345, pooled ECE=0.314.

## Component Ablations

- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.585, ablation AUROC=0.521, delta=0.064, FDR q=0.399 (point_estimate_component_gain).
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.585, ablation AUROC=0.555, delta=0.030, FDR q=0.540 (point_estimate_component_gain).
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.585, ablation AUROC=0.431, delta=0.154, FDR q=0.399 (point_estimate_component_gain).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.373, ablation AUROC=0.504, delta=-0.131, FDR q=0.023 (component_not_supported).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.373, ablation AUROC=0.467, delta=-0.095, FDR q=0.009 (component_not_supported).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.373, ablation AUROC=0.503, delta=-0.130, FDR q=0.066 (component_not_supported).

## Claim Boundary

Treat component gains as supported only when paired bootstrap deltas are positive and FDR-supported in the pre-specified stratum. Point-estimate gains can motivate mechanistic interpretation but should not be worded as statistically proven superiority.
