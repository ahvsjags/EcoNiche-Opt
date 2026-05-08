# Word-Full EcoNiche Graph Ablation

This report tests whether the Word-spec components add value over ablations: removing interaction edges, removing signed gene directions, or removing the biological objective from inner model selection.

## Full Model Snapshot

- clinical_benefit / melanoma_anti_pd1_primary: n=160, pooled AUROC=0.594, mean fold AUROC=0.578, pooled ECE=0.267.
- clinical_benefit / melanoma_binary_response_stress: n=29, pooled AUROC=0.460, mean fold AUROC=0.449, pooled ECE=0.397.
- clinical_benefit / melanoma_core_high_evidence: n=117, pooled AUROC=0.633, mean fold AUROC=0.630, pooled ECE=0.268.
- clinical_benefit / melanoma_core_plus_phs000452: n=238, pooled AUROC=0.608, mean fold AUROC=0.616, pooled ECE=0.246.
- clinical_benefit / melanoma_recist_supported_primary: n=131, pooled AUROC=0.634, mean fold AUROC=0.639, pooled ECE=0.245.
- clinical_benefit / pan_cancer_response_all: n=554, pooled AUROC=0.626, mean fold AUROC=0.622, pooled ECE=0.222.
- clinical_benefit / pan_cancer_without_secondary: n=369, pooled AUROC=0.621, mean fold AUROC=0.596, pooled ECE=0.227.
- clinical_benefit / secondary_confounded_transfer: n=185, pooled AUROC=0.632, mean fold AUROC=0.714, pooled ECE=0.228.
- primary_recist / melanoma_anti_pd1_primary: n=160, pooled AUROC=0.662, mean fold AUROC=0.616, pooled ECE=0.233.
- primary_recist / melanoma_binary_response_stress: n=29, pooled AUROC=0.460, mean fold AUROC=0.449, pooled ECE=0.397.
- primary_recist / melanoma_core_high_evidence: n=117, pooled AUROC=0.704, mean fold AUROC=0.697, pooled ECE=0.201.
- primary_recist / melanoma_core_plus_phs000452: n=238, pooled AUROC=0.631, mean fold AUROC=0.658, pooled ECE=0.262.
- primary_recist / melanoma_recist_supported_primary: n=131, pooled AUROC=0.681, mean fold AUROC=0.664, pooled ECE=0.221.
- primary_recist / pan_cancer_response_all: n=554, pooled AUROC=0.645, mean fold AUROC=0.632, pooled ECE=0.299.
- primary_recist / pan_cancer_without_secondary: n=369, pooled AUROC=0.647, mean fold AUROC=0.611, pooled ECE=0.274.
- primary_recist / secondary_confounded_transfer: n=185, pooled AUROC=0.616, mean fold AUROC=0.703, pooled ECE=0.265.
- strict_recist / melanoma_anti_pd1_primary: n=137, pooled AUROC=0.640, mean fold AUROC=0.612, pooled ECE=0.246.
- strict_recist / melanoma_binary_response_stress: n=29, pooled AUROC=0.460, mean fold AUROC=0.449, pooled ECE=0.397.
- strict_recist / melanoma_core_high_evidence: n=95, pooled AUROC=0.708, mean fold AUROC=0.696, pooled ECE=0.203.
- strict_recist / melanoma_core_plus_phs000452: n=200, pooled AUROC=0.637, mean fold AUROC=0.658, pooled ECE=0.252.
- strict_recist / melanoma_recist_supported_primary: n=108, pooled AUROC=0.692, mean fold AUROC=0.674, pooled ECE=0.201.
- strict_recist / pan_cancer_response_all: n=480, pooled AUROC=0.654, mean fold AUROC=0.646, pooled ECE=0.256.
- strict_recist / pan_cancer_without_secondary: n=316, pooled AUROC=0.657, mean fold AUROC=0.628, pooled ECE=0.247.
- strict_recist / secondary_confounded_transfer: n=164, pooled AUROC=0.636, mean fold AUROC=0.720, pooled ECE=0.247.

## Component Ablations

- primary_recist / melanoma_anti_pd1_primary vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.662, ablation AUROC=0.661, delta=0.001, FDR q=0.836 (point_estimate_component_gain).
- primary_recist / melanoma_anti_pd1_primary vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.662, ablation AUROC=0.663, delta=-0.000, FDR q=0.833 (component_not_supported).
- primary_recist / melanoma_anti_pd1_primary vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.662, ablation AUROC=0.665, delta=-0.003, FDR q=0.836 (component_not_supported).
- primary_recist / melanoma_binary_response_stress vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.460, ablation AUROC=0.480, delta=-0.020, FDR q=0.854 (component_not_supported).
- primary_recist / melanoma_binary_response_stress vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.460, ablation AUROC=0.460, delta=0.000, FDR q=1.000 (component_not_supported).
- primary_recist / melanoma_binary_response_stress vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.460, ablation AUROC=0.455, delta=0.005, FDR q=1.000 (point_estimate_component_gain).
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.704, ablation AUROC=0.702, delta=0.002, FDR q=0.747 (point_estimate_component_gain).
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.704, ablation AUROC=0.704, delta=0.000, FDR q=0.908 (point_estimate_component_gain).
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.704, ablation AUROC=0.733, delta=-0.029, FDR q=0.018 (component_not_supported).
- primary_recist / melanoma_core_plus_phs000452 vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.631, ablation AUROC=0.631, delta=0.000, FDR q=0.868 (point_estimate_component_gain).
- primary_recist / melanoma_core_plus_phs000452 vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.631, ablation AUROC=0.631, delta=0.000, FDR q=0.978 (component_not_supported).
- primary_recist / melanoma_core_plus_phs000452 vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.631, ablation AUROC=0.638, delta=-0.007, FDR q=0.116 (component_not_supported).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.681, ablation AUROC=0.682, delta=-0.001, FDR q=0.760 (component_not_supported).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.681, ablation AUROC=0.683, delta=-0.002, FDR q=0.185 (component_not_supported).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.681, ablation AUROC=0.703, delta=-0.023, FDR q=0.115 (component_not_supported).
- primary_recist / pan_cancer_response_all vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.645, ablation AUROC=0.642, delta=0.003, FDR q=0.118 (point_estimate_component_gain).
- primary_recist / pan_cancer_response_all vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.645, ablation AUROC=0.645, delta=-0.000, FDR q=0.234 (component_not_supported).
- primary_recist / pan_cancer_response_all vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.645, ablation AUROC=0.650, delta=-0.005, FDR q=0.118 (component_not_supported).
- primary_recist / pan_cancer_without_secondary vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.647, ablation AUROC=0.643, delta=0.004, FDR q=0.115 (point_estimate_component_gain).
- primary_recist / pan_cancer_without_secondary vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.647, ablation AUROC=0.647, delta=-0.000, FDR q=0.593 (component_not_supported).
- primary_recist / pan_cancer_without_secondary vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.647, ablation AUROC=0.639, delta=0.008, FDR q=0.405 (point_estimate_component_gain).
- primary_recist / secondary_confounded_transfer vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.616, ablation AUROC=0.616, delta=-0.000, FDR q=1.000 (component_not_supported).
- primary_recist / secondary_confounded_transfer vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.616, ablation AUROC=0.616, delta=0.000, FDR q=1.000 (point_estimate_component_gain).
- primary_recist / secondary_confounded_transfer vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.616, ablation AUROC=0.618, delta=-0.002, FDR q=1.000 (component_not_supported).

## Claim Boundary

Treat component gains as supported only when paired bootstrap deltas are positive and FDR-supported in the pre-specified stratum. Point-estimate gains can motivate mechanistic interpretation but should not be worded as statistically proven superiority.
