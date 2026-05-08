# Word-Full EcoNiche Graph Ablation

This report tests whether the Word-spec components add value over ablations: removing interaction edges, removing signed gene directions, or removing the biological objective from inner model selection.

## Full Model Snapshot

- clinical_benefit / melanoma_anti_pd1_primary: n=160, pooled AUROC=0.455, mean fold AUROC=0.428, pooled ECE=0.231.
- clinical_benefit / melanoma_binary_response_stress: n=29, pooled AUROC=0.596, mean fold AUROC=0.556, pooled ECE=0.219.
- clinical_benefit / melanoma_core_high_evidence: n=117, pooled AUROC=0.517, mean fold AUROC=0.478, pooled ECE=0.239.
- clinical_benefit / melanoma_core_plus_phs000452: n=238, pooled AUROC=0.509, mean fold AUROC=0.469, pooled ECE=0.114.
- clinical_benefit / melanoma_recist_supported_primary: n=131, pooled AUROC=0.515, mean fold AUROC=0.488, pooled ECE=0.199.
- clinical_benefit / pan_cancer_response_all: n=554, pooled AUROC=0.568, mean fold AUROC=0.621, pooled ECE=0.073.
- clinical_benefit / pan_cancer_without_secondary: n=369, pooled AUROC=0.529, mean fold AUROC=0.548, pooled ECE=0.083.
- clinical_benefit / secondary_confounded_transfer: n=185, pooled AUROC=0.580, mean fold AUROC=0.614, pooled ECE=0.035.
- primary_recist / melanoma_anti_pd1_primary: n=160, pooled AUROC=0.506, mean fold AUROC=0.456, pooled ECE=0.187.
- primary_recist / melanoma_binary_response_stress: n=29, pooled AUROC=0.596, mean fold AUROC=0.556, pooled ECE=0.219.
- primary_recist / melanoma_core_high_evidence: n=117, pooled AUROC=0.559, mean fold AUROC=0.546, pooled ECE=0.335.
- primary_recist / melanoma_core_plus_phs000452: n=238, pooled AUROC=0.477, mean fold AUROC=0.475, pooled ECE=0.192.
- primary_recist / melanoma_recist_supported_primary: n=131, pooled AUROC=0.501, mean fold AUROC=0.551, pooled ECE=0.159.
- primary_recist / pan_cancer_response_all: n=554, pooled AUROC=0.613, mean fold AUROC=0.577, pooled ECE=0.156.
- primary_recist / pan_cancer_without_secondary: n=369, pooled AUROC=0.566, mean fold AUROC=0.541, pooled ECE=0.225.
- primary_recist / secondary_confounded_transfer: n=185, pooled AUROC=0.595, mean fold AUROC=0.658, pooled ECE=0.189.
- strict_recist / melanoma_anti_pd1_primary: n=137, pooled AUROC=0.512, mean fold AUROC=0.422, pooled ECE=0.205.
- strict_recist / melanoma_binary_response_stress: n=29, pooled AUROC=0.596, mean fold AUROC=0.556, pooled ECE=0.219.
- strict_recist / melanoma_core_high_evidence: n=95, pooled AUROC=0.575, mean fold AUROC=0.548, pooled ECE=0.249.
- strict_recist / melanoma_core_plus_phs000452: n=200, pooled AUROC=0.495, mean fold AUROC=0.460, pooled ECE=0.139.
- strict_recist / melanoma_recist_supported_primary: n=108, pooled AUROC=0.569, mean fold AUROC=0.566, pooled ECE=0.197.
- strict_recist / pan_cancer_response_all: n=480, pooled AUROC=0.618, mean fold AUROC=0.614, pooled ECE=0.110.
- strict_recist / pan_cancer_without_secondary: n=316, pooled AUROC=0.559, mean fold AUROC=0.545, pooled ECE=0.170.
- strict_recist / secondary_confounded_transfer: n=164, pooled AUROC=0.574, mean fold AUROC=0.624, pooled ECE=0.180.

## Component Ablations

- primary_recist / melanoma_anti_pd1_primary vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.506, ablation AUROC=0.506, delta=0.000, FDR q=1.000 (component_not_supported).
- primary_recist / melanoma_anti_pd1_primary vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.506, ablation AUROC=0.584, delta=-0.078, FDR q=0.040 (component_not_supported).
- primary_recist / melanoma_anti_pd1_primary vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.506, ablation AUROC=0.471, delta=0.035, FDR q=0.383 (point_estimate_component_gain).
- primary_recist / melanoma_binary_response_stress vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.596, ablation AUROC=0.596, delta=0.000, FDR q=1.000 (component_not_supported).
- primary_recist / melanoma_binary_response_stress vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.596, ablation AUROC=0.677, delta=-0.081, FDR q=0.713 (component_not_supported).
- primary_recist / melanoma_binary_response_stress vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.596, ablation AUROC=0.540, delta=0.056, FDR q=0.756 (point_estimate_component_gain).
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.559, ablation AUROC=0.559, delta=0.000, FDR q=1.000 (component_not_supported).
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.559, ablation AUROC=0.603, delta=-0.044, FDR q=0.216 (component_not_supported).
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.559, ablation AUROC=0.553, delta=0.006, FDR q=0.974 (point_estimate_component_gain).
- primary_recist / melanoma_core_plus_phs000452 vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.477, ablation AUROC=0.477, delta=0.000, FDR q=1.000 (component_not_supported).
- primary_recist / melanoma_core_plus_phs000452 vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.477, ablation AUROC=0.543, delta=-0.066, FDR q=0.007 (component_not_supported).
- primary_recist / melanoma_core_plus_phs000452 vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.477, ablation AUROC=0.471, delta=0.007, FDR q=0.858 (point_estimate_component_gain).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.501, ablation AUROC=0.501, delta=0.000, FDR q=1.000 (component_not_supported).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.501, ablation AUROC=0.547, delta=-0.046, FDR q=0.221 (component_not_supported).
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.501, ablation AUROC=0.546, delta=-0.045, FDR q=0.383 (component_not_supported).
- primary_recist / pan_cancer_response_all vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.613, ablation AUROC=0.613, delta=0.000, FDR q=1.000 (component_not_supported).
- primary_recist / pan_cancer_response_all vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.613, ablation AUROC=0.643, delta=-0.030, FDR q=0.012 (component_not_supported).
- primary_recist / pan_cancer_response_all vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.613, ablation AUROC=0.608, delta=0.004, FDR q=0.906 (point_estimate_component_gain).
- primary_recist / pan_cancer_without_secondary vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.566, ablation AUROC=0.566, delta=0.000, FDR q=1.000 (component_not_supported).
- primary_recist / pan_cancer_without_secondary vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.566, ablation AUROC=0.591, delta=-0.025, FDR q=0.219 (component_not_supported).
- primary_recist / pan_cancer_without_secondary vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.566, ablation AUROC=0.590, delta=-0.024, FDR q=0.271 (component_not_supported).
- primary_recist / secondary_confounded_transfer vs EcoNiche-Opt-WordNoBioObjective: full AUROC=0.595, ablation AUROC=0.595, delta=0.000, FDR q=1.000 (component_not_supported).
- primary_recist / secondary_confounded_transfer vs EcoNiche-Opt-WordNoInteraction: full AUROC=0.595, ablation AUROC=0.590, delta=0.006, FDR q=0.885 (point_estimate_component_gain).
- primary_recist / secondary_confounded_transfer vs EcoNiche-Opt-WordUnsignedGraph: full AUROC=0.595, ablation AUROC=0.589, delta=0.006, FDR q=0.885 (point_estimate_component_gain).

## Claim Boundary

Treat component gains as supported only when paired bootstrap deltas are positive and FDR-supported in the pre-specified stratum. Point-estimate gains can motivate mechanistic interpretation but should not be worded as statistically proven superiority.
