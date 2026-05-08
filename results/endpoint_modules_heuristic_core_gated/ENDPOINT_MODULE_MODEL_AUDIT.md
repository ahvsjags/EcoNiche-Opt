# Endpoint-Stratified Module Model Audit

This audit separates endpoint definitions, cancer/therapy strata, the Word-spec signed-rank ecological graph model, module-level priors, strong immune signatures, calibration, and decision-curve outputs.

## Endpoint Definitions

- strict_recist: CR/PR/MR/R/DCB vs PD/NR/NDB; SD is excluded.
- primary_recist: CR/PR/MR/R/DCB vs SD/PD/NR/NDB; this is the conservative primary endpoint.
- clinical_benefit: CR/PR/MR/SD/R/DCB vs PD/NR/NDB.

## Main Result Snapshot

- primary_recist / melanoma_core_high_evidence: EcoNiche-Opt-HeuristicEcology pooled AUROC=0.687, mean fold AUROC=0.642, ECE=0.285; best comparator=EcoNiche-Opt-ModulePriorFixed AUROC=0.705, delta=-0.018.
- primary_recist / melanoma_recist_supported_primary: EcoNiche-Opt-HeuristicEcology pooled AUROC=0.521, mean fold AUROC=0.537, ECE=0.217; best comparator=EcoNiche-Opt-ModulePriorFixed AUROC=0.685, delta=-0.164.

## Strong Signature Claim Gate

- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-ModulePriorFixed: target AUROC=0.687, baseline AUROC=0.705, bootstrap delta=-0.023, 95% CI [-0.150, 0.104], FDR q=0.936; target is not above this comparator.
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordNoInteraction: target AUROC=0.687, baseline AUROC=0.701, bootstrap delta=-0.015, 95% CI [-0.040, 0.011], FDR q=0.550; target is not above this comparator.
- primary_recist / melanoma_core_high_evidence vs CXCL9: target AUROC=0.687, baseline AUROC=0.666, bootstrap delta=0.017, 95% CI [-0.113, 0.152], FDR q=0.936; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs TIG: target AUROC=0.687, baseline AUROC=0.666, bootstrap delta=0.016, 95% CI [-0.114, 0.150], FDR q=0.936; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs IFNG: target AUROC=0.687, baseline AUROC=0.664, bootstrap delta=0.018, 95% CI [-0.117, 0.160], FDR q=0.936; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs APM: target AUROC=0.687, baseline AUROC=0.664, bootstrap delta=0.018, 95% CI [-0.117, 0.157], FDR q=0.936; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs TIDE_dysfunction: target AUROC=0.687, baseline AUROC=0.655, bootstrap delta=0.027, 95% CI [-0.110, 0.163], FDR q=0.936; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-ImmuneComposite: target AUROC=0.687, baseline AUROC=0.632, bootstrap delta=0.051, 95% CI [-0.081, 0.194], FDR q=0.936; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs CYT: target AUROC=0.687, baseline AUROC=0.630, bootstrap delta=0.052, 95% CI [-0.083, 0.186], FDR q=0.936; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs PDCD1LG2: target AUROC=0.687, baseline AUROC=0.591, bootstrap delta=0.092, 95% CI [-0.047, 0.241], FDR q=0.528; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs IPRES: target AUROC=0.687, baseline AUROC=0.567, bootstrap delta=0.117, 95% CI [-0.018, 0.253], FDR q=0.317; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs TIDE_exclusion: target AUROC=0.687, baseline AUROC=0.544, bootstrap delta=0.138, 95% CI [-0.003, 0.287], FDR q=0.252; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordUnsignedGraph: target AUROC=0.687, baseline AUROC=0.471, bootstrap delta=0.213, 95% CI [0.101, 0.328], FDR q=0.000; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordNoBioObjective: target AUROC=0.687, baseline AUROC=0.355, bootstrap delta=0.327, 95% CI [0.175, 0.472], FDR q=0.018; target is above this comparator.
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-ModulePriorFixed: target AUROC=0.521, baseline AUROC=0.685, bootstrap delta=-0.161, 95% CI [-0.249, -0.070], FDR q=0.036; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs CXCL9: target AUROC=0.521, baseline AUROC=0.649, bootstrap delta=-0.126, 95% CI [-0.230, -0.023], FDR q=0.058; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs TIDE_dysfunction: target AUROC=0.521, baseline AUROC=0.646, bootstrap delta=-0.122, 95% CI [-0.224, -0.015], FDR q=0.062; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs TIG: target AUROC=0.521, baseline AUROC=0.644, bootstrap delta=-0.121, 95% CI [-0.214, -0.020], FDR q=0.058; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs APM: target AUROC=0.521, baseline AUROC=0.642, bootstrap delta=-0.118, 95% CI [-0.242, 0.001], FDR q=0.094; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs IFNG: target AUROC=0.521, baseline AUROC=0.640, bootstrap delta=-0.116, 95% CI [-0.216, -0.015], FDR q=0.060; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs CYT: target AUROC=0.521, baseline AUROC=0.624, bootstrap delta=-0.101, 95% CI [-0.194, -0.001], FDR q=0.094; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-ImmuneComposite: target AUROC=0.521, baseline AUROC=0.620, bootstrap delta=-0.097, 95% CI [-0.207, 0.011], FDR q=0.144; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs PDCD1LG2: target AUROC=0.521, baseline AUROC=0.589, bootstrap delta=-0.068, 95% CI [-0.179, 0.051], FDR q=0.399; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordNoBioObjective: target AUROC=0.521, baseline AUROC=0.575, bootstrap delta=-0.053, 95% CI [-0.172, 0.062], FDR q=0.414; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs IPRES: target AUROC=0.521, baseline AUROC=0.573, bootstrap delta=-0.053, 95% CI [-0.187, 0.066], FDR q=0.487; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs TIDE_exclusion: target AUROC=0.521, baseline AUROC=0.539, bootstrap delta=-0.018, 95% CI [-0.149, 0.109], FDR q=0.784; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordNoInteraction: target AUROC=0.521, baseline AUROC=0.506, bootstrap delta=0.017, 95% CI [-0.073, 0.111], FDR q=0.784; target is above this comparator.
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordUnsignedGraph: target AUROC=0.521, baseline AUROC=0.453, bootstrap delta=0.069, 95% CI [-0.043, 0.172], FDR q=0.288; target is above this comparator.

## Label Sensitivity Audit

- clinical_benefit: used=554, dropped=0, responders=274, nonresponders=280.
- primary_recist: used=554, dropped=0, responders=200, nonresponders=354.
- strict_recist: used=480, dropped=74, responders=200, nonresponders=280.

## Interpretation Guardrail

Do not claim superiority over all existing models unless the paired strong-signature comparisons are positive and FDR-supported in the pre-specified primary stratum. The Word-spec graph terms should be claimed as component gains only where the ablation table supports them.
