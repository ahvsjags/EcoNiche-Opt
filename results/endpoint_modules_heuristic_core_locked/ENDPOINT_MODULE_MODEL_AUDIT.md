# Endpoint-Stratified Module Model Audit

This audit separates endpoint definitions, cancer/therapy strata, the Word-spec signed-rank ecological graph model, module-level priors, strong immune signatures, calibration, and decision-curve outputs.

## Endpoint Definitions

- strict_recist: CR/PR/MR/R/DCB vs PD/NR/NDB; SD is excluded.
- primary_recist: CR/PR/MR/R/DCB vs SD/PD/NR/NDB; this is the conservative primary endpoint.
- clinical_benefit: CR/PR/MR/SD/R/DCB vs PD/NR/NDB.

## Main Result Snapshot

- primary_recist / melanoma_core_high_evidence: EcoNiche-Opt-HeuristicEcology pooled AUROC=0.705, mean fold AUROC=0.712, ECE=0.235; best comparator=EcoNiche-Opt-ModulePriorFixed AUROC=0.705, delta=0.000.
- primary_recist / melanoma_recist_supported_primary: EcoNiche-Opt-HeuristicEcology pooled AUROC=0.685, mean fold AUROC=0.664, ECE=0.239; best comparator=EcoNiche-Opt-ModulePriorFixed AUROC=0.685, delta=0.000.

## Strong Signature Claim Gate

- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-ModulePriorFixed: target AUROC=0.705, baseline AUROC=0.705, bootstrap delta=0.000, 95% CI [0.000, 0.000], FDR q=1.000; target is not above this comparator.
- primary_recist / melanoma_core_high_evidence vs CXCL9: target AUROC=0.705, baseline AUROC=0.666, bootstrap delta=0.040, 95% CI [-0.011, 0.093], FDR q=0.153; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs TIG: target AUROC=0.705, baseline AUROC=0.666, bootstrap delta=0.040, 95% CI [-0.001, 0.083], FDR q=0.111; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs IFNG: target AUROC=0.705, baseline AUROC=0.664, bootstrap delta=0.042, 95% CI [-0.006, 0.091], FDR q=0.121; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs APM: target AUROC=0.705, baseline AUROC=0.664, bootstrap delta=0.041, 95% CI [-0.029, 0.114], FDR q=0.260; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs TIDE_dysfunction: target AUROC=0.705, baseline AUROC=0.655, bootstrap delta=0.051, 95% CI [-0.003, 0.105], FDR q=0.111; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-ImmuneComposite: target AUROC=0.705, baseline AUROC=0.632, bootstrap delta=0.075, 95% CI [0.015, 0.135], FDR q=0.041; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs CYT: target AUROC=0.705, baseline AUROC=0.630, bootstrap delta=0.076, 95% CI [0.023, 0.133], FDR q=0.036; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordNoBioObjective: target AUROC=0.705, baseline AUROC=0.609, bootstrap delta=0.102, 95% CI [-0.027, 0.236], FDR q=0.137; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs PDCD1LG2: target AUROC=0.705, baseline AUROC=0.591, bootstrap delta=0.116, 95% CI [0.030, 0.199], FDR q=0.036; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs IPRES: target AUROC=0.705, baseline AUROC=0.567, bootstrap delta=0.140, 95% CI [-0.007, 0.292], FDR q=0.111; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordUnsignedGraph: target AUROC=0.705, baseline AUROC=0.549, bootstrap delta=0.161, 95% CI [0.035, 0.282], FDR q=0.041; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordNoInteraction: target AUROC=0.705, baseline AUROC=0.545, bootstrap delta=0.162, 95% CI [0.044, 0.283], FDR q=0.036; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs TIDE_exclusion: target AUROC=0.705, baseline AUROC=0.544, bootstrap delta=0.161, 95% CI [0.007, 0.310], FDR q=0.095; target is above this comparator.
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-ModulePriorFixed: target AUROC=0.685, baseline AUROC=0.685, bootstrap delta=0.000, 95% CI [0.000, 0.000], FDR q=1.000; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs CXCL9: target AUROC=0.685, baseline AUROC=0.649, bootstrap delta=0.035, 95% CI [-0.015, 0.087], FDR q=0.209; target is above this comparator.
- primary_recist / melanoma_recist_supported_primary vs TIDE_dysfunction: target AUROC=0.685, baseline AUROC=0.646, bootstrap delta=0.039, 95% CI [-0.009, 0.088], FDR q=0.139; target is above this comparator.
- primary_recist / melanoma_recist_supported_primary vs TIG: target AUROC=0.685, baseline AUROC=0.644, bootstrap delta=0.040, 95% CI [0.000, 0.086], FDR q=0.090; target is above this comparator.
- primary_recist / melanoma_recist_supported_primary vs APM: target AUROC=0.685, baseline AUROC=0.642, bootstrap delta=0.043, 95% CI [-0.025, 0.113], FDR q=0.229; target is above this comparator.
- primary_recist / melanoma_recist_supported_primary vs IFNG: target AUROC=0.685, baseline AUROC=0.640, bootstrap delta=0.045, 95% CI [-0.002, 0.094], FDR q=0.101; target is above this comparator.
- primary_recist / melanoma_recist_supported_primary vs CYT: target AUROC=0.685, baseline AUROC=0.624, bootstrap delta=0.060, 95% CI [0.009, 0.113], FDR q=0.036; target is above this comparator.
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-ImmuneComposite: target AUROC=0.685, baseline AUROC=0.620, bootstrap delta=0.064, 95% CI [0.004, 0.125], FDR q=0.082; target is above this comparator.
- primary_recist / melanoma_recist_supported_primary vs PDCD1LG2: target AUROC=0.685, baseline AUROC=0.589, bootstrap delta=0.093, 95% CI [0.013, 0.177], FDR q=0.036; target is above this comparator.
- primary_recist / melanoma_recist_supported_primary vs IPRES: target AUROC=0.685, baseline AUROC=0.573, bootstrap delta=0.108, 95% CI [-0.024, 0.239], FDR q=0.123; target is above this comparator.
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordNoInteraction: target AUROC=0.685, baseline AUROC=0.550, bootstrap delta=0.133, 95% CI [0.039, 0.229], FDR q=0.012; target is above this comparator.
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordNoBioObjective: target AUROC=0.685, baseline AUROC=0.548, bootstrap delta=0.133, 95% CI [0.004, 0.266], FDR q=0.088; target is above this comparator.
- primary_recist / melanoma_recist_supported_primary vs TIDE_exclusion: target AUROC=0.685, baseline AUROC=0.539, bootstrap delta=0.143, 95% CI [0.004, 0.277], FDR q=0.088; target is above this comparator.
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordUnsignedGraph: target AUROC=0.685, baseline AUROC=0.536, bootstrap delta=0.147, 95% CI [0.045, 0.254], FDR q=0.000; target is above this comparator.

## Label Sensitivity Audit

- clinical_benefit: used=554, dropped=0, responders=274, nonresponders=280.
- primary_recist: used=554, dropped=0, responders=200, nonresponders=354.
- strict_recist: used=480, dropped=74, responders=200, nonresponders=280.

## Interpretation Guardrail

Do not claim superiority over all existing models unless the paired strong-signature comparisons are positive and FDR-supported in the pre-specified primary stratum. The Word-spec graph terms should be claimed as component gains only where the ablation table supports them.
