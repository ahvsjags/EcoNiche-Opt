# Endpoint-Stratified Module Model Audit

This audit separates endpoint definitions, cancer/therapy strata, the Word-spec signed-rank ecological graph model, module-level priors, strong immune signatures, calibration, and decision-curve outputs.

## Endpoint Definitions

- strict_recist: CR/PR/MR/R/DCB vs PD/NR/NDB; SD is excluded.
- primary_recist: CR/PR/MR/R/DCB vs SD/PD/NR/NDB; this is the conservative primary endpoint.
- clinical_benefit: CR/PR/MR/SD/R/DCB vs PD/NR/NDB.

## Main Result Snapshot

- primary_recist / melanoma_core_high_evidence: EcoNiche-Opt-HeuristicEcology pooled AUROC=0.600, mean fold AUROC=0.579, ECE=0.341; best comparator=EcoNiche-Opt-ModulePriorFixed AUROC=0.705, delta=-0.105.
- primary_recist / melanoma_recist_supported_primary: EcoNiche-Opt-HeuristicEcology pooled AUROC=0.621, mean fold AUROC=0.595, ECE=0.237; best comparator=EcoNiche-Opt-ModulePriorFixed AUROC=0.685, delta=-0.063.

## Strong Signature Claim Gate

- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-ModulePriorFixed: target AUROC=0.600, baseline AUROC=0.705, bootstrap delta=-0.108, 95% CI [-0.227, 0.005], FDR q=0.369; target is not above this comparator.
- primary_recist / melanoma_core_high_evidence vs CXCL9: target AUROC=0.600, baseline AUROC=0.666, bootstrap delta=-0.068, 95% CI [-0.194, 0.052], FDR q=0.468; target is not above this comparator.
- primary_recist / melanoma_core_high_evidence vs TIG: target AUROC=0.600, baseline AUROC=0.666, bootstrap delta=-0.068, 95% CI [-0.189, 0.052], FDR q=0.468; target is not above this comparator.
- primary_recist / melanoma_core_high_evidence vs IFNG: target AUROC=0.600, baseline AUROC=0.664, bootstrap delta=-0.066, 95% CI [-0.191, 0.056], FDR q=0.468; target is not above this comparator.
- primary_recist / melanoma_core_high_evidence vs APM: target AUROC=0.600, baseline AUROC=0.664, bootstrap delta=-0.066, 95% CI [-0.195, 0.061], FDR q=0.468; target is not above this comparator.
- primary_recist / melanoma_core_high_evidence vs TIDE_dysfunction: target AUROC=0.600, baseline AUROC=0.655, bootstrap delta=-0.057, 95% CI [-0.186, 0.065], FDR q=0.532; target is not above this comparator.
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-ImmuneComposite: target AUROC=0.600, baseline AUROC=0.632, bootstrap delta=-0.033, 95% CI [-0.163, 0.091], FDR q=0.682; target is not above this comparator.
- primary_recist / melanoma_core_high_evidence vs CYT: target AUROC=0.600, baseline AUROC=0.630, bootstrap delta=-0.032, 95% CI [-0.148, 0.094], FDR q=0.682; target is not above this comparator.
- primary_recist / melanoma_core_high_evidence vs PDCD1LG2: target AUROC=0.600, baseline AUROC=0.591, bootstrap delta=0.008, 95% CI [-0.130, 0.144], FDR q=0.900; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs IPRES: target AUROC=0.600, baseline AUROC=0.567, bootstrap delta=0.032, 95% CI [-0.121, 0.184], FDR q=0.722; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordNoInteraction: target AUROC=0.600, baseline AUROC=0.556, bootstrap delta=0.046, 95% CI [-0.014, 0.110], FDR q=0.468; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs TIDE_exclusion: target AUROC=0.600, baseline AUROC=0.544, bootstrap delta=0.053, 95% CI [-0.094, 0.206], FDR q=0.620; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordNoBioObjective: target AUROC=0.600, baseline AUROC=0.475, bootstrap delta=0.123, 95% CI [-0.012, 0.278], FDR q=0.369; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordUnsignedGraph: target AUROC=0.600, baseline AUROC=0.355, bootstrap delta=0.243, 95% CI [0.076, 0.412], FDR q=0.036; target is above this comparator.
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-ModulePriorFixed: target AUROC=0.621, baseline AUROC=0.685, bootstrap delta=-0.062, 95% CI [-0.139, 0.015], FDR q=0.720; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs CXCL9: target AUROC=0.621, baseline AUROC=0.649, bootstrap delta=-0.027, 95% CI [-0.126, 0.068], FDR q=0.897; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs TIDE_dysfunction: target AUROC=0.621, baseline AUROC=0.646, bootstrap delta=-0.023, 95% CI [-0.128, 0.076], FDR q=0.897; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs TIG: target AUROC=0.621, baseline AUROC=0.644, bootstrap delta=-0.022, 95% CI [-0.120, 0.073], FDR q=0.897; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs APM: target AUROC=0.621, baseline AUROC=0.642, bootstrap delta=-0.019, 95% CI [-0.129, 0.089], FDR q=0.914; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs IFNG: target AUROC=0.621, baseline AUROC=0.640, bootstrap delta=-0.017, 95% CI [-0.114, 0.077], FDR q=0.914; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs CYT: target AUROC=0.621, baseline AUROC=0.624, bootstrap delta=-0.002, 95% CI [-0.105, 0.099], FDR q=0.956; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-ImmuneComposite: target AUROC=0.621, baseline AUROC=0.620, bootstrap delta=0.002, 95% CI [-0.108, 0.110], FDR q=0.956; target is above this comparator.
- primary_recist / melanoma_recist_supported_primary vs PDCD1LG2: target AUROC=0.621, baseline AUROC=0.589, bootstrap delta=0.031, 95% CI [-0.094, 0.157], FDR q=0.897; target is above this comparator.
- primary_recist / melanoma_recist_supported_primary vs IPRES: target AUROC=0.621, baseline AUROC=0.573, bootstrap delta=0.046, 95% CI [-0.077, 0.160], FDR q=0.897; target is above this comparator.
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordNoInteraction: target AUROC=0.621, baseline AUROC=0.557, bootstrap delta=0.063, 95% CI [-0.031, 0.162], FDR q=0.749; target is above this comparator.
- primary_recist / melanoma_recist_supported_primary vs TIDE_exclusion: target AUROC=0.621, baseline AUROC=0.539, bootstrap delta=0.081, 95% CI [-0.041, 0.199], FDR q=0.749; target is above this comparator.
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordUnsignedGraph: target AUROC=0.621, baseline AUROC=0.445, bootstrap delta=0.175, 95% CI [0.062, 0.288], FDR q=0.018; target is above this comparator.
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordNoBioObjective: target AUROC=0.621, baseline AUROC=0.438, bootstrap delta=0.182, 95% CI [0.080, 0.277], FDR q=0.000; target is above this comparator.

## Label Sensitivity Audit

- clinical_benefit: used=554, dropped=0, responders=274, nonresponders=280.
- primary_recist: used=554, dropped=0, responders=200, nonresponders=354.
- strict_recist: used=480, dropped=74, responders=200, nonresponders=280.

## Interpretation Guardrail

Do not claim superiority over all existing models unless the paired strong-signature comparisons are positive and FDR-supported in the pre-specified primary stratum. The Word-spec graph terms should be claimed as component gains only where the ablation table supports them.
