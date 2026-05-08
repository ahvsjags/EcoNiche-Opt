# Endpoint-Stratified Module Model Audit

This audit separates endpoint definitions, cancer/therapy strata, the Word-spec signed-rank ecological graph model, module-level priors, strong immune signatures, calibration, and decision-curve outputs.

## Endpoint Definitions

- strict_recist: CR/PR/MR/R/DCB vs PD/NR/NDB; SD is excluded.
- primary_recist: CR/PR/MR/R/DCB vs SD/PD/NR/NDB; this is the conservative primary endpoint.
- clinical_benefit: CR/PR/MR/SD/R/DCB vs PD/NR/NDB.

## Main Result Snapshot

- primary_recist / melanoma_core_high_evidence: EcoNiche-Opt-HeuristicEcology pooled AUROC=0.585, mean fold AUROC=0.466, ECE=0.277; best comparator=EcoNiche-Opt-ModulePriorFixed AUROC=0.705, delta=-0.120.
- primary_recist / melanoma_recist_supported_primary: EcoNiche-Opt-HeuristicEcology pooled AUROC=0.373, mean fold AUROC=0.345, ECE=0.314; best comparator=EcoNiche-Opt-ModulePriorFixed AUROC=0.685, delta=-0.312.

## Strong Signature Claim Gate

- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-ModulePriorFixed: target AUROC=0.585, baseline AUROC=0.705, bootstrap delta=-0.126, 95% CI [-0.257, -0.004], FDR q=0.399; target is not above this comparator.
- primary_recist / melanoma_core_high_evidence vs CXCL9: target AUROC=0.585, baseline AUROC=0.666, bootstrap delta=-0.086, 95% CI [-0.222, 0.048], FDR q=0.399; target is not above this comparator.
- primary_recist / melanoma_core_high_evidence vs TIG: target AUROC=0.585, baseline AUROC=0.666, bootstrap delta=-0.086, 95% CI [-0.222, 0.044], FDR q=0.399; target is not above this comparator.
- primary_recist / melanoma_core_high_evidence vs IFNG: target AUROC=0.585, baseline AUROC=0.664, bootstrap delta=-0.084, 95% CI [-0.226, 0.049], FDR q=0.399; target is not above this comparator.
- primary_recist / melanoma_core_high_evidence vs APM: target AUROC=0.585, baseline AUROC=0.664, bootstrap delta=-0.084, 95% CI [-0.223, 0.052], FDR q=0.399; target is not above this comparator.
- primary_recist / melanoma_core_high_evidence vs TIDE_dysfunction: target AUROC=0.585, baseline AUROC=0.655, bootstrap delta=-0.075, 95% CI [-0.211, 0.052], FDR q=0.399; target is not above this comparator.
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-ImmuneComposite: target AUROC=0.585, baseline AUROC=0.632, bootstrap delta=-0.051, 95% CI [-0.192, 0.080], FDR q=0.581; target is not above this comparator.
- primary_recist / melanoma_core_high_evidence vs CYT: target AUROC=0.585, baseline AUROC=0.630, bootstrap delta=-0.050, 95% CI [-0.180, 0.072], FDR q=0.581; target is not above this comparator.
- primary_recist / melanoma_core_high_evidence vs PDCD1LG2: target AUROC=0.585, baseline AUROC=0.591, bootstrap delta=-0.010, 95% CI [-0.145, 0.125], FDR q=0.866; target is not above this comparator.
- primary_recist / melanoma_core_high_evidence vs IPRES: target AUROC=0.585, baseline AUROC=0.567, bootstrap delta=0.014, 95% CI [-0.128, 0.158], FDR q=0.866; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordNoInteraction: target AUROC=0.585, baseline AUROC=0.555, bootstrap delta=0.029, 95% CI [-0.033, 0.096], FDR q=0.540; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs TIDE_exclusion: target AUROC=0.585, baseline AUROC=0.544, bootstrap delta=0.035, 95% CI [-0.112, 0.187], FDR q=0.738; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordNoBioObjective: target AUROC=0.585, baseline AUROC=0.521, bootstrap delta=0.062, 95% CI [-0.031, 0.153], FDR q=0.399; target is above this comparator.
- primary_recist / melanoma_core_high_evidence vs EcoNiche-Opt-WordUnsignedGraph: target AUROC=0.585, baseline AUROC=0.431, bootstrap delta=0.147, 95% CI [-0.013, 0.309], FDR q=0.399; target is above this comparator.
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-ModulePriorFixed: target AUROC=0.373, baseline AUROC=0.685, bootstrap delta=-0.310, 95% CI [-0.434, -0.186], FDR q=0.000; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs CXCL9: target AUROC=0.373, baseline AUROC=0.649, bootstrap delta=-0.275, 95% CI [-0.406, -0.144], FDR q=0.000; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs TIDE_dysfunction: target AUROC=0.373, baseline AUROC=0.646, bootstrap delta=-0.271, 95% CI [-0.408, -0.145], FDR q=0.000; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs TIG: target AUROC=0.373, baseline AUROC=0.644, bootstrap delta=-0.270, 95% CI [-0.398, -0.140], FDR q=0.000; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs APM: target AUROC=0.373, baseline AUROC=0.642, bootstrap delta=-0.267, 95% CI [-0.416, -0.126], FDR q=0.000; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs IFNG: target AUROC=0.373, baseline AUROC=0.640, bootstrap delta=-0.265, 95% CI [-0.392, -0.136], FDR q=0.000; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs CYT: target AUROC=0.373, baseline AUROC=0.624, bootstrap delta=-0.250, 95% CI [-0.381, -0.112], FDR q=0.000; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-ImmuneComposite: target AUROC=0.373, baseline AUROC=0.620, bootstrap delta=-0.246, 95% CI [-0.378, -0.112], FDR q=0.000; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs PDCD1LG2: target AUROC=0.373, baseline AUROC=0.589, bootstrap delta=-0.217, 95% CI [-0.354, -0.081], FDR q=0.005; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs IPRES: target AUROC=0.373, baseline AUROC=0.573, bootstrap delta=-0.202, 95% CI [-0.315, -0.085], FDR q=0.003; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs TIDE_exclusion: target AUROC=0.373, baseline AUROC=0.539, bootstrap delta=-0.167, 95% CI [-0.286, -0.045], FDR q=0.005; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordNoBioObjective: target AUROC=0.373, baseline AUROC=0.504, bootstrap delta=-0.129, 95% CI [-0.246, -0.016], FDR q=0.023; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordUnsignedGraph: target AUROC=0.373, baseline AUROC=0.503, bootstrap delta=-0.130, 95% CI [-0.265, 0.008], FDR q=0.066; target is not above this comparator.
- primary_recist / melanoma_recist_supported_primary vs EcoNiche-Opt-WordNoInteraction: target AUROC=0.373, baseline AUROC=0.467, bootstrap delta=-0.094, 95% CI [-0.164, -0.024], FDR q=0.009; target is not above this comparator.

## Label Sensitivity Audit

- clinical_benefit: used=554, dropped=0, responders=274, nonresponders=280.
- primary_recist: used=554, dropped=0, responders=200, nonresponders=354.
- strict_recist: used=480, dropped=74, responders=200, nonresponders=280.

## Interpretation Guardrail

Do not claim superiority over all existing models unless the paired strong-signature comparisons are positive and FDR-supported in the pre-specified primary stratum. The Word-spec graph terms should be claimed as component gains only where the ablation table supports them.
