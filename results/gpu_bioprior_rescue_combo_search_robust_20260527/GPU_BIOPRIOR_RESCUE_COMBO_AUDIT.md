# GPU biological-prior rescue-combo audit

Device: cuda (NVIDIA GeForce RTX 4070 Laptop GPU).
Selected candidate: `0.80*base+0.20*rz__PLA2G2D` under `lipid_pi3k` prior and `robust_only` transform policy.
Selection used primary melanoma LODO only; strict external labels were used only for locked scoring.

Primary AUROC=0.777, AUPRC=0.665, balanced accuracy=0.659.
Strict external AUROC=0.713, AUPRC=0.726, balanced accuracy=0.629, ECE=0.128.
Family gate: mean AUROC=0.567, delta=0.146, q=0.000, claim=strict_external_family_FDR_supported_numeric_target_met.