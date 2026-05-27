# GPU biological-prior rescue-combo audit

Device: cuda (NVIDIA GeForce RTX 4070 Laptop GPU).
Selected candidate: `0.80*base+0.20*rz__CD247` under `immune` prior and `robust_only` transform policy.
Selection used primary melanoma LODO only; strict external labels were used only for locked scoring.

Primary AUROC=0.778, AUPRC=0.661, balanced accuracy=0.693.
Strict external AUROC=0.690, AUPRC=0.705, balanced accuracy=0.588, ECE=0.130.
Family gate: mean AUROC=0.567, delta=0.124, q=0.011, claim=family_point_estimate_only.