# GPU biological-prior rescue-combo audit

Device: cuda (NVIDIA GeForce RTX 4070 Laptop GPU).
Selected candidate: `0.80*base+0.20*pct__PLA2G2D` under `immune` prior and `all` transform policy.
Selection used primary melanoma LODO only; strict external labels were used only for locked scoring.

Primary AUROC=0.783, AUPRC=0.673, balanced accuracy=0.683.
Strict external AUROC=0.699, AUPRC=0.714, balanced accuracy=0.627, ECE=0.092.
Family gate: mean AUROC=0.567, delta=0.132, q=0.004, claim=family_point_estimate_only.