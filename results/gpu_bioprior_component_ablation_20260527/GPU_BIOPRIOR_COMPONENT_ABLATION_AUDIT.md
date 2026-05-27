# GPU biological-prior component ablation

Device: cuda (NVIDIA GeForce RTX 4070 Laptop GPU).
Target candidate: `0.80*base+0.20*rz__PLA2G2D`.
Ablated candidate: `base_rescue_robust`.

This audit evaluates the frozen PLA2G2D lipid/PI3K rescue component against the base MAP4K1-TBX3/AXL rescue axis.
It is an evaluation-only component ablation; strict external labels were not used for candidate or weight selection.

Strict external AUROC changes from 0.686 to 0.713; AUPRC changes from 0.689 to 0.726.