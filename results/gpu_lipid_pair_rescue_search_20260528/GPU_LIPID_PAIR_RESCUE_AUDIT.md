# GPU lipid/PI3K pair rescue audit

Device: cuda (NVIDIA GeForce RTX 4070 Laptop GPU).
Selected candidate: `0.80*base+0.20*(0.65*pct__PLA2G2D+0.35*rz__PIK3CD)`.
Candidate and weight selection used primary melanoma LODO only; strict external and cBioPortal labels were used only for locked scoring.

Primary AUROC=0.788, AUPRC=0.678, balanced accuracy=0.742.
Strict current external AUROC=0.699, AUPRC=0.718, balanced accuracy=0.620, ECE=0.077.
cBioPortal Liu/DFCI AUROC=0.672, AUPRC=0.678, ECE=0.054, family q=0.051.