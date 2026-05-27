# GPU lipid/PI3K pair rescue audit

Device: cuda (NVIDIA GeForce RTX 4070 Laptop GPU).
Selected candidate: `0.35*base+0.65*(0.50*pct__PLA2G2D+0.50*rz__PIK3CD)`.
Candidate and weight selection used primary melanoma LODO only; strict external and cBioPortal labels were used only for locked scoring.

Primary AUROC=0.777, AUPRC=0.697, balanced accuracy=0.628.
Strict current external AUROC=0.709, AUPRC=0.737, balanced accuracy=0.633, ECE=0.078.
cBioPortal Liu/DFCI AUROC=0.704, AUPRC=0.710, ECE=0.105, family q=0.001.