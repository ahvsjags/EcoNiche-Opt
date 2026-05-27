# GPU lipid/PI3K pair rescue audit

Device: cuda (NVIDIA GeForce RTX 4070 Laptop GPU).
Selected candidate: `0.35*base+0.65*(0.35*pct__PLA2G2D+0.65*rz__PIK3CD)`.
Candidate and weight selection used primary melanoma LODO only; strict external and cBioPortal labels were used only for locked scoring.

Primary AUROC=0.772, AUPRC=0.699, balanced accuracy=0.686.
Strict current external AUROC=0.720, AUPRC=0.746, balanced accuracy=0.645, ECE=0.080.
cBioPortal Liu/DFCI AUROC=0.701, AUPRC=0.705, ECE=0.077, family q=0.000.