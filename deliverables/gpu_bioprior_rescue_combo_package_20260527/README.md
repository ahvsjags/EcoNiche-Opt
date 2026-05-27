# GPU Biological-Prior Rescue Combo

This package freezes the GPU-audited lipid/PI3K biological-prior rescue combo for locked melanoma external scoring.

## Locked Candidate

- Candidate: `0.80*base+0.20*rz__PLA2G2D`.
- Genes: MAP4K1, TBX3, AXL, PLA2G2D.
- Threshold: 0.370894.
- Selection boundary: candidate_and_weight_selected_by_primary_lodo_only_with_biological_prior.

## Current Evidence

- Primary melanoma LODO AUROC: 0.777.
- Strict melanoma external AUROC: 0.713.
- Strict melanoma external AUPRC: 0.726.
- Family-comparison FDR q: 0.000.
- GPU device: NVIDIA GeForce RTX 4070 Laptop GPU.

External cohorts must be scored without changing genes, weights, transform policy, threshold, or calibration.
