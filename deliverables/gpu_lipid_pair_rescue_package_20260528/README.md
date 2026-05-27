# EcoNiche-Opt GPU lipid/PI3K pair rescue package

Release tag: `v0.3.4-gpu-lipid-pair-rescue-20260528`.

This package freezes the primary-LODO-selected component-dominant lipid/PI3K rescue score.
The scoring rule uses MAP4K1/TBX3/AXL as the base tumor-immune balance axis and a locked PLA2G2D/PIK3CD pair component.
External labels are not used for gene, weight, threshold, or calibration selection.

Selected candidate: `0.35*base+0.65*(0.35*pct__PLA2G2D+0.65*rz__PIK3CD)`.
Locked threshold: `0.597088384628`.