# Tumor-Immune Balance Audit

This audit adds a literature-prior MAP4K1-TBX3 tumor-immune balance axis and related ecological balance variants.
The strict external cohorts are never used for training, feature selection, thresholding, calibration, or model selection.

## Primary LODO

- `EcoNiche-Opt-MAP4K1APM-TBX3IPRESBalance`: AUROC=0.771, AUPRC=0.653, BA=0.689
- `EcoNiche-Opt-MAP4K1APM-TBX3IPRESBlend50`: AUROC=0.745, AUPRC=0.602, BA=0.681
- `EcoNiche-Opt-TumorImmuneBalancePairBlend50`: AUROC=0.735, AUPRC=0.630, BA=0.635
- `EcoNiche-Opt-CYTAPM-IPRESBlend50`: AUROC=0.724, AUPRC=0.575, BA=0.665
- `EcoNiche-Opt-TumorImmuneBalancePair`: AUROC=0.720, AUPRC=0.619, BA=0.684

## Strict External

- `EcoNiche-Opt-TumorImmuneBalancePair`: AUROC=0.653, AUPRC=0.661, BA=0.595
- `EcoNiche-Opt-TumorImmuneBalancePairBlend50`: AUROC=0.628, AUPRC=0.646, BA=0.555
- `EcoNiche-Opt-MAP4K1APM-TBX3IPRESBalance`: AUROC=0.609, AUPRC=0.648, BA=0.577
- `EcoNiche-Opt-MAP4K1APM-TBX3IPRESBlend50`: AUROC=0.579, AUPRC=0.617, BA=0.565
- `EcoNiche-Opt-CYTAPM-IPRESBlend50`: AUROC=0.551, AUPRC=0.560, BA=0.544