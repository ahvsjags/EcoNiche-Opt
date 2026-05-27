# Aligned Interaction Edge Audit

This registered analysis tests whether response-resistance module interaction edges improve the current module-prior model under training-only selection. It does not use strict external labels for feature selection, thresholding, calibration or regularization selection.

## Primary LODO

- melanoma_core_high_evidence: edge AUROC=0.647 versus baseline AUROC=0.705; delta AUROC=-0.058; delta ECE=-0.110; delta fold-AUROC SD=-0.015; claim=interaction_edge_calibration_and_stability_tradeoff.
- melanoma_recist_supported_primary: edge AUROC=0.616 versus baseline AUROC=0.685; delta AUROC=-0.069; delta ECE=-0.132; delta fold-AUROC SD=-0.015; claim=interaction_edge_calibration_and_stability_tradeoff.

## Strict external

- strict_pd1_like_external: edge AUROC=0.542 versus baseline AUROC=0.572; delta AUROC=-0.030; delta ECE=-0.191; claim=interaction_edge_calibration_tradeoff.

## Interpretation

Interaction edges are not currently supported as a discrimination-improving component of the locked melanoma predictor. They may be retained as an ecological interpretation and calibration/stability diagnostic only where the audit shows lower ECE or lower fold variability.