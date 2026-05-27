# Ecological Polarity Candidate Audit

This audit tests a predeclared family of immune-effector, antigen-presentation, dedifferentiation, IPRES, and stromal polarity scores.
Candidate choice is made from primary melanoma LODO evidence only. Current strict external stress rows are diagnostic and cannot define the locked model.

- `high_evidence_primary` `primary_selected_candidate` `MAP4K1_minus_TBX3_AXL` weight=1.25: primary AUROC=0.784; strict external AUROC=0.656; boundary=selected_by_primary_lodo_only_not_by_external
- `high_evidence_primary` `current_external_stress_best` `MAP4K1_minus_TBX3_AXL` weight=1.00: primary AUROC=0.779; strict external AUROC=0.679; boundary=current_external_stress_screen_not_a_locked_selection_claim
- `expanded_primary_with_mgh` `primary_selected_candidate` `MAP4K1_minus_TBX3_AXL` weight=1.25: primary AUROC=0.729; strict external AUROC=0.656; boundary=selected_by_primary_lodo_only_not_by_external
- `expanded_primary_with_mgh` `current_external_stress_best` `MAP4K1_minus_TBX3_AXL` weight=1.00: primary AUROC=0.708; strict external AUROC=0.679; boundary=current_external_stress_screen_not_a_locked_selection_claim

## Claim Boundary

A candidate can support the strict external target only if the primary-selected row reaches AUROC >=0.70 on strict external scoring. Otherwise it remains a negative optimization audit.
