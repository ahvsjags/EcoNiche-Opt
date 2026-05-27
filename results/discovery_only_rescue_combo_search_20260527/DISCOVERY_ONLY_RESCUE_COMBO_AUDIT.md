# Discovery-only Rescue Combo Search

This audit screens genes and rescue-head combinations using primary melanoma LODO only. Strict external labels are used only after the primary-selected candidate is locked.

- Selected candidate: `0.65*base+0.35*z__SLC6A12`.
- Primary LODO: AUROC=0.845, AUPRC=0.753, balanced accuracy=0.737.
- Strict external GSE145996+PHS000452: AUROC=0.625, AUPRC=0.656, balanced accuracy=0.571, ECE=0.110.
- Eight-signature family comparison: family mean AUROC=0.567, delta=0.058, q=0.228, claim=family_point_estimate_only.

Claim boundary: this is a no-external-label feature-selection/thresholding audit. Because it was added after prior external failure-mode work, it should be frozen as the next locked melanoma rescue-combo candidate and confirmed on any newly obtained independent controlled cohort.