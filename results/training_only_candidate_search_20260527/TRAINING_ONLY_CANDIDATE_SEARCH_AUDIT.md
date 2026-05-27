# Training-only Melanoma Candidate Search Audit

This registered analysis asks whether a leakage-safe candidate-selection rule can improve the primary melanoma and strict external evidence without using external labels for feature selection, thresholding, calibration, or model selection.

## Selection Rule

- Candidate families are restricted to fixed module composites and sparse module-level logistic models.
- For primary LODO, each holdout fold selects a candidate from the remaining training cohorts only.
- For strict external testing, the selected candidate is chosen from GSE91061, GSE78220 and PRJEB23709_PD1_PRE only; GSE145996 and PHS000452_LIU_LIKE_PRE labels are read only after the candidate and threshold rule are fixed.

## Primary LODO Summary

- melanoma_core_high_evidence: selected=module_prior_composite; n=117; AUROC=0.705; AUPRC=0.558; balanced accuracy=0.632; ECE=0.235.
- melanoma_recist_supported_primary: selected=module_prior_composite; n=131; AUROC=0.685; AUPRC=0.580; balanced accuracy=0.613; ECE=0.239.

## Strict External Summary

- GSE145996+PHS000452_LIU_LIKE_PRE: selected=module_prior_composite; n=116; AUROC=0.572; AUPRC=0.600; balanced accuracy=0.531; ECE=0.256.

## Interpretation

The training-only search does not reach the strict external AUROC >=0.70 target. This rules out a simple no-leakage candidate-selection fix and prioritizes either new independent melanoma tumor-tissue data or a materially different training-only representation.

## Discovery-only External Candidate Ranking

The top discovery-only candidate was module_prior_composite with inner mean AUROC=0.710, inner mean AUPRC=0.682, and selection score=0.776.