# PD1-Like External Rescue Audit

This analysis addresses the weak strict melanoma PD1-like external performance for GSE145996 and PHS000452_LIU_LIKE_PRE.

Important boundary: the locked primary model remains the valid external-validation model. The transfer head is a secondary model-development rescue introduced after seeing the weakness, so it needs a fresh independent locked validation cohort before it can support a primary superiority claim.

## Strict RECIST Correction

Minor response (MR) is excluded from strict RECIST rather than counted as CR/PR response. This is an endpoint-rule correction, not model tuning.

## Pooled External Performance

- EcoNiche-Opt-HeuristicEcology-LockedPanel: AUROC=0.571, balanced_accuracy=0.542, ECE=0.261, n=116 (locked_primary).
- EcoNiche-Opt-PD1LikeTransferHead: AUROC=0.595, balanced_accuracy=0.512, ECE=0.145, n=116 (secondary_model_development_rescue).

## Fixed 0.5 Threshold Sensitivity

- EcoNiche-Opt-HeuristicEcology-LockedPanel: balanced_accuracy=0.531, sensitivity=0.429, specificity=0.633.
- EcoNiche-Opt-PD1LikeTransferHead: balanced_accuracy=0.627, sensitivity=0.554, specificity=0.700.

## Discovery LODO Selection

- EcoNiche-Opt-HeuristicEcology-LockedPanel: inner mean AUROC=0.710, inner min AUROC=0.648, status=locked_primary.
- EcoNiche-Opt-PD1LikeTransferHead: inner mean AUROC=0.660, inner min AUROC=0.544, status=secondary_model_development_rescue.

## Claim Boundary

Allowed: report the transfer head as a transparent secondary rescue that improves the weak PD1-like stress cohorts.
Not allowed: replace the locked external-validation claim with the transfer-head result unless a new independent external cohort validates it after this model is frozen.
