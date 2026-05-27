# Strict Melanoma PD1-like External Claim Gate

Strict melanoma PD1-like external gate: modest_point_estimate_only for strict_recist; n=116, AUROC=0.571, family mean AUROC=0.542, delta=0.030, q=0.309. This supports locked refitting-free external scoring, not a high-strength clinical validation claim.

## What this gate checks

- Strict external evidence is restricted to GSE145996 and PHS000452_LIU_LIKE_PRE.
- The primary endpoint for this gate is strict RECIST.
- The gate reads locked external outputs only; it does not refit, recalibrate, select thresholds, or change labels.
- A primary external superiority claim requires AUROC >= 0.70 and FDR-supported improvement over the predeclared signature family.

## Evidence summary

- pd1_like_transfer_head_pooled_strict_recist: secondary_model_development_rescue_not_primary_claim; endpoint=strict_recist; cohort_set=GSE145996+PHS000452_LIU_LIKE_PRE; n=116; AUROC=0.595; q=NA.
- strict_family_clinical_benefit: modest_point_estimate_only; endpoint=clinical_benefit; cohort_set=GSE145996+PHS000452_LIU_LIKE_PRE; n=135; AUROC=0.586; q=0.309.
- strict_family_primary_recist: modest_point_estimate_only; endpoint=primary_recist; cohort_set=GSE145996+PHS000452_LIU_LIKE_PRE; n=135; AUROC=0.531; q=0.540.
- strict_family_strict_recist: modest_point_estimate_only; endpoint=strict_recist; cohort_set=GSE145996+PHS000452_LIU_LIKE_PRE; n=116; AUROC=0.571; q=0.309.
- strict_cohort_clinical_benefit_GSE145996: cohort_support_moderate; endpoint=clinical_benefit; cohort_set=GSE145996; n=14; AUROC=0.622; q=NA.
- strict_cohort_clinical_benefit_PHS000452_LIU_LIKE_PRE: cohort_support_weak_or_modest; endpoint=clinical_benefit; cohort_set=PHS000452_LIU_LIKE_PRE; n=121; AUROC=0.583; q=NA.
- strict_cohort_primary_recist_GSE145996: cohort_support_weak_or_modest; endpoint=primary_recist; cohort_set=GSE145996; n=14; AUROC=0.521; q=NA.
- strict_cohort_primary_recist_PHS000452_LIU_LIKE_PRE: cohort_support_weak_or_modest; endpoint=primary_recist; cohort_set=PHS000452_LIU_LIKE_PRE; n=121; AUROC=0.535; q=NA.
- strict_cohort_strict_recist_GSE145996: cohort_support_weak_or_modest; endpoint=strict_recist; cohort_set=GSE145996; n=13; AUROC=0.575; q=NA.
- strict_cohort_strict_recist_PHS000452_LIU_LIKE_PRE: cohort_support_weak_or_modest; endpoint=strict_recist; cohort_set=PHS000452_LIU_LIKE_PRE; n=103; AUROC=0.573; q=NA.

## Allowed wording

The manuscript may state that strict melanoma PD1-like external cohorts were scored with a locked, refitting-free rule and showed modest point-estimate support. It should include cohort set, endpoint, n, AUROC, family mean AUROC, FDR q value, and calibration metrics.

## Blocked wording

Do not state that the strict melanoma external validation is clinically strong, AUROC >=0.70, FDR-supported, or prospectively validated unless a newly independent tumor-tissue cohort proves those claims after the score is frozen.

## Required next evidence

The next decisive dataset is an independent pretreatment melanoma tumor-tissue cohort treated with anti-PD1/anti-PD1-based therapy, RECIST CR/PR/SD/PD labels, and RNA-seq/NanoString/qPCR measurement of the locked panel, ideally n>=50-100.
