# Claim and External Validation Audit

Generated for the 2026-05-08 claim-strengthening pass. Evidence comes from registered scripts and result tables, not hand-edited result files.

## What Was Strengthened

1. FDR claim strength was upgraded from mostly point-estimate superiority to a predeclared eight-signature family test.
2. Locked external and clinical-assay panel transfer validation now has a reproducible runner and status file instead of a RESULT_PENDING placeholder.
3. A frozen prospective qPCR/NanoString validation package was generated for future independent clinical samples.

## Primary Melanoma Claim Evidence

Source: `results/claim_strengthening/strong_signature_family_omnibus.tsv`

- `primary_recist / melanoma_core_high_evidence`: EcoNiche-Opt-HeuristicEcology AUROC 0.705 versus eight-signature mean AUROC 0.632, delta 0.073, two-sided FDR q=0.002.
- `primary_recist / melanoma_recist_supported_primary`: EcoNiche-Opt-HeuristicEcology AUROC 0.685 versus eight-signature mean AUROC 0.620, delta 0.064, two-sided FDR q=0.002.

Allowed claim: EcoNiche-Opt significantly outperforms the predeclared eight-signature family in the primary melanoma benchmark strata.

Not allowed claim: EcoNiche-Opt significantly outperforms every individual signature or all existing models.

## Locked External and Panel Validation

Source: `results/locked_external_panel_validation/LOCKED_EXTERNAL_PANEL_VALIDATION_AUDIT.md`

Discovery-only threshold cohorts:

- `GSE91061`
- `GSE78220`
- `PRJEB23709_PD1_PRE`

Locked external/panel cohorts:

- `GSE145996`
- `PHS000452_LIU_LIKE_PRE`
- `PRJEB23709_COMBO_PRE`
- `GSE93157`
- `GSE140901`

External family-level results:

- `clinical_benefit / all_locked_external_and_panel`: target AUROC 0.604 versus eight-signature mean AUROC 0.566, two-sided FDR q=0.027.
- `strict_recist / all_locked_external_and_panel`: target AUROC 0.610 versus eight-signature mean AUROC 0.570, two-sided FDR q=0.031 after correcting strict RECIST to exclude minor response (MR).
- `primary_recist / all_locked_external_and_panel`: target AUROC 0.577 versus eight-signature mean AUROC 0.549, one-sided FDR q=0.035 and two-sided FDR q=0.071.

External per-cohort highlights:

- `PRJEB23709_COMBO_PRE`: primary RECIST AUROC 0.740 and strict RECIST AUROC 0.762.
- `GSE140901` NanoString transfer: primary RECIST AUROC 0.759 and strict RECIST AUROC 0.833.
- `GSE145996` and `PHS000452_LIU_LIKE_PRE` remain modest on locked strict PD-1-like external validation. After strict RECIST correction, the locked pooled AUROC is 0.571.

## PD1-Like Rescue Analysis

Source: `results/pd1_like_external_rescue/PD1_LIKE_EXTERNAL_RESCUE_AUDIT.md`

The strict RECIST rule was corrected so minor response (MR) is excluded rather than counted as CR/PR response. A secondary `EcoNiche-Opt-PD1LikeTransferHead` was added as a transparent model-development rescue for the weak PD1-like stress cohorts.

- Locked primary model on `GSE145996 + PHS000452_LIU_LIKE_PRE`: AUROC 0.571, balanced accuracy 0.542, ECE 0.261.
- Secondary transfer head on the same stress set: AUROC 0.595, discovery-threshold balanced accuracy 0.512, ECE 0.145.
- Secondary transfer head with fixed probability threshold 0.5: balanced accuracy 0.627, sensitivity 0.554, specificity 0.700.

Allowed claim: the secondary transfer head improves the weak PD1-like stress cohorts and gives a better-calibrated rescue analysis.

Not allowed claim: the transfer head is a new locked external-validation success unless it is frozen and validated on a fresh independent cohort.

## Prospective Validation Boundary

Source: `deliverables/prospective_validation/`

The package contains:

- `locked_scoring_spec.json`
- `locked_panel_genes.tsv`
- `assay_sample_manifest_template.tsv`
- `clinical_annotation_template.tsv`
- `prospective_validation_protocol.md`
- `statistical_analysis_plan.md`

Allowed claim: a prospective-ready locked assay package is provided.

Not allowed claim: prospective clinical validation has been completed.

## Verification

- `python -m py_compile scripts/analysis/run_claim_strengthening.py scripts/analysis/run_locked_external_panel_validation.py scripts/reporting/make_prospective_validation_package.py scripts/benchmark/run_locked_external.py`
- `python -m pytest tests/test_claim_gate.py tests/test_endpoint_modules.py -q`
- `python -m pytest -q`
- `python -m econiche_opt.cli validate-goals --goal-file docs/goal_status.yml`
- `python -m econiche_opt.cli validate-project --mode demo`
