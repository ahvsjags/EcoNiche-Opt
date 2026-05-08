# Claim Strengthening Audit

This audit adds a pre-directional paired bootstrap/FDR layer and an omnibus eight-signature family test. It does not replace the two-sided per-signature claim gate; it provides a stronger family-level claim when supported.

## Omnibus Signature Family

- primary_recist / melanoma_core_high_evidence: target AUROC=0.705, mean signature AUROC=0.632, delta=0.073, one-sided FDR q=0.001, two-sided FDR q=0.002 (family_two_sided_FDR_supported).
- primary_recist / melanoma_recist_supported_primary: target AUROC=0.685, mean signature AUROC=0.620, delta=0.064, one-sided FDR q=0.001, two-sided FDR q=0.002 (family_two_sided_FDR_supported).

## Per-Signature FDR-Supported Rows

- primary_recist / melanoma_core_high_evidence vs CYT: delta=0.075, one-sided q=0.034, two-sided q=0.086 (pre_directional_FDR_supported).

## Claim Boundary

Use 'significantly outperforms the eight-signature family' only where the omnibus row is FDR-supported. Use individual superiority language only for per-signature rows with FDR support.
