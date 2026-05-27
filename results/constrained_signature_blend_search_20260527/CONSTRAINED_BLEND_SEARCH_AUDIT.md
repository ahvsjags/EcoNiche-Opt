# Constrained Signature Blend Search Audit

Predeclared constrained blend candidates tested: 60. Candidate selection uses only inner LODO training folds or discovery-only inner LODO for external groups.

- primary_recist / melanoma_core_high_evidence: AUROC=0.562, AUPRC=0.435, balanced accuracy=0.521, ECE=0.149, selected=module_plus_APM_IPRES,module_plus_CXCL9_TIDE_dysfunction,module_plus_TIDE_dysfunction_IPRES.
- primary_recist / melanoma_recist_supported_primary: AUROC=0.570, AUPRC=0.495, balanced accuracy=0.587, ECE=0.107, selected=module_plus_CXCL9_IPRES,module_plus_CYT_IPRES,module_plus_CYT_alpha0.8,module_plus_TIDE_dysfunction_IPRES.
- strict_recist / strict_cbio_liu_plus_gse145996: AUROC=0.559, AUPRC=0.569, balanced accuracy=0.541, ECE=0.043, selected=module_plus_CXCL9_IPRES.
- strict_recist / strict_melanoma_pd1_like_external: AUROC=0.557, AUPRC=0.566, balanced accuracy=0.523, ECE=0.064, selected=module_plus_CXCL9_IPRES.