# Manual ICB Response Curation Audit

- Sample-level evidence table: `data\metadata\manual_response_curation.tsv` (360 labeled samples before PRJEB23709 split)
- PRJEB23709 evidence map: `data\external\PRJEB23709\gide_prjeb23709_response_map.tsv` (91 ENA runs; 73 pretreatment samples modeled)
- phs000452 TIGER stress-test map: `data\processed\bulk\PHS000452_LIU_LIKE_PRE.metadata.tsv` (121 Patient-like anti-PD1 samples modeled as secondary/stress-test)
- Cohort audit table: `results\curation\manual_curation_audit.tsv` (17 cohorts)
- Newly recovered cohorts: GSE140901, GSE145996, GSE165252, GSE67501, GSE93157, PRJEB23709, phs000452
- Usable modeled response cohorts after curation/stress-test: 14

## Newly recovered cohorts

| Accession | Samples | Genes | Labels | Decision | Evidence |
|---|---:|---:|---|---|---|
| GSE145996 | 14 | 19734 | `{"0": 8, "1": 6}` | usable_melanoma_response | MDPI Cancers supplementary workbook: SuppTable1 Patient Number/Best Response plus SuppTable2 Tumor ID number/RNASEQ; expression columns matched by MB tumor ID |
| PRJEB23709_PD1_PRE | 41 | 47708 | `{"0": 19, "1": 22}` | usable_melanoma_core_response | TIGER processed expression plus Gide/ENA response map; ENA run accession matched to sample_title/patient/timepoint; pretreatment anti-PD-1 monotherapy only |
| PRJEB23709_COMBO_PRE | 32 | 47708 | `{"0": 21, "1": 11}` | usable_secondary_combo_response | TIGER processed expression plus Gide/ENA response map; pretreatment anti-PD-1 plus anti-CTLA-4 kept separate from melanoma monotherapy primary |
| PHS000452_LIU_LIKE_PRE | 121 | 47708 | `{"0": 50, "1": 71}` | secondary_stress_test_only | TIGER Melanoma-phs000452 processed expression/clinical; Patient-like anti-PD1 rows retained, 32 IPI-like rows excluded from primary; baseline/timepoint is inferred from TIGER dataset context |
| GSE93157 | 65 | 770 | `{"0.0": 20, "1.0": 45}` | usable_pan_cancer_response | characteristics_ch1 best.resp/drug/biopsy; expression columns matched by GEO sample order |
| GSE67501 | 11 | 20819 | `{"0.0": 4, "1.0": 7}` | usable_pan_cancer_response | description (RCC sample ID); characteristics_ch1 response to anti-PD-1 nivolumab; GEO series design states pre-treatment tumors |
| GSE140901 | 24 | 784 | `{"0.0": 6, "1.0": 18}` | usable_pan_cancer_response | characteristics_ch1 best_response/clinical_benefit_response; title Sxx matched to expression column Sxx |
| GSE165252 | 32 | 61475 | `{"0.0": 12, "1.0": 20}` | usable_secondary_dynamic_confounded | characteristics_ch1 response; title sample_N baseline/on_treatment/resection; expression columns matched by GEO sample order |

## Held out / excluded cohorts

| Accession | Decision | Reason |
|---|---|---|
| GSE244982 | exclude_primary_response | Registry and sample metadata describe acquired resistance/progression mechanism groups rather than primary binary ICB response. |
| GSE121810 | exclude_primary_response | GEO metadata describes neoadjuvant/adjuvant pembrolizumab treatment-arm and survival setting in recurrent glioblastoma, not per-sample RECIST response. |
| GSE183924 | survival_only_not_response | Clinical workbook contains relapse/RFS/OS fields but no CR/PR/SD/PD or DCB/NDB response endpoint. |
| GSE123728 | missing_expression_file | No public expression matrix is present locally for this accession. |
| GSE165745 | missing_expression_file | No public expression matrix is present locally for this accession. |
| GSE122220 | missing_expression_file | No public expression matrix is present locally for this accession. |

## Label rule

The response endpoint is harmonized as `0=responder/clinical benefit` and `1=non-responder/no benefit`. For the primary RECIST-style endpoint, `SD` is conservatively grouped with non-response unless an explicit durable-clinical-benefit endpoint is used.
