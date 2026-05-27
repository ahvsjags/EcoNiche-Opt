# Clinical Partner Intake Checklist

Use this checklist before accepting a new independent cohort for EcoNiche-Opt validation.

## Required Cohort Definition

- Disease: melanoma.
- Specimen: pretreatment tumor tissue, not blood, plasma, serum, or PBMC.
- Treatment: anti-PD-1 or anti-PD-1-based immune-checkpoint blockade.
- Endpoint evidence: RECIST 1.1 CR/PR/SD/PD with scan dates, or a prospectively defined DCB/NDB endpoint.
- Preferred sample size: 50-100 independent patients for a credible first validation; smaller sets are pilot/feasibility only.

## Files to Request From the Clinical Team

- De-identified sample manifest using `assay_sample_manifest_template.tsv`.
- De-identified clinical annotation using `clinical_annotation_template.tsv`.
- Expression matrix with samples as rows and HGNC gene symbols as columns.
- Pathology QC summary: tumor content, necrosis, macrodissection, FFPE/fresh-frozen status, RNA QC.
- Response evidence: source eCRF, tumor board export, RECIST worksheet, or annotated clinical Excel.
- Data dictionary explaining every non-standard field.

## Pre-Scoring Gate

Do not score the cohort until subject IDs, sample IDs, baseline status, therapy dates, assay platform, response labels, and exclusion decisions are complete. Validation labels must not be used for feature selection, calibration, threshold tuning, or model refitting.
