# Controlled External Validation Access Package

This package supports the strict melanoma external-validation target for EcoNiche-Opt.
It does not contain controlled data and does not count any controlled cohort as validation evidence.

## Why This Package Exists

The current public strict melanoma external AUROC is below the top-tier target. The realistic path to a stronger locked external test is to obtain independent controlled melanoma checkpoint-blockade tumor RNA-seq cohorts, then score them without refitting.

## Locked Boundary

Controlled external cohort labels must never enter model training, feature selection, threshold selection, calibration, or candidate selection. They may only be used after the model, genes, weights, thresholds, and scoring script are frozen.

## Targets

- `EGAS00001001552` (EGA): Lee_Rizos_PD1_melanoma_EGA; request: https://ega-archive.org/studies/EGAS00001001552
- `phs001919` (dbGaP): Abril_Rodriguez_PD1_melanoma_dbGaP; request: https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001919.v1.p1
- `phs002683` (dbGaP): MGH_Hacohen_checkpoint_melanoma_dbGaP; request: https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002683.v1.p1

## Expected Local Import Files

- `expression.tsv`: sample-by-gene expression matrix.
- `assay_sample_manifest.tsv`: file/sample/patient/platform/timing manifest.
- `clinical_annotation.tsv`: patient-level response, treatment, timing, and provenance table.

## Locked Scoring Command

```bash
python -m econiche_opt.cli score-locked-validation \
  --package-dir deliverables/prospective_validation \
  --expression data/raw/controlled/<accession>/expression.tsv \
  --sample-manifest data/raw/controlled/<accession>/assay_sample_manifest.tsv \
  --clinical-annotation data/raw/controlled/<accession>/clinical_annotation.tsv \
  --out-dir results/controlled_external_validation/<accession>
```

## Claim Rule

If expression, baseline timing, or response provenance is incomplete, report `RESULT_PENDING`. Do not impute controlled external validation.
