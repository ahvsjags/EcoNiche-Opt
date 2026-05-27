# dbGaP Access Request Checklist

Purpose: obtain controlled melanoma checkpoint-blockade tumor RNA-seq and phenotype files for locked external validation of EcoNiche-Opt.

Required safeguards:

- Use approved institutional access only.
- Store controlled files outside Git and public archives.
- Import only through the documented schema in this package.
- Keep the cohort fully locked: no training, feature selection, thresholding, calibration, or model selection.

Targets:

- `phs001919`: Abril_Rodriguez_PD1_melanoma_dbGaP (https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001919.v1.p1)
- `phs002683`: MGH_Hacohen_checkpoint_melanoma_dbGaP (https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002683.v1.p1)

After access is approved:

1. Place raw controlled files under `data/raw/controlled/<accession>/` locally only.
2. Fill `controlled_assay_sample_manifest_template.tsv` and `controlled_clinical_annotation_template.tsv`.
3. Run the registered preprocessing/import path and then the locked external scoring command.
4. Record missing labels or missing baseline evidence as `RESULT_PENDING`, not as negative or imputed validation evidence.