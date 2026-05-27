# EGA Access Request Checklist

Purpose: obtain controlled melanoma checkpoint-blockade tumor RNA-seq and phenotype files for locked external validation of EcoNiche-Opt.

Required safeguards:

- Use approved institutional access only.
- Store controlled files outside Git and public archives.
- Import only through the documented schema in this package.
- Keep the cohort fully locked: no training, feature selection, thresholding, calibration, or model selection.

Targets:

- `EGAS00001001552`: Lee_Rizos_PD1_melanoma_EGA (https://ega-archive.org/studies/EGAS00001001552)

After access is approved:

1. Place raw controlled files under `data/raw/controlled/<accession>/` locally only.
2. Fill `controlled_assay_sample_manifest_template.tsv` and `controlled_clinical_annotation_template.tsv`.
3. Run the registered preprocessing/import path and then the locked external scoring command.
4. Record missing labels or missing baseline evidence as `RESULT_PENDING`, not as negative or imputed validation evidence.