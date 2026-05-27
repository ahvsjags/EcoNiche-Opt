# Tumor-Tissue Assay SOP for Independent EcoNiche-Opt Validation

## Specimen Requirement

EcoNiche-Opt is currently locked for pretreatment tumor-tissue transcriptomes. Blood-derived samples cannot directly validate this model because the score was trained and calibrated on tumor tissue and includes tumor-microenvironment states such as antigen presentation, myeloid suppression, stromal exclusion, TRM/TLS, and T/NK effector activity.

## Tissue and Pathology QC

1. Confirm melanoma diagnosis and sampling date relative to first ICB dose.
2. Prefer pretreatment FFPE or fresh-frozen tumor tissue with pathology review.
3. Record tumor content, necrosis percentage, macrodissection status, specimen block/slide ID, and anatomic site.
4. Extract RNA under the local certified workflow and record RIN or DV200, input amount, and assay batch.
5. Quantify the locked panel genes in `locked_panel_genes.tsv`; failed genes must be reported, not silently replaced.

## Expression Matrix

Rows must be sample IDs and columns must be HGNC gene symbols. Values should be normalized within the assay workflow and comparable across samples within the validation freeze. If the platform exports genes-by-samples, transpose before scoring and retain the raw export in the local audit archive.

## Locked Scoring

Use `locked_scoring_spec.json`, verify the SHA256 hash in `locked_scoring_spec.sha256`, score all QC-passing samples once, and apply the endpoint-specific locked threshold. Do not retrain or recalibrate on the independent validation labels.
