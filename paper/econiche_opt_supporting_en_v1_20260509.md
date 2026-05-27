# Supplementary Information

## EcoNiche-Opt: a locked immune-ecology transcriptomic score for multicohort prediction of immune checkpoint blockade response

Pengyuan Xu1,3, Guang Yang2,3, Moyan Li3*

1 Department of Materials Science and Engineering, Monash University, Clayton, VIC 3800, Australia.

2 School of Economics and Management, China University of Mining and Technology, Xuzhou 221116, Jiangsu, China.

3 Hong Kong University of Science and Technology (Guangzhou), Guangzhou 510000, China.

*Correspondence: moyanli@hkust-gz.edu.cn

## Inventory

- Supplementary Tables 1-24
- Supplementary References
- Supplementary Figures 1-10
- Supplementary Methods

## Supplementary Tables

**Supplementary Table 1 | Data registry roles.**  
Lists all registered data sources, cohort roles, cancer types, treatment types, assay platforms, data-use categories, and pipeline status. This table supports Figure 1a-b.

**Supplementary Table 2 | Data access status.**  
Records the access status of each data source, including public, controlled, and unknown access. This table supports Figure 1d and ensures that controlled data are not replaced by simulated outputs.

**Supplementary Table 3 | Expression QC.**  
Reports n_samples, n_genes, source file, and QC status for each expression matrix. This table supports Figure 1b and Supplementary Figure 2.

**Supplementary Table 4 | Manual curation evidence.**  
Records response labels, sample IDs, patient IDs, baseline/on-treatment status, treatment regimen, and evidence source. This table supports the main-text description of manual response-label curation.

**Supplementary Table 5 | External cohort availability assessment.**  
Records availability, missing data, and access boundaries for candidate external cohorts. This table supports the external-validation design in Figures 1 and 5.

**Supplementary Table 6 | Endpoint label sensitivity.**  
Lists n_used, n_dropped, responders, and nonresponders for each cohort under strict RECIST, primary RECIST, and clinical-benefit endpoints. This table supports Figure 4a and Supplementary Figure 5.

**Supplementary Table 7 | Module gene sets.**  
Lists EcoNiche-Opt prior modules, gene symbols, module weights, and score directions. This table supports Figure 2a, Figure 6d, and Figure 7b.

**Supplementary Table 8 | Ecological edges.**  
Lists optimized ecological interaction edges, including source_state, target_state, gene_a, gene_b, and edge_class. This table supports Figure 2d and Figure 6b.

**Supplementary Table 9 | Optimizer history.**  
Records optimizer objective, best score, best AUROC, best AUPRC, best ECE, and model variant. This table supports Figure 2e and Supplementary Figure 6.

**Supplementary Table 10 | Melanoma benchmark summary.**  
Summarizes pooled AUROC, AUPRC, balanced accuracy, ECE, and fold statistics for EcoNiche-Opt and baseline signatures in the primary melanoma benchmark. This table supports Figure 3.

**Supplementary Table 11 | Signature family FDR.**  
Reports the strong-signature family omnibus test, including target AUROC, family mean AUROC, mean delta, confidence interval, P value, and FDR q. This table supports Figure 3b and Figure 3f.

**Supplementary Table 12 | LODO metrics.**  
Lists leave-one-dataset-out metrics by holdout cohort, endpoint, model, and stratum. This table supports Figure 3d, Figure 4e, and Supplementary Figure 3.

**Supplementary Table 13 | Aligned locked-panel ablation.**  
Reports ablations matched to the locked module-panel scoring family used for the main melanoma result, including no-response-module, no-resistance-module, unsigned-direction, IFNG-only, equal-weight response-module, and calibrated variants. The table reports full AUROC, ablation AUROC, delta AUROC, ECE, confidence interval, FDR q value, and claim level. This table supports Figure 4b.

**Supplementary Table 14 | Decision curve.**  
Lists threshold probability, model net benefit, treat-all net benefit, and treat-none net benefit. This table supports Figure 3e and Figure 4c.

**Supplementary Table 15 | Locked external metrics.**  
Reports AUROC, AUPRC, balanced accuracy, ECE, threshold, and calibration metrics for each endpoint and model in locked external and panel-transfer cohorts. This table supports Figure 5b-c and Supplementary Figure 7.

**Supplementary Table 16 | NanoString panel transfer.**  
Reports n_module_genes, n_available_genes, coverage_fraction, and available genes for each module in NanoString panel-transfer cohorts. This table supports Figure 5d.

**Supplementary Table 17 | PD1-like rescue.**  
Reports AUROC, AUPRC, balanced accuracy, ECE, and threshold sensitivity for locked panel scoring and PD1LikeTransferHead in GSE145996, PHS000452_LIU_LIKE_PRE, and the pooled stress cohort. This table supports Figure 5e-f and Supplementary Figure 8.

**Supplementary Table 18 | Single-cell enrichment.**  
Reports cell_type, state, mean, median, and count for single-cell module localization. This table supports Figure 6a and Supplementary Figure 9.

**Supplementary Table 19 | Perturbation hypotheses.**  
Lists perturbation_name, target_gene, mechanism, reversal_score, depmap_score, DGIdb evidence, priority_score, and status. This table supports Figure 6e.

**Supplementary Table 20 | Algorithmic reproducibility and independent-scoring boundary.**  
Summarizes the reproducible analysis boundary from public ICB cohort curation, endpoint harmonization, signed-rank ecological modules, interaction features, heuristic optimization, threshold locking, external scoring, coverage QC, signature-family benchmarking, and the FDR-aware claim gate to open scoring code. This table clarifies which information is learned in the discovery layer, which information is read only during external evaluation, and how EcoNiche-Opt becomes an independently scorable method without using external labels for model selection.

**Supplementary Table 21 | GPU biological-prior rescue combo.**  
Reports the frozen GPU lipid/PI3K rescue-combo candidate, biological prior, transform policy, GPU device, primary melanoma LODO metrics, strict melanoma external metrics, family-level delta, FDR q value, claim level, and no-leakage selection boundary. This table supports Figure 5e.

**Supplementary Table 22 | GPU rescue component ablation.**  
Compares the selected GPU lipid/PI3K rescue combo with the ablated MAP4K1/TBX3/AXL base rescue axis in primary melanoma LODO and strict melanoma external settings. It reports AUROC, AUPRC, balanced accuracy, paired bootstrap intervals, P values, FDR q values, and claim levels. This table supports Figure 5e and the component-ablation claim boundary.

**Supplementary Table 23 | cBioPortal GPU rescue external validation.**  
Reports cBioPortal Liu/DFCI and cBioPortal Liu plus GSE145996 strict RECIST metrics for the frozen GPU lipid/PI3K rescue combo, including AUROC, AUPRC, balanced accuracy, ECE, Brier score, sensitivity, specificity, PPV, NPV, calibration slope/intercept, and no-leakage scoring boundary. This table supports Figure 5f.

**Supplementary Table 24 | GPU lipid/PI3K pair rescue.**  
Reports the component-dominant PLA2G2D/PIK3CD lipid/PI3K pair rescue selected by primary melanoma LODO with a balanced-accuracy guardrail. It includes the locked candidate formula, GPU device, primary AUROC/AUPRC/balanced accuracy, strict melanoma external metrics, cBioPortal Liu/DFCI metrics, family-level delta, FDR q value, claim level, and no-leakage selection boundary. This table supports Figure 5e-f.

## Supplementary References

1. Pardoll DM. The blockade of immune checkpoints in cancer immunotherapy. Nat Rev Cancer. 2012;12:252-264. doi:10.1038/nrc3239.
2. Eisenhauer EA, Therasse P, Bogaerts J, et al. New response evaluation criteria in solid tumours: revised RECIST guideline (version 1.1). Eur J Cancer. 2009;45:228-247. doi:10.1016/j.ejca.2008.10.026.
3. Ayers M, Lunceford J, Nebozhyn M, et al. IFN-gamma-related mRNA profile predicts clinical response to PD-1 blockade. J Clin Invest. 2017;127:2930-2940. doi:10.1172/JCI91190.
4. Jiang P, Gu S, Pan D, et al. Signatures of T cell dysfunction and exclusion predict cancer immunotherapy response. Nat Med. 2018;24:1550-1558. doi:10.1038/s41591-018-0136-1.
5. Rooney MS, Shukla SA, Wu CJ, Getz G, Hacohen N. Molecular and genetic properties of tumors associated with local immune cytolytic activity. Cell. 2015;160:48-61. doi:10.1016/j.cell.2014.12.033.
6. Hugo W, Zaretsky JM, Sun L, et al. Genomic and transcriptomic features of response to anti-PD-1 therapy in metastatic melanoma. Cell. 2016;165:35-44. doi:10.1016/j.cell.2016.02.065.
7. Riaz N, Havel JJ, Makarov V, et al. Tumor and microenvironment evolution during immunotherapy with nivolumab. Cell. 2017;171:934-949.e16. doi:10.1016/j.cell.2017.09.028.
8. Auslander N, Zhang G, Lee JS, et al. Robust prediction of response to immune checkpoint blockade therapy in metastatic melanoma. Nat Med. 2018;24:1545-1549. doi:10.1038/s41591-018-0157-9.
9. Gide TN, Quek C, Menzies AM, et al. Distinct immune cell populations define response to anti-PD-1 monotherapy and anti-PD-1/anti-CTLA-4 combined therapy. Cancer Cell. 2019;35:238-255.e6. doi:10.1016/j.ccell.2019.01.003.
10. Liu D, Schilling B, Liu D, et al. Integrative molecular and clinical modeling of clinical outcomes to PD1 blockade in patients with metastatic melanoma. Nat Med. 2019;25:1916-1927. doi:10.1038/s41591-019-0654-5.
11. Jerby-Arnon L, Shah P, Cuoco MS, et al. A cancer cell program promotes T cell exclusion and resistance to checkpoint blockade. Cell. 2018;175:984-997.e24. doi:10.1016/j.cell.2018.09.006.
12. Du K, Wei S, Wei Z, et al. Pathway signatures derived from on-treatment tumor specimens predict response to anti-PD1 blockade in metastatic melanoma. Nat Commun. 2021;12:6023. doi:10.1038/s41467-021-26299-4.
13. Prat A, Navarro A, Pare L, et al. Immune-related gene expression profiling after PD-1 blockade in non-small cell lung carcinoma, head and neck squamous cell carcinoma, and melanoma. Cancer Res. 2017;77:3540-3550. doi:10.1158/0008-5472.CAN-16-3556.
14. Rose TL, Weir WH, Mayhew GM, et al. Fibroblast growth factor receptor 3 alterations and response to immune checkpoint inhibition in metastatic urothelial cancer: a real world experience. Br J Cancer. 2021;125:1251-1260. doi:10.1038/s41416-021-01488-6.
15. Hwang S, Kwon AY, Jeong JY, et al. Immune gene signatures for predicting durable clinical benefit of anti-PD-1 immunotherapy in patients with non-small cell lung cancer. Sci Rep. 2020;10:643. doi:10.1038/s41598-019-57218-9.
16. Mariathasan S, Turley SJ, Nickles D, et al. TGF-beta attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. Nature. 2018;554:544-548. doi:10.1038/nature25501.
17. van den Ende T, de Clercq NC, van Berge Henegouwen MI, et al. Neoadjuvant chemoradiotherapy combined with atezolizumab for resectable esophageal adenocarcinoma: a single-arm phase II feasibility trial (PERFECT). Clin Cancer Res. 2021;27:3351-3359. doi:10.1158/1078-0432.CCR-20-4443.
18. Thorsson V, Gibbs DL, Brown SD, et al. The immune landscape of cancer. Immunity. 2018;48:812-830.e14. doi:10.1016/j.immuni.2018.03.023.
19. Tirosh I, Izar B, Prakadan SM, et al. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. Science. 2016;352:189-196. doi:10.1126/science.aad0501.
20. Sade-Feldman M, Yizhak K, Bjorgaard SL, et al. Defining T cell states associated with response to checkpoint immunotherapy in melanoma. Cell. 2019;176:404. doi:10.1016/j.cell.2018.12.034.
21. Collins GS, Reitsma JB, Altman DG, Moons KGM. Transparent reporting of a multivariable prediction model for individual prognosis or diagnosis (TRIPOD): the TRIPOD statement. Ann Intern Med. 2015;162:55-63. doi:10.7326/M14-0697.
22. Collins GS, Moons KGM, Dhiman P, et al. TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods. BMJ. 2024;385:e078378. doi:10.1136/bmj-2023-078378.
23. DeLong ER, DeLong DM, Clarke-Pearson DL. Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach. Biometrics. 1988;44:837-845. doi:10.2307/2531595.
24. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B. 1995;57:289-300. doi:10.1111/j.2517-6161.1995.tb02031.x.
25. Vickers AJ, Elkin EB. Decision curve analysis: a novel method for evaluating prediction models. Med Decis Making. 2006;26:565-574. doi:10.1177/0272989X06295361.
26. Guo C, Pleiss G, Sun Y, Weinberger KQ. On calibration of modern neural networks. Proceedings of the 34th International Conference on Machine Learning. 2017;70:1321-1330. https://proceedings.mlr.press/v70/guo17a.html.
27. Subramanian A, Narayan R, Corsello SM, et al. A next generation Connectivity Map: L1000 platform and the first 1,000,000 profiles. Cell. 2017;171:1437-1452.e17. doi:10.1016/j.cell.2017.10.049.
28. Freshour SL, Kiwala S, Cotto KC, et al. Integration of the Drug-Gene Interaction Database (DGIdb 4.0) with open crowdsource efforts. Nucleic Acids Res. 2021;49:D1144-D1151. doi:10.1093/nar/gkaa1084.
29. Tsherniak A, Vazquez F, Montgomery PG, et al. Defining a cancer dependency map. Cell. 2017;170:564-576.e16. doi:10.1016/j.cell.2017.06.010.

## Supplementary Figures

**Supplementary Figure 1 | Cohort curation and data access.**  
Extends the cohort-curation component of Figure 1 by showing sample size, data access status, cohort role, and inclusion boundary for each registered data source. This figure corresponds to Supplementary Tables 1-5 and documents the traceable evidence for cohort registration, access-status assessment, manual curation, and external cohort availability.

**Supplementary Figure 2 | Platform and gene coverage QC.**  
Extends Figure 1f by showing gene coverage across expression platforms, cohorts, and ecological modules. The figure identifies which cohorts can fully compute IFN/T inflamed, CD8 cytotoxic, checkpoint/exhaustion, antigen presentation, myeloid suppression, stromal exclusion, and TRM/TLS modules and which cohorts show panel gene loss. It corresponds to Supplementary Tables 3, 7, and 16.

**Supplementary Figure 3 | Discovery benchmark detailed performance.**  
Extends the melanoma benchmark in Figure 3 by showing discovery/LODO AUROC, AUPRC, balanced accuracy, and holdout cohort heterogeneity. It supports the detailed numerical interpretation of primary melanoma performance and corresponds to Supplementary Tables 10 and 12.

**Supplementary Figure 4 | Full baseline comparison.**  
Shows EcoNiche-Opt against the broader baseline set, including expanded signatures or model variants beyond the strong-signature family. It supports the conclusion that EcoNiche-Opt's advantage is concentrated in family-level comparison and ecological-module structure and corresponds to Supplementary Tables 10-12.

**Supplementary Figure 5 | Endpoint sensitivity.**  
Extends Figure 4a by showing how strict RECIST, primary RECIST, and clinical-benefit endpoints affect sample inclusion, response balance, and model evaluation. It corresponds to Supplementary Table 6 and demonstrates why endpoint harmonization is required.

**Supplementary Figure 6 | Aligned ablation and optimizer diagnostics.**  
Extends Figures 2 and 4b by showing aligned locked-panel component ablation together with optimizer diagnostics for the ecological graph search space. It corresponds to Supplementary Tables 8, 9, and 13.

**Supplementary Figure 7 | Expanded locked external validation.**  
Extends Figure 5b-c by showing performance of locked external and panel-transfer cohorts under primary RECIST, strict RECIST, and clinical-benefit endpoints. It corresponds to Supplementary Table 15.

**Supplementary Figure 8 | PD1-like and GPU biological-prior rescue.**  
Extends Figure 5e-f by showing strict RECIST PD1-like stress analysis in GSE145996 and PHS000452_LIU_LIKE_PRE. It includes the original locked-panel stress analysis, the frozen GPU lipid/PI3K rescue combo, the balanced-accuracy-guarded PLA2G2D/PIK3CD pair rescue, component ablation against the MAP4K1/TBX3/AXL base axis, and cBioPortal Liu/DFCI source cross-check evidence. It corresponds to Supplementary Tables 17 and 21-24.

**Supplementary Figure 9 | Single-cell and ecological mechanism.**  
Extends Figure 6a-b by showing module-score distributions across single-cell cell types and the interaction structure among ecological states. It corresponds to Supplementary Table 18 and supports the biological interpretability of the model.

**Supplementary Figure 10 | Reproducibility path for EcoNiche-Opt.**  
Extends Figure 7 by showing how cohort curation, locked scoring rules, sample-level traceability, prespecified endpoint analysis, external scoring boundaries, and open reproducibility interfaces support independent rescoring. It corresponds to Supplementary Table 20.

## Supplementary Methods

### Cohort registration and access control

All data sources first underwent cohort registration and then proceeded to expression QC, label curation, and model evaluation. Controlled or unknown-access sources retain their access status and are not replaced with simulated data. Analyses depending on inaccessible data are reported as access-restricted rather than filled with substitutes.

### Label traceability

Response-label curation requires each sample_id to be traceable to subject_id, treatment, timepoint, response_raw, and evidence source. CR/PR/SD/PD, R/NR, and DCB/NDB labels are mapped to three endpoints while preserving the original response_raw annotation.

### Leakage-safe evaluation

Training, feature selection, gene-direction estimation, thresholding, calibration, and model selection use training data only. Holdout, locked external, and panel-transfer cohorts do not participate in model selection.

### Claim gate

Superiority claims are used only when paired bootstrap or DeLong-compatible comparisons pass FDR correction. Family-level and individual-signature claims are reported separately. Perturbation rankings are hypothesis-only outputs.

### Reporting and reproducibility

The study follows prediction-model reporting logic for data sources, sample inclusion, endpoint definitions, model formulas, training/validation splits, performance metrics, calibration, decision-curve analysis, external validation, code availability, and data availability. The open reproducibility layer documents how the same frozen algorithm preserves module scores, interaction features, thresholds, and probability outputs across Python/R environments, allowing readers to distinguish discovery training, external scoring, and performance evaluation. Figures 1-7, Supplementary Figures 1-10, and Supplementary Tables 1-24 form the reporting evidence map.
