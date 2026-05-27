# Secondary Melanoma External Sensitivity

This audit scores fixed melanoma rescue-head candidates on secondary public melanoma cohorts without refitting, recalibration, feature selection, or threshold selection on external labels.

## current_strict_pd1_like

- `EcoNiche-Opt-MAP4K1-TBX3AXL-RobustZStressOnly`: cohorts=GSE145996,PHS000452_LIU_LIKE_PRE; AUROC=0.686; AUPRC=0.690; BA=0.624; boundary=current_external_stress_screen_not_locked_selection
- `EcoNiche-Opt-MAP4K1-TBX3AXL-ZScoreStressOnly`: cohorts=GSE145996,PHS000452_LIU_LIKE_PRE; AUROC=0.679; AUPRC=0.684; BA=0.602; boundary=current_external_stress_screen_not_locked_selection
- `EcoNiche-Opt-MAP4K1-TBX3AXL-PercentilePrimarySelected`: cohorts=GSE145996,PHS000452_LIU_LIKE_PRE; AUROC=0.656; AUPRC=0.599; BA=0.627; boundary=primary_selected_locked_candidate
- `EcoNiche-Opt-MAP4K1-TBX3-PercentileFixed`: cohorts=GSE145996,PHS000452_LIU_LIKE_PRE; AUROC=0.653; AUPRC=0.661; BA=0.595; boundary=literature_prior_fixed_pair
- `EcoNiche-Opt-MAP4K1-TBX3AXL-ZScoreStressOnly`: cohorts=GSE145996,PHS000452_LIU_LIKE_PRE; AUROC=0.643; AUPRC=0.592; BA=0.572; boundary=current_external_stress_screen_not_locked_selection
- `EcoNiche-Opt-MAP4K1-TBX3AXL-RobustZStressOnly`: cohorts=GSE145996,PHS000452_LIU_LIKE_PRE; AUROC=0.638; AUPRC=0.577; BA=0.631; boundary=current_external_stress_screen_not_locked_selection
- `EcoNiche-Opt-MAP4K1-TBX3AXL-PercentilePrimarySelected`: cohorts=GSE145996,PHS000452_LIU_LIKE_PRE; AUROC=0.636; AUPRC=0.538; BA=0.626; boundary=primary_selected_locked_candidate
- `EcoNiche-Opt-MAP4K1-TBX3-PercentileFixed`: cohorts=GSE145996,PHS000452_LIU_LIKE_PRE; AUROC=0.615; AUPRC=0.561; BA=0.573; boundary=literature_prior_fixed_pair

## expanded_public_melanoma

- `EcoNiche-Opt-MAP4K1-TBX3-PercentileFixed`: cohorts=GSE115821,GSE145996,GSE168204,PHS000452_LIU_LIKE_PRE; AUROC=0.660; AUPRC=0.626; BA=0.586; boundary=literature_prior_fixed_pair
- `EcoNiche-Opt-MAP4K1-TBX3AXL-PercentilePrimarySelected`: cohorts=GSE115821,GSE145996,GSE168204,PHS000452_LIU_LIKE_PRE; AUROC=0.644; AUPRC=0.579; BA=0.613; boundary=primary_selected_locked_candidate
- `EcoNiche-Opt-MAP4K1-TBX3AXL-RobustZStressOnly`: cohorts=GSE115821,GSE145996,GSE168204,PHS000452_LIU_LIKE_PRE; AUROC=0.634; AUPRC=0.602; BA=0.590; boundary=current_external_stress_screen_not_locked_selection
- `EcoNiche-Opt-MAP4K1-TBX3AXL-PercentilePrimarySelected`: cohorts=GSE115821,GSE145996,GSE168204,PHS000452_LIU_LIKE_PRE; AUROC=0.627; AUPRC=0.528; BA=0.608; boundary=primary_selected_locked_candidate
- `EcoNiche-Opt-MAP4K1-TBX3-PercentileFixed`: cohorts=GSE115821,GSE145996,GSE168204,PHS000452_LIU_LIKE_PRE; AUROC=0.626; AUPRC=0.550; BA=0.568; boundary=literature_prior_fixed_pair
- `EcoNiche-Opt-MAP4K1-TBX3AXL-ZScoreStressOnly`: cohorts=GSE115821,GSE145996,GSE168204,PHS000452_LIU_LIKE_PRE; AUROC=0.624; AUPRC=0.594; BA=0.565; boundary=current_external_stress_screen_not_locked_selection
- `EcoNiche-Opt-MAP4K1-TBX3AXL-RobustZStressOnly`: cohorts=GSE115821,GSE145996,GSE168204,PHS000452_LIU_LIKE_PRE; AUROC=0.607; AUPRC=0.532; BA=0.606; boundary=current_external_stress_screen_not_locked_selection
- `EcoNiche-Opt-MAP4K1-TBX3AXL-ZScoreStressOnly`: cohorts=GSE115821,GSE145996,GSE168204,PHS000452_LIU_LIKE_PRE; AUROC=0.605; AUPRC=0.541; BA=0.547; boundary=current_external_stress_screen_not_locked_selection

## low_n_array_sensitivity

- `EcoNiche-Opt-MAP4K1-TBX3-PercentileFixed`: cohorts=GSE122220; AUROC=0.875; AUPRC=0.804; BA=0.775; boundary=literature_prior_fixed_pair
- `EcoNiche-Opt-MAP4K1-TBX3AXL-PercentilePrimarySelected`: cohorts=GSE122220; AUROC=0.850; AUPRC=0.804; BA=0.800; boundary=primary_selected_locked_candidate
- `EcoNiche-Opt-MAP4K1-TBX3AXL-ZScoreStressOnly`: cohorts=GSE122220; AUROC=0.800; AUPRC=0.679; BA=0.800; boundary=current_external_stress_screen_not_locked_selection
- `EcoNiche-Opt-MAP4K1-TBX3AXL-RobustZStressOnly`: cohorts=GSE122220; AUROC=0.800; AUPRC=0.679; BA=0.700; boundary=current_external_stress_screen_not_locked_selection
- `EcoNiche-Opt-MAP4K1-TBX3-PercentileFixed`: cohorts=GSE122220; AUROC=0.792; AUPRC=0.608; BA=0.708; boundary=literature_prior_fixed_pair
- `EcoNiche-Opt-MAP4K1-TBX3AXL-PercentilePrimarySelected`: cohorts=GSE122220; AUROC=0.708; AUPRC=0.567; BA=0.750; boundary=primary_selected_locked_candidate
- `EcoNiche-Opt-MAP4K1-TBX3AXL-ZScoreStressOnly`: cohorts=GSE122220; AUROC=0.667; AUPRC=0.525; BA=0.750; boundary=current_external_stress_screen_not_locked_selection
- `EcoNiche-Opt-MAP4K1-TBX3AXL-RobustZStressOnly`: cohorts=GSE122220; AUROC=0.667; AUPRC=0.525; BA=0.667; boundary=current_external_stress_screen_not_locked_selection

## public_melanoma_nontraining_with_combo

- `EcoNiche-Opt-MAP4K1-TBX3AXL-RobustZStressOnly`: cohorts=GSE115821,GSE145996,GSE168204,PHS000452_LIU_LIKE_PRE,PRJEB23709_COMBO_PRE; AUROC=0.671; AUPRC=0.680; BA=0.612; boundary=current_external_stress_screen_not_locked_selection
- `EcoNiche-Opt-MAP4K1-TBX3-PercentileFixed`: cohorts=GSE115821,GSE145996,GSE168204,PHS000452_LIU_LIKE_PRE,PRJEB23709_COMBO_PRE; AUROC=0.667; AUPRC=0.673; BA=0.588; boundary=literature_prior_fixed_pair
- `EcoNiche-Opt-MAP4K1-TBX3AXL-ZScoreStressOnly`: cohorts=GSE115821,GSE145996,GSE168204,PHS000452_LIU_LIKE_PRE,PRJEB23709_COMBO_PRE; AUROC=0.658; AUPRC=0.672; BA=0.587; boundary=current_external_stress_screen_not_locked_selection
- `EcoNiche-Opt-MAP4K1-TBX3-PercentileFixed`: cohorts=GSE115821,GSE145996,GSE168204,PHS000452_LIU_LIKE_PRE,PRJEB23709_COMBO_PRE; AUROC=0.628; AUPRC=0.571; BA=0.564; boundary=literature_prior_fixed_pair
- `EcoNiche-Opt-MAP4K1-TBX3AXL-PercentilePrimarySelected`: cohorts=GSE115821,GSE145996,GSE168204,PHS000452_LIU_LIKE_PRE,PRJEB23709_COMBO_PRE; AUROC=0.626; AUPRC=0.619; BA=0.578; boundary=primary_selected_locked_candidate
- `EcoNiche-Opt-MAP4K1-TBX3AXL-ZScoreStressOnly`: cohorts=GSE115821,GSE145996,GSE168204,PHS000452_LIU_LIKE_PRE,PRJEB23709_COMBO_PRE; AUROC=0.603; AUPRC=0.557; BA=0.549; boundary=current_external_stress_screen_not_locked_selection
- `EcoNiche-Opt-MAP4K1-TBX3AXL-RobustZStressOnly`: cohorts=GSE115821,GSE145996,GSE168204,PHS000452_LIU_LIKE_PRE,PRJEB23709_COMBO_PRE; AUROC=0.600; AUPRC=0.546; BA=0.598; boundary=current_external_stress_screen_not_locked_selection
- `EcoNiche-Opt-MAP4K1-TBX3AXL-PercentilePrimarySelected`: cohorts=GSE115821,GSE145996,GSE168204,PHS000452_LIU_LIKE_PRE,PRJEB23709_COMBO_PRE; AUROC=0.588; AUPRC=0.526; BA=0.586; boundary=primary_selected_locked_candidate

## public_melanoma_nontraining_with_combo_and_array

- `EcoNiche-Opt-MAP4K1-TBX3-PercentileFixed`: cohorts=GSE115821,GSE122220,GSE145996,GSE168204,PHS000452_LIU_LIKE_PRE,PRJEB23709_COMBO_PRE; AUROC=0.683; AUPRC=0.681; BA=0.596; boundary=literature_prior_fixed_pair
- `EcoNiche-Opt-MAP4K1-TBX3AXL-RobustZStressOnly`: cohorts=GSE115821,GSE122220,GSE145996,GSE168204,PHS000452_LIU_LIKE_PRE,PRJEB23709_COMBO_PRE; AUROC=0.677; AUPRC=0.670; BA=0.616; boundary=current_external_stress_screen_not_locked_selection
- `EcoNiche-Opt-MAP4K1-TBX3AXL-ZScoreStressOnly`: cohorts=GSE115821,GSE122220,GSE145996,GSE168204,PHS000452_LIU_LIKE_PRE,PRJEB23709_COMBO_PRE; AUROC=0.669; AUPRC=0.665; BA=0.596; boundary=current_external_stress_screen_not_locked_selection
- `EcoNiche-Opt-MAP4K1-TBX3-PercentileFixed`: cohorts=GSE115821,GSE122220,GSE145996,GSE168204,PHS000452_LIU_LIKE_PRE,PRJEB23709_COMBO_PRE; AUROC=0.642; AUPRC=0.573; BA=0.569; boundary=literature_prior_fixed_pair
- `EcoNiche-Opt-MAP4K1-TBX3AXL-PercentilePrimarySelected`: cohorts=GSE115821,GSE122220,GSE145996,GSE168204,PHS000452_LIU_LIKE_PRE,PRJEB23709_COMBO_PRE; AUROC=0.641; AUPRC=0.627; BA=0.589; boundary=primary_selected_locked_candidate
- `EcoNiche-Opt-MAP4K1-TBX3AXL-ZScoreStressOnly`: cohorts=GSE115821,GSE122220,GSE145996,GSE168204,PHS000452_LIU_LIKE_PRE,PRJEB23709_COMBO_PRE; AUROC=0.613; AUPRC=0.550; BA=0.557; boundary=current_external_stress_screen_not_locked_selection
- `EcoNiche-Opt-MAP4K1-TBX3AXL-RobustZStressOnly`: cohorts=GSE115821,GSE122220,GSE145996,GSE168204,PHS000452_LIU_LIKE_PRE,PRJEB23709_COMBO_PRE; AUROC=0.609; AUPRC=0.542; BA=0.600; boundary=current_external_stress_screen_not_locked_selection
- `EcoNiche-Opt-MAP4K1-TBX3AXL-PercentilePrimarySelected`: cohorts=GSE115821,GSE122220,GSE145996,GSE168204,PHS000452_LIU_LIKE_PRE,PRJEB23709_COMBO_PRE; AUROC=0.599; AUPRC=0.523; BA=0.593; boundary=primary_selected_locked_candidate

## secondary_small_melanoma

- `EcoNiche-Opt-MAP4K1-TBX3-PercentileFixed`: cohorts=GSE115821,GSE168204; AUROC=0.657; AUPRC=0.554; BA=0.535; boundary=literature_prior_fixed_pair
- `EcoNiche-Opt-MAP4K1-TBX3-PercentileFixed`: cohorts=GSE115821,GSE168204; AUROC=0.657; AUPRC=0.554; BA=0.535; boundary=literature_prior_fixed_pair
- `EcoNiche-Opt-MAP4K1-TBX3AXL-RobustZStressOnly`: cohorts=GSE115821,GSE168204; AUROC=0.533; AUPRC=0.464; BA=0.475; boundary=current_external_stress_screen_not_locked_selection
- `EcoNiche-Opt-MAP4K1-TBX3AXL-RobustZStressOnly`: cohorts=GSE115821,GSE168204; AUROC=0.533; AUPRC=0.464; BA=0.503; boundary=current_external_stress_screen_not_locked_selection
- `EcoNiche-Opt-MAP4K1-TBX3AXL-PercentilePrimarySelected`: cohorts=GSE115821,GSE168204; AUROC=0.530; AUPRC=0.476; BA=0.540; boundary=primary_selected_locked_candidate
- `EcoNiche-Opt-MAP4K1-TBX3AXL-PercentilePrimarySelected`: cohorts=GSE115821,GSE168204; AUROC=0.530; AUPRC=0.476; BA=0.505; boundary=primary_selected_locked_candidate
- `EcoNiche-Opt-MAP4K1-TBX3AXL-ZScoreStressOnly`: cohorts=GSE115821,GSE168204; AUROC=0.513; AUPRC=0.473; BA=0.439; boundary=current_external_stress_screen_not_locked_selection
- `EcoNiche-Opt-MAP4K1-TBX3AXL-ZScoreStressOnly`: cohorts=GSE115821,GSE168204; AUROC=0.513; AUPRC=0.473; BA=0.439; boundary=current_external_stress_screen_not_locked_selection

## Interpretation

The low-n GSE122220 microarray sensitivity set can show favorable point estimates, but it is not a strict bulk RNA-seq external validation cohort. Strict-compatible public melanoma sets still do not close the AUROC >=0.70 target, so these results support sensitivity reporting and reinforce the need for newly obtained controlled independent melanoma tumor-tissue validation.
