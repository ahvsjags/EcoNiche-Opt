# Clinical Interpretability Audit

This audit is derived from registered prediction tables only. It does not retrain the model, choose features, alter thresholds, or use locked external labels for model selection.

- High-score enrichment contexts: 22.
- Subgroup metric rows: 92.
- Threshold operating-point rows: 88.
- Calibration-bin rows: 220.
- Primary high-evidence melanoma high-score response rate: 0.481 versus 0.270 low-score; Fisher two-sided q=0.238.
- Strict melanoma PD1-like external high-score response rate: 0.527 versus 0.443 low-score; Fisher two-sided q=0.67.

Files:
- high_score_enrichment.tsv
- subgroup_metrics.tsv
- threshold_operating_points.tsv
- calibration_bins.tsv