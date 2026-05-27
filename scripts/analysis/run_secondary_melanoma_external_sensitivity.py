from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import balanced_accuracy_score

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.io import load_processed_bulk
from econiche.metrics import compute_binary_metrics
from econiche_opt.model.endpoint_modules import _concat, prepare_endpoint_data, select_threshold


DISCOVERY_COHORTS = ["GSE91061", "GSE78220", "PRJEB23709_PD1_PRE"]
EXTERNAL_SETS: dict[str, list[str]] = {
    "current_strict_pd1_like": ["GSE145996", "PHS000452_LIU_LIKE_PRE"],
    "secondary_small_melanoma": ["GSE115821", "GSE168204"],
    "low_n_array_sensitivity": ["GSE122220"],
    "expanded_public_melanoma": ["GSE115821", "GSE168204", "GSE145996", "PHS000452_LIU_LIKE_PRE"],
    "public_melanoma_nontraining_with_combo": [
        "GSE115821",
        "GSE168204",
        "GSE145996",
        "PHS000452_LIU_LIKE_PRE",
        "PRJEB23709_COMBO_PRE",
    ],
    "public_melanoma_nontraining_with_combo_and_array": [
        "GSE115821",
        "GSE168204",
        "GSE145996",
        "PHS000452_LIU_LIKE_PRE",
        "PRJEB23709_COMBO_PRE",
        "GSE122220",
    ],
}
ENDPOINTS = ["primary_recist", "strict_recist"]

CANDIDATES: list[dict[str, object]] = [
    {
        "model_name": "EcoNiche-Opt-MAP4K1-TBX3AXL-PercentilePrimarySelected",
        "axis": "MAP4K1_minus_TBX3_AXL",
        "transform": "cohort_gene_percentile",
        "positive_genes": ["MAP4K1"],
        "negative_genes": ["TBX3", "AXL"],
        "negative_weight": 1.25,
        "claim_boundary": "primary_selected_locked_candidate",
    },
    {
        "model_name": "EcoNiche-Opt-MAP4K1-TBX3AXL-RobustZStressOnly",
        "axis": "MAP4K1_minus_TBX3_AXL",
        "transform": "cohort_robust_zscore",
        "positive_genes": ["MAP4K1"],
        "negative_genes": ["TBX3", "AXL"],
        "negative_weight": 1.0,
        "claim_boundary": "current_external_stress_screen_not_locked_selection",
    },
    {
        "model_name": "EcoNiche-Opt-MAP4K1-TBX3AXL-ZScoreStressOnly",
        "axis": "MAP4K1_minus_TBX3_AXL",
        "transform": "cohort_zscore",
        "positive_genes": ["MAP4K1"],
        "negative_genes": ["TBX3", "AXL"],
        "negative_weight": 1.0,
        "claim_boundary": "current_external_stress_screen_not_locked_selection",
    },
    {
        "model_name": "EcoNiche-Opt-MAP4K1-TBX3-PercentileFixed",
        "axis": "MAP4K1_minus_TBX3",
        "transform": "cohort_gene_percentile",
        "positive_genes": ["MAP4K1"],
        "negative_genes": ["TBX3"],
        "negative_weight": 1.0,
        "claim_boundary": "literature_prior_fixed_pair",
    },
]


def _all_genes() -> list[str]:
    return sorted(
        {
            gene
            for spec in CANDIDATES
            for gene in [*spec["positive_genes"], *spec["negative_genes"]]
        }
    )


def _transform_values(X: pd.DataFrame, method: str, genes: list[str]) -> pd.DataFrame:
    values = X[[gene for gene in genes if gene in X.columns]].apply(pd.to_numeric, errors="coerce")
    if method == "cohort_gene_percentile":
        return values.rank(axis=0, pct=True).fillna(0.5)
    if method == "cohort_zscore":
        return ((values - values.mean()) / (values.std(ddof=0) + 1e-6)).fillna(0.0)
    if method == "cohort_robust_zscore":
        median = values.median()
        mad = (values - median).abs().median() + 1e-6
        return ((values - median) / (1.4826 * mad)).clip(-5.0, 5.0).fillna(0.0)
    raise ValueError(f"Unsupported transform: {method}")


def _candidate_score(X: pd.DataFrame, spec: dict[str, object]) -> pd.Series:
    transformed = _transform_values(X, str(spec["transform"]), _all_genes())
    positives = [gene for gene in spec["positive_genes"] if gene in transformed.columns]
    negatives = [gene for gene in spec["negative_genes"] if gene in transformed.columns]
    score = pd.Series(0.0, index=X.index)
    if positives:
        score = score + transformed[positives].mean(axis=1)
    if negatives:
        score = score - float(spec["negative_weight"]) * transformed[negatives].mean(axis=1)
    if score.max() > score.min():
        score = (score - score.min()) / (score.max() - score.min())
    else:
        score = pd.Series(0.5, index=X.index)
    return score.astype(float)


def _metric_row(y: pd.Series, p: pd.Series, threshold: float) -> dict[str, object]:
    metrics = compute_binary_metrics(y.astype(int), p.astype(float), threshold=threshold)
    pred = (p.astype(float) >= threshold).astype(int)
    metrics["balanced_accuracy"] = float(balanced_accuracy_score(y.astype(int), pred))
    metrics["threshold"] = float(threshold)
    metrics["n_samples"] = int(len(y))
    metrics["n_responders"] = int(y.astype(int).sum())
    metrics["n_nonresponders"] = int((y.astype(int) == 0).sum())
    metrics["response_prevalence"] = float(y.astype(int).mean())
    metrics["AUPRC_minus_prevalence"] = float(metrics["AUPRC"] - metrics["response_prevalence"])
    return metrics


def run_sensitivity(
    X_by_cohort: dict[str, pd.DataFrame],
    metadata_by_cohort: dict[str, pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    all_cohorts = sorted({cohort for cohorts in EXTERNAL_SETS.values() for cohort in cohorts} | set(DISCOVERY_COHORTS))
    metric_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    for endpoint in ENDPOINTS:
        endpoint_data = prepare_endpoint_data(X_by_cohort, metadata_by_cohort, all_cohorts, endpoint)
        for spec in CANDIDATES:
            scores = {
                cohort: _candidate_score(X, spec)
                for cohort, X in endpoint_data.X_by_cohort.items()
            }
            train = [cohort for cohort in DISCOVERY_COHORTS if cohort in scores and cohort in endpoint_data.y_response_by_cohort]
            if not train:
                continue
            y_train = _concat(endpoint_data.y_response_by_cohort, train).astype(int)
            p_train = _concat(scores, train).astype(float)
            threshold = select_threshold(y_train.to_numpy(dtype=int), p_train.to_numpy(dtype=float))
            for external_set, cohorts in EXTERNAL_SETS.items():
                parts = []
                for cohort in cohorts:
                    if cohort not in endpoint_data.y_response_by_cohort or cohort not in scores:
                        continue
                    y = endpoint_data.y_response_by_cohort[cohort].astype(int)
                    if y.nunique() < 2:
                        continue
                    p = scores[cohort].reindex(y.index).astype(float)
                    for sample_id, prob in p.items():
                        parts.append(
                            {
                                "endpoint": endpoint,
                                "external_set": external_set,
                                "cohort": cohort,
                                "sample_id": sample_id,
                                "true_response_label": int(y.loc[sample_id]),
                                "response_probability": float(prob),
                                "threshold": float(threshold),
                                "model_name": spec["model_name"],
                                "axis": spec["axis"],
                                "transform": spec["transform"],
                                "negative_weight": spec["negative_weight"],
                                "claim_boundary": spec["claim_boundary"],
                                "selection_boundary": "discovery_only_threshold_locked_external_scoring",
                            }
                        )
                if not parts:
                    continue
                pred = pd.DataFrame(parts)
                metric_rows.append(
                    {
                        "endpoint": endpoint,
                        "external_set": external_set,
                        "model_name": spec["model_name"],
                        "axis": spec["axis"],
                        "transform": spec["transform"],
                        "negative_weight": spec["negative_weight"],
                        "claim_boundary": spec["claim_boundary"],
                        "cohorts": ",".join(sorted(pred["cohort"].unique())),
                        **_metric_row(
                            pred["true_response_label"].astype(int),
                            pred["response_probability"].astype(float),
                            float(threshold),
                        ),
                    }
                )
                prediction_rows.extend(parts)
    return pd.DataFrame(metric_rows), pd.DataFrame(prediction_rows)


def write_markdown(metrics: pd.DataFrame, out_md: Path) -> None:
    lines = [
        "# Secondary Melanoma External Sensitivity",
        "",
        "This audit scores fixed melanoma rescue-head candidates on secondary public melanoma cohorts without refitting, recalibration, feature selection, or threshold selection on external labels.",
        "",
    ]
    for external_set, frame in metrics.groupby("external_set"):
        lines.extend([f"## {external_set}", ""])
        for _, row in frame.sort_values("AUROC", ascending=False).iterrows():
            lines.append(
                "- `{}`: cohorts={}; AUROC={:.3f}; AUPRC={:.3f}; BA={:.3f}; boundary={}".format(
                    row["model_name"],
                    row["cohorts"],
                    float(row["AUROC"]),
                    float(row["AUPRC"]),
                    float(row["balanced_accuracy"]),
                    row["claim_boundary"],
                )
            )
        lines.append("")
    lines.extend(
        [
            "## Interpretation",
            "",
            "The low-n GSE122220 microarray sensitivity set can show favorable point estimates, but it is not a strict bulk RNA-seq external validation cohort. Strict-compatible public melanoma sets still do not close the AUROC >=0.70 target, so these results support sensitivity reporting and reinforce the need for newly obtained controlled independent melanoma tumor-tissue validation.",
            "",
        ]
    )
    out_md.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--out", default="results/secondary_melanoma_external_sensitivity_20260527")
    args = parser.parse_args()

    out = ROOT / args.out
    out.mkdir(parents=True, exist_ok=True)
    X_by_cohort, _, metadata_by_cohort = load_processed_bulk(ROOT / args.processed_dir)
    metrics, predictions = run_sensitivity(X_by_cohort, metadata_by_cohort)
    metrics.to_csv(out / "secondary_melanoma_external_metrics.tsv", sep="\t", index=False)
    predictions.to_csv(out / "secondary_melanoma_external_predictions.tsv", sep="\t", index=False)
    write_markdown(metrics, out / "SECONDARY_MELANOMA_EXTERNAL_SENSITIVITY.md")
    best = metrics.sort_values("AUROC", ascending=False).head(5)
    print(json.dumps(best[["external_set", "model_name", "AUROC", "AUPRC", "balanced_accuracy"]].to_dict("records"), ensure_ascii=False))
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
