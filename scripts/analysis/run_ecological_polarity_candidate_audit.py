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


DISCOVERY_SETS: dict[str, list[str]] = {
    "high_evidence_primary": ["GSE91061", "GSE78220", "PRJEB23709_PD1_PRE"],
    "expanded_primary_with_mgh": ["GSE91061", "GSE78220", "PRJEB23709_PD1_PRE", "GSE115821", "GSE168204"],
}
STRICT_EXTERNAL_COHORTS = ["GSE145996", "PHS000452_LIU_LIKE_PRE"]
PRIMARY_ENDPOINT = "primary_recist"
STRICT_ENDPOINT = "strict_recist"
PRIMARY_STRATUM = "melanoma_core_high_evidence"
STRICT_STRATUM = "strict_melanoma_pd1_like_external"

CANDIDATE_AXES: dict[str, tuple[list[str], list[str]]] = {
    "MAP4K1_minus_TBX3": (["MAP4K1"], ["TBX3"]),
    "MAP4K1_minus_TBX3_AXL": (["MAP4K1"], ["TBX3", "AXL"]),
    "MAP4K1_minus_TBX3_IPRES": (
        ["MAP4K1"],
        ["TBX3", "AXL", "VIM", "ZEB1", "FN1", "VEGFA", "WNT5A", "MITF", "SOX10"],
    ),
    "MAP4K1_IFN_minus_TBX3_AXL": (
        ["MAP4K1", "IFNG", "CXCL9", "CXCL10", "STAT1", "IRF1", "GBP1"],
        ["TBX3", "AXL"],
    ),
    "MAP4K1_APM_minus_TBX3_AXL": (
        ["MAP4K1", "HLA-A", "B2M", "TAP1", "TAP2"],
        ["TBX3", "AXL"],
    ),
    "MAP4K1_EFF_minus_TBX3_AXL": (
        ["MAP4K1", "CD8A", "CD8B", "GZMB", "PRF1", "NKG7"],
        ["TBX3", "AXL"],
    ),
    "IMMUNE_ALL_minus_TBX3_AXL": (
        [
            "MAP4K1",
            "IFNG",
            "CXCL9",
            "CXCL10",
            "STAT1",
            "IRF1",
            "GBP1",
            "CD8A",
            "CD8B",
            "GZMB",
            "PRF1",
            "NKG7",
            "HLA-A",
            "B2M",
            "TAP1",
            "TAP2",
        ],
        ["TBX3", "AXL"],
    ),
    "IMMUNE_ALL_minus_IPRES_STROMAL": (
        [
            "MAP4K1",
            "IFNG",
            "CXCL9",
            "CXCL10",
            "STAT1",
            "IRF1",
            "GBP1",
            "CD8A",
            "CD8B",
            "GZMB",
            "PRF1",
            "NKG7",
            "HLA-A",
            "B2M",
            "TAP1",
            "TAP2",
        ],
        [
            "TBX3",
            "AXL",
            "VIM",
            "ZEB1",
            "FN1",
            "VEGFA",
            "WNT5A",
            "MITF",
            "SOX10",
            "COL1A1",
            "TGFBI",
            "TGFB1",
            "FAP",
            "ACTA2",
            "MMP2",
        ],
    ),
}
TRANSFORMS = ["cohort_gene_percentile", "cohort_zscore", "raw_log"]
NEGATIVE_WEIGHTS = [0.75, 1.0, 1.25]


def _all_candidate_genes() -> list[str]:
    return sorted({gene for pair in CANDIDATE_AXES.values() for genes in pair for gene in genes})


def _precompute_transforms(X_by_cohort: dict[str, pd.DataFrame]) -> dict[str, dict[str, pd.DataFrame]]:
    genes = _all_candidate_genes()
    transforms: dict[str, dict[str, pd.DataFrame]] = {method: {} for method in TRANSFORMS}
    for cohort, X in X_by_cohort.items():
        available = [gene for gene in genes if gene in X.columns]
        values = X[available].apply(pd.to_numeric, errors="coerce") if available else pd.DataFrame(index=X.index)
        transforms["cohort_gene_percentile"][cohort] = values.rank(axis=0, pct=True).fillna(0.5)
        transforms["cohort_zscore"][cohort] = ((values - values.mean()) / (values.std(ddof=0) + 1e-6)).fillna(0.0)
        transforms["raw_log"][cohort] = np.log2(values.clip(lower=1e-3)).fillna(0.0)
    return transforms


def score_axis(
    transformed: pd.DataFrame,
    positive_genes: list[str],
    negative_genes: list[str],
    negative_weight: float,
) -> tuple[pd.Series, dict[str, object]]:
    positives = [gene for gene in positive_genes if gene in transformed.columns]
    negatives = [gene for gene in negative_genes if gene in transformed.columns]
    score = pd.Series(0.0, index=transformed.index)
    if positives:
        score = score + transformed[positives].mean(axis=1)
    if negatives:
        score = score - negative_weight * transformed[negatives].mean(axis=1)
    if score.max() > score.min():
        score = (score - score.min()) / (score.max() - score.min())
    else:
        score = pd.Series(0.5, index=transformed.index)
    return score.astype(float), {
        "n_positive_genes_available": len(positives),
        "n_negative_genes_available": len(negatives),
        "positive_genes_available": ",".join(positives),
        "negative_genes_available": ",".join(negatives),
    }


def _metric_row(y: pd.Series, p: pd.Series, thresholds: pd.Series | float) -> dict[str, object]:
    if isinstance(thresholds, pd.Series):
        pred = (p.astype(float) >= thresholds.astype(float)).astype(int)
        threshold = float(thresholds.median())
        threshold_min = float(thresholds.min())
        threshold_max = float(thresholds.max())
    else:
        pred = (p.astype(float) >= float(thresholds)).astype(int)
        threshold = threshold_min = threshold_max = float(thresholds)
    metrics = compute_binary_metrics(y.astype(int), p.astype(float), threshold=threshold)
    metrics["balanced_accuracy"] = float(balanced_accuracy_score(y.astype(int), pred))
    metrics["threshold_median"] = threshold
    metrics["threshold_min"] = threshold_min
    metrics["threshold_max"] = threshold_max
    metrics["n_samples"] = int(len(y))
    metrics["n_responders"] = int(y.astype(int).sum())
    metrics["n_nonresponders"] = int((y.astype(int) == 0).sum())
    metrics["response_prevalence"] = float(y.astype(int).mean())
    metrics["AUPRC_minus_prevalence"] = float(metrics["AUPRC"] - metrics["response_prevalence"])
    metrics["selection_score"] = float(
        metrics["AUROC"] + 0.25 * metrics["AUPRC"] + 0.10 * metrics["balanced_accuracy"] - 0.10 * metrics["ECE"]
    )
    return metrics


def _candidate_scores(
    transforms: dict[str, dict[str, pd.DataFrame]],
    cohorts: list[str],
    axis: str,
    method: str,
    negative_weight: float,
) -> tuple[dict[str, pd.Series], pd.DataFrame]:
    positive, negative = CANDIDATE_AXES[axis]
    scores: dict[str, pd.Series] = {}
    coverage_rows = []
    for cohort in cohorts:
        if cohort not in transforms[method]:
            continue
        score, coverage = score_axis(transforms[method][cohort], positive, negative, negative_weight)
        scores[cohort] = score
        coverage_rows.append(
            {
                "cohort": cohort,
                "axis": axis,
                "transform": method,
                "negative_weight": negative_weight,
                "n_samples": len(score),
                **coverage,
            }
        )
    return scores, pd.DataFrame(coverage_rows)


def _evaluate_primary_candidate(
    endpoint_data,
    transforms: dict[str, dict[str, pd.DataFrame]],
    discovery_cohorts: list[str],
    discovery_set: str,
    axis: str,
    method: str,
    negative_weight: float,
) -> tuple[dict[str, object] | None, list[dict[str, object]]]:
    scores, _ = _candidate_scores(transforms, discovery_cohorts, axis, method, negative_weight)
    parts: list[dict[str, object]] = []
    for holdout in discovery_cohorts:
        if holdout not in scores or holdout not in endpoint_data.y_response_by_cohort:
            continue
        train = [cohort for cohort in discovery_cohorts if cohort != holdout and cohort in scores]
        if len(train) < 2:
            continue
        y_train = _concat(endpoint_data.y_response_by_cohort, train).astype(int)
        y_test = endpoint_data.y_response_by_cohort[holdout].astype(int)
        if y_train.nunique() < 2 or y_test.nunique() < 2:
            continue
        p_train = _concat(scores, train).astype(float)
        threshold = select_threshold(y_train.to_numpy(dtype=int), p_train.to_numpy(dtype=float))
        p_test = scores[holdout].reindex(y_test.index).astype(float)
        for sample_id, prob in p_test.items():
            parts.append(
                {
                    "endpoint": PRIMARY_ENDPOINT,
                    "stratum": PRIMARY_STRATUM,
                    "discovery_set": discovery_set,
                    "cohort": holdout,
                    "sample_id": sample_id,
                    "true_response_label": int(y_test.loc[sample_id]),
                    "response_probability": float(prob),
                    "threshold": float(threshold),
                    "axis": axis,
                    "transform": method,
                    "negative_weight": negative_weight,
                    "selection_boundary": "outer_lodo_training_threshold_only",
                }
            )
    if not parts:
        return None, []
    predictions = pd.DataFrame(parts)
    metrics = _metric_row(
        predictions["true_response_label"].astype(int),
        predictions["response_probability"].astype(float),
        predictions["threshold"].astype(float),
    )
    return {
        "endpoint": PRIMARY_ENDPOINT,
        "stratum": PRIMARY_STRATUM,
        "discovery_set": discovery_set,
        "axis": axis,
        "transform": method,
        "negative_weight": negative_weight,
        "selection_boundary": "primary_lodo_only_no_external_selection",
        **metrics,
    }, parts


def _evaluate_strict_external_candidate(
    endpoint_data,
    transforms: dict[str, dict[str, pd.DataFrame]],
    discovery_cohorts: list[str],
    discovery_set: str,
    axis: str,
    method: str,
    negative_weight: float,
) -> tuple[dict[str, object] | None, list[dict[str, object]]]:
    cohorts = [*discovery_cohorts, *STRICT_EXTERNAL_COHORTS]
    scores, _ = _candidate_scores(transforms, cohorts, axis, method, negative_weight)
    train = [cohort for cohort in discovery_cohorts if cohort in scores and cohort in endpoint_data.y_response_by_cohort]
    if len(train) < 2:
        return None, []
    y_train = _concat(endpoint_data.y_response_by_cohort, train).astype(int)
    p_train = _concat(scores, train).astype(float)
    threshold = select_threshold(y_train.to_numpy(dtype=int), p_train.to_numpy(dtype=float))
    parts: list[dict[str, object]] = []
    for cohort in STRICT_EXTERNAL_COHORTS:
        if cohort not in scores or cohort not in endpoint_data.y_response_by_cohort:
            continue
        y_test = endpoint_data.y_response_by_cohort[cohort].astype(int)
        if y_test.nunique() < 2:
            continue
        p_test = scores[cohort].reindex(y_test.index).astype(float)
        for sample_id, prob in p_test.items():
            parts.append(
                {
                    "endpoint": STRICT_ENDPOINT,
                    "stratum": STRICT_STRATUM,
                    "discovery_set": discovery_set,
                    "cohort": cohort,
                    "sample_id": sample_id,
                    "true_response_label": int(y_test.loc[sample_id]),
                    "response_probability": float(prob),
                    "threshold": float(threshold),
                    "axis": axis,
                    "transform": method,
                    "negative_weight": negative_weight,
                    "selection_boundary": "discovery_only_threshold_locked_external_scoring",
                }
            )
    if not parts:
        return None, []
    predictions = pd.DataFrame(parts)
    metrics = _metric_row(
        predictions["true_response_label"].astype(int),
        predictions["response_probability"].astype(float),
        float(threshold),
    )
    return {
        "endpoint": STRICT_ENDPOINT,
        "stratum": STRICT_STRATUM,
        "discovery_set": discovery_set,
        "axis": axis,
        "transform": method,
        "negative_weight": negative_weight,
        "selection_boundary": "discovery_only_no_external_selection",
        **metrics,
    }, parts


def run_audit(
    X_by_cohort: dict[str, pd.DataFrame],
    metadata_by_cohort: dict[str, pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    primary_rows: list[dict[str, object]] = []
    external_rows: list[dict[str, object]] = []
    primary_predictions: list[dict[str, object]] = []
    external_predictions: list[dict[str, object]] = []
    selection_rows: list[dict[str, object]] = []

    for discovery_set, discovery_cohorts in DISCOVERY_SETS.items():
        primary_data = prepare_endpoint_data(X_by_cohort, metadata_by_cohort, discovery_cohorts, PRIMARY_ENDPOINT)
        strict_data = prepare_endpoint_data(
            X_by_cohort,
            metadata_by_cohort,
            [*discovery_cohorts, *STRICT_EXTERNAL_COHORTS],
            STRICT_ENDPOINT,
        )
        primary_transforms = _precompute_transforms(primary_data.X_by_cohort)
        strict_transforms = _precompute_transforms(strict_data.X_by_cohort)
        for axis in CANDIDATE_AXES:
            for method in TRANSFORMS:
                for negative_weight in NEGATIVE_WEIGHTS:
                    primary_row, primary_parts = _evaluate_primary_candidate(
                        primary_data,
                        primary_transforms,
                        discovery_cohorts,
                        discovery_set,
                        axis,
                        method,
                        negative_weight,
                    )
                    external_row, external_parts = _evaluate_strict_external_candidate(
                        strict_data,
                        strict_transforms,
                        discovery_cohorts,
                        discovery_set,
                        axis,
                        method,
                        negative_weight,
                    )
                    if primary_row is not None:
                        primary_rows.append(primary_row)
                        primary_predictions.extend(primary_parts)
                    if external_row is not None:
                        external_rows.append(external_row)
                        external_predictions.extend(external_parts)
        primary_frame = pd.DataFrame([row for row in primary_rows if row["discovery_set"] == discovery_set])
        external_frame = pd.DataFrame([row for row in external_rows if row["discovery_set"] == discovery_set])
        if not primary_frame.empty:
            selected = primary_frame.sort_values(
                ["selection_score", "AUROC", "AUPRC", "balanced_accuracy"],
                ascending=False,
            ).iloc[0]
            external_match = external_frame[
                external_frame["axis"].astype(str).eq(str(selected["axis"]))
                & external_frame["transform"].astype(str).eq(str(selected["transform"]))
                & np.isclose(external_frame["negative_weight"].astype(float), float(selected["negative_weight"]))
            ]
            stress = external_frame.sort_values(["AUROC", "AUPRC", "balanced_accuracy"], ascending=False).iloc[0]
            selection_rows.append(
                {
                    "discovery_set": discovery_set,
                    "selection_id": "primary_selected_candidate",
                    "axis": selected["axis"],
                    "transform": selected["transform"],
                    "negative_weight": float(selected["negative_weight"]),
                    "claim_boundary": "selected_by_primary_lodo_only_not_by_external",
                    "primary_AUROC": float(selected["AUROC"]),
                    "primary_AUPRC": float(selected["AUPRC"]),
                    "primary_balanced_accuracy": float(selected["balanced_accuracy"]),
                    "primary_ECE": float(selected["ECE"]),
                    "strict_external_AUROC": float(external_match.iloc[0]["AUROC"]) if not external_match.empty else np.nan,
                    "strict_external_AUPRC": float(external_match.iloc[0]["AUPRC"]) if not external_match.empty else np.nan,
                    "strict_external_balanced_accuracy": float(external_match.iloc[0]["balanced_accuracy"])
                    if not external_match.empty
                    else np.nan,
                    "strict_external_ECE": float(external_match.iloc[0]["ECE"]) if not external_match.empty else np.nan,
                }
            )
            selection_rows.append(
                {
                    "discovery_set": discovery_set,
                    "selection_id": "current_external_stress_best",
                    "axis": stress["axis"],
                    "transform": stress["transform"],
                    "negative_weight": float(stress["negative_weight"]),
                    "claim_boundary": "current_external_stress_screen_not_a_locked_selection_claim",
                    "primary_AUROC": float(
                        primary_frame[
                            primary_frame["axis"].astype(str).eq(str(stress["axis"]))
                            & primary_frame["transform"].astype(str).eq(str(stress["transform"]))
                            & np.isclose(primary_frame["negative_weight"].astype(float), float(stress["negative_weight"]))
                        ].iloc[0]["AUROC"]
                    ),
                    "primary_AUPRC": float(
                        primary_frame[
                            primary_frame["axis"].astype(str).eq(str(stress["axis"]))
                            & primary_frame["transform"].astype(str).eq(str(stress["transform"]))
                            & np.isclose(primary_frame["negative_weight"].astype(float), float(stress["negative_weight"]))
                        ].iloc[0]["AUPRC"]
                    ),
                    "primary_balanced_accuracy": float(
                        primary_frame[
                            primary_frame["axis"].astype(str).eq(str(stress["axis"]))
                            & primary_frame["transform"].astype(str).eq(str(stress["transform"]))
                            & np.isclose(primary_frame["negative_weight"].astype(float), float(stress["negative_weight"]))
                        ].iloc[0]["balanced_accuracy"]
                    ),
                    "primary_ECE": float(
                        primary_frame[
                            primary_frame["axis"].astype(str).eq(str(stress["axis"]))
                            & primary_frame["transform"].astype(str).eq(str(stress["transform"]))
                            & np.isclose(primary_frame["negative_weight"].astype(float), float(stress["negative_weight"]))
                        ].iloc[0]["ECE"]
                    ),
                    "strict_external_AUROC": float(stress["AUROC"]),
                    "strict_external_AUPRC": float(stress["AUPRC"]),
                    "strict_external_balanced_accuracy": float(stress["balanced_accuracy"]),
                    "strict_external_ECE": float(stress["ECE"]),
                }
            )

    return (
        pd.DataFrame(primary_rows),
        pd.DataFrame(external_rows),
        pd.DataFrame(primary_predictions),
        pd.DataFrame(external_predictions),
        pd.DataFrame(selection_rows),
    )


def write_markdown(selection: pd.DataFrame, out_md: Path) -> None:
    lines = [
        "# Ecological Polarity Candidate Audit",
        "",
        "This audit tests a predeclared family of immune-effector, antigen-presentation, dedifferentiation, IPRES, and stromal polarity scores.",
        "Candidate choice is made from primary melanoma LODO evidence only. Current strict external stress rows are diagnostic and cannot define the locked model.",
        "",
    ]
    if selection.empty:
        lines.append("No candidate-selection rows were produced.")
    else:
        for _, row in selection.iterrows():
            lines.append(
                "- `{}` `{}` `{}` weight={:.2f}: primary AUROC={:.3f}; strict external AUROC={:.3f}; boundary={}".format(
                    row["discovery_set"],
                    row["selection_id"],
                    row["axis"],
                    float(row["negative_weight"]),
                    float(row["primary_AUROC"]),
                    float(row["strict_external_AUROC"]),
                    row["claim_boundary"],
                )
            )
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "A candidate can support the strict external target only if the primary-selected row reaches AUROC >=0.70 on strict external scoring. Otherwise it remains a negative optimization audit.",
            "",
        ]
    )
    out_md.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--out", default="results/ecological_polarity_candidate_audit_20260527")
    args = parser.parse_args()

    out = ROOT / args.out
    out.mkdir(parents=True, exist_ok=True)
    X_by_cohort, _, metadata_by_cohort = load_processed_bulk(ROOT / args.processed_dir)
    primary, external, primary_predictions, external_predictions, selection = run_audit(X_by_cohort, metadata_by_cohort)
    primary.to_csv(out / "ecological_polarity_primary_summary.tsv", sep="\t", index=False)
    external.to_csv(out / "ecological_polarity_strict_external_summary.tsv", sep="\t", index=False)
    primary_predictions.to_csv(out / "ecological_polarity_primary_predictions.tsv", sep="\t", index=False)
    external_predictions.to_csv(out / "ecological_polarity_strict_external_predictions.tsv", sep="\t", index=False)
    selection.to_csv(out / "ecological_polarity_selection.tsv", sep="\t", index=False)
    write_markdown(selection, out / "ECOLOGICAL_POLARITY_CANDIDATE_AUDIT.md")
    best = selection[selection["selection_id"].astype(str).eq("primary_selected_candidate")]
    print(json.dumps(best.to_dict("records"), ensure_ascii=False))
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
