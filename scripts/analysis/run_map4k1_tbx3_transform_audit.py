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


PRIMARY_COHORTS = ["GSE91061", "GSE78220", "PRJEB23709_PD1_PRE"]
STRICT_EXTERNAL_COHORTS = ["GSE145996", "PHS000452_LIU_LIKE_PRE"]

AXES: dict[str, tuple[list[str], list[str], str]] = {
    "MAP4K1_minus_TBX3": (["MAP4K1"], ["TBX3"], "literature_pair_tumor_immune_balance"),
    "MAP4K1_minus_TBX3_AXL": (["MAP4K1"], ["TBX3", "AXL"], "literature_pair_plus_ipres_dedifferentiation_axis"),
    "MAP4K1_minus_TBX3_PAK4": (["MAP4K1"], ["TBX3", "PAK4"], "literature_pair_plus_wnt_pak4_resistance_axis"),
    "MAP4K1_CXCL9_minus_TBX3": (["MAP4K1", "CXCL9"], ["TBX3"], "tcell_ifn_extension_of_pair"),
    "MAP4K1_APM_minus_TBX3": (["MAP4K1", "HLA-A", "B2M", "TAP1"], ["TBX3"], "antigen_presentation_extension_of_pair"),
}

TRANSFORMS = ["cohort_gene_percentile", "cohort_zscore", "cohort_robust_zscore", "sample_gene_rank", "raw_log_ratio"]


def _axis_genes() -> list[str]:
    return sorted({gene for pos, neg, _ in AXES.values() for gene in [*pos, *neg]})


def transform_axis_score(X: pd.DataFrame, positive: list[str], negative: list[str], method: str) -> tuple[pd.Series, dict[str, object]]:
    genes = [gene for gene in _axis_genes() if gene in X.columns]
    values = X[genes].apply(pd.to_numeric, errors="coerce") if genes else pd.DataFrame(index=X.index)
    pos = [gene for gene in positive if gene in values.columns]
    neg = [gene for gene in negative if gene in values.columns]
    if method == "cohort_gene_percentile":
        transformed = values.rank(axis=0, pct=True).fillna(0.5)
    elif method == "cohort_zscore":
        transformed = ((values - values.mean()) / (values.std(ddof=0) + 1e-6)).fillna(0.0)
    elif method == "cohort_robust_zscore":
        median = values.median()
        mad = (values - median).abs().median() + 1e-6
        transformed = ((values - median) / (1.4826 * mad)).clip(-5.0, 5.0).fillna(0.0)
    elif method == "sample_gene_rank":
        transformed = values.rank(axis=1, pct=True).fillna(0.5)
    elif method == "raw_log_ratio":
        transformed = np.log2(values.clip(lower=1e-3)).fillna(0.0)
    else:
        raise ValueError(f"Unsupported transform: {method}")
    score = pd.Series(0.0, index=X.index)
    if pos:
        score = score + transformed[pos].mean(axis=1)
    if neg:
        score = score - transformed[neg].mean(axis=1)
    if score.max() > score.min():
        score = (score - score.min()) / (score.max() - score.min())
    else:
        score = pd.Series(0.5, index=X.index)
    return score.astype(float), {
        "n_positive_genes_available": len(pos),
        "n_negative_genes_available": len(neg),
        "positive_genes_available": ",".join(pos),
        "negative_genes_available": ",".join(neg),
    }


def build_scores(X_by_cohort: dict[str, pd.DataFrame]) -> tuple[dict[str, dict[str, dict[str, pd.Series]]], pd.DataFrame]:
    scores: dict[str, dict[str, dict[str, pd.Series]]] = {}
    coverage_rows = []
    for method in TRANSFORMS:
        scores[method] = {}
        for axis, (pos, neg, rationale) in AXES.items():
            scores[method][axis] = {}
            for cohort, X in X_by_cohort.items():
                score, coverage = transform_axis_score(X, pos, neg, method)
                scores[method][axis][cohort] = score
                coverage_rows.append(
                    {
                        "method": method,
                        "axis": axis,
                        "cohort": cohort,
                        "rationale": rationale,
                        "n_samples": len(score),
                        **coverage,
                    }
                )
    return scores, pd.DataFrame(coverage_rows)


def _metric_row(y: pd.Series, p: pd.Series, thresholds: pd.Series | float) -> dict[str, object]:
    if isinstance(thresholds, pd.Series):
        pred = (p.astype(float) >= thresholds.astype(float)).astype(int)
        threshold = float(thresholds.median())
    else:
        pred = (p.astype(float) >= float(thresholds)).astype(int)
        threshold = float(thresholds)
    metrics = compute_binary_metrics(y.astype(int), p.astype(float), threshold=threshold)
    metrics["balanced_accuracy"] = float(balanced_accuracy_score(y.astype(int), pred))
    metrics["threshold"] = threshold
    metrics["n_samples"] = int(len(y))
    metrics["n_responders"] = int(y.astype(int).sum())
    metrics["n_nonresponders"] = int((y.astype(int) == 0).sum())
    metrics["response_prevalence"] = float(y.astype(int).mean())
    metrics["AUPRC_minus_prevalence"] = float(metrics["AUPRC"] - metrics["response_prevalence"])
    return metrics


def evaluate_primary(scores: dict[str, dict[str, dict[str, pd.Series]]], y_by: dict[str, pd.Series]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    pred_rows = []
    for method in TRANSFORMS:
        for axis in AXES:
            parts = []
            for holdout in PRIMARY_COHORTS:
                if holdout not in y_by:
                    continue
                train = [cohort for cohort in PRIMARY_COHORTS if cohort != holdout and cohort in y_by]
                y_train = _concat(y_by, train).astype(int)
                p_train = _concat(scores[method][axis], train).astype(float)
                threshold = select_threshold(y_train.to_numpy(dtype=int), p_train.to_numpy(dtype=float))
                y_test = y_by[holdout].astype(int)
                p_test = scores[method][axis][holdout].astype(float)
                for sample_id, prob in p_test.items():
                    parts.append(
                        {
                            "endpoint": "primary_recist",
                            "stratum": "melanoma_core_high_evidence",
                            "cohort": holdout,
                            "sample_id": sample_id,
                            "true_response_label": int(y_test.loc[sample_id]),
                            "response_probability": float(prob),
                            "threshold": float(threshold),
                            "method": method,
                            "axis": axis,
                            "model_name": f"EcoNiche-Opt-{axis}-{method}",
                            "selection_boundary": "outer_lodo_training_threshold_only",
                        }
                    )
            pred = pd.DataFrame(parts)
            if pred.empty:
                continue
            row = {
                "endpoint": "primary_recist",
                "stratum": "melanoma_core_high_evidence",
                "method": method,
                "axis": axis,
                "model_name": f"EcoNiche-Opt-{axis}-{method}",
                "selection_boundary": "outer_lodo_training_threshold_only",
                **_metric_row(
                    pred["true_response_label"].astype(int),
                    pred["response_probability"].astype(float),
                    pred["threshold"].astype(float),
                ),
            }
            rows.append(row)
            pred_rows.extend(parts)
    return pd.DataFrame(rows), pd.DataFrame(pred_rows)


def evaluate_strict_external(scores: dict[str, dict[str, dict[str, pd.Series]]], y_by: dict[str, pd.Series]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    pred_rows = []
    for method in TRANSFORMS:
        for axis in AXES:
            discovery = [cohort for cohort in PRIMARY_COHORTS if cohort in y_by]
            y_train = _concat(y_by, discovery).astype(int)
            p_train = _concat(scores[method][axis], discovery).astype(float)
            threshold = select_threshold(y_train.to_numpy(dtype=int), p_train.to_numpy(dtype=float))
            parts = []
            for cohort in STRICT_EXTERNAL_COHORTS:
                if cohort not in y_by:
                    continue
                y_test = y_by[cohort].astype(int)
                p_test = scores[method][axis][cohort].astype(float)
                for sample_id, prob in p_test.items():
                    parts.append(
                        {
                            "endpoint": "strict_recist",
                            "stratum": "strict_melanoma_pd1_like_external",
                            "cohort": cohort,
                            "sample_id": sample_id,
                            "true_response_label": int(y_test.loc[sample_id]),
                            "response_probability": float(prob),
                            "threshold": float(threshold),
                            "method": method,
                            "axis": axis,
                            "model_name": f"EcoNiche-Opt-{axis}-{method}",
                            "selection_boundary": "discovery_only_threshold_current_external_stress",
                        }
                    )
            pred = pd.DataFrame(parts)
            if pred.empty:
                continue
            rows.append(
                {
                    "endpoint": "strict_recist",
                    "stratum": "strict_melanoma_pd1_like_external",
                    "method": method,
                    "axis": axis,
                    "model_name": f"EcoNiche-Opt-{axis}-{method}",
                    "selection_boundary": "discovery_only_threshold_current_external_stress",
                    **_metric_row(
                        pred["true_response_label"].astype(int),
                        pred["response_probability"].astype(float),
                        float(threshold),
                    ),
                }
            )
            pred_rows.extend(parts)
    return pd.DataFrame(rows), pd.DataFrame(pred_rows)


def write_selection(primary: pd.DataFrame, external: pd.DataFrame, out: Path) -> pd.DataFrame:
    primary_selected = primary.sort_values(["AUROC", "AUPRC", "balanced_accuracy"], ascending=False).iloc[0]
    external_match = external[
        external["method"].astype(str).eq(str(primary_selected["method"]))
        & external["axis"].astype(str).eq(str(primary_selected["axis"]))
    ].iloc[0]
    external_stress_best = external.sort_values(["AUROC", "AUPRC", "balanced_accuracy"], ascending=False).iloc[0]
    rows = [
        {
            "selection_id": "primary_selected_candidate",
            "claim_boundary": "selected_by_primary_lodo_only_not_by_external",
            "method": primary_selected["method"],
            "axis": primary_selected["axis"],
            "primary_AUROC": float(primary_selected["AUROC"]),
            "primary_AUPRC": float(primary_selected["AUPRC"]),
            "primary_balanced_accuracy": float(primary_selected["balanced_accuracy"]),
            "strict_external_AUROC": float(external_match["AUROC"]),
            "strict_external_AUPRC": float(external_match["AUPRC"]),
            "strict_external_balanced_accuracy": float(external_match["balanced_accuracy"]),
        },
        {
            "selection_id": "current_external_stress_best",
            "claim_boundary": "current_external_stress_screen_not_a_locked_selection_claim",
            "method": external_stress_best["method"],
            "axis": external_stress_best["axis"],
            "primary_AUROC": float(
                primary[
                    primary["method"].astype(str).eq(str(external_stress_best["method"]))
                    & primary["axis"].astype(str).eq(str(external_stress_best["axis"]))
                ].iloc[0]["AUROC"]
            ),
            "primary_AUPRC": float(
                primary[
                    primary["method"].astype(str).eq(str(external_stress_best["method"]))
                    & primary["axis"].astype(str).eq(str(external_stress_best["axis"]))
                ].iloc[0]["AUPRC"]
            ),
            "primary_balanced_accuracy": float(
                primary[
                    primary["method"].astype(str).eq(str(external_stress_best["method"]))
                    & primary["axis"].astype(str).eq(str(external_stress_best["axis"]))
                ].iloc[0]["balanced_accuracy"]
            ),
            "strict_external_AUROC": float(external_stress_best["AUROC"]),
            "strict_external_AUPRC": float(external_stress_best["AUPRC"]),
            "strict_external_balanced_accuracy": float(external_stress_best["balanced_accuracy"]),
        },
    ]
    selected = pd.DataFrame(rows)
    selected.to_csv(out / "map4k1_tbx3_transform_selection.tsv", sep="\t", index=False)
    return selected


def write_markdown(selected: pd.DataFrame, out_md: Path) -> None:
    lines = [
        "# MAP4K1-TBX3 Transform Audit",
        "",
        "This audit evaluates locked transformations of the MAP4K1-TBX3 tumor-immune balance axis and related literature-prior extensions.",
        "Primary selection is based on primary LODO only. The current external best row is reported as a stress screen, not as a locked external-selection claim.",
        "",
    ]
    for _, row in selected.iterrows():
        lines.append(
            "- **{}** `{}` `{}`: primary AUROC={:.3f}; strict external AUROC={:.3f}; boundary={}".format(
                row["selection_id"],
                row["method"],
                row["axis"],
                float(row["primary_AUROC"]),
                float(row["strict_external_AUROC"]),
                row["claim_boundary"],
            )
        )
    out_md.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--out", default="results/map4k1_tbx3_transform_audit_20260527")
    args = parser.parse_args()

    out = ROOT / args.out
    out.mkdir(parents=True, exist_ok=True)
    X_by_cohort, _, metadata_by_cohort = load_processed_bulk(ROOT / args.processed_dir)
    primary_data = prepare_endpoint_data(X_by_cohort, metadata_by_cohort, PRIMARY_COHORTS, "primary_recist")
    strict_data = prepare_endpoint_data(X_by_cohort, metadata_by_cohort, [*PRIMARY_COHORTS, *STRICT_EXTERNAL_COHORTS], "strict_recist")
    primary_scores, primary_coverage = build_scores(primary_data.X_by_cohort)
    strict_scores, strict_coverage = build_scores(strict_data.X_by_cohort)
    primary_summary, primary_predictions = evaluate_primary(primary_scores, primary_data.y_response_by_cohort)
    external_summary, external_predictions = evaluate_strict_external(strict_scores, strict_data.y_response_by_cohort)
    primary_summary.to_csv(out / "map4k1_tbx3_transform_primary_summary.tsv", sep="\t", index=False)
    primary_predictions.to_csv(out / "map4k1_tbx3_transform_primary_predictions.tsv", sep="\t", index=False)
    external_summary.to_csv(out / "map4k1_tbx3_transform_strict_external_summary.tsv", sep="\t", index=False)
    external_predictions.to_csv(out / "map4k1_tbx3_transform_strict_external_predictions.tsv", sep="\t", index=False)
    pd.concat([primary_coverage.assign(context="primary"), strict_coverage.assign(context="strict_external")], ignore_index=True).to_csv(
        out / "map4k1_tbx3_transform_gene_coverage.tsv",
        sep="\t",
        index=False,
    )
    selected = write_selection(primary_summary, external_summary, out)
    write_markdown(selected, out / "MAP4K1_TBX3_TRANSFORM_AUDIT.md")
    print(json.dumps(selected.to_dict("records"), ensure_ascii=False))
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
