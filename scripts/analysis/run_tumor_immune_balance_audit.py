from __future__ import annotations

import argparse
import json
import math
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


PRIMARY_ENDPOINT = "primary_recist"
STRICT_ENDPOINT = "strict_recist"
PRIMARY_STRATUM = "melanoma_core_high_evidence"
MAIN_MODEL = "EcoNiche-Opt-HeuristicEcology"
LOCKED_EXTERNAL_MODEL = "EcoNiche-Opt-HeuristicEcology-LockedPanel"
PRIMARY_COHORTS = ["GSE91061", "GSE78220", "PRJEB23709_PD1_PRE"]
STRICT_EXTERNAL_COHORTS = ["GSE145996", "PHS000452_LIU_LIKE_PRE"]


BALANCE_AXES: dict[str, tuple[list[str], list[str]]] = {
    "MAP4K1_minus_TBX3": (["MAP4K1"], ["TBX3"]),
    "map4k1_apm_minus_tbx3_ipres": (["MAP4K1", "HLA-A", "B2M", "TAP1", "CXCL9"], ["TBX3", "AXL", "VIM", "ZEB1", "FN1"]),
    "cyt_apm_minus_ipres": (["GZMA", "PRF1", "HLA-A", "B2M", "TAP1", "CXCL9"], ["AXL", "VIM", "ZEB1", "COL1A1", "FN1"]),
    "immune_effector_minus_tumor": (
        ["MAP4K1", "CD8A", "CD8B", "GZMB", "PRF1", "NKG7", "CXCL9", "CXCL10", "IFNG"],
        ["TBX3", "AXL", "MITF", "SOX10", "WNT5A", "PAK4"],
    ),
}

BLEND_SPECS: list[dict[str, object]] = [
    {
        "model_name": "EcoNiche-Opt-TumorImmuneBalancePair",
        "axis": "MAP4K1_minus_TBX3",
        "balance_weight": 1.0,
        "main_weight": 0.0,
        "selection_boundary": "literature_prior_fixed_no_label_fit",
    },
    {
        "model_name": "EcoNiche-Opt-TumorImmuneBalancePairBlend50",
        "axis": "MAP4K1_minus_TBX3",
        "balance_weight": 0.5,
        "main_weight": 0.5,
        "selection_boundary": "literature_prior_fixed_50_50_with_main_score",
    },
    {
        "model_name": "EcoNiche-Opt-MAP4K1APM-TBX3IPRESBalance",
        "axis": "map4k1_apm_minus_tbx3_ipres",
        "balance_weight": 1.0,
        "main_weight": 0.0,
        "selection_boundary": "biological_axis_fixed_no_external_selection",
    },
    {
        "model_name": "EcoNiche-Opt-MAP4K1APM-TBX3IPRESBlend50",
        "axis": "map4k1_apm_minus_tbx3_ipres",
        "balance_weight": 0.5,
        "main_weight": 0.5,
        "selection_boundary": "biological_axis_fixed_50_50_with_main_score",
    },
    {
        "model_name": "EcoNiche-Opt-CYTAPM-IPRESBlend50",
        "axis": "cyt_apm_minus_ipres",
        "balance_weight": 0.5,
        "main_weight": 0.5,
        "selection_boundary": "biological_axis_fixed_50_50_with_main_score",
    },
]


def _rank_subset(X: pd.DataFrame) -> pd.DataFrame:
    genes = sorted({gene for pair in BALANCE_AXES.values() for genes in pair for gene in genes})
    available = [gene for gene in genes if gene in X.columns]
    if not available:
        return pd.DataFrame(index=X.index)
    return X[available].apply(pd.to_numeric, errors="coerce").rank(axis=0, pct=True).fillna(0.5)


def balance_score(ranked: pd.DataFrame, positive_genes: list[str], negative_genes: list[str]) -> tuple[pd.Series, dict[str, object]]:
    pos = [gene for gene in positive_genes if gene in ranked.columns]
    neg = [gene for gene in negative_genes if gene in ranked.columns]
    raw = pd.Series(0.0, index=ranked.index)
    if pos:
        raw = raw + ranked[pos].mean(axis=1)
    if neg:
        raw = raw - ranked[neg].mean(axis=1)
    if raw.max() > raw.min():
        scaled = (raw - raw.min()) / (raw.max() - raw.min())
    else:
        scaled = pd.Series(0.5, index=ranked.index)
    return scaled.astype(float), {
        "n_positive_genes_available": len(pos),
        "n_negative_genes_available": len(neg),
        "positive_genes_available": ",".join(pos),
        "negative_genes_available": ",".join(neg),
    }


def build_axis_scores(
    X_by_cohort: dict[str, pd.DataFrame],
) -> tuple[dict[str, dict[str, pd.Series]], pd.DataFrame]:
    scores: dict[str, dict[str, pd.Series]] = {axis: {} for axis in BALANCE_AXES}
    coverage_rows = []
    for cohort, X in X_by_cohort.items():
        ranked = _rank_subset(X)
        for axis, (pos, neg) in BALANCE_AXES.items():
            score, coverage = balance_score(ranked, pos, neg)
            scores[axis][cohort] = score
            coverage_rows.append({"axis": axis, "cohort": cohort, "n_samples": len(score), **coverage})
    return scores, pd.DataFrame(coverage_rows)


def _main_primary_scores(path: Path) -> dict[str, pd.Series]:
    frame = pd.read_csv(path, sep="\t")
    frame = frame[
        frame["endpoint"].astype(str).eq(PRIMARY_ENDPOINT)
        & frame["stratum"].astype(str).eq(PRIMARY_STRATUM)
        & frame["model_name"].astype(str).eq(MAIN_MODEL)
    ].copy()
    return {
        cohort: group.set_index("sample_id")["response_probability"].astype(float)
        for cohort, group in frame.groupby("cohort")
    }


def _main_external_scores(path: Path) -> dict[str, pd.Series]:
    frame = pd.read_csv(path, sep="\t")
    frame = frame[frame["model_name"].astype(str).eq(LOCKED_EXTERNAL_MODEL)].copy()
    return {
        cohort: group.set_index("sample_id")["response_probability"].astype(float)
        for cohort, group in frame.groupby("cohort")
    }


def _blend_score(balance: pd.Series, main: pd.Series | None, balance_weight: float, main_weight: float) -> pd.Series:
    aligned_main = pd.Series(0.5, index=balance.index) if main is None else main.reindex(balance.index).fillna(0.5)
    score = balance_weight * balance + main_weight * aligned_main
    if score.max() > score.min():
        score = (score - score.min()) / (score.max() - score.min())
    return score.astype(float)


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
    return metrics


def evaluate_primary(
    X_by_cohort: dict[str, pd.DataFrame],
    metadata_by_cohort: dict[str, pd.DataFrame],
    main_scores: dict[str, pd.Series],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    endpoint_data = prepare_endpoint_data(X_by_cohort, metadata_by_cohort, PRIMARY_COHORTS, PRIMARY_ENDPOINT)
    axis_scores, coverage = build_axis_scores(endpoint_data.X_by_cohort)
    rows = []
    prediction_rows = []
    for spec in BLEND_SPECS:
        parts = []
        for holdout in PRIMARY_COHORTS:
            if holdout not in endpoint_data.y_response_by_cohort:
                continue
            train = [cohort for cohort in PRIMARY_COHORTS if cohort != holdout and cohort in endpoint_data.y_response_by_cohort]
            axis = str(spec["axis"])
            train_scores = {
                cohort: _blend_score(
                    axis_scores[axis][cohort],
                    main_scores.get(cohort),
                    float(spec["balance_weight"]),
                    float(spec["main_weight"]),
                )
                for cohort in train
            }
            y_train = _concat(endpoint_data.y_response_by_cohort, train).astype(int)
            p_train = _concat(train_scores, train).astype(float)
            threshold = select_threshold(y_train.to_numpy(dtype=int), p_train.to_numpy(dtype=float))
            p_test = _blend_score(
                axis_scores[axis][holdout],
                main_scores.get(holdout),
                float(spec["balance_weight"]),
                float(spec["main_weight"]),
            )
            y_test = endpoint_data.y_response_by_cohort[holdout].astype(int)
            for sample_id, prob in p_test.items():
                parts.append(
                    {
                        "endpoint": PRIMARY_ENDPOINT,
                        "stratum": PRIMARY_STRATUM,
                        "cohort": holdout,
                        "sample_id": sample_id,
                        "true_response_label": int(y_test.loc[sample_id]),
                        "response_probability": float(prob),
                        "threshold": float(threshold),
                        "model_name": spec["model_name"],
                        "axis": axis,
                        "selection_boundary": spec["selection_boundary"],
                    }
                )
        pred = pd.DataFrame(parts)
        if pred.empty:
            continue
        y = pred["true_response_label"].astype(int)
        p = pred["response_probability"].astype(float)
        thresholds = pred["threshold"].astype(float)
        rows.append(
            {
                "endpoint": PRIMARY_ENDPOINT,
                "stratum": PRIMARY_STRATUM,
                "model_name": spec["model_name"],
                "axis": spec["axis"],
                "selection_boundary": spec["selection_boundary"],
                "balance_weight": spec["balance_weight"],
                "main_weight": spec["main_weight"],
                **_metric_row(y, p, thresholds),
            }
        )
        prediction_rows.extend(parts)
    return pd.DataFrame(rows), pd.DataFrame(prediction_rows), coverage


def evaluate_strict_external(
    X_by_cohort: dict[str, pd.DataFrame],
    metadata_by_cohort: dict[str, pd.DataFrame],
    primary_main_scores: dict[str, pd.Series],
    external_main_scores: dict[str, pd.Series],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cohorts = [*PRIMARY_COHORTS, *STRICT_EXTERNAL_COHORTS]
    endpoint_data = prepare_endpoint_data(X_by_cohort, metadata_by_cohort, cohorts, STRICT_ENDPOINT)
    axis_scores, coverage = build_axis_scores(endpoint_data.X_by_cohort)
    rows = []
    prediction_rows = []
    for spec in BLEND_SPECS:
        axis = str(spec["axis"])
        discovery_scores = {
            cohort: _blend_score(
                axis_scores[axis][cohort],
                primary_main_scores.get(cohort),
                float(spec["balance_weight"]),
                float(spec["main_weight"]),
            )
            for cohort in PRIMARY_COHORTS
            if cohort in endpoint_data.y_response_by_cohort
        }
        discovery = list(discovery_scores)
        y_train = _concat(endpoint_data.y_response_by_cohort, discovery).astype(int)
        p_train = _concat(discovery_scores, discovery).astype(float)
        threshold = select_threshold(y_train.to_numpy(dtype=int), p_train.to_numpy(dtype=float))
        parts = []
        for cohort in STRICT_EXTERNAL_COHORTS:
            if cohort not in endpoint_data.y_response_by_cohort:
                continue
            p_test = _blend_score(
                axis_scores[axis][cohort],
                external_main_scores.get(cohort),
                float(spec["balance_weight"]),
                float(spec["main_weight"]),
            )
            y_test = endpoint_data.y_response_by_cohort[cohort].astype(int)
            for sample_id, prob in p_test.items():
                parts.append(
                    {
                        "endpoint": STRICT_ENDPOINT,
                        "stratum": "strict_melanoma_pd1_like_external",
                        "cohort": cohort,
                        "sample_id": sample_id,
                        "true_response_label": int(y_test.loc[sample_id]),
                        "response_probability": float(prob),
                        "threshold": float(threshold),
                        "model_name": spec["model_name"],
                        "axis": axis,
                        "selection_boundary": "discovery_only_threshold_locked_external_scoring",
                    }
                )
        pred = pd.DataFrame(parts)
        if pred.empty:
            continue
        rows.append(
            {
                "endpoint": STRICT_ENDPOINT,
                "stratum": "strict_melanoma_pd1_like_external",
                "model_name": spec["model_name"],
                "axis": spec["axis"],
                "selection_boundary": "discovery_only_threshold_locked_external_scoring",
                "balance_weight": spec["balance_weight"],
                "main_weight": spec["main_weight"],
                **_metric_row(
                    pred["true_response_label"].astype(int),
                    pred["response_probability"].astype(float),
                    float(threshold),
                ),
            }
        )
        prediction_rows.extend(parts)
    return pd.DataFrame(rows), pd.DataFrame(prediction_rows), coverage


def write_markdown(primary: pd.DataFrame, external: pd.DataFrame, out_md: Path) -> None:
    lines = [
        "# Tumor-Immune Balance Audit",
        "",
        "This audit adds a literature-prior MAP4K1-TBX3 tumor-immune balance axis and related ecological balance variants.",
        "The strict external cohorts are never used for training, feature selection, thresholding, calibration, or model selection.",
        "",
        "## Primary LODO",
        "",
    ]
    for _, row in primary.sort_values("AUROC", ascending=False).iterrows():
        lines.append(f"- `{row['model_name']}`: AUROC={float(row['AUROC']):.3f}, AUPRC={float(row['AUPRC']):.3f}, BA={float(row['balanced_accuracy']):.3f}")
    lines.extend(["", "## Strict External", ""])
    for _, row in external.sort_values("AUROC", ascending=False).iterrows():
        lines.append(f"- `{row['model_name']}`: AUROC={float(row['AUROC']):.3f}, AUPRC={float(row['AUPRC']):.3f}, BA={float(row['balanced_accuracy']):.3f}")
    out_md.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument(
        "--primary-predictions",
        default="results/endpoint_modules_heuristic_deep_primary_20260519/endpoint_module_predictions.tsv",
    )
    parser.add_argument("--external-predictions", default="results/pd1_like_external_rescue/pd1_like_rescue_predictions.tsv")
    parser.add_argument("--out", default="results/tumor_immune_balance_audit_20260527")
    args = parser.parse_args()

    out = ROOT / args.out
    out.mkdir(parents=True, exist_ok=True)
    X_by_cohort, _, metadata_by_cohort = load_processed_bulk(ROOT / args.processed_dir)
    primary_main = _main_primary_scores(ROOT / args.primary_predictions)
    external_main = _main_external_scores(ROOT / args.external_predictions)
    primary_summary, primary_predictions, primary_coverage = evaluate_primary(X_by_cohort, metadata_by_cohort, primary_main)
    external_summary, external_predictions, external_coverage = evaluate_strict_external(
        X_by_cohort,
        metadata_by_cohort,
        primary_main,
        external_main,
    )
    primary_summary.to_csv(out / "tumor_immune_balance_primary_summary.tsv", sep="\t", index=False)
    primary_predictions.to_csv(out / "tumor_immune_balance_primary_predictions.tsv", sep="\t", index=False)
    external_summary.to_csv(out / "tumor_immune_balance_strict_external_summary.tsv", sep="\t", index=False)
    external_predictions.to_csv(out / "tumor_immune_balance_strict_external_predictions.tsv", sep="\t", index=False)
    pd.concat([primary_coverage.assign(context="primary"), external_coverage.assign(context="strict_external")], ignore_index=True).to_csv(
        out / "tumor_immune_balance_gene_coverage.tsv",
        sep="\t",
        index=False,
    )
    write_markdown(primary_summary, external_summary, out / "TUMOR_IMMUNE_BALANCE_AUDIT.md")
    best_primary = primary_summary.sort_values("AUROC", ascending=False).iloc[0].to_dict()
    best_external = external_summary.sort_values("AUROC", ascending=False).iloc[0].to_dict()
    print(
        json.dumps(
            {
                "best_primary": {
                    "model_name": best_primary["model_name"],
                    "AUROC": best_primary["AUROC"],
                    "balanced_accuracy": best_primary["balanced_accuracy"],
                },
                "best_external": {
                    "model_name": best_external["model_name"],
                    "AUROC": best_external["AUROC"],
                    "balanced_accuracy": best_external["balanced_accuracy"],
                },
            },
            ensure_ascii=False,
        )
    )
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
