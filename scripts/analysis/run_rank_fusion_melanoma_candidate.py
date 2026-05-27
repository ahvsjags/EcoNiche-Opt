from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.io import load_processed_bulk
from econiche.metrics import compute_binary_metrics
from econiche_opt.model.endpoint_modules import (
    _concat,
    build_fixed_scores_by_cohort,
    build_module_features_by_cohort,
    default_strata,
    prepare_endpoint_data,
    select_threshold,
)


PRIMARY_ENDPOINT = "primary_recist"
STRICT_ENDPOINT = "strict_recist"
PRIMARY_STRATA = ["melanoma_core_high_evidence", "melanoma_recist_supported_primary"]
DISCOVERY_COHORTS = ["GSE91061", "GSE78220", "PRJEB23709_PD1_PRE"]
STRICT_EXTERNAL_COHORTS = ["GSE145996", "PHS000452_LIU_LIKE_PRE"]


def rank_percentile(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    if values.notna().sum() <= 1:
        return pd.Series(0.5, index=series.index)
    return values.rank(method="average", pct=True).fillna(0.5).astype(float)


def candidate_specs() -> list[dict[str, object]]:
    specs: list[dict[str, object]] = [{"candidate": "rank_module_prior", "weights": {"EcoNiche-Opt-ModulePriorFixed": 1.0}}]
    for cyt_penalty in [-0.25, -0.375, -0.5, -0.625, -0.75, -1.0]:
        specs.append(
            {
                "candidate": f"rank_module_cytotoxic_penalty_{abs(cyt_penalty):.3f}",
                "weights": {"EcoNiche-Opt-ModulePriorFixed": 1.0, "CYT": cyt_penalty},
            }
        )
    for dysfunction_weight in [0.25, 0.5, 0.75]:
        for cyt_penalty in [-0.5, -0.75, -1.0, -1.25]:
            specs.append(
                {
                    "candidate": f"rank_module_dysfunction_{dysfunction_weight:.2f}_cytpen_{abs(cyt_penalty):.2f}",
                    "weights": {
                        "EcoNiche-Opt-ModulePriorFixed": 1.0,
                        "TIDE_dysfunction": dysfunction_weight,
                        "CYT": cyt_penalty,
                    },
                }
            )
    for exclusion_weight in [0.25, 0.5, 0.75]:
        specs.append(
            {
                "candidate": f"rank_module_exclusion_context_{exclusion_weight:.2f}",
                "weights": {"EcoNiche-Opt-ModulePriorFixed": 1.0, "TIDE_exclusion": exclusion_weight},
            }
        )
    return specs


def build_rank_fusion_scores(
    fixed_scores: dict[str, dict[str, pd.Series]],
    spec: dict[str, object],
) -> dict[str, pd.Series]:
    weights = spec["weights"]
    if not isinstance(weights, dict):
        raise TypeError("candidate weights must be a dict")
    out: dict[str, pd.Series] = {}
    for cohort, score_map in fixed_scores.items():
        fused = pd.Series(0.0, index=next(iter(score_map.values())).index)
        for name, weight in weights.items():
            if name not in score_map:
                fused = fused + 0.5 * float(weight)
            else:
                fused = fused + rank_percentile(score_map[name]) * float(weight)
        out[cohort] = rank_percentile(fused)
    return out


def _score_candidate_inner(
    spec: dict[str, object],
    train_cohorts: list[str],
    score_by_cohort: dict[str, pd.Series],
    y_by_cohort: dict[str, pd.Series],
) -> dict[str, object]:
    rows = []
    for holdout in train_cohorts:
        inner_train = [cohort for cohort in train_cohorts if cohort != holdout and cohort in score_by_cohort]
        if not inner_train or holdout not in score_by_cohort:
            continue
        y_train = _concat(y_by_cohort, inner_train).astype(int)
        y_test = y_by_cohort[holdout].astype(int)
        if y_train.nunique() < 2 or y_test.nunique() < 2:
            continue
        train_score = _concat(score_by_cohort, inner_train).astype(float)
        threshold = select_threshold(y_train.to_numpy(dtype=int), train_score.to_numpy(dtype=float))
        metrics = compute_binary_metrics(y_test, score_by_cohort[holdout].astype(float), threshold=threshold)
        rows.append(metrics)
    if not rows:
        return {
            "candidate": spec["candidate"],
            "selection_score": -np.inf,
            "inner_mean_AUROC": np.nan,
            "inner_mean_AUPRC": np.nan,
            "inner_mean_balanced_accuracy": np.nan,
            "inner_mean_ECE": np.nan,
        }
    frame = pd.DataFrame(rows)
    score = (
        frame["AUROC"].mean()
        - 0.20 * frame["AUROC"].std(ddof=0)
        + 0.10 * frame["AUPRC"].mean()
        + 0.10 * frame["balanced_accuracy"].mean()
        - 0.10 * frame["ECE"].mean()
    )
    return {
        "candidate": spec["candidate"],
        "selection_score": float(score),
        "inner_mean_AUROC": float(frame["AUROC"].mean()),
        "inner_mean_AUPRC": float(frame["AUPRC"].mean()),
        "inner_mean_balanced_accuracy": float(frame["balanced_accuracy"].mean()),
        "inner_mean_ECE": float(frame["ECE"].mean()),
        "inner_sd_AUROC": float(frame["AUROC"].std(ddof=0)),
    }


def _select_spec(
    specs: list[dict[str, object]],
    train_cohorts: list[str],
    fixed_scores: dict[str, dict[str, pd.Series]],
    y_by_cohort: dict[str, pd.Series],
) -> tuple[dict[str, object], pd.DataFrame]:
    rows = []
    scored: list[tuple[dict[str, object], dict[str, object]]] = []
    for spec in specs:
        score_by_cohort = build_rank_fusion_scores(fixed_scores, spec)
        row = _score_candidate_inner(spec, train_cohorts, score_by_cohort, y_by_cohort)
        rows.append({**row, "weights": json.dumps(spec["weights"], sort_keys=True)})
        scored.append((spec, row))
    finite = [(spec, row) for spec, row in scored if np.isfinite(float(row.get("selection_score", np.nan)))]
    best_spec = max(finite, key=lambda item: float(item[1]["selection_score"]))[0] if finite else specs[0]
    return best_spec, pd.DataFrame(rows)


def _summarize_predictions(predictions: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    rows = []
    if predictions.empty:
        return pd.DataFrame()
    for keys, frame in predictions.groupby(group_cols):
        if not isinstance(keys, tuple):
            keys = (keys,)
        y = frame["true_response_label"].astype(int)
        prob = frame["response_probability"].astype(float)
        threshold = float(frame["threshold"].median()) if "threshold" in frame else 0.5
        metrics = compute_binary_metrics(y, prob, threshold=threshold)
        rows.append(
            {
                **dict(zip(group_cols, keys)),
                **metrics,
                "n_samples": int(len(frame)),
                "n_cohorts": int(frame["cohort"].nunique()) if "cohort" in frame else 1,
                "n_responders": int(y.sum()),
                "n_nonresponders": int((y == 0).sum()),
                "response_prevalence": float(y.mean()),
                "AUPRC_minus_prevalence": float(metrics["AUPRC"] - y.mean()),
                "selected_candidates": ",".join(sorted(frame["selected_candidate"].astype(str).unique())),
                "median_threshold": threshold,
            }
        )
    return pd.DataFrame(rows)


def evaluate_primary(
    X_by_cohort: dict[str, pd.DataFrame],
    metadata_by_cohort: dict[str, pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    active = [cohort for cohort in X_by_cohort if not cohort.startswith("demo_cohort_")]
    strata = default_strata(active)
    specs = candidate_specs()
    metric_rows = []
    prediction_rows = []
    selection_rows = []
    for stratum in PRIMARY_STRATA:
        cohorts = strata[stratum]["cohorts"]
        endpoint_data = prepare_endpoint_data(X_by_cohort, metadata_by_cohort, cohorts, PRIMARY_ENDPOINT)
        module_features, _ = build_module_features_by_cohort(endpoint_data.X_by_cohort)
        fixed_scores = build_fixed_scores_by_cohort(endpoint_data.X_by_cohort, module_features)
        for holdout in strata[stratum]["holdouts"]:
            train_cohorts = [cohort for cohort in strata[stratum]["train_pool"] if cohort != holdout and cohort in fixed_scores]
            if holdout not in fixed_scores or len(train_cohorts) < 2:
                continue
            y_train = _concat(endpoint_data.y_response_by_cohort, train_cohorts).astype(int)
            y_test = endpoint_data.y_response_by_cohort[holdout].astype(int)
            if y_train.nunique() < 2 or y_test.nunique() < 2:
                continue
            best_spec, selection = _select_spec(specs, train_cohorts, fixed_scores, endpoint_data.y_response_by_cohort)
            selection.insert(0, "endpoint", PRIMARY_ENDPOINT)
            selection.insert(1, "stratum", stratum)
            selection.insert(2, "holdout", holdout)
            selection_rows.append(selection)
            score_by_cohort = build_rank_fusion_scores(fixed_scores, best_spec)
            train_score = _concat(score_by_cohort, train_cohorts).astype(float)
            threshold = select_threshold(y_train.to_numpy(dtype=int), train_score.to_numpy(dtype=float))
            test_score = score_by_cohort[holdout].astype(float)
            metrics = compute_binary_metrics(y_test, test_score, threshold=threshold)
            metric_rows.append(
                {
                    **metrics,
                    "endpoint": PRIMARY_ENDPOINT,
                    "stratum": stratum,
                    "cohort": holdout,
                    "selected_candidate": best_spec["candidate"],
                    "weights": json.dumps(best_spec["weights"], sort_keys=True),
                    "threshold": float(threshold),
                    "n_samples": int(len(y_test)),
                    "n_responders": int(y_test.sum()),
                    "n_nonresponders": int((y_test == 0).sum()),
                    "train_cohorts": ",".join(train_cohorts),
                    "selection_boundary": "outer_lodo_training_only_rank_fusion",
                }
            )
            for sample_id, score in test_score.items():
                prediction_rows.append(
                    {
                        "endpoint": PRIMARY_ENDPOINT,
                        "stratum": stratum,
                        "cohort": holdout,
                        "sample_id": sample_id,
                        "true_response_label": int(y_test.loc[sample_id]),
                        "response_probability": float(score),
                        "pred_response_label": int(float(score) >= threshold),
                        "selected_candidate": best_spec["candidate"],
                        "threshold": float(threshold),
                    }
                )
    selection_df = pd.concat(selection_rows, ignore_index=True) if selection_rows else pd.DataFrame()
    return pd.DataFrame(metric_rows), pd.DataFrame(prediction_rows), selection_df


def evaluate_strict_external(
    X_by_cohort: dict[str, pd.DataFrame],
    metadata_by_cohort: dict[str, pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    endpoint_data = prepare_endpoint_data(
        X_by_cohort,
        metadata_by_cohort,
        [*DISCOVERY_COHORTS, *STRICT_EXTERNAL_COHORTS],
        STRICT_ENDPOINT,
    )
    module_features, _ = build_module_features_by_cohort(endpoint_data.X_by_cohort)
    fixed_scores = build_fixed_scores_by_cohort(endpoint_data.X_by_cohort, module_features)
    train_cohorts = [cohort for cohort in DISCOVERY_COHORTS if cohort in fixed_scores]
    best_spec, selection = _select_spec(candidate_specs(), train_cohorts, fixed_scores, endpoint_data.y_response_by_cohort)
    selection.insert(0, "endpoint", STRICT_ENDPOINT)
    selection.insert(1, "stratum", "strict_melanoma_pd1_like_external")
    selection.insert(2, "holdout", "discovery_inner_lodo")
    score_by_cohort = build_rank_fusion_scores(fixed_scores, best_spec)
    y_train = _concat(endpoint_data.y_response_by_cohort, train_cohorts).astype(int)
    train_score = _concat(score_by_cohort, train_cohorts).astype(float)
    threshold = select_threshold(y_train.to_numpy(dtype=int), train_score.to_numpy(dtype=float))

    metric_rows = []
    prediction_rows = []
    for holdout in STRICT_EXTERNAL_COHORTS:
        if holdout not in score_by_cohort or holdout not in endpoint_data.y_response_by_cohort:
            continue
        y_test = endpoint_data.y_response_by_cohort[holdout].astype(int)
        if y_test.nunique() < 2:
            continue
        score = score_by_cohort[holdout].astype(float)
        metrics = compute_binary_metrics(y_test, score, threshold=threshold)
        metric_rows.append(
            {
                **metrics,
                "endpoint": STRICT_ENDPOINT,
                "stratum": "strict_melanoma_pd1_like_external",
                "cohort": holdout,
                "selected_candidate": best_spec["candidate"],
                "weights": json.dumps(best_spec["weights"], sort_keys=True),
                "threshold": float(threshold),
                "n_samples": int(len(y_test)),
                "n_responders": int(y_test.sum()),
                "n_nonresponders": int((y_test == 0).sum()),
                "train_cohorts": ",".join(train_cohorts),
                "selection_boundary": "discovery_only_rank_fusion_no_external_selection",
            }
        )
        for sample_id, value in score.items():
            prediction_rows.append(
                {
                    "endpoint": STRICT_ENDPOINT,
                    "stratum": "strict_melanoma_pd1_like_external",
                    "cohort": holdout,
                    "sample_id": sample_id,
                    "true_response_label": int(y_test.loc[sample_id]),
                    "response_probability": float(value),
                    "pred_response_label": int(float(value) >= threshold),
                    "selected_candidate": best_spec["candidate"],
                    "threshold": float(threshold),
                }
            )
    return pd.DataFrame(metric_rows), pd.DataFrame(prediction_rows), selection


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--out", default="results/rank_fusion_melanoma_candidate_20260527")
    args = parser.parse_args()

    X_by_cohort, _, metadata_by_cohort = load_processed_bulk(ROOT / args.processed_dir)
    out = ROOT / args.out
    out.mkdir(parents=True, exist_ok=True)

    primary_metrics, primary_predictions, primary_selection = evaluate_primary(X_by_cohort, metadata_by_cohort)
    external_metrics, external_predictions, external_selection = evaluate_strict_external(X_by_cohort, metadata_by_cohort)
    primary_summary = _summarize_predictions(primary_predictions, ["endpoint", "stratum"])
    external_summary = _summarize_predictions(external_predictions, ["endpoint", "stratum"])

    primary_metrics.to_csv(out / "rank_fusion_primary_fold_metrics.tsv", sep="\t", index=False)
    primary_predictions.to_csv(out / "rank_fusion_primary_predictions.tsv", sep="\t", index=False)
    primary_selection.to_csv(out / "rank_fusion_primary_selection.tsv", sep="\t", index=False)
    primary_summary.to_csv(out / "rank_fusion_primary_summary.tsv", sep="\t", index=False)
    external_metrics.to_csv(out / "rank_fusion_strict_external_fold_metrics.tsv", sep="\t", index=False)
    external_predictions.to_csv(out / "rank_fusion_strict_external_predictions.tsv", sep="\t", index=False)
    external_selection.to_csv(out / "rank_fusion_strict_external_selection.tsv", sep="\t", index=False)
    external_summary.to_csv(out / "rank_fusion_strict_external_summary.tsv", sep="\t", index=False)

    manifest = pd.DataFrame(
        [
            {"artifact": path.name, "path": str(path), "n_bytes": path.stat().st_size}
            for path in sorted(out.glob("*.tsv"))
        ]
    )
    manifest.to_csv(out / "rank_fusion_output_manifest.tsv", sep="\t", index=False)
    print(json.dumps({"out": str(out), "primary_rows": len(primary_predictions), "external_rows": len(external_predictions)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
