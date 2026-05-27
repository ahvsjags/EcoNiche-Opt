from __future__ import annotations

import argparse
import sys
import warnings
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
    _candidate_specs,
    _concat,
    _predict_candidate,
    _score_inner_candidate,
    build_module_features_by_cohort,
    default_strata,
    prepare_endpoint_data,
    select_threshold,
)


DISCOVERY_COHORTS = ["GSE91061", "GSE78220", "PRJEB23709_PD1_PRE"]
STRICT_EXTERNAL_COHORTS = ["GSE145996", "PHS000452_LIU_LIKE_PRE"]
PRIMARY_ENDPOINT = "primary_recist"
STRICT_ENDPOINT = "strict_recist"
PRIMARY_STRATA = ["melanoma_core_high_evidence", "melanoma_recist_supported_primary"]


def _finite_selection_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    return [row for row in rows if np.isfinite(float(row.get("selection_score", np.nan)))]


def _selection_specs() -> list[dict[str, object]]:
    return [spec for spec in _candidate_specs(optimize_word=False) if str(spec.get("kind")) != "word_ecology"]


def _best_spec(inner_scores: list[dict[str, object]], specs: list[dict[str, object]]) -> tuple[dict[str, object], dict[str, object]]:
    valid = _finite_selection_rows(inner_scores)
    best = max(valid, key=lambda row: float(row["selection_score"])) if valid else {"candidate": "module_prior_composite"}
    spec = next(spec for spec in specs if str(spec["candidate"]) == str(best["candidate"]))
    return best, spec


def _evaluate_lodo(
    X_by_cohort: dict[str, pd.DataFrame],
    metadata_by_cohort: dict[str, pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    active = [cohort for cohort in X_by_cohort if not cohort.startswith("demo_cohort_")]
    strata = default_strata(active)
    specs = _selection_specs()
    metric_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    inner_rows: list[dict[str, object]] = []

    for stratum in PRIMARY_STRATA:
        cohorts = strata[stratum]["cohorts"]
        endpoint_data = prepare_endpoint_data(X_by_cohort, metadata_by_cohort, cohorts, PRIMARY_ENDPOINT)
        module_features, _ = build_module_features_by_cohort(endpoint_data.X_by_cohort)
        for holdout in strata[stratum]["holdouts"]:
            if holdout not in module_features:
                continue
            train_cohorts = [cohort for cohort in strata[stratum]["train_pool"] if cohort != holdout and cohort in module_features]
            if len(train_cohorts) < 2:
                continue
            y_train = _concat(endpoint_data.y_response_by_cohort, train_cohorts).astype(int)
            y_test = endpoint_data.y_response_by_cohort[holdout].astype(int)
            if y_train.nunique() < 2 or y_test.nunique() < 2:
                continue
            inner_scores = [
                _score_inner_candidate(spec, train_cohorts, module_features, endpoint_data.y_response_by_cohort)
                for spec in specs
            ]
            for row in inner_scores:
                inner_rows.append(
                    {
                        "endpoint": PRIMARY_ENDPOINT,
                        "stratum": stratum,
                        "holdout": holdout,
                        "candidate": row.get("candidate"),
                        "inner_mean_AUROC": row.get("inner_mean_AUROC"),
                        "inner_mean_AUPRC": row.get("inner_mean_AUPRC"),
                        "inner_mean_balanced_accuracy": row.get("inner_mean_balanced_accuracy"),
                        "inner_mean_ECE": row.get("inner_mean_ECE"),
                        "selection_score": row.get("selection_score"),
                    }
                )
            best, spec = _best_spec(inner_scores, specs)
            train_prob, test_prob, _ = _predict_candidate(
                spec,
                train_cohorts,
                holdout,
                module_features,
                endpoint_data.y_response_by_cohort,
            )
            threshold = select_threshold(y_train.to_numpy(dtype=int), np.asarray(train_prob, dtype=float))
            metrics = compute_binary_metrics(y_test, test_prob, threshold=threshold)
            metric_rows.append(
                {
                    **metrics,
                    "endpoint": PRIMARY_ENDPOINT,
                    "stratum": stratum,
                    "cohort": holdout,
                    "selected_candidate": best.get("candidate"),
                    "n_samples": int(len(y_test)),
                    "n_responders": int(y_test.sum()),
                    "n_nonresponders": int((y_test == 0).sum()),
                    "response_prevalence": float(y_test.mean()),
                    "threshold": float(threshold),
                    "train_cohorts": ",".join(train_cohorts),
                    "selection_boundary": "training_only_inner_lodo",
                }
            )
            for sample_id, prob in zip(y_test.index, test_prob):
                prediction_rows.append(
                    {
                        "endpoint": PRIMARY_ENDPOINT,
                        "stratum": stratum,
                        "cohort": holdout,
                        "sample_id": sample_id,
                        "true_response_label": int(y_test.loc[sample_id]),
                        "response_probability": float(prob),
                        "selected_candidate": best.get("candidate"),
                        "threshold": float(threshold),
                    }
                )
    return pd.DataFrame(metric_rows), pd.DataFrame(prediction_rows), pd.DataFrame(inner_rows)


def _summarize_lodo_predictions(predictions: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if predictions.empty:
        return pd.DataFrame()
    for (endpoint, stratum), frame in predictions.groupby(["endpoint", "stratum"]):
        y = frame["true_response_label"].astype(int)
        prob = frame["response_probability"].astype(float)
        metrics = compute_binary_metrics(y, prob)
        rows.append(
            {
                **metrics,
                "endpoint": endpoint,
                "stratum": stratum,
                "n_samples": int(len(frame)),
                "n_cohorts": int(frame["cohort"].nunique()),
                "n_responders": int(y.sum()),
                "n_nonresponders": int((y == 0).sum()),
                "response_prevalence": float(y.mean()),
                "AUPRC_minus_prevalence": float(metrics["AUPRC"] - y.mean()),
                "selected_candidates": ",".join(sorted(frame["selected_candidate"].astype(str).unique())),
                "selection_boundary": "training_only_inner_lodo",
            }
        )
    return pd.DataFrame(rows).sort_values(["endpoint", "stratum"]).reset_index(drop=True)


def _evaluate_strict_external(
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
    train_cohorts = [cohort for cohort in DISCOVERY_COHORTS if cohort in module_features]
    specs = _selection_specs()
    inner_scores = [
        _score_inner_candidate(spec, train_cohorts, module_features, endpoint_data.y_response_by_cohort)
        for spec in specs
    ]
    best, spec = _best_spec(inner_scores, specs)
    y_train = _concat(endpoint_data.y_response_by_cohort, train_cohorts).astype(int)

    metric_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    for holdout in STRICT_EXTERNAL_COHORTS:
        if holdout not in module_features:
            continue
        y_test = endpoint_data.y_response_by_cohort[holdout].astype(int)
        if y_test.nunique() < 2:
            continue
        train_prob, test_prob, _ = _predict_candidate(
            spec,
            train_cohorts,
            holdout,
            module_features,
            endpoint_data.y_response_by_cohort,
        )
        threshold = select_threshold(y_train.to_numpy(dtype=int), np.asarray(train_prob, dtype=float))
        metrics = compute_binary_metrics(y_test, test_prob, threshold=threshold)
        metric_rows.append(
            {
                **metrics,
                "endpoint": STRICT_ENDPOINT,
                "cohort": holdout,
                "selected_candidate": best.get("candidate"),
                "n_samples": int(len(y_test)),
                "n_responders": int(y_test.sum()),
                "n_nonresponders": int((y_test == 0).sum()),
                "response_prevalence": float(y_test.mean()),
                "threshold": float(threshold),
                "train_cohorts": ",".join(train_cohorts),
                "selection_boundary": "discovery_only_inner_lodo_no_external_selection",
            }
        )
        for sample_id, prob in zip(y_test.index, test_prob):
            prediction_rows.append(
                {
                    "endpoint": STRICT_ENDPOINT,
                    "cohort": holdout,
                    "sample_id": sample_id,
                    "true_response_label": int(y_test.loc[sample_id]),
                    "response_probability": float(prob),
                    "selected_candidate": best.get("candidate"),
                    "threshold": float(threshold),
                }
            )
    inner = pd.DataFrame(inner_scores)
    if not inner.empty:
        inner.insert(0, "endpoint", STRICT_ENDPOINT)
        inner.insert(1, "selection_set", ",".join(train_cohorts))
        inner["selected_candidate"] = best.get("candidate")
    return pd.DataFrame(metric_rows), pd.DataFrame(prediction_rows), inner


def _summarize_external_predictions(predictions: pd.DataFrame) -> pd.DataFrame:
    if predictions.empty:
        return pd.DataFrame()
    y = predictions["true_response_label"].astype(int)
    prob = predictions["response_probability"].astype(float)
    metrics = compute_binary_metrics(y, prob)
    return pd.DataFrame(
        [
            {
                **metrics,
                "endpoint": STRICT_ENDPOINT,
                "cohort_set": "+".join(STRICT_EXTERNAL_COHORTS),
                "selected_candidate": ",".join(sorted(predictions["selected_candidate"].astype(str).unique())),
                "n_samples": int(len(predictions)),
                "n_cohorts": int(predictions["cohort"].nunique()),
                "n_responders": int(y.sum()),
                "n_nonresponders": int((y == 0).sum()),
                "response_prevalence": float(y.mean()),
                "AUPRC_minus_prevalence": float(metrics["AUPRC"] - y.mean()),
                "selection_boundary": "discovery_only_inner_lodo_no_external_selection",
            }
        ]
    )


def _write_audit(out_dir: Path, lodo_summary: pd.DataFrame, external_summary: pd.DataFrame, external_inner: pd.DataFrame) -> None:
    lines = [
        "# Training-only Melanoma Candidate Search Audit",
        "",
        "This registered analysis asks whether a leakage-safe candidate-selection rule can improve the primary melanoma and strict external evidence without using external labels for feature selection, thresholding, calibration, or model selection.",
        "",
        "## Selection Rule",
        "",
        "- Candidate families are restricted to fixed module composites and sparse module-level logistic models.",
        "- For primary LODO, each holdout fold selects a candidate from the remaining training cohorts only.",
        "- For strict external testing, the selected candidate is chosen from GSE91061, GSE78220 and PRJEB23709_PD1_PRE only; GSE145996 and PHS000452_LIU_LIKE_PRE labels are read only after the candidate and threshold rule are fixed.",
        "",
        "## Primary LODO Summary",
        "",
    ]
    if lodo_summary.empty:
        lines.append("No LODO rows were produced.")
    else:
        for _, row in lodo_summary.iterrows():
            lines.append(
                f"- {row['stratum']}: selected={row['selected_candidates']}; n={int(row['n_samples'])}; "
                f"AUROC={row['AUROC']:.3f}; AUPRC={row['AUPRC']:.3f}; "
                f"balanced accuracy={row['balanced_accuracy']:.3f}; ECE={row['ECE']:.3f}."
            )
    lines.extend(["", "## Strict External Summary", ""])
    if external_summary.empty:
        lines.append("No strict external rows were produced.")
    else:
        row = external_summary.iloc[0]
        lines.append(
            f"- {row['cohort_set']}: selected={row['selected_candidate']}; n={int(row['n_samples'])}; "
            f"AUROC={row['AUROC']:.3f}; AUPRC={row['AUPRC']:.3f}; "
            f"balanced accuracy={row['balanced_accuracy']:.3f}; ECE={row['ECE']:.3f}."
        )
    lines.extend(["", "## Interpretation", ""])
    if not external_summary.empty and float(external_summary.iloc[0]["AUROC"]) < 0.70:
        lines.append(
            "The training-only search does not reach the strict external AUROC >=0.70 target. "
            "This rules out a simple no-leakage candidate-selection fix and prioritizes either new independent melanoma tumor-tissue data or a materially different training-only representation."
        )
    else:
        lines.append("The training-only search reaches the strict external AUROC target and should be frozen before any further independent validation.")
    if not external_inner.empty:
        best = external_inner.sort_values("selection_score", ascending=False).iloc[0]
        lines.extend(
            [
                "",
                "## Discovery-only External Candidate Ranking",
                "",
                f"The top discovery-only candidate was {best['candidate']} with inner mean AUROC={best['inner_mean_AUROC']:.3f}, inner mean AUPRC={best['inner_mean_AUPRC']:.3f}, and selection score={best['selection_score']:.3f}.",
            ]
        )
    (out_dir / "TRAINING_ONLY_CANDIDATE_SEARCH_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--out", default="results/training_only_candidate_search_20260527")
    args = parser.parse_args()

    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=UserWarning)
    out_dir = ROOT / args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    X_by_cohort, _, metadata_by_cohort = load_processed_bulk(ROOT / args.processed_dir)
    lodo_metrics, lodo_predictions, lodo_inner = _evaluate_lodo(X_by_cohort, metadata_by_cohort)
    lodo_summary = _summarize_lodo_predictions(lodo_predictions)
    external_metrics, external_predictions, external_inner = _evaluate_strict_external(X_by_cohort, metadata_by_cohort)
    external_summary = _summarize_external_predictions(external_predictions)

    lodo_metrics.to_csv(out_dir / "training_only_lodo_metrics.tsv", sep="\t", index=False)
    lodo_predictions.to_csv(out_dir / "training_only_lodo_predictions.tsv", sep="\t", index=False)
    lodo_inner.to_csv(out_dir / "training_only_lodo_inner_selection.tsv", sep="\t", index=False)
    lodo_summary.to_csv(out_dir / "training_only_lodo_summary.tsv", sep="\t", index=False)
    external_metrics.to_csv(out_dir / "training_only_strict_external_metrics.tsv", sep="\t", index=False)
    external_predictions.to_csv(out_dir / "training_only_strict_external_predictions.tsv", sep="\t", index=False)
    external_inner.to_csv(out_dir / "training_only_strict_external_inner_selection.tsv", sep="\t", index=False)
    external_summary.to_csv(out_dir / "training_only_strict_external_summary.tsv", sep="\t", index=False)
    _write_audit(out_dir, lodo_summary, external_summary, external_inner)

    if external_summary.empty:
        print("Training-only search: strict external RESULT_PENDING")
    else:
        row = external_summary.iloc[0]
        print(
            "Training-only search strict external: "
            f"candidate={row['selected_candidate']}; n={int(row['n_samples'])}; "
            f"AUROC={row['AUROC']:.3f}; AUPRC={row['AUPRC']:.3f}; "
            f"balanced_accuracy={row['balanced_accuracy']:.3f}; ECE={row['ECE']:.3f}"
        )
    print(f"Wrote {out_dir}")


if __name__ == "__main__":
    main()
