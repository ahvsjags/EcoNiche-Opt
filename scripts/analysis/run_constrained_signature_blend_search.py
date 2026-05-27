from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.io import load_processed_bulk
from econiche.metrics import compute_binary_metrics
from econiche_opt.model.endpoint_modules import (
    STRONG_BASELINES,
    build_fixed_scores_by_cohort,
    build_module_features_by_cohort,
    default_strata,
    endpoint_label_series,
    prepare_endpoint_data,
    select_threshold,
    sigmoid,
)


PRIMARY_ENDPOINT = "primary_recist"
STRICT_ENDPOINT = "strict_recist"
PRIMARY_STRATA = ["melanoma_core_high_evidence", "melanoma_recist_supported_primary"]
DISCOVERY_COHORTS = ["GSE91061", "GSE78220", "PRJEB23709_PD1_PRE"]
STRICT_EXTERNAL_COHORTS = ["GSE145996", "PHS000452_LIU_LIKE_PRE"]
CBIO_EXTERNAL_COHORTS = ["GSE145996", "CBIO_LIU_DFCI_2019_PRE"]
TARGET_BLEND_MODEL = "EcoNiche-Opt-ConstrainedBlend"
BASE_FEATURES = [
    "EcoNiche-Opt-ModulePriorFixed",
    "IFNG",
    "CXCL9",
    "TIG",
    "TIDE_dysfunction",
    "TIDE_exclusion",
    "CYT",
    "APM",
    "IPRES",
]


def merge_processed_dirs(*processed_dirs: Path) -> tuple[dict[str, pd.DataFrame], dict[str, pd.Series], dict[str, pd.DataFrame]]:
    merged_X: dict[str, pd.DataFrame] = {}
    merged_y: dict[str, pd.Series] = {}
    merged_meta: dict[str, pd.DataFrame] = {}
    for processed_dir in processed_dirs:
        if not processed_dir.exists():
            continue
        X, y, meta = load_processed_bulk(processed_dir)
        merged_X.update(X)
        merged_y.update(y)
        merged_meta.update(meta)
    return merged_X, merged_y, merged_meta


def candidate_weight_specs() -> list[dict[str, object]]:
    specs: list[dict[str, object]] = []

    def add(name: str, weights: dict[str, float]) -> None:
        total = float(sum(max(0.0, value) for value in weights.values()))
        if total <= 0:
            return
        normalized = {key: float(value) / total for key, value in weights.items() if float(value) > 0}
        specs.append({"candidate": name, "weights": normalized})

    for feature in BASE_FEATURES:
        add(f"single__{feature}", {feature: 1.0})
    immune_family = ["IFNG", "CXCL9", "TIG", "TIDE_dysfunction", "CYT", "APM"]
    eight_family = ["IFNG", "CXCL9", "TIG", "TIDE_dysfunction", "TIDE_exclusion", "CYT", "APM", "IPRES"]
    add("immune_family_mean", {feature: 1.0 for feature in immune_family})
    add("eight_signature_family_mean", {feature: 1.0 for feature in eight_family})
    for alpha in [0.50, 0.65, 0.75, 0.85]:
        add(
            f"module_plus_immune_family_alpha{alpha:g}",
            {"EcoNiche-Opt-ModulePriorFixed": alpha, **{feature: (1 - alpha) / len(immune_family) for feature in immune_family}},
        )
        add(
            f"module_plus_eight_family_alpha{alpha:g}",
            {"EcoNiche-Opt-ModulePriorFixed": alpha, **{feature: (1 - alpha) / len(eight_family) for feature in eight_family}},
        )
    for feature in ["IFNG", "CXCL9", "TIG", "TIDE_dysfunction", "CYT", "APM", "IPRES", "TIDE_exclusion"]:
        for alpha in [0.50, 0.65, 0.80]:
            add(f"module_plus_{feature}_alpha{alpha:g}", {"EcoNiche-Opt-ModulePriorFixed": alpha, feature: 1 - alpha})
    pair_features = ["CXCL9", "TIG", "TIDE_dysfunction", "CYT", "APM", "IPRES"]
    for i, left in enumerate(pair_features):
        for right in pair_features[i + 1 :]:
            add(
                f"module_plus_{left}_{right}",
                {"EcoNiche-Opt-ModulePriorFixed": 0.60, left: 0.20, right: 0.20},
            )
    # A small set of biology-driven blends that remain interpretable.
    add("ifn_apm_effector_resistance_balance", {"IFNG": 0.25, "CXCL9": 0.20, "APM": 0.25, "CYT": 0.15, "IPRES": 0.15})
    add(
        "module_ifn_apm_cyt",
        {"EcoNiche-Opt-ModulePriorFixed": 0.55, "IFNG": 0.15, "CXCL9": 0.10, "APM": 0.10, "CYT": 0.10},
    )
    dedup: dict[str, dict[str, object]] = {}
    for spec in specs:
        key = tuple(sorted((feature, round(weight, 6)) for feature, weight in spec["weights"].items()))
        dedup[str(key)] = spec
    return list(dedup.values())


def build_score_tables(X_by_cohort: dict[str, pd.DataFrame]) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    module_features, coverage = build_module_features_by_cohort(X_by_cohort)
    fixed = build_fixed_scores_by_cohort(
        X_by_cohort,
        module_features,
        baselines=sorted(set([*STRONG_BASELINES, *BASE_FEATURES])),
    )
    tables: dict[str, pd.DataFrame] = {}
    for cohort, scores in fixed.items():
        frame = pd.DataFrame({feature: scores[feature] for feature in BASE_FEATURES if feature in scores})
        tables[cohort] = frame.fillna(0.0)
    return tables, coverage


def blend_score(score_table: pd.DataFrame, weights: dict[str, float]) -> pd.Series:
    score = pd.Series(0.0, index=score_table.index)
    for feature, weight in weights.items():
        if feature in score_table.columns:
            score = score + float(weight) * score_table[feature].astype(float)
    return score.fillna(0.0)


def fit_monotone_platt(score: pd.Series, y: pd.Series) -> LogisticRegression | None:
    common = score.index.intersection(y.index)
    if len(common) < 8 or y.loc[common].nunique() < 2:
        return None
    model = LogisticRegression(C=1.0, class_weight="balanced", solver="lbfgs", max_iter=5000, random_state=20260527)
    model.fit(score.loc[common].astype(float).to_numpy().reshape(-1, 1), y.loc[common].astype(int).to_numpy())
    return model if float(model.coef_[0, 0]) > 0 else None


def predict_prob(score: pd.Series, calibrator: LogisticRegression | None) -> pd.Series:
    if calibrator is None:
        return pd.Series(sigmoid(score), index=score.index)
    values = calibrator.predict_proba(score.astype(float).to_numpy().reshape(-1, 1))[:, 1]
    return pd.Series(values, index=score.index)


def concat_series(series_by_cohort: dict[str, pd.Series], cohorts: list[str]) -> pd.Series:
    return pd.concat([series_by_cohort[cohort] for cohort in cohorts], axis=0)


def score_inner_candidate(
    spec: dict[str, object],
    train_cohorts: list[str],
    score_tables: dict[str, pd.DataFrame],
    y_by_cohort: dict[str, pd.Series],
) -> dict[str, object]:
    rows = []
    if len(train_cohorts) < 2:
        return {"candidate": spec["candidate"], "selection_score": -math.inf}
    for holdout in train_cohorts:
        inner_train = [cohort for cohort in train_cohorts if cohort != holdout]
        if not inner_train:
            continue
        y_train = pd.concat([y_by_cohort[cohort] for cohort in inner_train]).astype(int)
        y_test = y_by_cohort[holdout].astype(int)
        train_score = pd.concat(
            [blend_score(score_tables[cohort].reindex(y_by_cohort[cohort].index), spec["weights"]) for cohort in inner_train]
        ).reindex(y_train.index)
        test_score = blend_score(score_tables[holdout].reindex(y_test.index), spec["weights"])
        calibrator = fit_monotone_platt(train_score, y_train)
        train_prob = predict_prob(train_score, calibrator)
        test_prob = predict_prob(test_score, calibrator)
        threshold = select_threshold(y_train.to_numpy(dtype=int), train_prob.to_numpy(dtype=float))
        rows.append(compute_binary_metrics(y_test, test_prob, threshold=threshold))
    if not rows:
        return {"candidate": spec["candidate"], "selection_score": -math.inf}
    frame = pd.DataFrame(rows)
    selection_score = (
        frame["AUROC"].mean()
        - 0.45 * frame["AUROC"].std(ddof=0)
        + 0.10 * frame["AUPRC"].mean()
        + 0.10 * frame["balanced_accuracy"].mean()
        - 0.15 * frame["ECE"].mean()
    )
    return {
        "candidate": spec["candidate"],
        "inner_mean_AUROC": float(frame["AUROC"].mean()),
        "inner_sd_AUROC": float(frame["AUROC"].std(ddof=0)),
        "inner_mean_AUPRC": float(frame["AUPRC"].mean()),
        "inner_mean_balanced_accuracy": float(frame["balanced_accuracy"].mean()),
        "inner_mean_ECE": float(frame["ECE"].mean()),
        "selection_score": float(selection_score),
    }


def best_spec(
    specs: list[dict[str, object]],
    train_cohorts: list[str],
    score_tables: dict[str, pd.DataFrame],
    y_by_cohort: dict[str, pd.Series],
) -> tuple[dict[str, object], pd.DataFrame]:
    rows = [score_inner_candidate(spec, train_cohorts, score_tables, y_by_cohort) for spec in specs]
    selection = pd.DataFrame(rows).sort_values("selection_score", ascending=False).reset_index(drop=True)
    valid = selection[np.isfinite(selection["selection_score"].astype(float))]
    best_name = valid.iloc[0]["candidate"] if not valid.empty else "single__EcoNiche-Opt-ModulePriorFixed"
    spec = next(candidate for candidate in specs if candidate["candidate"] == best_name)
    return spec, selection


def predict_holdout(
    spec: dict[str, object],
    train_cohorts: list[str],
    holdout: str,
    score_tables: dict[str, pd.DataFrame],
    y_by_cohort: dict[str, pd.Series],
) -> tuple[pd.Series, pd.Series, float, LogisticRegression | None]:
    y_train = pd.concat([y_by_cohort[cohort] for cohort in train_cohorts]).astype(int)
    y_test = y_by_cohort[holdout].astype(int)
    train_score = pd.concat(
        [blend_score(score_tables[cohort].reindex(y_by_cohort[cohort].index), spec["weights"]) for cohort in train_cohorts]
    ).reindex(y_train.index)
    test_score = blend_score(score_tables[holdout].reindex(y_test.index), spec["weights"])
    calibrator = fit_monotone_platt(train_score, y_train)
    train_prob = predict_prob(train_score, calibrator)
    test_prob = predict_prob(test_score, calibrator)
    threshold = select_threshold(y_train.to_numpy(dtype=int), train_prob.to_numpy(dtype=float))
    return y_test, test_prob, float(threshold), calibrator


def evaluate_primary_lodo(
    X_by_cohort: dict[str, pd.DataFrame],
    metadata_by_cohort: dict[str, pd.DataFrame],
    score_tables: dict[str, pd.DataFrame],
    specs: list[dict[str, object]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    active = [cohort for cohort in X_by_cohort if not cohort.startswith("demo_cohort_")]
    strata = default_strata(active)
    metric_rows = []
    prediction_rows = []
    selection_rows = []
    for stratum in PRIMARY_STRATA:
        endpoint_data = prepare_endpoint_data(X_by_cohort, metadata_by_cohort, strata[stratum]["cohorts"], PRIMARY_ENDPOINT)
        y_by_cohort = endpoint_data.y_response_by_cohort
        for holdout in strata[stratum]["holdouts"]:
            if holdout not in y_by_cohort:
                continue
            train_cohorts = [cohort for cohort in strata[stratum]["train_pool"] if cohort != holdout and cohort in y_by_cohort]
            if len(train_cohorts) < 2:
                continue
            best, selection = best_spec(specs, train_cohorts, score_tables, y_by_cohort)
            selection["endpoint"] = PRIMARY_ENDPOINT
            selection["stratum"] = stratum
            selection["holdout"] = holdout
            selection_rows.append(selection)
            y_test, prob, threshold, _ = predict_holdout(best, train_cohorts, holdout, score_tables, y_by_cohort)
            metrics = compute_binary_metrics(y_test, prob, threshold=threshold)
            metric_rows.append(
                {
                    **metrics,
                    "endpoint": PRIMARY_ENDPOINT,
                    "stratum": stratum,
                    "cohort": holdout,
                    "model_name": TARGET_BLEND_MODEL,
                    "selected_candidate": best["candidate"],
                    "weights_json": json.dumps(best["weights"], sort_keys=True),
                    "n_samples": int(len(y_test)),
                    "n_responders": int(y_test.sum()),
                    "n_nonresponders": int((y_test == 0).sum()),
                    "threshold": threshold,
                    "selection_boundary": "inner_lodo_training_only",
                }
            )
            for sample_id in y_test.index:
                prediction_rows.append(
                    {
                        "endpoint": PRIMARY_ENDPOINT,
                        "stratum": stratum,
                        "cohort": holdout,
                        "model_name": TARGET_BLEND_MODEL,
                        "selected_candidate": best["candidate"],
                        "sample_id": sample_id,
                        "true_response_label": int(y_test.loc[sample_id]),
                        "response_probability": float(prob.loc[sample_id]),
                        "threshold": threshold,
                    }
                )
    return pd.DataFrame(metric_rows), pd.DataFrame(prediction_rows), pd.concat(selection_rows, ignore_index=True) if selection_rows else pd.DataFrame()


def endpoint_labels(metadata_by_cohort: dict[str, pd.DataFrame], cohorts: list[str], endpoint: str) -> dict[str, pd.Series]:
    labels: dict[str, pd.Series] = {}
    for cohort in cohorts:
        meta = metadata_by_cohort.get(cohort)
        if meta is None or "response_raw" not in meta.columns:
            continue
        y = endpoint_label_series(meta["response_raw"], endpoint).dropna().astype(int)
        if len(y) >= 4 and y.nunique() == 2:
            labels[cohort] = y
    return labels


def evaluate_external_group(
    group_id: str,
    cohorts: list[str],
    endpoint: str,
    metadata_by_cohort: dict[str, pd.DataFrame],
    score_tables: dict[str, pd.DataFrame],
    specs: list[dict[str, object]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    y_by_cohort = endpoint_labels(metadata_by_cohort, sorted(set(DISCOVERY_COHORTS + cohorts)), endpoint)
    train_cohorts = [cohort for cohort in DISCOVERY_COHORTS if cohort in y_by_cohort]
    best, selection = best_spec(specs, train_cohorts, score_tables, y_by_cohort)
    selection["endpoint"] = endpoint
    selection["stratum"] = group_id
    selection["holdout"] = "+".join(cohorts)
    metric_rows = []
    prediction_rows = []
    for cohort in cohorts:
        if cohort not in y_by_cohort:
            continue
        y_test, prob, threshold, _ = predict_holdout(best, train_cohorts, cohort, score_tables, y_by_cohort)
        metrics = compute_binary_metrics(y_test, prob, threshold=threshold)
        metric_rows.append(
            {
                **metrics,
                "endpoint": endpoint,
                "stratum": group_id,
                "cohort": cohort,
                "model_name": TARGET_BLEND_MODEL,
                "selected_candidate": best["candidate"],
                "weights_json": json.dumps(best["weights"], sort_keys=True),
                "n_samples": int(len(y_test)),
                "n_responders": int(y_test.sum()),
                "n_nonresponders": int((y_test == 0).sum()),
                "threshold": threshold,
                "selection_boundary": "discovery_only_inner_lodo_no_external_selection",
            }
        )
        for sample_id in y_test.index:
            prediction_rows.append(
                {
                    "endpoint": endpoint,
                    "stratum": group_id,
                    "cohort": cohort,
                    "model_name": TARGET_BLEND_MODEL,
                    "selected_candidate": best["candidate"],
                    "sample_id": sample_id,
                    "true_response_label": int(y_test.loc[sample_id]),
                    "response_probability": float(prob.loc[sample_id]),
                    "threshold": threshold,
                }
            )
    return pd.DataFrame(metric_rows), pd.DataFrame(prediction_rows), selection


def summarize_predictions(predictions: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if predictions.empty:
        return pd.DataFrame()
    for (endpoint, stratum, model_name), frame in predictions.groupby(["endpoint", "stratum", "model_name"]):
        y = frame["true_response_label"].astype(int)
        p = frame["response_probability"].astype(float)
        if y.nunique() < 2:
            continue
        metrics = compute_binary_metrics(y, p)
        rows.append(
            {
                **metrics,
                "endpoint": endpoint,
                "stratum": stratum,
                "model_name": model_name,
                "n_samples": int(len(frame)),
                "n_cohorts": int(frame["cohort"].nunique()),
                "n_responders": int(y.sum()),
                "n_nonresponders": int((y == 0).sum()),
                "response_prevalence": float(y.mean()),
                "AUPRC_minus_prevalence": float(metrics["AUPRC"] - y.mean()),
                "selected_candidates": ",".join(sorted(frame["selected_candidate"].astype(str).unique())),
            }
        )
    return pd.DataFrame(rows).sort_values(["endpoint", "stratum"]).reset_index(drop=True)


def write_audit(out_dir: Path, summary: pd.DataFrame, n_candidates: int) -> None:
    lines = [
        "# Constrained Signature Blend Search Audit",
        "",
        f"Predeclared constrained blend candidates tested: {n_candidates}. Candidate selection uses only inner LODO training folds or discovery-only inner LODO for external groups.",
        "",
    ]
    if summary.empty:
        lines.append("No summary rows were produced.")
    else:
        for _, row in summary.iterrows():
            lines.append(
                f"- {row['endpoint']} / {row['stratum']}: AUROC={row['AUROC']:.3f}, AUPRC={row['AUPRC']:.3f}, "
                f"balanced accuracy={row['balanced_accuracy']:.3f}, ECE={row['ECE']:.3f}, selected={row['selected_candidates']}."
            )
    (out_dir / "CONSTRAINED_BLEND_SEARCH_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--cbio-dir", default="data/processed/cbioportal_melanoma")
    parser.add_argument("--out", default="results/constrained_signature_blend_search_20260527")
    args = parser.parse_args()

    out_dir = ROOT / args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    X_by_cohort, _, metadata_by_cohort = merge_processed_dirs(ROOT / args.processed_dir, ROOT / args.cbio_dir)
    score_tables, coverage = build_score_tables(X_by_cohort)
    specs = candidate_weight_specs()
    pd.DataFrame(
        [{"candidate": spec["candidate"], "weights_json": json.dumps(spec["weights"], sort_keys=True)} for spec in specs]
    ).to_csv(out_dir / "constrained_blend_candidates.tsv", sep="\t", index=False)
    coverage.to_csv(out_dir / "constrained_blend_module_coverage.tsv", sep="\t", index=False)

    primary_metrics, primary_predictions, primary_selection = evaluate_primary_lodo(X_by_cohort, metadata_by_cohort, score_tables, specs)
    external_results = []
    external_predictions = []
    external_selection = []
    for group_id, cohorts in {
        "strict_melanoma_pd1_like_external": STRICT_EXTERNAL_COHORTS,
        "strict_cbio_liu_plus_gse145996": CBIO_EXTERNAL_COHORTS,
    }.items():
        metrics, predictions, selection = evaluate_external_group(group_id, cohorts, STRICT_ENDPOINT, metadata_by_cohort, score_tables, specs)
        external_results.append(metrics)
        external_predictions.append(predictions)
        external_selection.append(selection)
    metrics = pd.concat([primary_metrics, *external_results], ignore_index=True)
    predictions = pd.concat([primary_predictions, *external_predictions], ignore_index=True)
    selection = pd.concat([primary_selection, *external_selection], ignore_index=True)
    summary = summarize_predictions(predictions)
    metrics.to_csv(out_dir / "constrained_blend_metrics.tsv", sep="\t", index=False)
    predictions.to_csv(out_dir / "constrained_blend_predictions.tsv", sep="\t", index=False)
    selection.to_csv(out_dir / "constrained_blend_inner_selection.tsv", sep="\t", index=False)
    summary.to_csv(out_dir / "constrained_blend_summary.tsv", sep="\t", index=False)
    write_audit(out_dir, summary, len(specs))
    print(f"Wrote constrained blend search outputs to {out_dir}")


if __name__ == "__main__":
    main()
