from __future__ import annotations

import argparse
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
from econiche.statistics import benjamini_hochberg, paired_bootstrap_delta
from econiche_opt.model.endpoint_modules import (
    MODULE_PRIOR_WEIGHTS,
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
BIO_MODEL = "EcoNiche-Opt-BioObjectivePanelSearch"
NO_BIO_MODEL = "EcoNiche-Opt-NoBioObjectivePanelSearch"


def _active_real_cohorts(X_by_cohort: dict[str, pd.DataFrame], metadata_by_cohort: dict[str, pd.DataFrame]) -> list[str]:
    cohorts = []
    for cohort in sorted(X_by_cohort):
        if cohort.startswith("demo_cohort_"):
            continue
        meta = metadata_by_cohort.get(cohort)
        if meta is not None and "response_raw" in meta.columns and meta["response_raw"].notna().any():
            cohorts.append(cohort)
    return cohorts


def panel_weight_candidates() -> pd.DataFrame:
    modules = list(MODULE_PRIOR_WEIGHTS)
    prior = np.asarray([MODULE_PRIOR_WEIGHTS[module] for module in modules], dtype=float)
    rows: list[dict[str, object]] = []

    def add(name: str, weights: list[float] | np.ndarray, family: str) -> None:
        row: dict[str, object] = {"candidate": name, "candidate_family": family}
        for module, weight in zip(modules, weights):
            row[module] = float(weight)
        rows.append(row)

    add("bio_prior", prior, "prior")
    add("bio_prior_strong_suppression", [1.0, 0.5, 0.25, 0.5, -1.0, -1.0, 0.25], "prior_perturbation")
    add("bio_prior_no_exhaustion", [1.0, 0.5, 0.0, 0.5, -0.5, -0.5, 0.25], "prior_perturbation")
    add("response_equal_resistance_negative", [1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 1.0], "ecology_symmetric")
    add("response_equal_no_resistance", [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0], "response_only")
    add("ifn_only", [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], "single_axis")
    add("ifn_apm_cyto", [1.0, 0.7, 0.0, 0.7, 0.0, 0.0, 0.3], "response_axis")
    add("apm_cyto_resistance", [0.0, 0.7, 0.0, 1.0, -0.7, -0.7, 0.0], "response_resistance_axis")
    add("exhaustion_penalty", [1.0, 0.5, -0.5, 0.5, -0.5, -0.5, 0.25], "prior_perturbation")
    for idx, module in enumerate(modules):
        weights = prior.copy()
        weights[idx] = 0.0
        add(f"drop_{module}", weights, "single_module_drop")
    candidates = pd.DataFrame(rows)
    candidates["weight_l1"] = candidates[modules].abs().sum(axis=1)
    candidates["n_nonzero_modules"] = (candidates[modules].abs() > 0).sum(axis=1)
    return candidates


def biological_alignment_terms(weights: pd.Series) -> dict[str, float]:
    modules = list(MODULE_PRIOR_WEIGHTS)
    w = weights.reindex(modules).astype(float).to_numpy()
    prior = np.asarray([MODULE_PRIOR_WEIGHTS[module] for module in modules], dtype=float)
    if np.linalg.norm(w) > 0 and np.linalg.norm(prior) > 0:
        prior_cosine = float(np.dot(w / np.linalg.norm(w), prior / np.linalg.norm(prior)))
    else:
        prior_cosine = 0.0
    sign_agreement = float(np.mean(np.sign(w) == np.sign(prior)))
    compactness = float(np.count_nonzero(w) / len(w) <= 0.8)
    bio_bonus = 0.08 * prior_cosine + 0.04 * sign_agreement + 0.02 * compactness
    return {
        "bio_prior_cosine": prior_cosine,
        "bio_sign_agreement": sign_agreement,
        "bio_compactness_bonus": compactness,
        "bio_objective_bonus": float(bio_bonus),
    }


def _module_score(features: pd.DataFrame, weights: pd.Series) -> pd.Series:
    modules = list(MODULE_PRIOR_WEIGHTS)
    matrix = features.reindex(columns=modules).fillna(0.0).to_numpy(dtype=float)
    vector = pd.to_numeric(weights.reindex(modules), errors="coerce").fillna(0.0).to_numpy(dtype=float)
    return pd.Series(matrix @ vector, index=features.index)


def _fit_monotone_platt(score: pd.Series, y: pd.Series) -> LogisticRegression | None:
    if len(y) < 8 or y.nunique() < 2:
        return None
    model = LogisticRegression(C=1.0, class_weight="balanced", solver="lbfgs", max_iter=5000, random_state=20260527)
    model.fit(score.astype(float).to_numpy().reshape(-1, 1), y.astype(int).to_numpy())
    return model if float(model.coef_[0, 0]) > 0 else None


def _predict_prob(score: pd.Series, calibrator: LogisticRegression | None) -> pd.Series:
    if calibrator is None:
        return pd.Series(sigmoid(score), index=score.index)
    prob = calibrator.predict_proba(score.astype(float).to_numpy().reshape(-1, 1))[:, 1]
    return pd.Series(prob, index=score.index)


def _concat_series(series_by_cohort: dict[str, pd.Series], cohorts: list[str]) -> pd.Series:
    return pd.concat([series_by_cohort[cohort] for cohort in cohorts]).astype(int)


def _inner_candidate_metrics(
    candidate: pd.Series,
    train_cohorts: list[str],
    module_features: dict[str, pd.DataFrame],
    y_by_cohort: dict[str, pd.Series],
) -> pd.DataFrame:
    rows: list[dict[str, float]] = []
    for inner_holdout in train_cohorts:
        inner_train = [cohort for cohort in train_cohorts if cohort != inner_holdout]
        if not inner_train:
            continue
        y_train = _concat_series(y_by_cohort, inner_train)
        y_test = y_by_cohort[inner_holdout].astype(int)
        if y_train.nunique() < 2 or y_test.nunique() < 2:
            continue
        train_score = pd.concat(
            [_module_score(module_features[cohort].reindex(y_by_cohort[cohort].index), candidate) for cohort in inner_train]
        ).reindex(y_train.index)
        test_score = _module_score(module_features[inner_holdout].reindex(y_test.index), candidate)
        calibrator = _fit_monotone_platt(train_score, y_train)
        train_prob = _predict_prob(train_score, calibrator)
        test_prob = _predict_prob(test_score, calibrator)
        threshold = select_threshold(y_train.to_numpy(dtype=int), train_prob.to_numpy(dtype=float))
        rows.append(compute_binary_metrics(y_test, test_prob, threshold=threshold))
    return pd.DataFrame(rows)


def score_candidate(
    candidate: pd.Series,
    train_cohorts: list[str],
    module_features: dict[str, pd.DataFrame],
    y_by_cohort: dict[str, pd.Series],
    mode: str,
) -> dict[str, object]:
    metrics = _inner_candidate_metrics(candidate, train_cohorts, module_features, y_by_cohort)
    if metrics.empty:
        return {
            "candidate": candidate["candidate"],
            "candidate_family": candidate["candidate_family"],
            "inner_mean_AUROC": math.nan,
            "inner_sd_AUROC": math.nan,
            "inner_mean_AUPRC": math.nan,
            "inner_mean_balanced_accuracy": math.nan,
            "inner_mean_ECE": math.nan,
            "performance_score": -math.inf,
            "selection_score": -math.inf,
        }
    perf = float(
        metrics["AUROC"].mean()
        - 0.50 * metrics["AUROC"].std(ddof=0)
        + 0.10 * metrics["AUPRC"].mean()
        + 0.05 * metrics["balanced_accuracy"].mean()
        - 0.15 * metrics["ECE"].mean()
    )
    bio_terms = biological_alignment_terms(candidate)
    bio_bonus = float(bio_terms["bio_objective_bonus"]) if mode == "bio_objective" else 0.0
    return {
        "candidate": candidate["candidate"],
        "candidate_family": candidate["candidate_family"],
        "selection_mode": mode,
        "inner_mean_AUROC": float(metrics["AUROC"].mean()),
        "inner_sd_AUROC": float(metrics["AUROC"].std(ddof=0)),
        "inner_mean_AUPRC": float(metrics["AUPRC"].mean()),
        "inner_mean_balanced_accuracy": float(metrics["balanced_accuracy"].mean()),
        "inner_mean_ECE": float(metrics["ECE"].mean()),
        "performance_score": perf,
        **bio_terms,
        "selection_score": float(perf + bio_bonus),
    }


def select_candidate(
    candidates: pd.DataFrame,
    train_cohorts: list[str],
    module_features: dict[str, pd.DataFrame],
    y_by_cohort: dict[str, pd.Series],
    mode: str,
) -> tuple[pd.Series, pd.DataFrame]:
    rows = [score_candidate(row, train_cohorts, module_features, y_by_cohort, mode) for _, row in candidates.iterrows()]
    scores = pd.DataFrame(rows)
    finite = scores[np.isfinite(scores["selection_score"].astype(float))]
    if finite.empty:
        chosen_name = "bio_prior"
    else:
        chosen_name = str(finite.sort_values("selection_score", ascending=False).iloc[0]["candidate"])
    chosen = candidates[candidates["candidate"].astype(str).eq(chosen_name)].iloc[0]
    return chosen, scores


def _fit_evaluate_candidate(
    candidate: pd.Series,
    train_cohorts: list[str],
    test_cohort: str,
    module_features: dict[str, pd.DataFrame],
    y_by_cohort: dict[str, pd.Series],
) -> tuple[pd.Series, pd.Series, float, pd.Series | None]:
    y_train = _concat_series(y_by_cohort, train_cohorts)
    y_test = y_by_cohort[test_cohort].astype(int)
    train_score = pd.concat(
        [_module_score(module_features[cohort].reindex(y_by_cohort[cohort].index), candidate) for cohort in train_cohorts]
    ).reindex(y_train.index)
    test_score = _module_score(module_features[test_cohort].reindex(y_test.index), candidate)
    calibrator = _fit_monotone_platt(train_score, y_train)
    train_prob = _predict_prob(train_score, calibrator)
    test_prob = _predict_prob(test_score, calibrator)
    threshold = select_threshold(y_train.to_numpy(dtype=int), train_prob.to_numpy(dtype=float))
    return y_test, test_prob, float(threshold), train_prob


def _evaluation_rows(
    endpoint: str,
    stratum: str,
    holdout: str,
    train_cohorts: list[str],
    mode: str,
    candidate: pd.Series,
    y_test: pd.Series,
    prob: pd.Series,
    threshold: float,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    model_name = BIO_MODEL if mode == "bio_objective" else NO_BIO_MODEL
    metrics = compute_binary_metrics(y_test, prob, threshold=threshold)
    bio_terms = biological_alignment_terms(candidate)
    metric_row = {
        **metrics,
        "endpoint": endpoint,
        "stratum": stratum,
        "cohort": holdout,
        "model_name": model_name,
        "selection_mode": mode,
        "selected_candidate": candidate["candidate"],
        "candidate_family": candidate["candidate_family"],
        "threshold": threshold,
        "train_cohorts": ",".join(train_cohorts),
        "n_samples": int(len(y_test)),
        "n_responders": int(y_test.sum()),
        "n_nonresponders": int((y_test == 0).sum()),
        **bio_terms,
    }
    prediction_rows = [
        {
            "endpoint": endpoint,
            "stratum": stratum,
            "cohort": holdout,
            "model_name": model_name,
            "selection_mode": mode,
            "selected_candidate": candidate["candidate"],
            "sample_id": sample_id,
            "true_response_label": int(y_test.loc[sample_id]),
            "response_probability": float(prob.loc[sample_id]),
            "threshold": threshold,
            "pred_response_label": int(prob.loc[sample_id] >= threshold),
        }
        for sample_id in y_test.index
    ]
    return metric_row, prediction_rows


def evaluate_lodo(
    X_by_cohort: dict[str, pd.DataFrame],
    metadata_by_cohort: dict[str, pd.DataFrame],
    candidates: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    active = _active_real_cohorts(X_by_cohort, metadata_by_cohort)
    endpoint_data = prepare_endpoint_data(X_by_cohort, metadata_by_cohort, active, PRIMARY_ENDPOINT)
    module_features, _ = build_module_features_by_cohort(endpoint_data.X_by_cohort)
    strata = default_strata(endpoint_data.X_by_cohort)
    metrics_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    selection_rows: list[dict[str, object]] = []
    for stratum in PRIMARY_STRATA:
        spec = strata.get(stratum, {})
        train_pool = [cohort for cohort in spec.get("train_pool", []) if cohort in endpoint_data.y_response_by_cohort]
        holdouts = [cohort for cohort in spec.get("holdouts", []) if cohort in endpoint_data.y_response_by_cohort]
        for holdout in holdouts:
            train_cohorts = [cohort for cohort in train_pool if cohort != holdout and cohort in endpoint_data.y_response_by_cohort]
            if len(train_cohorts) < 2:
                continue
            y_train = _concat_series(endpoint_data.y_response_by_cohort, train_cohorts)
            y_test = endpoint_data.y_response_by_cohort[holdout].astype(int)
            if y_train.nunique() < 2 or y_test.nunique() < 2:
                continue
            for mode in ["bio_objective", "no_bio_objective"]:
                chosen, scores = select_candidate(candidates, train_cohorts, module_features, endpoint_data.y_response_by_cohort, mode)
                scores.insert(0, "endpoint", PRIMARY_ENDPOINT)
                scores.insert(1, "stratum", stratum)
                scores.insert(2, "holdout", holdout)
                scores.insert(3, "train_cohorts", ",".join(train_cohorts))
                scores["selected_candidate"] = chosen["candidate"]
                selection_rows.extend(scores.to_dict("records"))
                y_eval, prob, threshold, _ = _fit_evaluate_candidate(
                    chosen, train_cohorts, holdout, module_features, endpoint_data.y_response_by_cohort
                )
                metric_row, pred_rows = _evaluation_rows(
                    PRIMARY_ENDPOINT,
                    stratum,
                    holdout,
                    train_cohorts,
                    mode,
                    chosen,
                    y_eval,
                    prob,
                    threshold,
                )
                metrics_rows.append(metric_row)
                prediction_rows.extend(pred_rows)
    return pd.DataFrame(metrics_rows), pd.DataFrame(prediction_rows), pd.DataFrame(selection_rows)


def evaluate_strict_external(
    X_by_cohort: dict[str, pd.DataFrame],
    metadata_by_cohort: dict[str, pd.DataFrame],
    candidates: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cohorts = [*DISCOVERY_COHORTS, *STRICT_EXTERNAL_COHORTS]
    endpoint_data = prepare_endpoint_data(X_by_cohort, metadata_by_cohort, cohorts, STRICT_ENDPOINT)
    module_features, _ = build_module_features_by_cohort(endpoint_data.X_by_cohort)
    train_cohorts = [cohort for cohort in DISCOVERY_COHORTS if cohort in endpoint_data.y_response_by_cohort]
    metrics_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    selection_rows: list[dict[str, object]] = []
    for mode in ["bio_objective", "no_bio_objective"]:
        chosen, scores = select_candidate(candidates, train_cohorts, module_features, endpoint_data.y_response_by_cohort, mode)
        scores.insert(0, "endpoint", STRICT_ENDPOINT)
        scores.insert(1, "stratum", "strict_melanoma_pd1_like_external")
        scores.insert(2, "holdout", "+".join(STRICT_EXTERNAL_COHORTS))
        scores.insert(3, "train_cohorts", ",".join(train_cohorts))
        scores["selected_candidate"] = chosen["candidate"]
        selection_rows.extend(scores.to_dict("records"))
        for holdout in STRICT_EXTERNAL_COHORTS:
            if holdout not in endpoint_data.y_response_by_cohort or holdout not in module_features:
                continue
            y_test = endpoint_data.y_response_by_cohort[holdout].astype(int)
            if y_test.nunique() < 2:
                continue
            y_eval, prob, threshold, _ = _fit_evaluate_candidate(
                chosen, train_cohorts, holdout, module_features, endpoint_data.y_response_by_cohort
            )
            metric_row, pred_rows = _evaluation_rows(
                STRICT_ENDPOINT,
                "strict_melanoma_pd1_like_external",
                holdout,
                train_cohorts,
                mode,
                chosen,
                y_eval,
                prob,
                threshold,
            )
            metric_row["selection_boundary"] = "discovery_only_inner_lodo_no_external_selection"
            metrics_rows.append(metric_row)
            prediction_rows.extend(pred_rows)
    return pd.DataFrame(metrics_rows), pd.DataFrame(prediction_rows), pd.DataFrame(selection_rows)


def summarize_predictions(predictions: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if predictions.empty:
        return pd.DataFrame()
    for (endpoint, stratum, model_name, mode), frame in predictions.groupby(["endpoint", "stratum", "model_name", "selection_mode"]):
        y = frame["true_response_label"].astype(int)
        prob = frame["response_probability"].astype(float)
        metrics = compute_binary_metrics(y, prob)
        fold_metrics = []
        for _, fold in frame.groupby("cohort"):
            if fold["true_response_label"].nunique() == 2:
                fold_metrics.append(compute_binary_metrics(fold["true_response_label"].astype(int), fold["response_probability"].astype(float)))
        fold_frame = pd.DataFrame(fold_metrics)
        rows.append(
            {
                "endpoint": endpoint,
                "stratum": stratum,
                "model_name": model_name,
                "selection_mode": mode,
                "n_samples": int(len(frame)),
                "n_cohorts": int(frame["cohort"].nunique()),
                "n_responders": int(y.sum()),
                "n_nonresponders": int((y == 0).sum()),
                "pooled_AUROC": metrics["AUROC"],
                "pooled_AUPRC": metrics["AUPRC"],
                "pooled_balanced_accuracy": metrics["balanced_accuracy"],
                "pooled_Brier": metrics["Brier"],
                "pooled_ECE": metrics["ECE"],
                "fold_AUROC_sd": float(fold_frame["AUROC"].std(ddof=0)) if not fold_frame.empty else math.nan,
                "selected_candidates": ",".join(sorted(frame["selected_candidate"].astype(str).unique())),
            }
        )
    return pd.DataFrame(rows).sort_values(["endpoint", "stratum", "selection_mode"])


def paired_comparison(predictions: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if predictions.empty:
        return pd.DataFrame()
    for (endpoint, stratum), frame in predictions.groupby(["endpoint", "stratum"]):
        bio = frame[frame["model_name"].eq(BIO_MODEL)]
        no_bio = frame[frame["model_name"].eq(NO_BIO_MODEL)]
        merged = bio.merge(
            no_bio[["cohort", "sample_id", "true_response_label", "response_probability"]],
            on=["cohort", "sample_id", "true_response_label"],
            suffixes=("_bio", "_no_bio"),
        )
        if len(merged) < 8 or merged["true_response_label"].nunique() < 2:
            continue
        y = merged["true_response_label"].astype(int)
        bio_prob = merged["response_probability_bio"].astype(float)
        no_bio_prob = merged["response_probability_no_bio"].astype(float)
        bio_metrics = compute_binary_metrics(y, bio_prob)
        no_bio_metrics = compute_binary_metrics(y, no_bio_prob)
        fold_deltas: list[float] = []
        for _, fold in merged.groupby("cohort"):
            if fold["true_response_label"].nunique() == 2:
                bm = compute_binary_metrics(fold["true_response_label"].astype(int), fold["response_probability_bio"].astype(float))
                nm = compute_binary_metrics(fold["true_response_label"].astype(int), fold["response_probability_no_bio"].astype(float))
                fold_deltas.append(float(bm["AUROC"] - nm["AUROC"]))
        stats = paired_bootstrap_delta(y, bio_prob, no_bio_prob, n_bootstrap=1000, random_state=20260527)
        rows.append(
            {
                "endpoint": endpoint,
                "stratum": stratum,
                "target_model": BIO_MODEL,
                "ablation_model": NO_BIO_MODEL,
                "n_samples": int(len(merged)),
                "target_AUROC": bio_metrics["AUROC"],
                "ablation_AUROC": no_bio_metrics["AUROC"],
                "delta_AUROC": float(bio_metrics["AUROC"] - no_bio_metrics["AUROC"]),
                "target_AUPRC": bio_metrics["AUPRC"],
                "ablation_AUPRC": no_bio_metrics["AUPRC"],
                "delta_AUPRC": float(bio_metrics["AUPRC"] - no_bio_metrics["AUPRC"]),
                "target_balanced_accuracy": bio_metrics["balanced_accuracy"],
                "ablation_balanced_accuracy": no_bio_metrics["balanced_accuracy"],
                "delta_balanced_accuracy": float(bio_metrics["balanced_accuracy"] - no_bio_metrics["balanced_accuracy"]),
                "target_ECE": bio_metrics["ECE"],
                "ablation_ECE": no_bio_metrics["ECE"],
                "delta_ECE": float(bio_metrics["ECE"] - no_bio_metrics["ECE"]),
                "fold_delta_AUROC_sd": float(np.std(fold_deltas, ddof=0)) if fold_deltas else math.nan,
                **stats,
            }
        )
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result["fdr_q"] = benjamini_hochberg(result["p_value"].fillna(1.0))
    result["claim_level"] = np.select(
        [
            (result["delta_AUROC"] > 0) & (result["fdr_q"] <= 0.05),
            result["delta_AUROC"] > 0,
            (result["delta_balanced_accuracy"] > 0) | (result["delta_ECE"] < 0),
        ],
        [
            "FDR_supported_biological_objective_gain",
            "point_estimate_biological_objective_gain",
            "biological_objective_calibration_or_threshold_tradeoff",
        ],
        default="biological_objective_not_supported_for_this_context",
    )
    return result


def write_audit(out_dir: Path, summary: pd.DataFrame, comparison: pd.DataFrame) -> None:
    lines = [
        "# Aligned Biological Objective Audit",
        "",
        "This audit tests whether adding a biological-prior term to training-only panel-weight selection changes primary melanoma and strict external performance. Candidate selection, thresholding, and calibration use training or discovery cohorts only; locked external labels are used only after the selected candidate and threshold policy are fixed.",
        "",
        "## Model Summary",
        "",
    ]
    for _, row in summary.iterrows():
        lines.append(
            f"- {row['endpoint']} / {row['stratum']} / {row['model_name']}: n={int(row['n_samples'])}, "
            f"AUROC={float(row['pooled_AUROC']):.3f}, AUPRC={float(row['pooled_AUPRC']):.3f}, "
            f"BA={float(row['pooled_balanced_accuracy']):.3f}, ECE={float(row['pooled_ECE']):.3f}; "
            f"selected={row['selected_candidates']}."
        )
    lines.extend(["", "## Paired Bio Objective Comparison", ""])
    for _, row in comparison.iterrows():
        lines.append(
            f"- {row['endpoint']} / {row['stratum']}: delta AUROC={float(row['delta_AUROC']):.3f}, "
            f"delta BA={float(row['delta_balanced_accuracy']):.3f}, delta ECE={float(row['delta_ECE']):.3f}, "
            f"95% bootstrap CI [{float(row['ci_low']):.3f}, {float(row['ci_high']):.3f}], "
            f"FDR q={float(row['fdr_q']):.3f}; {row['claim_level']}."
        )
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "Use biological-objective performance language only for contexts with positive delta AUROC and FDR support. Otherwise, restrict claims to point-estimate, calibration, threshold-operation, or stability tradeoffs reflected in the table.",
        ]
    )
    (out_dir / "ALIGNED_BIOLOGICAL_OBJECTIVE_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")


def run(processed_dir: Path, out_dir: Path) -> None:
    X_by_cohort, _, metadata_by_cohort = load_processed_bulk(processed_dir)
    candidates = panel_weight_candidates()
    lodo_metrics, lodo_predictions, lodo_selection = evaluate_lodo(X_by_cohort, metadata_by_cohort, candidates)
    external_metrics, external_predictions, external_selection = evaluate_strict_external(X_by_cohort, metadata_by_cohort, candidates)
    predictions = pd.concat([lodo_predictions, external_predictions], ignore_index=True, sort=False)
    metrics = pd.concat([lodo_metrics, external_metrics], ignore_index=True, sort=False)
    selection = pd.concat([lodo_selection, external_selection], ignore_index=True, sort=False)
    summary = summarize_predictions(predictions)
    comparison = paired_comparison(predictions)

    out_dir.mkdir(parents=True, exist_ok=True)
    candidates.to_csv(out_dir / "aligned_biological_objective_candidates.tsv", sep="\t", index=False)
    metrics.to_csv(out_dir / "aligned_biological_objective_metrics.tsv", sep="\t", index=False)
    predictions.to_csv(out_dir / "aligned_biological_objective_predictions.tsv", sep="\t", index=False)
    selection.to_csv(out_dir / "aligned_biological_objective_inner_selection.tsv", sep="\t", index=False)
    summary.to_csv(out_dir / "aligned_biological_objective_summary.tsv", sep="\t", index=False)
    comparison.to_csv(out_dir / "aligned_biological_objective_comparison.tsv", sep="\t", index=False)
    write_audit(out_dir, summary, comparison)
    print(f"Wrote aligned biological objective audit to {out_dir}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--out", default="results/aligned_biological_objective_audit_20260527")
    args = parser.parse_args()
    run(ROOT / args.processed_dir, ROOT / args.out)


if __name__ == "__main__":
    main()
