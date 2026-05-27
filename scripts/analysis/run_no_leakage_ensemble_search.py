from __future__ import annotations

import argparse
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.io import load_processed_bulk
from econiche.metrics import compute_binary_metrics
from econiche.statistics import benjamini_hochberg, paired_bootstrap_delta
from econiche_opt.model.endpoint_modules import (
    STRONG_BASELINES,
    build_fixed_scores_by_cohort,
    build_module_features_by_cohort,
    default_strata,
    prepare_endpoint_data,
    select_threshold,
)

from scripts.analysis.run_aligned_biological_objective_audit import (
    _module_score,
    panel_weight_candidates,
)
from scripts.analysis.run_aligned_interaction_edge_audit import build_edge_features


PRIMARY_ENDPOINT = "primary_recist"
STRICT_ENDPOINT = "strict_recist"
PRIMARY_STRATA = ["melanoma_core_high_evidence", "melanoma_recist_supported_primary"]
DISCOVERY_COHORTS = ["GSE91061", "GSE78220", "PRJEB23709_PD1_PRE"]
STRICT_EXTERNAL_COHORTS = ["GSE145996", "PHS000452_LIU_LIKE_PRE"]
ENSEMBLE_MODEL = "EcoNiche-Opt-NoLeakageStackedEnsemble"
BASELINE_MODEL = "EcoNiche-Opt-ModulePriorFixed"
EIGHT_SIGNATURES = ["IFNG", "CXCL9", "TIG", "TIDE_dysfunction", "TIDE_exclusion", "CYT", "APM", "IPRES"]


def _active_real_cohorts(X_by_cohort: dict[str, pd.DataFrame], metadata_by_cohort: dict[str, pd.DataFrame]) -> list[str]:
    cohorts = []
    for cohort in sorted(X_by_cohort):
        if cohort.startswith("demo_cohort_"):
            continue
        meta = metadata_by_cohort.get(cohort)
        if meta is not None and "response_raw" in meta.columns and meta["response_raw"].notna().any():
            cohorts.append(cohort)
    return cohorts


def build_feature_matrix_by_cohort(
    X_by_cohort: dict[str, pd.DataFrame],
    module_features_by_cohort: dict[str, pd.DataFrame],
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    fixed_scores = build_fixed_scores_by_cohort(
        X_by_cohort,
        module_features_by_cohort,
        baselines=sorted(set([*STRONG_BASELINES, *EIGHT_SIGNATURES])),
    )
    bio_candidates = panel_weight_candidates()
    feature_by_cohort: dict[str, pd.DataFrame] = {}
    coverage_rows: list[dict[str, object]] = []
    for cohort, module_features in module_features_by_cohort.items():
        frames = []
        module_frame = module_features.add_prefix("module__")
        frames.append(module_frame)
        fixed_frame = pd.DataFrame(
            {f"signature__{name}": score.reindex(module_features.index).fillna(0.0) for name, score in fixed_scores.get(cohort, {}).items()},
            index=module_features.index,
        )
        frames.append(fixed_frame)
        bio_frame = pd.DataFrame(
            {
                f"bio_candidate__{row['candidate']}": _module_score(module_features, row).reindex(module_features.index).fillna(0.0)
                for _, row in bio_candidates.iterrows()
            },
            index=module_features.index,
        )
        frames.append(bio_frame)
        edge_frame = build_edge_features(module_features).add_prefix("edge_feature__")
        frames.append(edge_frame)
        feature_by_cohort[cohort] = pd.concat(frames, axis=1).fillna(0.0)
        for family, frame in [
            ("module", module_frame),
            ("signature", fixed_frame),
            ("bio_candidate", bio_frame),
            ("edge", edge_frame),
        ]:
            coverage_rows.append(
                {
                    "cohort": cohort,
                    "feature_family": family,
                    "n_features": int(frame.shape[1]),
                    "n_samples": int(frame.shape[0]),
                }
            )
    return feature_by_cohort, pd.DataFrame(coverage_rows)


def feature_family_columns(features: pd.DataFrame, family: str) -> list[str]:
    prefixes = {
        "module": ["module__"],
        "signature": ["signature__"],
        "compact": [
            "module__",
            "signature__EcoNiche-Opt-ModulePriorFixed",
            "signature__IFNG",
            "signature__CXCL9",
            "signature__TIG",
            "signature__APM",
            "signature__TIDE_dysfunction",
            "signature__TIDE_exclusion",
            "signature__IPRES",
        ],
        "module_signature": ["module__", "signature__"],
        "module_signature_bio": ["module__", "signature__", "bio_candidate__"],
        "module_signature_edge": ["module__", "signature__", "edge_feature__"],
        "all": ["module__", "signature__", "bio_candidate__", "edge_feature__"],
    }[family]
    if family == "compact":
        cols = [
            col
            for col in features.columns
            if col.startswith("module__")
            or col
            in {
                "signature__EcoNiche-Opt-ModulePriorFixed",
                "signature__IFNG",
                "signature__CXCL9",
                "signature__TIG",
                "signature__APM",
                "signature__TIDE_dysfunction",
                "signature__TIDE_exclusion",
                "signature__IPRES",
            }
        ]
    else:
        cols = [col for col in features.columns if any(col.startswith(prefix) for prefix in prefixes)]
    return sorted(cols)


def ensemble_specs() -> list[dict[str, object]]:
    specs: list[dict[str, object]] = []
    feature_families = ["module", "signature", "compact", "module_signature", "module_signature_bio", "module_signature_edge", "all"]
    for family in feature_families:
        for penalty in ["l1", "l2"]:
            for c_value in [0.03, 0.1, 0.3, 1.0]:
                specs.append(
                    {
                        "candidate": f"{family}_{penalty}_C{c_value:g}",
                        "feature_family": family,
                        "penalty": penalty,
                        "C": float(c_value),
                    }
                )
    return specs


def _concat_features(features_by_cohort: dict[str, pd.DataFrame], cohorts: list[str], columns: list[str]) -> pd.DataFrame:
    return pd.concat([features_by_cohort[cohort].reindex(columns=columns).fillna(0.0) for cohort in cohorts], axis=0)


def _concat_labels(y_by_cohort: dict[str, pd.Series], cohorts: list[str]) -> pd.Series:
    return pd.concat([y_by_cohort[cohort] for cohort in cohorts], axis=0).astype(int)


def _fit_model(X: pd.DataFrame, y: pd.Series, spec: dict[str, object]) -> Pipeline:
    solver = "liblinear" if spec["penalty"] == "l1" else "lbfgs"
    model = Pipeline(
        [
            ("scale", StandardScaler()),
            (
                "logistic",
                LogisticRegression(
                    C=float(spec["C"]),
                    penalty=str(spec["penalty"]),
                    solver=solver,
                    class_weight="balanced",
                    max_iter=5000,
                    random_state=20260527,
                ),
            ),
        ]
    )
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning, module="sklearn")
        warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")
        model.fit(X.astype(float), y.astype(int))
    return model


def _evaluate_spec_inner(
    spec: dict[str, object],
    train_cohorts: list[str],
    features_by_cohort: dict[str, pd.DataFrame],
    y_by_cohort: dict[str, pd.Series],
) -> dict[str, object]:
    rows: list[dict[str, float]] = []
    all_columns = feature_family_columns(next(iter(features_by_cohort.values())), str(spec["feature_family"]))
    if not all_columns:
        return {"candidate": spec["candidate"], "selection_score": -math.inf}
    for inner_holdout in train_cohorts:
        inner_train = [cohort for cohort in train_cohorts if cohort != inner_holdout]
        if not inner_train:
            continue
        y_train = _concat_labels(y_by_cohort, inner_train)
        y_test = y_by_cohort[inner_holdout].astype(int)
        if y_train.nunique() < 2 or y_test.nunique() < 2:
            continue
        X_train = _concat_features(features_by_cohort, inner_train, all_columns).reindex(y_train.index)
        X_test = features_by_cohort[inner_holdout].reindex(y_test.index).reindex(columns=all_columns).fillna(0.0)
        try:
            model = _fit_model(X_train, y_train, spec)
        except Exception:
            continue
        p_train = model.predict_proba(X_train.astype(float))[:, 1]
        p_test = model.predict_proba(X_test.astype(float))[:, 1]
        threshold = select_threshold(y_train.to_numpy(dtype=int), p_train)
        rows.append(compute_binary_metrics(y_test, p_test, threshold=threshold))
    if not rows:
        return {"candidate": spec["candidate"], "selection_score": -math.inf}
    frame = pd.DataFrame(rows)
    selection_score = float(
        frame["AUROC"].mean()
        - 0.50 * frame["AUROC"].std(ddof=0)
        + 0.10 * frame["AUPRC"].mean()
        + 0.05 * frame["balanced_accuracy"].mean()
        - 0.15 * frame["ECE"].mean()
        - 0.001 * len(all_columns)
    )
    return {
        "candidate": spec["candidate"],
        "feature_family": spec["feature_family"],
        "penalty": spec["penalty"],
        "C": spec["C"],
        "n_features": len(all_columns),
        "inner_mean_AUROC": float(frame["AUROC"].mean()),
        "inner_sd_AUROC": float(frame["AUROC"].std(ddof=0)),
        "inner_mean_AUPRC": float(frame["AUPRC"].mean()),
        "inner_mean_balanced_accuracy": float(frame["balanced_accuracy"].mean()),
        "inner_mean_ECE": float(frame["ECE"].mean()),
        "selection_score": selection_score,
    }


def select_spec(
    specs: list[dict[str, object]],
    train_cohorts: list[str],
    features_by_cohort: dict[str, pd.DataFrame],
    y_by_cohort: dict[str, pd.Series],
) -> tuple[dict[str, object], pd.DataFrame]:
    rows = [_evaluate_spec_inner(spec, train_cohorts, features_by_cohort, y_by_cohort) for spec in specs]
    selection = pd.DataFrame(rows)
    finite = selection[np.isfinite(selection["selection_score"].astype(float))]
    if finite.empty:
        chosen = next(spec for spec in specs if spec["candidate"] == "compact_l2_C0.1")
    else:
        chosen_name = str(finite.sort_values("selection_score", ascending=False).iloc[0]["candidate"])
        chosen = next(spec for spec in specs if spec["candidate"] == chosen_name)
    selection["selected_candidate"] = chosen["candidate"]
    return chosen, selection


def fit_predict_holdout(
    spec: dict[str, object],
    train_cohorts: list[str],
    holdout: str,
    features_by_cohort: dict[str, pd.DataFrame],
    y_by_cohort: dict[str, pd.Series],
) -> tuple[pd.Series, pd.Series, float, list[str]]:
    columns = feature_family_columns(next(iter(features_by_cohort.values())), str(spec["feature_family"]))
    y_train = _concat_labels(y_by_cohort, train_cohorts)
    y_test = y_by_cohort[holdout].astype(int)
    X_train = _concat_features(features_by_cohort, train_cohorts, columns).reindex(y_train.index)
    X_test = features_by_cohort[holdout].reindex(y_test.index).reindex(columns=columns).fillna(0.0)
    model = _fit_model(X_train, y_train, spec)
    p_train = model.predict_proba(X_train.astype(float))[:, 1]
    p_test = pd.Series(model.predict_proba(X_test.astype(float))[:, 1], index=y_test.index)
    threshold = select_threshold(y_train.to_numpy(dtype=int), p_train)
    return y_test, p_test, float(threshold), columns


def evaluate_lodo(
    X_by_cohort: dict[str, pd.DataFrame],
    metadata_by_cohort: dict[str, pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    active = _active_real_cohorts(X_by_cohort, metadata_by_cohort)
    endpoint_data = prepare_endpoint_data(X_by_cohort, metadata_by_cohort, active, PRIMARY_ENDPOINT)
    module_features, module_coverage = build_module_features_by_cohort(endpoint_data.X_by_cohort)
    features, feature_coverage = build_feature_matrix_by_cohort(endpoint_data.X_by_cohort, module_features)
    strata = default_strata(endpoint_data.X_by_cohort)
    specs = ensemble_specs()
    metric_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    selection_rows: list[dict[str, object]] = []
    for stratum in PRIMARY_STRATA:
        spec_frame = strata.get(stratum, {})
        train_pool = [cohort for cohort in spec_frame.get("train_pool", []) if cohort in endpoint_data.y_response_by_cohort]
        holdouts = [cohort for cohort in spec_frame.get("holdouts", []) if cohort in endpoint_data.y_response_by_cohort]
        for holdout in holdouts:
            train_cohorts = [cohort for cohort in train_pool if cohort != holdout and cohort in endpoint_data.y_response_by_cohort]
            if len(train_cohorts) < 2:
                continue
            y_train = _concat_labels(endpoint_data.y_response_by_cohort, train_cohorts)
            y_test = endpoint_data.y_response_by_cohort[holdout].astype(int)
            if y_train.nunique() < 2 or y_test.nunique() < 2:
                continue
            chosen, selection = select_spec(specs, train_cohorts, features, endpoint_data.y_response_by_cohort)
            selection.insert(0, "endpoint", PRIMARY_ENDPOINT)
            selection.insert(1, "stratum", stratum)
            selection.insert(2, "holdout", holdout)
            selection.insert(3, "train_cohorts", ",".join(train_cohorts))
            selection_rows.extend(selection.to_dict("records"))
            y_eval, prob, threshold, columns = fit_predict_holdout(
                chosen, train_cohorts, holdout, features, endpoint_data.y_response_by_cohort
            )
            metrics = compute_binary_metrics(y_eval, prob, threshold=threshold)
            metric_rows.append(
                {
                    **metrics,
                    "endpoint": PRIMARY_ENDPOINT,
                    "stratum": stratum,
                    "cohort": holdout,
                    "model_name": ENSEMBLE_MODEL,
                    "selected_candidate": chosen["candidate"],
                    "feature_family": chosen["feature_family"],
                    "penalty": chosen["penalty"],
                    "C": chosen["C"],
                    "n_features": len(columns),
                    "n_samples": int(len(y_eval)),
                    "n_responders": int(y_eval.sum()),
                    "n_nonresponders": int((y_eval == 0).sum()),
                    "threshold": threshold,
                    "train_cohorts": ",".join(train_cohorts),
                    "selection_boundary": "training_only_nested_lodo",
                }
            )
            for sample_id in y_eval.index:
                prediction_rows.append(
                    {
                        "endpoint": PRIMARY_ENDPOINT,
                        "stratum": stratum,
                        "cohort": holdout,
                        "sample_id": sample_id,
                        "model_name": ENSEMBLE_MODEL,
                        "true_response_label": int(y_eval.loc[sample_id]),
                        "response_probability": float(prob.loc[sample_id]),
                        "threshold": threshold,
                        "pred_response_label": int(prob.loc[sample_id] >= threshold),
                        "selected_candidate": chosen["candidate"],
                        "selection_boundary": "training_only_nested_lodo",
                    }
                )
    coverage = pd.concat([module_coverage.assign(feature_family="module_gene_coverage"), feature_coverage], ignore_index=True, sort=False)
    return pd.DataFrame(metric_rows), pd.DataFrame(prediction_rows), pd.DataFrame(selection_rows), coverage


def evaluate_strict_external(
    X_by_cohort: dict[str, pd.DataFrame],
    metadata_by_cohort: dict[str, pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    endpoint_data = prepare_endpoint_data(X_by_cohort, metadata_by_cohort, [*DISCOVERY_COHORTS, *STRICT_EXTERNAL_COHORTS], STRICT_ENDPOINT)
    module_features, _ = build_module_features_by_cohort(endpoint_data.X_by_cohort)
    features, _ = build_feature_matrix_by_cohort(endpoint_data.X_by_cohort, module_features)
    specs = ensemble_specs()
    train_cohorts = [cohort for cohort in DISCOVERY_COHORTS if cohort in endpoint_data.y_response_by_cohort]
    chosen, selection = select_spec(specs, train_cohorts, features, endpoint_data.y_response_by_cohort)
    selection.insert(0, "endpoint", STRICT_ENDPOINT)
    selection.insert(1, "stratum", "strict_melanoma_pd1_like_external")
    selection.insert(2, "holdout", "+".join(STRICT_EXTERNAL_COHORTS))
    selection.insert(3, "train_cohorts", ",".join(train_cohorts))
    metric_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    for holdout in STRICT_EXTERNAL_COHORTS:
        if holdout not in endpoint_data.y_response_by_cohort:
            continue
        y_test = endpoint_data.y_response_by_cohort[holdout].astype(int)
        if y_test.nunique() < 2:
            continue
        y_eval, prob, threshold, columns = fit_predict_holdout(
            chosen, train_cohorts, holdout, features, endpoint_data.y_response_by_cohort
        )
        metrics = compute_binary_metrics(y_eval, prob, threshold=threshold)
        metric_rows.append(
            {
                **metrics,
                "endpoint": STRICT_ENDPOINT,
                "stratum": "strict_melanoma_pd1_like_external",
                "cohort": holdout,
                "model_name": ENSEMBLE_MODEL,
                "selected_candidate": chosen["candidate"],
                "feature_family": chosen["feature_family"],
                "penalty": chosen["penalty"],
                "C": chosen["C"],
                "n_features": len(columns),
                "n_samples": int(len(y_eval)),
                "n_responders": int(y_eval.sum()),
                "n_nonresponders": int((y_eval == 0).sum()),
                "threshold": threshold,
                "train_cohorts": ",".join(train_cohorts),
                "selection_boundary": "discovery_only_nested_lodo_no_external_selection",
            }
        )
        for sample_id in y_eval.index:
            prediction_rows.append(
                {
                    "endpoint": STRICT_ENDPOINT,
                    "stratum": "strict_melanoma_pd1_like_external",
                    "cohort": holdout,
                    "sample_id": sample_id,
                    "model_name": ENSEMBLE_MODEL,
                    "true_response_label": int(y_eval.loc[sample_id]),
                    "response_probability": float(prob.loc[sample_id]),
                    "threshold": threshold,
                    "pred_response_label": int(prob.loc[sample_id] >= threshold),
                    "selected_candidate": chosen["candidate"],
                    "selection_boundary": "discovery_only_nested_lodo_no_external_selection",
                }
            )
    return pd.DataFrame(metric_rows), pd.DataFrame(prediction_rows), selection


def summarize_predictions(predictions: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if predictions.empty:
        return pd.DataFrame()
    for (endpoint, stratum), frame in predictions.groupby(["endpoint", "stratum"]):
        y = frame["true_response_label"].astype(int)
        prob = frame["response_probability"].astype(float)
        threshold = frame["threshold"].astype(float)
        metrics = compute_binary_metrics(y, prob, threshold=float(threshold.median()))
        fold_metrics = []
        for _, fold in frame.groupby("cohort"):
            if fold["true_response_label"].nunique() == 2:
                fold_metrics.append(compute_binary_metrics(fold["true_response_label"].astype(int), fold["response_probability"].astype(float)))
        fold_frame = pd.DataFrame(fold_metrics)
        rows.append(
            {
                **metrics,
                "endpoint": endpoint,
                "stratum": stratum,
                "model_name": ENSEMBLE_MODEL,
                "n_samples": int(len(frame)),
                "n_cohorts": int(frame["cohort"].nunique()),
                "n_responders": int(y.sum()),
                "n_nonresponders": int((y == 0).sum()),
                "response_prevalence": float(y.mean()),
                "AUPRC_minus_prevalence": float(metrics["AUPRC"] - y.mean()),
                "fold_AUROC_sd": float(fold_frame["AUROC"].std(ddof=0)) if not fold_frame.empty else math.nan,
                "selected_candidates": ",".join(sorted(frame["selected_candidate"].astype(str).unique())),
                "selection_boundary": ",".join(sorted(frame["selection_boundary"].astype(str).unique())),
            }
        )
    return pd.DataFrame(rows).sort_values(["endpoint", "stratum"]).reset_index(drop=True)


def compare_to_baseline(predictions: pd.DataFrame, baseline_predictions: pd.DataFrame) -> pd.DataFrame:
    if predictions.empty or baseline_predictions.empty:
        return pd.DataFrame()
    rows: list[dict[str, object]] = []
    for (endpoint, stratum), frame in predictions.groupby(["endpoint", "stratum"]):
        baseline = baseline_predictions[
            baseline_predictions["endpoint"].astype(str).eq(endpoint)
            & baseline_predictions["stratum"].astype(str).eq(stratum)
            & baseline_predictions["model_name"].astype(str).eq(BASELINE_MODEL)
        ]
        if baseline.empty:
            continue
        merged = frame.merge(
            baseline[["cohort", "sample_id", "true_response_label", "response_probability"]],
            on=["cohort", "sample_id", "true_response_label"],
            suffixes=("_ensemble", "_baseline"),
        )
        if len(merged) < 8 or merged["true_response_label"].nunique() < 2:
            continue
        y = merged["true_response_label"].astype(int)
        ensemble_prob = merged["response_probability_ensemble"].astype(float)
        baseline_prob = merged["response_probability_baseline"].astype(float)
        ensemble_metrics = compute_binary_metrics(y, ensemble_prob)
        baseline_metrics = compute_binary_metrics(y, baseline_prob)
        stats = paired_bootstrap_delta(y, ensemble_prob, baseline_prob, n_bootstrap=1000, random_state=20260527)
        rows.append(
            {
                "endpoint": endpoint,
                "stratum": stratum,
                "target_model": ENSEMBLE_MODEL,
                "baseline_model": BASELINE_MODEL,
                "n_samples": int(len(merged)),
                "target_AUROC": ensemble_metrics["AUROC"],
                "baseline_AUROC": baseline_metrics["AUROC"],
                "delta_AUROC": float(ensemble_metrics["AUROC"] - baseline_metrics["AUROC"]),
                "target_AUPRC": ensemble_metrics["AUPRC"],
                "baseline_AUPRC": baseline_metrics["AUPRC"],
                "delta_AUPRC": float(ensemble_metrics["AUPRC"] - baseline_metrics["AUPRC"]),
                "target_ECE": ensemble_metrics["ECE"],
                "baseline_ECE": baseline_metrics["ECE"],
                "delta_ECE": float(ensemble_metrics["ECE"] - baseline_metrics["ECE"]),
                **stats,
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["fdr_q"] = benjamini_hochberg(result["p_value"].fillna(1.0))
        result["claim_level"] = np.where(
            (result["delta_AUROC"] > 0) & (result["fdr_q"] <= 0.05),
            "FDR_supported_ensemble_gain",
            np.where(result["delta_AUROC"] > 0, "point_estimate_ensemble_gain", "ensemble_not_supported"),
        )
    return result


def write_audit(out_dir: Path, summary: pd.DataFrame, comparison: pd.DataFrame) -> None:
    lines = [
        "# No-leakage Ensemble Search Audit",
        "",
        "This registered search tests whether a training-only stacked ensemble can improve primary melanoma and strict external evidence. Feature family, regularization, threshold, and calibration are selected only within training folds or discovery cohorts. Strict external labels are used only for final scoring.",
        "",
        "## Summary",
        "",
    ]
    if summary.empty:
        lines.append("No ensemble summary rows were produced.")
    else:
        for _, row in summary.iterrows():
            lines.append(
                f"- {row['endpoint']} / {row['stratum']}: n={int(row['n_samples'])}, "
                f"AUROC={float(row['AUROC']):.3f}, AUPRC={float(row['AUPRC']):.3f}, "
                f"BA={float(row['balanced_accuracy']):.3f}, ECE={float(row['ECE']):.3f}; "
                f"selected={row['selected_candidates']}."
            )
    lines.extend(["", "## Comparison To ModulePriorFixed", ""])
    if comparison.empty:
        lines.append("No paired comparison rows were produced.")
    else:
        for _, row in comparison.iterrows():
            lines.append(
                f"- {row['endpoint']} / {row['stratum']}: delta AUROC={float(row['delta_AUROC']):.3f}, "
                f"delta AUPRC={float(row['delta_AUPRC']):.3f}, delta ECE={float(row['delta_ECE']):.3f}, "
                f"95% CI [{float(row['ci_low']):.3f}, {float(row['ci_high']):.3f}], "
                f"FDR q={float(row['fdr_q']):.3f}; {row['claim_level']}."
            )
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "Only promote the ensemble as the primary model if it improves primary LODO and does not degrade strict external validation. Otherwise, retain it as a negative optimization audit.",
        ]
    )
    (out_dir / "NO_LEAKAGE_ENSEMBLE_SEARCH_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")


def run(processed_dir: Path, out_dir: Path) -> None:
    X_by_cohort, _, metadata_by_cohort = load_processed_bulk(processed_dir)
    lodo_metrics, lodo_predictions, lodo_selection, coverage = evaluate_lodo(X_by_cohort, metadata_by_cohort)
    external_metrics, external_predictions, external_selection = evaluate_strict_external(X_by_cohort, metadata_by_cohort)
    predictions = pd.concat([lodo_predictions, external_predictions], ignore_index=True, sort=False)
    metrics = pd.concat([lodo_metrics, external_metrics], ignore_index=True, sort=False)
    selection = pd.concat([lodo_selection, external_selection], ignore_index=True, sort=False)
    summary = summarize_predictions(predictions)

    baseline_predictions_path = ROOT / "results" / "endpoint_modules_heuristic_deep_primary_20260519" / "endpoint_module_predictions.tsv"
    baseline_predictions = pd.read_csv(baseline_predictions_path, sep="\t") if baseline_predictions_path.exists() else pd.DataFrame()
    comparison = compare_to_baseline(predictions, baseline_predictions)

    out_dir.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(out_dir / "no_leakage_ensemble_metrics.tsv", sep="\t", index=False)
    predictions.to_csv(out_dir / "no_leakage_ensemble_predictions.tsv", sep="\t", index=False)
    selection.to_csv(out_dir / "no_leakage_ensemble_inner_selection.tsv", sep="\t", index=False)
    summary.to_csv(out_dir / "no_leakage_ensemble_summary.tsv", sep="\t", index=False)
    comparison.to_csv(out_dir / "no_leakage_ensemble_vs_module_prior.tsv", sep="\t", index=False)
    coverage.to_csv(out_dir / "no_leakage_ensemble_feature_coverage.tsv", sep="\t", index=False)
    write_audit(out_dir, summary, comparison)
    print(f"Wrote no-leakage ensemble search to {out_dir}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--out", default="results/no_leakage_ensemble_search_20260527")
    args = parser.parse_args()
    run(ROOT / args.processed_dir, ROOT / args.out)


if __name__ == "__main__":
    main()
