from __future__ import annotations

import argparse
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
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.io import load_processed_bulk
from econiche.metrics import compute_binary_metrics
from econiche_opt.model.endpoint_modules import (
    build_module_features_by_cohort,
    default_strata,
    module_prior_score,
    prepare_endpoint_data,
    select_threshold,
    sigmoid,
)


PRIMARY_ENDPOINT = "primary_recist"
STRICT_ENDPOINT = "strict_recist"
PRIMARY_STRATA = ["melanoma_core_high_evidence", "melanoma_recist_supported_primary"]
DISCOVERY_COHORTS = ["GSE91061", "GSE78220", "PRJEB23709_PD1_PRE"]
STRICT_EXTERNAL_COHORTS = ["GSE145996", "PHS000452_LIU_LIKE_PRE"]
RESPONSE_MODULES = ["ifn_t_cell_inflamed", "cytotoxic_cd8", "exhaustion_checkpoint", "antigen_presentation", "trm_tls"]
RESISTANCE_MODULES = ["myeloid_suppression", "stromal_exclusion"]
EDGE_C_VALUES = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0]


def build_edge_features(features: pd.DataFrame) -> pd.DataFrame:
    columns = [column for column in [*RESPONSE_MODULES, *RESISTANCE_MODULES] if column in features.columns]
    out = features[columns].copy()
    for response in RESPONSE_MODULES:
        for resistance in RESISTANCE_MODULES:
            if response in features.columns and resistance in features.columns:
                out[f"edge__{response}__x__{resistance}"] = features[response] * features[resistance]
                out[f"contrast__{response}__minus__{resistance}"] = features[response] - features[resistance]
    synergy_pairs = [
        ("ifn_t_cell_inflamed", "antigen_presentation"),
        ("ifn_t_cell_inflamed", "cytotoxic_cd8"),
        ("cytotoxic_cd8", "exhaustion_checkpoint"),
    ]
    for left, right in synergy_pairs:
        if left in features.columns and right in features.columns:
            out[f"synergy__{left}__x__{right}"] = features[left] * features[right]
    return out.fillna(0.0)


def _fit_edge_model(X_train: pd.DataFrame, y_train: pd.Series, X_test: pd.DataFrame, c_value: float) -> tuple[np.ndarray, np.ndarray]:
    model = Pipeline(
        [
            ("scale", StandardScaler()),
            (
                "logistic",
                LogisticRegression(
                    C=float(c_value),
                    class_weight="balanced",
                    solver="lbfgs",
                    max_iter=10000,
                    random_state=20260527,
                ),
            ),
        ]
    )
    model.fit(X_train.astype(float), y_train.astype(int))
    return model.predict_proba(X_train.astype(float))[:, 1], model.predict_proba(X_test.astype(float))[:, 1]


def _concat_series(y_by_cohort: dict[str, pd.Series], cohorts: list[str]) -> pd.Series:
    return pd.concat([y_by_cohort[cohort] for cohort in cohorts]).astype(int)


def _inner_select_c(module_features: dict[str, pd.DataFrame], y_by_cohort: dict[str, pd.Series], train_cohorts: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for c_value in EDGE_C_VALUES:
        metrics_rows: list[dict[str, float]] = []
        for inner_holdout in train_cohorts:
            inner_train = [cohort for cohort in train_cohorts if cohort != inner_holdout]
            if not inner_train:
                continue
            y_train = _concat_series(y_by_cohort, inner_train)
            y_test = y_by_cohort[inner_holdout].astype(int)
            if y_train.nunique() < 2 or y_test.nunique() < 2:
                continue
            X_train = pd.concat([build_edge_features(module_features[cohort]) for cohort in inner_train])
            X_test = build_edge_features(module_features[inner_holdout])
            train_prob, test_prob = _fit_edge_model(X_train, y_train, X_test, c_value)
            threshold = select_threshold(y_train.to_numpy(dtype=int), train_prob)
            metrics_rows.append(compute_binary_metrics(y_test, test_prob, threshold=threshold))
        if not metrics_rows:
            continue
        frame = pd.DataFrame(metrics_rows)
        score = float(
            frame["AUROC"].mean()
            - 0.50 * frame["AUROC"].std(ddof=0)
            + 0.05 * frame["AUPRC"].mean()
            + 0.05 * frame["balanced_accuracy"].mean()
            - 0.15 * frame["ECE"].mean()
        )
        rows.append(
            {
                "C": float(c_value),
                "inner_mean_AUROC": float(frame["AUROC"].mean()),
                "inner_sd_AUROC": float(frame["AUROC"].std(ddof=0)),
                "inner_mean_AUPRC": float(frame["AUPRC"].mean()),
                "inner_mean_balanced_accuracy": float(frame["balanced_accuracy"].mean()),
                "inner_mean_ECE": float(frame["ECE"].mean()),
                "selection_score": score,
            }
        )
    return pd.DataFrame(rows).sort_values("selection_score", ascending=False).reset_index(drop=True)


def _score_module_prior(
    module_features: dict[str, pd.DataFrame],
    y_by_cohort: dict[str, pd.Series],
    train_cohorts: list[str],
    holdout: str,
) -> tuple[np.ndarray, np.ndarray, float]:
    train_score = pd.concat([module_prior_score(module_features[cohort]) for cohort in train_cohorts])
    test_score = module_prior_score(module_features[holdout])
    train_prob = sigmoid(train_score)
    test_prob = sigmoid(test_score)
    y_train = _concat_series(y_by_cohort, train_cohorts)
    threshold = select_threshold(y_train.to_numpy(dtype=int), train_prob)
    return train_prob, test_prob, threshold


def _evaluate_primary_lodo(
    X_by_cohort: dict[str, pd.DataFrame],
    metadata_by_cohort: dict[str, pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    active = [cohort for cohort in X_by_cohort if not cohort.startswith("demo_cohort_")]
    strata = default_strata(active)
    metric_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    inner_rows: list[dict[str, object]] = []
    for stratum in PRIMARY_STRATA:
        endpoint_data = prepare_endpoint_data(X_by_cohort, metadata_by_cohort, strata[stratum]["cohorts"], PRIMARY_ENDPOINT)
        module_features, _ = build_module_features_by_cohort(endpoint_data.X_by_cohort)
        for holdout in strata[stratum]["holdouts"]:
            if holdout not in module_features:
                continue
            train_cohorts = [cohort for cohort in strata[stratum]["train_pool"] if cohort != holdout and cohort in module_features]
            if len(train_cohorts) < 2:
                continue
            y_train = _concat_series(endpoint_data.y_response_by_cohort, train_cohorts)
            y_test = endpoint_data.y_response_by_cohort[holdout].astype(int)
            if y_train.nunique() < 2 or y_test.nunique() < 2:
                continue
            selection = _inner_select_c(module_features, endpoint_data.y_response_by_cohort, train_cohorts)
            if selection.empty:
                continue
            selection.insert(0, "endpoint", PRIMARY_ENDPOINT)
            selection.insert(1, "stratum", stratum)
            selection.insert(2, "holdout", holdout)
            inner_rows.extend(selection.to_dict("records"))
            c_value = float(selection.iloc[0]["C"])

            X_train = pd.concat([build_edge_features(module_features[cohort]) for cohort in train_cohorts])
            X_test = build_edge_features(module_features[holdout])
            train_prob, test_prob = _fit_edge_model(X_train, y_train, X_test, c_value)
            threshold = select_threshold(y_train.to_numpy(dtype=int), train_prob)
            variants = {
                "EcoNiche-Opt-InteractionEdgeLogistic": (test_prob, threshold, c_value),
                "EcoNiche-Opt-ModulePriorComposite": _score_module_prior(
                    module_features, endpoint_data.y_response_by_cohort, train_cohorts, holdout
                )[1:],
            }
            for model_name, values in variants.items():
                if model_name == "EcoNiche-Opt-InteractionEdgeLogistic":
                    prob, model_threshold, model_c = values
                else:
                    prob, model_threshold = values
                    model_c = np.nan
                metrics = compute_binary_metrics(y_test, prob, threshold=float(model_threshold))
                metric_rows.append(
                    {
                        **metrics,
                        "endpoint": PRIMARY_ENDPOINT,
                        "stratum": stratum,
                        "cohort": holdout,
                        "model_name": model_name,
                        "C": model_c,
                        "threshold": float(model_threshold),
                        "n_samples": int(len(y_test)),
                        "n_responders": int(y_test.sum()),
                        "n_nonresponders": int((y_test == 0).sum()),
                        "train_cohorts": ",".join(train_cohorts),
                        "selection_boundary": "training_only_inner_lodo",
                    }
                )
                for sample_id, sample_prob in zip(y_test.index, prob):
                    prediction_rows.append(
                        {
                            "endpoint": PRIMARY_ENDPOINT,
                            "stratum": stratum,
                            "cohort": holdout,
                            "model_name": model_name,
                            "sample_id": sample_id,
                            "true_response_label": int(y_test.loc[sample_id]),
                            "response_probability": float(sample_prob),
                            "threshold": float(model_threshold),
                        }
                    )
    return pd.DataFrame(metric_rows), pd.DataFrame(prediction_rows), pd.DataFrame(inner_rows)


def _summarize_predictions(predictions: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if predictions.empty:
        return pd.DataFrame()
    for keys, frame in predictions.groupby(["endpoint", "stratum", "model_name"], dropna=False):
        endpoint, stratum, model_name = keys
        y = frame["true_response_label"].astype(int)
        prob = frame["response_probability"].astype(float)
        metrics = compute_binary_metrics(y, prob)
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
            }
        )
    return pd.DataFrame(rows).sort_values(["endpoint", "stratum", "model_name"]).reset_index(drop=True)


def _compare_models(summary: pd.DataFrame, metrics: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if summary.empty:
        return pd.DataFrame()
    for (endpoint, stratum), frame in summary.groupby(["endpoint", "stratum"], dropna=False):
        edge = frame[frame["model_name"].eq("EcoNiche-Opt-InteractionEdgeLogistic")]
        base = frame[frame["model_name"].eq("EcoNiche-Opt-ModulePriorComposite")]
        if edge.empty or base.empty:
            continue
        edge_row = edge.iloc[0]
        base_row = base.iloc[0]
        fold = metrics[metrics["endpoint"].eq(endpoint) & metrics["stratum"].eq(stratum)]
        fold_sd = fold.groupby("model_name")["AUROC"].std(ddof=0).to_dict()
        rows.append(
            {
                "endpoint": endpoint,
                "stratum": stratum,
                "edge_model": edge_row["model_name"],
                "baseline_model": base_row["model_name"],
                "n_samples": int(edge_row["n_samples"]),
                "edge_AUROC": float(edge_row["AUROC"]),
                "baseline_AUROC": float(base_row["AUROC"]),
                "delta_AUROC": float(edge_row["AUROC"] - base_row["AUROC"]),
                "edge_AUPRC": float(edge_row["AUPRC"]),
                "baseline_AUPRC": float(base_row["AUPRC"]),
                "delta_AUPRC": float(edge_row["AUPRC"] - base_row["AUPRC"]),
                "edge_ECE": float(edge_row["ECE"]),
                "baseline_ECE": float(base_row["ECE"]),
                "delta_ECE": float(edge_row["ECE"] - base_row["ECE"]),
                "edge_fold_AUROC_sd": float(fold_sd.get("EcoNiche-Opt-InteractionEdgeLogistic", np.nan)),
                "baseline_fold_AUROC_sd": float(fold_sd.get("EcoNiche-Opt-ModulePriorComposite", np.nan)),
                "delta_fold_AUROC_sd": float(
                    fold_sd.get("EcoNiche-Opt-InteractionEdgeLogistic", np.nan)
                    - fold_sd.get("EcoNiche-Opt-ModulePriorComposite", np.nan)
                ),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["claim_level"] = np.where(
        out["delta_AUROC"] > 0,
        "interaction_edge_discrimination_gain",
        np.where(
            (out["delta_ECE"] < 0) & (out["delta_fold_AUROC_sd"] <= 0),
            "interaction_edge_calibration_and_stability_tradeoff",
            np.where(out["delta_ECE"] < 0, "interaction_edge_calibration_tradeoff", "interaction_edge_not_supported"),
        ),
    )
    return out


def _evaluate_strict_external(
    X_by_cohort: dict[str, pd.DataFrame],
    metadata_by_cohort: dict[str, pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    endpoint_data = prepare_endpoint_data(X_by_cohort, metadata_by_cohort, [*DISCOVERY_COHORTS, *STRICT_EXTERNAL_COHORTS], STRICT_ENDPOINT)
    module_features, _ = build_module_features_by_cohort(endpoint_data.X_by_cohort)
    train_cohorts = [cohort for cohort in DISCOVERY_COHORTS if cohort in module_features]
    selection = _inner_select_c(module_features, endpoint_data.y_response_by_cohort, train_cohorts)
    if selection.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    c_value = float(selection.iloc[0]["C"])
    y_train = _concat_series(endpoint_data.y_response_by_cohort, train_cohorts)
    X_train = pd.concat([build_edge_features(module_features[cohort]) for cohort in train_cohorts])
    metric_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    for holdout in STRICT_EXTERNAL_COHORTS:
        if holdout not in module_features:
            continue
        y_test = endpoint_data.y_response_by_cohort[holdout].astype(int)
        if y_test.nunique() < 2:
            continue
        X_test = build_edge_features(module_features[holdout])
        train_prob, edge_prob = _fit_edge_model(X_train, y_train, X_test, c_value)
        edge_threshold = select_threshold(y_train.to_numpy(dtype=int), train_prob)
        _, base_prob, base_threshold = _score_module_prior(module_features, endpoint_data.y_response_by_cohort, train_cohorts, holdout)
        variants = {
            "EcoNiche-Opt-InteractionEdgeLogistic": (edge_prob, edge_threshold, c_value),
            "EcoNiche-Opt-ModulePriorComposite": (base_prob, base_threshold, np.nan),
        }
        for model_name, (prob, threshold, model_c) in variants.items():
            metrics = compute_binary_metrics(y_test, prob, threshold=float(threshold))
            metric_rows.append(
                {
                    **metrics,
                    "endpoint": STRICT_ENDPOINT,
                    "stratum": "strict_pd1_like_external",
                    "cohort": holdout,
                    "model_name": model_name,
                    "C": model_c,
                    "threshold": float(threshold),
                    "n_samples": int(len(y_test)),
                    "n_responders": int(y_test.sum()),
                    "n_nonresponders": int((y_test == 0).sum()),
                    "train_cohorts": ",".join(train_cohorts),
                    "selection_boundary": "discovery_only_inner_lodo_no_external_selection",
                }
            )
            for sample_id, sample_prob in zip(y_test.index, prob):
                prediction_rows.append(
                    {
                        "endpoint": STRICT_ENDPOINT,
                        "stratum": "strict_pd1_like_external",
                        "cohort": holdout,
                        "model_name": model_name,
                        "sample_id": sample_id,
                        "true_response_label": int(y_test.loc[sample_id]),
                        "response_probability": float(sample_prob),
                        "threshold": float(threshold),
                    }
                )
    selection.insert(0, "endpoint", STRICT_ENDPOINT)
    selection.insert(1, "selection_set", ",".join(train_cohorts))
    return pd.DataFrame(metric_rows), pd.DataFrame(prediction_rows), selection


def _write_audit(out_dir: Path, comparisons: pd.DataFrame, strict_comparisons: pd.DataFrame) -> None:
    lines = [
        "# Aligned Interaction Edge Audit",
        "",
        "This registered analysis tests whether response-resistance module interaction edges improve the current module-prior model under training-only selection. It does not use strict external labels for feature selection, thresholding, calibration or regularization selection.",
        "",
        "## Primary LODO",
        "",
    ]
    if comparisons.empty:
        lines.append("No primary comparison rows were produced.")
    else:
        for _, row in comparisons.iterrows():
            lines.append(
                f"- {row['stratum']}: edge AUROC={row['edge_AUROC']:.3f} versus baseline AUROC={row['baseline_AUROC']:.3f}; "
                f"delta AUROC={row['delta_AUROC']:.3f}; delta ECE={row['delta_ECE']:.3f}; "
                f"delta fold-AUROC SD={row['delta_fold_AUROC_sd']:.3f}; claim={row['claim_level']}."
            )
    lines.extend(["", "## Strict external", ""])
    if strict_comparisons.empty:
        lines.append("No strict external comparison rows were produced.")
    else:
        row = strict_comparisons.iloc[0]
        lines.append(
            f"- {row['stratum']}: edge AUROC={row['edge_AUROC']:.3f} versus baseline AUROC={row['baseline_AUROC']:.3f}; "
            f"delta AUROC={row['delta_AUROC']:.3f}; delta ECE={row['delta_ECE']:.3f}; claim={row['claim_level']}."
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "Interaction edges are not currently supported as a discrimination-improving component of the locked melanoma predictor. They may be retained as an ecological interpretation and calibration/stability diagnostic only where the audit shows lower ECE or lower fold variability.",
        ]
    )
    (out_dir / "ALIGNED_INTERACTION_EDGE_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--out", default="results/aligned_interaction_edge_audit_20260527")
    args = parser.parse_args()

    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=UserWarning)
    out_dir = ROOT / args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    X_by_cohort, _, metadata_by_cohort = load_processed_bulk(ROOT / args.processed_dir)
    primary_metrics, primary_predictions, primary_inner = _evaluate_primary_lodo(X_by_cohort, metadata_by_cohort)
    primary_summary = _summarize_predictions(primary_predictions)
    primary_comparisons = _compare_models(primary_summary, primary_metrics)
    strict_metrics, strict_predictions, strict_inner = _evaluate_strict_external(X_by_cohort, metadata_by_cohort)
    strict_summary = _summarize_predictions(strict_predictions)
    strict_comparisons = _compare_models(strict_summary, strict_metrics)

    primary_metrics.to_csv(out_dir / "aligned_interaction_edge_lodo_metrics.tsv", sep="\t", index=False)
    primary_predictions.to_csv(out_dir / "aligned_interaction_edge_lodo_predictions.tsv", sep="\t", index=False)
    primary_inner.to_csv(out_dir / "aligned_interaction_edge_lodo_inner_selection.tsv", sep="\t", index=False)
    primary_summary.to_csv(out_dir / "aligned_interaction_edge_lodo_summary.tsv", sep="\t", index=False)
    primary_comparisons.to_csv(out_dir / "aligned_interaction_edge_lodo_comparison.tsv", sep="\t", index=False)
    strict_metrics.to_csv(out_dir / "aligned_interaction_edge_strict_external_metrics.tsv", sep="\t", index=False)
    strict_predictions.to_csv(out_dir / "aligned_interaction_edge_strict_external_predictions.tsv", sep="\t", index=False)
    strict_inner.to_csv(out_dir / "aligned_interaction_edge_strict_external_inner_selection.tsv", sep="\t", index=False)
    strict_summary.to_csv(out_dir / "aligned_interaction_edge_strict_external_summary.tsv", sep="\t", index=False)
    strict_comparisons.to_csv(out_dir / "aligned_interaction_edge_strict_external_comparison.tsv", sep="\t", index=False)
    _write_audit(out_dir, primary_comparisons, strict_comparisons)

    if primary_comparisons.empty:
        print("Aligned interaction edge audit: no primary rows")
    else:
        core = primary_comparisons[primary_comparisons["stratum"].eq("melanoma_core_high_evidence")].iloc[0]
        print(
            "Aligned interaction edge core: "
            f"edge_AUROC={core['edge_AUROC']:.3f}; baseline_AUROC={core['baseline_AUROC']:.3f}; "
            f"delta_ECE={core['delta_ECE']:.3f}; claim={core['claim_level']}"
        )
    if not strict_comparisons.empty:
        row = strict_comparisons.iloc[0]
        print(
            "Aligned interaction edge strict external: "
            f"edge_AUROC={row['edge_AUROC']:.3f}; baseline_AUROC={row['baseline_AUROC']:.3f}; "
            f"delta_ECE={row['delta_ECE']:.3f}; claim={row['claim_level']}"
        )
    print(f"Wrote {out_dir}")


if __name__ == "__main__":
    main()
