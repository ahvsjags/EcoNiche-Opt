from __future__ import annotations

import argparse
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
from econiche.normalize import rank_gaussian_normalize
from econiche.statistics import benjamini_hochberg, paired_bootstrap_delta
from econiche_opt.model.endpoint_modules import (
    MODULE_GENE_SETS,
    MODULE_PRIOR_WEIGHTS,
    build_module_features_by_cohort,
    default_strata,
    endpoint_label_series,
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
CBIO_EXTERNAL_COHORTS = ["GSE145996", "CBIO_LIU_DFCI_2019_PRE"]
BASELINE_MODEL = "EcoNiche-Opt-RawModulePrior"


VARIANT_SPECS = [
    {
        "model_name": "EcoNiche-Opt-RankModulePrior",
        "rank": True,
        "gene_signed": False,
        "weight_mode": "prior",
        "calibrate": False,
    },
    {
        "model_name": "EcoNiche-Opt-RankModulePriorCalibrated",
        "rank": True,
        "gene_signed": False,
        "weight_mode": "prior",
        "calibrate": True,
    },
    {
        "model_name": "EcoNiche-Opt-SignedRankResponseModules",
        "rank": True,
        "gene_signed": True,
        "weight_mode": "absolute_prior",
        "calibrate": False,
    },
    {
        "model_name": "EcoNiche-Opt-SignedRankResponseModulesCalibrated",
        "rank": True,
        "gene_signed": True,
        "weight_mode": "absolute_prior",
        "calibrate": True,
    },
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


def module_genes() -> list[str]:
    genes: set[str] = set()
    for values in MODULE_GENE_SETS.values():
        genes.update(values)
    return sorted(genes)


def rank_by_cohort(X_by_cohort: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    genes = module_genes()
    ranked = {}
    for cohort, X in X_by_cohort.items():
        available = [gene for gene in genes if gene in X.columns]
        ranked[cohort] = rank_gaussian_normalize(X[available].astype(float)) if available else pd.DataFrame(index=X.index)
    return ranked


def estimate_gene_directions(
    ranked_by_cohort: dict[str, pd.DataFrame],
    y_by_cohort: dict[str, pd.Series],
    train_cohorts: list[str],
) -> dict[str, int]:
    genes = module_genes()
    frames = []
    labels = []
    for cohort in train_cohorts:
        if cohort not in ranked_by_cohort or cohort not in y_by_cohort:
            continue
        X = ranked_by_cohort[cohort]
        y = y_by_cohort[cohort].astype(float)
        available = [gene for gene in genes if gene in X.columns]
        if available:
            frames.append(X.reindex(y.index).loc[:, available])
            labels.append(y)
    if not frames:
        return {gene: 1 for gene in genes}
    X_train = pd.concat(frames, axis=0)
    y_train = pd.concat(labels, axis=0).reindex(X_train.index)
    directions: dict[str, int] = {}
    for gene in genes:
        if gene not in X_train.columns:
            directions[gene] = 1
            continue
        values = pd.to_numeric(X_train[gene], errors="coerce")
        valid = values.notna() & y_train.notna()
        if int(valid.sum()) < 6 or values.loc[valid].nunique() < 2 or y_train.loc[valid].nunique() < 2:
            directions[gene] = 1
            continue
        corr = values.loc[valid].corr(y_train.loc[valid])
        directions[gene] = -1 if pd.notna(corr) and corr < 0 else 1
    return directions


def signed_rank_module_features(
    ranked_X: pd.DataFrame,
    directions: dict[str, int] | None = None,
    signed: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    directions = directions or {}
    feature_rows: dict[str, pd.Series] = {}
    coverage_rows: list[dict[str, object]] = []
    for module, genes in MODULE_GENE_SETS.items():
        available = [gene for gene in genes if gene in ranked_X.columns]
        if not available:
            feature_rows[module] = pd.Series(0.0, index=ranked_X.index)
        else:
            values = ranked_X[available].astype(float).copy()
            if signed:
                for gene in available:
                    values[gene] = values[gene] * int(directions.get(gene, 1))
            feature_rows[module] = values.sum(axis=1) / np.sqrt(len(available))
        coverage_rows.append(
            {
                "module": module,
                "n_genes_defined": len(genes),
                "n_genes_available": len(available),
                "genes_available": ",".join(available),
                "signed": bool(signed),
            }
        )
    return pd.DataFrame(feature_rows, index=ranked_X.index).fillna(0.0), pd.DataFrame(coverage_rows)


def score_modules(features: pd.DataFrame, weight_mode: str) -> pd.Series:
    if weight_mode == "prior":
        weights = MODULE_PRIOR_WEIGHTS
    elif weight_mode == "absolute_prior":
        weights = {key: abs(value) for key, value in MODULE_PRIOR_WEIGHTS.items()}
    else:
        raise ValueError(f"Unknown weight_mode: {weight_mode}")
    score = pd.Series(0.0, index=features.index)
    for module, weight in weights.items():
        if module in features.columns:
            score = score + float(weight) * features[module]
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
    prob = calibrator.predict_proba(score.astype(float).to_numpy().reshape(-1, 1))[:, 1]
    return pd.Series(prob, index=score.index)


def raw_baseline_scores(X_by_cohort: dict[str, pd.DataFrame]) -> dict[str, pd.Series]:
    module_features, _ = build_module_features_by_cohort(X_by_cohort)
    return {cohort: module_prior_score(features) for cohort, features in module_features.items()}


def make_train_test_scores(
    spec: dict[str, object],
    train_cohorts: list[str],
    test_cohort: str,
    ranked_by_cohort: dict[str, pd.DataFrame],
    y_by_cohort: dict[str, pd.Series],
    raw_scores: dict[str, pd.Series],
) -> tuple[pd.Series, pd.Series, pd.DataFrame]:
    if not spec.get("rank", False):
        train_score = pd.concat([raw_scores[cohort].reindex(y_by_cohort[cohort].index) for cohort in train_cohorts])
        test_score = raw_scores[test_cohort].reindex(y_by_cohort[test_cohort].index)
        return train_score, test_score, pd.DataFrame()
    directions = estimate_gene_directions(ranked_by_cohort, y_by_cohort, train_cohorts) if spec.get("gene_signed") else {}
    train_parts = []
    coverage = []
    for cohort in train_cohorts:
        frame, cov = signed_rank_module_features(
            ranked_by_cohort[cohort].reindex(y_by_cohort[cohort].index),
            directions=directions,
            signed=bool(spec.get("gene_signed")),
        )
        train_parts.append(score_modules(frame, str(spec["weight_mode"])).reindex(y_by_cohort[cohort].index))
        cov.insert(0, "cohort", cohort)
        coverage.append(cov)
    test_frame, test_cov = signed_rank_module_features(
        ranked_by_cohort[test_cohort].reindex(y_by_cohort[test_cohort].index),
        directions=directions,
        signed=bool(spec.get("gene_signed")),
    )
    test_cov.insert(0, "cohort", test_cohort)
    coverage.append(test_cov)
    train_score = pd.concat(train_parts)
    test_score = score_modules(test_frame, str(spec["weight_mode"])).reindex(y_by_cohort[test_cohort].index)
    coverage_df = pd.concat(coverage, ignore_index=True)
    coverage_df["model_name"] = spec["model_name"]
    return train_score, test_score, coverage_df


def evaluate_fold(
    spec: dict[str, object],
    train_cohorts: list[str],
    test_cohort: str,
    ranked_by_cohort: dict[str, pd.DataFrame],
    y_by_cohort: dict[str, pd.Series],
    raw_scores: dict[str, pd.Series],
) -> tuple[dict[str, object], list[dict[str, object]], pd.DataFrame]:
    y_train = pd.concat([y_by_cohort[cohort] for cohort in train_cohorts]).astype(int)
    y_test = y_by_cohort[test_cohort].astype(int)
    train_score, test_score, coverage = make_train_test_scores(spec, train_cohorts, test_cohort, ranked_by_cohort, y_by_cohort, raw_scores)
    train_score = train_score.reindex(y_train.index)
    test_score = test_score.reindex(y_test.index)
    calibrator = fit_monotone_platt(train_score, y_train) if spec.get("calibrate") else None
    train_prob = predict_prob(train_score, calibrator)
    test_prob = predict_prob(test_score, calibrator)
    threshold = select_threshold(y_train.to_numpy(dtype=int), train_prob.to_numpy(dtype=float))
    metrics = compute_binary_metrics(y_test, test_prob, threshold=threshold)
    metric_row = {
        **metrics,
        "model_name": spec["model_name"],
        "cohort": test_cohort,
        "n_samples": int(len(y_test)),
        "n_responders": int(y_test.sum()),
        "n_nonresponders": int((y_test == 0).sum()),
        "threshold": float(threshold),
        "calibration": "training_only_platt" if calibrator is not None else "raw_sigmoid",
        "train_cohorts": ",".join(train_cohorts),
        "selection_boundary": "training_only_no_external_labels",
    }
    prediction_rows = [
        {
            "model_name": spec["model_name"],
            "cohort": test_cohort,
            "sample_id": sample_id,
            "true_response_label": int(y_test.loc[sample_id]),
            "response_probability": float(test_prob.loc[sample_id]),
            "threshold": float(threshold),
        }
        for sample_id in y_test.index
    ]
    return metric_row, prediction_rows, coverage


def evaluate_primary_lodo(
    X_by_cohort: dict[str, pd.DataFrame],
    metadata_by_cohort: dict[str, pd.DataFrame],
    ranked_by_cohort: dict[str, pd.DataFrame],
    raw_scores: dict[str, pd.Series],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    active = [cohort for cohort in X_by_cohort if not cohort.startswith("demo_cohort_")]
    strata = default_strata(active)
    metric_rows = []
    prediction_rows = []
    coverage_rows = []
    for stratum in PRIMARY_STRATA:
        endpoint_data = prepare_endpoint_data(X_by_cohort, metadata_by_cohort, strata[stratum]["cohorts"], PRIMARY_ENDPOINT)
        y_by_cohort = endpoint_data.y_response_by_cohort
        specs = [{"model_name": BASELINE_MODEL, "rank": False, "calibrate": False}, *VARIANT_SPECS]
        for holdout in strata[stratum]["holdouts"]:
            if holdout not in y_by_cohort:
                continue
            train_cohorts = [cohort for cohort in strata[stratum]["train_pool"] if cohort != holdout and cohort in y_by_cohort]
            if len(train_cohorts) < 2:
                continue
            for spec in specs:
                metric, predictions, coverage = evaluate_fold(spec, train_cohorts, holdout, ranked_by_cohort, y_by_cohort, raw_scores)
                metric.update({"endpoint": PRIMARY_ENDPOINT, "stratum": stratum})
                metric_rows.append(metric)
                for row in predictions:
                    row.update({"endpoint": PRIMARY_ENDPOINT, "stratum": stratum})
                    prediction_rows.append(row)
                if not coverage.empty:
                    coverage["endpoint"] = PRIMARY_ENDPOINT
                    coverage["stratum"] = stratum
                    coverage["holdout"] = holdout
                    coverage_rows.append(coverage)
    return (
        pd.DataFrame(metric_rows),
        pd.DataFrame(prediction_rows),
        pd.concat(coverage_rows, ignore_index=True) if coverage_rows else pd.DataFrame(),
    )


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
    X_by_cohort: dict[str, pd.DataFrame],
    metadata_by_cohort: dict[str, pd.DataFrame],
    ranked_by_cohort: dict[str, pd.DataFrame],
    raw_scores: dict[str, pd.Series],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    needed = sorted(set(DISCOVERY_COHORTS + cohorts))
    y_by_cohort = endpoint_labels(metadata_by_cohort, needed, STRICT_ENDPOINT)
    train_cohorts = [cohort for cohort in DISCOVERY_COHORTS if cohort in y_by_cohort]
    specs = [{"model_name": BASELINE_MODEL, "rank": False, "calibrate": False}, *VARIANT_SPECS]
    metric_rows = []
    prediction_rows = []
    coverage_rows = []
    for holdout in cohorts:
        if holdout not in y_by_cohort:
            continue
        for spec in specs:
            metric, predictions, coverage = evaluate_fold(spec, train_cohorts, holdout, ranked_by_cohort, y_by_cohort, raw_scores)
            metric.update({"endpoint": STRICT_ENDPOINT, "stratum": group_id})
            metric_rows.append(metric)
            for row in predictions:
                row.update({"endpoint": STRICT_ENDPOINT, "stratum": group_id})
                prediction_rows.append(row)
            if not coverage.empty:
                coverage["endpoint"] = STRICT_ENDPOINT
                coverage["stratum"] = group_id
                coverage["holdout"] = holdout
                coverage_rows.append(coverage)
    return (
        pd.DataFrame(metric_rows),
        pd.DataFrame(prediction_rows),
        pd.concat(coverage_rows, ignore_index=True) if coverage_rows else pd.DataFrame(),
    )


def summarize_predictions(predictions: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if predictions.empty:
        return pd.DataFrame()
    for (endpoint, stratum, model), frame in predictions.groupby(["endpoint", "stratum", "model_name"]):
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
                "model_name": model,
                "n_samples": int(len(frame)),
                "n_cohorts": int(frame["cohort"].nunique()),
                "n_responders": int(y.sum()),
                "n_nonresponders": int((y == 0).sum()),
                "response_prevalence": float(y.mean()),
                "AUPRC_minus_prevalence": float(metrics["AUPRC"] - y.mean()),
            }
        )
    return pd.DataFrame(rows).sort_values(["endpoint", "stratum", "model_name"]).reset_index(drop=True)


def paired_comparisons(predictions: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (endpoint, stratum), frame in predictions.groupby(["endpoint", "stratum"]):
        baseline = frame[frame["model_name"].eq(BASELINE_MODEL)]
        if baseline.empty:
            continue
        for model, target in frame.groupby("model_name"):
            if model == BASELINE_MODEL:
                continue
            merged = target.merge(
                baseline[["cohort", "sample_id", "true_response_label", "response_probability"]],
                on=["cohort", "sample_id", "true_response_label"],
                suffixes=("_target", "_baseline"),
            )
            if len(merged) < 8 or merged["true_response_label"].nunique() < 2:
                continue
            y = merged["true_response_label"].astype(int)
            target_prob = merged["response_probability_target"].astype(float)
            baseline_prob = merged["response_probability_baseline"].astype(float)
            target_metrics = compute_binary_metrics(y, target_prob)
            baseline_metrics = compute_binary_metrics(y, baseline_prob)
            stats = paired_bootstrap_delta(y, target_prob, baseline_prob, n_bootstrap=1000, random_state=20260527)
            rows.append(
                {
                    "endpoint": endpoint,
                    "stratum": stratum,
                    "target_model": model,
                    "baseline_model": BASELINE_MODEL,
                    "n_samples": int(len(merged)),
                    "target_AUROC": float(target_metrics["AUROC"]),
                    "baseline_AUROC": float(baseline_metrics["AUROC"]),
                    "delta_AUROC": float(target_metrics["AUROC"] - baseline_metrics["AUROC"]),
                    "target_balanced_accuracy": float(target_metrics["balanced_accuracy"]),
                    "baseline_balanced_accuracy": float(baseline_metrics["balanced_accuracy"]),
                    "delta_balanced_accuracy": float(target_metrics["balanced_accuracy"] - baseline_metrics["balanced_accuracy"]),
                    "target_ECE": float(target_metrics["ECE"]),
                    "baseline_ECE": float(baseline_metrics["ECE"]),
                    "delta_ECE": float(target_metrics["ECE"] - baseline_metrics["ECE"]),
                    **stats,
                }
            )
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result["fdr_q"] = 1.0
    for _, idx in result.groupby(["endpoint", "stratum"]).groups.items():
        result.loc[idx, "fdr_q"] = benjamini_hochberg(result.loc[idx, "p_value"].fillna(1.0))
    result["claim_level"] = np.where(
        (result["delta_AUROC"] > 0) & (result["fdr_q"] <= 0.05),
        "FDR_supported_AUROC_gain",
        np.where(
            (result["delta_AUROC"] > 0)
            | (result["delta_balanced_accuracy"] > 0)
            | (result["delta_ECE"] < 0),
            "point_estimate_or_calibration_tradeoff",
            "not_supported",
        ),
    )
    return result.sort_values(["endpoint", "stratum", "target_model"]).reset_index(drop=True)


def write_audit(out_dir: Path, summary: pd.DataFrame, comparisons: pd.DataFrame) -> None:
    lines = [
        "# Signed-rank Module Audit",
        "",
        "This audit tests whether sample-wise rank Gaussian module scoring and training-only gene-direction signing improve the current raw module-prior score. External cohorts are used only for evaluation.",
        "",
        "## Summary",
        "",
    ]
    if summary.empty:
        lines.append("No summary rows were produced.")
    else:
        for _, row in summary.iterrows():
            lines.append(
                f"- {row['endpoint']} / {row['stratum']} / {row['model_name']}: AUROC={row['AUROC']:.3f}, "
                f"AUPRC={row['AUPRC']:.3f}, balanced accuracy={row['balanced_accuracy']:.3f}, ECE={row['ECE']:.3f}."
            )
    lines.extend(["", "## Baseline Comparisons", ""])
    if comparisons.empty:
        lines.append("No paired comparisons were produced.")
    else:
        for _, row in comparisons.iterrows():
            lines.append(
                f"- {row['endpoint']} / {row['stratum']} / {row['target_model']}: dAUROC={row['delta_AUROC']:.3f}, "
                f"dBA={row['delta_balanced_accuracy']:.3f}, dECE={row['delta_ECE']:.3f}, q={row['fdr_q']:.3f} ({row['claim_level']})."
            )
    (out_dir / "SIGNED_RANK_MODULE_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--cbio-dir", default="data/processed/cbioportal_melanoma")
    parser.add_argument("--out", default="results/signed_rank_module_audit_20260527")
    args = parser.parse_args()

    out_dir = ROOT / args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    X_by_cohort, _, metadata_by_cohort = merge_processed_dirs(ROOT / args.processed_dir, ROOT / args.cbio_dir)
    ranked = rank_by_cohort(X_by_cohort)
    raw_scores = raw_baseline_scores(X_by_cohort)

    primary_metrics, primary_predictions, primary_coverage = evaluate_primary_lodo(
        X_by_cohort, metadata_by_cohort, ranked, raw_scores
    )
    external_frames = []
    prediction_frames = [primary_predictions]
    coverage_frames = [primary_coverage]
    for group_id, cohorts in {
        "strict_melanoma_pd1_like_external": STRICT_EXTERNAL_COHORTS,
        "strict_cbio_liu_plus_gse145996": CBIO_EXTERNAL_COHORTS,
    }.items():
        metrics, predictions, coverage = evaluate_external_group(
            group_id, cohorts, X_by_cohort, metadata_by_cohort, ranked, raw_scores
        )
        external_frames.append(metrics)
        prediction_frames.append(predictions)
        coverage_frames.append(coverage)
    metrics = pd.concat([primary_metrics, *external_frames], ignore_index=True)
    predictions = pd.concat(prediction_frames, ignore_index=True)
    coverage = pd.concat([frame for frame in coverage_frames if not frame.empty], ignore_index=True)
    summary = summarize_predictions(predictions)
    comparisons = paired_comparisons(predictions)

    metrics.to_csv(out_dir / "signed_rank_module_metrics.tsv", sep="\t", index=False)
    predictions.to_csv(out_dir / "signed_rank_module_predictions.tsv", sep="\t", index=False)
    coverage.to_csv(out_dir / "signed_rank_module_coverage.tsv", sep="\t", index=False)
    summary.to_csv(out_dir / "signed_rank_module_summary.tsv", sep="\t", index=False)
    comparisons.to_csv(out_dir / "signed_rank_module_comparisons.tsv", sep="\t", index=False)
    write_audit(out_dir, summary, comparisons)
    print(f"Wrote signed-rank module audit outputs to {out_dir}")


if __name__ == "__main__":
    main()
