from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.io import load_processed_bulk
from econiche.metrics import compute_binary_metrics
from econiche.statistics import benjamini_hochberg
from econiche_opt.model.endpoint_modules import (
    STRONG_BASELINES,
    build_fixed_scores_by_cohort,
    build_module_features_by_cohort,
    endpoint_label_series,
    select_threshold,
    sigmoid,
)


DISCOVERY_COHORTS = ["GSE91061", "GSE78220", "PRJEB23709_PD1_PRE"]
TARGET_MODEL = "EcoNiche-Opt-HeuristicEcology-LockedPanel"
MODULE_PRIOR_MODEL = "EcoNiche-Opt-ModulePriorFixed"
EIGHT_SIGNATURES = ["IFNG", "CXCL9", "TIG", "TIDE_dysfunction", "TIDE_exclusion", "CYT", "APM", "IPRES"]
MODEL_NAMES = [TARGET_MODEL, *EIGHT_SIGNATURES]

EVALUATION_GROUPS: dict[str, dict[str, object]] = {
    "cbio_liu_dfci_only": {
        "cohorts": ["CBIO_LIU_DFCI_2019_PRE"],
        "claim_status": "independent_cbioportal_liu_external",
        "notes": "cBioPortal original Liu/DFCI 2019 pretreatment anti-PD1 RNA-seq; use as an independent-source replacement/cross-check for TIGER phs000452, not together with duplicate Liu versions.",
    },
    "strict_cbio_liu_plus_gse145996": {
        "cohorts": ["CBIO_LIU_DFCI_2019_PRE", "GSE145996"],
        "claim_status": "strict_melanoma_pd1_like_external_cbioportal",
        "notes": "Strict external set combining cBioPortal Liu/DFCI and GSE145996; no discovery, threshold, calibration, or feature-selection labels are taken from these cohorts.",
    },
    "current_tiger_phs_plus_gse145996": {
        "cohorts": ["PHS000452_LIU_LIKE_PRE", "GSE145996"],
        "claim_status": "current_tiger_strict_external_reference",
        "notes": "Existing TIGER phs000452-derived strict external reference; included to quantify whether cBioPortal processing changes the locked score.",
    },
    "cbio_iatlas_liu_duplicate_crosscheck": {
        "cohorts": ["CBIO_IATLAS_LIU_2019_PRE"],
        "claim_status": "duplicate_source_crosscheck_not_independent",
        "notes": "iAtlas Liu harmonized cohort shares source patients with Liu/DFCI; use only as a duplicate-source processing cross-check.",
    },
    "cbio_discovery_overlap_crosscheck": {
        "cohorts": ["CBIO_IATLAS_GIDE_2019_PRE", "CBIO_IATLAS_RIAZ_2017_PRE", "CBIO_IATLAS_HUGO_2016_PRE"],
        "claim_status": "discovery_overlap_crosscheck_not_external",
        "notes": "cBioPortal versions of Gide, Riaz and Hugo overlap discovery cohorts and are not independent external validation.",
    },
}


def merge_processed_dirs(*processed_dirs: Path) -> tuple[dict[str, pd.DataFrame], dict[str, pd.Series], dict[str, pd.DataFrame]]:
    merged_X: dict[str, pd.DataFrame] = {}
    merged_y: dict[str, pd.Series] = {}
    merged_meta: dict[str, pd.DataFrame] = {}
    for processed_dir in processed_dirs:
        X, y, meta = load_processed_bulk(processed_dir)
        merged_X.update(X)
        merged_y.update(y)
        merged_meta.update(meta)
    return merged_X, merged_y, merged_meta


def labels_for_endpoint(metadata_by_cohort: dict[str, pd.DataFrame], cohorts: list[str], endpoint: str) -> dict[str, pd.Series]:
    labels: dict[str, pd.Series] = {}
    for cohort in cohorts:
        meta = metadata_by_cohort.get(cohort)
        if meta is None or "response_raw" not in meta.columns:
            continue
        y = endpoint_label_series(meta["response_raw"], endpoint).dropna().astype(int)
        if len(y) >= 4 and y.nunique() == 2:
            labels[cohort] = y
    return labels


def fit_monotone_platt(score: pd.Series, y: pd.Series) -> LogisticRegression | None:
    common = score.index.intersection(y.index)
    if len(common) < 8 or y.loc[common].nunique() < 2:
        return None
    model = LogisticRegression(C=1.0, class_weight="balanced", solver="lbfgs", max_iter=5000, random_state=20260527)
    model.fit(score.loc[common].astype(float).to_numpy().reshape(-1, 1), y.loc[common].astype(int).to_numpy())
    if float(model.coef_[0, 0]) <= 0:
        return None
    return model


def predict_with_calibrator(score: pd.Series, calibrator: LogisticRegression | None) -> pd.Series:
    if calibrator is None:
        return pd.Series(sigmoid(score), index=score.index)
    prob = calibrator.predict_proba(score.astype(float).to_numpy().reshape(-1, 1))[:, 1]
    return pd.Series(prob, index=score.index)


def build_scores_by_cohort(X_by_cohort: dict[str, pd.DataFrame]) -> tuple[dict[str, dict[str, pd.Series]], pd.DataFrame]:
    module_features, coverage = build_module_features_by_cohort(X_by_cohort)
    fixed = build_fixed_scores_by_cohort(
        X_by_cohort,
        module_features,
        baselines=sorted(set([*STRONG_BASELINES, *EIGHT_SIGNATURES])),
    )
    out: dict[str, dict[str, pd.Series]] = {}
    for cohort, scores in fixed.items():
        cohort_scores = {TARGET_MODEL: scores[MODULE_PRIOR_MODEL]}
        for model in EIGHT_SIGNATURES:
            if model in scores:
                cohort_scores[model] = scores[model]
        out[cohort] = cohort_scores
    return out, coverage


def fit_discovery_thresholds(
    scores_by_cohort: dict[str, dict[str, pd.Series]],
    labels_by_cohort: dict[str, pd.Series],
) -> dict[str, dict[str, object]]:
    thresholds: dict[str, dict[str, object]] = {}
    for model in MODEL_NAMES:
        score_parts = []
        y_parts = []
        train_cohorts = []
        for cohort in DISCOVERY_COHORTS:
            if cohort not in scores_by_cohort or cohort not in labels_by_cohort or model not in scores_by_cohort[cohort]:
                continue
            y = labels_by_cohort[cohort].astype(int)
            score = scores_by_cohort[cohort][model].reindex(y.index).dropna()
            common = y.index.intersection(score.index)
            if len(common) < 4:
                continue
            score_parts.append(score.loc[common])
            y_parts.append(y.loc[common])
            train_cohorts.append(cohort)
        if not score_parts:
            continue
        score_train = pd.concat(score_parts)
        y_train = pd.concat(y_parts).reindex(score_train.index).astype(int)
        calibrator = fit_monotone_platt(score_train, y_train)
        prob_train = predict_with_calibrator(score_train, calibrator)
        threshold = select_threshold(y_train.to_numpy(dtype=int), prob_train.to_numpy(dtype=float))
        train_metrics = compute_binary_metrics(y_train, prob_train, threshold=threshold)
        thresholds[model] = {
            "threshold": float(threshold),
            "calibrator": calibrator,
            "threshold_training_cohorts": ",".join(train_cohorts),
            "training_AUROC": float(train_metrics["AUROC"]),
            "training_balanced_accuracy": float(train_metrics["balanced_accuracy"]),
            "calibration": "discovery_only_monotone_platt" if calibrator is not None else "raw_sigmoid",
        }
    return thresholds


def make_predictions(
    endpoint: str,
    group_id: str,
    cohorts: list[str],
    scores_by_cohort: dict[str, dict[str, pd.Series]],
    labels_by_cohort: dict[str, pd.Series],
    metadata_by_cohort: dict[str, pd.DataFrame],
    thresholds: dict[str, dict[str, object]],
    claim_status: str,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for cohort in cohorts:
        if cohort not in labels_by_cohort or cohort not in scores_by_cohort:
            continue
        y = labels_by_cohort[cohort].astype(int)
        meta = metadata_by_cohort[cohort].reindex(y.index)
        for model in MODEL_NAMES:
            if model not in thresholds or model not in scores_by_cohort[cohort]:
                continue
            score = scores_by_cohort[cohort][model].reindex(y.index).dropna()
            common = y.index.intersection(score.index)
            if len(common) < 4:
                continue
            prob = predict_with_calibrator(score.loc[common], thresholds[model]["calibrator"])
            threshold = float(thresholds[model]["threshold"])
            for sample_id in common:
                rows.append(
                    {
                        "endpoint": endpoint,
                        "group_id": group_id,
                        "cohort": cohort,
                        "claim_status": claim_status,
                        "model_name": model,
                        "sample_id": sample_id,
                        "true_response_label": int(y.loc[sample_id]),
                        "response_raw": meta.loc[sample_id, "response_raw"] if "response_raw" in meta.columns else "",
                        "response_probability": float(prob.loc[sample_id]),
                        "threshold": threshold,
                        "predicted_response_label": int(prob.loc[sample_id] >= threshold),
                        "threshold_training_cohorts": thresholds[model]["threshold_training_cohorts"],
                        "calibration": thresholds[model]["calibration"],
                    }
                )
    return pd.DataFrame(rows)


def summarize_metrics(predictions: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if predictions.empty:
        return pd.DataFrame()
    group_cols = ["endpoint", "group_id", "claim_status", "model_name"]
    for keys, frame in predictions.groupby(group_cols):
        y = frame["true_response_label"].astype(int)
        if len(frame) < 8 or y.nunique() < 2:
            continue
        p = frame["response_probability"].astype(float)
        threshold = float(frame["threshold"].iloc[0])
        metrics = compute_binary_metrics(y, p, threshold=threshold)
        rows.append(
            {
                **dict(zip(group_cols, keys)),
                "n_samples": int(len(frame)),
                "n_cohorts": int(frame["cohort"].nunique()),
                "n_responders": int(y.sum()),
                "n_nonresponders": int((y == 0).sum()),
                "response_prevalence": float(y.mean()),
                **metrics,
            }
        )
    return pd.DataFrame(rows).sort_values(["endpoint", "group_id", "model_name"]).reset_index(drop=True)


def _bootstrap_family_delta(y: pd.Series, target: pd.Series, baselines: pd.DataFrame, n_bootstrap: int = 2000) -> dict[str, float]:
    y_arr = y.to_numpy(dtype=int)
    target_arr = target.to_numpy(dtype=float)
    baseline_arr = baselines.to_numpy(dtype=float)
    rng = np.random.default_rng(20260527)
    deltas: list[float] = []
    for _ in range(n_bootstrap):
        idx = rng.integers(0, len(y_arr), len(y_arr))
        if len(np.unique(y_arr[idx])) < 2:
            continue
        target_auc = roc_auc_score(y_arr[idx], target_arr[idx])
        baseline_auc = np.mean([roc_auc_score(y_arr[idx], baseline_arr[idx, j]) for j in range(baseline_arr.shape[1])])
        deltas.append(float(target_auc - baseline_auc))
    if not deltas:
        return {}
    arr = np.asarray(deltas)
    return {
        "delta_vs_family_mean": float(arr.mean()),
        "ci_low": float(np.quantile(arr, 0.025)),
        "ci_high": float(np.quantile(arr, 0.975)),
        "one_sided_p": float((arr <= 0).mean()),
        "two_sided_p": float(min(1.0, 2.0 * min((arr <= 0).mean(), (arr >= 0).mean()))),
    }


def family_comparison(predictions: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (endpoint, group_id, claim_status), frame in predictions.groupby(["endpoint", "group_id", "claim_status"]):
        target = frame[frame["model_name"].eq(TARGET_MODEL)]
        if target.empty:
            continue
        target_key = target["cohort"].astype(str) + "::" + target["sample_id"].astype(str)
        target_prob = pd.Series(target["response_probability"].to_numpy(dtype=float), index=target_key)
        y = pd.Series(target["true_response_label"].to_numpy(dtype=int), index=target_key)
        baseline_series = {}
        for baseline in EIGHT_SIGNATURES:
            base = frame[frame["model_name"].eq(baseline)]
            if base.empty:
                continue
            key = base["cohort"].astype(str) + "::" + base["sample_id"].astype(str)
            baseline_series[baseline] = pd.Series(base["response_probability"].to_numpy(dtype=float), index=key)
        baseline_frame = pd.DataFrame(baseline_series).dropna(axis=0)
        common = baseline_frame.index.intersection(target_prob.index).intersection(y.index)
        if len(common) < 8 or y.loc[common].nunique() < 2 or baseline_frame.shape[1] < 4:
            continue
        target_auc = float(roc_auc_score(y.loc[common], target_prob.loc[common]))
        baseline_aucs = {name: float(roc_auc_score(y.loc[common], baseline_frame.loc[common, name])) for name in baseline_frame.columns}
        stats = _bootstrap_family_delta(y.loc[common], target_prob.loc[common], baseline_frame.loc[common])
        rows.append(
            {
                "endpoint": endpoint,
                "group_id": group_id,
                "claim_status": claim_status,
                "target_model": TARGET_MODEL,
                "n_samples": int(len(common)),
                "n_signatures": int(baseline_frame.shape[1]),
                "target_AUROC": target_auc,
                "family_mean_AUROC": float(np.mean(list(baseline_aucs.values()))),
                "best_signature": max(baseline_aucs, key=baseline_aucs.get),
                "best_signature_AUROC": float(max(baseline_aucs.values())),
                **stats,
            }
        )
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result["two_sided_fdr_q"] = benjamini_hochberg(result["two_sided_p"].fillna(1.0))
    result["claim_level"] = np.where(
        (result["delta_vs_family_mean"] > 0) & (result["two_sided_fdr_q"] <= 0.05),
        "family_two_sided_FDR_supported",
        np.where(result["delta_vs_family_mean"] > 0, "family_point_estimate_only", "family_not_superior"),
    )
    return result.sort_values(["endpoint", "group_id"]).reset_index(drop=True)


def write_audit(out_dir: Path, metrics: pd.DataFrame, comparison: pd.DataFrame, coverage: pd.DataFrame) -> None:
    lines = [
        "# cBioPortal Melanoma External Validation Audit",
        "",
        "This audit adds public cBioPortal melanoma ICB expression profiles as a real-data source. Discovery-only cohorts are fixed to GSE91061, GSE78220 and PRJEB23709_PD1_PRE for calibration and thresholding. cBioPortal external labels are used only for evaluation.",
        "",
        "## Group Results",
        "",
    ]
    target = metrics[metrics["model_name"].eq(TARGET_MODEL)] if not metrics.empty else pd.DataFrame()
    if target.empty:
        lines.append("No target-model metrics were produced.")
    else:
        for _, row in target.iterrows():
            lines.append(
                f"- {row['endpoint']} / {row['group_id']}: n={int(row['n_samples'])}, AUROC={row['AUROC']:.3f}, "
                f"AUPRC={row['AUPRC']:.3f}, balanced accuracy={row['balanced_accuracy']:.3f}, ECE={row['ECE']:.3f} "
                f"({row['claim_status']})."
            )
    lines.extend(["", "## Family Claim Gate", ""])
    if comparison.empty:
        lines.append("No eight-signature family comparison was produced.")
    else:
        for _, row in comparison.iterrows():
            lines.append(
                f"- {row['endpoint']} / {row['group_id']}: target AUROC={row['target_AUROC']:.3f}, "
                f"family mean AUROC={row['family_mean_AUROC']:.3f}, delta={row['delta_vs_family_mean']:.3f}, "
                f"q={row['two_sided_fdr_q']:.3f}, best signature={row['best_signature']} ({row['claim_level']})."
            )
    lines.extend(["", "## Coverage", ""])
    if not coverage.empty:
        cohort_coverage = coverage.groupby("cohort")["n_genes_available"].sum().reset_index()
        for _, row in cohort_coverage.iterrows():
            lines.append(f"- {row['cohort']}: total available module-gene hits={int(row['n_genes_available'])}.")
    (out_dir / "CBIOPORTAL_EXTERNAL_VALIDATION_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--cbio-dir", default="data/processed/cbioportal_melanoma")
    parser.add_argument("--out", default="results/cbioportal_melanoma_external_validation_20260527")
    parser.add_argument("--endpoints", nargs="*", default=["strict_recist", "primary_recist"])
    args = parser.parse_args()

    out_dir = ROOT / args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    X_by_cohort, _, metadata_by_cohort = merge_processed_dirs(ROOT / args.processed_dir, ROOT / args.cbio_dir)
    scores_by_cohort, coverage = build_scores_by_cohort(X_by_cohort)
    coverage.to_csv(out_dir / "cbioportal_module_coverage.tsv", sep="\t", index=False)

    all_predictions = []
    all_metrics = []
    all_comparisons = []
    for endpoint in args.endpoints:
        needed_cohorts = sorted(set(DISCOVERY_COHORTS + [cohort for group in EVALUATION_GROUPS.values() for cohort in group["cohorts"]]))
        labels_by_cohort = labels_for_endpoint(metadata_by_cohort, needed_cohorts, endpoint)
        thresholds = fit_discovery_thresholds(scores_by_cohort, labels_by_cohort)
        for group_id, group in EVALUATION_GROUPS.items():
            predictions = make_predictions(
                endpoint,
                group_id,
                list(group["cohorts"]),
                scores_by_cohort,
                labels_by_cohort,
                metadata_by_cohort,
                thresholds,
                str(group["claim_status"]),
            )
            if not predictions.empty:
                predictions["group_notes"] = str(group["notes"])
                all_predictions.append(predictions)
                metrics = summarize_metrics(predictions)
                if not metrics.empty:
                    metrics["group_notes"] = str(group["notes"])
                    all_metrics.append(metrics)
                comparison = family_comparison(predictions)
                if not comparison.empty:
                    comparison["group_notes"] = str(group["notes"])
                    all_comparisons.append(comparison)
    predictions_df = pd.concat(all_predictions, ignore_index=True) if all_predictions else pd.DataFrame()
    metrics_df = pd.concat(all_metrics, ignore_index=True) if all_metrics else pd.DataFrame()
    comparison_df = pd.concat(all_comparisons, ignore_index=True) if all_comparisons else pd.DataFrame()
    predictions_df.to_csv(out_dir / "cbioportal_external_predictions.tsv", sep="\t", index=False)
    metrics_df.to_csv(out_dir / "cbioportal_external_metrics.tsv", sep="\t", index=False)
    comparison_df.to_csv(out_dir / "cbioportal_external_family_comparison.tsv", sep="\t", index=False)
    write_audit(out_dir, metrics_df, comparison_df, coverage)
    print(f"Wrote cBioPortal external validation outputs to {out_dir}")


if __name__ == "__main__":
    main()
