from __future__ import annotations

import argparse
import copy
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.io import load_processed_bulk
from econiche.metrics import compute_binary_metrics
from econiche.normalize import rank_gaussian_normalize
from econiche_opt.model.endpoint_modules import (
    MODULE_GENE_SETS,
    WORD_STATE_GENE_SETS,
    build_fixed_scores_by_cohort,
    build_module_features_by_cohort,
    endpoint_label_series,
    select_threshold,
    sigmoid,
)


DISCOVERY_COHORTS = ["GSE91061", "GSE78220", "PRJEB23709_PD1_PRE"]
PD1_LIKE_EXTERNAL_COHORTS = ["GSE145996", "PHS000452_LIU_LIKE_PRE"]
TARGET_LOCKED = "EcoNiche-Opt-HeuristicEcology-LockedPanel"
TRANSFER_HEAD = "EcoNiche-Opt-PD1LikeTransferHead"
BASELINE_SIGNATURES = ["IFNG", "CXCL9", "TIG", "TIDE_dysfunction", "APM", "CYT", "IPRES", "TIDE_exclusion"]


def _gene_universe(X_by_cohort: dict[str, pd.DataFrame], cohorts: list[str]) -> list[str]:
    genes: set[str] = set()
    for module_genes in MODULE_GENE_SETS.values():
        genes.update(module_genes)
    for state_genes in WORD_STATE_GENE_SETS.values():
        genes.update(state_genes)
    return [gene for gene in sorted(genes) if all(gene in X_by_cohort[cohort].columns for cohort in cohorts)]


def _labels(metadata_by_cohort: dict[str, pd.DataFrame], cohorts: list[str]) -> dict[str, pd.Series]:
    labels: dict[str, pd.Series] = {}
    for cohort in cohorts:
        y = endpoint_label_series(metadata_by_cohort[cohort]["response_raw"], "strict_recist").dropna().astype(int)
        if len(y) >= 4 and y.nunique() == 2:
            labels[cohort] = y
    return labels


def _build_tables(
    X_by_cohort: dict[str, pd.DataFrame],
    labels_by_cohort: dict[str, pd.Series],
    module_features_by_cohort: dict[str, pd.DataFrame],
    fixed_scores_by_cohort: dict[str, dict[str, pd.Series]],
    cohorts: list[str],
) -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame], dict[str, pd.Series]]:
    genes = _gene_universe(X_by_cohort, cohorts)
    ranked = {cohort: rank_gaussian_normalize(X_by_cohort[cohort][genes].astype(float)) for cohort in cohorts}
    module_tables: dict[str, pd.DataFrame] = {}
    gene_tables: dict[str, pd.DataFrame] = {}
    locked_scores: dict[str, pd.Series] = {}
    for cohort in cohorts:
        y = labels_by_cohort[cohort]
        module_tables[cohort] = module_features_by_cohort[cohort].reindex(y.index).fillna(0.0)
        gene_tables[cohort] = ranked[cohort].reindex(y.index).fillna(0.0)
        locked_scores[cohort] = pd.Series(
            sigmoid(fixed_scores_by_cohort[cohort]["EcoNiche-Opt-ModulePriorFixed"].reindex(y.index)),
            index=y.index,
        )
    return module_tables, gene_tables, locked_scores


def _fit_transfer_components(
    train_cohorts: list[str],
    module_tables: dict[str, pd.DataFrame],
    gene_tables: dict[str, pd.DataFrame],
    labels_by_cohort: dict[str, pd.Series],
) -> tuple[ExtraTreesClassifier, object]:
    X_gene = pd.concat([gene_tables[cohort] for cohort in train_cohorts])
    X_module = pd.concat([module_tables[cohort] for cohort in train_cohorts])
    y = pd.concat([labels_by_cohort[cohort].reindex(gene_tables[cohort].index) for cohort in train_cohorts])
    gene_model = ExtraTreesClassifier(
        n_estimators=500,
        max_depth=3,
        min_samples_leaf=3,
        class_weight="balanced",
        random_state=20260508,
        n_jobs=-1,
    )
    module_model = make_pipeline(
        SelectKBest(f_classif, k=min(5, X_module.shape[1])),
        StandardScaler(),
        LogisticRegression(C=0.1, solver="liblinear", class_weight="balanced", random_state=20260508),
    )
    gene_model.fit(X_gene, y)
    module_model.fit(X_module, y)
    return gene_model, module_model


def _predict_transfer(
    gene_model: ExtraTreesClassifier,
    module_model: object,
    cohort: str,
    module_tables: dict[str, pd.DataFrame],
    gene_tables: dict[str, pd.DataFrame],
) -> pd.Series:
    gene_prob = gene_model.predict_proba(gene_tables[cohort])[:, 1]
    module_prob = module_model.predict_proba(module_tables[cohort])[:, 1]
    prob = 0.25 * gene_prob + 0.75 * module_prob
    return pd.Series(prob, index=module_tables[cohort].index)


def _inner_lodo(
    module_tables: dict[str, pd.DataFrame],
    gene_tables: dict[str, pd.DataFrame],
    labels_by_cohort: dict[str, pd.Series],
    locked_scores: dict[str, pd.Series],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for model_name in [TARGET_LOCKED, TRANSFER_HEAD]:
        fold_aucs: list[float] = []
        for holdout in DISCOVERY_COHORTS:
            train = [cohort for cohort in DISCOVERY_COHORTS if cohort != holdout]
            y = labels_by_cohort[holdout]
            if model_name == TARGET_LOCKED:
                prob = locked_scores[holdout].reindex(y.index)
            else:
                gene_model, module_model = _fit_transfer_components(train, module_tables, gene_tables, labels_by_cohort)
                prob = _predict_transfer(gene_model, module_model, holdout, module_tables, gene_tables).reindex(y.index)
            fold_aucs.append(float(roc_auc_score(y, prob)))
        rows.append(
            {
                "endpoint": "strict_recist",
                "model_name": model_name,
                "selection_scope": "discovery_lodo_only",
                "inner_mean_AUROC": float(np.mean(fold_aucs)),
                "inner_min_AUROC": float(np.min(fold_aucs)),
                "inner_sd_AUROC": float(np.std(fold_aucs)),
                "selected_as_primary_locked": model_name == TARGET_LOCKED,
                "claim_status": "locked_primary" if model_name == TARGET_LOCKED else "secondary_model_development_rescue",
            }
        )
    return pd.DataFrame(rows)


def _metrics_and_predictions(
    labels_by_cohort: dict[str, pd.Series],
    locked_scores: dict[str, pd.Series],
    transfer_scores: dict[str, pd.Series],
    thresholds: dict[str, float],
    cohorts: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    metric_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    score_by_model = {TARGET_LOCKED: locked_scores, TRANSFER_HEAD: transfer_scores}
    for model_name, scores in score_by_model.items():
        pooled_y_parts = []
        pooled_p_parts = []
        for cohort in cohorts:
            y = labels_by_cohort[cohort]
            prob = scores[cohort].reindex(y.index)
            metrics = compute_binary_metrics(y, prob, threshold=thresholds[model_name])
            metric_rows.append(
                {
                    "endpoint": "strict_recist",
                    "cohort": cohort,
                    "model_name": model_name,
                    "analysis_type": "pd1_like_external_rescue",
                    "claim_status": "locked_primary" if model_name == TARGET_LOCKED else "secondary_model_development_rescue",
                    "n_samples": int(len(y)),
                    "n_responders": int(y.sum()),
                    "n_nonresponders": int((y == 0).sum()),
                    "threshold": thresholds[model_name],
                    **metrics,
                }
            )
            pooled_y_parts.append(y)
            pooled_p_parts.append(prob)
            for sample_id in y.index:
                prediction_rows.append(
                    {
                        "endpoint": "strict_recist",
                        "cohort": cohort,
                        "model_name": model_name,
                        "sample_id": sample_id,
                        "true_response_label": int(y.loc[sample_id]),
                        "response_probability": float(prob.loc[sample_id]),
                        "predicted_response_label": int(prob.loc[sample_id] >= thresholds[model_name]),
                        "claim_status": "locked_primary" if model_name == TARGET_LOCKED else "secondary_model_development_rescue",
                    }
                )
        pooled_y = pd.concat(pooled_y_parts)
        pooled_prob = pd.concat(pooled_p_parts).reindex(pooled_y.index)
        pooled_metrics = compute_binary_metrics(pooled_y, pooled_prob, threshold=thresholds[model_name])
        metric_rows.append(
            {
                "endpoint": "strict_recist",
                "cohort": "GSE145996+PHS000452_LIU_LIKE_PRE",
                "model_name": model_name,
                "analysis_type": "pd1_like_external_rescue_pooled",
                "claim_status": "locked_primary" if model_name == TARGET_LOCKED else "secondary_model_development_rescue",
                "n_samples": int(len(pooled_y)),
                "n_responders": int(pooled_y.sum()),
                "n_nonresponders": int((pooled_y == 0).sum()),
                "threshold": thresholds[model_name],
                **pooled_metrics,
            }
        )
    return pd.DataFrame(metric_rows), pd.DataFrame(prediction_rows)


def _threshold_sensitivity(
    labels_by_cohort: dict[str, pd.Series],
    locked_scores: dict[str, pd.Series],
    transfer_scores: dict[str, pd.Series],
    training_thresholds: dict[str, float],
    cohorts: list[str],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    score_by_model = {TARGET_LOCKED: locked_scores, TRANSFER_HEAD: transfer_scores}
    for model_name, scores in score_by_model.items():
        threshold_policies = {
            "discovery_youden": training_thresholds[model_name],
            "fixed_probability_0_5": 0.5,
        }
        for policy, threshold in threshold_policies.items():
            for cohort_group, group_cohorts in {
                "GSE145996": ["GSE145996"],
                "PHS000452_LIU_LIKE_PRE": ["PHS000452_LIU_LIKE_PRE"],
                "GSE145996+PHS000452_LIU_LIKE_PRE": cohorts,
            }.items():
                y = pd.concat([labels_by_cohort[cohort] for cohort in group_cohorts if cohort in labels_by_cohort])
                prob = pd.concat([scores[cohort].reindex(labels_by_cohort[cohort].index) for cohort in group_cohorts if cohort in labels_by_cohort])
                prob = prob.reindex(y.index)
                metrics = compute_binary_metrics(y, prob, threshold=float(threshold))
                rows.append(
                    {
                        "endpoint": "strict_recist",
                        "cohort": cohort_group,
                        "model_name": model_name,
                        "threshold_policy": policy,
                        "threshold": float(threshold),
                        "claim_status": "locked_primary"
                        if model_name == TARGET_LOCKED
                        else "secondary_model_development_rescue",
                        "n_samples": int(len(y)),
                        "n_responders": int(y.sum()),
                        "n_nonresponders": int((y == 0).sum()),
                        **metrics,
                    }
                )
    return pd.DataFrame(rows)


def run(processed_dir: Path, out_dir: Path) -> None:
    cohorts = DISCOVERY_COHORTS + PD1_LIKE_EXTERNAL_COHORTS
    X_by_cohort, _, metadata_by_cohort = load_processed_bulk(processed_dir)
    X_by_cohort = {cohort: X_by_cohort[cohort] for cohort in cohorts if cohort in X_by_cohort}
    cohorts = [cohort for cohort in cohorts if cohort in X_by_cohort]
    labels_by_cohort = _labels(metadata_by_cohort, cohorts)
    cohorts = [cohort for cohort in cohorts if cohort in labels_by_cohort]
    module_features_by_cohort, _ = build_module_features_by_cohort(X_by_cohort)
    fixed_scores_by_cohort = build_fixed_scores_by_cohort(X_by_cohort, module_features_by_cohort, baselines=BASELINE_SIGNATURES)
    module_tables, gene_tables, locked_scores = _build_tables(
        X_by_cohort,
        labels_by_cohort,
        module_features_by_cohort,
        fixed_scores_by_cohort,
        cohorts,
    )
    inner = _inner_lodo(module_tables, gene_tables, labels_by_cohort, locked_scores)
    gene_model, module_model = _fit_transfer_components(DISCOVERY_COHORTS, module_tables, gene_tables, labels_by_cohort)
    transfer_scores = {
        cohort: _predict_transfer(gene_model, module_model, cohort, module_tables, gene_tables) for cohort in cohorts
    }
    thresholds = {
        TARGET_LOCKED: select_threshold(
            pd.concat([labels_by_cohort[cohort] for cohort in DISCOVERY_COHORTS]).to_numpy(dtype=int),
            pd.concat([locked_scores[cohort] for cohort in DISCOVERY_COHORTS]).to_numpy(dtype=float),
        ),
        TRANSFER_HEAD: select_threshold(
            pd.concat([labels_by_cohort[cohort] for cohort in DISCOVERY_COHORTS]).to_numpy(dtype=int),
            pd.concat([transfer_scores[cohort] for cohort in DISCOVERY_COHORTS]).to_numpy(dtype=float),
        ),
    }
    external_cohorts = [cohort for cohort in PD1_LIKE_EXTERNAL_COHORTS if cohort in labels_by_cohort]
    metrics, predictions = _metrics_and_predictions(labels_by_cohort, locked_scores, transfer_scores, thresholds, external_cohorts)
    threshold_sensitivity = _threshold_sensitivity(labels_by_cohort, locked_scores, transfer_scores, thresholds, external_cohorts)
    out_dir.mkdir(parents=True, exist_ok=True)
    inner.to_csv(out_dir / "pd1_like_rescue_candidate_selection.tsv", sep="\t", index=False)
    metrics.to_csv(out_dir / "pd1_like_rescue_metrics.tsv", sep="\t", index=False)
    predictions.to_csv(out_dir / "pd1_like_rescue_predictions.tsv", sep="\t", index=False)
    threshold_sensitivity.to_csv(out_dir / "pd1_like_rescue_threshold_sensitivity.tsv", sep="\t", index=False)
    pooled = metrics[metrics["cohort"] == "GSE145996+PHS000452_LIU_LIKE_PRE"].set_index("model_name")
    lines = [
        "# PD1-Like External Rescue Audit",
        "",
        "This analysis addresses the weak strict melanoma PD1-like external performance for GSE145996 and PHS000452_LIU_LIKE_PRE.",
        "",
        "Important boundary: the locked primary model remains the valid external-validation model. The transfer head is a secondary model-development rescue introduced after seeing the weakness, so it needs a fresh independent locked validation cohort before it can support a primary superiority claim.",
        "",
        "## Strict RECIST Correction",
        "",
        "Minor response (MR) is excluded from strict RECIST rather than counted as CR/PR response. This is an endpoint-rule correction, not model tuning.",
        "",
        "## Pooled External Performance",
        "",
    ]
    for model_name in [TARGET_LOCKED, TRANSFER_HEAD]:
        if model_name in pooled.index:
            row = pooled.loc[model_name]
            lines.append(
                f"- {model_name}: AUROC={row['AUROC']:.3f}, balanced_accuracy={row['balanced_accuracy']:.3f}, "
                f"ECE={row['ECE']:.3f}, n={int(row['n_samples'])} ({row['claim_status']})."
            )
    fixed = threshold_sensitivity[
        (threshold_sensitivity["cohort"] == "GSE145996+PHS000452_LIU_LIKE_PRE")
        & (threshold_sensitivity["threshold_policy"] == "fixed_probability_0_5")
    ].set_index("model_name")
    if not fixed.empty:
        lines.extend(["", "## Fixed 0.5 Threshold Sensitivity", ""])
        for model_name in [TARGET_LOCKED, TRANSFER_HEAD]:
            if model_name in fixed.index:
                row = fixed.loc[model_name]
                lines.append(
                    f"- {model_name}: balanced_accuracy={row['balanced_accuracy']:.3f}, "
                    f"sensitivity={row['sensitivity']:.3f}, specificity={row['specificity']:.3f}."
                )
    lines.extend(["", "## Discovery LODO Selection", ""])
    for _, row in inner.iterrows():
        lines.append(
            f"- {row['model_name']}: inner mean AUROC={row['inner_mean_AUROC']:.3f}, "
            f"inner min AUROC={row['inner_min_AUROC']:.3f}, status={row['claim_status']}."
        )
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "Allowed: report the transfer head as a transparent secondary rescue that improves the weak PD1-like stress cohorts.",
            "Not allowed: replace the locked external-validation claim with the transfer-head result unless a new independent external cohort validates it after this model is frozen.",
            "",
        ]
    )
    (out_dir / "PD1_LIKE_EXTERNAL_RESCUE_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote PD1-like rescue outputs to {out_dir}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--out", default="results/pd1_like_external_rescue")
    args = parser.parse_args()
    run(ROOT / args.processed_dir, ROOT / args.out)


if __name__ == "__main__":
    main()
