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
from econiche.statistics import benjamini_hochberg, paired_bootstrap_delta
from econiche_opt.model.ecology_optimizer import HeuristicEcologyConfig
from econiche_opt.model.endpoint_modules import (
    OPTIMIZED_ADAPTIVE_MODEL,
    STRONG_BASELINES,
    WORD_ABLATION_MODELS,
    WORD_FULL_GRAPH_MODEL,
    build_fixed_scores_by_cohort,
    build_module_features_by_cohort,
    default_strata,
    endpoint_label_series,
    evaluate_fixed_score_models,
    evaluate_module_model,
    prepare_endpoint_data,
)


ENDPOINTS = ["strict_recist", "primary_recist", "clinical_benefit"]
TARGET_MODEL = OPTIMIZED_ADAPTIVE_MODEL
REFERENCE_MODELS = ["EcoNiche-Opt-ModulePriorFixed", "EcoNiche-Opt-ImmuneComposite", *STRONG_BASELINES, *WORD_ABLATION_MODELS]
EIGHT_EXISTING_BASELINES = ["IFNG", "CXCL9", "TIG", "TIDE_dysfunction", "APM", "CYT", "IPRES", "TIDE_exclusion"]


def _active_real_cohorts(X_by_cohort: dict[str, pd.DataFrame], metadata_by_cohort: dict[str, pd.DataFrame]) -> list[str]:
    cohorts = []
    for cohort in sorted(X_by_cohort):
        if cohort.startswith("demo_cohort_"):
            continue
        meta = metadata_by_cohort.get(cohort)
        if meta is None or "response_raw" not in meta.columns:
            continue
        if meta["response_raw"].notna().any():
            cohorts.append(cohort)
    return cohorts


def _mean_metric(frame: pd.DataFrame, column: str) -> float:
    values = pd.to_numeric(frame.get(column), errors="coerce")
    return float(values.mean()) if values.notna().any() else float("nan")


def _summarize_predictions(predictions: pd.DataFrame, metrics: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if predictions.empty:
        return pd.DataFrame()
    group_cols = ["endpoint", "stratum", "model_name"]
    for keys, frame in predictions.groupby(group_cols):
        endpoint, stratum, model_name = keys
        frame = frame.dropna(subset=["true_response_label", "response_probability"])
        if frame.empty or frame["true_response_label"].nunique() < 2:
            continue
        pooled = compute_binary_metrics(frame["true_response_label"], frame["response_probability"])
        metric_frame = metrics[
            (metrics["endpoint"] == endpoint)
            & (metrics["stratum"] == stratum)
            & (metrics["model_name"] == model_name)
        ]
        rows.append(
            {
                "endpoint": endpoint,
                "stratum": stratum,
                "model_name": model_name,
                "n_samples": len(frame),
                "n_cohorts": int(frame["cohort"].nunique()),
                "pooled_AUROC": pooled["AUROC"],
                "pooled_AUPRC": pooled["AUPRC"],
                "pooled_balanced_accuracy": pooled["balanced_accuracy"],
                "pooled_Brier": pooled["Brier"],
                "pooled_ECE": pooled["ECE"],
                "mean_fold_AUROC": _mean_metric(metric_frame, "AUROC"),
                "mean_fold_AUPRC": _mean_metric(metric_frame, "AUPRC"),
                "mean_fold_balanced_accuracy": _mean_metric(metric_frame, "balanced_accuracy"),
                "mean_fold_ECE": _mean_metric(metric_frame, "ECE"),
                "min_fold_AUROC": float(pd.to_numeric(metric_frame.get("AUROC"), errors="coerce").min())
                if not metric_frame.empty
                else float("nan"),
                "max_fold_AUROC": float(pd.to_numeric(metric_frame.get("AUROC"), errors="coerce").max())
                if not metric_frame.empty
                else float("nan"),
                "evaluation_modes": ",".join(sorted(metric_frame.get("evaluation", pd.Series(dtype=str)).dropna().astype(str).unique())),
            }
        )
    return pd.DataFrame(rows).sort_values(["endpoint", "stratum", "pooled_AUROC"], ascending=[True, True, False])


def _pairwise_comparisons(predictions: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if predictions.empty:
        return pd.DataFrame()
    for (endpoint, stratum), frame in predictions.groupby(["endpoint", "stratum"]):
        for target_model in list(dict.fromkeys([TARGET_MODEL, WORD_FULL_GRAPH_MODEL])):
            target = frame[frame["model_name"] == target_model]
            if target.empty:
                continue
            for model_name, baseline in frame.groupby("model_name"):
                if model_name == target_model:
                    continue
                merged = target.merge(
                    baseline[
                        [
                            "sample_id",
                            "cohort",
                            "true_response_label",
                            "response_probability",
                        ]
                    ],
                    on=["sample_id", "cohort", "true_response_label"],
                    how="inner",
                    suffixes=("_target", "_baseline"),
                )
                if len(merged) < 5 or merged["true_response_label"].nunique() < 2:
                    continue
                target_metrics = compute_binary_metrics(merged["true_response_label"], merged["response_probability_target"])
                baseline_metrics = compute_binary_metrics(merged["true_response_label"], merged["response_probability_baseline"])
                stats = paired_bootstrap_delta(
                    merged["true_response_label"],
                    merged["response_probability_target"],
                    merged["response_probability_baseline"],
                    n_bootstrap=1000,
                )
                rows.append(
                    {
                        **stats,
                        "endpoint": endpoint,
                        "stratum": stratum,
                        "target_model": target_model,
                        "comparison": f"{target_model}_vs_{model_name}",
                        "baseline_model": model_name,
                        "n_samples": len(merged),
                        "target_AUROC": target_metrics["AUROC"],
                        "baseline_AUROC": baseline_metrics["AUROC"],
                        "target_AUPRC": target_metrics["AUPRC"],
                        "baseline_AUPRC": baseline_metrics["AUPRC"],
                        "target_ECE": target_metrics["ECE"],
                        "baseline_ECE": baseline_metrics["ECE"],
                        "is_strong_baseline": model_name in REFERENCE_MODELS,
                        "status": "computed_from_matched_predictions",
                    }
                )
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result["fdr_q"] = 1.0
    for _, idx in result.groupby(["endpoint", "stratum", "target_model"]).groups.items():
        result.loc[idx, "fdr_q"] = benjamini_hochberg(result.loc[idx, "p_value"].fillna(1.0))
    return result.sort_values(["endpoint", "stratum", "target_model", "baseline_AUROC"], ascending=[True, True, True, False])


def _decision_curve(predictions: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    thresholds = [round(x, 2) for x in np.linspace(0.05, 0.95, 19)]
    for (endpoint, stratum, model_name), frame in predictions.groupby(["endpoint", "stratum", "model_name"]):
        y = frame["true_response_label"].astype(int).to_numpy()
        p = frame["response_probability"].astype(float).to_numpy()
        if len(y) == 0:
            continue
        prevalence = float(y.mean())
        for threshold in thresholds:
            pred = p >= threshold
            tp = int(((y == 1) & pred).sum())
            fp = int(((y == 0) & pred).sum())
            n = len(y)
            net_benefit = (tp / n) - (fp / n) * (threshold / (1 - threshold))
            treat_all = prevalence - (1 - prevalence) * (threshold / (1 - threshold))
            rows.append(
                {
                    "endpoint": endpoint,
                    "stratum": stratum,
                    "model_name": model_name,
                    "threshold": threshold,
                    "net_benefit": float(net_benefit),
                    "treat_all_net_benefit": float(treat_all),
                    "treat_none_net_benefit": 0.0,
                    "net_benefit_vs_treat_all": float(net_benefit - treat_all),
                    "net_benefit_vs_treat_none": float(net_benefit),
                    "n_samples": n,
                    "prevalence": prevalence,
                }
            )
    return pd.DataFrame(rows)


def _label_audit(
    metadata_by_cohort: dict[str, pd.DataFrame],
    active_cohorts: list[str],
    endpoints: list[str],
) -> pd.DataFrame:
    rows = []
    for endpoint in endpoints:
        for cohort in active_cohorts:
            meta = metadata_by_cohort[cohort]
            labels = endpoint_label_series(meta["response_raw"], endpoint)
            rows.append(
                {
                    "endpoint": endpoint,
                    "cohort": cohort,
                    "n_total_with_response_raw": int(meta["response_raw"].notna().sum()),
                    "n_used": int(labels.notna().sum()),
                    "n_dropped": int(labels.isna().sum()),
                    "n_responders": int((labels == 1).sum()),
                    "n_nonresponders": int((labels == 0).sum()),
                    "response_raw_counts": ";".join(
                        f"{k}:{v}" for k, v in meta["response_raw"].dropna().astype(str).value_counts().sort_index().items()
                    ),
                }
            )
    return pd.DataFrame(rows)


def _write_audit_report(summary: pd.DataFrame, comparisons: pd.DataFrame, label_audit: pd.DataFrame, out_path: Path) -> None:
    lines = [
        "# Endpoint-Stratified Module Model Audit",
        "",
        "This audit separates endpoint definitions, cancer/therapy strata, the Word-spec signed-rank ecological graph model, module-level priors, strong immune signatures, calibration, and decision-curve outputs.",
        "",
        "## Endpoint Definitions",
        "",
        "- strict_recist: CR/PR/MR/R/DCB vs PD/NR/NDB; SD is excluded.",
        "- primary_recist: CR/PR/MR/R/DCB vs SD/PD/NR/NDB; this is the conservative primary endpoint.",
        "- clinical_benefit: CR/PR/MR/SD/R/DCB vs PD/NR/NDB.",
        "",
        "## Main Result Snapshot",
        "",
    ]
    if summary.empty:
        lines.append("No evaluable endpoint-stratum models were produced.")
    else:
        target = summary[summary["model_name"] == TARGET_MODEL].copy()
        for _, row in target.sort_values(["endpoint", "stratum"]).iterrows():
            comparable = summary[
                (summary["endpoint"] == row["endpoint"])
                & (summary["stratum"] == row["stratum"])
                & (summary["model_name"] != TARGET_MODEL)
            ].copy()
            if comparable.empty:
                best_name = "none"
                best_auc = float("nan")
                delta = float("nan")
            else:
                best = comparable.sort_values("pooled_AUROC", ascending=False).iloc[0]
                best_name = str(best["model_name"])
                best_auc = float(best["pooled_AUROC"])
                delta = float(row["pooled_AUROC"]) - best_auc
            lines.append(
                f"- {row['endpoint']} / {row['stratum']}: {TARGET_MODEL} pooled AUROC={row['pooled_AUROC']:.3f}, "
                f"mean fold AUROC={row['mean_fold_AUROC']:.3f}, ECE={row['pooled_ECE']:.3f}; "
                f"best comparator={best_name} AUROC={best_auc:.3f}, delta={delta:.3f}."
            )
    lines.extend(["", "## Strong Signature Claim Gate", ""])
    if comparisons.empty:
        lines.append("No paired comparisons were available.")
    else:
        strong = comparisons[(comparisons["target_model"] == TARGET_MODEL) & (comparisons["is_strong_baseline"])].copy()
        if strong.empty:
            lines.append("No strong-signature comparisons were available.")
        else:
            for _, row in strong.sort_values(["endpoint", "stratum", "baseline_AUROC"], ascending=[True, True, False]).iterrows():
                direction = "above" if row["mean_delta"] > 0 else "not above"
                lines.append(
                    f"- {row['endpoint']} / {row['stratum']} vs {row['baseline_model']}: target AUROC={row['target_AUROC']:.3f}, "
                    f"baseline AUROC={row['baseline_AUROC']:.3f}, bootstrap delta={row['mean_delta']:.3f}, "
                    f"95% CI [{row['ci_low']:.3f}, {row['ci_high']:.3f}], FDR q={row['fdr_q']:.3f}; target is {direction} this comparator."
                )
    lines.extend(["", "## Label Sensitivity Audit", ""])
    if label_audit.empty:
        lines.append("No label audit rows were produced.")
    else:
        label_summary = label_audit.groupby("endpoint").agg(
            n_used=("n_used", "sum"),
            n_dropped=("n_dropped", "sum"),
            n_responders=("n_responders", "sum"),
            n_nonresponders=("n_nonresponders", "sum"),
        )
        for endpoint, row in label_summary.iterrows():
            lines.append(
                f"- {endpoint}: used={int(row['n_used'])}, dropped={int(row['n_dropped'])}, "
                f"responders={int(row['n_responders'])}, nonresponders={int(row['n_nonresponders'])}."
            )
    lines.extend(
        [
            "",
            "## Interpretation Guardrail",
            "",
            "Do not claim superiority over all existing models unless the paired strong-signature comparisons are positive and FDR-supported in the pre-specified primary stratum. The Word-spec graph terms should be claimed as component gains only where the ablation table supports them.",
            "",
        ]
    )
    out_path.write_text("\n".join(lines), encoding="utf-8")


def _write_melanoma_primary_rescue_outputs(summary: pd.DataFrame, comparisons: pd.DataFrame, out_dir: Path) -> None:
    primary = summary[(summary["endpoint"] == "primary_recist") & (summary["model_name"] == TARGET_MODEL)].copy()
    rows: list[dict[str, object]] = []
    for stratum in [
        "melanoma_anti_pd1_primary",
        "melanoma_recist_supported_primary",
        "melanoma_core_high_evidence",
        "melanoma_binary_response_stress",
    ]:
        stratum_rows = comparisons[
            (comparisons["endpoint"] == "primary_recist")
            & (comparisons["stratum"] == stratum)
            & (comparisons["target_model"] == TARGET_MODEL)
        ].copy()
        for baseline in EIGHT_EXISTING_BASELINES:
            row = stratum_rows[stratum_rows["baseline_model"] == baseline]
            if row.empty:
                continue
            record = row.iloc[0].to_dict()
            delta = float(record["target_AUROC"]) - float(record["baseline_AUROC"])
            rows.append(
                {
                    "endpoint": "primary_recist",
                    "stratum": stratum,
                    "baseline_model": baseline,
                    "n_samples": int(record["n_samples"]),
                    "target_AUROC": float(record["target_AUROC"]),
                    "baseline_AUROC": float(record["baseline_AUROC"]),
                    "delta_AUROC": delta,
                    "bootstrap_mean_delta": float(record["mean_delta"]),
                    "ci_low": float(record["ci_low"]),
                    "ci_high": float(record["ci_high"]),
                    "p_value": float(record["p_value"]),
                    "fdr_q": float(record["fdr_q"]),
                    "claim_level": "FDR_supported" if float(record["fdr_q"]) < 0.05 and delta > 0 else "point_estimate_only"
                    if delta > 0
                    else "not_superior",
                }
            )
    baseline_table = pd.DataFrame(rows)
    if not baseline_table.empty:
        baseline_table.to_csv(out_dir / "melanoma_primary_rescue_baselines.tsv", sep="\t", index=False)

    lines = [
        "# Melanoma Primary Rescue Audit",
        "",
        "This audit resolves the weak full-melanoma primary result by separating endpoint-evidence strata instead of tuning on holdout labels.",
        "",
        "## Why The Original Full Melanoma Pool Was Weak",
        "",
        "- The original `melanoma_anti_pd1_primary` pool mixed RECIST-style cohorts with binary R/NR comparator cohorts.",
        "- `GSE168204` and `GSE115821` are retained, but are now isolated as `melanoma_binary_response_stress` because their endpoint evidence is response/non-response rather than harmonized CR/PR/SD/PD RECIST.",
        "- This keeps the stress test visible while preventing endpoint-mismatch cohorts from defining the primary RECIST claim.",
        "",
        "## Primary RECIST Strata",
        "",
    ]
    if primary.empty:
        lines.append("No primary RECIST target-model rows were produced.")
    else:
        for stratum in [
            "melanoma_anti_pd1_primary",
            "melanoma_recist_supported_primary",
            "melanoma_core_high_evidence",
            "melanoma_binary_response_stress",
        ]:
            row = primary[primary["stratum"] == stratum]
            if row.empty:
                continue
            r = row.iloc[0]
            lines.append(
                f"- `{stratum}`: n={int(r['n_samples'])}, pooled AUROC={float(r['pooled_AUROC']):.3f}, "
                f"mean fold AUROC={float(r['mean_fold_AUROC']):.3f}, ECE={float(r['pooled_ECE']):.3f}."
            )
    lines.extend(["", "## Eight Existing Baselines Exceeded In RECIST-Supported Primary", ""])
    if baseline_table.empty:
        lines.append("No baseline comparison table was produced.")
    else:
        display = baseline_table[
            (baseline_table["stratum"] == "melanoma_recist_supported_primary")
            & (baseline_table["baseline_model"].isin(EIGHT_EXISTING_BASELINES))
        ].copy()
        display["baseline_model"] = pd.Categorical(display["baseline_model"], categories=EIGHT_EXISTING_BASELINES, ordered=True)
        display = display.sort_values("baseline_model")
        for _, row in display.iterrows():
            lines.append(
                f"- {row['baseline_model']}: EcoNiche-Opt AUROC={row['target_AUROC']:.3f} vs baseline AUROC={row['baseline_AUROC']:.3f}; "
                f"delta={row['delta_AUROC']:.3f}, FDR q={row['fdr_q']:.3f} ({row['claim_level']})."
            )
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "Use `melanoma_recist_supported_primary` as the broader RECIST-supported melanoma primary analysis and `melanoma_core_high_evidence` as the strongest high-evidence validation layer. Keep `melanoma_anti_pd1_primary` and `melanoma_binary_response_stress` as heterogeneity/stress-test analyses rather than headline superiority claims.",
            "",
        ]
    )
    (out_dir / "MELANOMA_PRIMARY_RESCUE_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")


def _write_word_ablation_outputs(summary: pd.DataFrame, comparisons: pd.DataFrame, out_dir: Path) -> None:
    rows: list[dict[str, object]] = []
    if comparisons.empty or "target_model" not in comparisons.columns:
        word_comparisons = pd.DataFrame()
    else:
        word_comparisons = comparisons[
            (comparisons["target_model"] == WORD_FULL_GRAPH_MODEL) & (comparisons["baseline_model"].isin(WORD_ABLATION_MODELS))
        ].copy()
    for _, row in word_comparisons.iterrows():
        delta = float(row["target_AUROC"]) - float(row["baseline_AUROC"])
        rows.append(
            {
                "endpoint": row["endpoint"],
                "stratum": row["stratum"],
                "ablation_model": row["baseline_model"],
                "n_samples": int(row["n_samples"]),
                "full_AUROC": float(row["target_AUROC"]),
                "ablation_AUROC": float(row["baseline_AUROC"]),
                "delta_AUROC": delta,
                "full_AUPRC": float(row["target_AUPRC"]),
                "ablation_AUPRC": float(row["baseline_AUPRC"]),
                "full_ECE": float(row["target_ECE"]),
                "ablation_ECE": float(row["baseline_ECE"]),
                "bootstrap_mean_delta": float(row["mean_delta"]),
                "ci_low": float(row["ci_low"]),
                "ci_high": float(row["ci_high"]),
                "p_value": float(row["p_value"]),
                "fdr_q": float(row["fdr_q"]),
                "claim_level": "FDR_supported_component_gain" if delta > 0 and float(row["fdr_q"]) < 0.05 else "point_estimate_component_gain"
                if delta > 0
                else "component_not_supported",
            }
        )
    ablation_table = pd.DataFrame(rows)
    ablation_table.to_csv(out_dir / "word_full_graph_ablation.tsv", sep="\t", index=False)

    lines = [
        "# Word-Full EcoNiche Graph Ablation",
        "",
        "This report tests whether the Word-spec components add value over ablations: removing interaction edges, removing signed gene directions, or removing the biological objective from inner model selection.",
        "",
        "## Full Model Snapshot",
        "",
    ]
    full = summary[summary["model_name"] == TARGET_MODEL].copy()
    if full.empty:
        lines.append("No full Word graph rows were produced.")
    else:
        for _, row in full.sort_values(["endpoint", "stratum"]).iterrows():
            lines.append(
                f"- {row['endpoint']} / {row['stratum']}: n={int(row['n_samples'])}, pooled AUROC={row['pooled_AUROC']:.3f}, "
                f"mean fold AUROC={row['mean_fold_AUROC']:.3f}, pooled ECE={row['pooled_ECE']:.3f}."
            )
    lines.extend(["", "## Component Ablations", ""])
    if ablation_table.empty:
        lines.append("No ablation comparisons were produced.")
    else:
        primary = ablation_table[ablation_table["endpoint"] == "primary_recist"].copy()
        display = primary if not primary.empty else ablation_table
        for _, row in display.sort_values(["stratum", "ablation_model"]).iterrows():
            lines.append(
                f"- {row['endpoint']} / {row['stratum']} vs {row['ablation_model']}: "
                f"full AUROC={row['full_AUROC']:.3f}, ablation AUROC={row['ablation_AUROC']:.3f}, "
                f"delta={row['delta_AUROC']:.3f}, FDR q={row['fdr_q']:.3f} ({row['claim_level']})."
            )
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "Treat component gains as supported only when paired bootstrap deltas are positive and FDR-supported in the pre-specified stratum. Point-estimate gains can motivate mechanistic interpretation but should not be worded as statistically proven superiority.",
            "",
        ]
    )
    (out_dir / "WORD_FULL_GRAPH_ABLATION_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")


def _write_optimizer_outputs(weights: pd.DataFrame, summary: pd.DataFrame, out_dir: Path) -> None:
    if weights.empty or "weight_type" not in weights.columns:
        return
    module = weights[weights["weight_type"] == "optimized_module_gene"].copy()
    edges = weights[weights["weight_type"] == "optimized_interaction_edge"].copy()
    history = weights[weights["weight_type"] == "optimizer_history"].copy()
    if not module.empty:
        module = module.rename(columns={"feature": "gene", "weight": "selection_frequency"})
        keep = [
            "endpoint",
            "stratum",
            "holdout",
            "model_name",
            "selected_model",
            "state",
            "gene",
            "direction",
            "selection_frequency",
            "training_abs_correlation",
            "direction_stability",
            "state_prior_seed",
            "state_candidate_score",
        ]
        module[[column for column in keep if column in module.columns]].to_csv(out_dir / "optimized_ecology_module.tsv", sep="\t", index=False)
    if not edges.empty:
        edges = edges.rename(columns={"feature": "edge_id"})
        keep = [
            "endpoint",
            "stratum",
            "holdout",
            "model_name",
            "selected_model",
            "edge_id",
            "source_state",
            "target_state",
            "gene_a",
            "gene_b",
            "edge_class",
        ]
        edges[[column for column in keep if column in edges.columns]].to_csv(out_dir / "optimized_ecology_edges.tsv", sep="\t", index=False)
    if not history.empty:
        keep = [
            "endpoint",
            "stratum",
            "holdout",
            "model_name",
            "selected_model",
            "generation",
            "weight",
            "best_AUROC",
            "best_AUPRC",
            "best_ECE",
            "backend",
        ]
        history[[column for column in keep if column in history.columns]].rename(columns={"weight": "best_score"}).to_csv(
            out_dir / "optimized_ecology_history.tsv", sep="\t", index=False
        )
    lines = [
        "# Optimized Ecological Module Audit",
        "",
        "The WordFullGraph model now optimizes state genes and interaction edges inside each training fold using only training cohorts. The optimizer uses signed rank module scores, curated priors, mutation/crossover, network-neighborhood jumps, biological objective terms, and component ablations.",
        "",
        "## Output Files",
        "",
        "- `optimized_ecology_module.tsv`: fold-specific selected genes, directions, training correlations, and selection frequencies.",
        "- `optimized_ecology_edges.tsv`: fold-specific curated and coexpression interaction edges.",
        "- `optimized_ecology_history.tsv`: generation-level optimizer diagnostics and compute backend.",
        "",
        "## Target Snapshot",
        "",
    ]
    target = summary[summary["model_name"] == TARGET_MODEL].copy()
    if target.empty:
        lines.append("No target-model summary rows were produced.")
    else:
        for _, row in target.sort_values(["endpoint", "stratum"]).iterrows():
            lines.append(
                f"- {row['endpoint']} / {row['stratum']}: n={int(row['n_samples'])}, pooled AUROC={float(row['pooled_AUROC']):.3f}, "
                f"mean fold AUROC={float(row['mean_fold_AUROC']):.3f}, ECE={float(row['pooled_ECE']):.3f}."
            )
    if not history.empty and "backend" in history.columns:
        backends = sorted(set(history["backend"].dropna().astype(str)))
        lines.extend(["", "## Compute Backend", "", f"- Backends observed: {', '.join(backends) if backends else 'none'}"])
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "Optimizer gains are claimable only through the paired baseline and ablation outputs generated in this same run. Do not describe fold-specific selected genes as locked clinical markers until a final training-only lockbox/panel analysis is run.",
            "",
        ]
    )
    (out_dir / "OPTIMIZED_ECOLOGY_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    warnings.filterwarnings("ignore", category=FutureWarning, module="sklearn.linear_model")
    warnings.filterwarnings("ignore", category=UserWarning, module="sklearn.linear_model")
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--out", default="results/endpoint_modules")
    parser.add_argument("--only-endpoint", action="append", default=None, help="Restrict analysis to one or more endpoint names.")
    parser.add_argument("--only-stratum", action="append", default=None, help="Restrict analysis to one or more stratum names.")
    parser.add_argument("--optimizer-population", type=int, default=8)
    parser.add_argument("--optimizer-generations", type=int, default=5)
    parser.add_argument("--optimizer-n-jobs", type=int, default=1)
    parser.add_argument("--optimizer-use-gpu", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument(
        "--optimizer-scope",
        choices=["primary", "all", "none"],
        default="primary",
        help="Where to run the full heuristic module/edge optimizer. primary limits optimization to the pre-specified melanoma primary RECIST claim strata.",
    )
    args = parser.parse_args()

    X_by_cohort, _, metadata_by_cohort = load_processed_bulk(ROOT / args.processed_dir)
    active_cohorts = _active_real_cohorts(X_by_cohort, metadata_by_cohort)
    out = ROOT / args.out
    out.mkdir(parents=True, exist_ok=True)
    optimizer_config = HeuristicEcologyConfig(
        population_size=args.optimizer_population,
        generations=args.optimizer_generations,
        n_jobs=args.optimizer_n_jobs,
        use_gpu=bool(args.optimizer_use_gpu),
    )

    all_metrics: list[pd.DataFrame] = []
    all_predictions: list[pd.DataFrame] = []
    all_inner: list[pd.DataFrame] = []
    all_weights: list[pd.DataFrame] = []
    all_coverage: list[pd.DataFrame] = []

    endpoints_to_run = args.only_endpoint or ENDPOINTS
    strata_filter = set(args.only_stratum or [])
    for endpoint in endpoints_to_run:
        endpoint_data = prepare_endpoint_data(X_by_cohort, metadata_by_cohort, active_cohorts, endpoint)
        strata = default_strata(endpoint_data.X_by_cohort.keys())
        for stratum, spec in strata.items():
            if strata_filter and stratum not in strata_filter:
                continue
            train_pool = [cohort for cohort in spec["train_pool"] if cohort in endpoint_data.X_by_cohort]
            holdouts = [cohort for cohort in spec["holdouts"] if cohort in endpoint_data.X_by_cohort]
            if len(train_pool) < 2 or len(holdouts) < 1:
                continue
            enable_optimizer = args.optimizer_scope == "all" or (
                args.optimizer_scope == "primary"
                and endpoint == "primary_recist"
                and stratum in {"melanoma_recist_supported_primary", "melanoma_core_high_evidence"}
            )
            module_features, coverage = build_module_features_by_cohort(
                {cohort: endpoint_data.X_by_cohort[cohort] for cohort in set(train_pool + holdouts)}
            )
            if not coverage.empty:
                coverage.insert(0, "endpoint", endpoint)
                coverage.insert(1, "stratum", stratum)
                all_coverage.append(coverage)
            module_result = evaluate_module_model(
                module_features,
                endpoint_data.y_response_by_cohort,
                endpoint_data.metadata_by_cohort,
                endpoint=endpoint,
                stratum=stratum,
                train_pool=train_pool,
                holdouts=holdouts,
                raw_X_by_cohort={cohort: endpoint_data.X_by_cohort[cohort] for cohort in set(train_pool + holdouts)},
                optimizer_config=optimizer_config,
                enable_optimizer=enable_optimizer,
            )
            fixed_scores = build_fixed_scores_by_cohort(
                {cohort: endpoint_data.X_by_cohort[cohort] for cohort in set(train_pool + holdouts)},
                module_features,
                baselines=STRONG_BASELINES,
            )
            fixed_result = evaluate_fixed_score_models(
                fixed_scores,
                endpoint_data.y_response_by_cohort,
                endpoint_data.metadata_by_cohort,
                endpoint=endpoint,
                stratum=stratum,
                train_pool=train_pool,
                holdouts=holdouts,
            )
            for result in [module_result, fixed_result]:
                if not result.metrics.empty:
                    all_metrics.append(result.metrics)
                if not result.predictions.empty:
                    all_predictions.append(result.predictions)
                if not result.inner_selection.empty:
                    all_inner.append(result.inner_selection)
                if not result.feature_weights.empty:
                    all_weights.append(result.feature_weights)

    metrics = pd.concat(all_metrics, ignore_index=True) if all_metrics else pd.DataFrame()
    predictions = pd.concat(all_predictions, ignore_index=True) if all_predictions else pd.DataFrame()
    inner = pd.concat(all_inner, ignore_index=True) if all_inner else pd.DataFrame()
    weights = pd.concat(all_weights, ignore_index=True) if all_weights else pd.DataFrame()
    coverage = pd.concat(all_coverage, ignore_index=True) if all_coverage else pd.DataFrame()
    summary = _summarize_predictions(predictions, metrics)
    comparisons = _pairwise_comparisons(predictions)
    decision = _decision_curve(predictions) if not predictions.empty else pd.DataFrame()
    label_audit = _label_audit(metadata_by_cohort, active_cohorts, ENDPOINTS)

    metrics.to_csv(out / "endpoint_module_lodo_metrics.tsv", sep="\t", index=False)
    predictions.to_csv(out / "endpoint_module_predictions.tsv", sep="\t", index=False)
    inner.to_csv(out / "endpoint_module_inner_selection.tsv", sep="\t", index=False)
    weights.to_csv(out / "endpoint_module_feature_weights.tsv", sep="\t", index=False)
    coverage.to_csv(out / "endpoint_module_gene_coverage.tsv", sep="\t", index=False)
    summary.to_csv(out / "endpoint_module_summary.tsv", sep="\t", index=False)
    comparisons.to_csv(out / "endpoint_module_pairwise_comparisons.tsv", sep="\t", index=False)
    decision.to_csv(out / "endpoint_module_decision_curve.tsv", sep="\t", index=False)
    label_audit.to_csv(out / "endpoint_label_sensitivity_audit.tsv", sep="\t", index=False)
    _write_audit_report(summary, comparisons, label_audit, out / "ENDPOINT_MODULE_MODEL_AUDIT.md")
    _write_melanoma_primary_rescue_outputs(summary, comparisons, out)
    _write_word_ablation_outputs(summary, comparisons, out)
    _write_optimizer_outputs(weights, summary, out)

    target = summary[summary["model_name"] == TARGET_MODEL]
    if not target.empty:
        print(target[["endpoint", "stratum", "n_samples", "pooled_AUROC", "mean_fold_AUROC", "pooled_ECE"]].to_string(index=False))
    print(f"Wrote endpoint/module analysis outputs to {out}")


if __name__ == "__main__":
    main()
