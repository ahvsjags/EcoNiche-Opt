from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import average_precision_score, balanced_accuracy_score, roc_auc_score


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.statistics import benjamini_hochberg  # noqa: E402
from scripts.analysis.run_gpu_bioprior_rescue_combo_search import (  # noqa: E402
    BASE_GENES,
    LIPID_PI3K_PRIOR_GENES,
    PRIMARY_COHORTS,
    STRICT_EXTERNAL_COHORTS,
    build_transforms,
    candidate_specs,
    external_score,
    labels_for_endpoint,
    load_bulk,
    primary_lodo,
)


def _selected_candidate(selection_path: Path) -> str:
    if not selection_path.exists():
        raise FileNotFoundError(selection_path)
    frame = pd.read_csv(selection_path, sep="\t")
    rows = frame[frame["selection_id"].astype(str).str.contains("gpu_", na=False)]
    if rows.empty:
        raise ValueError("Selection table lacks gpu_* row")
    row = rows.sort_values("strict_external_AUROC", ascending=False).iloc[0]
    if str(row["selection_boundary"]) != "candidate_and_weight_selected_by_primary_lodo_only_with_biological_prior":
        raise ValueError("Selected GPU rescue combo is not primary-LODO-only")
    return str(row["candidate"])


def _paired_metric_rows(
    predictions: pd.DataFrame,
    summary: pd.DataFrame,
    context: str,
    target_candidate: str,
    ablated_candidate: str,
    n_bootstrap: int,
) -> list[dict[str, object]]:
    pred = predictions[predictions["candidate"].astype(str).isin([target_candidate, ablated_candidate])].copy()
    pred["key"] = pred["cohort"].astype(str) + "::" + pred["sample_id"].astype(str)
    y = pred.drop_duplicates("key").set_index("key")["true_response_label"].astype(int)
    wide = pred.pivot_table(index="key", columns="candidate", values="response_probability", aggfunc="first").dropna()
    y = y.reindex(wide.index).astype(int)
    thresholds = (
        pred.drop_duplicates(["candidate", "key"])
        .groupby("candidate")["threshold"]
        .median()
        .to_dict()
    )
    target = wide[target_candidate].to_numpy(dtype=float)
    ablated = wide[ablated_candidate].to_numpy(dtype=float)
    labels = y.to_numpy(dtype=int)

    rng = np.random.default_rng(20260527)
    auc_deltas: list[float] = []
    auprc_deltas: list[float] = []
    for _ in range(n_bootstrap):
        idx = rng.integers(0, len(labels), len(labels))
        yy = labels[idx]
        if len(np.unique(yy)) < 2:
            continue
        auc_deltas.append(float(roc_auc_score(yy, target[idx]) - roc_auc_score(yy, ablated[idx])))
        auprc_deltas.append(float(average_precision_score(yy, target[idx]) - average_precision_score(yy, ablated[idx])))
    auc_arr = np.asarray(auc_deltas, dtype=float)
    auprc_arr = np.asarray(auprc_deltas, dtype=float)

    target_row = summary[summary["candidate"].astype(str).eq(target_candidate)].iloc[0]
    ablated_row = summary[summary["candidate"].astype(str).eq(ablated_candidate)].iloc[0]
    target_pred = target >= float(thresholds[target_candidate])
    ablated_pred = ablated >= float(thresholds[ablated_candidate])
    rows = [
        {
            "context": context,
            "comparison_id": "selected_gpu_combo_vs_base_rescue",
            "target_candidate": target_candidate,
            "ablated_candidate": ablated_candidate,
            "metric": "AUROC",
            "target_value": float(target_row["AUROC"]),
            "ablated_value": float(ablated_row["AUROC"]),
            "delta": float(target_row["AUROC"] - ablated_row["AUROC"]),
            "bootstrap_mean_delta": float(auc_arr.mean()),
            "ci_low": float(np.quantile(auc_arr, 0.025)),
            "ci_high": float(np.quantile(auc_arr, 0.975)),
            "two_sided_p": float(min(1.0, 2.0 * min((auc_arr <= 0).mean(), (auc_arr >= 0).mean()))),
            "n_samples": int(len(labels)),
            "claim_boundary": "component_ablation_evaluation_not_used_for_candidate_selection",
        },
        {
            "context": context,
            "comparison_id": "selected_gpu_combo_vs_base_rescue",
            "target_candidate": target_candidate,
            "ablated_candidate": ablated_candidate,
            "metric": "AUPRC",
            "target_value": float(target_row["AUPRC"]),
            "ablated_value": float(ablated_row["AUPRC"]),
            "delta": float(target_row["AUPRC"] - ablated_row["AUPRC"]),
            "bootstrap_mean_delta": float(auprc_arr.mean()),
            "ci_low": float(np.quantile(auprc_arr, 0.025)),
            "ci_high": float(np.quantile(auprc_arr, 0.975)),
            "two_sided_p": float(min(1.0, 2.0 * min((auprc_arr <= 0).mean(), (auprc_arr >= 0).mean()))),
            "n_samples": int(len(labels)),
            "claim_boundary": "component_ablation_evaluation_not_used_for_candidate_selection",
        },
        {
            "context": context,
            "comparison_id": "selected_gpu_combo_vs_base_rescue",
            "target_candidate": target_candidate,
            "ablated_candidate": ablated_candidate,
            "metric": "balanced_accuracy",
            "target_value": float(balanced_accuracy_score(labels, target_pred.astype(int))),
            "ablated_value": float(balanced_accuracy_score(labels, ablated_pred.astype(int))),
            "delta": float(
                balanced_accuracy_score(labels, target_pred.astype(int))
                - balanced_accuracy_score(labels, ablated_pred.astype(int))
            ),
            "bootstrap_mean_delta": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "two_sided_p": np.nan,
            "n_samples": int(len(labels)),
            "claim_boundary": "threshold_metric_reported_without_bootstrap_claim",
        },
    ]
    return rows


def run(processed_dir: Path, selection_path: Path, out_dir: Path, n_bootstrap: int) -> dict[str, str]:
    target_candidate = _selected_candidate(selection_path)
    ablated_candidate = "base_rescue_robust"
    X_by_cohort, metadata_by_cohort = load_bulk(processed_dir)
    primary_y = labels_for_endpoint(X_by_cohort, metadata_by_cohort, PRIMARY_COHORTS, "primary_recist")
    strict_y = labels_for_endpoint(
        X_by_cohort,
        metadata_by_cohort,
        [*PRIMARY_COHORTS, *STRICT_EXTERNAL_COHORTS],
        "strict_recist",
    )
    genes = sorted(set(BASE_GENES + LIPID_PI3K_PRIOR_GENES))
    transforms = build_transforms(X_by_cohort, primary_y, strict_y, genes)
    specs = [
        spec
        for spec in candidate_specs("lipid_pi3k", transform_policy="robust_only")
        if str(spec["candidate"]) in {target_candidate, ablated_candidate}
    ]
    if {str(spec["candidate"]) for spec in specs} != {target_candidate, ablated_candidate}:
        raise ValueError("Could not recover selected and ablated candidate specifications")
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    primary_summary, primary_predictions = primary_lodo(transforms, primary_y, specs, device)
    external_summary, external_predictions = external_score(transforms, primary_y, strict_y, specs, device)

    rows = []
    rows.extend(
        _paired_metric_rows(
            primary_predictions,
            primary_summary,
            "primary_melanoma_lodo",
            target_candidate,
            ablated_candidate,
            n_bootstrap,
        )
    )
    rows.extend(
        _paired_metric_rows(
            external_predictions,
            external_summary,
            "strict_melanoma_external",
            target_candidate,
            ablated_candidate,
            n_bootstrap,
        )
    )
    ablation = pd.DataFrame(rows)
    mask = ablation["two_sided_p"].notna()
    ablation.loc[mask, "two_sided_fdr_q"] = benjamini_hochberg(ablation.loc[mask, "two_sided_p"])
    ablation.loc[~mask, "two_sided_fdr_q"] = np.nan
    ablation["claim_level"] = np.where(
        (ablation["delta"] > 0) & (ablation["two_sided_fdr_q"].fillna(1.0) <= 0.05),
        "FDR_supported_component_gain",
        np.where(ablation["delta"] > 0, "point_estimate_component_gain", "component_tradeoff_or_no_gain"),
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "ablation": out_dir / "gpu_bioprior_component_ablation.tsv",
        "primary_summary": out_dir / "gpu_bioprior_component_primary_summary.tsv",
        "primary_predictions": out_dir / "gpu_bioprior_component_primary_predictions.tsv",
        "external_summary": out_dir / "gpu_bioprior_component_external_summary.tsv",
        "external_predictions": out_dir / "gpu_bioprior_component_external_predictions.tsv",
        "audit": out_dir / "GPU_BIOPRIOR_COMPONENT_ABLATION_AUDIT.md",
    }
    ablation.to_csv(outputs["ablation"], sep="\t", index=False)
    primary_summary.to_csv(outputs["primary_summary"], sep="\t", index=False)
    primary_predictions.to_csv(outputs["primary_predictions"], sep="\t", index=False)
    external_summary.to_csv(outputs["external_summary"], sep="\t", index=False)
    external_predictions.to_csv(outputs["external_predictions"], sep="\t", index=False)
    selected_external = external_summary[external_summary["candidate"].astype(str).eq(target_candidate)].iloc[0]
    base_external = external_summary[external_summary["candidate"].astype(str).eq(ablated_candidate)].iloc[0]
    lines = [
        "# GPU biological-prior component ablation",
        "",
        f"Device: {device} ({torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'cpu'}).",
        f"Target candidate: `{target_candidate}`.",
        f"Ablated candidate: `{ablated_candidate}`.",
        "",
        "This audit evaluates the frozen PLA2G2D lipid/PI3K rescue component against the base MAP4K1-TBX3/AXL rescue axis.",
        "It is an evaluation-only component ablation; strict external labels were not used for candidate or weight selection.",
        "",
        "Strict external AUROC changes from {:.3f} to {:.3f}; AUPRC changes from {:.3f} to {:.3f}.".format(
            float(base_external["AUROC"]),
            float(selected_external["AUROC"]),
            float(base_external["AUPRC"]),
            float(selected_external["AUPRC"]),
        ),
    ]
    outputs["audit"].write_text("\n".join(lines), encoding="utf-8")
    return {key: str(value) for key, value in outputs.items()}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument(
        "--selection",
        default="results/gpu_bioprior_rescue_combo_search_robust_20260527/gpu_bioprior_rescue_combo_selection.tsv",
    )
    parser.add_argument("--out", default="results/gpu_bioprior_component_ablation_20260527")
    parser.add_argument("--bootstrap", type=int, default=2000)
    args = parser.parse_args()
    outputs = run(ROOT / args.processed_dir, ROOT / args.selection, ROOT / args.out, args.bootstrap)
    print(json.dumps(outputs, ensure_ascii=False))


if __name__ == "__main__":
    main()
