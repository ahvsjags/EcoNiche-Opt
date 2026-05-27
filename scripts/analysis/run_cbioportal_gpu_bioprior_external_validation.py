from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.baselines import BASELINE_SIGNATURES, signature_score  # noqa: E402
from econiche.metrics import compute_binary_metrics  # noqa: E402
from econiche.statistics import benjamini_hochberg  # noqa: E402
from econiche_opt.model.endpoint_modules import endpoint_label_series  # noqa: E402


TARGET_GENES = ["MAP4K1", "TBX3", "AXL", "PLA2G2D"]
GROUPS = {
    "cbio_liu_dfci_only": ["CBIO_LIU_DFCI_2019_PRE"],
    "strict_cbio_liu_plus_gse145996": ["CBIO_LIU_DFCI_2019_PRE", "GSE145996"],
}
EIGHT_SIGNATURES = ["IFNG", "CXCL9", "TIG", "TIDE_dysfunction", "TIDE_exclusion", "CYT", "APM", "IPRES"]


def _minmax(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce").fillna(0.0)
    return (values - values.min()) / (values.max() - values.min() + 1e-9)


def _zscore(values: pd.DataFrame) -> pd.DataFrame:
    return ((values - values.mean()) / (values.std(ddof=0) + 1e-6)).fillna(0.0)


def _robust_zscore(values: pd.DataFrame) -> pd.DataFrame:
    median = values.median()
    mad = (values - median).abs().median() + 1e-6
    return ((values - median) / (1.4826 * mad)).clip(-5.0, 5.0).fillna(0.0)


def score_gpu_combo(X: pd.DataFrame, spec: dict[str, object]) -> pd.Series:
    missing = [gene for gene in TARGET_GENES if gene not in X.columns]
    if missing:
        raise ValueError(f"Cannot score GPU bioprior combo; missing genes: {missing}")
    locked = spec["locked_score"]
    weight_base = float(locked["weight_base"])
    weight_component = float(locked["weight_component"])
    sign = float(locked["component_training_direction_sign"])
    component_gene = str(locked["component_gene"])
    values = X[TARGET_GENES].apply(pd.to_numeric, errors="coerce")
    rz = _robust_zscore(values)
    z = _zscore(values)
    base_rz = _minmax(rz["MAP4K1"] - rz[["TBX3", "AXL"]].mean(axis=1))
    base_z = _minmax(z["MAP4K1"] - z[["TBX3", "AXL"]].mean(axis=1))
    base = _minmax(0.95 * base_rz + 0.05 * base_z)
    component = _minmax(sign * rz[component_gene])
    score = _minmax(weight_base * base + weight_component * component)
    score.name = "response_probability"
    return score


def _read_cohort(cohort: str, bulk_dir: Path, cbio_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    base = cbio_dir / cohort if cohort.startswith("CBIO_") else bulk_dir / cohort
    expr = pd.read_csv(base.with_suffix(".expr.tsv"), sep="\t", index_col=0)
    metadata = pd.read_csv(base.with_suffix(".metadata.tsv"), sep="\t")
    return expr, metadata


def _cohort_predictions(cohort: str, bulk_dir: Path, cbio_dir: Path, spec: dict[str, object]) -> tuple[pd.DataFrame, dict[str, object]]:
    expr, metadata = _read_cohort(cohort, bulk_dir, cbio_dir)
    labels = endpoint_label_series(metadata.set_index("sample_id")["response_raw"], "strict_recist").dropna().astype(int)
    common = expr.index.intersection(labels.index)
    available = [gene for gene in TARGET_GENES if gene in expr.columns]
    coverage = {
        "cohort": cohort,
        "n_target_genes": len(TARGET_GENES),
        "n_target_genes_available": len(available),
        "available_genes": ",".join(available),
        "missing_genes": ",".join(gene for gene in TARGET_GENES if gene not in available),
        "status": "ready" if len(available) == len(TARGET_GENES) else "missing_target_genes",
    }
    if coverage["status"] != "ready" or len(common) < 8 or labels.loc[common].nunique() < 2:
        return pd.DataFrame(), coverage
    score = score_gpu_combo(expr.loc[common], spec)
    threshold = float(spec["locked_score"]["locked_threshold"])
    predictions = pd.DataFrame(
        {
            "cohort": cohort,
            "sample_id": common.astype(str),
            "true_response_label": labels.loc[common].to_numpy(dtype=int),
            "response_probability": score.reindex(common).to_numpy(dtype=float),
            "threshold": threshold,
            "model_name": spec["rescue_combo_id"],
        }
    )
    return predictions, coverage


def _baseline_predictions(cohort: str, bulk_dir: Path, cbio_dir: Path) -> pd.DataFrame:
    expr, metadata = _read_cohort(cohort, bulk_dir, cbio_dir)
    labels = endpoint_label_series(metadata.set_index("sample_id")["response_raw"], "strict_recist").dropna().astype(int)
    common = expr.index.intersection(labels.index)
    rows = []
    for model_name in EIGHT_SIGNATURES:
        genes = BASELINE_SIGNATURES.get(model_name, [model_name])
        raw = signature_score(expr.loc[common], genes)
        score = _zscore(pd.DataFrame({"score": raw}))["score"]
        for sample_id, value in score.items():
            rows.append(
                {
                    "cohort": cohort,
                    "sample_id": str(sample_id),
                    "true_response_label": int(labels.loc[sample_id]),
                    "response_probability": float(value),
                    "model_name": model_name,
                }
            )
    return pd.DataFrame(rows)


def _family_gate(target_predictions: pd.DataFrame, baseline_predictions: pd.DataFrame, n_bootstrap: int) -> pd.DataFrame:
    rows = []
    rng = np.random.default_rng(20260527)
    for group_id, cohorts in GROUPS.items():
        target = target_predictions[target_predictions["cohort"].astype(str).isin(cohorts)].copy()
        baseline = baseline_predictions[baseline_predictions["cohort"].astype(str).isin(cohorts)].copy()
        if target.empty or baseline.empty:
            continue
        target["key"] = target["cohort"].astype(str) + "::" + target["sample_id"].astype(str)
        baseline["key"] = baseline["cohort"].astype(str) + "::" + baseline["sample_id"].astype(str)
        y = target.drop_duplicates("key").set_index("key")["true_response_label"].astype(int)
        target_score = target.drop_duplicates("key").set_index("key")["response_probability"].astype(float)
        wide = baseline.pivot_table(index="key", columns="model_name", values="response_probability", aggfunc="first")
        common = y.index.intersection(target_score.index).intersection(wide.dropna().index)
        y = y.loc[common]
        target_score = target_score.loc[common]
        wide = wide.loc[common]
        if len(common) < 8 or y.nunique() < 2:
            continue
        baseline_aucs = {col: float(roc_auc_score(y, wide[col])) for col in wide.columns}
        deltas = []
        for _ in range(n_bootstrap):
            idx = rng.integers(0, len(y), len(y))
            yy = y.to_numpy(dtype=int)[idx]
            if len(np.unique(yy)) < 2:
                continue
            tt = target_score.to_numpy(dtype=float)[idx]
            bb = wide.to_numpy(dtype=float)[idx]
            deltas.append(float(roc_auc_score(yy, tt) - np.mean([roc_auc_score(yy, bb[:, col]) for col in range(bb.shape[1])])))
        arr = np.asarray(deltas)
        rows.append(
            {
                "group_id": group_id,
                "target_model": str(target["model_name"].iloc[0]),
                "baseline_family": "eight_strong_signatures",
                "n_samples": int(len(common)),
                "n_signatures": int(wide.shape[1]),
                "target_AUROC": float(roc_auc_score(y, target_score)),
                "family_mean_AUROC": float(np.mean(list(baseline_aucs.values()))),
                "best_signature": max(baseline_aucs, key=baseline_aucs.get),
                "best_signature_AUROC": float(max(baseline_aucs.values())),
                "delta_vs_family_mean": float(arr.mean()),
                "ci_low": float(np.quantile(arr, 0.025)),
                "ci_high": float(np.quantile(arr, 0.975)),
                "two_sided_p": float(min(1.0, 2.0 * min((arr <= 0).mean(), (arr >= 0).mean()))),
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["two_sided_fdr_q"] = benjamini_hochberg(result["two_sided_p"])
        result["claim_level"] = np.where(
            (result["target_AUROC"] >= 0.70) & (result["delta_vs_family_mean"] > 0) & (result["two_sided_fdr_q"] <= 0.05),
            "strict_external_family_FDR_supported_numeric_target_met",
            np.where(result["delta_vs_family_mean"] > 0, "family_point_estimate_only", "family_not_superior"),
        )
    return result


def run(bulk_dir: Path, cbio_dir: Path, package_dir: Path, out_dir: Path, n_bootstrap: int) -> dict[str, str]:
    spec = json.loads((package_dir / "gpu_bioprior_rescue_combo_scoring_spec.json").read_text(encoding="utf-8"))
    prediction_frames = []
    coverage_rows = []
    baseline_frames = []
    for cohort in ["CBIO_LIU_DFCI_2019_PRE", "GSE145996"]:
        predictions, coverage = _cohort_predictions(cohort, bulk_dir, cbio_dir, spec)
        coverage_rows.append(coverage)
        if not predictions.empty:
            prediction_frames.append(predictions)
            baseline_frames.append(_baseline_predictions(cohort, bulk_dir, cbio_dir))
    predictions = pd.concat(prediction_frames, ignore_index=True) if prediction_frames else pd.DataFrame()
    baseline_predictions = pd.concat(baseline_frames, ignore_index=True) if baseline_frames else pd.DataFrame()
    metric_rows = []
    for group_id, cohorts in GROUPS.items():
        frame = predictions[predictions["cohort"].astype(str).isin(cohorts)]
        if frame.empty:
            continue
        metrics = compute_binary_metrics(
            frame["true_response_label"].astype(int),
            frame["response_probability"].astype(float),
            threshold=float(frame["threshold"].iloc[0]),
        )
        metric_rows.append(
            {
                "group_id": group_id,
                "model_name": spec["rescue_combo_id"],
                "endpoint": "strict_recist",
                "n_samples": int(len(frame)),
                "n_responders": int(frame["true_response_label"].astype(int).sum()),
                "n_nonresponders": int((frame["true_response_label"].astype(int) == 0).sum()),
                **metrics,
                "selection_boundary": "locked_gpu_bioprior_package_no_cbio_label_fit",
            }
        )
    metrics = pd.DataFrame(metric_rows)
    family = _family_gate(predictions, baseline_predictions, n_bootstrap)
    coverage = pd.DataFrame(coverage_rows)

    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "metrics": out_dir / "cbioportal_gpu_bioprior_external_metrics.tsv",
        "predictions": out_dir / "cbioportal_gpu_bioprior_external_predictions.tsv",
        "baseline_predictions": out_dir / "cbioportal_gpu_bioprior_external_baseline_predictions.tsv",
        "family": out_dir / "cbioportal_gpu_bioprior_external_family_comparison.tsv",
        "coverage": out_dir / "cbioportal_gpu_bioprior_gene_coverage.tsv",
        "audit": out_dir / "CBIOPORTAL_GPU_BIOPRIOR_EXTERNAL_AUDIT.md",
    }
    metrics.to_csv(outputs["metrics"], sep="\t", index=False)
    predictions.to_csv(outputs["predictions"], sep="\t", index=False)
    baseline_predictions.to_csv(outputs["baseline_predictions"], sep="\t", index=False)
    family.to_csv(outputs["family"], sep="\t", index=False)
    coverage.to_csv(outputs["coverage"], sep="\t", index=False)
    liu = metrics[metrics["group_id"].eq("cbio_liu_dfci_only")]
    pooled = metrics[metrics["group_id"].eq("strict_cbio_liu_plus_gse145996")]
    lines = [
        "# cBioPortal GPU biological-prior external validation",
        "",
        "The frozen GPU lipid/PI3K rescue-combo package was scored on cBioPortal Liu/DFCI without using cBioPortal labels for gene, weight, threshold, or calibration selection.",
        "",
    ]
    if not liu.empty:
        lines.append(
            "CBIO_LIU_DFCI_2019_PRE AUROC={:.3f}, AUPRC={:.3f}, ECE={:.3f}.".format(
                float(liu.iloc[0]["AUROC"]), float(liu.iloc[0]["AUPRC"]), float(liu.iloc[0]["ECE"])
            )
        )
    if not pooled.empty:
        lines.append(
            "CBIO_LIU_DFCI_2019_PRE plus GSE145996 AUROC={:.3f}, AUPRC={:.3f}, ECE={:.3f}.".format(
                float(pooled.iloc[0]["AUROC"]), float(pooled.iloc[0]["AUPRC"]), float(pooled.iloc[0]["ECE"])
            )
        )
    outputs["audit"].write_text("\n".join(lines), encoding="utf-8")
    return {key: str(value) for key, value in outputs.items()}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bulk-dir", default="data/processed/bulk")
    parser.add_argument("--cbio-dir", default="data/processed/cbioportal_melanoma")
    parser.add_argument("--package-dir", default="deliverables/gpu_bioprior_rescue_combo_package_20260527")
    parser.add_argument("--out", default="results/cbioportal_gpu_bioprior_external_validation_20260527")
    parser.add_argument("--bootstrap", type=int, default=2000)
    args = parser.parse_args()
    outputs = run(ROOT / args.bulk_dir, ROOT / args.cbio_dir, ROOT / args.package_dir, ROOT / args.out, args.bootstrap)
    print(json.dumps(outputs, ensure_ascii=False))


if __name__ == "__main__":
    main()
