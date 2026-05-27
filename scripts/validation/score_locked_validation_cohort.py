from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.metrics import compute_binary_metrics, decision_curve
from econiche_opt.deploy import read_expression_table
from econiche_opt.model.endpoint_modules import build_module_features, module_prior_score, sigmoid

ENDPOINT_LABEL_COLUMNS = {
    "primary_recist": "primary_recist_label",
    "strict_recist": "strict_recist_label",
    "clinical_benefit": "clinical_benefit_label",
}


def _read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    return pd.read_csv(path, sep=sep)


def _load_spec(package_dir: Path) -> dict[str, object]:
    spec_path = package_dir / "locked_scoring_spec.json"
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    hash_path = package_dir / "locked_scoring_spec.sha256"
    if hash_path.exists():
        expected = hash_path.read_text(encoding="utf-8").split()[0]
        observed = hashlib.sha256(spec_path.read_bytes()).hexdigest()
        if observed != expected:
            raise ValueError(f"locked_scoring_spec.json SHA256 mismatch: observed {observed}, expected {expected}")
    return spec


def _endpoint_thresholds(spec: dict[str, object]) -> dict[str, dict[str, float | str]]:
    out: dict[str, dict[str, float | str]] = {}
    for item in spec.get("endpoint_thresholds", []):
        if not isinstance(item, dict) or "endpoint" not in item:
            continue
        endpoint = str(item["endpoint"])
        out[endpoint] = {
            "threshold": float(item.get("threshold", 0.5)),
            "calibration": str(item.get("calibration", "raw_sigmoid")),
            "calibration_coef": float(item.get("calibration_coef", 1.0)),
            "calibration_intercept": float(item.get("calibration_intercept", 0.0)),
        }
    return out


def _audit_manifest(sample_manifest: pd.DataFrame, expression: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    required = [
        "sample_id",
        "sample_source",
        "baseline_status",
        "therapy",
        "qc_pass",
        "locked_validation_use_flag",
    ]
    for col in required:
        rows.append({"check": f"sample_manifest_column:{col}", "is_valid": col in sample_manifest.columns, "detail": ""})
    if "sample_id" not in sample_manifest.columns:
        return pd.DataFrame(rows)
    sample_ids = sample_manifest["sample_id"].astype(str)
    expression_ids = set(expression.index.astype(str))
    missing_expr = sorted(set(sample_ids) - expression_ids)
    rows.append({"check": "sample_ids_have_expression", "is_valid": not missing_expr, "detail": ",".join(missing_expr[:20])})
    if "sample_source" in sample_manifest.columns:
        values = sample_manifest["sample_source"].fillna("").astype(str).str.lower()
        rows.append({"check": "tumor_tissue_only", "is_valid": bool(values.str.contains("tumor").all()), "detail": ";".join(sorted(values.unique()))})
    if "baseline_status" in sample_manifest.columns:
        values = sample_manifest["baseline_status"].fillna("").astype(str).str.lower()
        rows.append({"check": "pretreatment_status", "is_valid": bool(values.str.contains("pretreatment|baseline|pre").all()), "detail": ";".join(sorted(values.unique()))})
    if "therapy" in sample_manifest.columns:
        values = sample_manifest["therapy"].fillna("").astype(str).str.lower()
        rows.append({"check": "anti_pd1_based_therapy", "is_valid": bool(values.str.contains("pd-1|pd1").all()), "detail": ";".join(sorted(values.unique()))})
    return pd.DataFrame(rows)


def _panel_gene_coverage(panel: pd.DataFrame, expression: pd.DataFrame) -> pd.DataFrame:
    available = set(expression.columns.astype(str))
    rows = []
    for gene, sub in panel.groupby("gene_symbol", sort=True):
        rows.append(
            {
                "gene_symbol": gene,
                "modules": ";".join(sorted(sub["module"].astype(str).unique())),
                "available": gene in available,
                "n_nonmissing_samples": int(expression[gene].notna().sum()) if gene in expression.columns else 0,
            }
        )
    return pd.DataFrame(rows)


def _probability_from_raw(raw_score: pd.Series, threshold_info: dict[str, float | str]) -> pd.Series:
    coef = float(threshold_info.get("calibration_coef", 1.0))
    intercept = float(threshold_info.get("calibration_intercept", 0.0))
    calibrated_logit = coef * raw_score.astype(float) + intercept
    return pd.Series(sigmoid(calibrated_logit), index=raw_score.index)


def _score_rows(module_features: pd.DataFrame, raw_score: pd.Series, thresholds: dict[str, dict[str, float | str]]) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for endpoint, info in thresholds.items():
        prob = _probability_from_raw(raw_score, info)
        threshold = float(info["threshold"])
        frame = pd.DataFrame(
            {
                "sample_id": raw_score.index.astype(str),
                "endpoint": endpoint,
                "raw_module_prior_score": raw_score.values,
                "response_probability": prob.values,
                "locked_threshold": threshold,
                "predicted_label": (prob >= threshold).astype(int).values,
                "calibration": info.get("calibration", "raw_sigmoid"),
                "calibration_coef": float(info.get("calibration_coef", 1.0)),
                "calibration_intercept": float(info.get("calibration_intercept", 0.0)),
            }
        )
        rows.append(frame)
    out = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    for module in module_features.columns:
        lookup = module_features[module].to_dict()
        out[f"module_score__{module}"] = out["sample_id"].map(lookup)
    return out


def _metrics(scores: pd.DataFrame, clinical: pd.DataFrame | None) -> tuple[pd.DataFrame, pd.DataFrame]:
    if clinical is None or clinical.empty or "sample_id" not in clinical.columns:
        return pd.DataFrame(), pd.DataFrame()
    clinical = clinical.copy()
    clinical["sample_id"] = clinical["sample_id"].astype(str)
    rows = []
    curve_rows = []
    for endpoint, label_col in ENDPOINT_LABEL_COLUMNS.items():
        if label_col not in clinical.columns:
            continue
        y = pd.to_numeric(clinical.set_index("sample_id")[label_col], errors="coerce").dropna().astype(int)
        subset = scores[scores["endpoint"] == endpoint].set_index("sample_id")
        common = y.index.intersection(subset.index)
        if len(common) < 4 or y.loc[common].nunique() < 2:
            rows.append({"endpoint": endpoint, "status": "RESULT_PENDING", "reason": "fewer than four labeled samples or one response class"})
            continue
        prob = subset.loc[common, "response_probability"].astype(float)
        threshold = float(subset.loc[common, "locked_threshold"].iloc[0])
        metrics = compute_binary_metrics(y.loc[common], prob, threshold=threshold)
        rows.append(
            {
                "endpoint": endpoint,
                "status": "completed",
                "n_samples": int(len(common)),
                "n_responders": int(y.loc[common].sum()),
                "n_nonresponders": int((y.loc[common] == 0).sum()),
                "locked_threshold": threshold,
                **metrics,
            }
        )
        curve = decision_curve(y.loc[common], prob).assign(endpoint=endpoint)
        curve_rows.append(curve)
    return pd.DataFrame(rows), pd.concat(curve_rows, ignore_index=True) if curve_rows else pd.DataFrame()


def score_locked_validation_cohort(
    package_dir: str | Path,
    expression_path: str | Path,
    sample_manifest_path: str | Path,
    out_dir: str | Path,
    clinical_annotation_path: str | Path | None = None,
    transpose: bool = False,
) -> dict[str, Path]:
    package_dir = Path(package_dir)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    spec = _load_spec(package_dir)
    panel = pd.read_csv(package_dir / "locked_panel_genes.tsv", sep="\t")
    expression = read_expression_table(expression_path, transpose=transpose)
    sample_manifest = _read_table(sample_manifest_path)
    sample_manifest["sample_id"] = sample_manifest["sample_id"].astype(str)

    audit = _audit_manifest(sample_manifest, expression)
    audit.to_csv(out / "locked_validation_manifest_audit.tsv", sep="\t", index=False)

    eligible = sample_manifest["sample_id"].astype(str)
    if "locked_validation_use_flag" in sample_manifest.columns:
        flag = sample_manifest["locked_validation_use_flag"].fillna("").astype(str).str.lower()
        use_mask = flag.isin(["1", "true", "yes", "y", "use", "include", "locked"])
        if use_mask.any():
            eligible = sample_manifest.loc[use_mask, "sample_id"].astype(str)
    common = [sample for sample in eligible if sample in expression.index]
    if len(common) < 2:
        raise ValueError("At least two validation samples with expression are required for locked panel scoring")
    expression = expression.loc[common].copy()

    module_features, module_coverage = build_module_features(expression)
    raw_score = module_prior_score(module_features)
    thresholds = _endpoint_thresholds(spec)
    if not thresholds:
        raise ValueError("locked_scoring_spec.json does not contain endpoint_thresholds")
    scores = _score_rows(module_features, raw_score, thresholds)

    clinical = _read_table(clinical_annotation_path) if clinical_annotation_path else None
    metric_table, curve_table = _metrics(scores, clinical)

    gene_coverage = _panel_gene_coverage(panel, expression)
    module_coverage = module_coverage.assign(
        coverage_fraction=lambda x: x["n_genes_available"].astype(float) / x["n_genes_defined"].replace(0, np.nan).astype(float)
    )
    module_features.insert(0, "sample_id", module_features.index.astype(str))

    outputs = {
        "scores": out / "locked_validation_scores.tsv",
        "module_scores": out / "locked_validation_module_scores.tsv",
        "module_coverage": out / "locked_validation_module_coverage.tsv",
        "gene_coverage": out / "locked_validation_gene_coverage.tsv",
        "manifest_audit": out / "locked_validation_manifest_audit.tsv",
        "metrics": out / "locked_validation_metrics.tsv",
        "decision_curve": out / "locked_validation_decision_curve.tsv",
    }
    scores.to_csv(outputs["scores"], sep="\t", index=False)
    module_features.to_csv(outputs["module_scores"], sep="\t", index=False)
    module_coverage.to_csv(outputs["module_coverage"], sep="\t", index=False)
    gene_coverage.to_csv(outputs["gene_coverage"], sep="\t", index=False)
    metric_table.to_csv(outputs["metrics"], sep="\t", index=False)
    curve_table.to_csv(outputs["decision_curve"], sep="\t", index=False)

    manifest = pd.DataFrame(
        [
            {"artifact": key, "path": str(path), "status": "written", "n_bytes": path.stat().st_size if path.exists() else 0}
            for key, path in outputs.items()
        ]
    )
    manifest.to_csv(out / "locked_validation_output_manifest.tsv", sep="\t", index=False)
    outputs["output_manifest"] = out / "locked_validation_output_manifest.tsv"
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description="Score an independent cohort with the frozen EcoNiche-Opt validation package.")
    parser.add_argument("--package-dir", default="deliverables/prospective_validation")
    parser.add_argument("--expression", required=True)
    parser.add_argument("--sample-manifest", required=True)
    parser.add_argument("--clinical-annotation")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--transpose", action="store_true")
    args = parser.parse_args()
    outputs = score_locked_validation_cohort(
        package_dir=args.package_dir,
        expression_path=args.expression,
        sample_manifest_path=args.sample_manifest,
        clinical_annotation_path=args.clinical_annotation,
        out_dir=args.out_dir,
        transpose=args.transpose,
    )
    print(pd.DataFrame([{"artifact": key, "path": str(path)} for key, path in outputs.items()]).to_string(index=False))


if __name__ == "__main__":
    main()
