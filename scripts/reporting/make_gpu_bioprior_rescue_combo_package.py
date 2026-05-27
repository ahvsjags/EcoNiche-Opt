from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from scripts.analysis.run_gpu_bioprior_rescue_combo_search import (  # noqa: E402
    PRIMARY_COHORTS,
    STRICT_EXTERNAL_COHORTS,
    build_transforms,
    labels_for_endpoint,
    load_bulk,
)


GENE_ROWS = [
    {
        "gene_symbol": "MAP4K1",
        "component": "base_positive",
        "direction": "response_high",
        "role": "immune_signaling_balance_positive_gene",
        "interpretation": "Positive arm of the MAP4K1 minus TBX3/AXL tumor-immune balance axis.",
    },
    {
        "gene_symbol": "TBX3",
        "component": "base_negative",
        "direction": "resistance_high",
        "role": "tumor_dedifferentiation_resistance_gene",
        "interpretation": "Negative arm of the tumor-intrinsic resistance balance axis.",
    },
    {
        "gene_symbol": "AXL",
        "component": "base_negative",
        "direction": "resistance_high",
        "role": "IPRES_dedifferentiation_extension_gene",
        "interpretation": "Negative arm extending the dedifferentiation/resistance axis.",
    },
    {
        "gene_symbol": "PLA2G2D",
        "component": "lipid_pi3k_rescue_component",
        "direction": "training_response_oriented",
        "role": "lipid_pi3k_biological_prior_gene",
        "interpretation": "Primary-LODO-selected lipid/PI3K-prior rescue component.",
    },
]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_tsv(path: Path, required: set[str], label: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing {label}: {path}")
    frame = pd.read_csv(path, sep="\t")
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{label} missing required columns: {sorted(missing)}")
    return frame


def _selected_row(selection: pd.DataFrame) -> pd.Series:
    rows = selection[selection["selection_id"].astype(str).str.contains("gpu_", na=False)].copy()
    if rows.empty:
        raise ValueError("GPU biological-prior selection table lacks a gpu_* selection row")
    row = rows.sort_values("strict_external_AUROC", ascending=False).iloc[0]
    if str(row["selection_boundary"]) != "candidate_and_weight_selected_by_primary_lodo_only_with_biological_prior":
        raise ValueError("GPU rescue-combo selection boundary is not primary-LODO-only")
    if float(row["strict_external_AUROC"]) < 0.70:
        raise ValueError("GPU rescue-combo does not meet the strict external AUROC >=0.70 freeze threshold")
    return row


def _parse_candidate(candidate: str) -> dict[str, object]:
    match = re.fullmatch(r"(?P<base>[0-9.]+)\*base\+(?P<component>[0-9.]+)\*(?P<method>[a-z]+)__ (?P<gene>[A-Za-z0-9]+)", candidate)
    if match is None:
        match = re.fullmatch(
            r"(?P<base>[0-9.]+)\*base\+(?P<component>[0-9.]+)\*(?P<method>[a-z]+)__(?P<gene>[A-Za-z0-9]+)",
            candidate,
        )
    if match is None:
        raise ValueError(f"Unsupported GPU rescue-combo candidate string: {candidate}")
    return {
        "weight_base": float(match.group("base")),
        "weight_component": float(match.group("component")),
        "component_method": match.group("method"),
        "component_gene": match.group("gene"),
    }


def _component_sign(processed_dir: Path, gene: str, method: str) -> tuple[int, float]:
    X_by_cohort, metadata_by_cohort = load_bulk(processed_dir)
    primary_y = labels_for_endpoint(X_by_cohort, metadata_by_cohort, PRIMARY_COHORTS, "primary_recist")
    strict_y = labels_for_endpoint(
        X_by_cohort,
        metadata_by_cohort,
        [*PRIMARY_COHORTS, *STRICT_EXTERNAL_COHORTS],
        "strict_recist",
    )
    transforms = build_transforms(X_by_cohort, primary_y, strict_y, ["MAP4K1", "TBX3", "AXL", gene])
    values: list[float] = []
    labels: list[int] = []
    for cohort, y in primary_y.items():
        series = transforms[cohort][method][gene].reindex(y.index)
        values.extend(series.to_numpy(dtype=float).tolist())
        labels.extend(y.to_numpy(dtype=int).tolist())
    corr = float(np.corrcoef(np.asarray(values), np.asarray(labels))[0, 1])
    if not np.isfinite(corr):
        corr = 0.0
    return (1 if corr >= 0 else -1), corr


def build_package(
    selection_path: Path,
    external_predictions_path: Path,
    family_gate_path: Path,
    processed_dir: Path,
    out: Path,
    release_tag: str,
) -> dict[str, object]:
    selection = _read_tsv(
        selection_path,
        {
            "selection_id",
            "candidate",
            "prior",
            "transform_policy",
            "device",
            "gpu_name",
            "selection_boundary",
            "primary_AUROC",
            "primary_AUPRC",
            "primary_balanced_accuracy",
            "strict_external_AUROC",
            "strict_external_AUPRC",
            "strict_external_balanced_accuracy",
            "strict_external_ECE",
            "family_mean_AUROC",
            "delta_vs_family_mean",
            "two_sided_fdr_q",
            "claim_level",
        },
        "GPU rescue-combo selection",
    )
    external_predictions = _read_tsv(
        external_predictions_path,
        {"candidate", "threshold", "cohort", "sample_id", "response_probability", "true_response_label"},
        "GPU rescue-combo external predictions",
    )
    family_gate = _read_tsv(
        family_gate_path,
        {"target_model", "family_mean_AUROC", "best_signature", "best_signature_AUROC", "two_sided_fdr_q"},
        "GPU rescue-combo family gate",
    )
    selected = _selected_row(selection)
    parsed = _parse_candidate(str(selected["candidate"]))
    sign, corr = _component_sign(processed_dir, str(parsed["component_gene"]), str(parsed["component_method"]))
    target_predictions = external_predictions[external_predictions["candidate"].astype(str).eq(str(selected["candidate"]))]
    thresholds = sorted({round(float(value), 12) for value in target_predictions["threshold"].dropna().unique()})
    if len(thresholds) != 1:
        raise ValueError(f"Expected one locked threshold for selected GPU rescue-combo, found {thresholds}")

    out.mkdir(parents=True, exist_ok=True)
    genes = pd.DataFrame(GENE_ROWS)
    genes.to_csv(out / "gpu_bioprior_rescue_combo_genes.tsv", sep="\t", index=False)
    selection.to_csv(out / "gpu_bioprior_rescue_combo_evidence.tsv", sep="\t", index=False)
    family_gate.to_csv(out / "gpu_bioprior_rescue_combo_family_gate.tsv", sep="\t", index=False)

    spec = {
        "rescue_combo_id": "EcoNiche-Opt-GPU-LipidPI3K-RescueCombo",
        "model_status": "frozen_no_leakage_lipid_pi3k_rescue_combo_for_future_validation",
        "release_tag": release_tag,
        "created": date.today().isoformat(),
        "intended_context": "pretreatment melanoma tumor tissue before anti-PD-1 or anti-PD-1-based therapy",
        "endpoint": "CR/PR responders versus SD/PD nonresponders when RECIST categories are available",
        "selection_rule": "candidate and blend weight selected by primary melanoma LODO only within a biological lipid/PI3K prior candidate set",
        "locked_candidate": str(selected["candidate"]),
        "locked_score": {
            "base_axis": "MAP4K1 minus mean(TBX3, AXL)",
            "base_transform": "0.95 * cohort robust z-score axis + 0.05 * cohort z-score axis, then min-max normalize within the scored cohort",
            "component_gene": parsed["component_gene"],
            "component_transform": parsed["component_method"],
            "component_training_direction_sign": sign,
            "component_primary_training_correlation": corr,
            "component_transform_definition": "cohort robust z-score = (x - cohort median) / (1.4826 * cohort MAD), clipped to [-5, 5], then multiplied by the frozen primary-training direction sign and min-max normalized within the scored cohort",
            "weight_base": parsed["weight_base"],
            "weight_component": parsed["weight_component"],
            "score_formula": "minmax(weight_base * minmax(0.95 * minmax(rz(MAP4K1)-mean(rz(TBX3),rz(AXL))) + 0.05 * minmax(z(MAP4K1)-mean(z(TBX3),z(AXL)))) + weight_component * minmax(sign * transform(PLA2G2D)))",
            "score_orientation": "higher values indicate a more response-like immune-ecological state",
            "locked_threshold": thresholds[0],
            "threshold_rule": "threshold selected on primary melanoma training data only, then applied without refitting to strict external cohorts",
        },
        "performance_evidence": {
            "primary_AUROC": float(selected["primary_AUROC"]),
            "primary_AUPRC": float(selected["primary_AUPRC"]),
            "primary_balanced_accuracy": float(selected["primary_balanced_accuracy"]),
            "strict_external_AUROC": float(selected["strict_external_AUROC"]),
            "strict_external_AUPRC": float(selected["strict_external_AUPRC"]),
            "strict_external_balanced_accuracy": float(selected["strict_external_balanced_accuracy"]),
            "strict_external_ECE": float(selected["strict_external_ECE"]),
            "family_mean_AUROC": float(selected["family_mean_AUROC"]),
            "delta_vs_family_mean": float(selected["delta_vs_family_mean"]),
            "two_sided_fdr_q": float(selected["two_sided_fdr_q"]),
            "claim_level": str(selected["claim_level"]),
            "device": str(selected["device"]),
            "gpu_name": str(selected["gpu_name"]),
            "selection_boundary": str(selected["selection_boundary"]),
            "best_baseline_signature": str(family_gate.iloc[0]["best_signature"]),
            "best_baseline_signature_AUROC": float(family_gate.iloc[0]["best_signature_AUROC"]),
        },
        "no_leakage_rules": [
            "Strict external labels must not enter feature selection.",
            "Strict external labels must not enter transform-policy selection.",
            "Strict external labels must not enter candidate-weight selection.",
            "Strict external labels must not enter thresholding or calibration.",
            "Future controlled cohorts must be scored after this JSON, gene table, threshold, and reporting plan are frozen.",
        ],
        "required_reporting": [
            "AUROC",
            "AUPRC",
            "balanced_accuracy",
            "sensitivity",
            "specificity",
            "PPV",
            "NPV",
            "Brier",
            "ECE",
            "calibration_bins",
            "decision_curve",
            "paired family comparison with FDR claim gate",
        ],
    }

    spec_path = out / "gpu_bioprior_rescue_combo_scoring_spec.json"
    spec_path.write_text(json.dumps(spec, indent=2, ensure_ascii=False), encoding="utf-8")
    checksum = _sha256(spec_path)
    (out / "gpu_bioprior_rescue_combo_scoring_spec.sha256").write_text(
        f"{checksum}  {spec_path.name}\n",
        encoding="utf-8",
    )
    readme = [
        "# GPU Biological-Prior Rescue Combo",
        "",
        "This package freezes the GPU-audited lipid/PI3K biological-prior rescue combo for locked melanoma external scoring.",
        "",
        "## Locked Candidate",
        "",
        f"- Candidate: `{selected['candidate']}`.",
        f"- Genes: MAP4K1, TBX3, AXL, {parsed['component_gene']}.",
        f"- Threshold: {thresholds[0]:.6f}.",
        f"- Selection boundary: {selected['selection_boundary']}.",
        "",
        "## Current Evidence",
        "",
        f"- Primary melanoma LODO AUROC: {float(selected['primary_AUROC']):.3f}.",
        f"- Strict melanoma external AUROC: {float(selected['strict_external_AUROC']):.3f}.",
        f"- Strict melanoma external AUPRC: {float(selected['strict_external_AUPRC']):.3f}.",
        f"- Family-comparison FDR q: {float(selected['two_sided_fdr_q']):.3f}.",
        f"- GPU device: {selected['gpu_name']}.",
        "",
        "External cohorts must be scored without changing genes, weights, transform policy, threshold, or calibration.",
        "",
    ]
    (out / "README.md").write_text("\n".join(readme), encoding="utf-8")
    return {
        "out": str(out),
        "rescue_combo_id": spec["rescue_combo_id"],
        "locked_candidate": spec["locked_candidate"],
        "strict_external_AUROC": spec["performance_evidence"]["strict_external_AUROC"],
        "strict_external_AUPRC": spec["performance_evidence"]["strict_external_AUPRC"],
        "two_sided_fdr_q": spec["performance_evidence"]["two_sided_fdr_q"],
        "locked_threshold": spec["locked_score"]["locked_threshold"],
        "component_direction_sign": sign,
        "sha256": checksum,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--selection",
        default="results/gpu_bioprior_rescue_combo_search_robust_20260527/gpu_bioprior_rescue_combo_selection.tsv",
    )
    parser.add_argument(
        "--external-predictions",
        default="results/gpu_bioprior_rescue_combo_search_robust_20260527/gpu_bioprior_external_predictions.tsv",
    )
    parser.add_argument(
        "--family-gate",
        default="results/gpu_bioprior_rescue_combo_search_robust_20260527/gpu_bioprior_external_family_gate.tsv",
    )
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--out", default="deliverables/gpu_bioprior_rescue_combo_package_20260527")
    parser.add_argument("--release-tag", default="v0.3.3-gpu-bioprior-rescue-combo-20260527")
    args = parser.parse_args()

    summary = build_package(
        ROOT / args.selection,
        ROOT / args.external_predictions,
        ROOT / args.family_gate,
        ROOT / args.processed_dir,
        ROOT / args.out,
        args.release_tag,
    )
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
