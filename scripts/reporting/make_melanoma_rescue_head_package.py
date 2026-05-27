from __future__ import annotations

import argparse
import hashlib
import json
from datetime import date
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]


GENE_ROWS = [
    {
        "gene_symbol": "MAP4K1",
        "direction": "response_high",
        "role": "tumor_immune_balance_positive_gene",
        "interpretation": "immune-signaling axis gene retained by primary melanoma LODO selection",
    },
    {
        "gene_symbol": "TBX3",
        "direction": "resistance_high",
        "role": "tumor_dedifferentiation_resistance_gene",
        "interpretation": "tumor-intrinsic resistance axis gene retained by primary melanoma LODO selection",
    },
    {
        "gene_symbol": "AXL",
        "direction": "resistance_high",
        "role": "IPRES_dedifferentiation_extension_gene",
        "interpretation": "resistance-extension gene retained by primary melanoma LODO selection",
    },
]


def _read_selection(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing MAP4K1-TBX3 selection table: {path}")
    selection = pd.read_csv(path, sep="\t")
    required = {
        "selection_id",
        "claim_boundary",
        "method",
        "axis",
        "primary_AUROC",
        "primary_AUPRC",
        "primary_balanced_accuracy",
        "strict_external_AUROC",
        "strict_external_AUPRC",
        "strict_external_balanced_accuracy",
    }
    missing = required - set(selection.columns)
    if missing:
        raise ValueError(f"Selection table missing required columns: {sorted(missing)}")
    return selection


def _primary_selected(selection: pd.DataFrame) -> pd.Series:
    rows = selection[selection["selection_id"].astype(str).eq("primary_selected_candidate")]
    if rows.empty:
        raise ValueError("Selection table lacks primary_selected_candidate row")
    row = rows.iloc[0]
    if str(row["claim_boundary"]) != "selected_by_primary_lodo_only_not_by_external":
        raise ValueError("Primary-selected rescue head must be selected by primary LODO only")
    return row


def _stress_best(selection: pd.DataFrame) -> pd.Series | None:
    rows = selection[selection["selection_id"].astype(str).eq("current_external_stress_best")]
    if rows.empty:
        return None
    return rows.iloc[0]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_package(selection_path: Path, out: Path, release_tag: str) -> dict[str, object]:
    selection = _read_selection(selection_path)
    primary = _primary_selected(selection)
    stress = _stress_best(selection)
    out.mkdir(parents=True, exist_ok=True)

    genes = pd.DataFrame(GENE_ROWS)
    genes.to_csv(out / "melanoma_rescue_head_genes.tsv", sep="\t", index=False)
    selection.to_csv(out / "melanoma_rescue_head_evidence.tsv", sep="\t", index=False)

    spec = {
        "rescue_head_id": "EcoNiche-Opt-Melanoma-MAP4K1-TBX3-AXL",
        "model_status": "frozen_extension_for_future_locked_external_validation",
        "release_tag": release_tag,
        "created": date.today().isoformat(),
        "intended_context": "pretreatment melanoma tumor tissue before anti-PD-1 or anti-PD-1-based therapy",
        "endpoint": "CR/PR responders versus SD/PD nonresponders when RECIST categories are available",
        "selection_rule": "choose the MAP4K1-TBX3/AXL transform using primary melanoma LODO only; do not use current or future external labels for feature, transform, threshold, or calibration selection",
        "locked_transform": {
            "method": str(primary["method"]),
            "axis": str(primary["axis"]),
            "positive_genes": ["MAP4K1"],
            "negative_genes": ["TBX3", "AXL"],
            "transform_definition": "For each external cohort, rank each locked gene across samples to percentiles, average positive-gene percentiles, subtract the average negative-gene percentiles, then min-max scale the resulting sample score within the scored cohort.",
            "score_orientation": "higher values indicate a more response-like tumor-immune balance",
        },
        "primary_development_evidence": {
            "selection_boundary": str(primary["claim_boundary"]),
            "primary_AUROC": float(primary["primary_AUROC"]),
            "primary_AUPRC": float(primary["primary_AUPRC"]),
            "primary_balanced_accuracy": float(primary["primary_balanced_accuracy"]),
            "current_strict_external_AUROC_for_primary_selected_head": float(primary["strict_external_AUROC"]),
            "current_strict_external_AUPRC_for_primary_selected_head": float(primary["strict_external_AUPRC"]),
            "current_strict_external_balanced_accuracy_for_primary_selected_head": float(
                primary["strict_external_balanced_accuracy"]
            ),
        },
        "stress_screen_not_for_locked_selection": None,
        "no_leakage_rules": [
            "External labels must not enter feature selection.",
            "External labels must not enter transform selection.",
            "External labels must not enter thresholding or calibration.",
            "Current external stress-screen rows are reported only for weakness diagnosis and cannot define the locked rescue head.",
            "Controlled or unavailable external data remain ACCESS_RESTRICTED or RESULT_PENDING until imported through a registered pipeline.",
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
            "paired baseline comparison with FDR claim gate",
        ],
    }
    if stress is not None:
        spec["stress_screen_not_for_locked_selection"] = {
            "selection_id": str(stress["selection_id"]),
            "claim_boundary": str(stress["claim_boundary"]),
            "method": str(stress["method"]),
            "axis": str(stress["axis"]),
            "primary_AUROC": float(stress["primary_AUROC"]),
            "strict_external_AUROC": float(stress["strict_external_AUROC"]),
        }

    spec_path = out / "melanoma_rescue_head_scoring_spec.json"
    spec_path.write_text(json.dumps(spec, indent=2, ensure_ascii=False), encoding="utf-8")
    (out / "melanoma_rescue_head_scoring_spec.sha256").write_text(
        f"{_sha256(spec_path)}  {spec_path.name}\n",
        encoding="utf-8",
    )

    readme = [
        "# Melanoma MAP4K1-TBX3/AXL Rescue Head",
        "",
        "This package freezes the primary-selected melanoma rescue head for future locked external validation.",
        "It is an extension evidence package, not a replacement for the 62-gene EcoNiche-Opt locked panel.",
        "",
        "## Locked Rule",
        "",
        "The frozen rule is the primary-LODO-selected `cohort_gene_percentile` MAP4K1 minus TBX3/AXL axis.",
        "External cohorts may be scored only after the transform, genes, endpoint, and reporting plan are fixed.",
        "",
        "## Current Evidence Boundary",
        "",
        f"- Primary melanoma LODO AUROC: {float(primary['primary_AUROC']):.3f}.",
        f"- Current strict melanoma external AUROC for this primary-selected head: {float(primary['strict_external_AUROC']):.3f}.",
        "- The current external stress-screen row is diagnostic only and must not be used as a locked selection claim.",
        "",
        "## Claim Boundary",
        "",
        "Until newly obtained controlled external cohorts are scored with this frozen package, this artifact supports model-development and validation-readiness claims only. It does not close the strict external AUROC >=0.70 target.",
        "",
    ]
    (out / "README.md").write_text("\n".join(readme), encoding="utf-8")

    return {
        "out": str(out),
        "rescue_head_id": spec["rescue_head_id"],
        "locked_method": spec["locked_transform"]["method"],
        "locked_axis": spec["locked_transform"]["axis"],
        "primary_AUROC": spec["primary_development_evidence"]["primary_AUROC"],
        "current_strict_external_AUROC": spec["primary_development_evidence"][
            "current_strict_external_AUROC_for_primary_selected_head"
        ],
        "sha256": _sha256(spec_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--selection",
        default="results/map4k1_tbx3_transform_audit_20260527/map4k1_tbx3_transform_selection.tsv",
    )
    parser.add_argument("--out", default="deliverables/melanoma_rescue_head_package_20260527")
    parser.add_argument("--release-tag", default="v0.3.2-melanoma-rescue-head-20260527")
    args = parser.parse_args()

    summary = build_package(ROOT / args.selection, ROOT / args.out, args.release_tag)
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
