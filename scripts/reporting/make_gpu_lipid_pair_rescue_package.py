from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.io import load_processed_bulk
from econiche_opt.model.endpoint_modules import endpoint_label_series


PRIMARY_COHORTS = ["GSE91061", "GSE78220", "PRJEB23709_PD1_PRE"]
GENES = ["MAP4K1", "TBX3", "AXL", "PLA2G2D", "PIK3CD"]
PATTERN = re.compile(
    r"(?P<weight_base>[0-9.]+)\*base\+(?P<weight_pair>[0-9.]+)\*"
    r"\((?P<mix1>[0-9.]+)\*(?P<method1>rz|z|pct)__PLA2G2D\+"
    r"(?P<mix2>[0-9.]+)\*(?P<method2>rz|z|pct)__PIK3CD\)"
)


def _direction(processed_dir: Path, method: str, gene: str) -> float:
    X_by_cohort, _, metadata_by_cohort = load_processed_bulk(processed_dir)
    values = []
    labels = []
    for cohort in PRIMARY_COHORTS:
        metadata = metadata_by_cohort[cohort]
        if "sample_id" in metadata.columns:
            metadata = metadata.set_index("sample_id")
        y = endpoint_label_series(metadata["response_raw"], "primary_recist").dropna().astype(int)
        common = X_by_cohort[cohort].index.intersection(y.index)
        x = X_by_cohort[cohort].loc[common, [gene]].apply(pd.to_numeric, errors="coerce")
        if method == "z":
            score = ((x[gene] - x[gene].mean()) / (x[gene].std(ddof=0) + 1e-6)).fillna(0.0)
        elif method == "rz":
            median = x[gene].median()
            mad = (x[gene] - median).abs().median() + 1e-6
            score = ((x[gene] - median) / (1.4826 * mad)).clip(-5.0, 5.0).fillna(0.0)
        else:
            score = x[gene].rank(pct=True).fillna(0.5)
        values.append(score)
        labels.append(y.loc[common])
    v = pd.concat(values)
    y = pd.concat(labels).reindex(v.index).astype(int)
    return 1.0 if float(v[y == 1].mean() - v[y == 0].mean()) >= 0.0 else -1.0


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run(
    selection_path: Path,
    external_metrics_path: Path,
    external_predictions_path: Path,
    family_gate_path: Path,
    processed_dir: Path,
    out_dir: Path,
    release_tag: str,
) -> dict[str, Path]:
    selection = pd.read_csv(selection_path, sep="\t").iloc[0]
    metrics = pd.read_csv(external_metrics_path, sep="\t")
    predictions = pd.read_csv(external_predictions_path, sep="\t")
    family = pd.read_csv(family_gate_path, sep="\t")
    candidate = str(selection["candidate"])
    match = PATTERN.fullmatch(candidate)
    if not match:
        raise ValueError(f"Unsupported lipid-pair candidate format: {candidate}")
    parsed = match.groupdict()
    method1 = parsed["method1"]
    method2 = parsed["method2"]
    direction1 = _direction(processed_dir, method1, "PLA2G2D")
    direction2 = _direction(processed_dir, method2, "PIK3CD")
    threshold = float(predictions["threshold"].dropna().iloc[0])
    out_dir.mkdir(parents=True, exist_ok=True)

    gene_roles = pd.DataFrame(
        [
            {"gene_symbol": "MAP4K1", "role": "base_response_axis_positive", "transform": "rz,z", "locked_in_model": True},
            {"gene_symbol": "TBX3", "role": "base_resistance_axis_negative", "transform": "rz,z", "locked_in_model": True},
            {"gene_symbol": "AXL", "role": "base_resistance_axis_negative", "transform": "rz,z", "locked_in_model": True},
            {"gene_symbol": "PLA2G2D", "role": "lipid_immunoregulatory_component", "transform": method1, "locked_in_model": True},
            {"gene_symbol": "PIK3CD", "role": "PI3K_delta_immune_component", "transform": method2, "locked_in_model": True},
        ]
    )
    evidence = pd.concat(
        [
            selection.to_frame().T.assign(evidence_type="selection"),
            metrics.assign(evidence_type="external_metrics"),
            family.assign(evidence_type="family_gate"),
        ],
        ignore_index=True,
        sort=False,
    )
    spec = {
        "model_id": "EcoNiche-Opt-GPU-LipidPI3K-PairRescue",
        "release_tag": release_tag,
        "candidate": candidate,
        "selection_boundary": selection.get(
            "selection_boundary",
            "pair_candidate_and_weights_selected_by_primary_lodo_only_lipid_pi3k_prior_no_external_labels",
        ),
        "training_endpoint": "primary_recist",
        "training_cohorts": PRIMARY_COHORTS,
        "locked_score": {
            "base_genes": {"positive": "MAP4K1", "negative_mean": ["TBX3", "AXL"]},
            "base_formula": "minmax(0.95*robust_z(MAP4K1-mean(TBX3,AXL))+0.05*z(MAP4K1-mean(TBX3,AXL)))",
            "weight_base": float(parsed["weight_base"]),
            "weight_pair_component": float(parsed["weight_pair"]),
            "pair_components": [
                {
                    "gene": "PLA2G2D",
                    "method": method1,
                    "mix_weight": float(parsed["mix1"]),
                    "training_direction_sign": direction1,
                },
                {
                    "gene": "PIK3CD",
                    "method": method2,
                    "mix_weight": float(parsed["mix2"]),
                    "training_direction_sign": direction2,
                },
            ],
            "locked_threshold": threshold,
            "normalization": "cohort-wise z, robust-z, or within-cohort percentile as specified; min-max normalization is applied inside each scoring cohort without labels",
        },
        "claim_boundary": {
            "external_or_holdout_labels_used_for_training": False,
            "external_or_holdout_labels_used_for_feature_selection": False,
            "external_or_holdout_labels_used_for_threshold_or_calibration": False,
            "claim_language": "family-level FDR-supported improvement only where the paired family gate passes",
        },
    }
    outputs = {
        "genes": out_dir / "gpu_lipid_pair_rescue_genes.tsv",
        "evidence": out_dir / "gpu_lipid_pair_rescue_evidence.tsv",
        "family_gate": out_dir / "gpu_lipid_pair_rescue_family_gate.tsv",
        "spec": out_dir / "gpu_lipid_pair_rescue_scoring_spec.json",
        "checksum": out_dir / "gpu_lipid_pair_rescue_scoring_spec.sha256",
        "readme": out_dir / "README.md",
    }
    gene_roles.to_csv(outputs["genes"], sep="\t", index=False)
    evidence.to_csv(outputs["evidence"], sep="\t", index=False)
    family.to_csv(outputs["family_gate"], sep="\t", index=False)
    outputs["spec"].write_text(json.dumps(spec, indent=2, ensure_ascii=False), encoding="utf-8")
    outputs["checksum"].write_text(_sha256(outputs["spec"]) + "\n", encoding="utf-8")
    readme = [
        "# EcoNiche-Opt GPU lipid/PI3K pair rescue package",
        "",
        f"Release tag: `{release_tag}`.",
        "",
        "This package freezes the primary-LODO-selected component-dominant lipid/PI3K rescue score.",
        "The scoring rule uses MAP4K1/TBX3/AXL as the base tumor-immune balance axis and a locked PLA2G2D/PIK3CD pair component.",
        "External labels are not used for gene, weight, threshold, or calibration selection.",
        "",
        f"Selected candidate: `{candidate}`.",
        f"Locked threshold: `{threshold:.12g}`.",
    ]
    outputs["readme"].write_text("\n".join(readme), encoding="utf-8")
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selection", default="results/gpu_lipid_pair_rescue_component_dominant_ba_guardrail_20260528/gpu_lipid_pair_rescue_selection.tsv")
    parser.add_argument("--external-metrics", default="results/gpu_lipid_pair_rescue_component_dominant_ba_guardrail_20260528/gpu_lipid_pair_external_metrics.tsv")
    parser.add_argument("--external-predictions", default="results/gpu_lipid_pair_rescue_component_dominant_ba_guardrail_20260528/gpu_lipid_pair_external_predictions.tsv")
    parser.add_argument("--family-gate", default="results/gpu_lipid_pair_rescue_component_dominant_ba_guardrail_20260528/gpu_lipid_pair_external_family_gate.tsv")
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--out", default="deliverables/gpu_lipid_pair_rescue_package_20260528")
    parser.add_argument("--release-tag", default="v0.3.4-gpu-lipid-pair-rescue-20260528")
    args = parser.parse_args()
    outputs = run(
        ROOT / args.selection,
        ROOT / args.external_metrics,
        ROOT / args.external_predictions,
        ROOT / args.family_gate,
        ROOT / args.processed_dir,
        ROOT / args.out,
        args.release_tag,
    )
    print(json.dumps({key: str(value) for key, value in outputs.items()}, ensure_ascii=False))


if __name__ == "__main__":
    main()
