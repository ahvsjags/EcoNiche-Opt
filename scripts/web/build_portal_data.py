from __future__ import annotations

"""Build the data bundle used by the static EcoNiche-Opt research portal.

The portal is intentionally a generated view of registered repository files. It
does not contain hand-entered clinical results or replacement data.
"""

import argparse
import csv
import json
import shutil
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "deliverables" / "prospective_validation"
BMC = ROOT / "paper" / "BMC_Bioinformatics_revision_20260808"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def number(value: str | None) -> float | int | None:
    if value is None or value == "":
        return None
    parsed = float(value)
    return int(parsed) if parsed.is_integer() else parsed


def build_bundle() -> dict[str, Any]:
    spec = read_json(PACKAGE / "locked_scoring_spec.json")
    panel_rows = read_tsv(PACKAGE / "locked_panel_genes.tsv")
    modules: list[dict[str, Any]] = []
    for module in spec.get("module_weights", {}):
        rows = [row for row in panel_rows if row.get("module") == module]
        modules.append(
            {
                "id": module,
                "label": module.replace("_", " ").title(),
                "weight": spec["module_weights"][module],
                "direction": rows[0].get("score_direction", "response_high") if rows else "response_high",
                "genes": [row["gene_symbol"] for row in rows],
            }
        )

    review_rows = read_tsv(BMC / "reviewer_comment_resolution_matrix.tsv")
    audit_rows = read_tsv(BMC / "bmc_revision_pre_submission_audit.tsv")
    figure_rows = read_tsv(ROOT / "figures" / "article" / "figure_manifest.tsv")

    benchmark_rows = read_tsv(ROOT / "results" / "endpoint_modules_heuristic_core_locked_gpu" / "endpoint_module_summary.tsv")
    benchmark = []
    for row in benchmark_rows:
        if row.get("model_name") != "EcoNiche-Opt-HeuristicEcology":
            continue
        benchmark.append(
            {
                "endpoint": row.get("endpoint", ""),
                "stratum": row.get("stratum", ""),
                "n_samples": number(row.get("n_samples")),
                "n_cohorts": number(row.get("n_cohorts")),
                "pooled_AUROC": number(row.get("pooled_AUROC")),
                "pooled_AUPRC": number(row.get("pooled_AUPRC")),
                "pooled_balanced_accuracy": number(row.get("pooled_balanced_accuracy")),
                "pooled_ECE": number(row.get("pooled_ECE")),
                "evaluation_modes": row.get("evaluation_modes", ""),
            }
        )

    external_rows = read_tsv(ROOT / "results" / "locked_external_panel_validation_calibrated_20260519" / "locked_external_metrics.tsv")
    external = []
    seen: set[tuple[str, str]] = set()
    for row in external_rows:
        if row.get("model_name") != "EcoNiche-Opt-HeuristicEcology-LockedPanel":
            continue
        key = (row.get("endpoint", ""), row.get("cohort", ""))
        if key in seen:
            continue
        seen.add(key)
        external.append(
            {
                "endpoint": row.get("endpoint", ""),
                "cohort": row.get("cohort", ""),
                "analysis_type": row.get("analysis_type", ""),
                "n_samples": number(row.get("n_samples")),
                "n_responders": number(row.get("n_responders")),
                "n_nonresponders": number(row.get("n_nonresponders")),
                "AUROC": number(row.get("AUROC")),
                "AUPRC": number(row.get("AUPRC")),
                "balanced_accuracy": number(row.get("balanced_accuracy")),
                "Brier": number(row.get("Brier")),
                "ECE": number(row.get("ECE")),
                "threshold_source": row.get("threshold_source", ""),
            }
        )

    audit_pass = sum(row.get("is_valid", "").lower() == "true" for row in audit_rows)
    return {
        "generated_by": "scripts/web/build_portal_data.py",
        "source_files": [
            "deliverables/prospective_validation/locked_scoring_spec.json",
            "deliverables/prospective_validation/locked_panel_genes.tsv",
            "paper/BMC_Bioinformatics_revision_20260808/reviewer_comment_resolution_matrix.tsv",
            "paper/BMC_Bioinformatics_revision_20260808/bmc_revision_pre_submission_audit.tsv",
            "figures/article/figure_manifest.tsv",
            "results/endpoint_modules_heuristic_core_locked_gpu/endpoint_module_summary.tsv",
            "results/locked_external_panel_validation_calibrated_20260519/locked_external_metrics.tsv",
        ],
        "model": {
            "name": spec.get("model_name"),
            "status": spec.get("model_status"),
            "release_tag": spec.get("release_tag"),
            "context": spec.get("intended_validation_context"),
            "unsupported_context": spec.get("unsupported_context_without_new_training"),
            "input_format": spec.get("input_matrix_format"),
            "transform": spec.get("input_transform"),
            "formula": spec.get("score_formula"),
            "panel_unique_genes": spec.get("panel_unique_genes"),
            "panel_gene_module_rows": spec.get("panel_gene_module_rows"),
            "threshold_source": spec.get("training_threshold_source"),
            "allowed_claim": spec.get("allowed_primary_claim"),
            "forbidden_claim": spec.get("forbidden_claim"),
        },
        "modules": modules,
        "endpoints": spec.get("endpoint_thresholds", []),
        "comparators": spec.get("predeclared_comparators", []),
        "required_metrics": spec.get("required_primary_report_metrics", []),
        "review": {
            "total": len(review_rows),
            "resolved": sum(row.get("status", "").lower() == "resolved" for row in review_rows),
            "rows": review_rows,
        },
        "audit": {
            "total": len(audit_rows),
            "passed": audit_pass,
            "failed": len(audit_rows) - audit_pass,
            "rows": audit_rows,
        },
        "figures": {
            "count": len(figure_rows),
            "rows": figure_rows,
        },
        "benchmark": benchmark,
        "external": external,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="web/data/portal_manifest.json")
    args = parser.parse_args()
    output = ROOT / args.out
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(build_bundle(), indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    shutil.copyfile(PACKAGE / "locked_scoring_spec.json", output.parent / "locked_scoring_spec.json")
    shutil.copyfile(PACKAGE / "locked_panel_genes.tsv", output.parent / "locked_panel_genes.tsv")
    print(output)


if __name__ == "__main__":
    main()
