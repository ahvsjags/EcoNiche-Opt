from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.io import load_processed_bulk
from econiche.priors import load_priors, make_default_cell_state_priors

METRIC_COLUMNS = ["AUROC", "AUPRC", "balanced_accuracy", "MCC", "F1", "ECE", "Brier"]


def _load_config(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def _write_config(base: dict, updates: dict, out: Path) -> Path:
    data = dict(base)
    data.update(updates)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    return out


def _write_zero_priors(base_priors: Path, processed_dir: Path, out: Path) -> Path:
    if base_priors.exists():
        priors = load_priors(str(base_priors))
        zero = priors * 0.0
    else:
        X_by_cohort, _, _ = load_processed_bulk(processed_dir)
        common = sorted(set.intersection(*(set(frame.columns) for frame in X_by_cohort.values()))) if X_by_cohort else []
        zero = make_default_cell_state_priors(common) * 0.0
    out.parent.mkdir(parents=True, exist_ok=True)
    zero.to_csv(out, sep="\t")
    return out


def _filtered_processed_dir(processed_dir: Path, out_dir: Path, demo: bool) -> tuple[Path, int]:
    X_by_cohort, _, _ = load_processed_bulk(processed_dir)
    if demo:
        cohorts = sorted(cohort for cohort in X_by_cohort if cohort.startswith("demo_cohort_"))
    else:
        cohorts = sorted(cohort for cohort in X_by_cohort if not cohort.startswith("demo_cohort_"))
    filtered = out_dir / "processed_inputs"
    if filtered.exists():
        shutil.rmtree(filtered)
    filtered.mkdir(parents=True, exist_ok=True)
    for cohort in cohorts:
        for suffix in [".expr.tsv", ".metadata.tsv"]:
            src = processed_dir / f"{cohort}{suffix}"
            if src.exists():
                shutil.copyfile(src, filtered / src.name)
    common = processed_dir / "common_genes_primary.txt"
    if common.exists():
        shutil.copyfile(common, filtered / common.name)
    return filtered, len(cohorts)


def _summarize_metrics(variant: str, metrics_path: Path, status: str, reason: str = "") -> list[dict[str, object]]:
    if status != "computed_from_pipeline" or not metrics_path.exists():
        return [
            {
                "ablation": variant,
                "metric": metric,
                "estimate": pd.NA,
                "sd": pd.NA,
                "n_cohorts": 0,
                "status": status,
                "reason": reason,
            }
            for metric in METRIC_COLUMNS
        ]
    metrics = pd.read_csv(metrics_path, sep="\t")
    rows = []
    for metric in METRIC_COLUMNS:
        values = pd.to_numeric(metrics.get(metric), errors="coerce")
        rows.append(
            {
                "ablation": variant,
                "metric": metric,
                "estimate": float(values.mean()) if values.notna().any() else pd.NA,
                "sd": float(values.std(ddof=0)) if values.notna().any() else pd.NA,
                "n_cohorts": int(values.notna().sum()),
                "status": status,
                "reason": reason,
            }
        )
    return rows


def _run_variant(
    variant: dict,
    base_config: dict,
    config_dir: Path,
    processed_dir: Path,
    priors_path: Path,
    out_dir: Path,
    demo: bool,
) -> dict[str, object]:
    variant_name = variant["name"]
    variant_out = out_dir / "runs" / variant_name
    config_path = _write_config(base_config, variant.get("config", {}), config_dir / f"{variant_name}.yml")
    active_priors = priors_path
    if variant.get("zero_priors"):
        active_priors = _write_zero_priors(priors_path, processed_dir, config_dir / f"{variant_name}.zero_priors.tsv")
    command = [
        sys.executable,
        str(ROOT / "scripts/model/run_econiche.py"),
        "--config",
        str(config_path.relative_to(ROOT)),
        "--processed-dir",
        str(processed_dir.relative_to(ROOT)),
        "--priors",
        str(active_priors.relative_to(ROOT)),
        "--out",
        str(variant_out.relative_to(ROOT)),
    ]
    if demo:
        command.append("--demo")
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    status = "computed_from_pipeline" if result.returncode == 0 else "FAILED"
    reason = "" if result.returncode == 0 else (result.stderr or result.stdout or "variant failed").strip().splitlines()[-1]
    return {
        "ablation": variant_name,
        "status": status,
        "reason": reason,
        "command": " ".join(command),
        "out_dir": str(variant_out.relative_to(ROOT)),
        "config": json.dumps(variant.get("config", {}), sort_keys=True),
        "zero_priors": bool(variant.get("zero_priors")),
        "metrics_path": variant_out / "lodo_metrics.tsv",
    }


def _write_delta(metrics: pd.DataFrame, out: Path) -> None:
    estimates = metrics[metrics["status"] == "computed_from_pipeline"].pivot_table(
        index="ablation", columns="metric", values="estimate", aggfunc="first"
    )
    rows = []
    if "full_model" in estimates.index:
        baseline = estimates.loc["full_model"]
        for variant, values in estimates.iterrows():
            for metric in estimates.columns:
                rows.append(
                    {
                        "ablation": variant,
                        "metric": metric,
                        "delta_vs_full_model": values.get(metric, pd.NA) - baseline.get(metric, pd.NA),
                    }
                )
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--execute-real", action="store_true")
    parser.add_argument("--config", default="config/model_config.yml")
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--priors", default="data/priors/cell_state_priors.tsv")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    out_dir = ROOT / (args.out or ("results/demo_ablation" if args.demo else "results/real_ablation"))
    out_dir.mkdir(parents=True, exist_ok=True)
    processed_dir, n_cohorts = _filtered_processed_dir(ROOT / args.processed_dir, out_dir, args.demo)
    should_execute = args.demo or args.execute_real

    variants = [
        {"name": "full_model", "config": {}},
        {"name": "zero_cell_state_priors", "config": {}, "zero_priors": True},
        {
            "name": "no_biological_objective_terms",
            "config": {"w_cell_specificity": 0.0, "w_pathway": 0.0, "w_network": 0.0, "w_lr": 0.0},
        },
        {"name": "no_size_penalty", "config": {"w_size": 0.0}},
        {"name": "compact_module", "config": {"max_genes_per_state": 3}},
    ]

    manifest_rows = []
    metric_rows = []
    if should_execute and n_cohorts >= 2:
        base_config = _load_config(ROOT / args.config)
        config_dir = out_dir / "configs"
        for variant in variants:
            manifest = _run_variant(
                variant,
                base_config,
                config_dir,
                processed_dir,
                ROOT / args.priors,
                out_dir,
                args.demo,
            )
            manifest_rows.append({k: v for k, v in manifest.items() if k != "metrics_path"})
            metric_rows.extend(
                _summarize_metrics(
                    variant["name"],
                    manifest["metrics_path"],
                    manifest["status"],
                    manifest["reason"],
                )
            )
    else:
        reason = (
            "real ablation rerun requires --execute-real and at least two curated non-demo processed cohorts"
            if not args.demo
            else "need at least two demo cohorts"
        )
        for variant in variants:
            manifest_rows.append(
                {
                    "ablation": variant["name"],
                    "status": "RESULT_PENDING",
                    "reason": reason,
                    "command": "",
                    "out_dir": "",
                    "config": json.dumps(variant.get("config", {}), sort_keys=True),
                    "zero_priors": bool(variant.get("zero_priors")),
                }
            )
            metric_rows.extend(_summarize_metrics(variant["name"], Path(), "RESULT_PENDING", reason))

    metrics = pd.DataFrame(metric_rows)
    manifest = pd.DataFrame(manifest_rows)
    metrics.to_csv(out_dir / "ablation_metrics.tsv", sep="\t", index=False)
    manifest.to_csv(out_dir / "ablation_manifest.tsv", sep="\t", index=False)
    _write_delta(metrics, out_dir / "ablation_delta_vs_full.tsv")
    print(f"Wrote {out_dir / 'ablation_metrics.tsv'}")


if __name__ == "__main__":
    main()
