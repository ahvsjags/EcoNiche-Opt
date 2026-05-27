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


def _run_setting(
    setting: dict,
    base_config: dict,
    config_dir: Path,
    processed_dir: Path,
    priors_path: Path,
    out_dir: Path,
    demo: bool,
) -> dict[str, object]:
    setting_name = setting["name"]
    setting_out = out_dir / "runs" / setting_name
    config_path = _write_config(base_config, setting.get("config", {}), config_dir / f"{setting_name}.yml")
    command = [
        sys.executable,
        str(ROOT / "scripts/model/run_econiche.py"),
        "--config",
        str(config_path.relative_to(ROOT)),
        "--processed-dir",
        str(processed_dir.relative_to(ROOT)),
        "--priors",
        str(priors_path.relative_to(ROOT)),
        "--out",
        str(setting_out.relative_to(ROOT)),
    ]
    if demo:
        command.append("--demo")
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    status = "computed_from_pipeline" if result.returncode == 0 else "FAILED"
    reason = "" if result.returncode == 0 else (result.stderr or result.stdout or "setting failed").strip().splitlines()[-1]
    return {
        "setting": setting_name,
        "status": status,
        "reason": reason,
        "command": " ".join(command),
        "out_dir": str(setting_out.relative_to(ROOT)),
        "config": json.dumps(setting.get("config", {}), sort_keys=True),
        "metrics_path": setting_out / "lodo_metrics.tsv",
    }


def _summarize_metrics(setting: str, metrics_path: Path, status: str, reason: str = "") -> list[dict[str, object]]:
    if status != "computed_from_pipeline" or not metrics_path.exists():
        return [
            {
                "setting": setting,
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
                "setting": setting,
                "metric": metric,
                "estimate": float(values.mean()) if values.notna().any() else pd.NA,
                "sd": float(values.std(ddof=0)) if values.notna().any() else pd.NA,
                "n_cohorts": int(values.notna().sum()),
                "status": status,
                "reason": reason,
            }
        )
    return rows


def _write_stability(metrics: pd.DataFrame, out: Path) -> None:
    computed = metrics[metrics["status"] == "computed_from_pipeline"].copy()
    rows = []
    for metric, group in computed.groupby("metric"):
        values = pd.to_numeric(group["estimate"], errors="coerce")
        rows.append(
            {
                "metric": metric,
                "n_settings": int(values.notna().sum()),
                "estimate_mean": float(values.mean()) if values.notna().any() else pd.NA,
                "estimate_sd": float(values.std(ddof=0)) if values.notna().any() else pd.NA,
                "estimate_min": float(values.min()) if values.notna().any() else pd.NA,
                "estimate_max": float(values.max()) if values.notna().any() else pd.NA,
                "status": "computed_from_pipeline" if values.notna().any() else "RESULT_PENDING",
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

    out_dir = ROOT / (args.out or ("results/demo_sensitivity" if args.demo else "results/real_sensitivity"))
    out_dir.mkdir(parents=True, exist_ok=True)
    processed_dir, n_cohorts = _filtered_processed_dir(ROOT / args.processed_dir, out_dir, args.demo)
    should_execute = args.demo or args.execute_real

    settings = [
        {"name": "baseline_seed_42", "config": {"random_state": 42}},
        {"name": "seed_7", "config": {"random_state": 7}},
        {"name": "seed_99", "config": {"random_state": 99}},
        {"name": "compact_3_genes", "config": {"max_genes_per_state": 3}},
        {"name": "expanded_min_6_genes", "config": {"min_genes_per_state": 6, "max_genes_per_state": 25}},
        {"name": "high_robustness_penalty", "config": {"robust_rho": 1.0}},
        {"name": "no_calibration_penalty", "config": {"w_ece": 0.0}},
    ]

    manifest_rows = []
    metric_rows = []
    if should_execute and n_cohorts >= 2:
        base_config = _load_config(ROOT / args.config)
        config_dir = out_dir / "configs"
        for setting in settings:
            manifest = _run_setting(
                setting,
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
                    setting["name"],
                    manifest["metrics_path"],
                    manifest["status"],
                    manifest["reason"],
                )
            )
    else:
        reason = (
            "real sensitivity rerun requires --execute-real and at least two curated non-demo processed cohorts"
            if not args.demo
            else "need at least two demo cohorts"
        )
        for setting in settings:
            manifest_rows.append(
                {
                    "setting": setting["name"],
                    "status": "RESULT_PENDING",
                    "reason": reason,
                    "command": "",
                    "out_dir": "",
                    "config": json.dumps(setting.get("config", {}), sort_keys=True),
                }
            )
            metric_rows.extend(_summarize_metrics(setting["name"], Path(), "RESULT_PENDING", reason))

    metrics = pd.DataFrame(metric_rows)
    manifest = pd.DataFrame(manifest_rows)
    metrics.to_csv(out_dir / "sensitivity_metrics.tsv", sep="\t", index=False)
    manifest.to_csv(out_dir / "sensitivity_manifest.tsv", sep="\t", index=False)
    _write_stability(metrics, out_dir / "sensitivity_stability.tsv")
    print(f"Wrote {out_dir / 'sensitivity_metrics.tsv'}")


if __name__ == "__main__":
    main()
