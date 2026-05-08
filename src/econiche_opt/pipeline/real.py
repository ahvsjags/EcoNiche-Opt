from __future__ import annotations

import json
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from econiche_opt.data.download_plan import write_download_plan
from econiche_opt.data.registry import audit_accession, load_registry, validate_registry
from econiche_opt.validation.schema import validate_result_schema


@dataclass(frozen=True)
class RealPipelineConfig:
    root: Path
    registry: Path
    results_dir: Path
    out_dir: Path
    execute_download: bool = False
    execute_preprocess: bool = False
    execute_training: bool = False
    execute_secondary: bool = True
    strict_existing_results: bool = True


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _rel(root: Path, path: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def _stage(
    name: str,
    status: str,
    command: str = "",
    output: str = "",
    reason: str = "",
    started_utc: str | None = None,
    finished_utc: str | None = None,
) -> dict[str, object]:
    return {
        "stage": name,
        "status": status,
        "command": command,
        "output": output,
        "reason": reason,
        "started_utc": started_utc or _utc_now(),
        "finished_utc": finished_utc or _utc_now(),
    }


def _run(root: Path, name: str, command: list[str], output: str = "") -> dict[str, object]:
    started = _utc_now()
    command_text = " ".join(command)
    result = subprocess.run(command, cwd=root, text=True, capture_output=True)
    finished = _utc_now()
    if result.returncode == 0:
        return _stage(name, "PASS", command_text, output, "", started, finished)
    reason = (result.stderr or result.stdout or f"exit code {result.returncode}").strip().splitlines()[-1]
    return _stage(name, "FAILED", command_text, output, reason, started, finished)


def _validate_existing_results(cfg: RealPipelineConfig) -> dict[str, object]:
    started = _utc_now()
    report = validate_result_schema(cfg.results_dir)
    out = cfg.out_dir / "real_result_schema_validation.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(out, sep="\t", index=False)
    failed = report[~report["is_valid"]]
    if failed.empty:
        status = "PASS"
        reason = ""
    elif cfg.strict_existing_results:
        status = "FAILED"
        reason = "required real result schema files are missing or invalid"
    else:
        status = "RESULT_PENDING"
        reason = "real result schema files are missing or invalid"
    return _stage("validate_existing_real_results", status, "", _rel(cfg.root, out), reason, started, _utc_now())


def _write_registry_reports(cfg: RealPipelineConfig) -> list[dict[str, object]]:
    rows = []
    registry = load_registry(cfg.registry)
    validation = validate_registry(registry)
    validation_out = cfg.out_dir / "real_registry_validation.tsv"
    validation.to_csv(validation_out, sep="\t", index=False)
    rows.append(
        _stage(
            "validate_registry",
            "PASS" if validation["is_valid"].all() else "FAILED",
            "",
            _rel(cfg.root, validation_out),
            "" if validation["is_valid"].all() else "registry validation failed",
        )
    )

    audit = audit_accession(registry)
    audit_out = cfg.out_dir / "real_data_access_audit.tsv"
    audit.to_csv(audit_out, sep="\t", index=False)
    rows.append(_stage("audit_data_access", "PASS", "", _rel(cfg.root, audit_out)))

    plan_out = cfg.out_dir / "real_download_plan.tsv"
    write_download_plan(cfg.registry, plan_out)
    rows.append(_stage("download_plan", "PASS", "", _rel(cfg.root, plan_out), "dry-run plan only"))
    return rows


def _optional_stage(
    cfg: RealPipelineConfig,
    name: str,
    command: list[str],
    enabled: bool,
    reason: str,
    output: str = "",
) -> dict[str, object]:
    if enabled:
        return _run(cfg.root, name, command, output=output)
    return _stage(name, "SKIPPED", " ".join(command), output, reason)


def run_real_pipeline(cfg: RealPipelineConfig) -> pd.DataFrame:
    cfg.out_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    rows.extend(_write_registry_reports(cfg))

    rows.append(
        _optional_stage(
            cfg,
            "download_public_metadata",
            [sys.executable, "scripts/download/download_geo.py", "--registry", str(cfg.registry), "--metadata-only"],
            cfg.execute_download,
            "not requested; dry-run plan recorded instead",
            "data/metadata/geo_download_summary.tsv",
        )
    )

    for name, command in [
        ("preprocess_build_metadata", [sys.executable, "scripts/preprocess/build_metadata.py"]),
        ("preprocess_harmonize_labels", [sys.executable, "scripts/preprocess/harmonize_labels.py"]),
        ("preprocess_deduplicate_patients", [sys.executable, "scripts/preprocess/deduplicate_patients.py"]),
        ("preprocess_bulk", [sys.executable, "scripts/preprocess/preprocess_bulk.py"]),
    ]:
        rows.append(
            _optional_stage(
                cfg,
                name,
                command,
                cfg.execute_preprocess,
                "not requested; use --execute-preprocess after raw/metadata QC",
            )
        )

    training_commands = [
        ("train_econiche_real", [sys.executable, "scripts/model/run_econiche.py", "--config", "config/model_config.yml"]),
        ("score_baselines_real", [sys.executable, "scripts/baselines/run_baselines.py", "--config", "config/baselines.yml"]),
        ("train_ml_baselines_real", [sys.executable, "scripts/baselines/train_ml_baselines.py"]),
        ("response_composite_real", [sys.executable, "scripts/model/run_response_composite.py"]),
        ("run_lodo_real", [sys.executable, "scripts/benchmark/run_lodo.py"]),
        ("bootstrap_compare_real", [sys.executable, "scripts/benchmark/bootstrap_compare.py"]),
        ("calibration_real", [sys.executable, "scripts/benchmark/calibration.py"]),
        ("decision_curve_real", [sys.executable, "scripts/benchmark/decision_curve.py"]),
    ]
    for name, command in training_commands:
        rows.append(
            _optional_stage(
                cfg,
                name,
                command,
                cfg.execute_training,
                "not requested; existing real results are validated without overwriting",
            )
        )

    rows.append(_validate_existing_results(cfg))

    rows.append(
        _run(
            cfg.root,
            "public_cohort_smoke_tests",
            [
                sys.executable,
                "scripts/analysis/run_public_smoke_tests.py",
                "--out",
                "results/audit/public_cohort_smoke_tests.tsv",
            ],
            "results/audit/public_cohort_smoke_tests.tsv",
        )
    )
    rows.append(
        _run(
            cfg.root,
            "therapy_timepoint_stratification",
            [
                sys.executable,
                "scripts/analysis/run_stratification.py",
                "--out",
                "results/real_stratification/therapy_timepoint_stratification.tsv",
            ],
            "results/real_stratification/therapy_timepoint_stratification.tsv",
        )
    )

    secondary_commands = [
        ("xena_tcga_manifest", [sys.executable, "scripts/download/download_xena.py"]),
        ("deconvolution_baselines", [sys.executable, "scripts/preprocess/run_deconvolution.py", "--out-dir", "results/real"]),
        ("pancancer_transfer", [sys.executable, "scripts/benchmark/run_pan_cancer_transfer.py"]),
        ("survival_validation", [sys.executable, "scripts/benchmark/survival_analysis.py"]),
        ("single_cell_mapping", [sys.executable, "scripts/single_cell/map_modules_scrna.py"]),
        ("perturbation_prioritization", [sys.executable, "scripts/analysis/run_perturbation_prioritization.py"]),
        (
            "real_reproducibility_report",
            [sys.executable, "scripts/reporting/make_reproducibility_report.py", "--out", "paper/real_reproducibility_report.md"],
        ),
        ("real_manuscript", [sys.executable, "scripts/reporting/make_manuscript.py", "--out", "paper/real_manuscript.md"]),
    ]
    for name, command in secondary_commands:
        rows.append(
            _optional_stage(
                cfg,
                name,
                command,
                cfg.execute_secondary,
                "not requested; secondary validation skipped",
            )
        )

    manifest = pd.DataFrame(rows)
    manifest_out = cfg.out_dir / "pipeline_run_manifest.tsv"
    manifest.to_csv(manifest_out, sep="\t", index=False)
    summary = {
        "generated_utc": _utc_now(),
        "registry": _rel(cfg.root, cfg.registry),
        "results_dir": _rel(cfg.root, cfg.results_dir),
        "out_dir": _rel(cfg.root, cfg.out_dir),
        "status_counts": manifest["status"].value_counts().to_dict(),
        "failed_stages": manifest.loc[manifest["status"] == "FAILED", "stage"].tolist(),
    }
    (cfg.out_dir / "pipeline_run_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


def assert_real_pipeline_ok(manifest: pd.DataFrame) -> None:
    failed = manifest[manifest["status"] == "FAILED"]
    if not failed.empty:
        raise SystemExit("Real pipeline failed:\n" + failed[["stage", "reason"]].to_string(index=False))
