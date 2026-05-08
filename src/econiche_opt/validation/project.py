from __future__ import annotations

from pathlib import Path

import pandas as pd

from econiche_opt.data.registry import load_registry, validate_registry
from econiche_opt.reporting.citations import validate_source_registry
from econiche_opt.validation.goals import validate_goal_status
from econiche_opt.validation.schema import validate_result_schema


def _check_exists(path: Path) -> dict[str, object]:
    return {"check": f"exists:{path.as_posix()}", "is_valid": path.exists(), "detail": ""}


def validate_project(root: str | Path = ".", mode: str = "demo") -> pd.DataFrame:
    root_path = Path(root)
    rows = [
        _check_exists(root_path / "README.md"),
        _check_exists(root_path / "AGENTS.md"),
        _check_exists(root_path / "config" / "data_registry.yml"),
        _check_exists(root_path / "docs" / "goal_status.yml"),
        _check_exists(root_path / "docs" / "reproducibility" / "no_fabrication_policy.md"),
    ]

    goal_path = root_path / "docs" / "goal_status.yml"
    if goal_path.exists():
        goal_report = validate_goal_status(goal_path, demo_mode=(mode == "demo"))
        rows.append(
            {
                "check": "goal_status",
                "is_valid": bool(goal_report["is_valid"].all()),
                "detail": f"{int(goal_report['is_valid'].sum())}/{len(goal_report)} valid",
            }
        )

    registry_path = root_path / "config" / "data_registry.yml"
    if registry_path.exists():
        registry_report = validate_registry(load_registry(registry_path))
        rows.append(
            {
                "check": "registry",
                "is_valid": bool(registry_report["is_valid"].all()),
                "detail": f"{int(registry_report['is_valid'].sum())}/{len(registry_report)} valid",
            }
        )

    results_dir = root_path / "results" / ("demo" if mode == "demo" else "real")
    schema_report = validate_result_schema(results_dir, demo=(mode == "demo"))
    rows.append(
        {
            "check": "result_schema",
            "is_valid": bool(schema_report["is_valid"].all()),
            "detail": f"{int(schema_report['is_valid'].sum())}/{len(schema_report)} valid",
        }
    )

    source_path = root_path / "config" / "source_registry.yml"
    if source_path.exists():
        source_report = validate_source_registry(source_path)
        rows.append(
            {
                "check": "source_registry",
                "is_valid": bool(source_report["is_valid"].all()),
                "detail": f"{int(source_report['is_valid'].sum())}/{len(source_report)} valid",
            }
        )

    return pd.DataFrame(rows)


def assert_project_valid(root: str | Path = ".", mode: str = "demo") -> None:
    report = validate_project(root, mode=mode)
    failed = report[~report["is_valid"]]
    if not failed.empty:
        raise SystemExit("Project validation failed:\n" + failed.to_string(index=False))
