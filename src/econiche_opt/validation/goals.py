from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import yaml

EXPECTED_GOALS = [f"GOAL-{idx:03d}" for idx in range(81)]
ALLOWED_STATUSES = {
    "pending",
    "completed",
    "interface_completed",
    "blocked",
    "blocked_by_access",
    "failed",
    "result_pending",
    "unavailable_with_reason",
}


def load_goal_status(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def validate_goal_status(path: str | Path, demo_mode: bool = True) -> pd.DataFrame:
    data = load_goal_status(path)
    goals = data.get("goals", {})
    rows = []
    for goal_id in EXPECTED_GOALS:
        goal = goals.get(goal_id)
        if goal is None:
            rows.append({"goal_id": goal_id, "is_valid": False, "issue": "missing_goal"})
            continue
        status = str(goal.get("status", "")).strip()
        priority = str(goal.get("priority", "")).strip()
        issues = []
        if not status:
            issues.append("missing_status")
        elif status not in ALLOWED_STATUSES:
            issues.append("invalid_status")
        if priority in {"P0", "P1"} and not status:
            issues.append("p0_p1_missing_status")
        if status == "completed" and not (goal.get("files_changed") or goal.get("tests_run")):
            issues.append("completed_missing_evidence")
        if status in {"blocked", "blocked_by_access", "unavailable_with_reason"} and not goal.get("blocking_issues"):
            issues.append("blocked_missing_reason")
        if (not demo_mode) and status == "interface_completed":
            issues.append("interface_completed_only_allowed_in_demo_mode")
        rows.append({"goal_id": goal_id, "is_valid": not issues, "issue": ",".join(issues)})
    return pd.DataFrame(rows)


def assert_goal_status_valid(path: str | Path, demo_mode: bool = True) -> None:
    report = validate_goal_status(path, demo_mode=demo_mode)
    failed = report[~report["is_valid"]]
    if not failed.empty:
        raise SystemExit("Goal status validation failed:\n" + failed.to_string(index=False))
