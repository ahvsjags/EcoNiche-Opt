from pathlib import Path

import yaml

from econiche_opt.validation.goals import EXPECTED_GOALS, validate_goal_status


def test_goal_validator_requires_all_goals(tmp_path: Path):
    path = tmp_path / "goals.yml"
    path.write_text(yaml.safe_dump({"goals": {"GOAL-000": {"priority": "P0", "status": "completed"}}}), encoding="utf-8")
    report = validate_goal_status(path)
    assert len(report) == len(EXPECTED_GOALS)
    assert not report["is_valid"].all()


def test_goal_validator_accepts_completed_evidence(tmp_path: Path):
    goals = {
        goal_id: {
            "title": goal_id,
            "priority": "P0" if goal_id == "GOAL-000" else "P2",
            "status": "completed",
            "files_changed": ["README.md"],
            "tests_run": [],
            "blocking_issues": [],
            "notes": "",
        }
        for goal_id in EXPECTED_GOALS
    }
    path = tmp_path / "goals.yml"
    path.write_text(yaml.safe_dump({"goals": goals}), encoding="utf-8")
    report = validate_goal_status(path)
    assert report["is_valid"].all()
