from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_prospective_validation_package_is_locked_and_auditable() -> None:
    subprocess.run(
        [sys.executable, "scripts/validation/validate_prospective_package.py", "--package-dir", "deliverables/prospective_validation"],
        cwd=ROOT,
        check=True,
    )
