from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_submission_package_audit_passes(tmp_path: Path) -> None:
    out = tmp_path / "submission_readiness_audit.tsv"
    subprocess.run(
        [
            sys.executable,
            "scripts/validation/audit_submission_package.py",
            "--manuscript",
            "paper/Journal of Translational Medicine投稿/EcoNiche-Opt_JTM_Main_Manuscript.md",
            "--jtm-dir",
            "paper/Journal of Translational Medicine投稿",
            "--table-dir",
            "tables/article",
            "--source-data",
            "paper/Journal of Translational Medicine投稿/Additional_file_2_Source_Data.xlsx",
            "--out",
            str(out),
        ],
        cwd=ROOT,
        check=True,
    )
    assert out.exists()
