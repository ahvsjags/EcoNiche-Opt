from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


if __name__ == "__main__":
    subprocess.run([sys.executable, str(ROOT / "scripts/model/run_econiche.py"), "--config", "config/model_config.yml"], cwd=ROOT, check=True)
