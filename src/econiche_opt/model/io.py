from __future__ import annotations

from pathlib import Path
from typing import Any

import joblib


def save_model(model: Any, path: str | Path) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(model, out)
    return out


def load_model(path: str | Path) -> Any:
    return joblib.load(Path(path))
