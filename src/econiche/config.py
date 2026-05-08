from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from econiche.module import EcoNicheConfig


def load_yaml(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def load_model_config(path: str | Path) -> EcoNicheConfig:
    data = load_yaml(path)
    allowed = set(EcoNicheConfig.__dataclass_fields__)
    return EcoNicheConfig(**{key: value for key, value in data.items() if key in allowed})
