from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import yaml

REQUIRED_SOURCE_FIELDS = {"id", "type", "title", "source", "access_status"}


def load_source_registry(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {"sources": []}


def validate_source_registry(path: str | Path) -> pd.DataFrame:
    registry = load_source_registry(path)
    rows = []
    for source in registry.get("sources", []):
        missing = sorted(field for field in REQUIRED_SOURCE_FIELDS if not source.get(field))
        rows.append(
            {
                "id": source.get("id", "UNKNOWN"),
                "is_valid": not missing,
                "missing_fields": ",".join(missing),
            }
        )
    return pd.DataFrame(rows)
