from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd

CORE_PREDICTION_COLUMNS = {"sample_id", "patient_id", "cohort", "pred_prob"}
CORE_METRIC_COLUMNS = {"cohort", "model_name"}


def validate_tsv_columns(path: str | Path, required_columns: Iterable[str]) -> pd.DataFrame:
    tsv = Path(path)
    if not tsv.exists():
        return pd.DataFrame([{"path": str(tsv), "is_valid": False, "missing_columns": "FILE_NOT_FOUND"}])
    frame = pd.read_csv(tsv, sep="\t", nrows=5)
    missing = sorted(set(required_columns) - set(frame.columns))
    return pd.DataFrame([{"path": str(tsv), "is_valid": not missing, "missing_columns": ",".join(missing)}])


def validate_result_schema(results_dir: str | Path, demo: bool = False) -> pd.DataFrame:
    root = Path(results_dir)
    checks = [
        validate_tsv_columns(root / "lodo_predictions.tsv", CORE_PREDICTION_COLUMNS),
        validate_tsv_columns(root / "lodo_metrics.tsv", CORE_METRIC_COLUMNS),
        validate_tsv_columns(root / "econiche_module.tsv", {"state", "gene", "direction"}),
    ]
    return pd.concat(checks, ignore_index=True)


def assert_result_schema_valid(results_dir: str | Path, demo: bool = False) -> None:
    report = validate_result_schema(results_dir, demo=demo)
    failed = report[~report["is_valid"]]
    if not failed.empty:
        raise SystemExit("Result schema validation failed:\n" + failed.to_string(index=False))
