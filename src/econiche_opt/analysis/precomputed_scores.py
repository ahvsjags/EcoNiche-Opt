from __future__ import annotations

from pathlib import Path

import pandas as pd

REQUIRED_COLUMNS = {"sample_id", "model", "score"}


def import_precomputed_scores(path: str | Path, out: str | Path | None = None) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    missing = sorted(REQUIRED_COLUMNS - set(frame.columns))
    if missing:
        raise ValueError(f"Missing required precomputed score columns: {', '.join(missing)}")
    frame = frame.copy()
    frame["source"] = Path(path).name
    if out is not None:
        out_path = Path(out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(out_path, sep="\t", index=False)
    return frame
