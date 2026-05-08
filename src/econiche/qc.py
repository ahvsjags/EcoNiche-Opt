from __future__ import annotations

import hashlib
from typing import Any

import pandas as pd


def conservative_patient_id(row: pd.Series) -> str:
    parts = [
        str(row.get("accession", row.get("cohort", ""))),
        str(row.get("source_name", "")),
        str(row.get("title", "")),
        str(row.get("sample_id", "")),
    ]
    digest = hashlib.sha1("|".join(parts).encode("utf-8")).hexdigest()[:10]
    return f"conservative_{digest}"


def ensure_patient_id(metadata: pd.DataFrame) -> pd.DataFrame:
    out = metadata.copy()
    if "patient_id" not in out.columns:
        out["patient_id"] = pd.NA
    missing = out["patient_id"].isna() | (out["patient_id"].astype(str).str.strip() == "")
    out.loc[missing, "patient_id"] = out.loc[missing].apply(conservative_patient_id, axis=1)
    out["patient_id_confidence"] = "high"
    out.loc[missing, "patient_id_confidence"] = "low"
    return out


def _normalized_timepoint(value: Any) -> str:
    return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")


def split_by_timepoint_priority(metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    meta = ensure_patient_id(metadata)
    timepoint = meta.get("timepoint", pd.Series(["pretreatment"] * len(meta), index=meta.index)).map(_normalized_timepoint)
    progression_mask = timepoint.str.contains("progression|acquired|resistance", regex=True, na=False)
    secondary_mask = timepoint.str.contains("on_treatment|ontreatment|during|post", regex=True, na=False) & ~progression_mask
    primary_mask = timepoint.str.contains("pretreatment|baseline|pre_treatment|pre", regex=True, na=False)
    if not primary_mask.any():
        primary_mask = ~(secondary_mask | progression_mask)

    primary = meta.loc[primary_mask].copy()
    if not primary.empty:
        primary = primary.sort_values(["patient_id", "sample_id"] if "sample_id" in primary.columns else ["patient_id"])
        primary = primary.drop_duplicates("patient_id", keep="first")
    secondary = meta.loc[secondary_mask].copy()
    progression = meta.loc[progression_mask].copy()
    return primary.reset_index(drop=True), secondary.reset_index(drop=True), progression.reset_index(drop=True)


def check_patient_leakage(folds: dict[str, dict[str, pd.DataFrame]]) -> pd.DataFrame:
    rows = []
    for fold_name, fold in folds.items():
        train_ids = set(fold["train"].get("patient_id", pd.Series(dtype=str)).dropna().astype(str))
        test_ids = set(fold["test"].get("patient_id", pd.Series(dtype=str)).dropna().astype(str))
        overlap = sorted(train_ids & test_ids)
        rows.append(
            {
                "fold": fold_name,
                "ok": len(overlap) == 0,
                "overlap_count": len(overlap),
                "overlapping_patient_ids": ",".join(overlap),
            }
        )
    return pd.DataFrame(rows)
