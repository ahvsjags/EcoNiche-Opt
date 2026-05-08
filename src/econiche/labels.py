from __future__ import annotations

import re
from typing import Any

import pandas as pd


RESPONDER = {
    "cr",
    "complete response",
    "complete_response",
    "pr",
    "partial response",
    "partial_response",
    "responder",
    "response",
    "r",
    "prcr",
    "pr/cr",
    "pr cr",
}
NON_RESPONDER = {
    "pd",
    "progressive disease",
    "progressive_disease",
    "progression",
    "non-responder",
    "non responder",
    "non_response",
    "nonresponse",
    "nr",
    "no response",
    "no_response",
    "non responder",
}
STABLE = {"sd", "stable disease", "stable_disease"}
DCB = {"dcb", "durable clinical benefit", "durable_clinical_benefit"}
NDB = {"ndb", "no durable clinical benefit", "non durable clinical benefit", "non_durable_clinical_benefit"}


def _clean_label(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip().lower()
    text = re.sub(r"[\s/]+", " ", text)
    text = text.replace("_", " ")
    return text


def harmonize_response(value: Any, endpoint: str = "primary_recist") -> int | None:
    """Map response labels to 0=responder/sensitive and 1=non-responder/resistant."""
    text = _clean_label(value)
    if not text:
        return None

    if endpoint == "clinical_benefit":
        if text in DCB:
            return 0
        if text in NDB:
            return 1

    if text in DCB:
        return 0
    if text in NDB:
        return 1
    if text in RESPONDER:
        return 0
    if text in NON_RESPONDER:
        return 1
    if text in STABLE:
        return None if endpoint == "strict_recist" else 1
    return None


def harmonize_metadata(
    metadata: pd.DataFrame,
    endpoint: str = "primary_recist",
    source_column: str = "response_raw",
    label_column: str = "label",
) -> pd.DataFrame:
    out = metadata.copy()
    if source_column not in out.columns:
        out[label_column] = pd.NA
        out["label_status"] = "missing_response_column"
        return out
    out[label_column] = out[source_column].map(lambda value: harmonize_response(value, endpoint=endpoint))
    out["label_status"] = out[label_column].map(lambda value: "parsed" if pd.notna(value) else "needs_manual_curation")
    return out
