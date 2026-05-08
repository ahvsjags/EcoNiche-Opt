from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path
from typing import Any

import pandas as pd


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.name.endswith(".gz") else path.open("r", encoding="utf-8", errors="replace")


def _clean(value: str) -> str:
    return value.strip().strip('"')


def parse_series_matrix_metadata(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    sample_fields: dict[str, list[str]] = {}
    characteristics: list[list[str]] = []
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            row = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
            key = row[0].replace("!Sample_", "")
            values = [_clean(value) for value in row[1:]]
            if key == "characteristics_ch1":
                characteristics.append(values)
            else:
                sample_fields[key] = values
    n = len(sample_fields.get("geo_accession", sample_fields.get("title", [])))
    rows = []
    for idx in range(n):
        row: dict[str, Any] = {"matrix_file": path.name}
        for key, values in sample_fields.items():
            row[key] = values[idx] if idx < len(values) else pd.NA
        chars = []
        for values in characteristics:
            if idx < len(values):
                value = values[idx]
                chars.append(value)
                if ":" in value:
                    key, val = value.split(":", 1)
                    norm_key = re.sub(r"[^a-z0-9]+", "_", key.strip().lower()).strip("_")
                    row[norm_key] = val.strip()
        row["characteristics_ch1"] = "|".join(chars)
        rows.append(row)
    return pd.DataFrame(rows)


def infer_patient_id(row: pd.Series) -> str:
    for key in ["patient_id", "patient", "patient id", "patient_id_raw", "samples", "tumor_of_origin"]:
        normalized = key.replace(" ", "_")
        if normalized in row and pd.notna(row[normalized]) and str(row[normalized]).strip():
            return str(row[normalized]).strip()
    title = str(row.get("title", ""))
    match = re.search(r"(Pt\d+|Patient\s*\d+|sample[_\s-]*\d+|BACI\d+|D\d+|N\d+|RCC-\d+|Mel\d+)", title, flags=re.I)
    if match:
        value = match.group(1)
        if re.match(r"Patient\s*\d+", value, flags=re.I):
            return re.sub(r"\s+", "", value)
        return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")
    description = str(row.get("description", ""))
    match = re.search(r"(RCC-\d+|S\d+|SL\d+|Pt\d+)", description, flags=re.I)
    if match:
        return match.group(1)
    sample = str(row.get("geo_accession", row.get("sample_id", "")))
    return sample


def infer_response(row: pd.Series) -> str | None:
    explicit_values = []
    for key, value in row.items():
        key_text = str(key).lower()
        if any(token in key_text for token in ["response", "benefit", "recist", "best_resp", "best_response", "io_response"]) and pd.notna(value):
            explicit_values.append(str(value))
    for value in explicit_values:
        text = value.lower().strip()
        if text in {"r", "responder", "response"}:
            return "R"
        if text in {"nr", "non-responder", "non responder", "no_response", "no response"}:
            return "NR"
        if text in {"cr", "complete response"}:
            return "CR"
        if text in {"pr", "partial response", "prcr", "pr/cr"}:
            return "PR"
        if text in {"pd", "progressive disease"}:
            return "PD"
        if text in {"sd", "stable disease"}:
            return "SD"
        if text in {"dcb", "durable clinical benefit"}:
            return "DCB"
        if text in {"ndb", "no durable clinical benefit"}:
            return "NDB"
    texts = explicit_values + [str(row.get("source_name_ch1", "")), str(row.get("title", ""))]
    text = " | ".join(texts).lower()
    if "complete response" in text or re.search(r"\bcr\b", text):
        return "CR"
    if "partial response" in text or "prcr" in text or re.search(r"\bpr\b", text):
        return "PR"
    if "progressive disease" in text or re.search(r"\bpd\b", text):
        return "PD"
    if "stable disease" in text or re.search(r"\bsd\b", text):
        return "SD"
    if "no_response" in text or "no response" in text or re.search(r"\bnr\b", text):
        return "NR"
    if "dcb" in text or "durable clinical benefit" in text:
        return "DCB"
    if "ndb" in text or "no durable" in text:
        return "NDB"
    return None


def infer_timepoint(row: pd.Series) -> str:
    visit = str(row.get("visit_pre_or_on_treatment", "")).strip().lower()
    if visit in {"pre", "baseline", "pre-treatment", "pretreatment"}:
        return "pretreatment"
    if visit in {"on", "on-treatment", "on_treatment"}:
        return "on_treatment"
    title = str(row.get("title", "")).lower()
    if re.search(r"(^|[_\s.-])pre([_\s.-]|$)", title):
        return "pretreatment"
    if re.search(r"(^|[_\s.-])on([_\s.-]|$)", title):
        return "on_treatment"
    text = " | ".join(str(row.get(key, "")) for key in ["title", "source_name_ch1", "characteristics_ch1", "visit_pre_or_on_treatment", "treatment_state", "timepoint", "treatment_group", "biopsy", "description"]).lower()
    if any(token in text for token in ["progression", "acquired"]):
        return "progression"
    if any(token in text for token in ["on treatment", "on_treatment", "ontx", "_on_", "post.treatment", "post treatment", "on immune checkpoint"]):
        return "on_treatment"
    if any(token in text for token in ["pre", "baseline", "treatment.naive", "treatment naive"]):
        return "pretreatment"
    return "unknown"


def harmonized_geo_metadata(raw: pd.DataFrame, accession: str, platform: str | None = None) -> pd.DataFrame:
    out = raw.copy()
    out["sample_id"] = out.get("geo_accession", out.get("title", pd.Series(range(len(out))).astype(str)))
    out["cohort"] = accession
    out["accession"] = accession
    out["platform"] = platform or out.get("platform_id", pd.Series(["unknown"] * len(out))).astype(str)
    out["title"] = out.get("title", pd.Series([""] * len(out))).astype(str)
    out["source_name"] = out.get("source_name_ch1", pd.Series([""] * len(out))).astype(str)
    out["patient_id_raw"] = out.apply(infer_patient_id, axis=1)
    out["patient_id"] = out["patient_id_raw"]
    out["response_raw"] = out.apply(infer_response, axis=1)
    if accession == "GSE136961":
        title = out["title"].astype(str).str.upper()
        out.loc[title.str.startswith("D"), "response_raw"] = "DCB"
        out.loc[title.str.startswith("N"), "response_raw"] = "NDB"
    if accession == "GSE93157" and "best_resp" in out.columns:
        out["response_raw"] = out["best_resp"]
    if accession == "GSE140901" and "best_response" in out.columns:
        out["response_raw"] = out["best_response"]
    if accession == "GSE67501" and "description" in out.columns:
        out["patient_id"] = out["description"].where(out["description"].notna(), out["patient_id"])
        out["patient_id_raw"] = out["patient_id"]
    out["timepoint"] = out.apply(infer_timepoint, axis=1)
    if accession in {"GSE67501", "GSE93157", "GSE140901", "GSE136961"}:
        out["timepoint"] = "pretreatment"
    return out
