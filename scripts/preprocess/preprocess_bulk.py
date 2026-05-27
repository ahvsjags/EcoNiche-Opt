from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.expression import (
    choose_expression_file,
    clean_expression_matrix,
    download_geo_platform_soft,
    download_ncbi_gene_info,
    load_geo_probe_symbol_map,
    load_entrez_symbol_map,
    load_or_query_mygene_ensembl_symbol_map,
    read_table_matrix,
)
from econiche.normalize import log2_if_needed
from econiche.registry import load_registry, normalize_access_status


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", default="data/raw")
    parser.add_argument("--registry", default="config/data_registry.yml")
    parser.add_argument("--metadata", default="data/metadata/metadata_harmonized.tsv")
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--report", default="tables/expression_qc_report.tsv")
    args = parser.parse_args()
    processed = ROOT / args.processed_dir
    processed.mkdir(parents=True, exist_ok=True)
    registry = load_registry(ROOT / args.registry)
    metadata_path = ROOT / args.metadata
    metadata_all = pd.read_csv(metadata_path, sep="\t") if metadata_path.exists() else pd.DataFrame()
    gene_info = download_ncbi_gene_info(ROOT / "data/priors/Homo_sapiens.gene_info.gz")
    entrez_map = load_entrez_symbol_map(gene_info)
    rows = []
    for cohort in registry.get("cohorts", []):
        accession = cohort["accession"]
        if not str(accession).startswith("GSE") or normalize_access_status(cohort.get("access")) != "public":
            continue
        if "scRNA" in str(cohort.get("platform", "")):
            continue
        expr_path = choose_expression_file(accession, ROOT / args.raw_dir)
        if expr_path is None:
            rows.append({"cohort": accession, "n_samples": 0, "n_genes": 0, "status": "missing_expression_file", "source_file": ""})
            continue
        try:
            if expr_path.suffix == ".xlsx":
                df = read_table_matrix(expr_path)
            elif expr_path.name.endswith(".csv.gz") or expr_path.suffix == ".csv":
                df = read_table_matrix(expr_path, sep=",")
            else:
                df = read_table_matrix(expr_path, sep="\t")
            gene_map = dict(entrez_map)
            if any(str(value).startswith("ENSG") for value in df.index[: min(200, len(df.index))]):
                cache = ROOT / f"data/priors/{accession}_ensembl_symbol_map.tsv"
                gene_map.update(load_or_query_mygene_ensembl_symbol_map(df.index, cache))
            if accession == "GSE67501":
                platform_path = download_geo_platform_soft("GPL14951", ROOT / "data/external/geo_platforms")
                gene_map.update(load_geo_probe_symbol_map(platform_path))
            X = clean_expression_matrix(df, entrez_map=gene_map)
            if X.empty or X.shape[1] < 10:
                rows.append({"cohort": accession, "n_samples": X.shape[0], "n_genes": X.shape[1], "status": "unusable_expression_shape", "source_file": str(expr_path)})
                continue
            X = log2_if_needed(X)
            if accession == "GSE145996":
                meta = load_gse145996_supplement(ROOT, X.index)
            else:
                meta = metadata_all[metadata_all.get("accession", pd.Series(dtype=str)).astype(str) == accession].copy()
            if meta.empty:
                meta = pd.DataFrame({"sample_id": X.index, "patient_id": X.index, "cohort": accession, "accession": accession, "label": pd.NA})
            meta = align_metadata_to_expression(meta, X.index, accession)
            if "label" in meta.columns:
                keep = meta["label"].notna()
                X = X.loc[meta.loc[keep, "sample_id"]]
                meta = meta.loc[keep].copy()
            if "timepoint" in meta.columns and meta["timepoint"].notna().any():
                primary_mask = meta["timepoint"].astype(str).str.lower().isin(["pretreatment", "baseline", "pre_treatment"])
                if primary_mask.any():
                    X = X.loc[meta.loc[primary_mask, "sample_id"]]
                    meta = meta.loc[primary_mask].copy()
            if "patient_id" in meta.columns and meta["patient_id"].notna().any():
                meta = meta.drop_duplicates("patient_id", keep="first").copy()
                X = X.loc[meta["sample_id"]]
            if len(meta) < 4 or meta.get("label", pd.Series(dtype=float)).nunique(dropna=True) < 2:
                status = "metadata_needs_manual_curation"
            else:
                status = "processed"
            X.to_csv(processed / f"{accession}.expr.tsv", sep="\t")
            meta.to_csv(processed / f"{accession}.metadata.tsv", sep="\t", index=False)
            rows.append({"cohort": accession, "n_samples": X.shape[0], "n_genes": X.shape[1], "status": status, "source_file": str(expr_path)})
        except Exception as exc:
            rows.append({"cohort": accession, "n_samples": 0, "n_genes": 0, "status": "preprocess_failed", "source_file": str(expr_path), "reason": str(exc)})
    report = pd.DataFrame(rows)
    out = ROOT / args.report
    out.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(out, sep="\t", index=False)
    print(f"Wrote expression QC report to {out}")


def normalize_sample_token(value: object) -> str:
    text = str(value)
    text = re.sub(r"\.bam$", "", text)
    text = text.replace(".baseline", "_Pre").replace(".OnTx", "_On")
    text = text.replace("-", "_").replace(".", "_")
    return re.sub(r"[^A-Za-z0-9]+", "", text).lower()


def _response_to_nonresponse_label(response: object) -> object:
    text = str(response).strip().upper()
    if text in {"CR", "PR", "R", "DCB"}:
        return 0
    if text in {"SD", "PD", "NR", "NDB"}:
        return 1
    return pd.NA


def load_gse145996_supplement(root: Path, expression_index) -> pd.DataFrame:
    supplement_root = root / "data/external/GSE145996"
    xlsx_files = sorted(supplement_root.rglob("*.xlsx"))
    if not xlsx_files:
        return pd.DataFrame()
    xlsx = xlsx_files[0]
    clinical = pd.read_excel(xlsx, sheet_name="SuppTable1", header=1)
    specimens = pd.read_excel(xlsx, sheet_name="SuppTable2", header=1)
    clinical = clinical.rename(columns=lambda value: str(value).strip())
    specimens = specimens.rename(columns=lambda value: str(value).strip())
    clinical = clinical.dropna(subset=["Patient Number"]).copy()
    specimens = specimens.dropna(subset=["Patient Number", "Tumor ID number"]).copy()
    specimens["RNASEQ"] = pd.to_numeric(specimens.get("RNASEQ"), errors="coerce").fillna(0).astype(int)
    specimens = specimens[specimens["RNASEQ"] == 1].copy()
    merged = specimens.merge(
        clinical[["Patient Number", "Drug", "Best Response"]],
        on="Patient Number",
        how="left",
        validate="many_to_one",
    )
    merged["tumor_token"] = merged["Tumor ID number"].astype(str).map(normalize_sample_token)
    by_tumor = {row["tumor_token"]: row for _, row in merged.iterrows()}

    rows = []
    for sample in expression_index:
        sample_text = str(sample)
        match = re.search(r"(MB[_-]?\d+)", sample_text, flags=re.I)
        tumor_token = normalize_sample_token(match.group(1)) if match else normalize_sample_token(sample_text)
        row = by_tumor.get(tumor_token)
        if row is None:
            rows.append(
                {
                    "sample_id": sample,
                    "patient_id": sample,
                    "cohort": "GSE145996",
                    "accession": "GSE145996",
                    "label": pd.NA,
                    "match_method": "supplement_tumor_id_unmatched",
                    "evidence_source_file": str(xlsx),
                }
            )
            continue
        response = row.get("Best Response", pd.NA)
        rows.append(
            {
                "sample_id": sample,
                "expression_sample_id": sample,
                "patient_id": row.get("Patient Number", pd.NA),
                "patient_id_raw": row.get("Patient Number", pd.NA),
                "tumor_id": row.get("Tumor ID number", pd.NA),
                "geo_accession": pd.NA,
                "cohort": "GSE145996",
                "accession": "GSE145996",
                "response_raw": response,
                "label": _response_to_nonresponse_label(response),
                "timepoint": "pretreatment",
                "therapy": row.get("Drug", "anti-PD1"),
                "treatment": row.get("Drug", "anti-PD1"),
                "match_method": "supplement_tumor_id",
                "evidence_source_file": str(xlsx),
                "evidence_field": "SuppTable1 Patient Number/Best Response plus SuppTable2 Tumor ID number/RNASEQ",
                "characteristics_ch1": f"patient={row.get('Patient Number', '')}|tumor_id={row.get('Tumor ID number', '')}|best_response={response}|rnaseq=1",
            }
        )
    return pd.DataFrame(rows)


def align_metadata_to_expression(meta: pd.DataFrame, expression_index, accession: str) -> pd.DataFrame:
    meta = meta.copy()
    if "sample_id" not in meta.columns:
        meta["sample_id"] = meta.get("geo_accession", meta.get("title", pd.Series(dtype=str)))
    if accession in {"GSE165252", "GSE93157"} and len(meta) == len(expression_index):
        ordered = meta.reset_index(drop=True).copy()
        ordered["sample_id"] = list(expression_index)
        ordered["expression_sample_id"] = ordered["sample_id"]
        ordered["cohort"] = accession
        ordered["accession"] = accession
        ordered["match_method"] = "geo_order"
        if "patient_id" not in ordered.columns or ordered["patient_id"].isna().all():
            ordered["patient_id"] = ordered.get("patient_id_raw", ordered["sample_id"])
        return ordered
    candidates = []
    for _, row in meta.iterrows():
        values = [
            row.get("sample_id"),
            row.get("geo_accession"),
            row.get("title"),
            row.get("description"),
            row.get("patient_id"),
            row.get("patient_id_raw"),
        ]
        for value in values:
            if pd.notna(value):
                candidates.append((normalize_sample_token(value), row))
    matched = []
    used = set()
    candidate_map = {}
    for key, row in candidates:
        candidate_map.setdefault(key, row)
    for sample in expression_index:
        norm = normalize_sample_token(sample)
        row = None
        if norm in candidate_map:
            row = candidate_map[norm].copy()
        else:
            for key, candidate in sorted(candidate_map.items(), key=lambda item: len(item[0]), reverse=True):
                if key and (key in norm or norm in key):
                    row = candidate.copy()
                    break
        if row is None:
            row = pd.Series({"sample_id": sample, "patient_id": sample, "cohort": accession, "accession": accession, "label": pd.NA})
        row["sample_id"] = sample
        row["expression_sample_id"] = sample
        row["cohort"] = accession
        row["accession"] = accession
        if pd.isna(row.get("match_method", pd.NA)) or not str(row.get("match_method", "")).strip():
            row["match_method"] = "token"
        if pd.isna(row.get("patient_id", pd.NA)):
            row["patient_id"] = row.get("patient_id_raw", sample)
        matched.append(row)
        used.add(sample)
    return pd.DataFrame(matched).reset_index(drop=True)


if __name__ == "__main__":
    main()
