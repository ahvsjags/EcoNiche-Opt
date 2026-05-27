from __future__ import annotations

import argparse
import gzip
import io
import json
import re
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
GSE122220_MATRIX_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE122nnn/GSE122220/matrix/GSE122220_series_matrix.txt.gz"
GPL10558_ANNOT_URL = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz"


def _strip(value: object) -> str:
    return str(value).strip().strip('"')


def _fields(line: str) -> list[str]:
    return [_strip(part) for part in line.split("\t")]


def _first_metadata(lines: list[str], key: str) -> list[str]:
    prefix = f"!{key}\t"
    for line in lines:
        if line.startswith(prefix):
            return _fields(line)[1:]
    return []


def _metadata_rows(lines: list[str], key: str) -> list[list[str]]:
    prefix = f"!{key}\t"
    return [_fields(line)[1:] for line in lines if line.startswith(prefix)]


def _parse_characteristics(rows: list[list[str]], n_samples: int) -> list[dict[str, str]]:
    parsed = [dict() for _ in range(n_samples)]
    for row in rows:
        for idx, value in enumerate(row[:n_samples]):
            if ":" not in value:
                continue
            field, observed = value.split(":", 1)
            key = re.sub(r"[^a-z0-9]+", "_", field.strip().lower()).strip("_")
            parsed[idx][key] = observed.strip()
    return parsed


def _characteristic_value(parsed: dict[str, str], prefix: str) -> str:
    prefix = prefix.lower()
    for key, value in parsed.items():
        if key.startswith(prefix):
            return value
    return ""


def _normalize_gene_symbol(symbol: object) -> str:
    text = _strip(symbol)
    if not text or text.lower() == "nan":
        return ""
    text = re.split(r"\s*///\s*|;", text)[0].strip()
    text = text.split(",")[0].strip()
    return text.upper()


def _normalize_response(response: object) -> str:
    token = _strip(response).upper()
    if token in {"CR", "PR", "SD", "PD"}:
        return token
    return token


def _response_label(response_raw: str) -> float:
    token = _normalize_response(response_raw)
    if token in {"CR", "PR"}:
        return 1.0
    if token in {"SD", "PD"}:
        return 0.0
    return np.nan


def parse_gpl10558_annotation(text: str) -> pd.DataFrame:
    lines = text.splitlines()
    start = lines.index("!platform_table_begin") + 1
    end = lines.index("!platform_table_end") if "!platform_table_end" in lines else len(lines)
    table = pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t", dtype=str)
    if "ID" not in table.columns or "Gene symbol" not in table.columns:
        raise ValueError("GPL10558 annotation must contain ID and Gene symbol columns")
    annot = table[["ID", "Gene symbol"]].rename(columns={"ID": "probe_id", "Gene symbol": "gene_symbol"})
    annot["probe_id"] = annot["probe_id"].map(_strip)
    annot["gene_symbol"] = annot["gene_symbol"].map(_normalize_gene_symbol)
    annot = annot[annot["probe_id"].ne("") & annot["gene_symbol"].ne("")]
    return annot.drop_duplicates("probe_id").reset_index(drop=True)


def parse_gse122220_series_matrix(text: str, annotation: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    lines = text.splitlines()
    start = lines.index("!series_matrix_table_begin") + 1
    end = lines.index("!series_matrix_table_end") if "!series_matrix_table_end" in lines else len(lines)
    raw_matrix = pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t")
    raw_matrix = raw_matrix.rename(columns={raw_matrix.columns[0]: "probe_id"})
    raw_matrix["probe_id"] = raw_matrix["probe_id"].map(_strip)
    sample_columns = [col for col in raw_matrix.columns if col != "probe_id"]

    mapped = raw_matrix.merge(annotation, on="probe_id", how="inner")
    numeric = mapped[sample_columns].apply(pd.to_numeric, errors="coerce")
    if float(np.nanmedian(numeric.to_numpy())) > 50:
        numeric = np.log2(numeric + 1.0)
    numeric["gene_symbol"] = mapped["gene_symbol"].to_numpy()
    expression = numeric.groupby("gene_symbol", sort=True)[sample_columns].mean().T
    expression.index = [_strip(idx) for idx in expression.index]

    titles = _first_metadata(lines, "Sample_title")
    accessions = _first_metadata(lines, "Sample_geo_accession")
    source_names = _first_metadata(lines, "Sample_source_name_ch1")
    descriptions = _metadata_rows(lines, "Sample_description")
    characteristics = _parse_characteristics(_metadata_rows(lines, "Sample_characteristics_ch1"), len(accessions))

    rows: list[dict[str, object]] = []
    for idx, sample_id in enumerate(accessions):
        title = titles[idx] if idx < len(titles) else sample_id
        parsed = characteristics[idx] if idx < len(characteristics) else {}
        treatment_raw = _characteristic_value(parsed, "treatment")
        response_raw = _normalize_response(_characteristic_value(parsed, "response"))
        therapy = "anti-CTLA-4 plus anti-PD-1 combination" if "IPI" in treatment_raw.upper() else "anti-PD-1 monotherapy"
        description = " | ".join(row[idx] for row in descriptions if idx < len(row))
        rows.append(
            {
                "sample_id": sample_id,
                "patient_id": title,
                "cohort": "GSE122220",
                "accession": "GSE122220",
                "title": title,
                "cancer_type": "melanoma",
                "source_name": source_names[idx] if idx < len(source_names) else "Melanoma tumor biopsy",
                "sample_source": "tumor_tissue",
                "specimen_type": "tumor_biopsy",
                "platform": "Illumina_HumanHT12_V4_expression_beadchip",
                "therapy": therapy,
                "treatment_raw": treatment_raw,
                "previous_ipilimumab": _characteristic_value(parsed, "previous_ipilumimab"),
                "timepoint": "pretreatment",
                "baseline_status": "pretreatment_before_checkpoint_inhibitor",
                "response_raw": response_raw,
                "label": _response_label(response_raw),
                "label_status": "GEO_sample_characteristics_CR_PR_SD_PD; low_n_array_platform_sensitivity_only",
                "age": _characteristic_value(parsed, "age"),
                "gender": _characteristic_value(parsed, "gender"),
                "description": description,
                "expression_sample_id": sample_id,
                "match_method": "GEO_series_matrix_sample_accession",
            }
        )
    metadata = pd.DataFrame(rows)
    expression = expression.reindex(metadata["sample_id"].astype(str))
    return expression, metadata


def _download_text(url: str, raw_out: Path) -> str:
    raw_out.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=90) as response:
        payload = response.read()
    raw_out.write_bytes(payload)
    return gzip.decompress(payload).decode("utf-8", errors="replace")


def _strict_label(response_raw: str) -> str:
    token = _normalize_response(response_raw)
    if token in {"CR", "PR"}:
        return "1"
    if token == "PD":
        return "0"
    return ""


def _clinical_benefit_label(response_raw: str) -> str:
    token = _normalize_response(response_raw)
    if token in {"CR", "PR", "SD"}:
        return "1"
    if token == "PD":
        return "0"
    return ""


def _write_validation_templates(metadata: pd.DataFrame, out_dir: Path) -> tuple[Path, Path]:
    manifest_rows = []
    clinical_rows = []
    for _, row in metadata.iterrows():
        primary_label = "" if pd.isna(row["label"]) else int(row["label"])
        manifest_rows.append(
            {
                "site_id": "GEO",
                "subject_id": row["patient_id"],
                "sample_id": row["sample_id"],
                "patient_id_at_site": row["patient_id"],
                "specimen_id": row["sample_id"],
                "collection_timepoint": "baseline_pre_treatment",
                "collection_date": "",
                "days_before_icb_start": "",
                "cancer_type": "melanoma",
                "disease_stage": "metastatic",
                "sample_source": "tumor_tissue",
                "specimen_type": "tumor_biopsy",
                "anatomic_site": "",
                "baseline_status": "pretreatment_before_checkpoint_inhibitor",
                "therapy": row["therapy"],
                "line_of_therapy": "",
                "icb_start_date": "",
                "assay_platform": "Illumina HumanHT-12 V4.0 expression beadchip",
                "rna_input_ng": "",
                "rna_dv200_or_rin": "",
                "tumor_content_percent": "",
                "necrosis_percent": "",
                "macrodissection_performed": "",
                "housekeeping_qc_pass": "",
                "panel_gene_coverage_percent": "",
                "qc_pass": "TRUE",
                "locked_validation_use_flag": "include_low_n_array_sensitivity",
                "exclusion_reason": "",
                "data_freeze_id": "GSE122220_public_20260527",
            }
        )
        clinical_rows.append(
            {
                "subject_id": row["patient_id"],
                "sample_id": row["sample_id"],
                "icb_regimen": row["therapy"],
                "icb_start_date": "",
                "baseline_scan_date": "",
                "response_raw": row["response_raw"],
                "best_overall_response_date": "",
                "progression_date": "",
                "last_follow_up_date": "",
                "dcbr_6mo": "",
                "recist_version": "GEO sample characteristics list CR/PR/SD/PD; RECIST version not stated",
                "primary_recist_label": primary_label,
                "strict_recist_label": _strict_label(row["response_raw"]),
                "clinical_benefit_label": _clinical_benefit_label(row["response_raw"]),
                "response_assessor": "",
                "label_source_document": "GSE122220 GEO series matrix sample characteristics",
                "source_page_or_record_id": row["sample_id"],
                "curation_notes": "Pretreatment melanoma tumor biopsy expression array; low-n platform sensitivity evidence only, not counted as strict independent bulk RNA-seq external validation.",
            }
        )
    manifest = out_dir / "GSE122220.assay_sample_manifest.tsv"
    clinical = out_dir / "GSE122220.clinical_annotation.tsv"
    pd.DataFrame(manifest_rows).to_csv(manifest, sep="\t", index=False)
    pd.DataFrame(clinical_rows).to_csv(clinical, sep="\t", index=False)
    return manifest, clinical


def _write_qc(expression: pd.DataFrame, metadata: pd.DataFrame, panel_genes: Path, out: Path) -> None:
    panel = pd.read_csv(panel_genes, sep="\t")
    available = set(expression.columns.astype(str))
    locked = set(panel["gene_symbol"].astype(str))
    overlap = sorted(available & locked)
    labels = pd.to_numeric(metadata["label"], errors="coerce")
    qc = pd.DataFrame(
        [
            {
                "cohort": "GSE122220",
                "n_samples": int(len(expression)),
                "n_responders": int((labels == 1).sum()),
                "n_nonresponders": int((labels == 0).sum()),
                "response_raw_counts": ";".join(
                    f"{key}:{int(value)}" for key, value in metadata["response_raw"].value_counts().sort_index().items()
                ),
                "n_genes": int(expression.shape[1]),
                "locked_panel_genes": int(len(locked)),
                "locked_panel_overlap": int(len(overlap)),
                "locked_panel_overlap_fraction": len(overlap) / len(locked) if locked else np.nan,
                "overlap_genes": ",".join(overlap),
                "primary_use": "low_n_array_platform_sensitivity_not_strict_bulk_external",
            }
        ]
    )
    qc.to_csv(out, sep="\t", index=False)


def preprocess_gse122220(
    out_dir: Path,
    raw_dir: Path,
    panel_genes: Path,
    matrix_url: str = GSE122220_MATRIX_URL,
    annot_url: str = GPL10558_ANNOT_URL,
) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)
    matrix_text = _download_text(matrix_url, raw_dir / "GSE122220_series_matrix.txt.gz")
    annotation_text = _download_text(annot_url, raw_dir / "GPL10558.annot.gz")
    annotation = parse_gpl10558_annotation(annotation_text)
    expression, metadata = parse_gse122220_series_matrix(matrix_text, annotation)

    expr_path = out_dir / "GSE122220.expr.tsv"
    meta_path = out_dir / "GSE122220.metadata.tsv"
    qc_path = ROOT / "tables" / "gse122220_array_qc.tsv"
    expression.to_csv(expr_path, sep="\t", index_label="sample_id")
    metadata.to_csv(meta_path, sep="\t", index=False)
    manifest, clinical = _write_validation_templates(metadata, out_dir)
    qc_path.parent.mkdir(parents=True, exist_ok=True)
    _write_qc(expression, metadata, panel_genes, qc_path)
    return {
        "expression": expr_path,
        "metadata": meta_path,
        "sample_manifest": manifest,
        "clinical_annotation": clinical,
        "qc": qc_path,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", default="data/processed/bulk")
    parser.add_argument("--raw-dir", default="data/raw/GSE122220")
    parser.add_argument("--panel-genes", default="deliverables/prospective_validation/locked_panel_genes.tsv")
    parser.add_argument("--matrix-url", default=GSE122220_MATRIX_URL)
    parser.add_argument("--annot-url", default=GPL10558_ANNOT_URL)
    args = parser.parse_args()

    outputs = preprocess_gse122220(
        out_dir=ROOT / args.out_dir,
        raw_dir=ROOT / args.raw_dir,
        panel_genes=ROOT / args.panel_genes,
        matrix_url=args.matrix_url,
        annot_url=args.annot_url,
    )
    print(json.dumps({key: str(path) for key, path in outputs.items()}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
