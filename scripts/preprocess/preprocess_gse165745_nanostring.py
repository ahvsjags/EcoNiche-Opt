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
GSE165745_MATRIX_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165745/matrix/GSE165745_series_matrix.txt.gz"


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
            parsed[idx][field.strip().lower().replace(" ", "_")] = observed.strip()
    return parsed


def _normalize_gene_symbol(symbol: object) -> str:
    text = _strip(symbol)
    aliases = {
        "CD8a": "CD8A",
        "Cd8a": "CD8A",
        "Arg1": "ARG1",
        "Tgfb1": "TGFB1",
    }
    return aliases.get(text, text.upper())


def parse_gse165745_series_matrix(text: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    lines = text.splitlines()
    start = lines.index("!series_matrix_table_begin") + 1
    end = lines.index("!series_matrix_table_end") if "!series_matrix_table_end" in lines else len(lines)
    matrix_text = "\n".join(lines[start:end])
    raw_matrix = pd.read_csv(io.StringIO(matrix_text), sep="\t")
    raw_matrix = raw_matrix.rename(columns={raw_matrix.columns[0]: "gene_symbol"})
    raw_matrix["gene_symbol"] = raw_matrix["gene_symbol"].map(_normalize_gene_symbol)
    sample_columns = [col for col in raw_matrix.columns if col != "gene_symbol"]
    expression = raw_matrix.set_index("gene_symbol")[sample_columns].apply(pd.to_numeric, errors="coerce")
    expression = np.log2(expression + 1.0).T
    expression.index = [_strip(idx) for idx in expression.index]
    expression = expression.T.groupby(level=0).mean().T

    titles = _first_metadata(lines, "Sample_title")
    accessions = _first_metadata(lines, "Sample_geo_accession")
    descriptions = _metadata_rows(lines, "Sample_description")
    data_processing = _first_metadata(lines, "Sample_data_processing")
    source_names = _first_metadata(lines, "Sample_source_name_ch1")
    characteristics = _parse_characteristics(_metadata_rows(lines, "Sample_characteristics_ch1"), len(accessions))

    rows: list[dict[str, object]] = []
    for idx, sample_id in enumerate(accessions):
        title = titles[idx] if idx < len(titles) else sample_id
        patient_id = title.split(":", 1)[0].strip()
        parsed = characteristics[idx] if idx < len(characteristics) else {}
        phenotype = parsed.get("phenotype", "")
        response_label = 1 if phenotype.lower() == "responder" else 0 if phenotype.lower() == "nonresponder" else np.nan
        description = " | ".join(row[idx] for row in descriptions if idx < len(row))
        rows.append(
            {
                "sample_id": sample_id,
                "patient_id": patient_id,
                "cohort": "GSE165745",
                "accession": "GSE165745",
                "title": title,
                "cancer_type": "melanoma",
                "source_name": source_names[idx] if idx < len(source_names) else "Melanoma",
                "sample_source": "tumor_tissue",
                "specimen_type": "FFPE",
                "platform": "NanoString_nCounter_Vantage_3D_Wnt_panel",
                "therapy": "anti-PD-1 monotherapy",
                "timepoint": "pretreatment",
                "baseline_status": "pretreatment_before_first_ICB_dose",
                "response_raw": phenotype,
                "label": response_label,
                "label_status": "source_binary_responder_nonresponder_not_recist",
                "age": parsed.get("age", ""),
                "gender": parsed.get("gender", ""),
                "disease_site": parsed.get("disease_site", ""),
                "cell_type": parsed.get("cell_type", ""),
                "description": description,
                "data_processing": data_processing[idx] if idx < len(data_processing) else "",
                "expression_sample_id": sample_id,
                "match_method": "GEO_series_matrix_sample_accession",
            }
        )
    metadata = pd.DataFrame(rows)
    expression = expression.reindex(metadata["sample_id"].astype(str))
    return expression, metadata


def _download_text(url: str, raw_out: Path) -> str:
    raw_out.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=60) as response:
        payload = response.read()
    raw_out.write_bytes(payload)
    return gzip.decompress(payload).decode("utf-8", errors="replace")


def _write_validation_templates(metadata: pd.DataFrame, out_dir: Path) -> tuple[Path, Path]:
    manifest_rows = []
    clinical_rows = []
    for _, row in metadata.iterrows():
        label = "" if pd.isna(row["label"]) else int(row["label"])
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
                "specimen_type": "FFPE",
                "anatomic_site": row["disease_site"],
                "baseline_status": "pretreatment_before_first_ICB_dose",
                "therapy": "anti-PD-1 monotherapy",
                "line_of_therapy": "",
                "icb_start_date": "",
                "assay_platform": "NanoString nCounter Vantage 3D Wnt Pathways Panel",
                "rna_input_ng": "",
                "rna_dv200_or_rin": "",
                "tumor_content_percent": "",
                "necrosis_percent": "",
                "macrodissection_performed": "",
                "housekeeping_qc_pass": "",
                "panel_gene_coverage_percent": "",
                "qc_pass": "TRUE",
                "locked_validation_use_flag": "include",
                "exclusion_reason": "",
                "data_freeze_id": "GSE165745_public_20260527",
            }
        )
        clinical_rows.append(
            {
                "subject_id": row["patient_id"],
                "sample_id": row["sample_id"],
                "icb_regimen": "anti-PD-1 monotherapy",
                "icb_start_date": "",
                "baseline_scan_date": "",
                "response_raw": row["response_raw"],
                "best_overall_response_date": "",
                "progression_date": "",
                "last_follow_up_date": "",
                "dcbr_6mo": "",
                "recist_version": "source binary responder/nonresponder; RECIST category not available in GEO record",
                "primary_recist_label": label,
                "strict_recist_label": "",
                "clinical_benefit_label": label,
                "response_assessor": "",
                "label_source_document": "GSE165745 GEO sample phenotype",
                "source_page_or_record_id": row["sample_id"],
                "curation_notes": "Responder/nonresponder label is used only for panel-transfer sensitivity; not counted as strict RECIST melanoma bulk RNA-seq external validation.",
            }
        )
    manifest = out_dir / "GSE165745.assay_sample_manifest.tsv"
    clinical = out_dir / "GSE165745.clinical_annotation.tsv"
    pd.DataFrame(manifest_rows).to_csv(manifest, sep="\t", index=False)
    pd.DataFrame(clinical_rows).to_csv(clinical, sep="\t", index=False)
    return manifest, clinical


def _write_qc(expression: pd.DataFrame, metadata: pd.DataFrame, panel_genes: Path, out: Path) -> None:
    panel = pd.read_csv(panel_genes, sep="\t")
    available = set(expression.columns.astype(str))
    locked = set(panel["gene_symbol"].astype(str))
    overlap = sorted(available & locked)
    qc = pd.DataFrame(
        [
            {
                "cohort": "GSE165745",
                "n_samples": int(len(expression)),
                "n_responders": int(pd.to_numeric(metadata["label"], errors="coerce").fillna(0).sum()),
                "n_nonresponders": int((pd.to_numeric(metadata["label"], errors="coerce") == 0).sum()),
                "n_genes": int(expression.shape[1]),
                "locked_panel_genes": int(len(locked)),
                "locked_panel_overlap": int(len(overlap)),
                "locked_panel_overlap_fraction": len(overlap) / len(locked) if locked else np.nan,
                "overlap_genes": ",".join(overlap),
                "primary_use": "panel_transfer_sensitivity_not_strict_bulk_external",
            }
        ]
    )
    qc.to_csv(out, sep="\t", index=False)


def preprocess_gse165745(
    out_dir: Path,
    raw_dir: Path,
    panel_genes: Path,
    matrix_url: str = GSE165745_MATRIX_URL,
) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)
    text = _download_text(matrix_url, raw_dir / "GSE165745_series_matrix.txt.gz")
    expression, metadata = parse_gse165745_series_matrix(text)
    expr_path = out_dir / "GSE165745.expr.tsv"
    meta_path = out_dir / "GSE165745.metadata.tsv"
    qc_path = ROOT / "tables" / "gse165745_panel_qc.tsv"
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
    parser.add_argument("--raw-dir", default="data/raw/GSE165745")
    parser.add_argument("--panel-genes", default="deliverables/prospective_validation/locked_panel_genes.tsv")
    parser.add_argument("--matrix-url", default=GSE165745_MATRIX_URL)
    args = parser.parse_args()

    outputs = preprocess_gse165745(
        out_dir=ROOT / args.out_dir,
        raw_dir=ROOT / args.raw_dir,
        panel_genes=ROOT / args.panel_genes,
        matrix_url=args.matrix_url,
    )
    print(json.dumps({key: str(path) for key, path in outputs.items()}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
