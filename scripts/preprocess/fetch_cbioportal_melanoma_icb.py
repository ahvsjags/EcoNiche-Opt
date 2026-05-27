from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.baselines import BASELINE_SIGNATURES
from econiche_opt.model.endpoint_modules import MODULE_GENE_SETS, WORD_INTERACTION_EDGES, WORD_STATE_GENE_SETS


CBIOPORTAL_BASE_URL = "https://www.cbioportal.org/api"
LOCKED_RESCUE_HEAD_GENES = ("MAP4K1", "TBX3", "AXL", "PLA2G2D", "PIK3CD")


@dataclass(frozen=True)
class CbioStudySpec:
    study_id: str
    molecular_profile_id: str
    sample_list_id: str
    output_cohort: str
    response_attribute: str
    treatment_attribute: str
    sample_treatment_attribute: str | None = None
    sample_treatment_keep: tuple[str, ...] = ("pre", "pretreatment", "baseline")
    response_source_level: str = "sample"
    treatment_source_level: str = "sample"
    notes: str = ""


STUDY_SPECS: dict[str, CbioStudySpec] = {
    "mel_dfci_2019": CbioStudySpec(
        study_id="mel_dfci_2019",
        molecular_profile_id="mel_dfci_2019_mrna_seq_tpm",
        sample_list_id="mel_dfci_2019_tpm",
        output_cohort="CBIO_LIU_DFCI_2019_PRE",
        response_attribute="BR",
        treatment_attribute="IO_THERAPY",
        response_source_level="patient",
        treatment_source_level="patient",
        notes="Original Liu/DFCI Nat Med 2019 anti-PD1 melanoma cBioPortal study; all samples are pretreatment by study definition.",
    ),
    "mel_iatlas_liu_2019": CbioStudySpec(
        study_id="mel_iatlas_liu_2019",
        molecular_profile_id="mel_iatlas_liu_2019_rna_seq_mrna",
        sample_list_id="mel_iatlas_liu_2019_rna_seq_mrna",
        output_cohort="CBIO_IATLAS_LIU_2019_PRE",
        response_attribute="RESPONSE",
        treatment_attribute="ICI_RX",
        sample_treatment_attribute="SAMPLE_TREATMENT",
        response_source_level="sample",
        treatment_source_level="patient",
        notes="iAtlas harmonized Liu cohort; duplicate-source cross-check rather than an independent Liu validation cohort.",
    ),
    "mel_iatlas_gide_2019": CbioStudySpec(
        study_id="mel_iatlas_gide_2019",
        molecular_profile_id="mel_iatlas_gide_2019_rna_seq_mrna",
        sample_list_id="mel_iatlas_gide_2019_rna_seq_mrna",
        output_cohort="CBIO_IATLAS_GIDE_2019_PRE",
        response_attribute="RESPONSE",
        treatment_attribute="ICI_RX",
        sample_treatment_attribute="SAMPLE_TREATMENT",
        response_source_level="sample",
        treatment_source_level="patient",
        notes="iAtlas harmonized Gide cohort; overlaps PRJEB23709 and should be used for source cross-check unless explicitly split from discovery.",
    ),
    "mel_iatlas_riaz_nivolumab_2017": CbioStudySpec(
        study_id="mel_iatlas_riaz_nivolumab_2017",
        molecular_profile_id="mel_iatlas_riaz_nivolumab_2017_rna_seq_mrna",
        sample_list_id="mel_iatlas_riaz_nivolumab_2017_rna_seq_mrna",
        output_cohort="CBIO_IATLAS_RIAZ_2017_PRE",
        response_attribute="RESPONSE",
        treatment_attribute="ICI_RX",
        sample_treatment_attribute="SAMPLE_TREATMENT",
        response_source_level="sample",
        treatment_source_level="patient",
        notes="iAtlas harmonized Riaz cohort; overlaps GSE91061 discovery data and is not independent external evidence.",
    ),
    "mel_iatlas_hugo_ucla_2016": CbioStudySpec(
        study_id="mel_iatlas_hugo_ucla_2016",
        molecular_profile_id="mel_iatlas_hugo_ucla_2016_rna_seq_mrna",
        sample_list_id="mel_iatlas_hugo_ucla_2016_rna_seq_mrna",
        output_cohort="CBIO_IATLAS_HUGO_2016_PRE",
        response_attribute="RESPONSE",
        treatment_attribute="ICI_RX",
        sample_treatment_attribute="SAMPLE_TREATMENT",
        response_source_level="sample",
        treatment_source_level="patient",
        notes="iAtlas harmonized Hugo/UCLA cohort; overlaps GSE78220 discovery data and is not independent external evidence.",
    ),
}


def requested_gene_symbols(extra_gene_files: Iterable[Path] = ()) -> list[str]:
    genes: set[str] = set()
    genes.update(LOCKED_RESCUE_HEAD_GENES)
    for gene_set in MODULE_GENE_SETS.values():
        genes.update(gene_set)
    for gene_set in WORD_STATE_GENE_SETS.values():
        genes.update(gene_set)
    for _, _, gene_a, gene_b, _ in WORD_INTERACTION_EDGES:
        genes.add(gene_a)
        genes.add(gene_b)
    for gene_set in BASELINE_SIGNATURES.values():
        genes.update(gene_set)
    for path in extra_gene_files:
        if not path.exists():
            continue
        frame = pd.read_csv(path, sep="\t")
        if "gene" in frame.columns:
            genes.update(frame["gene"].dropna().astype(str))
        elif "gene_symbol" in frame.columns:
            genes.update(frame["gene_symbol"].dropna().astype(str))
        elif frame.shape[1] >= 1:
            genes.update(frame.iloc[:, 0].dropna().astype(str))
    return sorted(gene.strip() for gene in genes if gene and str(gene).strip())


def _post_json(url: str, payload: object, params: dict[str, object] | None = None, timeout: int = 60) -> object:
    response = requests.post(url, params=params or {}, json=payload, timeout=timeout)
    response.raise_for_status()
    return response.json()


def _get_json(url: str, params: dict[str, object] | None = None, timeout: int = 60) -> object:
    response = requests.get(url, params=params or {}, timeout=timeout)
    response.raise_for_status()
    return response.json()


def resolve_genes(gene_symbols: list[str], base_url: str = CBIOPORTAL_BASE_URL) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for start in range(0, len(gene_symbols), 250):
        chunk = gene_symbols[start : start + 250]
        data = _post_json(
            f"{base_url}/genes/fetch",
            chunk,
            params={"geneIdType": "HUGO_GENE_SYMBOL", "projection": "SUMMARY"},
        )
        rows.extend(data)
    frame = pd.DataFrame(rows).drop_duplicates(subset=["hugoGeneSymbol", "entrezGeneId"])
    return frame.sort_values(["hugoGeneSymbol", "entrezGeneId"]).reset_index(drop=True)


def fetch_sample_ids(sample_list_id: str, base_url: str = CBIOPORTAL_BASE_URL) -> list[str]:
    data = _get_json(f"{base_url}/sample-lists/{sample_list_id}/sample-ids")
    return [str(sample_id) for sample_id in data]


def fetch_expression(
    molecular_profile_id: str,
    sample_list_id: str,
    genes: pd.DataFrame,
    base_url: str = CBIOPORTAL_BASE_URL,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    entrez_ids = genes["entrezGeneId"].dropna().astype(int).drop_duplicates().tolist()
    for start in range(0, len(entrez_ids), 120):
        chunk = entrez_ids[start : start + 120]
        payload = {"entrezGeneIds": chunk, "sampleListId": sample_list_id}
        data = _post_json(
            f"{base_url}/molecular-profiles/{molecular_profile_id}/molecular-data/fetch",
            payload,
            params={"projection": "DETAILED"},
            timeout=90,
        )
        rows.extend(data)
        time.sleep(0.05)
    raw = pd.DataFrame(rows)
    if raw.empty:
        return pd.DataFrame()
    raw["hugoGeneSymbol"] = raw["gene"].map(lambda item: item.get("hugoGeneSymbol") if isinstance(item, dict) else np.nan)
    raw["value"] = pd.to_numeric(raw["value"], errors="coerce")
    expr = raw.pivot_table(index="sampleId", columns="hugoGeneSymbol", values="value", aggfunc="first")
    expr = expr.reindex(columns=genes["hugoGeneSymbol"].drop_duplicates().tolist())
    expr.index.name = "sample_id"
    return expr.sort_index(axis=0).sort_index(axis=1)


def fetch_clinical_matrix(study_id: str, level: str, base_url: str = CBIOPORTAL_BASE_URL) -> pd.DataFrame:
    data = _get_json(
        f"{base_url}/studies/{study_id}/clinical-data",
        params={"clinicalDataType": level.upper(), "projection": "SUMMARY", "pageSize": 1000000},
        timeout=120,
    )
    frame = pd.DataFrame(data)
    if frame.empty:
        return pd.DataFrame()
    index_col = "sampleId" if level.upper() == "SAMPLE" else "patientId"
    matrix = frame.pivot_table(index=index_col, columns="clinicalAttributeId", values="value", aggfunc="first")
    matrix.index.name = "sample_id" if level.upper() == "SAMPLE" else "patient_id"
    return matrix.sort_index(axis=0).sort_index(axis=1)


def response_to_raw(value: object) -> str | None:
    if value is None or pd.isna(value):
        return None
    text = str(value).strip().lower().replace("_", " ").replace("-", " ")
    if not text:
        return None
    if text in {"cr", "complete response", "complete responder"}:
        return "CR"
    if text in {"pr", "partial response", "partial responder"}:
        return "PR"
    if text in {"sd", "stable disease"}:
        return "SD"
    if text in {"pd", "progressive disease", "progression"}:
        return "PD"
    if text in {"mixed response", "mr"}:
        return "MR"
    if text in {"responder", "response", "r"}:
        return "R"
    if text in {"non responder", "nonresponse", "non response", "non-responder", "nr"}:
        return "NR"
    if text in {"yes", "benefit", "clinical benefit", "durable clinical benefit"}:
        return "DCB"
    if text in {"no", "no benefit", "no clinical benefit", "nondurable clinical benefit"}:
        return "NDB"
    return str(value).strip()


def binary_label_from_response_raw(response_raw: object) -> float:
    token = response_to_raw(response_raw)
    if token in {"CR", "PR", "MR", "R", "DCB"}:
        return 1.0
    if token in {"SD", "PD", "NR", "NDB"}:
        return 0.0
    return np.nan


def _source_frame(sample_clinical: pd.DataFrame, patient_clinical: pd.DataFrame, source_level: str) -> pd.DataFrame:
    if source_level == "patient":
        return patient_clinical
    return sample_clinical


def build_metadata(
    spec: CbioStudySpec,
    expr: pd.DataFrame,
    sample_clinical: pd.DataFrame,
    patient_clinical: pd.DataFrame,
) -> pd.DataFrame:
    sample_meta = sample_clinical.reindex(expr.index).copy()
    sample_meta.insert(0, "sample_id", sample_meta.index.astype(str))
    if "PATIENT_ID" in sample_meta.columns:
        sample_meta["patient_id"] = sample_meta["PATIENT_ID"].astype(str)
    else:
        # cBioPortal sample IDs in these melanoma cohorts are one sample per patient
        # or use same numeric suffix; patient-level clinical joins below harden labels.
        sample_meta["patient_id"] = sample_meta["sample_id"].str.replace("Liu_Sample", "Patient", regex=False).str.replace("Sample", "Patient", regex=False)
    if spec.study_id == "mel_iatlas_liu_2019":
        sample_meta["patient_id"] = sample_meta["sample_id"].str.replace("Liu_Sample", "Patient", regex=False)
    elif spec.study_id.startswith("mel_iatlas"):
        # iAtlas sample identifiers such as PD02_Pre and Pt101_pre share patient
        # IDs with prefixes before the final timepoint suffix.
        sample_meta["patient_id"] = sample_meta["sample_id"].str.replace(r"_(Pre|pre|On|on|Prog|prog)$", "", regex=True)
    patient_meta = patient_clinical.reindex(sample_meta["patient_id"].astype(str)).copy()
    patient_meta.index = sample_meta.index

    response_frame = _source_frame(sample_meta, patient_meta, spec.response_source_level)
    treatment_frame = _source_frame(sample_meta, patient_meta, spec.treatment_source_level)
    response = response_frame.get(spec.response_attribute, pd.Series(index=expr.index, dtype=object)).reindex(expr.index)
    treatment = treatment_frame.get(spec.treatment_attribute, pd.Series(index=expr.index, dtype=object)).reindex(expr.index)
    response_raw = response.map(response_to_raw)
    label = response_raw.map(binary_label_from_response_raw)
    if spec.sample_treatment_attribute and spec.sample_treatment_attribute in sample_meta.columns:
        timepoint_raw = sample_meta[spec.sample_treatment_attribute].astype(str)
    else:
        timepoint_raw = pd.Series("Pre", index=expr.index)
    timepoint_norm = timepoint_raw.str.strip().str.lower()
    keep = timepoint_norm.isin({value.lower() for value in spec.sample_treatment_keep})

    metadata = pd.DataFrame(
        {
            "sample_id": sample_meta["sample_id"].astype(str),
            "patient_id": sample_meta["patient_id"].astype(str),
            "cohort": spec.output_cohort,
            "source_study_id": spec.study_id,
            "source_molecular_profile_id": spec.molecular_profile_id,
            "source_sample_list_id": spec.sample_list_id,
            "platform": "cBioPortal_RNA_seq_expression",
            "cancer_type": "melanoma",
            "treatment": treatment.fillna("").astype(str),
            "sample_treatment": timepoint_raw.astype(str),
            "timepoint": np.where(keep, "pretreatment", "non_baseline_or_unknown"),
            "response_source_attribute": spec.response_attribute,
            "response_raw_source": response.fillna("").astype(str),
            "response_raw": response_raw,
            "label": label,
            "endpoint_rule": "primary label uses CR/PR/MR/R/DCB=1 and SD/PD/NR/NDB=0; strict RECIST endpoint excludes MR and SD downstream",
            "source_notes": spec.notes,
        },
        index=expr.index,
    )
    metadata = metadata[keep & metadata["label"].notna()].copy()
    metadata["label"] = metadata["label"].astype(int)
    return metadata.reset_index(drop=True)


def fetch_and_write_study(spec: CbioStudySpec, genes: pd.DataFrame, out_dir: Path, raw_dir: Path) -> dict[str, object]:
    raw_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    sample_ids = fetch_sample_ids(spec.sample_list_id)
    expr = fetch_expression(spec.molecular_profile_id, spec.sample_list_id, genes)
    expr = expr.reindex(sample_ids).dropna(axis=0, how="all")
    sample_clinical = fetch_clinical_matrix(spec.study_id, "SAMPLE")
    patient_clinical = fetch_clinical_matrix(spec.study_id, "PATIENT")
    metadata = build_metadata(spec, expr, sample_clinical, patient_clinical)
    expr = expr.reindex(metadata["sample_id"]).dropna(axis=1, how="all")

    expr_path = out_dir / f"{spec.output_cohort}.expr.tsv"
    meta_path = out_dir / f"{spec.output_cohort}.metadata.tsv"
    expr.to_csv(expr_path, sep="\t")
    metadata.to_csv(meta_path, sep="\t", index=False)
    sample_clinical.to_csv(raw_dir / f"{spec.study_id}.sample_clinical.tsv", sep="\t")
    patient_clinical.to_csv(raw_dir / f"{spec.study_id}.patient_clinical.tsv", sep="\t")
    return {
        "study_id": spec.study_id,
        "cohort": spec.output_cohort,
        "molecular_profile_id": spec.molecular_profile_id,
        "sample_list_id": spec.sample_list_id,
        "n_samples_in_sample_list": len(sample_ids),
        "n_samples_written": int(expr.shape[0]),
        "n_genes_requested": int(len(genes)),
        "n_genes_written": int(expr.shape[1]),
        "n_responders": int(metadata["label"].sum()) if not metadata.empty else 0,
        "n_nonresponders": int((metadata["label"] == 0).sum()) if not metadata.empty else 0,
        "expr_path": str(expr_path.relative_to(ROOT)),
        "metadata_path": str(meta_path.relative_to(ROOT)),
        "notes": spec.notes,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch public cBioPortal melanoma ICB expression and labels.")
    parser.add_argument("--studies", nargs="*", default=list(STUDY_SPECS), choices=sorted(STUDY_SPECS))
    parser.add_argument("--out", default="data/processed/cbioportal_melanoma")
    parser.add_argument("--raw-out", default="data/raw/cbioportal_melanoma")
    parser.add_argument("--gene-file", action="append", default=[])
    parser.add_argument("--manifest", default="data/external/cbioportal_melanoma_manifest.tsv")
    args = parser.parse_args()

    out_dir = ROOT / args.out
    raw_dir = ROOT / args.raw_out
    manifest_path = ROOT / args.manifest
    gene_files = [ROOT / path for path in args.gene_file]
    symbols = requested_gene_symbols(gene_files)
    genes = resolve_genes(symbols)
    raw_dir.mkdir(parents=True, exist_ok=True)
    genes.to_csv(raw_dir / "resolved_gene_symbols.tsv", sep="\t", index=False)
    (raw_dir / "requested_gene_symbols.json").write_text(json.dumps(symbols, indent=2), encoding="utf-8")

    rows = []
    for study in args.studies:
        rows.append(fetch_and_write_study(STUDY_SPECS[study], genes, out_dir, raw_dir))
    manifest = pd.DataFrame(rows)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(manifest_path, sep="\t", index=False)
    print(f"Wrote {len(manifest)} cBioPortal melanoma cohorts to {out_dir}")
    print(f"Wrote manifest to {manifest_path}")


if __name__ == "__main__":
    main()
