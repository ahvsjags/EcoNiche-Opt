from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.expression import clean_expression_matrix
from econiche.normalize import log2_if_needed


def response_to_nonresponse_label(response: object, response_nr: object) -> object:
    token = str(response).strip().upper()
    nr = str(response_nr).strip().upper()
    if token in {"CR", "PR", "MR"} or nr == "R":
        return 0
    if token in {"SD", "PD"} or nr == "N":
        return 1
    return pd.NA


def load_tiger_expression(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    if "GENE_SYMBOL" not in raw.columns:
        raise ValueError(f"{path} does not contain a GENE_SYMBOL column")
    expression = clean_expression_matrix(raw.set_index("GENE_SYMBOL"))
    expression = log2_if_needed(expression)
    expression.index = expression.index.astype(str)
    return expression


def build_patient_like_cohort(expression: pd.DataFrame, clinical: pd.DataFrame, expression_path: Path, clinical_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    clinical = clinical.copy()
    clinical["sample_id"] = clinical["sample_id"].astype(str)
    clinical["patient_name"] = clinical["patient_name"].astype(str)
    cohort_meta = clinical[
        clinical["Therapy"].astype(str).str.upper().eq("ANTI-PD-1")
        & clinical["patient_name"].str.startswith("Patient")
    ].copy()
    if cohort_meta.empty:
        raise ValueError("No Patient-like anti-PD-1 melanoma samples found in TIGER phs000452 clinical table")
    missing = sorted(set(cohort_meta["sample_id"]) - set(expression.index))
    if missing:
        raise ValueError(f"Expression matrix is missing clinical samples: {missing[:5]}")
    cohort_meta = cohort_meta.sort_values(["patient_name", "sample_id"]).drop_duplicates("patient_name", keep="first")
    cohort_meta["label"] = [
        response_to_nonresponse_label(response, response_nr)
        for response, response_nr in zip(cohort_meta["response"], cohort_meta["response_NR"])
    ]
    cohort_meta = cohort_meta[cohort_meta["label"].notna()].copy()
    sample_ids = cohort_meta["sample_id"].astype(str).tolist()
    X = expression.loc[sample_ids].copy()
    meta = pd.DataFrame(
        {
            "sample_id": cohort_meta["sample_id"].astype(str),
            "expression_sample_id": cohort_meta["sample_id"].astype(str),
            "patient_id": cohort_meta["patient_name"].astype(str).str.replace(r"_T_[A-Z]+$", "", regex=True),
            "patient_id_raw": cohort_meta["patient_name"].astype(str),
            "cohort": "PHS000452_LIU_LIKE_PRE",
            "accession": "phs000452",
            "dataset_id": "Melanoma-phs000452",
            "response_raw": cohort_meta["response"].astype(str),
            "response_NR": cohort_meta["response_NR"].astype(str),
            "label": pd.to_numeric(cohort_meta["label"], errors="coerce").astype("Int64"),
            "label_definition": "0=CR/PR/MR or TIGER response_NR=R; 1=SD/PD or response_NR=N",
            "timepoint": "pretreatment_inferred_from_tiger_patient_like_antipd1_dataset",
            "therapy": "anti-PD-1",
            "treatment": "anti-PD-1",
            "tumor_type": cohort_meta["tumor_type"].astype(str),
            "seq_type": cohort_meta["seq_type"].astype(str),
            "m_stage": cohort_meta["M Stage"].astype(str),
            "overall_survival_days": cohort_meta["overall survival (days)"],
            "vital_status": cohort_meta["vital status"].astype(str),
            "sex": cohort_meta["Gender"].astype(str),
            "match_method": "tiger_sample_id_patient_like_subset",
            "evidence_source_file": str(clinical_path),
            "expression_source_file": str(expression_path),
            "evidence_field": "TIGER clinical sample_id/patient_name/Therapy/response/response_NR; Patient-like rows kept separate from IPI-like rows",
            "tiger_download_expression_url": "http://tiger.canceromics.org/tiger/Download/immunotherapy/expression/tsv/Melanoma-phs000452.Response.tsv",
            "tiger_download_clinical_url": "http://tiger.canceromics.org/tiger/Download/immunotherapy/clinical/tsv/Melanoma-phs000452.Response.tsv",
        }
    ).reset_index(drop=True)
    return X, meta


def main() -> None:
    parser = argparse.ArgumentParser(description="Preprocess TIGER Melanoma-phs000452 Patient-like anti-PD-1 samples.")
    parser.add_argument("--expression", default="data/external/TIGER_phs000452/Melanoma-phs000452.expression.tsv")
    parser.add_argument("--clinical", default="data/external/TIGER_phs000452/Melanoma-phs000452.clinical.tsv")
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--report", default="tables/phs000452_tiger_qc.tsv")
    args = parser.parse_args()

    expression_path = ROOT / args.expression
    clinical_path = ROOT / args.clinical
    processed_dir = ROOT / args.processed_dir
    processed_dir.mkdir(parents=True, exist_ok=True)

    expression = load_tiger_expression(expression_path)
    clinical = pd.read_csv(clinical_path, sep="\t")
    X, meta = build_patient_like_cohort(expression, clinical, expression_path, clinical_path)

    cohort = "PHS000452_LIU_LIKE_PRE"
    X.to_csv(processed_dir / f"{cohort}.expr.tsv", sep="\t")
    meta.to_csv(processed_dir / f"{cohort}.metadata.tsv", sep="\t", index=False)

    report = pd.DataFrame(
        [
            {
                "cohort": cohort,
                "n_samples": X.shape[0],
                "n_genes": X.shape[1],
                "n_patients": int(meta["patient_id"].nunique()),
                "n_responders_CR_PR_MR": int(meta["response_raw"].isin(["CR", "PR", "MR"]).sum()),
                "n_nonresponders_SD_PD": int(meta["response_raw"].isin(["SD", "PD"]).sum()),
                "response_raw_counts": ";".join(f"{k}:{v}" for k, v in meta["response_raw"].value_counts().sort_index().items()),
                "excluded_ipi_like_rows": int(clinical["patient_name"].astype(str).str.startswith("IPI").sum()),
                "status": "processed_patient_like_subset",
                "source_expression": str(expression_path),
                "source_clinical": str(clinical_path),
            }
        ]
    )
    report_path = ROOT / args.report
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(report_path, sep="\t", index=False)
    print(report.to_string(index=False))
    print(f"Wrote phs000452 TIGER QC report to {report_path}")


if __name__ == "__main__":
    main()
