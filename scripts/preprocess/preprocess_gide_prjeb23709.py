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


COHORT_DEFINITIONS = {
    "PRJEB23709_PD1_PRE": {
        "therapy_group": "anti-PD1_monotherapy",
        "therapy_label": "anti-PD-1",
    },
    "PRJEB23709_COMBO_PRE": {
        "therapy_group": "anti-PD1_plus_anti-CTLA4",
        "therapy_label": "anti-CTLA-4+anti-PD-1",
    },
}


def load_tiger_expression(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    if "GENE_SYMBOL" not in raw.columns:
        raise ValueError(f"{path} does not contain a GENE_SYMBOL column")
    raw = raw.set_index("GENE_SYMBOL")
    expression = clean_expression_matrix(raw)
    expression = log2_if_needed(expression)
    expression.index = expression.index.astype(str)
    return expression


def load_merged_clinical(clinical_path: Path, response_map_path: Path) -> pd.DataFrame:
    clinical = pd.read_csv(clinical_path, sep="\t")
    response_map = pd.read_csv(response_map_path, sep="\t")
    required_response_cols = {"run_accession", "patient_id", "therapy_group", "timepoint", "response_raw", "label"}
    missing = required_response_cols - set(response_map.columns)
    if missing:
        raise ValueError(f"{response_map_path} is missing required columns: {sorted(missing)}")
    merged = clinical.merge(
        response_map,
        left_on="sample_id",
        right_on="run_accession",
        how="inner",
        suffixes=("_tiger", "_curated"),
    )
    if len(merged) != len(response_map):
        raise ValueError(
            f"Only {len(merged)} of {len(response_map)} response-map rows matched TIGER clinical sample IDs"
        )
    mismatched = merged[
        merged["response"].astype(str).str.upper() != merged["response_raw"].astype(str).str.upper()
    ]
    if not mismatched.empty:
        examples = ",".join(mismatched["sample_id"].astype(str).head(5))
        raise ValueError(f"TIGER clinical and curated response map disagree for samples: {examples}")
    return merged


def build_cohort(
    expression: pd.DataFrame,
    merged: pd.DataFrame,
    cohort_name: str,
    therapy_group: str,
    therapy_label: str,
    expression_path: Path,
    clinical_path: Path,
    response_map_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cohort_meta = merged[
        (merged["timepoint"].astype(str).str.lower() == "pretreatment")
        & (merged["therapy_group"].astype(str) == therapy_group)
    ].copy()
    if cohort_meta.empty:
        raise ValueError(f"No pretreatment samples found for {therapy_group}")

    cohort_meta = cohort_meta.sort_values(["patient_id", "sample_id"]).drop_duplicates("patient_id", keep="first")
    sample_ids = cohort_meta["sample_id"].astype(str).tolist()
    missing_expr = sorted(set(sample_ids) - set(expression.index.astype(str)))
    if missing_expr:
        raise ValueError(f"{cohort_name} has samples without expression values: {missing_expr[:5]}")
    X = expression.loc[sample_ids].copy()

    meta = pd.DataFrame(
        {
            "sample_id": cohort_meta["sample_id"].astype(str),
            "expression_sample_id": cohort_meta["sample_id"].astype(str),
            "run_accession": cohort_meta["run_accession"].astype(str),
            "sample_accession": cohort_meta["sample_accession"].astype(str),
            "sample_title": cohort_meta["sample_title"].astype(str),
            "patient_id": cohort_meta["patient_id"].astype(str),
            "patient_id_raw": cohort_meta["patient_id"].astype(str),
            "patient_no": cohort_meta["patient_no"],
            "cohort": cohort_name,
            "accession": "PRJEB23709",
            "dataset_id": "Melanoma-PRJEB23709",
            "response_raw": cohort_meta["response_raw"].astype(str),
            "label": pd.to_numeric(cohort_meta["label"], errors="coerce").astype("Int64"),
            "label_definition": cohort_meta["label_definition"].astype(str),
            "timepoint": "pretreatment",
            "therapy_group": therapy_group,
            "therapy": therapy_label,
            "treatment": cohort_meta["treatment"].astype(str),
            "pfs_days": cohort_meta["pfs_days"],
            "os_days": cohort_meta["os_days"],
            "sex": cohort_meta["Gender"],
            "age_start": cohort_meta["age_start"],
            "tumor_type": cohort_meta["tumor_type"],
            "match_method": "tiger_sample_id_to_curated_ena_run_accession",
            "expression_source_file": str(expression_path),
            "clinical_source_file": str(clinical_path),
            "response_evidence_source_file": str(response_map_path),
            "response_evidence_original_tables": cohort_meta["source_clinical_table"].astype(str),
            "ena_evidence_source_file": cohort_meta["source_ena_table"].astype(str),
            "evidence_field": "run_accession/sample_title matched to Gide supplementary response tables and ENA metadata",
            "tiger_download_expression_url": "http://tiger.canceromics.org/tiger/Download/immunotherapy/expression/tsv/Melanoma-PRJEB23709.Response.tsv",
            "tiger_download_clinical_url": "http://tiger.canceromics.org/tiger/Download/immunotherapy/clinical/tsv/Melanoma-PRJEB23709.Response.tsv",
        }
    )
    meta = meta.reset_index(drop=True)
    return X, meta


def main() -> None:
    parser = argparse.ArgumentParser(description="Preprocess Gide/PRJEB23709 TIGER expression into PRE therapy-specific cohorts.")
    parser.add_argument("--expression", default="data/external/TIGER_PRJEB23709/Melanoma-PRJEB23709.expression.tsv")
    parser.add_argument("--clinical", default="data/external/TIGER_PRJEB23709/Melanoma-PRJEB23709.clinical.tsv")
    parser.add_argument("--response-map", default="data/external/PRJEB23709/gide_prjeb23709_response_map.tsv")
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--report", default="tables/gide_prjeb23709_pre_qc.tsv")
    args = parser.parse_args()

    expression_path = ROOT / args.expression
    clinical_path = ROOT / args.clinical
    response_map_path = ROOT / args.response_map
    processed_dir = ROOT / args.processed_dir
    processed_dir.mkdir(parents=True, exist_ok=True)

    expression = load_tiger_expression(expression_path)
    merged = load_merged_clinical(clinical_path, response_map_path)

    report_rows = []
    for cohort_name, definition in COHORT_DEFINITIONS.items():
        X, meta = build_cohort(
            expression=expression,
            merged=merged,
            cohort_name=cohort_name,
            therapy_group=definition["therapy_group"],
            therapy_label=definition["therapy_label"],
            expression_path=expression_path,
            clinical_path=clinical_path,
            response_map_path=response_map_path,
        )
        X.to_csv(processed_dir / f"{cohort_name}.expr.tsv", sep="\t")
        meta.to_csv(processed_dir / f"{cohort_name}.metadata.tsv", sep="\t", index=False)
        report_rows.append(
            {
                "cohort": cohort_name,
                "n_samples": X.shape[0],
                "n_genes": X.shape[1],
                "n_patients": int(meta["patient_id"].nunique()),
                "n_responders_CR_PR": int(meta["response_raw"].isin(["CR", "PR"]).sum()),
                "n_nonresponders_SD_PD": int(meta["response_raw"].isin(["SD", "PD"]).sum()),
                "response_raw_counts": ";".join(
                    f"{key}:{value}" for key, value in meta["response_raw"].value_counts().sort_index().items()
                ),
                "status": "processed",
                "source_expression": str(expression_path),
                "source_response_map": str(response_map_path),
            }
        )

    report = pd.DataFrame(report_rows)
    report_path = ROOT / args.report
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(report_path, sep="\t", index=False)
    print(report.to_string(index=False))
    print(f"Wrote PRJEB23709 PRE QC report to {report_path}")


if __name__ == "__main__":
    main()
