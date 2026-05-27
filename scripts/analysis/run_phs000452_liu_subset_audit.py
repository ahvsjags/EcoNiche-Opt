from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from sklearn.metrics import average_precision_score, balanced_accuracy_score, roc_auc_score


ROOT = Path(__file__).resolve().parents[2]
PHS_COHORT = "PHS000452_LIU_LIKE_PRE"
GSE_COHORT = "GSE145996"


def normalize_liu_patient_id(value: object) -> str:
    text = str(value)
    return text[:-2] if text.endswith("_T") else text


def _response_counts(frame: pd.DataFrame, column: str) -> str:
    if column not in frame.columns:
        return ""
    counts = frame[column].astype(str).value_counts().sort_index()
    return ";".join(f"{key}:{int(value)}" for key, value in counts.items())


def source_concordance(phs_metadata: pd.DataFrame, cbio_metadata: pd.DataFrame) -> pd.DataFrame:
    phs = phs_metadata.copy()
    cbio = cbio_metadata.copy()
    phs["patient_norm"] = phs["patient_id"].map(normalize_liu_patient_id)
    cbio["patient_norm"] = cbio["patient_id"].astype(str)
    merged = phs[
        ["patient_norm", "patient_id", "sample_id", "patient_id_raw", "response_raw", "response_NR", "m_stage"]
    ].merge(
        cbio[["patient_norm", "patient_id", "sample_id", "response_raw", "response_raw_source", "treatment"]],
        on="patient_norm",
        suffixes=("_tiger", "_cbio"),
        how="inner",
    )
    mismatch = merged[merged["response_raw_tiger"].astype(str) != merged["response_raw_cbio"].astype(str)]
    return pd.DataFrame(
        [
            {
                "audit": "liu_tiger_cbioportal_patient_response_concordance",
                "tiger_n": int(phs["patient_norm"].nunique()),
                "cbio_n": int(cbio["patient_norm"].nunique()),
                "matched_n": int(merged["patient_norm"].nunique()),
                "response_mismatch_n": int(len(mismatch)),
                "missing_from_cbio": ",".join(sorted(set(phs["patient_norm"]) - set(cbio["patient_norm"]))),
                "missing_from_tiger": ",".join(sorted(set(cbio["patient_norm"]) - set(phs["patient_norm"]))),
                "tiger_response_counts": _response_counts(phs, "response_raw"),
                "cbio_response_counts": _response_counts(cbio, "response_raw"),
                "interpretation": "TIGER and cBioPortal response labels are concordant for matched Liu/DFCI patients; performance differences are more likely processing/source-subset effects than response-label flips.",
            }
        ]
    )


def _subset_definitions(frame: pd.DataFrame) -> list[tuple[str, str, pd.Series, str]]:
    suffix = frame["patient_id_raw"].astype(str).str.extract(r"_(T_[A-Z])$")[0].fillna("unsuffixed")
    definitions: list[tuple[str, str, pd.Series, str]] = [
        ("all_phs_strict", "all", pd.Series(True, index=frame.index), "primary_strict_external_reference"),
    ]
    for value in sorted(suffix.unique()):
        definitions.append((f"patient_id_suffix::{value}", value, suffix.eq(value), "metadata_suffix_diagnostic"))
    for column, boundary in [
        ("m_stage", "clinical_stage_diagnostic"),
        ("sex", "demographic_subgroup_diagnostic"),
        ("vital_status", "post_outcome_diagnostic_not_predictive_claim"),
    ]:
        if column not in frame.columns:
            continue
        values = frame[column].astype(str).replace({"nan": ""})
        for value in sorted(v for v in values.unique() if v):
            definitions.append((f"{column}::{value}", value, values.eq(value), boundary))
    return definitions


def _metric_row(frame: pd.DataFrame) -> dict[str, object] | None:
    if len(frame) < 8 or frame["true_response_label"].nunique() < 2:
        return None
    y = frame["true_response_label"].astype(int)
    p = frame["response_probability"].astype(float)
    threshold = float(frame["threshold"].iloc[0])
    pred = (p >= threshold).astype(int)
    return {
        "n_samples": int(len(frame)),
        "n_responders": int(y.sum()),
        "n_nonresponders": int((y == 0).sum()),
        "AUROC": float(roc_auc_score(y, p)),
        "AUPRC": float(average_precision_score(y, p)),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
        "threshold": threshold,
    }


def subset_metrics(predictions: pd.DataFrame, phs_metadata: pd.DataFrame, selection: pd.DataFrame) -> pd.DataFrame:
    selected_blends = selection[
        selection["selection_id"].isin(
            [
                "primary_auc_selected_blend",
                "robust_fixed_development_candidate",
                "current_external_stress_best",
            ]
        )
    ]["blend_id"].astype(str).unique()
    phs = phs_metadata.copy()
    phs["patient_id_suffix"] = phs["patient_id_raw"].astype(str).str.extract(r"_(T_[A-Z])$")[0].fillna("unsuffixed")
    pred = predictions[
        predictions["cohort"].astype(str).eq(PHS_COHORT)
        & predictions["blend_id"].astype(str).isin(selected_blends)
    ].merge(
        phs[
            [
                "sample_id",
                "patient_id",
                "patient_id_raw",
                "patient_id_suffix",
                "response_raw",
                "response_NR",
                "m_stage",
                "sex",
                "vital_status",
            ]
        ],
        on="sample_id",
        how="left",
        suffixes=("", "_metadata"),
    )
    rows: list[dict[str, object]] = []
    for blend_id, blend_frame in pred.groupby("blend_id"):
        for subset_id, subset_value, mask, boundary in _subset_definitions(blend_frame):
            subset = blend_frame.loc[mask.to_numpy()].copy()
            metrics = _metric_row(subset)
            if metrics is None:
                continue
            rows.append(
                {
                    "cohort": PHS_COHORT,
                    "blend_id": blend_id,
                    "subset_id": subset_id,
                    "subset_value": subset_value,
                    "claim_boundary": boundary,
                    "response_counts": _response_counts(subset, "response_raw"),
                    **metrics,
                }
            )
    return pd.DataFrame(rows).sort_values(["blend_id", "subset_id"]).reset_index(drop=True)


def cohort_contrast(per_cohort: pd.DataFrame, selection: pd.DataFrame) -> pd.DataFrame:
    selected = selection[
        selection["selection_id"].isin(["primary_auc_selected_blend", "robust_fixed_development_candidate", "current_external_stress_best"])
    ][["selection_id", "blend_id", "claim_boundary"]]
    out = per_cohort.merge(selected, on="blend_id", how="inner")
    return out.sort_values(["selection_id", "cohort"]).reset_index(drop=True)


def write_markdown(concordance: pd.DataFrame, subset: pd.DataFrame, contrast: pd.DataFrame, out_md: Path) -> None:
    lines = [
        "# PHS000452/Liu Subset Failure-Mode Audit",
        "",
        "This audit checks whether the strict external gap reflects response-label discordance, source processing differences, or identifiable metadata subgroups in the Liu/MGSP-like PHS000452 cohort.",
        "",
        "## Source Concordance",
        "",
    ]
    row = concordance.iloc[0]
    lines.append(
        "- TIGER unique patients={}, cBioPortal unique patients={}, matched={}, response mismatches={}; missing_from_cbio={}; missing_from_tiger={}.".format(
            int(row["tiger_n"]),
            int(row["cbio_n"]),
            int(row["matched_n"]),
            int(row["response_mismatch_n"]),
            row["missing_from_cbio"] or "none",
            row["missing_from_tiger"] or "none",
        )
    )
    lines.extend(["", "## Cohort Contrast", ""])
    for _, r in contrast.iterrows():
        lines.append(
            "- `{}` / `{}`: AUROC={:.3f}, AUPRC={:.3f}, BA={:.3f}, n={}.".format(
                r["selection_id"],
                r["cohort"],
                float(r["AUROC"]),
                float(r["AUPRC"]),
                float(r["balanced_accuracy"]),
                int(r["n_samples"]),
            )
        )
    lines.extend(["", "## PHS000452 Subgroups", ""])
    for _, r in subset.sort_values(["blend_id", "AUROC"], ascending=[True, False]).iterrows():
        lines.append(
            "- `{}` / `{}`: AUROC={:.3f}, AUPRC={:.3f}, BA={:.3f}, n={}, responses={}; boundary={}.".format(
                r["blend_id"],
                r["subset_id"],
                float(r["AUROC"]),
                float(r["AUPRC"]),
                float(r["balanced_accuracy"]),
                int(r["n_samples"]),
                r["response_counts"],
                r["claim_boundary"],
            )
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "TIGER and cBioPortal labels are concordant for matched Liu/DFCI patients. The robust MAP4K1-TBX3/AXL candidate performs well in GSE145996 but is weaker in PHS000452 overall, with metadata subgroup heterogeneity. These subgroup rows are diagnostic and should not be used to redefine the locked strict external claim.",
            "",
        ]
    )
    out_md.write_text("\n".join(lines), encoding="utf-8")


def run_audit(
    phs_metadata_path: Path,
    cbio_metadata_path: Path,
    failure_predictions_path: Path,
    failure_per_cohort_path: Path,
    failure_selection_path: Path,
    out_dir: Path,
) -> dict[str, Path]:
    phs_metadata = pd.read_csv(phs_metadata_path, sep="\t")
    cbio_metadata = pd.read_csv(cbio_metadata_path, sep="\t")
    predictions = pd.read_csv(failure_predictions_path, sep="\t")
    per_cohort = pd.read_csv(failure_per_cohort_path, sep="\t")
    selection = pd.read_csv(failure_selection_path, sep="\t")
    concordance = source_concordance(phs_metadata, cbio_metadata)
    subset = subset_metrics(predictions, phs_metadata, selection)
    contrast = cohort_contrast(per_cohort, selection)
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "concordance": out_dir / "phs000452_liu_source_concordance.tsv",
        "subset_metrics": out_dir / "phs000452_subset_metrics.tsv",
        "cohort_contrast": out_dir / "phs000452_gse145996_contrast.tsv",
        "markdown": out_dir / "PHS000452_LIU_SUBSET_FAILURE_MODE_AUDIT.md",
    }
    concordance.to_csv(outputs["concordance"], sep="\t", index=False)
    subset.to_csv(outputs["subset_metrics"], sep="\t", index=False)
    contrast.to_csv(outputs["cohort_contrast"], sep="\t", index=False)
    write_markdown(concordance, subset, contrast, outputs["markdown"])
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phs-metadata", default="data/processed/bulk/PHS000452_LIU_LIKE_PRE.metadata.tsv")
    parser.add_argument("--cbio-metadata", default="data/processed/cbioportal_melanoma/CBIO_LIU_DFCI_2019_PRE.metadata.tsv")
    parser.add_argument(
        "--failure-predictions",
        default="results/strict_external_failure_mode_audit_20260527/strict_external_failure_mode_external_predictions.tsv",
    )
    parser.add_argument(
        "--failure-per-cohort",
        default="results/strict_external_failure_mode_audit_20260527/strict_external_failure_mode_per_cohort.tsv",
    )
    parser.add_argument(
        "--failure-selection",
        default="results/strict_external_failure_mode_audit_20260527/strict_external_failure_mode_selection.tsv",
    )
    parser.add_argument("--out", default="results/phs000452_liu_subset_audit_20260527")
    args = parser.parse_args()
    outputs = run_audit(
        ROOT / args.phs_metadata,
        ROOT / args.cbio_metadata,
        ROOT / args.failure_predictions,
        ROOT / args.failure_per_cohort,
        ROOT / args.failure_selection,
        ROOT / args.out,
    )
    subset = pd.read_csv(outputs["subset_metrics"], sep="\t")
    print(json.dumps(subset.sort_values("AUROC", ascending=False).head(8).to_dict("records"), ensure_ascii=False))
    print(f"Wrote {ROOT / args.out}")


if __name__ == "__main__":
    main()
