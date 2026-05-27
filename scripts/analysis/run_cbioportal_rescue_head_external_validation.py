from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche_opt.model.endpoint_modules import _concat, endpoint_label_series, select_threshold
from scripts.analysis.run_strict_external_failure_mode_audit import (
    METHODS,
    PRIMARY_COHORTS,
    TARGET_GENES,
    blend_score,
    blend_specs,
    build_method_scores,
    evaluate_primary_lodo,
    _metric_row,
)


CBIO_EXTERNAL_GROUPS: dict[str, dict[str, object]] = {
    "cbio_liu_dfci_only": {
        "cohorts": ["CBIO_LIU_DFCI_2019_PRE"],
        "claim_boundary": "independent_source_crosscheck_no_external_model_selection",
        "notes": "Original public cBioPortal Liu/DFCI melanoma anti-PD1 pretreatment source. It is a source-level cross-check for the Liu/MGSP-like external cohort.",
    },
    "strict_cbio_liu_plus_gse145996": {
        "cohorts": ["CBIO_LIU_DFCI_2019_PRE", "GSE145996"],
        "claim_boundary": "strict_external_crosscheck_no_refit",
        "notes": "Strict external cross-check combining cBioPortal Liu/DFCI and GSE145996. External labels are not used for training, feature selection, thresholding or calibration.",
    },
    "cbio_iatlas_liu_duplicate_crosscheck": {
        "cohorts": ["CBIO_IATLAS_LIU_2019_PRE"],
        "claim_boundary": "duplicate_source_crosscheck_not_independent",
        "notes": "iAtlas Liu harmonized profile shares source patients with Liu/DFCI and is retained only as a duplicate-source processing check.",
    },
}


def read_targeted_cohort(processed_dir: Path, cohort: str) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    expr_path = processed_dir / f"{cohort}.expr.tsv"
    meta_path = processed_dir / f"{cohort}.metadata.tsv"
    if not expr_path.exists() or not meta_path.exists():
        return None
    header = pd.read_csv(expr_path, sep="\t", nrows=0).columns.tolist()
    if not header:
        return None
    sample_col = header[0]
    available = [gene for gene in TARGET_GENES if gene in header]
    expression = pd.read_csv(expr_path, sep="\t", usecols=[sample_col, *available], index_col=0)
    metadata = pd.read_csv(meta_path, sep="\t").set_index("sample_id", drop=False).reindex(expression.index)
    metadata = metadata[metadata["label"].notna()]
    expression = expression.loc[metadata.index]
    return expression, metadata


def load_targeted_dirs(processed_dir: Path, cbio_dir: Path) -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame], pd.DataFrame]:
    cohorts = sorted(
        {
            *PRIMARY_COHORTS,
            "GSE145996",
            *[cohort for group in CBIO_EXTERNAL_GROUPS.values() for cohort in group["cohorts"]],
        }
    )
    rows: list[dict[str, object]] = []
    X_by_cohort: dict[str, pd.DataFrame] = {}
    metadata_by_cohort: dict[str, pd.DataFrame] = {}
    for cohort in cohorts:
        source_dir = cbio_dir if cohort.startswith("CBIO_") else processed_dir
        loaded = read_targeted_cohort(source_dir, cohort)
        if loaded is None:
            rows.append(
                {
                    "cohort": cohort,
                    "source_dir": str(source_dir.relative_to(ROOT)),
                    "n_samples": 0,
                    "n_target_genes_available": 0,
                    "target_genes_available": "",
                    "status": "missing_expression_or_metadata",
                }
            )
            continue
        expression, metadata = loaded
        available = [gene for gene in TARGET_GENES if gene in expression.columns]
        X_by_cohort[cohort] = expression
        metadata_by_cohort[cohort] = metadata
        rows.append(
            {
                "cohort": cohort,
                "source_dir": str(source_dir.relative_to(ROOT)),
                "n_samples": int(expression.shape[0]),
                "n_target_genes_available": int(len(available)),
                "target_genes_available": ",".join(available),
                "status": "ready" if len(available) == len(TARGET_GENES) else "partial_gene_coverage",
            }
        )
    return X_by_cohort, metadata_by_cohort, pd.DataFrame(rows)


def prepare_endpoint(
    X_by_cohort: dict[str, pd.DataFrame],
    metadata_by_cohort: dict[str, pd.DataFrame],
    cohorts: list[str],
    endpoint: str,
) -> tuple[dict[str, pd.DataFrame], dict[str, pd.Series]]:
    X_out: dict[str, pd.DataFrame] = {}
    y_out: dict[str, pd.Series] = {}
    for cohort in cohorts:
        if cohort not in X_by_cohort or cohort not in metadata_by_cohort:
            continue
        y = endpoint_label_series(metadata_by_cohort[cohort]["response_raw"], endpoint)
        mask = y.notna()
        if int(mask.sum()) < 4 or y.loc[mask].nunique() < 2:
            continue
        X_out[cohort] = X_by_cohort[cohort].loc[mask].copy()
        y_out[cohort] = y.loc[mask].astype(int)
    return X_out, y_out


def evaluate_external_groups(
    scores_by_method: dict[str, dict[str, pd.Series]],
    y_by_cohort: dict[str, pd.Series],
    specs: list[dict[str, object]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pooled_rows: list[dict[str, object]] = []
    per_cohort_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    train = [cohort for cohort in PRIMARY_COHORTS if cohort in y_by_cohort]
    y_train = _concat(y_by_cohort, train).astype(int)
    for spec in specs:
        p_train = _concat({cohort: blend_score(scores_by_method, cohort, spec) for cohort in train}, train).astype(float)
        threshold = select_threshold(y_train.to_numpy(dtype=int), p_train.to_numpy(dtype=float))
        for group_id, group in CBIO_EXTERNAL_GROUPS.items():
            parts: list[dict[str, object]] = []
            for cohort in group["cohorts"]:
                cohort = str(cohort)
                if cohort not in y_by_cohort:
                    continue
                y_test = y_by_cohort[cohort].astype(int)
                p_test = blend_score(scores_by_method, cohort, spec).reindex(y_test.index).astype(float)
                for sample_id, probability in p_test.items():
                    parts.append(
                        {
                            "endpoint": "strict_recist",
                            "group_id": group_id,
                            "cohort": cohort,
                            "sample_id": sample_id,
                            "true_response_label": int(y_test.loc[sample_id]),
                            "response_probability": float(probability),
                            "threshold": float(threshold),
                            "blend_id": spec["blend_id"],
                            "blend_type": spec["blend_type"],
                            "axis": "MAP4K1_minus_TBX3_AXL",
                            "claim_boundary": group["claim_boundary"],
                        }
                    )
            if not parts:
                continue
            pred = pd.DataFrame(parts)
            pooled_rows.append(
                {
                    "endpoint": "strict_recist",
                    "group_id": group_id,
                    "blend_id": spec["blend_id"],
                    "blend_type": spec["blend_type"],
                    "axis": "MAP4K1_minus_TBX3_AXL",
                    "claim_boundary": group["claim_boundary"],
                    "group_notes": group["notes"],
                    **_metric_row(
                        pred["true_response_label"].astype(int),
                        pred["response_probability"].astype(float),
                        float(threshold),
                    ),
                }
            )
            for cohort, cohort_pred in pred.groupby("cohort"):
                per_cohort_rows.append(
                    {
                        "endpoint": "strict_recist",
                        "group_id": group_id,
                        "cohort": cohort,
                        "blend_id": spec["blend_id"],
                        "blend_type": spec["blend_type"],
                        "axis": "MAP4K1_minus_TBX3_AXL",
                        "claim_boundary": group["claim_boundary"],
                        **_metric_row(
                            cohort_pred["true_response_label"].astype(int),
                            cohort_pred["response_probability"].astype(float),
                            float(threshold),
                        ),
                    }
                )
            prediction_rows.extend(parts)
    return pd.DataFrame(pooled_rows), pd.DataFrame(per_cohort_rows), pd.DataFrame(prediction_rows)


def build_selection(primary: pd.DataFrame, external: pd.DataFrame, per_cohort: pd.DataFrame) -> pd.DataFrame:
    if primary.empty or external.empty:
        return pd.DataFrame()
    primary_best = primary.sort_values(["AUROC", "AUPRC", "balanced_accuracy"], ascending=False).iloc[0]
    robust_blend = "0.05*cohort_zscore+0.95*cohort_robust_zscore"
    rows: list[dict[str, object]] = []
    for group_id, group_frame in external.groupby("group_id"):
        primary_match = group_frame[group_frame["blend_id"].astype(str).eq(str(primary_best["blend_id"]))]
        robust = group_frame[group_frame["blend_id"].astype(str).eq(robust_blend)]
        stress = group_frame.sort_values(["AUROC", "AUPRC", "balanced_accuracy"], ascending=False).iloc[0]
        selected = [
            ("primary_auc_selected_blend", primary_match.iloc[0] if not primary_match.empty else stress, "selected_by_primary_lodo_only_not_by_cbio_external"),
            ("robust_fixed_development_candidate", robust.iloc[0] if not robust.empty else stress, "fixed_robust_transform_candidate_no_cbio_external_label_fit"),
            ("current_cbio_external_stress_best", stress, "diagnostic_current_cbio_external_stress_screen_not_selection_claim"),
        ]
        for selection_id, row, boundary in selected:
            cohort_rows = per_cohort[
                per_cohort["group_id"].astype(str).eq(str(group_id))
                & per_cohort["blend_id"].astype(str).eq(str(row["blend_id"]))
            ]
            per_text = ";".join(
                f"{r['cohort']}:{float(r['AUROC']):.3f}"
                for _, r in cohort_rows.sort_values("cohort").iterrows()
            )
            primary_row = primary[primary["blend_id"].astype(str).eq(str(row["blend_id"]))]
            primary_row = primary_row.iloc[0] if not primary_row.empty else pd.Series(dtype=object)
            rows.append(
                {
                    "selection_id": selection_id,
                    "group_id": group_id,
                    "claim_boundary": boundary,
                    "blend_id": row["blend_id"],
                    "blend_type": row["blend_type"],
                    "primary_AUROC": float(primary_row.get("AUROC", float("nan"))),
                    "primary_AUPRC": float(primary_row.get("AUPRC", float("nan"))),
                    "primary_balanced_accuracy": float(primary_row.get("balanced_accuracy", float("nan"))),
                    "strict_external_AUROC": float(row["AUROC"]),
                    "strict_external_AUPRC": float(row["AUPRC"]),
                    "strict_external_balanced_accuracy": float(row["balanced_accuracy"]),
                    "strict_external_ECE": float(row.get("ECE", float("nan"))),
                    "strict_external_per_cohort_AUROC": per_text,
                }
            )
    return pd.DataFrame(rows)


def write_markdown(selection: pd.DataFrame, coverage: pd.DataFrame, out_md: Path) -> None:
    lines = [
        "# cBioPortal Rescue-Head External Validation Audit",
        "",
        "This audit checks whether the frozen MAP4K1-TBX3/AXL rescue head can be fairly recomputed on cBioPortal Liu/DFCI after explicitly requesting all three target genes.",
        "Primary-model selection and thresholding use only discovery cohorts. Rows marked as stress screens are diagnostic and cannot be used as locked external model selection.",
        "",
        "## Target-Gene Coverage",
        "",
    ]
    for _, row in coverage.sort_values("cohort").iterrows():
        lines.append(
            f"- {row['cohort']}: {int(row['n_target_genes_available'])}/{len(TARGET_GENES)} target genes "
            f"({row['target_genes_available'] or 'none'}); status={row['status']}."
        )
    lines.extend(["", "## Selection Summary", ""])
    if selection.empty:
        lines.append("No selection rows were produced.")
    else:
        for _, row in selection.iterrows():
            lines.append(
                "- `{}` / `{}`: blend={}; primary AUROC={:.3f}; strict external AUROC={:.3f}; "
                "BA={:.3f}; per-cohort {}; boundary={}".format(
                    row["selection_id"],
                    row["group_id"],
                    row["blend_id"],
                    float(row["primary_AUROC"]),
                    float(row["strict_external_AUROC"]),
                    float(row["strict_external_balanced_accuracy"]),
                    row["strict_external_per_cohort_AUROC"],
                    row["claim_boundary"],
                )
            )
    out_md.write_text("\n".join(lines), encoding="utf-8")


def run_audit(processed_dir: Path, cbio_dir: Path, out_dir: Path) -> dict[str, Path]:
    X_all, metadata_all, coverage = load_targeted_dirs(processed_dir, cbio_dir)
    primary_X, primary_y = prepare_endpoint(X_all, metadata_all, PRIMARY_COHORTS, "primary_recist")
    strict_X, strict_y = prepare_endpoint(
        X_all,
        metadata_all,
        sorted({*PRIMARY_COHORTS, "GSE145996", *[cohort for group in CBIO_EXTERNAL_GROUPS.values() for cohort in group["cohorts"]]}),
        "strict_recist",
    )
    primary_scores = build_method_scores(primary_X)
    strict_scores = build_method_scores(strict_X)
    specs = blend_specs()
    primary, primary_predictions = evaluate_primary_lodo(primary_scores, primary_y, specs)
    external, per_cohort, external_predictions = evaluate_external_groups(strict_scores, strict_y, specs)
    selection = build_selection(primary, external, per_cohort)
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "coverage": out_dir / "cbioportal_rescue_head_gene_coverage.tsv",
        "primary_summary": out_dir / "cbioportal_rescue_head_primary_summary.tsv",
        "external_summary": out_dir / "cbioportal_rescue_head_external_summary.tsv",
        "external_per_cohort": out_dir / "cbioportal_rescue_head_external_per_cohort.tsv",
        "selection": out_dir / "cbioportal_rescue_head_selection.tsv",
        "primary_predictions": out_dir / "cbioportal_rescue_head_primary_predictions.tsv",
        "external_predictions": out_dir / "cbioportal_rescue_head_external_predictions.tsv",
        "markdown": out_dir / "CBIOPORTAL_RESCUE_HEAD_EXTERNAL_AUDIT.md",
    }
    coverage.to_csv(outputs["coverage"], sep="\t", index=False)
    primary.to_csv(outputs["primary_summary"], sep="\t", index=False)
    external.to_csv(outputs["external_summary"], sep="\t", index=False)
    per_cohort.to_csv(outputs["external_per_cohort"], sep="\t", index=False)
    selection.to_csv(outputs["selection"], sep="\t", index=False)
    primary_predictions.to_csv(outputs["primary_predictions"], sep="\t", index=False)
    external_predictions.to_csv(outputs["external_predictions"], sep="\t", index=False)
    write_markdown(selection, coverage, outputs["markdown"])
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--cbio-dir", default="data/processed/cbioportal_melanoma")
    parser.add_argument("--out", default="results/cbioportal_rescue_head_external_validation_20260527")
    args = parser.parse_args()
    outputs = run_audit(ROOT / args.processed_dir, ROOT / args.cbio_dir, ROOT / args.out)
    selection = pd.read_csv(outputs["selection"], sep="\t")
    print(json.dumps(selection.to_dict("records"), ensure_ascii=False))
    print(f"Wrote {ROOT / args.out}")


if __name__ == "__main__":
    main()
