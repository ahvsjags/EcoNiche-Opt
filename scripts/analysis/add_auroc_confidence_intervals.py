from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn import metrics as sk_metrics


ROOT = Path(__file__).resolve().parents[2]
PRIMARY_PREDICTIONS = ROOT / "results" / "endpoint_modules_heuristic_core_locked_gpu" / "endpoint_module_predictions.tsv"
EXTERNAL_PREDICTIONS = (
    ROOT / "results" / "locked_external_panel_validation_calibrated_20260519" / "locked_external_predictions.tsv"
)
TABLE10 = ROOT / "tables" / "article" / "supp_table_10_melanoma_benchmark_summary.tsv"
TABLE11 = ROOT / "tables" / "article" / "supp_table_11_signature_family_fdr.tsv"
TABLE15 = ROOT / "tables" / "article" / "supp_table_15_locked_external_metrics.tsv"
EXTERNAL_FAMILY = (
    ROOT / "results" / "locked_external_panel_validation_calibrated_20260519" / "locked_external_signature_family_omnibus.tsv"
)


def bootstrap_auroc_ci(y_true: pd.Series, prob: pd.Series, n_bootstrap: int, seed: int) -> dict[str, float]:
    y = y_true.to_numpy(dtype=int)
    p = prob.to_numpy(dtype=float)
    if len(y) < 4 or len(np.unique(y)) < 2:
        return {"AUROC": np.nan, "AUROC_ci_low": np.nan, "AUROC_ci_high": np.nan, "AUROC_ci_bootstrap_n": 0}
    observed = float(sk_metrics.roc_auc_score(y, p))
    rng = np.random.default_rng(seed)
    values: list[float] = []
    for _ in range(n_bootstrap):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        values.append(float(sk_metrics.roc_auc_score(y[idx], p[idx])))
    if not values:
        return {"AUROC": observed, "AUROC_ci_low": np.nan, "AUROC_ci_high": np.nan, "AUROC_ci_bootstrap_n": 0}
    arr = np.asarray(values)
    return {
        "AUROC": observed,
        "AUROC_ci_low": float(np.quantile(arr, 0.025)),
        "AUROC_ci_high": float(np.quantile(arr, 0.975)),
        "AUROC_ci_bootstrap_n": int(len(arr)),
    }


def primary_ci(predictions: pd.DataFrame, n_bootstrap: int) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    group_cols = ["endpoint", "stratum", "model_name"]
    for key, frame in predictions.groupby(group_cols):
        endpoint, stratum, model_name = key
        y = frame["true_response_label"].astype(int)
        prob = frame["response_probability"].astype(float)
        stats = bootstrap_auroc_ci(y, prob, n_bootstrap, seed=20260527 + len(rows))
        rows.append(
            {
                "endpoint": endpoint,
                "stratum": stratum,
                "model_name": model_name,
                "n_samples_ci": int(len(frame)),
                "n_responders_ci": int(y.sum()),
                "n_nonresponders_ci": int((y == 0).sum()),
                **stats,
            }
        )
    return pd.DataFrame(rows)


def external_ci(predictions: pd.DataFrame, n_bootstrap: int) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    group_cols = ["endpoint", "cohort", "analysis_type", "model_name"]
    for key, frame in predictions.groupby(group_cols):
        endpoint, cohort, analysis_type, model_name = key
        y = frame["true_response_label"].astype(int)
        prob = frame["response_probability"].astype(float)
        stats = bootstrap_auroc_ci(y, prob, n_bootstrap, seed=20260627 + len(rows))
        rows.append(
            {
                "endpoint": endpoint,
                "cohort": cohort,
                "analysis_type": analysis_type,
                "model_name": model_name,
                "n_samples_ci": int(len(frame)),
                "n_responders_ci": int(y.sum()),
                "n_nonresponders_ci": int((y == 0).sum()),
                **stats,
            }
        )
    return pd.DataFrame(rows)


def external_family_target_ci(external_predictions: pd.DataFrame, family: pd.DataFrame, n_bootstrap: int) -> pd.DataFrame:
    family_masks = {
        "all_locked_external_and_panel": lambda frame: pd.Series(True, index=frame.index),
        "melanoma_external_and_panel": lambda frame: frame["cohort"].isin(
            ["GSE145996", "PHS000452_LIU_LIKE_PRE", "PRJEB23709_COMBO_PRE", "GSE93157"]
        ),
        "strict_pd1_like_external": lambda frame: frame["cohort"].isin(["GSE145996", "PHS000452_LIU_LIKE_PRE"]),
        "nanostring_panel_transfer": lambda frame: frame["cohort"].isin(["GSE93157", "GSE140901"]),
    }
    target_model = str(family["target_model"].iloc[0]) if not family.empty and "target_model" in family.columns else ""
    rows: list[dict[str, object]] = []
    for _, fam_row in family.iterrows():
        endpoint = str(fam_row["endpoint"])
        validation_family = str(fam_row["validation_family"])
        if validation_family not in family_masks:
            continue
        frame = external_predictions[
            external_predictions["endpoint"].astype(str).eq(endpoint)
            & external_predictions["model_name"].astype(str).eq(target_model)
        ].copy()
        frame = frame[family_masks[validation_family](frame)]
        if frame.empty:
            continue
        y = frame["true_response_label"].astype(int)
        prob = frame["response_probability"].astype(float)
        stats = bootstrap_auroc_ci(y, prob, n_bootstrap, seed=20260727 + len(rows))
        rows.append(
            {
                "endpoint": endpoint,
                "validation_family": validation_family,
                "target_model": target_model,
                "n_samples_target_ci": int(len(frame)),
                "target_AUROC_ci_low": stats["AUROC_ci_low"],
                "target_AUROC_ci_high": stats["AUROC_ci_high"],
                "target_AUROC_ci_bootstrap_n": stats["AUROC_ci_bootstrap_n"],
            }
        )
    return pd.DataFrame(rows)


def merge_table10(ci: pd.DataFrame) -> pd.DataFrame:
    table = pd.read_csv(TABLE10, sep="\t")
    table = table.drop(
        columns=[
            col
            for col in [
                "n_responders",
                "n_nonresponders",
                "pooled_AUROC_ci_low",
                "pooled_AUROC_ci_high",
                "pooled_AUROC_ci_bootstrap_n",
            ]
            if col in table.columns
        ]
    )
    merged = table.merge(
        ci[
            [
                "endpoint",
                "stratum",
                "model_name",
                "n_responders_ci",
                "n_nonresponders_ci",
                "AUROC_ci_low",
                "AUROC_ci_high",
                "AUROC_ci_bootstrap_n",
            ]
        ],
        on=["endpoint", "stratum", "model_name"],
        how="left",
    )
    return merged.rename(
        columns={
            "n_responders_ci": "n_responders",
            "n_nonresponders_ci": "n_nonresponders",
            "AUROC_ci_low": "pooled_AUROC_ci_low",
            "AUROC_ci_high": "pooled_AUROC_ci_high",
            "AUROC_ci_bootstrap_n": "pooled_AUROC_ci_bootstrap_n",
        }
    )


def merge_table11(ci: pd.DataFrame) -> pd.DataFrame:
    table = pd.read_csv(TABLE11, sep="\t")
    table = table.drop(
        columns=[
            col
            for col in [
                "target_AUROC_ci_low",
                "target_AUROC_ci_high",
                "target_AUROC_ci_bootstrap_n",
            ]
            if col in table.columns
        ]
    )
    target_ci = ci.rename(columns={"model_name": "target_model"})
    merged = table.merge(
        target_ci[
            [
                "endpoint",
                "stratum",
                "target_model",
                "AUROC_ci_low",
                "AUROC_ci_high",
                "AUROC_ci_bootstrap_n",
            ]
        ],
        on=["endpoint", "stratum", "target_model"],
        how="left",
    )
    return merged.rename(
        columns={
            "AUROC_ci_low": "target_AUROC_ci_low",
            "AUROC_ci_high": "target_AUROC_ci_high",
            "AUROC_ci_bootstrap_n": "target_AUROC_ci_bootstrap_n",
            "ci_low": "delta_AUROC_ci_low",
            "ci_high": "delta_AUROC_ci_high",
        }
    )


def merge_table15(ci: pd.DataFrame) -> pd.DataFrame:
    table = pd.read_csv(TABLE15, sep="\t")
    table = table.drop(
        columns=[col for col in ["AUROC_ci_low", "AUROC_ci_high", "AUROC_ci_bootstrap_n"] if col in table.columns]
    )
    merged = table.merge(
        ci[
            [
                "endpoint",
                "cohort",
                "analysis_type",
                "model_name",
                "AUROC_ci_low",
                "AUROC_ci_high",
                "AUROC_ci_bootstrap_n",
            ]
        ],
        on=["endpoint", "cohort", "analysis_type", "model_name"],
        how="left",
    )
    return merged


def write_audit(out_dir: Path, table10: pd.DataFrame, table11: pd.DataFrame, table15: pd.DataFrame, external_family: pd.DataFrame) -> pd.DataFrame:
    rows = [
        {
            "check": "table10_primary_benchmark_auroc_ci_columns",
            "is_valid": bool(
                {"pooled_AUROC_ci_low", "pooled_AUROC_ci_high"}.issubset(table10.columns)
                and table10["pooled_AUROC_ci_low"].notna().any()
            ),
            "detail": f"rows={len(table10)}",
        },
        {
            "check": "table11_signature_family_target_auroc_ci_columns",
            "is_valid": bool(
                {"target_AUROC_ci_low", "target_AUROC_ci_high", "delta_AUROC_ci_low", "delta_AUROC_ci_high"}.issubset(
                    table11.columns
                )
                and table11["target_AUROC_ci_low"].notna().all()
            ),
            "detail": f"rows={len(table11)}",
        },
        {
            "check": "table15_locked_external_auroc_ci_columns",
            "is_valid": bool(
                {"AUROC_ci_low", "AUROC_ci_high"}.issubset(table15.columns) and table15["AUROC_ci_low"].notna().any()
            ),
            "detail": f"rows={len(table15)}",
        },
        {
            "check": "external_family_target_auroc_ci_columns",
            "is_valid": bool(
                {"target_AUROC_ci_low", "target_AUROC_ci_high"}.issubset(external_family.columns)
                and external_family["target_AUROC_ci_low"].notna().any()
            ),
            "detail": f"rows={len(external_family)}",
        },
    ]
    report = pd.DataFrame(rows)
    report.to_csv(out_dir / "auroc_ci_audit.tsv", sep="\t", index=False)
    return report


def run(out_dir: Path, n_bootstrap: int) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    primary_predictions = pd.read_csv(PRIMARY_PREDICTIONS, sep="\t")
    external_predictions = pd.read_csv(EXTERNAL_PREDICTIONS, sep="\t")
    primary = primary_ci(primary_predictions, n_bootstrap)
    external = external_ci(external_predictions, n_bootstrap)
    primary.to_csv(out_dir / "primary_melanoma_auroc_ci.tsv", sep="\t", index=False)
    external.to_csv(out_dir / "locked_external_auroc_ci.tsv", sep="\t", index=False)

    table10 = merge_table10(primary)
    table11 = merge_table11(primary)
    table15 = merge_table15(external)
    table10.to_csv(TABLE10, sep="\t", index=False)
    table11.to_csv(TABLE11, sep="\t", index=False)
    table15.to_csv(TABLE15, sep="\t", index=False)

    family = pd.read_csv(EXTERNAL_FAMILY, sep="\t")
    family_target_ci = external_family_target_ci(external_predictions, family, n_bootstrap)
    if not family_target_ci.empty:
        family = family.merge(
            family_target_ci,
            on=["endpoint", "validation_family", "target_model"],
            how="left",
        )
    family = family.rename(columns={"ci_low": "delta_AUROC_ci_low", "ci_high": "delta_AUROC_ci_high"})
    family.to_csv(out_dir / "locked_external_signature_family_omnibus_with_ci.tsv", sep="\t", index=False)

    report = write_audit(out_dir, table10, table11, table15, family)
    print(report.to_string(index=False))
    if not bool(report["is_valid"].all()):
        raise SystemExit(1)
    print(f"Wrote AUROC CI outputs to {out_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Add bootstrap AUROC confidence intervals to article performance tables.")
    parser.add_argument("--out", default="results/performance_ci_audit_20260527")
    parser.add_argument("--n-bootstrap", type=int, default=1000)
    args = parser.parse_args()
    run(ROOT / args.out, args.n_bootstrap)


if __name__ == "__main__":
    main()
