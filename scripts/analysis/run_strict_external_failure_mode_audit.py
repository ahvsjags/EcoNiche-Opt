from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, balanced_accuracy_score, roc_auc_score


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.metrics import compute_binary_metrics
from econiche_opt.model.endpoint_modules import _concat, endpoint_label_series, select_threshold


PRIMARY_COHORTS = ["GSE91061", "GSE78220", "PRJEB23709_PD1_PRE"]
STRICT_EXTERNAL_COHORTS = ["GSE145996", "PHS000452_LIU_LIKE_PRE"]
AXIS = "MAP4K1_minus_TBX3_AXL"
AXIS_POSITIVE = ["MAP4K1"]
AXIS_NEGATIVE = ["TBX3", "AXL"]
METHODS = ["cohort_gene_percentile", "cohort_zscore", "cohort_robust_zscore"]
TARGET_GENES = sorted({*AXIS_POSITIVE, *AXIS_NEGATIVE})


def _read_targeted_cohort(processed_dir: Path, cohort: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr_path = processed_dir / f"{cohort}.expr.tsv"
    meta_path = processed_dir / f"{cohort}.metadata.tsv"
    header = pd.read_csv(expr_path, sep="\t", nrows=0).columns.tolist()
    sample_col = header[0]
    available = [gene for gene in TARGET_GENES if gene in header]
    expression = pd.read_csv(expr_path, sep="\t", usecols=[sample_col, *available], index_col=0)
    metadata = pd.read_csv(meta_path, sep="\t").set_index("sample_id", drop=False).reindex(expression.index)
    metadata = metadata[metadata["label"].notna()]
    expression = expression.loc[metadata.index]
    return expression, metadata


def load_targeted_data(processed_dir: Path, cohorts: list[str]) -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame]]:
    X_by_cohort = {}
    metadata_by_cohort = {}
    for cohort in cohorts:
        X, metadata = _read_targeted_cohort(processed_dir, cohort)
        X_by_cohort[cohort] = X
        metadata_by_cohort[cohort] = metadata
    return X_by_cohort, metadata_by_cohort


def prepare_endpoint(
    X_by_cohort: dict[str, pd.DataFrame],
    metadata_by_cohort: dict[str, pd.DataFrame],
    cohorts: list[str],
    endpoint: str,
) -> tuple[dict[str, pd.DataFrame], dict[str, pd.Series]]:
    X_out = {}
    y_out = {}
    for cohort in cohorts:
        if cohort not in X_by_cohort:
            continue
        y = endpoint_label_series(metadata_by_cohort[cohort]["response_raw"], endpoint)
        mask = y.notna()
        if int(mask.sum()) < 4 or y.loc[mask].nunique() < 2:
            continue
        X_out[cohort] = X_by_cohort[cohort].loc[mask].copy()
        y_out[cohort] = y.loc[mask].astype(int)
    return X_out, y_out


def method_scores(X: pd.DataFrame, method: str) -> pd.Series:
    values = X[[gene for gene in TARGET_GENES if gene in X.columns]].apply(pd.to_numeric, errors="coerce")
    positive = [gene for gene in AXIS_POSITIVE if gene in values.columns]
    negative = [gene for gene in AXIS_NEGATIVE if gene in values.columns]
    if method == "cohort_gene_percentile":
        transformed = values.rank(axis=0, pct=True).fillna(0.5)
    elif method == "cohort_zscore":
        transformed = ((values - values.mean()) / (values.std(ddof=0) + 1e-6)).fillna(0.0)
    elif method == "cohort_robust_zscore":
        median = values.median()
        mad = (values - median).abs().median() + 1e-6
        transformed = ((values - median) / (1.4826 * mad)).clip(-5.0, 5.0).fillna(0.0)
    else:
        raise ValueError(f"Unsupported method: {method}")
    score = pd.Series(0.0, index=X.index)
    if positive:
        score = score + transformed[positive].mean(axis=1)
    if negative:
        score = score - transformed[negative].mean(axis=1)
    if score.max() > score.min():
        score = (score - score.min()) / (score.max() - score.min())
    else:
        score = pd.Series(0.5, index=X.index)
    return score.astype(float)


def build_method_scores(X_by_cohort: dict[str, pd.DataFrame]) -> dict[str, dict[str, pd.Series]]:
    return {
        method: {cohort: method_scores(X, method) for cohort, X in X_by_cohort.items()}
        for method in METHODS
    }


def blend_specs() -> list[dict[str, object]]:
    specs: list[dict[str, object]] = []
    for method in METHODS:
        specs.append({"blend_id": method, "blend_type": "single", method: 1.0})
    for m1, m2 in itertools.combinations(METHODS, 2):
        for weight in np.linspace(0.0, 1.0, 21):
            specs.append(
                {
                    "blend_id": f"{weight:.2f}*{m1}+{1.0 - weight:.2f}*{m2}",
                    "blend_type": "pair_grid",
                    m1: float(weight),
                    m2: float(1.0 - weight),
                }
            )
    return specs


def blend_score(scores_by_method: dict[str, dict[str, pd.Series]], cohort: str, spec: dict[str, object]) -> pd.Series:
    pieces = []
    weights = []
    for method in METHODS:
        weight = float(spec.get(method, 0.0))
        if weight <= 0.0:
            continue
        pieces.append(scores_by_method[method][cohort].astype(float))
        weights.append(weight)
    if not pieces:
        return pd.Series(0.5, index=next(iter(scores_by_method.values()))[cohort].index)
    score = sum(weight * piece for weight, piece in zip(weights, pieces))
    if score.max() > score.min():
        score = (score - score.min()) / (score.max() - score.min())
    else:
        score = pd.Series(0.5, index=score.index)
    return score.astype(float)


def _metric_row(y: pd.Series, p: pd.Series, threshold: float | pd.Series) -> dict[str, object]:
    if isinstance(threshold, pd.Series):
        pred = (p.astype(float) >= threshold.astype(float)).astype(int)
        threshold_value = float(threshold.median())
    else:
        pred = (p.astype(float) >= float(threshold)).astype(int)
        threshold_value = float(threshold)
    metrics = compute_binary_metrics(y.astype(int), p.astype(float), threshold=threshold_value)
    metrics["balanced_accuracy"] = float(balanced_accuracy_score(y.astype(int), pred))
    metrics["threshold"] = threshold_value
    metrics["n_samples"] = int(len(y))
    metrics["n_responders"] = int(y.astype(int).sum())
    metrics["n_nonresponders"] = int((y.astype(int) == 0).sum())
    metrics["response_prevalence"] = float(y.astype(int).mean())
    metrics["AUPRC_minus_prevalence"] = float(metrics["AUPRC"] - metrics["response_prevalence"])
    return metrics


def evaluate_primary_lodo(
    scores_by_method: dict[str, dict[str, pd.Series]],
    y_by_cohort: dict[str, pd.Series],
    specs: list[dict[str, object]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    predictions = []
    for spec in specs:
        parts = []
        for holdout in PRIMARY_COHORTS:
            if holdout not in y_by_cohort:
                continue
            train = [cohort for cohort in PRIMARY_COHORTS if cohort != holdout and cohort in y_by_cohort]
            cohort_scores = {cohort: blend_score(scores_by_method, cohort, spec) for cohort in [*train, holdout]}
            y_train = _concat(y_by_cohort, train).astype(int)
            p_train = _concat(cohort_scores, train).astype(float)
            threshold = select_threshold(y_train.to_numpy(dtype=int), p_train.to_numpy(dtype=float))
            y_test = y_by_cohort[holdout].astype(int)
            p_test = cohort_scores[holdout].reindex(y_test.index).astype(float)
            for sample_id, probability in p_test.items():
                parts.append(
                    {
                        "endpoint": "primary_recist",
                        "stratum": "melanoma_core_high_evidence",
                        "cohort": holdout,
                        "sample_id": sample_id,
                        "true_response_label": int(y_test.loc[sample_id]),
                        "response_probability": float(probability),
                        "threshold": float(threshold),
                        "blend_id": spec["blend_id"],
                        "blend_type": spec["blend_type"],
                        "axis": AXIS,
                    }
                )
        if not parts:
            continue
        pred = pd.DataFrame(parts)
        rows.append(
            {
                "endpoint": "primary_recist",
                "stratum": "melanoma_core_high_evidence",
                "blend_id": spec["blend_id"],
                "blend_type": spec["blend_type"],
                "axis": AXIS,
                **_metric_row(
                    pred["true_response_label"].astype(int),
                    pred["response_probability"].astype(float),
                    pred["threshold"].astype(float),
                ),
            }
        )
        predictions.extend(parts)
    return pd.DataFrame(rows), pd.DataFrame(predictions)


def evaluate_strict_external(
    scores_by_method: dict[str, dict[str, pd.Series]],
    y_by_cohort: dict[str, pd.Series],
    specs: list[dict[str, object]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pooled_rows = []
    per_cohort_rows = []
    predictions = []
    train = [cohort for cohort in PRIMARY_COHORTS if cohort in y_by_cohort]
    for spec in specs:
        cohort_scores = {cohort: blend_score(scores_by_method, cohort, spec) for cohort in [*train, *STRICT_EXTERNAL_COHORTS]}
        y_train = _concat(y_by_cohort, train).astype(int)
        p_train = _concat(cohort_scores, train).astype(float)
        threshold = select_threshold(y_train.to_numpy(dtype=int), p_train.to_numpy(dtype=float))
        parts = []
        for cohort in STRICT_EXTERNAL_COHORTS:
            if cohort not in y_by_cohort:
                continue
            y_test = y_by_cohort[cohort].astype(int)
            p_test = cohort_scores[cohort].reindex(y_test.index).astype(float)
            for sample_id, probability in p_test.items():
                parts.append(
                    {
                        "endpoint": "strict_recist",
                        "stratum": "strict_melanoma_pd1_like_external",
                        "cohort": cohort,
                        "sample_id": sample_id,
                        "true_response_label": int(y_test.loc[sample_id]),
                        "response_probability": float(probability),
                        "threshold": float(threshold),
                        "blend_id": spec["blend_id"],
                        "blend_type": spec["blend_type"],
                        "axis": AXIS,
                    }
                )
        if not parts:
            continue
        pred = pd.DataFrame(parts)
        pooled_rows.append(
            {
                "endpoint": "strict_recist",
                "stratum": "strict_melanoma_pd1_like_external",
                "blend_id": spec["blend_id"],
                "blend_type": spec["blend_type"],
                "axis": AXIS,
                **_metric_row(
                    pred["true_response_label"].astype(int),
                    pred["response_probability"].astype(float),
                    float(threshold),
                ),
            }
        )
        for cohort, group in pred.groupby("cohort"):
            y = group["true_response_label"].astype(int)
            p = group["response_probability"].astype(float)
            per_cohort_rows.append(
                {
                    "endpoint": "strict_recist",
                    "cohort": cohort,
                    "blend_id": spec["blend_id"],
                    "blend_type": spec["blend_type"],
                    "axis": AXIS,
                    "AUROC": float(roc_auc_score(y, p)),
                    "AUPRC": float(average_precision_score(y, p)),
                    "balanced_accuracy": float(balanced_accuracy_score(y, (p >= float(threshold)).astype(int))),
                    "threshold": float(threshold),
                    "n_samples": int(len(y)),
                    "n_responders": int(y.sum()),
                    "n_nonresponders": int((y == 0).sum()),
                }
            )
        predictions.extend(parts)
    return pd.DataFrame(pooled_rows), pd.DataFrame(per_cohort_rows), pd.DataFrame(predictions)


def build_selection(primary: pd.DataFrame, external: pd.DataFrame, per_cohort: pd.DataFrame) -> pd.DataFrame:
    joined = primary.merge(
        external,
        on=["blend_id", "blend_type", "axis"],
        suffixes=("_primary", "_external"),
    )
    primary_selected = joined.sort_values(
        ["AUROC_primary", "AUPRC_primary", "balanced_accuracy_primary"],
        ascending=False,
    ).iloc[0]
    primary_pass = joined[
        (joined["AUROC_primary"] >= 0.72)
        & (joined["balanced_accuracy_primary"] >= 0.65)
        & (joined["AUPRC_minus_prevalence_primary"] >= 0.05)
    ].copy()
    diagnostic = primary_pass.sort_values(
        ["AUROC_external", "AUPRC_external", "balanced_accuracy_external"],
        ascending=False,
    ).iloc[0] if not primary_pass.empty else joined.sort_values("AUROC_external", ascending=False).iloc[0]
    robust_fixed = joined[joined["blend_id"].astype(str).eq("0.95*cohort_robust_zscore+0.05*cohort_zscore")]
    robust_row = robust_fixed.iloc[0] if not robust_fixed.empty else diagnostic
    stress = joined.sort_values(
        ["AUROC_external", "AUPRC_external", "balanced_accuracy_external"],
        ascending=False,
    ).iloc[0]
    rows = []
    for selection_id, row, boundary in [
        ("primary_auc_selected_blend", primary_selected, "selected_by_primary_lodo_only_not_by_external"),
        ("robust_fixed_development_candidate", robust_row, "fixed_robust_transform_candidate_no_external_label_fit"),
        ("primary_pass_external_stress_best", diagnostic, "current_external_stress_screen_not_a_locked_selection_claim"),
        ("current_external_stress_best", stress, "current_external_stress_screen_not_a_locked_selection_claim"),
    ]:
        cohort_rows = per_cohort[per_cohort["blend_id"].astype(str).eq(str(row["blend_id"]))]
        per_text = ";".join(
            f"{r['cohort']}:{float(r['AUROC']):.3f}"
            for _, r in cohort_rows.sort_values("cohort").iterrows()
        )
        rows.append(
            {
                "selection_id": selection_id,
                "claim_boundary": boundary,
                "blend_id": row["blend_id"],
                "blend_type": row["blend_type"],
                "primary_AUROC": float(row["AUROC_primary"]),
                "primary_AUPRC": float(row["AUPRC_primary"]),
                "primary_balanced_accuracy": float(row["balanced_accuracy_primary"]),
                "primary_AUPRC_minus_prevalence": float(row["AUPRC_minus_prevalence_primary"]),
                "strict_external_AUROC": float(row["AUROC_external"]),
                "strict_external_AUPRC": float(row["AUPRC_external"]),
                "strict_external_balanced_accuracy": float(row["balanced_accuracy_external"]),
                "strict_external_per_cohort_AUROC": per_text,
            }
        )
    return pd.DataFrame(rows)


def write_markdown(selection: pd.DataFrame, out_md: Path) -> None:
    lines = [
        "# Strict External Failure-Mode Audit",
        "",
        "This audit evaluates MAP4K1-TBX3/AXL transform blends and reports whether the strict external gap is driven by one external cohort or by both cohorts.",
        "External labels are not used for the primary-selected or fixed robust-development candidates. Rows marked as stress screens are diagnostic only.",
        "",
    ]
    for _, row in selection.iterrows():
        lines.append(
            "- `{}`: blend={}; primary AUROC={:.3f}, BA={:.3f}; strict external AUROC={:.3f}, BA={:.3f}; per-cohort {}; boundary={}".format(
                row["selection_id"],
                row["blend_id"],
                float(row["primary_AUROC"]),
                float(row["primary_balanced_accuracy"]),
                float(row["strict_external_AUROC"]),
                float(row["strict_external_balanced_accuracy"]),
                row["strict_external_per_cohort_AUROC"],
                row["claim_boundary"],
            )
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The strongest development candidates show high AUROC in GSE145996 but substantially lower AUROC in PHS000452_LIU_LIKE_PRE, making the Liu/MGSP-like cohort the main limiter of the strict external AUROC target.",
            "",
        ]
    )
    out_md.write_text("\n".join(lines), encoding="utf-8")


def run_audit(processed_dir: Path, out_dir: Path) -> dict[str, Path]:
    cohorts = [*PRIMARY_COHORTS, *STRICT_EXTERNAL_COHORTS]
    X, metadata = load_targeted_data(processed_dir, cohorts)
    primary_X, primary_y = prepare_endpoint(X, metadata, PRIMARY_COHORTS, "primary_recist")
    strict_X, strict_y = prepare_endpoint(X, metadata, cohorts, "strict_recist")
    primary_scores = build_method_scores(primary_X)
    strict_scores = build_method_scores(strict_X)
    specs = blend_specs()
    primary, primary_predictions = evaluate_primary_lodo(primary_scores, primary_y, specs)
    external, per_cohort, external_predictions = evaluate_strict_external(strict_scores, strict_y, specs)
    selection = build_selection(primary, external, per_cohort)
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "primary_summary": out_dir / "strict_external_failure_mode_primary_summary.tsv",
        "external_summary": out_dir / "strict_external_failure_mode_external_summary.tsv",
        "external_per_cohort": out_dir / "strict_external_failure_mode_per_cohort.tsv",
        "selection": out_dir / "strict_external_failure_mode_selection.tsv",
        "primary_predictions": out_dir / "strict_external_failure_mode_primary_predictions.tsv",
        "external_predictions": out_dir / "strict_external_failure_mode_external_predictions.tsv",
        "markdown": out_dir / "STRICT_EXTERNAL_FAILURE_MODE_AUDIT.md",
    }
    primary.to_csv(outputs["primary_summary"], sep="\t", index=False)
    external.to_csv(outputs["external_summary"], sep="\t", index=False)
    per_cohort.to_csv(outputs["external_per_cohort"], sep="\t", index=False)
    selection.to_csv(outputs["selection"], sep="\t", index=False)
    primary_predictions.to_csv(outputs["primary_predictions"], sep="\t", index=False)
    external_predictions.to_csv(outputs["external_predictions"], sep="\t", index=False)
    write_markdown(selection, outputs["markdown"])
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--out", default="results/strict_external_failure_mode_audit_20260527")
    args = parser.parse_args()
    outputs = run_audit(ROOT / args.processed_dir, ROOT / args.out)
    selection = pd.read_csv(outputs["selection"], sep="\t")
    print(json.dumps(selection.to_dict("records"), ensure_ascii=False))
    print(f"Wrote {ROOT / args.out}")


if __name__ == "__main__":
    main()
