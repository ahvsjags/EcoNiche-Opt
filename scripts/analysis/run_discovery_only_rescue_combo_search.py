from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, balanced_accuracy_score, roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.baselines import BASELINE_SIGNATURES, signature_score
from econiche.io import load_processed_bulk
from econiche.metrics import compute_binary_metrics
from econiche.statistics import benjamini_hochberg
from econiche_opt.model.endpoint_modules import endpoint_label_series, select_threshold


PRIMARY_COHORTS = ["GSE91061", "GSE78220", "PRJEB23709_PD1_PRE"]
STRICT_EXTERNAL_COHORTS = ["GSE145996", "PHS000452_LIU_LIKE_PRE"]
BASE_POSITIVE = ["MAP4K1"]
BASE_NEGATIVE = ["TBX3", "AXL"]
TRANSFORMS = ["pct", "z", "rz"]
WEIGHTS = [0.20, 0.35, 0.50, 0.65, 0.80]
EIGHT_SIGNATURES = ["IFNG", "CXCL9", "TIG", "TIDE_dysfunction", "TIDE_exclusion", "CYT", "APM", "IPRES"]
_BASE_SCORE_CACHE: dict[str, pd.Series] = {}
_DIRECTION_CACHE: dict[tuple[str, str, tuple[str, ...]], float] = {}
_COMPONENT_CACHE: dict[tuple[str, str, str, float], pd.Series] = {}


def _safe_gene_symbol(gene: object) -> bool:
    text = str(gene)
    return bool(re.match(r"^[A-Za-z0-9.-]+$", text)) and len(text) <= 20


def labels_for_endpoint(
    X_by_cohort: dict[str, pd.DataFrame],
    metadata_by_cohort: dict[str, pd.DataFrame],
    cohorts: list[str],
    endpoint: str,
) -> dict[str, pd.Series]:
    labels: dict[str, pd.Series] = {}
    for cohort in cohorts:
        y = endpoint_label_series(metadata_by_cohort[cohort]["response_raw"], endpoint)
        mask = y.notna()
        common = X_by_cohort[cohort].index.intersection(y[mask].index)
        if len(common) >= 8 and y.loc[common].nunique() == 2:
            labels[cohort] = y.loc[common].astype(int)
    return labels


def primary_gene_universe(X_by_cohort: dict[str, pd.DataFrame]) -> list[str]:
    genes = set(X_by_cohort[PRIMARY_COHORTS[0]].columns)
    for cohort in PRIMARY_COHORTS[1:]:
        genes &= set(X_by_cohort[cohort].columns)
    genes = {gene for gene in genes if _safe_gene_symbol(gene)}
    genes.update(BASE_POSITIVE)
    genes.update(BASE_NEGATIVE)
    return sorted(gene for gene in genes if gene)


def needed_samples(primary_y: dict[str, pd.Series], strict_y: dict[str, pd.Series], cohort: str) -> list[str]:
    samples: list[str] = []
    for source in [primary_y.get(cohort), strict_y.get(cohort)]:
        if source is not None:
            samples.extend(source.index.astype(str).tolist())
    return list(dict.fromkeys(samples))


def build_transforms(
    X_by_cohort: dict[str, pd.DataFrame],
    primary_y: dict[str, pd.Series],
    strict_y: dict[str, pd.Series],
    genes: list[str],
    cohorts: list[str],
) -> dict[str, dict[str, pd.DataFrame]]:
    transforms: dict[str, dict[str, pd.DataFrame]] = {}
    for cohort in cohorts:
        samples = needed_samples(primary_y, strict_y, cohort)
        available = [gene for gene in genes if gene in X_by_cohort[cohort].columns]
        values = X_by_cohort[cohort].loc[samples, available].apply(pd.to_numeric, errors="coerce")
        pct = values.rank(axis=0, pct=True).fillna(0.5)
        z = ((values - values.mean()) / (values.std(ddof=0) + 1e-6)).fillna(0.0)
        median = values.median()
        mad = (values - median).abs().median() + 1e-6
        rz = ((values - median) / (1.4826 * mad)).clip(-5.0, 5.0).fillna(0.0)
        transforms[cohort] = {"pct": pct.astype("float32"), "z": z.astype("float32"), "rz": rz.astype("float32")}
    return transforms


def _auc_safe(y: pd.Series | np.ndarray, p: pd.Series | np.ndarray) -> float:
    y_arr = np.asarray(y, dtype=int)
    if len(np.unique(y_arr)) < 2:
        return float("nan")
    return float(roc_auc_score(y_arr, np.asarray(p, dtype=float)))


def gene_screen(
    transforms: dict[str, dict[str, pd.DataFrame]],
    y_by_cohort: dict[str, pd.Series],
    genes: list[str],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for method in TRANSFORMS:
        pred_parts: list[np.ndarray] = []
        y_parts: list[int] = []
        threshold_parts: list[np.ndarray] = []
        available = [gene for gene in genes if all(gene in transforms[cohort][method].columns for cohort in PRIMARY_COHORTS)]
        for holdout in PRIMARY_COHORTS:
            train = [cohort for cohort in PRIMARY_COHORTS if cohort != holdout]
            X_train = pd.concat([transforms[cohort][method].loc[y_by_cohort[cohort].index, available] for cohort in train])
            y_train = pd.concat([y_by_cohort[cohort] for cohort in train]).reindex(X_train.index).astype(int)
            direction = X_train[y_train == 1].mean(axis=0) - X_train[y_train == 0].mean(axis=0)
            sign = np.where(direction.to_numpy(dtype=float) >= 0.0, 1.0, -1.0).astype("float32")
            train_scores = X_train.to_numpy(dtype="float32") * sign
            test_scores = transforms[holdout][method].loc[y_by_cohort[holdout].index, available].to_numpy(dtype="float32") * sign
            y_arr = y_train.to_numpy(dtype=int)
            thresholds = np.quantile(train_scores, np.linspace(0.05, 0.95, 19), axis=0)
            best_ba = np.full(train_scores.shape[1], -1.0)
            best_threshold = np.zeros(train_scores.shape[1])
            for row in thresholds:
                pred = (train_scores >= row).astype(int)
                sensitivity = ((pred == 1) & (y_arr[:, None] == 1)).sum(axis=0) / max(1, int((y_arr == 1).sum()))
                specificity = ((pred == 0) & (y_arr[:, None] == 0)).sum(axis=0) / max(1, int((y_arr == 0).sum()))
                balanced = (sensitivity + specificity) / 2.0
                update = balanced > best_ba
                best_ba[update] = balanced[update]
                best_threshold[update] = row[update]
            pred_parts.append(test_scores)
            y_parts.extend(y_by_cohort[holdout].astype(int).tolist())
            threshold_parts.append(best_threshold)
        prediction_matrix = np.vstack(pred_parts)
        thresholds = np.vstack([np.tile(threshold, (pred_parts[idx].shape[0], 1)) for idx, threshold in enumerate(threshold_parts)])
        y = np.asarray(y_parts, dtype=int)
        binary = (prediction_matrix >= thresholds).astype(int)
        for idx, gene in enumerate(available):
            p = prediction_matrix[:, idx]
            rows.append(
                {
                    "method": method,
                    "gene": gene,
                    "primary_AUROC": _auc_safe(y, p),
                    "primary_AUPRC": float(average_precision_score(y, p)),
                    "primary_balanced_accuracy": float(balanced_accuracy_score(y, binary[:, idx])),
                }
            )
    return pd.DataFrame(rows).sort_values(["primary_AUROC", "primary_AUPRC", "primary_balanced_accuracy"], ascending=False).reset_index(drop=True)


def base_rescue_score(transforms: dict[str, dict[str, pd.DataFrame]], cohort: str) -> pd.Series:
    if cohort in _BASE_SCORE_CACHE:
        return _BASE_SCORE_CACHE[cohort]
    pieces: list[pd.Series] = []
    for method, weight in [("z", 0.05), ("rz", 0.95)]:
        frame = transforms[cohort][method]
        positive = frame[[gene for gene in BASE_POSITIVE if gene in frame.columns]].mean(axis=1)
        negative = frame[[gene for gene in BASE_NEGATIVE if gene in frame.columns]].mean(axis=1)
        score = positive - negative
        score = (score - score.min()) / (score.max() - score.min() + 1e-9)
        pieces.append(weight * score)
    out = sum(pieces)
    score = ((out - out.min()) / (out.max() - out.min() + 1e-9)).astype(float)
    _BASE_SCORE_CACHE[cohort] = score
    return score


def orient_component(
    scorer: Callable[[str], pd.Series],
    train_cohorts: list[str],
    y_by_cohort: dict[str, pd.Series],
) -> float:
    values = pd.concat([scorer(cohort).reindex(y_by_cohort[cohort].index) for cohort in train_cohorts])
    labels = pd.concat([y_by_cohort[cohort] for cohort in train_cohorts]).reindex(values.index).astype(int)
    return 1.0 if float(values[labels == 1].mean() - values[labels == 0].mean()) >= 0.0 else -1.0


def combo_specs(top_genes: pd.DataFrame, top_k: int) -> list[dict[str, object]]:
    rows = [{"candidate": "base_rescue_robust", "weight_base": 1.0, "component_method": "", "component_gene": ""}]
    selected = top_genes.head(top_k)[["method", "gene"]].drop_duplicates()
    for _, row in selected.iterrows():
        for weight in WEIGHTS:
            rows.append(
                {
                    "candidate": f"{weight:.2f}*base+{1.0 - weight:.2f}*{row['method']}__{row['gene']}",
                    "weight_base": float(weight),
                    "component_method": str(row["method"]),
                    "component_gene": str(row["gene"]),
                }
            )
    return rows


def score_combo(
    transforms: dict[str, dict[str, pd.DataFrame]],
    cohort: str,
    spec: dict[str, object],
    train_cohorts: list[str],
    y_train_by_cohort: dict[str, pd.Series],
) -> pd.Series:
    base = base_rescue_score(transforms, cohort)
    method = str(spec.get("component_method", ""))
    gene = str(spec.get("component_gene", ""))
    if not method or not gene:
        return base
    if gene not in transforms[cohort][method].columns:
        return base

    def component(inner_cohort: str) -> pd.Series:
        return transforms[inner_cohort][method][gene].astype(float)

    direction_key = (method, gene, tuple(train_cohorts))
    if direction_key not in _DIRECTION_CACHE:
        _DIRECTION_CACHE[direction_key] = orient_component(component, train_cohorts, y_train_by_cohort)
    sign = _DIRECTION_CACHE[direction_key]
    component_key = (cohort, method, gene, sign)
    if component_key in _COMPONENT_CACHE:
        component_score = _COMPONENT_CACHE[component_key]
    else:
        component_score = component(cohort) * sign
        component_score = (component_score - component_score.min()) / (component_score.max() - component_score.min() + 1e-9)
        _COMPONENT_CACHE[component_key] = component_score
    score = float(spec["weight_base"]) * base + (1.0 - float(spec["weight_base"])) * component_score
    return ((score - score.min()) / (score.max() - score.min() + 1e-9)).astype(float)


def evaluate_primary_lodo(
    transforms: dict[str, dict[str, pd.DataFrame]],
    y_by_cohort: dict[str, pd.Series],
    specs: list[dict[str, object]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    for spec in specs:
        parts: list[dict[str, object]] = []
        for holdout in PRIMARY_COHORTS:
            train = [cohort for cohort in PRIMARY_COHORTS if cohort != holdout]
            y_train = pd.concat([y_by_cohort[cohort] for cohort in train]).astype(int)
            p_train = pd.concat(
                [score_combo(transforms, cohort, spec, train, y_by_cohort).reindex(y_by_cohort[cohort].index) for cohort in train]
            ).reindex(y_train.index)
            threshold = select_threshold(y_train.to_numpy(dtype=int), p_train.to_numpy(dtype=float))
            y_test = y_by_cohort[holdout].astype(int)
            p_test = score_combo(transforms, holdout, spec, train, y_by_cohort).reindex(y_test.index)
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
                        "candidate": spec["candidate"],
                        "weight_base": spec["weight_base"],
                        "component_method": spec.get("component_method", ""),
                        "component_gene": spec.get("component_gene", ""),
                        "selection_boundary": "selected_by_primary_lodo_only",
                    }
                )
        pred = pd.DataFrame(parts)
        if pred.empty or pred["true_response_label"].nunique() < 2:
            continue
        metrics = compute_binary_metrics(pred["true_response_label"].astype(int), pred["response_probability"].astype(float), threshold=0.5)
        pred_label = pred["response_probability"].astype(float) >= pred["threshold"].astype(float)
        metrics["balanced_accuracy"] = float(balanced_accuracy_score(pred["true_response_label"].astype(int), pred_label.astype(int)))
        rows.append(
            {
                "endpoint": "primary_recist",
                "stratum": "melanoma_core_high_evidence",
                "candidate": spec["candidate"],
                "weight_base": spec["weight_base"],
                "component_method": spec.get("component_method", ""),
                "component_gene": spec.get("component_gene", ""),
                "selection_boundary": "selected_by_primary_lodo_only",
                "n_samples": int(len(pred)),
                "n_responders": int(pred["true_response_label"].sum()),
                "n_nonresponders": int((pred["true_response_label"] == 0).sum()),
                **metrics,
            }
        )
        prediction_rows.extend(parts)
    return pd.DataFrame(rows).sort_values(["AUROC", "AUPRC", "balanced_accuracy"], ascending=False), pd.DataFrame(prediction_rows)


def evaluate_external(
    transforms: dict[str, dict[str, pd.DataFrame]],
    primary_y: dict[str, pd.Series],
    strict_y: dict[str, pd.Series],
    specs: list[dict[str, object]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    train = PRIMARY_COHORTS
    y_train = pd.concat([primary_y[cohort] for cohort in train]).astype(int)
    for spec in specs:
        p_train = pd.concat(
            [score_combo(transforms, cohort, spec, train, primary_y).reindex(primary_y[cohort].index) for cohort in train]
        ).reindex(y_train.index)
        threshold = select_threshold(y_train.to_numpy(dtype=int), p_train.to_numpy(dtype=float))
        parts: list[dict[str, object]] = []
        for cohort in STRICT_EXTERNAL_COHORTS:
            if cohort not in strict_y:
                continue
            y_test = strict_y[cohort].astype(int)
            p_test = score_combo(transforms, cohort, spec, train, primary_y).reindex(y_test.index)
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
                        "candidate": spec["candidate"],
                        "weight_base": spec["weight_base"],
                        "component_method": spec.get("component_method", ""),
                        "component_gene": spec.get("component_gene", ""),
                        "selection_boundary": "external_locked_scoring_after_primary_selection",
                    }
                )
        pred = pd.DataFrame(parts)
        if pred.empty or pred["true_response_label"].nunique() < 2:
            continue
        metrics = compute_binary_metrics(pred["true_response_label"].astype(int), pred["response_probability"].astype(float), threshold=float(threshold))
        rows.append(
            {
                "endpoint": "strict_recist",
                "stratum": "strict_melanoma_pd1_like_external",
                "candidate": spec["candidate"],
                "weight_base": spec["weight_base"],
                "component_method": spec.get("component_method", ""),
                "component_gene": spec.get("component_gene", ""),
                "selection_boundary": "external_locked_scoring_after_primary_selection",
                "n_samples": int(len(pred)),
                "n_responders": int(pred["true_response_label"].sum()),
                "n_nonresponders": int((pred["true_response_label"] == 0).sum()),
                **metrics,
            }
        )
        prediction_rows.extend(parts)
    return pd.DataFrame(rows).sort_values(["AUROC", "AUPRC", "balanced_accuracy"], ascending=False), pd.DataFrame(prediction_rows)


def baseline_predictions(
    X_by_cohort: dict[str, pd.DataFrame],
    primary_y: dict[str, pd.Series],
    strict_y: dict[str, pd.Series],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    train = PRIMARY_COHORTS
    for model_name in EIGHT_SIGNATURES:
        genes = BASELINE_SIGNATURES.get(model_name, [model_name])
        train_scores = []
        train_labels = []
        for cohort in train:
            score = signature_score(X_by_cohort[cohort], genes)
            train_scores.append(score.reindex(primary_y[cohort].index))
            train_labels.append(primary_y[cohort])
        p_train = pd.concat(train_scores)
        y_train = pd.concat(train_labels).reindex(p_train.index).astype(int)
        if p_train.notna().sum() < 8:
            continue
        score_mean = float(p_train.mean())
        score_std = float(p_train.std(ddof=0) + 1e-6)
        p_train_z = (p_train - score_mean) / score_std
        threshold = select_threshold(y_train.to_numpy(dtype=int), p_train_z.to_numpy(dtype=float))
        for cohort in STRICT_EXTERNAL_COHORTS:
            if cohort not in strict_y:
                continue
            raw = signature_score(X_by_cohort[cohort], genes).reindex(strict_y[cohort].index)
            score = (raw - float(raw.mean())) / float(raw.std(ddof=0) + 1e-6)
            for sample_id, probability in score.items():
                rows.append(
                    {
                        "endpoint": "strict_recist",
                        "stratum": "strict_melanoma_pd1_like_external",
                        "cohort": cohort,
                        "sample_id": sample_id,
                        "true_response_label": int(strict_y[cohort].loc[sample_id]),
                        "response_probability": float(probability),
                        "threshold": float(threshold),
                        "model_name": model_name,
                    }
                )
    return pd.DataFrame(rows)


def _bootstrap_delta(y: pd.Series, target: pd.Series, baselines: pd.DataFrame, n_bootstrap: int = 2000) -> dict[str, float]:
    y_arr = y.to_numpy(dtype=int)
    target_arr = target.to_numpy(dtype=float)
    baseline_arr = baselines.to_numpy(dtype=float)
    rng = np.random.default_rng(20260527)
    deltas = []
    for _ in range(n_bootstrap):
        idx = rng.integers(0, len(y_arr), len(y_arr))
        if len(np.unique(y_arr[idx])) < 2:
            continue
        target_auc = roc_auc_score(y_arr[idx], target_arr[idx])
        baseline_auc = np.mean([roc_auc_score(y_arr[idx], baseline_arr[idx, col]) for col in range(baseline_arr.shape[1])])
        deltas.append(float(target_auc - baseline_auc))
    arr = np.asarray(deltas, dtype=float)
    return {
        "delta_vs_family_mean": float(arr.mean()),
        "ci_low": float(np.quantile(arr, 0.025)),
        "ci_high": float(np.quantile(arr, 0.975)),
        "one_sided_p": float((arr <= 0.0).mean()),
        "two_sided_p": float(min(1.0, 2.0 * min((arr <= 0.0).mean(), (arr >= 0.0).mean()))),
    }


def family_comparison(selected_external_predictions: pd.DataFrame, baseline_external_predictions: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    target_key = selected_external_predictions["cohort"].astype(str) + "::" + selected_external_predictions["sample_id"].astype(str)
    target = pd.Series(selected_external_predictions["response_probability"].to_numpy(dtype=float), index=target_key)
    y = pd.Series(selected_external_predictions["true_response_label"].to_numpy(dtype=int), index=target_key)
    baseline_series: dict[str, pd.Series] = {}
    for model_name, frame in baseline_external_predictions.groupby("model_name"):
        key = frame["cohort"].astype(str) + "::" + frame["sample_id"].astype(str)
        baseline_series[str(model_name)] = pd.Series(frame["response_probability"].to_numpy(dtype=float), index=key)
    baseline_frame = pd.DataFrame(baseline_series).dropna(axis=0)
    common = baseline_frame.index.intersection(target.index).intersection(y.index)
    if len(common) >= 8 and y.loc[common].nunique() == 2 and baseline_frame.shape[1] >= 4:
        target_auc = float(roc_auc_score(y.loc[common], target.loc[common]))
        baseline_aucs = {name: float(roc_auc_score(y.loc[common], baseline_frame.loc[common, name])) for name in baseline_frame.columns}
        stats = _bootstrap_delta(y.loc[common], target.loc[common], baseline_frame.loc[common])
        rows.append(
            {
                "endpoint": "strict_recist",
                "stratum": "strict_melanoma_pd1_like_external",
                "target_model": str(selected_external_predictions["candidate"].iloc[0]),
                "baseline_family": "eight_strong_signatures",
                "n_samples": int(len(common)),
                "n_signatures": int(baseline_frame.shape[1]),
                "target_AUROC": target_auc,
                "family_mean_AUROC": float(np.mean(list(baseline_aucs.values()))),
                "best_signature": max(baseline_aucs, key=baseline_aucs.get),
                "best_signature_AUROC": float(max(baseline_aucs.values())),
                **stats,
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["two_sided_fdr_q"] = benjamini_hochberg(result["two_sided_p"].fillna(1.0))
        result["claim_level"] = np.where(
            (result["target_AUROC"] >= 0.70) & (result["delta_vs_family_mean"] > 0.0) & (result["two_sided_fdr_q"] <= 0.05),
            "strict_external_family_FDR_supported_numeric_target_met",
            np.where(result["delta_vs_family_mean"] > 0.0, "family_point_estimate_only", "family_not_superior"),
        )
    return result


def build_selection(primary: pd.DataFrame, external: pd.DataFrame, family: pd.DataFrame) -> pd.DataFrame:
    selected = primary.sort_values(["AUROC", "AUPRC", "balanced_accuracy"], ascending=False).iloc[0]
    external_row = external[external["candidate"].astype(str).eq(str(selected["candidate"]))].iloc[0]
    family_row = family.iloc[0] if not family.empty else pd.Series(dtype=object)
    return pd.DataFrame(
        [
            {
                "selection_id": "primary_lodo_selected_rescue_combo",
                "candidate": selected["candidate"],
                "component_gene": selected["component_gene"],
                "component_method": selected["component_method"],
                "weight_base": float(selected["weight_base"]),
                "selection_boundary": "gene_screen_and_combo_weight_selected_by_primary_lodo_only",
                "primary_AUROC": float(selected["AUROC"]),
                "primary_AUPRC": float(selected["AUPRC"]),
                "primary_balanced_accuracy": float(selected["balanced_accuracy"]),
                "strict_external_AUROC": float(external_row["AUROC"]),
                "strict_external_AUPRC": float(external_row["AUPRC"]),
                "strict_external_balanced_accuracy": float(external_row["balanced_accuracy"]),
                "strict_external_ECE": float(external_row["ECE"]),
                "family_mean_AUROC": float(family_row.get("family_mean_AUROC", np.nan)),
                "delta_vs_family_mean": float(family_row.get("delta_vs_family_mean", np.nan)),
                "two_sided_fdr_q": float(family_row.get("two_sided_fdr_q", np.nan)),
                "claim_level": str(family_row.get("claim_level", "family_comparison_missing")),
            }
        ]
    )


def write_markdown(selection: pd.DataFrame, out_md: Path) -> None:
    row = selection.iloc[0]
    lines = [
        "# Discovery-only Rescue Combo Search",
        "",
        "This audit screens genes and rescue-head combinations using primary melanoma LODO only. Strict external labels are used only after the primary-selected candidate is locked.",
        "",
        "- Selected candidate: `{}`.".format(row["candidate"]),
        "- Primary LODO: AUROC={:.3f}, AUPRC={:.3f}, balanced accuracy={:.3f}.".format(
            float(row["primary_AUROC"]),
            float(row["primary_AUPRC"]),
            float(row["primary_balanced_accuracy"]),
        ),
        "- Strict external GSE145996+PHS000452: AUROC={:.3f}, AUPRC={:.3f}, balanced accuracy={:.3f}, ECE={:.3f}.".format(
            float(row["strict_external_AUROC"]),
            float(row["strict_external_AUPRC"]),
            float(row["strict_external_balanced_accuracy"]),
            float(row["strict_external_ECE"]),
        ),
        "- Eight-signature family comparison: family mean AUROC={:.3f}, delta={:.3f}, q={:.3f}, claim={}.".format(
            float(row["family_mean_AUROC"]),
            float(row["delta_vs_family_mean"]),
            float(row["two_sided_fdr_q"]),
            row["claim_level"],
        ),
        "",
        "Claim boundary: this is a no-external-label feature-selection/thresholding audit. Because it was added after prior external failure-mode work, it should be frozen as the next locked melanoma rescue-combo candidate and confirmed on any newly obtained independent controlled cohort.",
    ]
    out_md.write_text("\n".join(lines), encoding="utf-8")


def run_search(processed_dir: Path, out_dir: Path, top_k: int) -> dict[str, Path]:
    _BASE_SCORE_CACHE.clear()
    _DIRECTION_CACHE.clear()
    _COMPONENT_CACHE.clear()
    X_by_cohort, _, metadata_by_cohort = load_processed_bulk(processed_dir)
    primary_y = labels_for_endpoint(X_by_cohort, metadata_by_cohort, PRIMARY_COHORTS, "primary_recist")
    strict_y = labels_for_endpoint(X_by_cohort, metadata_by_cohort, [*PRIMARY_COHORTS, *STRICT_EXTERNAL_COHORTS], "strict_recist")
    genes = primary_gene_universe(X_by_cohort)
    transforms = build_transforms(X_by_cohort, primary_y, strict_y, genes, [*PRIMARY_COHORTS, *STRICT_EXTERNAL_COHORTS])
    gene_summary = gene_screen(transforms, primary_y, genes)
    specs = combo_specs(gene_summary, top_k=top_k)
    primary_summary, primary_predictions = evaluate_primary_lodo(transforms, primary_y, specs)
    primary_selected = primary_summary.sort_values(["AUROC", "AUPRC", "balanced_accuracy"], ascending=False).head(1)
    selected_specs = [spec for spec in specs if str(spec["candidate"]) == str(primary_selected.iloc[0]["candidate"])]
    external_summary, external_predictions = evaluate_external(transforms, primary_y, strict_y, selected_specs)
    baseline_external = baseline_predictions(X_by_cohort, primary_y, strict_y)
    family = family_comparison(external_predictions, baseline_external)
    selection = build_selection(primary_summary, external_summary, family)
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "gene_summary": out_dir / "discovery_only_gene_screen_summary.tsv",
        "combo_primary_summary": out_dir / "discovery_only_combo_primary_summary.tsv",
        "primary_predictions": out_dir / "discovery_only_combo_primary_predictions.tsv",
        "external_summary": out_dir / "discovery_only_combo_external_summary.tsv",
        "external_predictions": out_dir / "discovery_only_combo_external_predictions.tsv",
        "baseline_external_predictions": out_dir / "discovery_only_combo_external_baseline_predictions.tsv",
        "family_comparison": out_dir / "discovery_only_combo_external_family_comparison.tsv",
        "selection": out_dir / "discovery_only_rescue_combo_selection.tsv",
        "markdown": out_dir / "DISCOVERY_ONLY_RESCUE_COMBO_AUDIT.md",
    }
    gene_summary.to_csv(outputs["gene_summary"], sep="\t", index=False)
    primary_summary.to_csv(outputs["combo_primary_summary"], sep="\t", index=False)
    primary_predictions.to_csv(outputs["primary_predictions"], sep="\t", index=False)
    external_summary.to_csv(outputs["external_summary"], sep="\t", index=False)
    external_predictions.to_csv(outputs["external_predictions"], sep="\t", index=False)
    baseline_external.to_csv(outputs["baseline_external_predictions"], sep="\t", index=False)
    family.to_csv(outputs["family_comparison"], sep="\t", index=False)
    selection.to_csv(outputs["selection"], sep="\t", index=False)
    write_markdown(selection, outputs["markdown"])
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--out", default="results/discovery_only_rescue_combo_search_20260527")
    parser.add_argument("--top-k", type=int, default=350)
    args = parser.parse_args()
    outputs = run_search(ROOT / args.processed_dir, ROOT / args.out, top_k=args.top_k)
    selection = pd.read_csv(outputs["selection"], sep="\t")
    print(json.dumps(selection.to_dict("records"), ensure_ascii=False))
    print(f"Wrote {ROOT / args.out}")


if __name__ == "__main__":
    main()
