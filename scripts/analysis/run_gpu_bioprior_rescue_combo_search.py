from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import average_precision_score, balanced_accuracy_score, roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.baselines import BASELINE_SIGNATURES, signature_score
from econiche.metrics import compute_binary_metrics
from econiche.statistics import benjamini_hochberg
from econiche_opt.model.endpoint_modules import endpoint_label_series, select_threshold


PRIMARY_COHORTS = ["GSE91061", "GSE78220", "PRJEB23709_PD1_PRE"]
STRICT_EXTERNAL_COHORTS = ["GSE145996", "PHS000452_LIU_LIKE_PRE"]
BASE_GENES = ["MAP4K1", "TBX3", "AXL"]
LIPID_PI3K_PRIOR_GENES = ["PLA2G2D", "PIK3CD"]
IMMUNE_PRIOR_GENES = ["PLA2G2D", "PIK3CD", "ICOS", "LCK", "CD247", "TIGIT", "SLAMF7", "MAP4K1"]
TRANSFORMS = ["rz", "z", "pct"]
WEIGHTS = [0.20, 0.35, 0.50, 0.65, 0.80]
EIGHT_SIGNATURES = ["IFNG", "CXCL9", "TIG", "TIDE_dysfunction", "TIDE_exclusion", "CYT", "APM", "IPRES"]


def labels_for_endpoint(X_by_cohort: dict[str, pd.DataFrame], metadata_by_cohort: dict[str, pd.DataFrame], cohorts: list[str], endpoint: str) -> dict[str, pd.Series]:
    out: dict[str, pd.Series] = {}
    for cohort in cohorts:
        y = endpoint_label_series(metadata_by_cohort[cohort]["response_raw"], endpoint)
        common = X_by_cohort[cohort].index.intersection(y[y.notna()].index)
        if len(common) >= 8 and y.loc[common].nunique() == 2:
            out[cohort] = y.loc[common].astype(int)
    return out


def load_bulk(processed_dir: Path) -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame]]:
    sys.path.insert(0, str(SRC))
    from econiche.io import load_processed_bulk

    X_by_cohort, _, metadata_by_cohort = load_processed_bulk(processed_dir)
    return X_by_cohort, metadata_by_cohort


def _needed_samples(primary_y: dict[str, pd.Series], strict_y: dict[str, pd.Series], cohort: str) -> list[str]:
    samples: list[str] = []
    for labels in [primary_y.get(cohort), strict_y.get(cohort)]:
        if labels is not None:
            samples.extend(labels.index.astype(str).tolist())
    return list(dict.fromkeys(samples))


def build_transforms(
    X_by_cohort: dict[str, pd.DataFrame],
    primary_y: dict[str, pd.Series],
    strict_y: dict[str, pd.Series],
    genes: list[str],
) -> dict[str, dict[str, pd.DataFrame]]:
    transforms: dict[str, dict[str, pd.DataFrame]] = {}
    for cohort in [*PRIMARY_COHORTS, *STRICT_EXTERNAL_COHORTS]:
        samples = _needed_samples(primary_y, strict_y, cohort)
        values = X_by_cohort[cohort].loc[samples, genes].apply(pd.to_numeric, errors="coerce")
        z = ((values - values.mean()) / (values.std(ddof=0) + 1e-6)).fillna(0.0)
        median = values.median()
        mad = (values - median).abs().median() + 1e-6
        rz = ((values - median) / (1.4826 * mad)).clip(-5.0, 5.0).fillna(0.0)
        pct = values.rank(axis=0, pct=True).fillna(0.5)
        transforms[cohort] = {"z": z.astype("float32"), "rz": rz.astype("float32"), "pct": pct.astype("float32")}
    return transforms


def _normalize_tensor(x: torch.Tensor) -> torch.Tensor:
    return (x - x.min(dim=0, keepdim=True).values) / (x.max(dim=0, keepdim=True).values - x.min(dim=0, keepdim=True).values + 1e-9)


def candidate_specs(prior: str, transform_policy: str = "all") -> list[dict[str, object]]:
    genes = LIPID_PI3K_PRIOR_GENES if prior == "lipid_pi3k" else IMMUNE_PRIOR_GENES
    transforms = ["rz"] if transform_policy == "robust_only" else TRANSFORMS
    specs = [{"candidate": "base_rescue_robust", "weight_base": 1.0, "method": "", "gene": "", "prior": prior}]
    for method in transforms:
        for gene in genes:
            for weight in WEIGHTS:
                specs.append(
                    {
                        "candidate": f"{weight:.2f}*base+{1.0 - weight:.2f}*{method}__{gene}",
                        "weight_base": weight,
                        "method": method,
                        "gene": gene,
                        "prior": prior,
                        "transform_policy": transform_policy,
                    }
                )
    return specs


def score_matrix(
    transforms: dict[str, dict[str, pd.DataFrame]],
    cohort: str,
    specs: list[dict[str, object]],
    train_cohorts: list[str],
    train_y: dict[str, pd.Series],
    device: torch.device,
) -> pd.DataFrame:
    rz = torch.tensor(transforms[cohort]["rz"][BASE_GENES].to_numpy(dtype=np.float32), device=device)
    z = torch.tensor(transforms[cohort]["z"][BASE_GENES].to_numpy(dtype=np.float32), device=device)
    base_rz = _normalize_tensor((rz[:, [0]] - rz[:, [1, 2]].mean(dim=1, keepdim=True)))
    base_z = _normalize_tensor((z[:, [0]] - z[:, [1, 2]].mean(dim=1, keepdim=True)))
    base = _normalize_tensor(0.95 * base_rz + 0.05 * base_z)
    columns = []
    parts = []
    direction_cache: dict[tuple[str, str], float] = {}
    for spec in specs:
        candidate = str(spec["candidate"])
        columns.append(candidate)
        method = str(spec["method"])
        gene = str(spec["gene"])
        if not method or not gene:
            parts.append(base)
            continue
        key = (method, gene)
        if key not in direction_cache:
            values = []
            labels = []
            for train_cohort in train_cohorts:
                values.append(transforms[train_cohort][method][gene].reindex(train_y[train_cohort].index))
                labels.append(train_y[train_cohort])
            v = pd.concat(values)
            y = pd.concat(labels).reindex(v.index).astype(int)
            direction_cache[key] = 1.0 if float(v[y == 1].mean() - v[y == 0].mean()) >= 0.0 else -1.0
        component = torch.tensor(transforms[cohort][method][[gene]].to_numpy(dtype=np.float32), device=device) * direction_cache[key]
        component = _normalize_tensor(component)
        score = _normalize_tensor(float(spec["weight_base"]) * base + (1.0 - float(spec["weight_base"])) * component)
        parts.append(score)
    matrix = torch.cat(parts, dim=1).detach().cpu().numpy()
    return pd.DataFrame(matrix, index=transforms[cohort]["rz"].index, columns=columns)


def primary_lodo(
    transforms: dict[str, dict[str, pd.DataFrame]],
    primary_y: dict[str, pd.Series],
    specs: list[dict[str, object]],
    device: torch.device,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pred_rows: list[pd.DataFrame] = []
    for holdout in PRIMARY_COHORTS:
        train = [cohort for cohort in PRIMARY_COHORTS if cohort != holdout]
        train_scores = []
        train_labels = []
        for cohort in train:
            mat = score_matrix(transforms, cohort, specs, train, primary_y, device).reindex(primary_y[cohort].index)
            train_scores.append(mat)
            train_labels.append(primary_y[cohort])
        X_train = pd.concat(train_scores)
        y_train = pd.concat(train_labels).reindex(X_train.index).astype(int)
        thresholds = {
            column: select_threshold(y_train.to_numpy(dtype=int), X_train[column].to_numpy(dtype=float))
            for column in X_train.columns
        }
        test = score_matrix(transforms, holdout, specs, train, primary_y, device).reindex(primary_y[holdout].index)
        for candidate in test.columns:
            frame = pd.DataFrame(
                {
                    "endpoint": "primary_recist",
                    "stratum": "melanoma_core_high_evidence",
                    "cohort": holdout,
                    "sample_id": test.index,
                    "true_response_label": primary_y[holdout].reindex(test.index).astype(int).to_numpy(),
                    "response_probability": test[candidate].to_numpy(dtype=float),
                    "threshold": float(thresholds[candidate]),
                    "candidate": candidate,
                }
            )
            pred_rows.append(frame)
    predictions = pd.concat(pred_rows, ignore_index=True)
    rows = []
    for candidate, frame in predictions.groupby("candidate"):
        y = frame["true_response_label"].astype(int)
        p = frame["response_probability"].astype(float)
        pred = p >= frame["threshold"].astype(float)
        rows.append(
            {
                "candidate": candidate,
                "AUROC": float(roc_auc_score(y, p)),
                "AUPRC": float(average_precision_score(y, p)),
                "balanced_accuracy": float(balanced_accuracy_score(y, pred.astype(int))),
                "n_samples": int(len(frame)),
                "n_responders": int(y.sum()),
                "n_nonresponders": int((y == 0).sum()),
            }
        )
    return pd.DataFrame(rows).sort_values(["AUROC", "AUPRC", "balanced_accuracy"], ascending=False), predictions


def external_score(
    transforms: dict[str, dict[str, pd.DataFrame]],
    primary_y: dict[str, pd.Series],
    strict_y: dict[str, pd.Series],
    specs: list[dict[str, object]],
    device: torch.device,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    train = PRIMARY_COHORTS
    train_scores = []
    train_labels = []
    for cohort in train:
        mat = score_matrix(transforms, cohort, specs, train, primary_y, device).reindex(primary_y[cohort].index)
        train_scores.append(mat)
        train_labels.append(primary_y[cohort])
    X_train = pd.concat(train_scores)
    y_train = pd.concat(train_labels).reindex(X_train.index).astype(int)
    thresholds = {
        column: select_threshold(y_train.to_numpy(dtype=int), X_train[column].to_numpy(dtype=float))
        for column in X_train.columns
    }
    frames = []
    for cohort in STRICT_EXTERNAL_COHORTS:
        test = score_matrix(transforms, cohort, specs, train, primary_y, device).reindex(strict_y[cohort].index)
        for candidate in test.columns:
            frames.append(
                pd.DataFrame(
                    {
                        "endpoint": "strict_recist",
                        "stratum": "strict_melanoma_pd1_like_external",
                        "cohort": cohort,
                        "sample_id": test.index,
                        "true_response_label": strict_y[cohort].reindex(test.index).astype(int).to_numpy(),
                        "response_probability": test[candidate].to_numpy(dtype=float),
                        "threshold": float(thresholds[candidate]),
                        "candidate": candidate,
                    }
                )
            )
    predictions = pd.concat(frames, ignore_index=True)
    rows = []
    for candidate, frame in predictions.groupby("candidate"):
        y = frame["true_response_label"].astype(int)
        p = frame["response_probability"].astype(float)
        metrics = compute_binary_metrics(y, p, threshold=float(frame["threshold"].iloc[0]))
        rows.append(
            {
                "candidate": candidate,
                "n_samples": int(len(frame)),
                "n_responders": int(y.sum()),
                "n_nonresponders": int((y == 0).sum()),
                **metrics,
            }
        )
    return pd.DataFrame(rows).sort_values(["AUROC", "AUPRC", "balanced_accuracy"], ascending=False), predictions


def baseline_external_predictions(X_by_cohort: dict[str, pd.DataFrame], primary_y: dict[str, pd.Series], strict_y: dict[str, pd.Series]) -> pd.DataFrame:
    rows = []
    for model_name in EIGHT_SIGNATURES:
        genes = BASELINE_SIGNATURES.get(model_name, [model_name])
        train_scores = []
        train_labels = []
        for cohort in PRIMARY_COHORTS:
            score = signature_score(X_by_cohort[cohort], genes).reindex(primary_y[cohort].index)
            train_scores.append(score)
            train_labels.append(primary_y[cohort])
        train = pd.concat(train_scores)
        y_train = pd.concat(train_labels).reindex(train.index).astype(int)
        train = (train - train.mean()) / (train.std(ddof=0) + 1e-6)
        threshold = select_threshold(y_train.to_numpy(dtype=int), train.to_numpy(dtype=float))
        for cohort in STRICT_EXTERNAL_COHORTS:
            raw = signature_score(X_by_cohort[cohort], genes).reindex(strict_y[cohort].index)
            score = (raw - raw.mean()) / (raw.std(ddof=0) + 1e-6)
            for sample_id, value in score.items():
                rows.append(
                    {
                        "cohort": cohort,
                        "sample_id": sample_id,
                        "model_name": model_name,
                        "true_response_label": int(strict_y[cohort].loc[sample_id]),
                        "response_probability": float(value),
                        "threshold": float(threshold),
                    }
                )
    return pd.DataFrame(rows)


def family_gate(target_predictions: pd.DataFrame, baseline_predictions: pd.DataFrame, n_bootstrap: int) -> pd.DataFrame:
    key = target_predictions["cohort"].astype(str) + "::" + target_predictions["sample_id"].astype(str)
    y = pd.Series(target_predictions["true_response_label"].to_numpy(dtype=int), index=key)
    target = pd.Series(target_predictions["response_probability"].to_numpy(dtype=float), index=key)
    baseline_series = {}
    for model_name, frame in baseline_predictions.groupby("model_name"):
        bkey = frame["cohort"].astype(str) + "::" + frame["sample_id"].astype(str)
        baseline_series[model_name] = pd.Series(frame["response_probability"].to_numpy(dtype=float), index=bkey)
    baselines = pd.DataFrame(baseline_series).dropna()
    common = y.index.intersection(target.index).intersection(baselines.index)
    y = y.loc[common]
    target = target.loc[common]
    baselines = baselines.loc[common]
    rng = np.random.default_rng(20260527)
    deltas = []
    for _ in range(n_bootstrap):
        idx = rng.integers(0, len(y), len(y))
        yy = y.to_numpy(dtype=int)[idx]
        if len(np.unique(yy)) < 2:
            continue
        tt = target.to_numpy(dtype=float)[idx]
        bb = baselines.to_numpy(dtype=float)[idx]
        deltas.append(float(roc_auc_score(yy, tt) - np.mean([roc_auc_score(yy, bb[:, col]) for col in range(bb.shape[1])])))
    arr = np.asarray(deltas)
    baseline_aucs = {col: float(roc_auc_score(y, baselines[col])) for col in baselines.columns}
    result = pd.DataFrame(
        [
            {
                "target_model": str(target_predictions["candidate"].iloc[0]),
                "baseline_family": "eight_strong_signatures",
                "n_samples": int(len(common)),
                "n_signatures": int(baselines.shape[1]),
                "target_AUROC": float(roc_auc_score(y, target)),
                "family_mean_AUROC": float(np.mean(list(baseline_aucs.values()))),
                "best_signature": max(baseline_aucs, key=baseline_aucs.get),
                "best_signature_AUROC": float(max(baseline_aucs.values())),
                "delta_vs_family_mean": float(arr.mean()),
                "ci_low": float(np.quantile(arr, 0.025)),
                "ci_high": float(np.quantile(arr, 0.975)),
                "one_sided_p": float((arr <= 0).mean()),
                "two_sided_p": float(min(1.0, 2.0 * min((arr <= 0).mean(), (arr >= 0).mean()))),
            }
        ]
    )
    result["two_sided_fdr_q"] = benjamini_hochberg(result["two_sided_p"])
    result["claim_level"] = np.where(
        (result["target_AUROC"] >= 0.70) & (result["delta_vs_family_mean"] > 0) & (result["two_sided_fdr_q"] <= 0.05),
        "strict_external_family_FDR_supported_numeric_target_met",
        np.where(result["delta_vs_family_mean"] > 0, "family_point_estimate_only", "family_not_superior"),
    )
    return result


def run_search(processed_dir: Path, out_dir: Path, prior: str, transform_policy: str, n_bootstrap: int) -> dict[str, Path]:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    X_by_cohort, metadata_by_cohort = load_bulk(processed_dir)
    primary_y = labels_for_endpoint(X_by_cohort, metadata_by_cohort, PRIMARY_COHORTS, "primary_recist")
    strict_y = labels_for_endpoint(X_by_cohort, metadata_by_cohort, [*PRIMARY_COHORTS, *STRICT_EXTERNAL_COHORTS], "strict_recist")
    genes = sorted(set(BASE_GENES + LIPID_PI3K_PRIOR_GENES + IMMUNE_PRIOR_GENES))
    transforms = build_transforms(X_by_cohort, primary_y, strict_y, genes)
    specs = candidate_specs(prior, transform_policy=transform_policy)
    primary, primary_predictions = primary_lodo(transforms, primary_y, specs, device)
    selected = primary.sort_values(["AUROC", "AUPRC", "balanced_accuracy"], ascending=False).iloc[0]
    selected_specs = [spec for spec in specs if str(spec["candidate"]) == str(selected["candidate"])]
    external, external_predictions = external_score(transforms, primary_y, strict_y, selected_specs, device)
    baseline_predictions = baseline_external_predictions(X_by_cohort, primary_y, strict_y)
    gate = family_gate(external_predictions, baseline_predictions, n_bootstrap=n_bootstrap)
    selection = pd.DataFrame(
        [
            {
                "selection_id": f"gpu_{prior}_primary_selected_rescue_combo",
                "candidate": selected["candidate"],
                "prior": prior,
                "transform_policy": transform_policy,
                "device": str(device),
                "gpu_name": torch.cuda.get_device_name(0) if torch.cuda.is_available() else "cpu",
                "selection_boundary": "candidate_and_weight_selected_by_primary_lodo_only_with_biological_prior",
                "primary_AUROC": float(selected["AUROC"]),
                "primary_AUPRC": float(selected["AUPRC"]),
                "primary_balanced_accuracy": float(selected["balanced_accuracy"]),
                "strict_external_AUROC": float(external.iloc[0]["AUROC"]),
                "strict_external_AUPRC": float(external.iloc[0]["AUPRC"]),
                "strict_external_balanced_accuracy": float(external.iloc[0]["balanced_accuracy"]),
                "strict_external_ECE": float(external.iloc[0]["ECE"]),
                "family_mean_AUROC": float(gate.iloc[0]["family_mean_AUROC"]),
                "delta_vs_family_mean": float(gate.iloc[0]["delta_vs_family_mean"]),
                "two_sided_fdr_q": float(gate.iloc[0]["two_sided_fdr_q"]),
                "claim_level": str(gate.iloc[0]["claim_level"]),
            }
        ]
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "primary_summary": out_dir / "gpu_bioprior_primary_summary.tsv",
        "primary_predictions": out_dir / "gpu_bioprior_primary_predictions.tsv",
        "external_summary": out_dir / "gpu_bioprior_external_summary.tsv",
        "external_predictions": out_dir / "gpu_bioprior_external_predictions.tsv",
        "baseline_predictions": out_dir / "gpu_bioprior_external_baseline_predictions.tsv",
        "family_gate": out_dir / "gpu_bioprior_external_family_gate.tsv",
        "selection": out_dir / "gpu_bioprior_rescue_combo_selection.tsv",
        "markdown": out_dir / "GPU_BIOPRIOR_RESCUE_COMBO_AUDIT.md",
    }
    primary.to_csv(outputs["primary_summary"], sep="\t", index=False)
    primary_predictions.to_csv(outputs["primary_predictions"], sep="\t", index=False)
    external.to_csv(outputs["external_summary"], sep="\t", index=False)
    external_predictions.to_csv(outputs["external_predictions"], sep="\t", index=False)
    baseline_predictions.to_csv(outputs["baseline_predictions"], sep="\t", index=False)
    gate.to_csv(outputs["family_gate"], sep="\t", index=False)
    selection.to_csv(outputs["selection"], sep="\t", index=False)
    lines = [
        "# GPU biological-prior rescue-combo audit",
        "",
        f"Device: {selection.iloc[0]['device']} ({selection.iloc[0]['gpu_name']}).",
        f"Selected candidate: `{selection.iloc[0]['candidate']}` under `{prior}` prior and `{transform_policy}` transform policy.",
        "Selection used primary melanoma LODO only; strict external labels were used only for locked scoring.",
        "",
        "Primary AUROC={:.3f}, AUPRC={:.3f}, balanced accuracy={:.3f}.".format(
            float(selection.iloc[0]["primary_AUROC"]),
            float(selection.iloc[0]["primary_AUPRC"]),
            float(selection.iloc[0]["primary_balanced_accuracy"]),
        ),
        "Strict external AUROC={:.3f}, AUPRC={:.3f}, balanced accuracy={:.3f}, ECE={:.3f}.".format(
            float(selection.iloc[0]["strict_external_AUROC"]),
            float(selection.iloc[0]["strict_external_AUPRC"]),
            float(selection.iloc[0]["strict_external_balanced_accuracy"]),
            float(selection.iloc[0]["strict_external_ECE"]),
        ),
        "Family gate: mean AUROC={:.3f}, delta={:.3f}, q={:.3f}, claim={}.".format(
            float(selection.iloc[0]["family_mean_AUROC"]),
            float(selection.iloc[0]["delta_vs_family_mean"]),
            float(selection.iloc[0]["two_sided_fdr_q"]),
            selection.iloc[0]["claim_level"],
        ),
    ]
    outputs["markdown"].write_text("\n".join(lines), encoding="utf-8")
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--out", default="results/gpu_bioprior_rescue_combo_search_20260527")
    parser.add_argument("--prior", choices=["lipid_pi3k", "immune"], default="lipid_pi3k")
    parser.add_argument("--transform-policy", choices=["all", "robust_only"], default="all")
    parser.add_argument("--bootstrap", type=int, default=2000)
    args = parser.parse_args()
    outputs = run_search(ROOT / args.processed_dir, ROOT / args.out, args.prior, args.transform_policy, args.bootstrap)
    selection = pd.read_csv(outputs["selection"], sep="\t")
    print(json.dumps(selection.to_dict("records"), ensure_ascii=False))
    print(f"Wrote {ROOT / args.out}")


if __name__ == "__main__":
    main()
