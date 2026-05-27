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
CBIO_COHORT = "CBIO_LIU_DFCI_2019_PRE"
BASE_GENES = ["MAP4K1", "TBX3", "AXL"]
PAIR_GENES = ["PLA2G2D", "PIK3CD"]
TRANSFORMS = ["rz", "z", "pct"]
BASE_WEIGHTS = [0.20, 0.35, 0.50, 0.65, 0.80]
PAIR_MIXES = [0.20, 0.35, 0.50, 0.65, 0.80]
EIGHT_SIGNATURES = ["IFNG", "CXCL9", "TIG", "TIDE_dysfunction", "TIDE_exclusion", "CYT", "APM", "IPRES"]


def _read_bulk(processed_dir: Path) -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame]]:
    from econiche.io import load_processed_bulk

    X_by_cohort, _, metadata_by_cohort = load_processed_bulk(processed_dir)
    return X_by_cohort, metadata_by_cohort


def _read_cbio(cbio_dir: Path, X_by_cohort: dict[str, pd.DataFrame], metadata_by_cohort: dict[str, pd.DataFrame]) -> None:
    X_by_cohort[CBIO_COHORT] = pd.read_csv(cbio_dir / f"{CBIO_COHORT}.expr.tsv", sep="\t", index_col=0)
    metadata_by_cohort[CBIO_COHORT] = pd.read_csv(cbio_dir / f"{CBIO_COHORT}.metadata.tsv", sep="\t")


def _labels(
    X_by_cohort: dict[str, pd.DataFrame],
    metadata_by_cohort: dict[str, pd.DataFrame],
    cohorts: list[str],
    endpoint: str,
) -> dict[str, pd.Series]:
    out: dict[str, pd.Series] = {}
    for cohort in cohorts:
        metadata = metadata_by_cohort[cohort]
        if "sample_id" in metadata.columns:
            metadata = metadata.set_index("sample_id")
        y = endpoint_label_series(metadata["response_raw"], endpoint).dropna().astype(int)
        common = X_by_cohort[cohort].index.intersection(y.index)
        if len(common) >= 8 and y.loc[common].nunique() == 2:
            out[cohort] = y.loc[common].astype(int)
    return out


def _needed_samples(primary_y: dict[str, pd.Series], strict_y: dict[str, pd.Series], cohort: str) -> list[str]:
    samples: list[str] = []
    for labels in [primary_y.get(cohort), strict_y.get(cohort)]:
        if labels is not None:
            samples.extend(labels.index.astype(str).tolist())
    return list(dict.fromkeys(samples))


def _build_transforms(
    X_by_cohort: dict[str, pd.DataFrame],
    primary_y: dict[str, pd.Series],
    strict_y: dict[str, pd.Series],
    genes: list[str],
) -> dict[str, dict[str, pd.DataFrame]]:
    transforms: dict[str, dict[str, pd.DataFrame]] = {}
    for cohort in [*PRIMARY_COHORTS, *STRICT_EXTERNAL_COHORTS, CBIO_COHORT]:
        samples = _needed_samples(primary_y, strict_y, cohort)
        values = X_by_cohort[cohort].loc[samples, genes].apply(pd.to_numeric, errors="coerce")
        z = ((values - values.mean()) / (values.std(ddof=0) + 1e-6)).fillna(0.0)
        median = values.median()
        mad = (values - median).abs().median() + 1e-6
        rz = ((values - median) / (1.4826 * mad)).clip(-5.0, 5.0).fillna(0.0)
        pct = values.rank(axis=0, pct=True).fillna(0.5)
        transforms[cohort] = {
            "z": z.astype("float32"),
            "rz": rz.astype("float32"),
            "pct": pct.astype("float32"),
        }
    return transforms


def _normalize_tensor(x: torch.Tensor) -> torch.Tensor:
    return (x - x.min(dim=0, keepdim=True).values) / (x.max(dim=0, keepdim=True).values - x.min(dim=0, keepdim=True).values + 1e-9)


def _candidate_specs(base_weight_policy: str) -> list[dict[str, object]]:
    specs: list[dict[str, object]] = []
    base_weights = [weight for weight in BASE_WEIGHTS if base_weight_policy == "all" or weight <= 0.35]
    for method1 in TRANSFORMS:
        for method2 in TRANSFORMS:
            for weight_base in base_weights:
                for mix in PAIR_MIXES:
                    specs.append(
                        {
                            "candidate": (
                                f"{weight_base:.2f}*base+{1.0 - weight_base:.2f}*"
                                f"({mix:.2f}*{method1}__PLA2G2D+{1.0 - mix:.2f}*{method2}__PIK3CD)"
                            ),
                            "method1": method1,
                            "method2": method2,
                            "weight_base": weight_base,
                            "mix": mix,
                        }
                    )
    return specs


def _direction(
    transforms: dict[str, dict[str, pd.DataFrame]],
    train_cohorts: list[str],
    train_y: dict[str, pd.Series],
    method: str,
    gene: str,
) -> float:
    values = []
    labels = []
    for cohort in train_cohorts:
        values.append(transforms[cohort][method][gene].reindex(train_y[cohort].index))
        labels.append(train_y[cohort])
    v = pd.concat(values)
    y = pd.concat(labels).reindex(v.index).astype(int)
    return 1.0 if float(v[y == 1].mean() - v[y == 0].mean()) >= 0.0 else -1.0


def _score_matrix(
    transforms: dict[str, dict[str, pd.DataFrame]],
    cohort: str,
    specs: list[dict[str, object]],
    train_cohorts: list[str],
    train_y: dict[str, pd.Series],
    device: torch.device,
) -> pd.DataFrame:
    rz = torch.tensor(transforms[cohort]["rz"][BASE_GENES].to_numpy(dtype=np.float32), device=device)
    z = torch.tensor(transforms[cohort]["z"][BASE_GENES].to_numpy(dtype=np.float32), device=device)
    base_rz = _normalize_tensor(rz[:, [0]] - rz[:, [1, 2]].mean(dim=1, keepdim=True))
    base_z = _normalize_tensor(z[:, [0]] - z[:, [1, 2]].mean(dim=1, keepdim=True))
    base = _normalize_tensor(0.95 * base_rz + 0.05 * base_z)
    direction_cache: dict[tuple[str, str], float] = {}
    columns: list[str] = []
    parts: list[torch.Tensor] = []
    for spec in specs:
        method1 = str(spec["method1"])
        method2 = str(spec["method2"])
        for method, gene in [(method1, "PLA2G2D"), (method2, "PIK3CD")]:
            key = (method, gene)
            if key not in direction_cache:
                direction_cache[key] = _direction(transforms, train_cohorts, train_y, method, gene)
        comp1 = torch.tensor(transforms[cohort][method1][["PLA2G2D"]].to_numpy(dtype=np.float32), device=device)
        comp2 = torch.tensor(transforms[cohort][method2][["PIK3CD"]].to_numpy(dtype=np.float32), device=device)
        comp1 = _normalize_tensor(direction_cache[(method1, "PLA2G2D")] * comp1)
        comp2 = _normalize_tensor(direction_cache[(method2, "PIK3CD")] * comp2)
        mix = float(spec["mix"])
        pair_component = _normalize_tensor(mix * comp1 + (1.0 - mix) * comp2)
        weight_base = float(spec["weight_base"])
        score = _normalize_tensor(weight_base * base + (1.0 - weight_base) * pair_component)
        columns.append(str(spec["candidate"]))
        parts.append(score)
    matrix = torch.cat(parts, dim=1).detach().cpu().numpy()
    return pd.DataFrame(matrix, index=transforms[cohort]["rz"].index, columns=columns)


def _primary_lodo(
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
            train_scores.append(_score_matrix(transforms, cohort, specs, train, primary_y, device).reindex(primary_y[cohort].index))
            train_labels.append(primary_y[cohort])
        X_train = pd.concat(train_scores)
        y_train = pd.concat(train_labels).reindex(X_train.index).astype(int)
        thresholds = {
            column: select_threshold(y_train.to_numpy(dtype=int), X_train[column].to_numpy(dtype=float))
            for column in X_train.columns
        }
        test = _score_matrix(transforms, holdout, specs, train, primary_y, device).reindex(primary_y[holdout].index)
        for candidate in test.columns:
            pred_rows.append(
                pd.DataFrame(
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
            )
    predictions = pd.concat(pred_rows, ignore_index=True)
    rows = []
    for candidate, frame in predictions.groupby("candidate"):
        y = frame["true_response_label"].astype(int)
        p = frame["response_probability"].astype(float)
        pred = p >= frame["threshold"].astype(float)
        fold_aucs = [roc_auc_score(g["true_response_label"].astype(int), g["response_probability"].astype(float)) for _, g in frame.groupby("cohort")]
        rows.append(
            {
                "candidate": candidate,
                "AUROC": float(roc_auc_score(y, p)),
                "AUPRC": float(average_precision_score(y, p)),
                "balanced_accuracy": float(balanced_accuracy_score(y, pred.astype(int))),
                "mean_fold_AUROC": float(np.mean(fold_aucs)),
                "min_fold_AUROC": float(np.min(fold_aucs)),
                "n_samples": int(len(frame)),
                "n_responders": int(y.sum()),
                "n_nonresponders": int((y == 0).sum()),
            }
        )
    return pd.DataFrame(rows).sort_values(["AUROC", "AUPRC", "balanced_accuracy"], ascending=False), predictions


def _score_locked_groups(
    transforms: dict[str, dict[str, pd.DataFrame]],
    primary_y: dict[str, pd.Series],
    strict_y: dict[str, pd.Series],
    selected_spec: list[dict[str, object]],
    device: torch.device,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    train = PRIMARY_COHORTS
    train_scores = []
    train_labels = []
    for cohort in train:
        train_scores.append(_score_matrix(transforms, cohort, selected_spec, train, primary_y, device).reindex(primary_y[cohort].index))
        train_labels.append(primary_y[cohort])
    X_train = pd.concat(train_scores)
    y_train = pd.concat(train_labels).reindex(X_train.index).astype(int)
    threshold = select_threshold(y_train.to_numpy(dtype=int), X_train.iloc[:, 0].to_numpy(dtype=float))
    frames = []
    for cohort in [*STRICT_EXTERNAL_COHORTS, CBIO_COHORT]:
        labels = strict_y[cohort]
        test = _score_matrix(transforms, cohort, selected_spec, train, primary_y, device).reindex(labels.index)
        frames.append(
            pd.DataFrame(
                {
                    "endpoint": "strict_recist",
                    "cohort": cohort,
                    "sample_id": test.index,
                    "true_response_label": labels.reindex(test.index).astype(int).to_numpy(),
                    "response_probability": test.iloc[:, 0].to_numpy(dtype=float),
                    "threshold": float(threshold),
                    "candidate": str(selected_spec[0]["candidate"]),
                }
            )
        )
    predictions = pd.concat(frames, ignore_index=True)
    groups = {
        "strict_current_gse145996_phs000452": STRICT_EXTERNAL_COHORTS,
        "cbio_liu_dfci_only": [CBIO_COHORT],
        "strict_cbio_liu_plus_gse145996": [CBIO_COHORT, "GSE145996"],
    }
    rows = []
    for group_id, cohorts in groups.items():
        frame = predictions[predictions["cohort"].astype(str).isin(cohorts)].copy()
        metrics = compute_binary_metrics(
            frame["true_response_label"].astype(int),
            frame["response_probability"].astype(float),
            threshold=float(frame["threshold"].iloc[0]),
        )
        rows.append(
            {
                "group_id": group_id,
                "endpoint": "strict_recist",
                "candidate": str(selected_spec[0]["candidate"]),
                "n_samples": int(len(frame)),
                "n_responders": int(frame["true_response_label"].astype(int).sum()),
                "n_nonresponders": int((frame["true_response_label"].astype(int) == 0).sum()),
                **metrics,
            }
        )
    return pd.DataFrame(rows), predictions


def _baseline_predictions(
    X_by_cohort: dict[str, pd.DataFrame],
    primary_y: dict[str, pd.Series],
    strict_y: dict[str, pd.Series],
    cohorts: list[str],
) -> pd.DataFrame:
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
        for cohort in cohorts:
            labels = strict_y[cohort]
            raw = signature_score(X_by_cohort[cohort], genes).reindex(labels.index)
            score = (raw - raw.mean()) / (raw.std(ddof=0) + 1e-6)
            for sample_id, value in score.items():
                rows.append(
                    {
                        "cohort": cohort,
                        "sample_id": str(sample_id),
                        "model_name": model_name,
                        "true_response_label": int(labels.loc[sample_id]),
                        "response_probability": float(value),
                        "threshold": float(threshold),
                    }
                )
    return pd.DataFrame(rows)


def _family_gate(target_predictions: pd.DataFrame, baseline_predictions: pd.DataFrame, n_bootstrap: int) -> pd.DataFrame:
    groups = {
        "strict_current_gse145996_phs000452": STRICT_EXTERNAL_COHORTS,
        "cbio_liu_dfci_only": [CBIO_COHORT],
        "strict_cbio_liu_plus_gse145996": [CBIO_COHORT, "GSE145996"],
    }
    rng = np.random.default_rng(20260528)
    rows = []
    for group_id, cohorts in groups.items():
        target_frame = target_predictions[target_predictions["cohort"].astype(str).isin(cohorts)].copy()
        baseline_frame = baseline_predictions[baseline_predictions["cohort"].astype(str).isin(cohorts)].copy()
        target_frame["key"] = target_frame["cohort"].astype(str) + "::" + target_frame["sample_id"].astype(str)
        baseline_frame["key"] = baseline_frame["cohort"].astype(str) + "::" + baseline_frame["sample_id"].astype(str)
        y = target_frame.drop_duplicates("key").set_index("key")["true_response_label"].astype(int)
        target = target_frame.drop_duplicates("key").set_index("key")["response_probability"].astype(float)
        wide = baseline_frame.pivot_table(index="key", columns="model_name", values="response_probability", aggfunc="first").dropna()
        common = y.index.intersection(target.index).intersection(wide.index)
        y = y.loc[common]
        target = target.loc[common]
        wide = wide.loc[common]
        deltas = []
        for _ in range(n_bootstrap):
            idx = rng.integers(0, len(y), len(y))
            yy = y.to_numpy(dtype=int)[idx]
            if len(np.unique(yy)) < 2:
                continue
            tt = target.to_numpy(dtype=float)[idx]
            bb = wide.to_numpy(dtype=float)[idx]
            deltas.append(float(roc_auc_score(yy, tt) - np.mean([roc_auc_score(yy, bb[:, col]) for col in range(bb.shape[1])])))
        arr = np.asarray(deltas)
        baseline_aucs = {column: float(roc_auc_score(y, wide[column])) for column in wide.columns}
        rows.append(
            {
                "group_id": group_id,
                "target_model": str(target_frame["candidate"].iloc[0]),
                "baseline_family": "eight_strong_signatures",
                "n_samples": int(len(common)),
                "n_signatures": int(wide.shape[1]),
                "target_AUROC": float(roc_auc_score(y, target)),
                "family_mean_AUROC": float(np.mean(list(baseline_aucs.values()))),
                "best_signature": max(baseline_aucs, key=baseline_aucs.get),
                "best_signature_AUROC": float(max(baseline_aucs.values())),
                "delta_vs_family_mean": float(arr.mean()),
                "ci_low": float(np.quantile(arr, 0.025)),
                "ci_high": float(np.quantile(arr, 0.975)),
                "two_sided_p": float(min(1.0, 2.0 * min((arr <= 0).mean(), (arr >= 0).mean()))),
            }
        )
    result = pd.DataFrame(rows)
    result["two_sided_fdr_q"] = benjamini_hochberg(result["two_sided_p"])
    result["claim_level"] = np.where(
        (result["target_AUROC"] >= 0.70) & (result["delta_vs_family_mean"] > 0) & (result["two_sided_fdr_q"] <= 0.05),
        "strict_external_family_FDR_supported_numeric_target_met",
        np.where(result["delta_vs_family_mean"] > 0, "family_point_estimate_only", "family_not_superior"),
    )
    return result


def run_search(
    bulk_dir: Path,
    cbio_dir: Path,
    out_dir: Path,
    n_bootstrap: int,
    base_weight_policy: str,
    selection_policy: str,
    min_primary_ba: float,
) -> dict[str, Path]:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    X_by_cohort, metadata_by_cohort = _read_bulk(bulk_dir)
    _read_cbio(cbio_dir, X_by_cohort, metadata_by_cohort)
    genes = BASE_GENES + PAIR_GENES
    missing = {
        cohort: [gene for gene in genes if gene not in X_by_cohort[cohort].columns]
        for cohort in [*PRIMARY_COHORTS, *STRICT_EXTERNAL_COHORTS, CBIO_COHORT]
    }
    blocking = {cohort: values for cohort, values in missing.items() if values}
    if blocking:
        raise ValueError(f"Missing lipid-pair genes: {blocking}")
    primary_y = _labels(X_by_cohort, metadata_by_cohort, PRIMARY_COHORTS, "primary_recist")
    strict_y = _labels(X_by_cohort, metadata_by_cohort, [*STRICT_EXTERNAL_COHORTS, CBIO_COHORT], "strict_recist")
    transforms = _build_transforms(X_by_cohort, primary_y, strict_y, genes)
    specs = _candidate_specs(base_weight_policy)
    primary, primary_predictions = _primary_lodo(transforms, primary_y, specs, device)
    selection_frame = primary
    if selection_policy == "ba_guardrail":
        guarded = primary[primary["balanced_accuracy"].astype(float) >= min_primary_ba].copy()
        if not guarded.empty:
            selection_frame = guarded.sort_values(["AUROC", "AUPRC", "balanced_accuracy"], ascending=False)
    selected = selection_frame.iloc[0]
    selected_spec = [spec for spec in specs if str(spec["candidate"]) == str(selected["candidate"])]
    external_metrics, external_predictions = _score_locked_groups(transforms, primary_y, strict_y, selected_spec, device)
    baseline_predictions = _baseline_predictions(X_by_cohort, primary_y, strict_y, [*STRICT_EXTERNAL_COHORTS, CBIO_COHORT])
    gate = _family_gate(external_predictions, baseline_predictions, n_bootstrap=n_bootstrap)
    strict_row = external_metrics[external_metrics["group_id"].eq("strict_current_gse145996_phs000452")].iloc[0]
    cbio_row = external_metrics[external_metrics["group_id"].eq("cbio_liu_dfci_only")].iloc[0]
    cbio_gate = gate[gate["group_id"].eq("cbio_liu_dfci_only")].iloc[0]
    selection = pd.DataFrame(
        [
            {
                "selection_id": f"gpu_lipid_pair_{base_weight_policy}_primary_selected_rescue",
                "candidate": str(selected["candidate"]),
                "prior": "lipid_pi3k_pair",
                "base_weight_policy": base_weight_policy,
                "selection_policy": selection_policy,
                "min_primary_balanced_accuracy_guardrail": min_primary_ba,
                "device": str(device),
                "gpu_name": torch.cuda.get_device_name(0) if torch.cuda.is_available() else "cpu",
                "selection_boundary": "pair_candidate_and_weights_selected_by_primary_lodo_only_lipid_pi3k_prior_no_external_labels",
                "primary_AUROC": float(selected["AUROC"]),
                "primary_AUPRC": float(selected["AUPRC"]),
                "primary_balanced_accuracy": float(selected["balanced_accuracy"]),
                "primary_min_fold_AUROC": float(selected["min_fold_AUROC"]),
                "strict_external_AUROC": float(strict_row["AUROC"]),
                "strict_external_AUPRC": float(strict_row["AUPRC"]),
                "strict_external_balanced_accuracy": float(strict_row["balanced_accuracy"]),
                "strict_external_ECE": float(strict_row["ECE"]),
                "cbio_liu_AUROC": float(cbio_row["AUROC"]),
                "cbio_liu_AUPRC": float(cbio_row["AUPRC"]),
                "cbio_liu_balanced_accuracy": float(cbio_row["balanced_accuracy"]),
                "cbio_liu_ECE": float(cbio_row["ECE"]),
                "cbio_family_mean_AUROC": float(cbio_gate["family_mean_AUROC"]),
                "cbio_delta_vs_family_mean": float(cbio_gate["delta_vs_family_mean"]),
                "cbio_two_sided_fdr_q": float(cbio_gate["two_sided_fdr_q"]),
                "cbio_claim_level": str(cbio_gate["claim_level"]),
            }
        ]
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "primary_summary": out_dir / "gpu_lipid_pair_primary_summary.tsv",
        "primary_predictions": out_dir / "gpu_lipid_pair_primary_predictions.tsv",
        "external_metrics": out_dir / "gpu_lipid_pair_external_metrics.tsv",
        "external_predictions": out_dir / "gpu_lipid_pair_external_predictions.tsv",
        "baseline_predictions": out_dir / "gpu_lipid_pair_external_baseline_predictions.tsv",
        "family_gate": out_dir / "gpu_lipid_pair_external_family_gate.tsv",
        "selection": out_dir / "gpu_lipid_pair_rescue_selection.tsv",
        "markdown": out_dir / "GPU_LIPID_PAIR_RESCUE_AUDIT.md",
    }
    primary.to_csv(outputs["primary_summary"], sep="\t", index=False)
    primary_predictions.to_csv(outputs["primary_predictions"], sep="\t", index=False)
    external_metrics.to_csv(outputs["external_metrics"], sep="\t", index=False)
    external_predictions.to_csv(outputs["external_predictions"], sep="\t", index=False)
    baseline_predictions.to_csv(outputs["baseline_predictions"], sep="\t", index=False)
    gate.to_csv(outputs["family_gate"], sep="\t", index=False)
    selection.to_csv(outputs["selection"], sep="\t", index=False)
    lines = [
        "# GPU lipid/PI3K pair rescue audit",
        "",
        f"Device: {selection.iloc[0]['device']} ({selection.iloc[0]['gpu_name']}).",
        f"Selected candidate: `{selection.iloc[0]['candidate']}`.",
        "Candidate and weight selection used primary melanoma LODO only; strict external and cBioPortal labels were used only for locked scoring.",
        "",
        "Primary AUROC={:.3f}, AUPRC={:.3f}, balanced accuracy={:.3f}.".format(
            float(selection.iloc[0]["primary_AUROC"]),
            float(selection.iloc[0]["primary_AUPRC"]),
            float(selection.iloc[0]["primary_balanced_accuracy"]),
        ),
        "Strict current external AUROC={:.3f}, AUPRC={:.3f}, balanced accuracy={:.3f}, ECE={:.3f}.".format(
            float(selection.iloc[0]["strict_external_AUROC"]),
            float(selection.iloc[0]["strict_external_AUPRC"]),
            float(selection.iloc[0]["strict_external_balanced_accuracy"]),
            float(selection.iloc[0]["strict_external_ECE"]),
        ),
        "cBioPortal Liu/DFCI AUROC={:.3f}, AUPRC={:.3f}, ECE={:.3f}, family q={:.3f}.".format(
            float(selection.iloc[0]["cbio_liu_AUROC"]),
            float(selection.iloc[0]["cbio_liu_AUPRC"]),
            float(selection.iloc[0]["cbio_liu_ECE"]),
            float(selection.iloc[0]["cbio_two_sided_fdr_q"]),
        ),
    ]
    outputs["markdown"].write_text("\n".join(lines), encoding="utf-8")
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bulk-dir", default="data/processed/bulk")
    parser.add_argument("--cbio-dir", default="data/processed/cbioportal_melanoma")
    parser.add_argument("--out", default="results/gpu_lipid_pair_rescue_search_20260528")
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--base-weight-policy", choices=["all", "component_dominant"], default="all")
    parser.add_argument("--selection-policy", choices=["auroc", "ba_guardrail"], default="auroc")
    parser.add_argument("--min-primary-ba", type=float, default=0.65)
    args = parser.parse_args()
    outputs = run_search(
        ROOT / args.bulk_dir,
        ROOT / args.cbio_dir,
        ROOT / args.out,
        n_bootstrap=args.bootstrap,
        base_weight_policy=args.base_weight_policy,
        selection_policy=args.selection_policy,
        min_primary_ba=args.min_primary_ba,
    )
    selection = pd.read_csv(outputs["selection"], sep="\t")
    print(json.dumps(selection.to_dict("records"), ensure_ascii=False))
    print(f"Wrote {ROOT / args.out}")


if __name__ == "__main__":
    main()
