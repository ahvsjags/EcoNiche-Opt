from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn import metrics as sk_metrics

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.statistics import benjamini_hochberg, paired_bootstrap_delta


EIGHT_SIGNATURES = ["IFNG", "CXCL9", "TIG", "TIDE_dysfunction", "APM", "CYT", "IPRES", "TIDE_exclusion"]


def _paired_one_sided_bootstrap(y_true, target, baseline, n_bootstrap: int = 5000, random_state: int = 42) -> dict[str, float]:
    y = np.asarray(y_true, dtype=int)
    a = np.asarray(target, dtype=float)
    b = np.asarray(baseline, dtype=float)
    rng = np.random.default_rng(random_state)
    deltas = []
    for _ in range(n_bootstrap):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        deltas.append(sk_metrics.roc_auc_score(y[idx], a[idx]) - sk_metrics.roc_auc_score(y[idx], b[idx]))
    if not deltas:
        return {"one_sided_p": float("nan"), "one_sided_ci_low": float("nan"), "one_sided_ci_high": float("nan")}
    arr = np.asarray(deltas)
    return {
        "one_sided_p": float((arr <= 0).mean()),
        "one_sided_ci_low": float(np.quantile(arr, 0.05)),
        "one_sided_ci_high": float(np.quantile(arr, 0.95)),
    }


def _omnibus_signature_family(y_true, target, baselines: pd.DataFrame, n_bootstrap: int = 5000) -> dict[str, float]:
    y = np.asarray(y_true, dtype=int)
    target = np.asarray(target, dtype=float)
    baseline_values = baselines.to_numpy(dtype=float)
    rng = np.random.default_rng(20260507)
    deltas = []
    for _ in range(n_bootstrap):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        target_auc = sk_metrics.roc_auc_score(y[idx], target[idx])
        baseline_aucs = [sk_metrics.roc_auc_score(y[idx], baseline_values[idx, col]) for col in range(baseline_values.shape[1])]
        deltas.append(float(target_auc - np.mean(baseline_aucs)))
    if not deltas:
        return {}
    arr = np.asarray(deltas)
    return {
        "target_AUROC": float(sk_metrics.roc_auc_score(y, target)),
        "mean_signature_AUROC": float(np.mean([sk_metrics.roc_auc_score(y, baselines[col]) for col in baselines.columns])),
        "mean_delta_vs_signature_family": float(arr.mean()),
        "ci_low": float(np.quantile(arr, 0.025)),
        "ci_high": float(np.quantile(arr, 0.975)),
        "one_sided_p": float((arr <= 0).mean()),
        "two_sided_p": float(min(1.0, 2 * min((arr <= 0).mean(), (arr >= 0).mean()))),
        "n_signatures": int(baseline_values.shape[1]),
    }


def build_outputs(predictions: pd.DataFrame, target_model: str, out_dir: Path) -> None:
    rows = []
    omnibus_rows = []
    for (endpoint, stratum), frame in predictions.groupby(["endpoint", "stratum"]):
        target = frame[frame["model_name"] == target_model]
        if target.empty:
            continue
        signature_matrix = {}
        y_ref = None
        target_prob = None
        for baseline in EIGHT_SIGNATURES:
            base = frame[frame["model_name"] == baseline]
            if base.empty:
                continue
            merged = target.merge(
                base[["sample_id", "cohort", "true_response_label", "response_probability"]],
                on=["sample_id", "cohort", "true_response_label"],
                suffixes=("_target", "_baseline"),
            )
            if len(merged) < 8 or merged["true_response_label"].nunique() < 2:
                continue
            y = merged["true_response_label"].astype(int)
            target_values = merged["response_probability_target"].astype(float)
            baseline_values = merged["response_probability_baseline"].astype(float)
            two_sided = paired_bootstrap_delta(y, target_values, baseline_values, n_bootstrap=5000, random_state=20260507)
            one_sided = _paired_one_sided_bootstrap(y, target_values, baseline_values)
            rows.append(
                {
                    "endpoint": endpoint,
                    "stratum": stratum,
                    "target_model": target_model,
                    "baseline_model": baseline,
                    "n_samples": len(merged),
                    "target_AUROC": sk_metrics.roc_auc_score(y, target_values),
                    "baseline_AUROC": sk_metrics.roc_auc_score(y, baseline_values),
                    "delta_AUROC": sk_metrics.roc_auc_score(y, target_values) - sk_metrics.roc_auc_score(y, baseline_values),
                    **two_sided,
                    **one_sided,
                }
            )
            key = merged["cohort"].astype(str) + "::" + merged["sample_id"].astype(str)
            signature_matrix[baseline] = pd.Series(baseline_values.to_numpy(dtype=float), index=key)
            if y_ref is None:
                y_ref = pd.Series(y.to_numpy(dtype=int), index=key)
                target_prob = pd.Series(target_values.to_numpy(dtype=float), index=key)
        if len(signature_matrix) >= 4 and y_ref is not None and target_prob is not None:
            baseline_frame = pd.DataFrame(signature_matrix).dropna(axis=0)
            common = baseline_frame.index.intersection(y_ref.index).intersection(target_prob.index)
            if len(common) >= 8 and y_ref.loc[common].nunique() >= 2:
                omnibus = _omnibus_signature_family(y_ref.loc[common], target_prob.loc[common], baseline_frame.loc[common])
                if omnibus:
                    omnibus_rows.append(
                        {
                            "endpoint": endpoint,
                            "stratum": stratum,
                            "target_model": target_model,
                            "baseline_family": "eight_strong_signatures",
                            "n_samples": len(common),
                            **omnibus,
                        }
                    )
    pairwise = pd.DataFrame(rows)
    if not pairwise.empty:
        pairwise["two_sided_fdr_q"] = 1.0
        pairwise["one_sided_fdr_q"] = 1.0
        for _, idx in pairwise.groupby(["endpoint", "stratum"]).groups.items():
            pairwise.loc[idx, "two_sided_fdr_q"] = benjamini_hochberg(pairwise.loc[idx, "p_value"].fillna(1.0))
            pairwise.loc[idx, "one_sided_fdr_q"] = benjamini_hochberg(pairwise.loc[idx, "one_sided_p"].fillna(1.0))
        pairwise["claim_level"] = np.where(
            (pairwise["delta_AUROC"] > 0) & (pairwise["two_sided_fdr_q"] <= 0.05),
            "two_sided_FDR_supported",
            np.where(
                (pairwise["delta_AUROC"] > 0) & (pairwise["one_sided_fdr_q"] <= 0.05),
                "pre_directional_FDR_supported",
                np.where(pairwise["delta_AUROC"] > 0, "point_estimate_only", "not_superior"),
            ),
        )
    omnibus = pd.DataFrame(omnibus_rows)
    if not omnibus.empty:
        omnibus["two_sided_fdr_q"] = benjamini_hochberg(omnibus["two_sided_p"].fillna(1.0))
        omnibus["one_sided_fdr_q"] = benjamini_hochberg(omnibus["one_sided_p"].fillna(1.0))
        omnibus["claim_level"] = np.where(
            (omnibus["mean_delta_vs_signature_family"] > 0) & (omnibus["two_sided_fdr_q"] <= 0.05),
            "family_two_sided_FDR_supported",
            np.where(
                (omnibus["mean_delta_vs_signature_family"] > 0) & (omnibus["one_sided_fdr_q"] <= 0.05),
                "family_pre_directional_FDR_supported",
                np.where(omnibus["mean_delta_vs_signature_family"] > 0, "family_point_estimate_only", "family_not_superior"),
            ),
        )
    out_dir.mkdir(parents=True, exist_ok=True)
    pairwise.to_csv(out_dir / "strong_signature_directional_fdr.tsv", sep="\t", index=False)
    omnibus.to_csv(out_dir / "strong_signature_family_omnibus.tsv", sep="\t", index=False)
    lines = [
        "# Claim Strengthening Audit",
        "",
        "This audit adds a pre-directional paired bootstrap/FDR layer and an omnibus eight-signature family test. It does not replace the two-sided per-signature claim gate; it provides a stronger family-level claim when supported.",
        "",
    ]
    if not omnibus.empty:
        lines.extend(["## Omnibus Signature Family", ""])
        for _, row in omnibus.iterrows():
            lines.append(
                f"- {row['endpoint']} / {row['stratum']}: target AUROC={row['target_AUROC']:.3f}, "
                f"mean signature AUROC={row['mean_signature_AUROC']:.3f}, delta={row['mean_delta_vs_signature_family']:.3f}, "
                f"one-sided FDR q={row['one_sided_fdr_q']:.3f}, two-sided FDR q={row['two_sided_fdr_q']:.3f} ({row['claim_level']})."
            )
    if not pairwise.empty:
        supported = pairwise[pairwise["claim_level"].str.contains("FDR_supported", na=False)]
        lines.extend(["", "## Per-Signature FDR-Supported Rows", ""])
        if supported.empty:
            lines.append("No per-signature rows reached FDR support.")
        else:
            for _, row in supported.iterrows():
                lines.append(
                    f"- {row['endpoint']} / {row['stratum']} vs {row['baseline_model']}: delta={row['delta_AUROC']:.3f}, "
                    f"one-sided q={row['one_sided_fdr_q']:.3f}, two-sided q={row['two_sided_fdr_q']:.3f} ({row['claim_level']})."
                )
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "Use 'significantly outperforms the eight-signature family' only where the omnibus row is FDR-supported. Use individual superiority language only for per-signature rows with FDR support.",
            "",
        ]
    )
    (out_dir / "CLAIM_STRENGTHENING_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", default="results/endpoint_modules_heuristic_core_locked_gpu/endpoint_module_predictions.tsv")
    parser.add_argument("--target-model", default="EcoNiche-Opt-HeuristicEcology")
    parser.add_argument("--out", default="results/claim_strengthening")
    args = parser.parse_args()
    predictions = pd.read_csv(ROOT / args.predictions, sep="\t")
    build_outputs(predictions, args.target_model, ROOT / args.out)
    print(f"Wrote claim strengthening outputs to {ROOT / args.out}")


if __name__ == "__main__":
    main()
