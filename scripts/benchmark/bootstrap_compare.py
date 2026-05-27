from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.statistics import benjamini_hochberg, paired_bootstrap_delta


def main() -> None:
    econiche_path = ROOT / "results/real/lodo_predictions.tsv"
    baseline_path = ROOT / "results/real/baseline_predictions.tsv"
    out = ROOT / "results/real/model_comparison_bootstrap.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    if not econiche_path.exists() or not baseline_path.exists():
        pd.DataFrame([{"comparison": "EcoNiche-Opt_vs_baseline", "status": "RESULT_PENDING", "reason": "missing predictions"}]).to_csv(out, sep="\t", index=False)
        return
    eco = pd.read_csv(econiche_path, sep="\t")
    base = pd.read_csv(baseline_path, sep="\t")
    rows = []
    for model_name, frame in base.groupby("model_name"):
        merged = eco.merge(frame, on="sample_id", suffixes=("_eco", "_base"))
        merged = merged.dropna(subset=["pred_prob_eco", "pred_prob_base", "true_label_eco"])
        if len(merged) < 5:
            continue
        stats = paired_bootstrap_delta(merged["true_label_eco"], merged["pred_prob_eco"], merged["pred_prob_base"], n_bootstrap=200)
        stats.update({"comparison": f"EcoNiche-Opt_vs_{model_name}", "metric": "AUROC"})
        rows.append(stats)
    result = pd.DataFrame(rows) if rows else pd.DataFrame([{"comparison": "EcoNiche-Opt_vs_baseline", "status": "RESULT_PENDING", "reason": "no comparable baseline rows"}])
    if "p_value" in result.columns:
        result["fdr_q"] = benjamini_hochberg(result["p_value"].fillna(1.0))
    result.to_csv(out, sep="\t", index=False)
    result.to_csv(ROOT / "results/real/model_comparison_fdr.tsv", sep="\t", index=False)
    print(f"Wrote bootstrap comparison to {out}")


if __name__ == "__main__":
    main()
