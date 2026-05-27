from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
from sklearn import metrics as sk_metrics

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche_opt.data.registry import load_registry


def _registry_lookup(registry_path: Path) -> dict[str, dict[str, str]]:
    registry = load_registry(registry_path)
    rows: dict[str, dict[str, str]] = {}
    for cohort in registry.get("cohorts", []):
        accession = str(cohort.get("accession", "UNKNOWN"))
        timepoints = cohort.get("timepoints", cohort.get("timepoint", "unknown"))
        if isinstance(timepoints, list):
            timepoint_text = ";".join(map(str, timepoints))
        else:
            timepoint_text = str(timepoints)
        rows[accession] = {
            "therapy": str(cohort.get("therapy", "unknown")),
            "timepoint": timepoint_text or "unknown",
            "cancer_type": str(cohort.get("cancer_type", "unknown")),
        }
    return rows


def _metric_row(frame: pd.DataFrame, stratum: str, value: str) -> dict[str, object]:
    labels = pd.to_numeric(frame.get("true_label"), errors="coerce")
    probs = pd.to_numeric(frame.get("pred_prob"), errors="coerce")
    mask = labels.notna() & probs.notna()
    labels = labels[mask].astype(int)
    probs = probs[mask]
    auroc = float(sk_metrics.roc_auc_score(labels, probs)) if labels.nunique() == 2 else pd.NA
    auprc = float(sk_metrics.average_precision_score(labels, probs)) if labels.nunique() == 2 else pd.NA
    return {
        "stratum": stratum,
        "value": value,
        "n_samples": int(mask.sum()),
        "n_responders": int((labels == 0).sum()),
        "n_nonresponders": int((labels == 1).sum()),
        "AUROC": auroc,
        "AUPRC": auprc,
        "status": "computed_from_pipeline" if mask.any() else "RESULT_PENDING",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--predictions", default=None)
    parser.add_argument("--registry", default="config/data_registry.yml")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    mode = "demo" if args.demo else "real"
    predictions_path = ROOT / (args.predictions or f"results/{mode}/lodo_predictions.tsv")
    out_path = ROOT / (args.out or f"results/{mode}_stratification/therapy_timepoint_stratification.tsv")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if not predictions_path.exists():
        pd.DataFrame(
            [{"stratum": "all", "value": "all", "n_samples": 0, "status": "RESULT_PENDING", "reason": "predictions_missing"}]
        ).to_csv(out_path, sep="\t", index=False)
        print(f"Wrote {out_path}")
        return

    predictions = pd.read_csv(predictions_path, sep="\t")
    lookup = _registry_lookup(ROOT / args.registry)
    predictions["therapy"] = predictions["cohort"].astype(str).map(lambda c: lookup.get(c, {}).get("therapy", "unknown"))
    predictions["timepoint"] = predictions["cohort"].astype(str).map(lambda c: lookup.get(c, {}).get("timepoint", "unknown"))
    predictions["cancer_type"] = predictions["cohort"].astype(str).map(lambda c: lookup.get(c, {}).get("cancer_type", "unknown"))

    rows = [_metric_row(predictions, "all", "all")]
    for column in ["cohort", "therapy", "timepoint", "cancer_type"]:
        for value, group in predictions.groupby(column, dropna=False):
            rows.append(_metric_row(group, column, str(value)))
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
