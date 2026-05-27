from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.config import load_model_config
from econiche.io import load_processed_bulk
from econiche.metrics import compute_binary_metrics
from econiche.model import EcoNicheOpt
from econiche.priors import load_priors, make_default_cell_state_priors
from econiche.registry import load_registry


def main() -> None:
    registry = load_registry(ROOT / "config/data_registry.yml")
    cancer_type = {cohort["accession"]: cohort.get("cancer_type", "") for cohort in registry.get("cohorts", [])}
    X_by_cohort, y_by_cohort, metadata_by_cohort = load_processed_bulk(ROOT / "data/processed/bulk_real")
    train = {cohort: X for cohort, X in X_by_cohort.items() if "melanoma" in str(cancer_type.get(cohort, "")).lower()}
    test = {cohort: X for cohort, X in X_by_cohort.items() if cohort not in train}
    out = ROOT / "results/real/pancancer_transfer_metrics.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    if len(train) < 2 or not test:
        pd.DataFrame(
            [
                {
                    "mode": "no_retraining_transfer",
                    "status": "RESULT_PENDING",
                    "reason": "Need at least two melanoma training cohorts and one processed pan-cancer cohort",
                }
            ]
        ).to_csv(out, sep="\t", index=False)
        return
    genes = sorted(set.intersection(*(set(X.columns) for X in X_by_cohort.values())))
    prior_path = ROOT / "data/priors/cell_state_priors_real.tsv"
    priors = load_priors(prior_path) if prior_path.exists() else make_default_cell_state_priors(genes)
    model = EcoNicheOpt(load_model_config(ROOT / "config/model_config.yml"), priors=priors)
    model.fit(
        train,
        {cohort: y_by_cohort[cohort] for cohort in train},
        {cohort: metadata_by_cohort[cohort] for cohort in train},
    )
    rows = []
    predictions = []
    for cohort, X in test.items():
        pred = model.score_samples(X, metadata_by_cohort[cohort])
        pred["true_label"] = y_by_cohort[cohort].values
        pred["mode"] = "no_retraining_transfer"
        pred["cohort"] = cohort
        predictions.append(pred)
        metrics = compute_binary_metrics(pred["true_label"], pred["pred_prob"])
        metrics.update(
            {
                "mode": "no_retraining_transfer",
                "cohort": cohort,
                "cancer_type": cancer_type.get(cohort),
                "n_samples": len(pred),
                "status": "available",
            }
        )
        rows.append(metrics)
    rows.extend(
        [
            {"mode": "cancer_type_intercept_only", "status": "RESULT_PENDING", "reason": "Calibration cohort not configured yet"},
            {"mode": "state_level_recalibration", "status": "RESULT_PENDING", "reason": "Cancer-type-specific calibration not configured yet"},
        ]
    )
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    pd.concat(predictions, ignore_index=True).to_csv(ROOT / "results/real/pancancer_transfer_predictions.tsv", sep="\t", index=False)
    direction = pd.DataFrame({"state": list(model.model_.feature_names_in_), "coefficient": model.model_.coef_[0], "mode": "melanoma_trained"})
    direction.to_csv(ROOT / "results/real/pancancer_state_direction.tsv", sep="\t", index=False)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
