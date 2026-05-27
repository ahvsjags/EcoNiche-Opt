from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.config import load_model_config
from econiche.io import load_processed_bulk
from econiche.model import EcoNicheOpt
from econiche.priors import load_priors, make_default_cell_state_priors
from econiche.utils import write_json


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/model_config.yml")
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--include-demo", action="store_true")
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--priors", default="data/priors/cell_state_priors.tsv")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    cfg = load_model_config(ROOT / args.config)
    X_by_cohort, y_by_cohort, metadata_by_cohort = load_processed_bulk(ROOT / args.processed_dir)
    if args.demo:
        X_by_cohort = {k: v for k, v in X_by_cohort.items() if k.startswith("demo_cohort_")}
        y_by_cohort = {k: v for k, v in y_by_cohort.items() if k in X_by_cohort}
        metadata_by_cohort = {k: v for k, v in metadata_by_cohort.items() if k in X_by_cohort}
    elif not args.include_demo:
        X_by_cohort = {k: v for k, v in X_by_cohort.items() if not k.startswith("demo_cohort_")}
        y_by_cohort = {k: v for k, v in y_by_cohort.items() if k in X_by_cohort}
        metadata_by_cohort = {k: v for k, v in metadata_by_cohort.items() if k in X_by_cohort}
    if len(X_by_cohort) < 2:
        raise SystemExit("Need at least two processed cohorts. Run scripts/make_demo_data.py first for demo mode.")

    prior_path = ROOT / args.priors
    if prior_path.exists():
        priors = load_priors(str(prior_path))
    else:
        genes = sorted(set.intersection(*(set(X.columns) for X in X_by_cohort.values())))
        priors = make_default_cell_state_priors(genes)
    model = EcoNicheOpt(cfg, priors=priors)
    result = model.fit(X_by_cohort, y_by_cohort, metadata_by_cohort)

    out_dir = ROOT / (args.out or ("results/demo" if args.demo else "results/real"))
    out_dir.mkdir(parents=True, exist_ok=True)
    result.module_table().to_csv(out_dir / "econiche_module.tsv", sep="\t", index=False)
    result.lodo_metrics.to_csv(out_dir / "lodo_metrics.tsv", sep="\t", index=False)
    result.predictions.to_csv(out_dir / "lodo_predictions.tsv", sep="\t", index=False)
    result.history.to_csv(out_dir / "objective_history.tsv", sep="\t", index=False)
    result.coefficients.to_csv(out_dir / "coefficients.tsv", sep="\t", index=False)
    write_json(result.objective_terms, out_dir / "objective_terms.json")
    print(f"Wrote EcoNiche-Opt outputs to {out_dir}")


if __name__ == "__main__":
    main()
