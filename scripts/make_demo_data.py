from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.demo import make_synthetic_data
from econiche.priors import make_default_cell_state_priors


def main() -> None:
    demo = make_synthetic_data(n_cohorts=4, n_per_cohort=48, random_state=42)
    processed = ROOT / "data" / "processed" / "bulk"
    priors_dir = ROOT / "data" / "priors"
    results = ROOT / "results" / "demo"
    processed.mkdir(parents=True, exist_ok=True)
    priors_dir.mkdir(parents=True, exist_ok=True)
    results.mkdir(parents=True, exist_ok=True)

    for cohort, X in demo["X_by_cohort"].items():
        X.to_csv(processed / f"{cohort}.expr.tsv", sep="\t")
        demo["metadata_by_cohort"][cohort].to_csv(processed / f"{cohort}.metadata.tsv", sep="\t", index=False)

    common_genes = sorted(set.intersection(*(set(X.columns) for X in demo["X_by_cohort"].values())))
    (processed / "common_genes_primary.txt").write_text("\n".join(common_genes) + "\n", encoding="utf-8")
    make_default_cell_state_priors(common_genes).to_csv(priors_dir / "cell_state_priors.tsv", sep="\t")
    manifest = {
        "cohorts": list(demo["X_by_cohort"]),
        "n_genes": len(common_genes),
        "planted_genes": demo["planted_genes"],
        "label_encoding": "0=responder/sensitive, 1=non-responder/resistant",
    }
    (results / "demo_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Wrote demo data for {len(demo['X_by_cohort'])} cohorts to {processed}")


if __name__ == "__main__":
    main()
