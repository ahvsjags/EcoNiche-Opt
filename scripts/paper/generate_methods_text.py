from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.registry import load_registry, write_registry_report


def main() -> None:
    registry = load_registry(ROOT / "config/data_registry.yml")
    audit_path = ROOT / "tables/dataset_access_audit.tsv"
    audit = pd.read_csv(audit_path, sep="\t") if audit_path.exists() else write_registry_report(registry, audit_path)
    text = f"""# Methods

## Dataset curation and access status

The dataset registry contains {len(audit)} cohorts. Access categories are audited before analysis, and controlled or unclear cohorts are not treated as available without verification.

## Response label harmonization

Response labels are encoded as 0 for responder or sensitive and 1 for non-responder or resistant. Primary RECIST maps CR/PR to responder and SD/PD to non-responder; strict RECIST excludes SD.

## Patient-level deduplication

All train/test splits use patient-level identifiers. When patient identifiers cannot be parsed, conservative IDs are constructed from accession, sample source, title, and sample tokens.

## Expression preprocessing

Bulk expression matrices are intersected on common genes and transformed with within-sample rank Gaussian normalization for cross-platform robustness.

## Construction of biological priors

Six cell-state prior columns are generated for tumor dedifferentiation, antigen presentation/MHC, T/NK effector activity, T-cell dysfunction, CAF/ECM exclusion, and myeloid suppression.

## EcoNiche-Opt mathematical formulation

EcoNiche-Opt scores a compact module for each state and combines state activities with logistic coefficients to estimate non-response probability.

## Optimization algorithm

The scaffold uses prior-weighted and correlation-aware module selection on training cohorts only. The objective reports AUROC lower-confidence-bound behavior, AUPRC, calibration error, biological prior support, and penalties.

## Baseline model scoring

Baseline signatures are emitted through a unified prediction schema. Missing signatures are labeled as unavailable_with_reason.

## Benchmark protocol

Primary benchmark evaluation uses leave-one-dataset-out folds. Gene directions, module optimization, coefficients, thresholds, and calibration are estimated without the holdout cohort.

## Statistical comparison

Paired bootstrap comparison and Benjamini-Hochberg FDR correction are provided. DeLong testing is left as a separate optional module.

## Single-cell validation

Single-cell scripts are mechanism-validation entrypoints only and must not use test cohort response labels to train priors.

## Survival analysis

Survival analysis is emitted only when OS/PFS columns are available. TCGA-style prognosis must not be described as ICB response prediction.

## Perturbation prioritization

Perturbation outputs are candidate axes and testable hypotheses, not validated treatments.

## Code and data availability

The demo pipeline is fully reproducible. Real public data retrieval is driven by the registry; controlled data require separate access approval.
"""
    out = ROOT / "paper/methods.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
