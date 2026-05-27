from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_prospective_validation_package_is_locked_and_auditable() -> None:
    subprocess.run(
        [sys.executable, "scripts/validation/validate_prospective_package.py", "--package-dir", "deliverables/prospective_validation"],
        cwd=ROOT,
        check=True,
    )


def test_locked_validation_scorer_outputs_scores_and_metrics(tmp_path: Path) -> None:
    panel = pd.read_csv(ROOT / "deliverables" / "prospective_validation" / "locked_panel_genes.tsv", sep="\t")
    genes = sorted(panel["gene_symbol"].unique())
    samples = [f"VAL{i:02d}" for i in range(8)]
    rng = np.random.default_rng(20260527)
    expression = pd.DataFrame(rng.normal(size=(len(samples), len(genes))), index=samples, columns=genes)
    expression_path = tmp_path / "expression.tsv"
    expression.to_csv(expression_path, sep="\t")

    manifest = pd.DataFrame(
        {
            "site_id": ["SITE001"] * len(samples),
            "subject_id": samples,
            "sample_id": samples,
            "sample_source": ["tumor_tissue"] * len(samples),
            "baseline_status": ["pretreatment_before_first_ICB_dose"] * len(samples),
            "therapy": ["anti-PD-1_or_anti-PD-1-based"] * len(samples),
            "qc_pass": [True] * len(samples),
            "locked_validation_use_flag": [True] * len(samples),
        }
    )
    manifest_path = tmp_path / "manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)

    clinical = pd.DataFrame(
        {
            "subject_id": samples,
            "sample_id": samples,
            "response_raw": ["CR", "PR", "SD", "PD", "CR", "PD", "PR", "SD"],
            "recist_version": ["RECIST 1.1"] * len(samples),
            "primary_recist_label": [1, 1, 0, 0, 1, 0, 1, 0],
            "strict_recist_label": [1, 1, "", 0, 1, 0, 1, ""],
            "clinical_benefit_label": [1, 1, 1, 0, 1, 0, 1, 1],
            "label_source_document": ["demo"] * len(samples),
            "source_page_or_record_id": [str(i) for i in range(len(samples))],
        }
    )
    clinical_path = tmp_path / "clinical.tsv"
    clinical.to_csv(clinical_path, sep="\t", index=False)

    out_dir = tmp_path / "locked_validation"
    subprocess.run(
        [
            sys.executable,
            "-m",
            "econiche_opt.cli",
            "score-locked-validation",
            "--package-dir",
            "deliverables/prospective_validation",
            "--expression",
            str(expression_path),
            "--sample-manifest",
            str(manifest_path),
            "--clinical-annotation",
            str(clinical_path),
            "--out-dir",
            str(out_dir),
        ],
        cwd=ROOT,
        check=True,
    )
    scores = pd.read_csv(out_dir / "locked_validation_scores.tsv", sep="\t")
    metrics = pd.read_csv(out_dir / "locked_validation_metrics.tsv", sep="\t")
    audit = pd.read_csv(out_dir / "locked_validation_manifest_audit.tsv", sep="\t")
    assert set(scores["endpoint"]) == {"primary_recist", "strict_recist", "clinical_benefit"}
    assert scores["response_probability"].between(0, 1).all()
    assert "AUROC" in metrics.columns
    assert audit["is_valid"].all()
