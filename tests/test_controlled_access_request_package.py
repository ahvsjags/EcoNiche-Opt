from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml

from scripts.reporting.make_controlled_access_request_package import build_package


def test_controlled_access_package_emits_targets_and_templates(tmp_path: Path):
    registry = {
        "cohorts": [
            {
                "accession": "phs001919",
                "name": "Abril",
                "layer": "A",
                "cancer_type": "melanoma",
                "therapy": "anti-PD1",
                "platform": "RNA-seq",
                "timepoints": ["controlled"],
                "endpoint": ["response"],
                "role": "high_value_restricted_external_candidate",
                "access": "DBGAP_ACCESS_RESTRICTED",
                "priority": "high",
                "download_script": "https://dbgap.example/phs001919",
                "preprocessing_script": "scripts/preprocess/preprocess_bulk.py",
                "uses": ["external_candidate", "melanoma_primary"],
            },
            {
                "accession": "GSE91061",
                "name": "Public",
                "layer": "A",
                "cancer_type": "melanoma",
                "therapy": "anti-PD1",
                "platform": "RNA-seq",
                "timepoints": ["pretreatment"],
                "endpoint": ["response"],
                "role": "discovery",
                "access": "public",
                "priority": "high",
                "download_script": "public",
                "preprocessing_script": "script",
                "uses": ["training"],
            },
        ]
    }
    registry_path = tmp_path / "registry.yml"
    registry_path.write_text(yaml.safe_dump(registry), encoding="utf-8")
    triage_path = tmp_path / "lead.tsv"
    pd.DataFrame(
        [
            {
                "lead_id": "ABRIL_RODRIGUEZ_PHS001919",
                "candidate_accession": "phs001919.v1.p1",
                "title": "Abril lead",
            }
        ]
    ).to_csv(triage_path, sep="\t", index=False)
    out = tmp_path / "pkg"

    targets = build_package(registry_path, triage_path, out)

    assert targets["accession"].tolist() == ["phs001919"]
    assert (out / "controlled_external_access_targets.tsv").exists()
    assert (out / "controlled_clinical_annotation_template.tsv").exists()
    assert (out / "controlled_assay_sample_manifest_template.tsv").exists()
    readme = (out / "README.md").read_text(encoding="utf-8")
    assert "Locked Boundary" in readme
    assert "RESULT_PENDING" in readme
