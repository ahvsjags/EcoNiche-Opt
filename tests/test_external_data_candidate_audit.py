from __future__ import annotations

from pathlib import Path

import yaml

from scripts.validation.audit_melanoma_external_data_candidates import build_candidate_audit


def test_restricted_ega_candidate_is_reported(tmp_path: Path):
    registry = {
        "cohorts": [
            {
                "accession": "EGAS00001001552",
                "name": "Lee_Rizos_PD1_melanoma_EGA",
                "layer": "A",
                "cancer_type": "melanoma",
                "therapy": "anti-PD1",
                "platform": "RNA-seq",
                "timepoints": ["pretreatment"],
                "endpoint": ["CR_PR_vs_SD_PD"],
                "role": "high_value_restricted_external_candidate",
                "access": "EGA_ACCESS_RESTRICTED",
                "priority": "high",
                "download_script": "https://ega-archive.org/studies/EGAS00001001552",
                "preprocessing_script": "scripts/preprocess/preprocess_bulk.py",
                "uses": ["external_candidate"],
                "notes": "controlled test row",
            }
        ]
    }
    path = tmp_path / "registry.yml"
    path.write_text(yaml.safe_dump(registry), encoding="utf-8")

    audit = build_candidate_audit(path)
    row = audit[audit["accession"].eq("EGAS00001001552")].iloc[0]

    assert row["normalized_access_status"] == "controlled"
    assert row["processed_status"] == "not_processed"
    assert "Request controlled access" in row["next_action"]
