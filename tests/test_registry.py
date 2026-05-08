from pathlib import Path

import yaml

from econiche.registry import audit_accession, load_registry, validate_registry, write_registry_report


def test_registry_validation_reports_missing_required_fields(tmp_path: Path):
    registry_path = tmp_path / "registry.yml"
    registry_path.write_text(
        yaml.safe_dump(
            {
                "cohorts": [
                    {
                        "accession": "GSE1",
                        "layer": "A",
                        "cancer_type": "melanoma",
                        "therapy": "anti-PD1",
                        "platform": "RNA-seq",
                        "timepoints": ["pretreatment"],
                        "endpoint": ["RECIST_response"],
                        "role": "discovery",
                        "access": "public",
                    },
                    {"accession": "CONTROLLED1", "access": "dbGaP_or_controlled_verify"},
                ]
            }
        ),
        encoding="utf-8",
    )

    registry = load_registry(registry_path)
    validation = validate_registry(registry)
    audit = audit_accession(registry)

    assert validation.loc[validation["accession"] == "GSE1", "is_valid"].item()
    assert not validation.loc[validation["accession"] == "CONTROLLED1", "is_valid"].item()
    assert audit.loc[audit["accession"] == "GSE1", "access_status"].item() == "public"
    assert audit.loc[audit["accession"] == "CONTROLLED1", "access_status"].item() == "controlled"


def test_write_registry_report_creates_audit_and_roles_tables(tmp_path: Path):
    registry = {
        "cohorts": [
            {
                "accession": "GSE1",
                "layer": "A",
                "cancer_type": "melanoma",
                "therapy": "anti-PD1",
                "platform": "RNA-seq",
                "timepoints": ["pretreatment"],
                "endpoint": ["RECIST_response"],
                "role": "discovery",
                "access": "public",
            }
        ]
    }

    audit_path = tmp_path / "dataset_access_audit.tsv"
    write_registry_report(registry, audit_path)

    assert audit_path.exists()
    assert (tmp_path / "dataset_roles.tsv").exists()
