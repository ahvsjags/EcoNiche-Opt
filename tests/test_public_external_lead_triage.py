from __future__ import annotations

from pathlib import Path

import yaml

from scripts.validation.audit_public_external_leads import build_lead_triage


def test_public_external_lead_triage_classifies_panel_controlled_and_ineligible(tmp_path: Path):
    registry = {
        "cohorts": [
            {"accession": "EGAS00001001552"},
            {"accession": "phs001919"},
            {"accession": "phs002683"},
            {"accession": "GSE123728"},
            {"accession": "GSE165745"},
            {"accession": "GSE122220"},
            {"accession": "IMvigor210"},
        ]
    }
    registry_path = tmp_path / "registry.yml"
    registry_path.write_text(yaml.safe_dump(registry), encoding="utf-8")

    triage = build_lead_triage(tmp_path, registry_path)
    by_lead = triage.set_index("lead_id")

    assert by_lead.loc["LEE_RIZOS_EGAS00001001552", "eligibility_status"] == "controlled_access_required"
    assert by_lead.loc["ABRIL_RODRIGUEZ_PHS001919", "strict_melanoma_primary_suitability"] == "potentially_high_value_after_access"
    assert by_lead.loc["MGH_HACOHEN_PHS002683", "eligibility_status"] == "controlled_access_required"
    assert by_lead.loc["GSE123728", "strict_melanoma_primary_suitability"] == "panel_transfer_not_bulk_primary"
    assert by_lead.loc["GSE165745", "eligibility_status"] == "panel_transfer_public"
    assert by_lead.loc["GSE122220", "eligibility_status"] == "low_n_array_public_processed"
    assert by_lead.loc["GSE122220", "strict_melanoma_primary_suitability"] == "low_n_platform_sensitivity_only"
    assert by_lead.loc["IMVIGOR210", "strict_melanoma_primary_suitability"] == "not_melanoma_primary_validation"
