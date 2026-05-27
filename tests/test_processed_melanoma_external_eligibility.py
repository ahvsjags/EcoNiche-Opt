from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml

from scripts.validation.audit_processed_melanoma_external_eligibility import build_audit


def test_processed_eligibility_classifies_known_strict_and_non_melanoma(tmp_path: Path):
    processed = tmp_path / "processed"
    processed.mkdir()
    pd.DataFrame({"label": [0, 1], "response_raw": ["PR", "PD"], "timepoint": ["pretreatment", "pretreatment"]}).to_csv(
        processed / "GSE145996.metadata.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame({"label": [0, 1], "response_raw": ["CR", "PD"], "timepoint": ["pretreatment", "pretreatment"]}).to_csv(
        processed / "GSE176307.metadata.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame({"label": [1, 0], "response_raw": ["PR", "PD"], "timepoint": ["pretreatment", "pretreatment"]}).to_csv(
        processed / "GSE122220.metadata.tsv",
        sep="\t",
        index=False,
    )
    registry = {
        "cohorts": [
            {
                "accession": "GSE145996",
                "cancer_type": "melanoma",
                "therapy": "anti-PD1",
                "platform": "RNA-seq",
                "timepoints": ["pretreatment"],
                "role": "external_small",
            },
            {
                "accession": "GSE176307",
                "cancer_type": "urothelial_cancer",
                "therapy": "ICB",
                "platform": "RNA-seq",
                "timepoints": ["pretreatment"],
                "role": "pan_cancer_external",
            },
            {
                "accession": "GSE122220",
                "cancer_type": "melanoma",
                "therapy": "ICB",
                "platform": "Illumina_HumanHT12_V4_expression_beadchip",
                "timepoints": ["pretreatment"],
                "role": "low_n_array_external_sensitivity",
            },
        ]
    }
    registry_path = tmp_path / "registry.yml"
    registry_path.write_text(yaml.safe_dump(registry), encoding="utf-8")

    audit = build_audit(processed, registry_path)
    statuses = audit.set_index("cohort")["eligibility_status"].to_dict()

    assert statuses["GSE145996"] == "strict_external_current"
    assert statuses["GSE176307"] == "not_melanoma_primary"
    assert statuses["GSE122220"] == "low_n_array_platform_sensitivity"
