from __future__ import annotations

import pandas as pd

from scripts.analysis.run_phs000452_liu_subset_audit import normalize_liu_patient_id, source_concordance


def test_normalize_liu_patient_id_removes_trailing_t_only():
    assert normalize_liu_patient_id("Patient192_T") == "Patient192"
    assert normalize_liu_patient_id("Patient100") == "Patient100"
    assert normalize_liu_patient_id("Patient100_T_P") == "Patient100_T_P"


def test_source_concordance_counts_matched_response_mismatches():
    phs = pd.DataFrame(
        {
            "patient_id": ["Patient1", "Patient2_T"],
            "sample_id": ["SRR1", "SRR2"],
            "patient_id_raw": ["Patient1_T_P", "Patient2_T"],
            "response_raw": ["PR", "PD"],
            "response_NR": ["R", "N"],
            "m_stage": ["M1c", "M1b"],
        }
    )
    cbio = pd.DataFrame(
        {
            "patient_id": ["Patient1", "Patient2"],
            "sample_id": ["Sample1", "Sample2"],
            "response_raw": ["PR", "SD"],
            "response_raw_source": ["Partial Response", "Stable Disease"],
            "treatment": ["Nivolumab", "MK3475"],
        }
    )

    audit = source_concordance(phs, cbio)

    assert audit.loc[0, "matched_n"] == 2
    assert audit.loc[0, "response_mismatch_n"] == 1
