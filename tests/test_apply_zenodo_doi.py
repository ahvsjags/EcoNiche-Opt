from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.reporting.apply_zenodo_doi import apply_doi, normalize_doi


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def test_normalize_zenodo_doi_accepts_doi_or_url():
    assert normalize_doi("10.5281/zenodo.12345678") == "10.5281/zenodo.12345678"
    assert normalize_doi("https://doi.org/10.5281/zenodo.12345678") == "10.5281/zenodo.12345678"


@pytest.mark.parametrize(
    "value",
    ["RESULT_PENDING", "10.5281/zenodo.RESULT_PENDING", "10.9999/zenodo.123", "https://example.org/10.5281/zenodo.123"],
)
def test_normalize_zenodo_doi_rejects_placeholders_and_non_zenodo(value: str):
    with pytest.raises(ValueError):
        normalize_doi(value)


def test_apply_zenodo_doi_updates_release_citation_files(tmp_path: Path):
    _write(
        tmp_path / "CITATION.cff",
        "cff-version: 1.2.0\ntitle: EcoNiche-Opt\nversion: 0.3.4\n",
    )
    _write(
        tmp_path / "README.md",
        "# EcoNiche-Opt\n\nNo DOI yet.\n",
    )
    _write(
        tmp_path / "DATA_RESULTS_FIGURES_UPLOAD_NOTES.md",
        "Release status: `v0.3.4-gpu-lipid-pair-rescue-20260528`. Zenodo metadata are prepared with\n"
        "`zenodo_doi=RESULT_PENDING`; no DOI should be cited until Zenodo or an\n"
        "institutional archive mints a real identifier.\n",
    )
    manifest = {
        "release_tag": "v0.3.4-gpu-lipid-pair-rescue-20260528",
        "repository_url": "https://github.com/ahvsjags/EcoNiche-Opt",
        "commit": "abc123",
        "zenodo_doi": "RESULT_PENDING",
        "doi_status": "metadata_prepared_no_doi_minted",
        "required_before_citation": ["Replace RESULT_PENDING with the minted DOI."],
    }
    _write(
        tmp_path / "deliverables/zenodo_release_metadata_20260527/zenodo_release_manifest.json",
        json.dumps(manifest),
    )
    _write(
        tmp_path / "deliverables/zenodo_manual_publication_steps_20260528.md",
        "Status: RESULT_PENDING until Zenodo mints a real DOI.\n"
        "After Zenodo publishes the record, copy the minted DOI. It must match the pattern `10.5281/zenodo.<record_id>`.\n"
        "The Zenodo DOI is not complete until the record is published by Zenodo. Until then, manuscript and citation files must retain `RESULT_PENDING` or omit the DOI.\n",
    )
    code_text = (
        "## Code availability\n\n"
        "The EcoNiche-Opt code repository is available at `https://github.com/ahvsjags/EcoNiche-Opt`. "
        "The manuscript version is `econiche-opt` v0.3.4, archived under release tag `v0.3.4-gpu-lipid-pair-rescue-20260528`; "
        "the release-specific source archive is available at `https://github.com/ahvsjags/EcoNiche-Opt/archive/refs/tags/v0.3.4-gpu-lipid-pair-rescue-20260528.zip`.\n"
    )
    _write(tmp_path / "paper/econiche_opt_manuscript_en_v1_20260509.md", code_text)
    _write(tmp_path / "paper/Journal of Translational Medicine投稿/EcoNiche-Opt_JTM_Main_Manuscript.md", code_text)
    _write(tmp_path / "paper/communications_medicine_submission/communications_medicine_submission_readiness.md", code_text)

    changed = apply_doi(tmp_path, "https://doi.org/10.5281/zenodo.12345678")

    assert Path("CITATION.cff") in changed
    assert "doi: 10.5281/zenodo.12345678" in (tmp_path / "CITATION.cff").read_text(encoding="utf-8")
    updated_manifest = json.loads(
        (tmp_path / "deliverables/zenodo_release_metadata_20260527/zenodo_release_manifest.json").read_text(encoding="utf-8")
    )
    assert updated_manifest["zenodo_doi"] == "10.5281/zenodo.12345678"
    assert updated_manifest["doi_status"] == "doi_minted"
    assert "https://doi.org/10.5281/zenodo.12345678" in (
        tmp_path / "paper/econiche_opt_manuscript_en_v1_20260509.md"
    ).read_text(encoding="utf-8")
