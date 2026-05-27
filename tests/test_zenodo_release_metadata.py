from __future__ import annotations

import json
from pathlib import Path

from scripts.reporting.make_zenodo_release_metadata import build_metadata


def test_zenodo_release_metadata_keeps_doi_pending(tmp_path: Path):
    out = tmp_path / "zenodo"

    manifest = build_metadata(out, "v-test", "https://github.com/example/repo")

    assert manifest["zenodo_doi"] == "RESULT_PENDING"
    assert manifest["doi_status"] == "metadata_prepared_no_doi_minted"
    assert (out / ".zenodo.json").exists()
    assert (out / "ZENODO_RELEASE_CHECKLIST.md").exists()
    metadata = json.loads((out / ".zenodo.json").read_text(encoding="utf-8"))
    assert metadata["upload_type"] == "software"
    assert metadata["related_identifiers"][0]["identifier"] == "https://github.com/example/repo"
