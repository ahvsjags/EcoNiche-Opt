from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.reporting.publish_zenodo_release import (
    ZenodoError,
    get_archive,
    load_token,
    read_metadata,
)


def test_read_metadata_requires_core_zenodo_fields(tmp_path: Path):
    path = tmp_path / ".zenodo.json"
    path.write_text(
        json.dumps(
            {
                "title": "EcoNiche-Opt",
                "upload_type": "software",
                "description": "Release archive",
                "creators": [{"name": "Xu, Pengyuan"}],
                "license": "MIT",
                "version": "0.3.4",
            }
        ),
        encoding="utf-8",
    )

    metadata = read_metadata(path)

    assert metadata["upload_type"] == "software"
    assert metadata["version"] == "0.3.4"


def test_read_metadata_rejects_incomplete_metadata(tmp_path: Path):
    path = tmp_path / ".zenodo.json"
    path.write_text(json.dumps({"title": "EcoNiche-Opt"}), encoding="utf-8")

    with pytest.raises(ValueError, match="Missing required Zenodo metadata fields"):
        read_metadata(path)


def test_load_token_uses_env_without_printing_value(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setenv("ZENODO_TEST_TOKEN", "secret-token")

    assert load_token("ZENODO_TEST_TOKEN") == "secret-token"


def test_load_token_requires_env_or_explicit_token(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.delenv("ZENODO_TEST_TOKEN", raising=False)

    with pytest.raises(ZenodoError, match="Missing Zenodo token"):
        load_token("ZENODO_TEST_TOKEN")


def test_get_archive_accepts_existing_zip(tmp_path: Path):
    archive = tmp_path / "release.zip"
    archive.write_bytes(b"zip")

    assert get_archive(archive, tmp_path / "missing-repo", "v-test", tmp_path) == archive.resolve()
