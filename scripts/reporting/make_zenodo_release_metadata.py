from __future__ import annotations

import argparse
import json
import subprocess
from datetime import date
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
ZENODO_DOI_RE = re.compile(r"^10\.5281/zenodo\.\d+$")


def _git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        return "UNAVAILABLE_LOCAL_GIT"


def _normalize_zenodo_doi(doi: str | None) -> str | None:
    if doi is None:
        return None
    value = doi.strip()
    for prefix in ("https://doi.org/", "http://doi.org/", "https://dx.doi.org/", "http://dx.doi.org/"):
        if value.lower().startswith(prefix):
            value = value[len(prefix) :]
            break
    value = value.rstrip("/")
    if not ZENODO_DOI_RE.fullmatch(value):
        raise ValueError(f"Invalid Zenodo DOI: {doi!r}")
    return value


def build_metadata(
    out: Path,
    release_tag: str,
    repository_url: str,
    commit_override: str | None = None,
    zenodo_doi: str | None = None,
) -> dict[str, object]:
    out.mkdir(parents=True, exist_ok=True)
    commit = commit_override or _git_commit()
    normalized_doi = _normalize_zenodo_doi(zenodo_doi)
    metadata = {
        "title": "EcoNiche-Opt: locked immune-ecology transcriptomic scoring framework",
        "upload_type": "software",
        "publication_date": date.today().isoformat(),
        "creators": [
            {
                "name": "Xu, Pengyuan",
                "affiliation": "Department of Materials Science and Engineering, Monash University",
            },
            {
                "name": "Yang, Guang",
                "affiliation": "School of Economics and Management, China University of Mining and Technology",
            },
            {
                "name": "Li, Moyan",
                "affiliation": "Hong Kong University of Science and Technology (Guangzhou)",
            },
        ],
        "description": (
            "EcoNiche-Opt provides reproducible code, locked scoring specifications, validation artifacts, "
            "and manuscript source data for multicohort immune-checkpoint blockade response benchmarking."
        ),
        "license": "MIT",
        "version": "0.3.4",
        "keywords": [
            "immune checkpoint blockade",
            "melanoma",
            "transcriptomics",
            "biomarker",
            "immune ecology",
            "machine learning",
            "reproducibility",
        ],
        "related_identifiers": [
            {
                "identifier": repository_url,
                "relation": "isSupplementTo",
                "scheme": "url",
            }
        ],
        "notes": (
            "Controlled-access or licensed datasets are not redistributed. Missing controlled outputs remain "
            "ACCESS_RESTRICTED or RESULT_PENDING until authorized files are imported through registered pipelines."
        ),
    }
    release_manifest = {
        "release_tag": release_tag,
        "repository_url": repository_url,
        "commit": commit,
        "zenodo_doi": normalized_doi or "RESULT_PENDING",
        "doi_status": "doi_minted" if normalized_doi else "metadata_prepared_no_doi_minted",
        "required_before_citation": (
            [
                "Zenodo DOI minted and recorded.",
                "Use the DOI URL in manuscript Code availability, CITATION.cff, README, and release notes.",
            ]
            if normalized_doi
            else [
                "Create GitHub release from the frozen tag.",
                "Archive the release with Zenodo or an institutional repository.",
                "Replace RESULT_PENDING with the minted DOI in CITATION.cff, README, manuscript, and Code availability.",
            ]
        ),
        "included_public_artifact_classes": [
            "package source code",
            "registered analysis scripts",
            "public metadata and manifests",
            "generated benchmark and validation tables",
            "manuscript figures and source data",
            "locked scoring specifications",
        ],
        "excluded_artifact_classes": [
            "controlled-access raw data",
            "licensed data",
            "oversized raw expression matrices not allowed for redistribution",
            "local caches",
        ],
    }
    if normalized_doi:
        release_manifest["zenodo_url"] = f"https://doi.org/{normalized_doi}"

    zenodo_path = out / ".zenodo.json"
    zenodo_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    manifest_path = out / "zenodo_release_manifest.json"
    manifest_path.write_text(json.dumps(release_manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    checklist = [
        "# Zenodo Release Checklist",
        "",
        f"- Release tag: `{release_tag}`",
        f"- Commit: `{commit}`",
        f"- Repository URL: {repository_url}",
        f"- DOI status: `{'doi_minted' if normalized_doi else 'RESULT_PENDING until Zenodo mints a DOI'}`",
        *( [f"- DOI: `{normalized_doi}`", f"- DOI URL: https://doi.org/{normalized_doi}"] if normalized_doi else [] ),
        "",
        "## Before manuscript submission",
        "",
        "1. Confirm the frozen tag is on GitHub.",
        "2. Confirm the GitHub release is public.",
        "3. Archive the release in Zenodo if no DOI is present.",
        "4. Use the real DOI in public-facing files after minting.",
        "5. Do not cite a placeholder DOI.",
        "",
    ]
    (out / "ZENODO_RELEASE_CHECKLIST.md").write_text("\n".join(checklist), encoding="utf-8")
    return release_manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="deliverables/zenodo_release_metadata_20260527")
    parser.add_argument("--release-tag", default="v0.3.4-gpu-lipid-pair-rescue-20260528")
    parser.add_argument("--repository-url", default="https://github.com/ahvsjags/EcoNiche-Opt")
    parser.add_argument("--commit", default=None)
    parser.add_argument("--zenodo-doi", default=None)
    args = parser.parse_args()

    manifest = build_metadata(ROOT / args.out, args.release_tag, args.repository_url, args.commit, args.zenodo_doi)
    print(json.dumps(manifest, ensure_ascii=False))


if __name__ == "__main__":
    main()
