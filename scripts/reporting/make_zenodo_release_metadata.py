from __future__ import annotations

import argparse
import json
import subprocess
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def _git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        return "UNAVAILABLE_LOCAL_GIT"


def build_metadata(out: Path, release_tag: str, repository_url: str) -> dict[str, object]:
    out.mkdir(parents=True, exist_ok=True)
    commit = _git_commit()
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
        "zenodo_doi": "RESULT_PENDING",
        "doi_status": "metadata_prepared_no_doi_minted",
        "required_before_citation": [
            "Create GitHub release from the frozen tag.",
            "Archive the release with Zenodo or an institutional repository.",
            "Replace RESULT_PENDING with the minted DOI in CITATION.cff, README, manuscript, and Code availability.",
        ],
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
        "- DOI status: `RESULT_PENDING` until Zenodo mints a DOI.",
        "",
        "## Before manuscript submission",
        "",
        "1. Push the frozen tag to GitHub.",
        "2. Create a GitHub release from the frozen tag.",
        "3. Archive the release in Zenodo.",
        "4. Replace `RESULT_PENDING` with the real DOI in public-facing files.",
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
    args = parser.parse_args()

    manifest = build_metadata(ROOT / args.out, args.release_tag, args.repository_url)
    print(json.dumps(manifest, ensure_ascii=False))


if __name__ == "__main__":
    main()
