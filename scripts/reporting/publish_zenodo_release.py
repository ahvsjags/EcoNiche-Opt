from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.reporting.apply_zenodo_doi import apply_doi, normalize_doi
DEFAULT_METADATA = Path("deliverables/zenodo_release_metadata_20260527/.zenodo.json")
DEFAULT_MANIFEST = Path("deliverables/zenodo_release_metadata_20260527/zenodo_release_manifest.json")
DEFAULT_RELEASE_REPO = Path("github_release/econiche-opt-package-20260508")
DEFAULT_RELEASE_TAG = "v0.3.4-gpu-lipid-pair-rescue-20260528"
DEFAULT_REPO_URL = "https://github.com/ahvsjags/EcoNiche-Opt"
ZENODO_API = "https://zenodo.org/api"
SANDBOX_API = "https://sandbox.zenodo.org/api"


class ZenodoError(RuntimeError):
    pass


def _json_request(
    method: str,
    url: str,
    token: str,
    payload: dict[str, Any] | None = None,
    expected: tuple[int, ...] = (200,),
) -> dict[str, Any]:
    data = json.dumps(payload or {}).encode("utf-8") if payload is not None else None
    headers = {"Authorization": f"Bearer {token}"}
    if payload is not None:
        headers["Content-Type"] = "application/json"
    request = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            body = response.read().decode("utf-8")
            if response.status not in expected:
                raise ZenodoError(f"Unexpected Zenodo status {response.status}: {body[:500]}")
            return json.loads(body) if body else {}
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        raise ZenodoError(f"Zenodo API {method} {url} failed with {exc.code}: {body[:1000]}") from exc


def _upload_file(bucket_url: str, archive: Path, token: str, expected: tuple[int, ...] = (200, 201)) -> dict[str, Any]:
    target = bucket_url.rstrip("/") + "/" + urllib.parse.quote(archive.name)
    headers = {"Authorization": f"Bearer {token}", "Content-Type": "application/zip"}
    request = urllib.request.Request(target, data=archive.read_bytes(), headers=headers, method="PUT")
    try:
        with urllib.request.urlopen(request, timeout=600) as response:
            body = response.read().decode("utf-8")
            if response.status not in expected:
                raise ZenodoError(f"Unexpected Zenodo upload status {response.status}: {body[:500]}")
            return json.loads(body) if body else {}
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        raise ZenodoError(f"Zenodo file upload failed with {exc.code}: {body[:1000]}") from exc


def read_metadata(path: Path) -> dict[str, Any]:
    metadata = json.loads(path.read_text(encoding="utf-8"))
    required = {"title", "upload_type", "description", "creators", "license", "version"}
    missing = sorted(required - set(metadata))
    if missing:
        raise ValueError(f"Missing required Zenodo metadata fields: {', '.join(missing)}")
    return metadata


def build_git_archive(release_repo: Path, release_tag: str, out_dir: Path) -> Path:
    if not (release_repo / ".git").exists():
        raise FileNotFoundError(f"Release repository is not a git repo: {release_repo}")
    out_dir.mkdir(parents=True, exist_ok=True)
    archive = out_dir / f"EcoNiche-Opt-{release_tag}.zip"
    subprocess.check_call(
        [
            "git",
            "-c",
            f"safe.directory={release_repo.as_posix()}",
            "archive",
            "--format=zip",
            "--output",
            str(archive),
            release_tag,
        ],
        cwd=release_repo,
    )
    return archive


def get_archive(
    archive: Path | None,
    release_repo: Path,
    release_tag: str,
    out_dir: Path,
) -> Path:
    if archive is not None:
        resolved = archive.resolve()
        if not resolved.exists():
            raise FileNotFoundError(f"Archive does not exist: {resolved}")
        return resolved
    return build_git_archive(release_repo, release_tag, out_dir)


def load_token(token_env: str, explicit_token: str | None = None) -> str:
    token = explicit_token or os.environ.get(token_env, "")
    token = token.strip()
    if not token:
        raise ZenodoError(
            f"Missing Zenodo token. Set {token_env} with scopes deposit:write and deposit:actions, "
            "or pass --token-env pointing to an environment variable that contains it."
        )
    return token


def publish_release(
    *,
    api_base: str,
    metadata: dict[str, Any],
    archive: Path,
    token: str,
    publish: bool,
) -> dict[str, Any]:
    create = _json_request(
        "POST",
        api_base.rstrip("/") + "/deposit/depositions",
        token,
        {"metadata": metadata},
        expected=(201,),
    )
    deposition_id = create["id"]
    bucket_url = create["links"]["bucket"]
    upload = _upload_file(bucket_url, archive, token)
    result: dict[str, Any] = {
        "status": "draft_uploaded",
        "deposition_id": deposition_id,
        "deposition_url": create.get("links", {}).get("html"),
        "uploaded_file": archive.name,
        "uploaded_size_bytes": archive.stat().st_size,
        "upload_response": {
            "key": upload.get("key"),
            "checksum": upload.get("checksum"),
            "size": upload.get("size"),
        },
        "reserved_doi": create.get("metadata", {}).get("prereserve_doi", {}).get("doi"),
    }
    if publish:
        published = _json_request(
            "POST",
            api_base.rstrip() + f"/deposit/depositions/{deposition_id}/actions/publish",
            token,
            None,
            expected=(202,),
        )
        doi = published.get("doi")
        result.update(
            {
                "status": "published",
                "doi": doi,
                "doi_url": published.get("doi_url"),
                "record_url": published.get("record_url"),
                "record_id": published.get("record_id"),
            }
        )
    return result


def write_status(path: Path, status: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Create/upload/publish the EcoNiche-Opt release on Zenodo.")
    parser.add_argument("--metadata", default=str(DEFAULT_METADATA))
    parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    parser.add_argument("--release-repo", default=str(DEFAULT_RELEASE_REPO))
    parser.add_argument("--release-tag", default=DEFAULT_RELEASE_TAG)
    parser.add_argument("--archive", default=None, help="Optional existing zip archive. Otherwise git archive is built from --release-repo.")
    parser.add_argument("--archive-out-dir", default="deliverables/zenodo_upload_archive_20260528")
    parser.add_argument("--out", default="deliverables/zenodo_api_publication_status_20260528.json")
    parser.add_argument("--token-env", default="ZENODO_TOKEN")
    parser.add_argument("--token", default=None, help="Optional token value. Prefer --token-env to avoid shell history exposure.")
    parser.add_argument("--sandbox", action="store_true", help="Use sandbox.zenodo.org API.")
    parser.add_argument("--execute", action="store_true", help="Actually create a Zenodo draft and upload the archive.")
    parser.add_argument("--publish", action="store_true", help="Publish the Zenodo draft. Requires --execute.")
    parser.add_argument("--apply-doi", action="store_true", help="After publish, apply the minted DOI to local citation/manuscript files.")
    args = parser.parse_args()

    root = ROOT
    metadata = read_metadata(root / args.metadata)
    archive = get_archive(
        Path(args.archive) if args.archive else None,
        (root / args.release_repo).resolve(),
        args.release_tag,
        root / args.archive_out_dir,
    )
    api_base = SANDBOX_API if args.sandbox else ZENODO_API
    status: dict[str, Any] = {
        "status": "dry_run",
        "api_base": api_base,
        "release_tag": args.release_tag,
        "metadata_file": args.metadata,
        "archive": str(archive),
        "archive_size_bytes": archive.stat().st_size,
        "publish_requested": bool(args.publish),
        "apply_doi_requested": bool(args.apply_doi),
        "token_env": args.token_env,
        "token_present": bool(args.token or os.environ.get(args.token_env)),
    }

    if args.publish and not args.execute:
        raise SystemExit("--publish requires --execute.")
    if args.apply_doi and not args.publish:
        raise SystemExit("--apply-doi requires --publish.")

    if args.execute:
        try:
            token = load_token(args.token_env, args.token)
            status = publish_release(
                api_base=api_base,
                metadata=metadata,
                archive=archive,
                token=token,
                publish=args.publish,
            )
            status.update(
                {
                    "api_base": api_base,
                    "release_tag": args.release_tag,
                    "metadata_file": args.metadata,
                    "archive": str(archive),
                }
            )
            if args.apply_doi:
                doi = status.get("doi")
                if not doi:
                    raise ZenodoError("Published response did not include a DOI; refusing to update citation files.")
                normalized = normalize_doi(str(doi))
                changed = apply_doi(root, normalized)
                status["applied_doi"] = normalized
                status["doi_updated_files"] = [str(path) for path in changed]
        except (ValueError, ZenodoError) as exc:
            raise SystemExit(str(exc)) from exc

    write_status(root / args.out, status)
    printable = {key: value for key, value in status.items() if key != "token"}
    print(json.dumps(printable, ensure_ascii=False))


if __name__ == "__main__":
    main()
