from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[2]

XENA_DATASETS = [
    {
        "dataset_id": "TCGA_SKCM_Xena_HiSeqV2",
        "accession": "TCGA_SKCM_Xena",
        "source": "UCSC_Xena_TCGA_legacy_hub",
        "kind": "expression",
        "url": "https://tcga.xenahubs.net/download/TCGA.SKCM.sampleMap/HiSeqV2.gz",
        "filename": "TCGA.SKCM.HiSeqV2.gz",
        "required": True,
    },
    {
        "dataset_id": "TCGA_SKCM_Xena_clinicalMatrix",
        "accession": "TCGA_SKCM_Xena",
        "source": "UCSC_Xena_TCGA_legacy_hub",
        "kind": "clinical",
        "url": "https://tcga.xenahubs.net/download/TCGA.SKCM.sampleMap/SKCM_clinicalMatrix",
        "filename": "TCGA.SKCM.clinicalMatrix.tsv",
        "required": True,
    },
    {
        "dataset_id": "GDC_TCGA_SKCM_STAR_TPM_Xena",
        "accession": "GDC_TCGA_SKCM",
        "source": "UCSC_Xena_GDC_hub",
        "kind": "expression",
        "url": "https://gdc.xenahubs.net/download/TCGA-SKCM.star_tpm.tsv.gz",
        "filename": "GDC_TCGA_SKCM.star_tpm.tsv.gz",
        "required": False,
    },
]


@dataclass(frozen=True)
class ProbeResult:
    status: str
    http_status: int | None
    content_length: int | None
    content_type: str
    reason: str


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def probe_url(url: str, timeout: int = 30) -> ProbeResult:
    try:
        response = requests.head(url, allow_redirects=True, timeout=timeout)
        if response.status_code in {403, 405}:
            response = requests.get(url, stream=True, timeout=timeout, headers={"Range": "bytes=0-1"})
        status = "AVAILABLE" if response.status_code < 400 else "UNAVAILABLE"
        length = response.headers.get("content-length")
        return ProbeResult(
            status=status,
            http_status=response.status_code,
            content_length=int(length) if length and length.isdigit() else None,
            content_type=response.headers.get("content-type", ""),
            reason="" if status == "AVAILABLE" else response.text[:180],
        )
    except requests.RequestException as exc:
        return ProbeResult("UNAVAILABLE", None, None, "", str(exc))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(url: str, out_path: Path, timeout: int = 60) -> tuple[str, str]:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_suffix(out_path.suffix + ".part")
    with requests.get(url, stream=True, timeout=timeout) as response:
        response.raise_for_status()
        with tmp_path.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
    tmp_path.replace(out_path)
    return "DOWNLOADED", _sha256(out_path)


def query_gdc_file_manifest(out_dir: Path, size: int = 1000, timeout: int = 60) -> dict[str, Path | int | str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {
        "filters": {
            "op": "and",
            "content": [
                {"op": "=", "content": {"field": "cases.project.project_id", "value": "TCGA-SKCM"}},
                {"op": "=", "content": {"field": "data_type", "value": "Gene Expression Quantification"}},
                {"op": "=", "content": {"field": "analysis.workflow_type", "value": "STAR - Counts"}},
            ],
        },
        "fields": "file_id,file_name,data_type,data_format,analysis.workflow_type,cases.submitter_id,cases.case_id,cases.samples.sample_type",
        "format": "JSON",
        "size": str(size),
    }
    response = requests.post("https://api.gdc.cancer.gov/files", json=payload, timeout=timeout)
    response.raise_for_status()
    data = response.json()
    json_out = out_dir / "GDC_TCGA_SKCM_STAR_counts_manifest.json"
    json_out.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")
    rows = []
    for hit in data.get("data", {}).get("hits", []):
        cases = hit.get("cases") or [{}]
        case = cases[0] if cases else {}
        samples = case.get("samples") or [{}]
        sample = samples[0] if samples else {}
        rows.append(
            {
                "file_id": hit.get("file_id") or hit.get("id"),
                "file_name": hit.get("file_name"),
                "data_type": hit.get("data_type"),
                "data_format": hit.get("data_format"),
                "workflow_type": (hit.get("analysis") or {}).get("workflow_type"),
                "case_submitter_id": case.get("submitter_id"),
                "case_id": case.get("case_id"),
                "sample_type": sample.get("sample_type"),
            }
        )
    tsv_out = out_dir / "GDC_TCGA_SKCM_STAR_counts_manifest.tsv"
    pd.DataFrame(rows).to_csv(tsv_out, sep="\t", index=False)
    return {
        "json_out": json_out,
        "tsv_out": tsv_out,
        "n_files": len(rows),
        "status": "PASS" if rows else "RESULT_PENDING",
    }


def build_xena_manifest(
    out_dir: Path,
    download: bool = False,
    include_gdc_expression: bool = True,
    max_download_mb: float = 200.0,
    timeout: int = 60,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    datasets = XENA_DATASETS if include_gdc_expression else XENA_DATASETS[:2]
    out_dir.mkdir(parents=True, exist_ok=True)
    for dataset in datasets:
        out_path = out_dir / str(dataset["filename"])
        probe = probe_url(str(dataset["url"]), timeout=timeout)
        status = probe.status
        sha = ""
        reason = probe.reason
        if out_path.exists() and out_path.stat().st_size > 0:
            status = "EXISTING"
            sha = _sha256(out_path)
            reason = ""
        elif download and probe.status == "AVAILABLE":
            if probe.content_length is not None and probe.content_length > max_download_mb * 1024 * 1024:
                status = "RESULT_PENDING"
                reason = f"remote file exceeds --max-download-mb ({max_download_mb})"
            else:
                try:
                    status, sha = download_file(str(dataset["url"]), out_path, timeout=timeout)
                    reason = ""
                except requests.RequestException as exc:
                    status = "FAILED"
                    reason = str(exc)
        rows.append(
            {
                **dataset,
                "local_path": str(out_path),
                "status": status,
                "http_status": probe.http_status,
                "content_length": probe.content_length,
                "bytes_local": out_path.stat().st_size if out_path.exists() else 0,
                "sha256": sha,
                "checked_utc": _utc_now(),
                "reason": reason,
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Download or audit public TCGA SKCM data from UCSC Xena and GDC.")
    parser.add_argument("--out-dir", default="data/raw/TCGA_SKCM_Xena")
    parser.add_argument("--manifest", default="results/real/xena_tcga_manifest.tsv")
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--skip-gdc-expression", action="store_true")
    parser.add_argument("--max-download-mb", type=float, default=200.0)
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    out_dir = ROOT / args.out_dir if not Path(args.out_dir).is_absolute() else Path(args.out_dir)
    manifest_path = ROOT / args.manifest if not Path(args.manifest).is_absolute() else Path(args.manifest)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    manifest = build_xena_manifest(
        out_dir=out_dir,
        download=args.download,
        include_gdc_expression=not args.skip_gdc_expression,
        max_download_mb=args.max_download_mb,
        timeout=args.timeout,
    )
    try:
        gdc_manifest = query_gdc_file_manifest(out_dir, timeout=args.timeout)
        manifest.loc[len(manifest)] = {
            "dataset_id": "GDC_TCGA_SKCM_STAR_counts_API_manifest",
            "accession": "GDC_TCGA_SKCM",
            "source": "GDC_API_files_endpoint",
            "kind": "manifest",
            "url": "https://api.gdc.cancer.gov/files",
            "filename": Path(gdc_manifest["tsv_out"]).name,
            "required": False,
            "local_path": str(gdc_manifest["tsv_out"]),
            "status": gdc_manifest["status"],
            "http_status": 200,
            "content_length": "",
            "bytes_local": Path(gdc_manifest["tsv_out"]).stat().st_size,
            "sha256": _sha256(Path(gdc_manifest["tsv_out"])),
            "checked_utc": _utc_now(),
            "reason": f"{gdc_manifest['n_files']} GDC STAR-count files listed",
        }
    except requests.RequestException as exc:
        manifest.loc[len(manifest)] = {
            "dataset_id": "GDC_TCGA_SKCM_STAR_counts_API_manifest",
            "accession": "GDC_TCGA_SKCM",
            "source": "GDC_API_files_endpoint",
            "kind": "manifest",
            "url": "https://api.gdc.cancer.gov/files",
            "filename": "GDC_TCGA_SKCM_STAR_counts_manifest.tsv",
            "required": False,
            "local_path": str(out_dir / "GDC_TCGA_SKCM_STAR_counts_manifest.tsv"),
            "status": "FAILED",
            "http_status": "",
            "content_length": "",
            "bytes_local": 0,
            "sha256": "",
            "checked_utc": _utc_now(),
            "reason": str(exc),
        }

    manifest.to_csv(manifest_path, sep="\t", index=False)
    print(manifest[["dataset_id", "status", "bytes_local", "reason"]].to_string(index=False))
    print(f"Wrote {manifest_path}")
    if args.strict:
        required = manifest[manifest["required"].astype(bool)]
        failed = required[~required["status"].isin(["AVAILABLE", "DOWNLOADED", "EXISTING", "PASS"])]
        if not failed.empty:
            print(failed[["dataset_id", "status", "reason"]].to_string(index=False), file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
