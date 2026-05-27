from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import requests
from bs4 import BeautifulSoup
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.registry import load_registry, normalize_access_status
from econiche_opt.data.download_plan import write_download_plan


GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"


def geo_prefix(accession: str) -> str:
    return accession[:-3] + "nnn"


def list_geo_ftp_files(accession: str, section: str) -> list[dict]:
    url = f"{GEO_FTP}/{geo_prefix(accession)}/{accession}/{section}/"
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    soup = BeautifulSoup(response.text, "html.parser")
    files = []
    for anchor in soup.find_all("a"):
        href = anchor.get("href", "")
        if not href or href == "../" or href.startswith("/") or href.startswith("http"):
            continue
        file_url = url + href
        size = None
        try:
            head = requests.head(file_url, allow_redirects=True, timeout=60)
            if head.headers.get("Content-Length"):
                size = int(head.headers["Content-Length"])
        except Exception:
            size = None
        files.append({"accession": accession, "section": section, "filename": href, "url": file_url, "remote_size": size})
    return files


def download_file(url: str, out: Path, remote_size: int | None = None) -> str:
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists() and remote_size and out.stat().st_size == remote_size:
        return "already_exists"
    with requests.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()
        total = int(response.headers.get("Content-Length", remote_size or 0))
        tmp = out.with_suffix(out.suffix + ".part")
        with tmp.open("wb") as handle, tqdm(total=total, unit="B", unit_scale=True, desc=out.name) as bar:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if not chunk:
                    continue
                handle.write(chunk)
                bar.update(len(chunk))
        tmp.replace(out)
    return "downloaded"


def extract_geo_metadata(accession: str, out_dir: Path, metadata_only: bool = True) -> pd.DataFrame:
    try:
        import GEOparse
    except Exception:
        return pd.DataFrame(
            [{"accession": accession, "sample_id": pd.NA, "status": "unavailable_with_reason", "reason": "GEOparse_not_installed"}]
        )
    try:
        gse = GEOparse.get_GEO(geo=accession, destdir=str(out_dir), how="full" if not metadata_only else "brief")
    except Exception as exc:
        return pd.DataFrame([{"accession": accession, "sample_id": pd.NA, "status": "download_failed", "reason": str(exc)}])
    rows = []
    for gsm_name, gsm in gse.gsms.items():
        metadata = gsm.metadata
        rows.append(
            {
                "sample_id": gsm_name,
                "patient_id_raw": pd.NA,
                "cohort": accession,
                "accession": accession,
                "platform": ";".join(metadata.get("platform_id", [])),
                "title": ";".join(metadata.get("title", [])),
                "source_name": ";".join(metadata.get("source_name_ch1", [])),
                "characteristics_ch1": "|".join(metadata.get("characteristics_ch1", [])),
                "therapy": pd.NA,
                "timepoint": pd.NA,
                "response_raw": pd.NA,
                "status": "metadata_downloaded",
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", default="config/data_registry.yml")
    parser.add_argument("--out-dir", default="data/raw")
    parser.add_argument("--metadata-only", action="store_true")
    parser.add_argument("--download-matrix", action="store_true", help="Download GEO series matrix files from NCBI FTP.")
    parser.add_argument("--download-supplementary", action="store_true", help="Download supplementary files from NCBI FTP.")
    parser.add_argument("--include-raw-tar", action="store_true", help="Include RAW.tar supplementary bundles when present.")
    parser.add_argument("--dry-run", action="store_true", help="Write a download plan without network access or file downloads.")
    args = parser.parse_args()
    if args.dry_run:
        plan = write_download_plan(ROOT / args.registry, ROOT / "data/metadata/download_dry_run.tsv")
        print(plan.to_string(index=False))
        print("Dry run only; no data downloaded.")
        return
    registry = load_registry(ROOT / args.registry)
    all_rows = []
    manifest_rows = []
    for cohort in registry.get("cohorts", []):
        accession = cohort["accession"]
        if not str(accession).startswith("GSE"):
            continue
        if normalize_access_status(cohort.get("access")) != "public":
            all_rows.append({"accession": accession, "status": "ACCESS_RESTRICTED_OR_UNKNOWN", "reason": cohort.get("access")})
            continue
        out_dir = ROOT / args.out_dir / accession
        out_dir.mkdir(parents=True, exist_ok=True)
        if args.metadata_only and not (args.download_matrix or args.download_supplementary):
            metadata = extract_geo_metadata(accession, out_dir, metadata_only=args.metadata_only)
            metadata.to_csv(out_dir / "metadata.tsv", sep="\t", index=False)
            all_rows.append(metadata)
        for section, should_download in [("matrix", args.download_matrix), ("suppl", args.download_supplementary)]:
            if not should_download:
                continue
            try:
                files = list_geo_ftp_files(accession, section)
            except Exception as exc:
                manifest_rows.append(
                    {
                        "accession": accession,
                        "section": section,
                        "filename": pd.NA,
                        "url": pd.NA,
                        "remote_size": pd.NA,
                        "local_path": pd.NA,
                        "status": "listing_failed",
                        "reason": str(exc),
                    }
                )
                continue
            for item in files:
                if item["filename"].endswith("RAW.tar") and not args.include_raw_tar:
                    status = "skipped_raw_tar"
                    local_path = ""
                else:
                    local = out_dir / section / item["filename"]
                    try:
                        status = download_file(item["url"], local, item["remote_size"])
                        local_path = str(local.relative_to(ROOT))
                    except Exception as exc:
                        status = "download_failed"
                        local_path = str(local.relative_to(ROOT))
                        item["reason"] = str(exc)
                manifest_rows.append(
                    {
                        **item,
                        "local_path": local_path,
                        "status": status,
                        "reason": item.get("reason", ""),
                    }
                )
    summary = (
        pd.concat([row if isinstance(row, pd.DataFrame) else pd.DataFrame([row]) for row in all_rows], ignore_index=True, sort=False)
        if all_rows
        else pd.DataFrame(columns=["accession", "sample_id", "status", "reason"])
    )
    summary_out = ROOT / "data/metadata/geo_download_summary.tsv"
    summary_out.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_out, sep="\t", index=False)
    manifest = pd.DataFrame(manifest_rows)
    if not manifest.empty:
        manifest.to_csv(ROOT / "data/metadata/geo_file_manifest.tsv", sep="\t", index=False)
    print(f"Wrote GEO download summary to {summary_out}")


if __name__ == "__main__":
    main()
