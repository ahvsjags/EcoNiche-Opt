from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche_opt.perturbation.reversal import build_reversal_signature  # noqa: E402

FIGSHARE_ARTICLE_ID = 25880521
FIGSHARE_API_URL = f"https://api.figshare.com/v2/articles/{FIGSHARE_ARTICLE_ID}"
NEEDED_FILES = {"CRISPRGeneEffect.csv", "Model.csv"}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch_figshare_manifest(timeout: int = 60) -> dict:
    response = requests.get(FIGSHARE_API_URL, timeout=timeout)
    response.raise_for_status()
    return response.json()


def _download(url: str, out: Path, timeout: int = 120) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".part")
    with requests.get(url, stream=True, timeout=timeout) as response:
        response.raise_for_status()
        with tmp.open("wb") as handle:
            for chunk in response.iter_content(1024 * 1024):
                if chunk:
                    handle.write(chunk)
    tmp.replace(out)


def ensure_depmap_files(raw_dir: Path, download: bool, max_download_mb: float, timeout: int) -> tuple[dict[str, Path], pd.DataFrame]:
    manifest = fetch_figshare_manifest(timeout=timeout)
    file_rows = []
    available: dict[str, Path] = {}
    files = {file_info["name"]: file_info for file_info in manifest.get("files", [])}
    for name in sorted(NEEDED_FILES):
        info = files.get(name)
        out = raw_dir / name
        if info is None:
            status = "RESULT_PENDING"
            reason = "file not listed in Figshare article"
            size = 0
            url = ""
        else:
            size = int(info.get("size") or 0)
            url = str(info.get("download_url") or "")
            if out.exists() and out.stat().st_size > 0:
                status = "EXISTING"
                reason = ""
                available[name] = out
            elif download and size <= max_download_mb * 1024 * 1024:
                _download(url, out, timeout=timeout)
                status = "DOWNLOADED"
                reason = ""
                available[name] = out
            elif download:
                status = "RESULT_PENDING"
                reason = f"remote file exceeds --max-download-mb ({max_download_mb})"
            else:
                status = "RESULT_PENDING"
                reason = "run with --download to cache the public DepMap file"
        file_rows.append(
            {
                "release_title": manifest.get("title", ""),
                "article_id": FIGSHARE_ARTICLE_ID,
                "file_name": name,
                "download_url": url,
                "remote_size_bytes": size,
                "local_path": str(out),
                "local_size_bytes": out.stat().st_size if out.exists() else 0,
                "sha256": _sha256(out) if out.exists() and out.stat().st_size > 0 else "",
                "status": status,
                "reason": reason,
            }
        )
    return available, pd.DataFrame(file_rows)


def _depmap_symbol(column: str) -> str:
    match = re.match(r"^(.+?)\s+\(\d+\)$", str(column))
    return (match.group(1) if match else str(column)).strip().upper()


def _model_id_column(frame: pd.DataFrame) -> str:
    for candidate in ["ModelID", "DepMap_ID", "depmap_id", "model_id"]:
        if candidate in frame.columns:
            return candidate
    return str(frame.columns[0])


def _melanoma_model_ids(model_path: Path) -> set[str]:
    models = pd.read_csv(model_path, dtype=str)
    model_col = _model_id_column(models)
    text_columns = [
        column
        for column in ["OncotreeLineage", "OncotreePrimaryDisease", "PrimaryDisease", "Lineage", "Subtype", "StrippedCellLineName"]
        if column in models.columns
    ]
    if not text_columns:
        return set(models[model_col].dropna().astype(str))
    text = models[text_columns].fillna("").agg(" ".join, axis=1).str.lower()
    return set(models.loc[text.str.contains("melanoma", regex=False), model_col].dropna().astype(str))


def score_depmap_dependencies(module: pd.DataFrame, gene_effect_path: Path, model_path: Path, top_genes: int) -> pd.DataFrame:
    signature = build_reversal_signature(module, top_genes=top_genes)
    genes = list(dict.fromkeys(signature["resistance_up"] + signature["resistance_down"]))
    if not genes:
        return pd.DataFrame()
    header = pd.read_csv(gene_effect_path, nrows=0)
    id_col = str(header.columns[0])
    column_by_gene = {_depmap_symbol(column): str(column) for column in header.columns[1:]}
    use_gene_columns = [column_by_gene[gene] for gene in genes if gene in column_by_gene]
    if not use_gene_columns:
        return pd.DataFrame(
            [
                {
                    "target_gene": gene,
                    "depmap_score": pd.NA,
                    "mean_gene_effect_all": pd.NA,
                    "mean_gene_effect_melanoma": pd.NA,
                    "n_models": 0,
                    "n_melanoma_models": 0,
                    "status": "RESULT_PENDING",
                    "reason": "gene not present in DepMap CRISPRGeneEffect columns",
                    "interpretation": "dependency support unavailable for this module gene",
                }
                for gene in genes
            ]
        )
    effects = pd.read_csv(gene_effect_path, usecols=[id_col, *use_gene_columns])
    effects[id_col] = effects[id_col].astype(str)
    melanoma_ids = _melanoma_model_ids(model_path)
    melanoma = effects[effects[id_col].isin(melanoma_ids)]
    rows = []
    reverse_column_map = {_depmap_symbol(column): column for column in use_gene_columns}
    for gene in genes:
        column = reverse_column_map.get(gene)
        if column is None:
            rows.append(
                {
                    "target_gene": gene,
                    "depmap_score": pd.NA,
                    "mean_gene_effect_all": pd.NA,
                    "mean_gene_effect_melanoma": pd.NA,
                    "n_models": 0,
                    "n_melanoma_models": 0,
                    "status": "RESULT_PENDING",
                    "reason": "gene not present in DepMap CRISPRGeneEffect columns",
                    "interpretation": "dependency support unavailable for this module gene",
                }
            )
            continue
        all_values = pd.to_numeric(effects[column], errors="coerce")
        melanoma_values = pd.to_numeric(melanoma[column], errors="coerce")
        mean_all = float(all_values.mean()) if all_values.notna().any() else pd.NA
        mean_melanoma = float(melanoma_values.mean()) if melanoma_values.notna().any() else pd.NA
        depmap_score = -mean_melanoma if pd.notna(mean_melanoma) else pd.NA
        rows.append(
            {
                "target_gene": gene,
                "depmap_score": depmap_score,
                "mean_gene_effect_all": mean_all,
                "mean_gene_effect_melanoma": mean_melanoma,
                "n_models": int(all_values.notna().sum()),
                "n_melanoma_models": int(melanoma_values.notna().sum()),
                "status": "DepMap_public_dependency_support",
                "reason": "",
                "interpretation": "higher depmap_score means stronger average melanoma dependency because DepMap gene-effect values are more negative for dependency",
            }
        )
    return pd.DataFrame(rows).sort_values("depmap_score", ascending=False, na_position="last")


def _pending(reason: str) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "target_gene": "RESULT_PENDING",
                "depmap_score": "",
                "mean_gene_effect_all": "",
                "mean_gene_effect_melanoma": "",
                "n_models": 0,
                "n_melanoma_models": 0,
                "status": "RESULT_PENDING",
                "reason": reason,
                "interpretation": "DepMap dependency support not computed",
            }
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Score EcoNiche module genes against public DepMap CRISPR dependency data.")
    parser.add_argument("--module", default="results/real/econiche_module.tsv")
    parser.add_argument("--raw-dir", default="data/raw/DepMap")
    parser.add_argument("--out", default="results/perturbation/depmap_targets.tsv")
    parser.add_argument("--manifest", default="results/perturbation/depmap_manifest.tsv")
    parser.add_argument("--top-genes", type=int, default=120)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--max-download-mb", type=float, default=600.0)
    parser.add_argument("--timeout", type=int, default=180)
    args = parser.parse_args()

    raw_dir = ROOT / args.raw_dir if not Path(args.raw_dir).is_absolute() else Path(args.raw_dir)
    out = ROOT / args.out if not Path(args.out).is_absolute() else Path(args.out)
    manifest_out = ROOT / args.manifest if not Path(args.manifest).is_absolute() else Path(args.manifest)
    out.parent.mkdir(parents=True, exist_ok=True)
    manifest_out.parent.mkdir(parents=True, exist_ok=True)

    try:
        files, manifest = ensure_depmap_files(raw_dir, args.download, args.max_download_mb, args.timeout)
    except requests.RequestException as exc:
        _pending(f"DepMap Figshare manifest/download failed: {exc}").to_csv(out, sep="\t", index=False)
        print(f"Wrote {out}")
        return 0
    manifest.to_csv(manifest_out, sep="\t", index=False)

    module_path = ROOT / args.module if not Path(args.module).is_absolute() else Path(args.module)
    if not module_path.exists():
        result = _pending(f"module table missing: {module_path}")
    elif not NEEDED_FILES.issubset(files):
        missing = sorted(NEEDED_FILES - set(files))
        result = _pending(f"DepMap files not cached: {', '.join(missing)}")
    else:
        module = pd.read_csv(module_path, sep="\t")
        result = score_depmap_dependencies(module, files["CRISPRGeneEffect.csv"], files["Model.csv"], args.top_genes)
        if result.empty:
            result = _pending("no module genes overlapped DepMap CRISPRGeneEffect columns")
    result.to_csv(out, sep="\t", index=False)
    print(f"Wrote {manifest_out}")
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
