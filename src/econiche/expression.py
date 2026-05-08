from __future__ import annotations

import gzip
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests


def download_ncbi_gene_info(out: str | Path) -> Path:
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists() and out.stat().st_size > 1_000_000:
        return out
    url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
    response = requests.get(url, timeout=120)
    response.raise_for_status()
    out.write_bytes(response.content)
    return out


def download_ncbi_gene2ensembl(out: str | Path) -> Path:
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists() and out.stat().st_size > 1_000_000:
        return out
    url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz"
    response = requests.get(url, timeout=240)
    response.raise_for_status()
    out.write_bytes(response.content)
    return out


def load_entrez_symbol_map(path: str | Path) -> dict[str, str]:
    path = Path(path)
    mapping: dict[str, str] = {}
    if not path.exists():
        return mapping
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2:
                mapping[parts[1]] = parts[2]
    return mapping


def load_ensembl_symbol_map(path: str | Path, entrez_map: dict[str, str]) -> dict[str, str]:
    path = Path(path)
    mapping: dict[str, str] = {}
    if not path.exists():
        return mapping
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0] != "9606":
                continue
            gene_id = parts[1]
            ensembl_gene = parts[2]
            symbol = entrez_map.get(gene_id)
            if ensembl_gene and symbol:
                mapping[ensembl_gene] = symbol
    return mapping


def load_or_query_mygene_ensembl_symbol_map(ensembl_ids, out: str | Path, batch_size: int = 1000) -> dict[str, str]:
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cleaned = sorted({re.sub(r"\.\d+$", "", str(value).strip()) for value in ensembl_ids if str(value).startswith("ENSG")})
    existing = pd.DataFrame(columns=["ensembl_gene_id", "symbol"])
    if out.exists() and out.stat().st_size > 100:
        cached = pd.read_csv(out, sep="\t")
        if {"ensembl_gene_id", "symbol"}.issubset(cached.columns):
            mapping = dict(zip(cached["ensembl_gene_id"].astype(str), cached["symbol"].astype(str)))
            if len(set(cleaned) - set(mapping)) < max(10, int(len(cleaned) * 0.05)):
                return {key: value for key, value in mapping.items() if value and value != "nan"}
            existing = cached[["ensembl_gene_id", "symbol"]].dropna().copy()
    done = set(existing["ensembl_gene_id"].astype(str)) if not existing.empty else set()
    rows: list[dict[str, str]] = existing.to_dict("records")
    missing = [gene for gene in cleaned if gene not in done]
    for start in range(0, len(missing), batch_size):
        batch = missing[start : start + batch_size]
        response = None
        for attempt in range(6):
            try:
                response = requests.post(
                    "https://mygene.info/v3/query",
                    data={
                        "q": ",".join(batch),
                        "scopes": "ensembl.gene",
                        "fields": "symbol",
                        "species": "human",
                        "size": 1,
                    },
                    timeout=120,
                )
                response.raise_for_status()
                break
            except requests.RequestException:
                if attempt == 5:
                    raise
                time.sleep(2 * (attempt + 1))
        if response is None:
            continue
        for item in response.json():
            query = str(item.get("query", ""))
            symbol = str(item.get("symbol", "")).strip()
            if query and symbol:
                rows.append({"ensembl_gene_id": query, "symbol": symbol})
        frame = pd.DataFrame(rows).drop_duplicates("ensembl_gene_id") if rows else pd.DataFrame(columns=["ensembl_gene_id", "symbol"])
        frame.to_csv(out, sep="\t", index=False)
    frame = pd.DataFrame(rows).drop_duplicates("ensembl_gene_id") if rows else pd.DataFrame(columns=["ensembl_gene_id", "symbol"])
    frame.to_csv(out, sep="\t", index=False)
    return dict(zip(frame["ensembl_gene_id"].astype(str), frame["symbol"].astype(str)))


def download_geo_platform_soft(platform_id: str, out_dir: str | Path) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{platform_id}_platform.soft"
    if out.exists() and out.stat().st_size > 100_000:
        return out
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={platform_id}&targ=self&form=text&view=full"
    response = requests.get(url, timeout=240)
    response.raise_for_status()
    out.write_bytes(response.content)
    return out


def load_geo_probe_symbol_map(path: str | Path) -> dict[str, str]:
    path = Path(path)
    mapping: dict[str, str] = {}
    if not path.exists():
        return mapping
    in_table = False
    columns: list[str] | None = None
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if not in_table or not line:
                continue
            parts = line.split("\t")
            if columns is None:
                columns = parts
                continue
            if not columns or len(parts) < len(columns):
                continue
            row = dict(zip(columns, parts))
            probe = row.get("ID", "").strip()
            symbol = row.get("Symbol", "").strip() or row.get("ILMN_Gene", "").strip()
            if probe and symbol and not symbol.upper().startswith("ERCC-"):
                mapping[probe] = symbol
    return mapping


def normalize_gene_index(index, entrez_map: dict[str, str] | None = None) -> list[str]:
    entrez_map = entrez_map or {}
    genes = []
    for value in index:
        text = str(value).strip().replace('"', "")
        text = re.sub(r"\.\d+$", "", text)
        if text in entrez_map:
            text = entrez_map[text]
        if "_" in text and not text.upper().startswith("ILMN_") and re.match(r"^[A-Za-z0-9.-]+_\d+$", text):
            text = text.split("_", 1)[0]
        genes.append(text)
    return genes


def collapse_duplicate_genes(expr: pd.DataFrame) -> pd.DataFrame:
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr = expr.loc[expr.index.notna()]
    expr = expr[~expr.index.astype(str).str.startswith("nan")]
    return expr.groupby(expr.index).median(numeric_only=True)


def read_table_matrix(path: str | Path, gene_column: str | int | None = 0, sep: str | None = None) -> pd.DataFrame:
    path = Path(path)
    if sep is None:
        sep = "," if path.suffix == ".csv" or ".csv" in path.name else "\t"
    if path.suffix == ".xlsx":
        df = pd.read_excel(path)
    else:
        skiprows = _count_non_matrix_header_rows(path, sep)
        df = pd.read_csv(path, sep=sep, comment="#", compression="infer", low_memory=False, skiprows=skiprows)
    if df.empty:
        return pd.DataFrame()
    if gene_column is None:
        df = df.set_index(df.columns[0])
    elif isinstance(gene_column, int):
        df = df.set_index(df.columns[gene_column])
    else:
        df = df.set_index(gene_column)
    return df


def _count_non_matrix_header_rows(path: Path, sep: str) -> int:
    open_func = gzip.open if path.name.endswith(".gz") else open
    with open_func(path, "rt", encoding="utf-8", errors="replace") as handle:
        for idx, line in enumerate(handle):
            stripped = line.lstrip().strip('"')
            if stripped.startswith("#"):
                continue
            cells = [cell.strip().strip('"') for cell in line.rstrip("\n\r").split(sep)]
            if any(cells):
                return idx
    return 0


def _looks_like_sample(value: object) -> bool:
    text = str(value)
    return bool(re.search(r"^(Pt|Patient|BACI|RS-|GSM|RCC-|AZ|TE|SL|D\d+|N\d+|MGH|BAC)", text, flags=re.I))


def orient_expression(df: pd.DataFrame, known_genes: set[str] | None = None) -> pd.DataFrame:
    known_genes = known_genes or set()
    index_hits = sum(1 for x in df.index[:200] if str(x).split(".")[0] in known_genes or re.match(r"^[A-Z0-9.-]+$", str(x)))
    column_hits = sum(1 for x in df.columns[:200] if str(x).split(".")[0] in known_genes or re.match(r"^[A-Z0-9.-]+$", str(x)))
    if column_hits > index_hits and df.shape[0] < df.shape[1]:
        return df.T
    return df


def clean_expression_matrix(df: pd.DataFrame, entrez_map: dict[str, str] | None = None) -> pd.DataFrame:
    df = df.copy()
    df = df.loc[:, ~pd.Index(df.columns).astype(str).str.match(r"^Unnamed")]
    drop_cols = [col for col in df.columns if str(col).lower() in {"chr", "start", "end", "strand", "length", "locus", "target", "entrez_gene_id", "ncbi_name", "ncbi_accession", "gene_function"}]
    df = df.drop(columns=drop_cols, errors="ignore")
    detection_cols = [col for col in df.columns if "detection pval" in str(col).lower()]
    df = df.drop(columns=detection_cols, errors="ignore")
    sample_like_index = sum(_looks_like_sample(value) for value in df.index[: min(100, len(df.index))])
    if df.shape[1] > df.shape[0] and sample_like_index >= max(1, min(5, len(df.index) // 2)):
        df.columns = normalize_gene_index(df.columns, entrez_map)
        numeric = df.apply(pd.to_numeric, errors="coerce")
        numeric = numeric.dropna(axis=1, how="all").dropna(axis=0, how="all")
        return numeric.T.groupby(level=0).median(numeric_only=True).T
    df.index = normalize_gene_index(df.index, entrez_map)
    df = collapse_duplicate_genes(df)
    numeric = df.apply(pd.to_numeric, errors="coerce")
    numeric = numeric.dropna(axis=1, how="all").dropna(axis=0, how="all")
    return numeric.T


def choose_expression_file(accession: str, raw_dir: str | Path) -> Path | None:
    suppl = Path(raw_dir) / accession / "suppl"
    if not suppl.exists():
        return None
    preferences = {
        "GSE91061": ["fpkm"],
        "GSE78220": ["PatientFPKM"],
        "GSE145996": ["FPKM"],
        "GSE168204": ["counts"],
        "GSE115821": ["counts"],
        "GSE93157": ["raw_data_values"],
        "GSE244982": ["ProcessedData"],
        "GSE136961": ["TPM"],
        "GSE176307": ["log_trans_normalized", "salmon_tpm"],
        "GSE67501": ["Non-normalized"],
        "GSE121810": ["HUGO"],
        "GSE140901": ["processed_data"],
        "GSE165252": ["vst", "norm"],
        "GSE183924": ["FPKM"],
    }
    files = sorted([p for p in suppl.iterdir() if p.is_file() and not p.name.endswith(".tar")])
    for token in preferences.get(accession, []):
        for path in files:
            if token.lower() in path.name.lower():
                return path
    return files[0] if files else None
