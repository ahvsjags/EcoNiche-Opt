from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


CELL_MARKER_SIGNATURES: dict[str, tuple[str, ...]] = {
    "cd8_t_effector": ("CD3D", "CD3E", "CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "NKG7", "CCL5", "CXCL13"),
    "checkpoint_t_cell": ("PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "ICOS", "TOX", "CXCL13"),
    "treg": ("FOXP3", "IL2RA", "IKZF2", "CCR8", "CTLA4", "TIGIT", "TNFRSF18"),
    "b_cell_plasma": ("MS4A1", "CD19", "CD79A", "CD79B", "MZB1", "JCHAIN", "TNFRSF17", "BANK1"),
    "macrophage_myeloid": ("CD68", "CD163", "CSF1R", "LST1", "LYZ", "AIF1", "C1QA", "C1QB", "C1QC", "FCGR3A"),
    "dendritic_apc": ("CD1C", "CLEC10A", "ITGAX", "FCER1A", "LAMP3", "CCR7", "HLA-DRA", "HLA-DPB1"),
    "fibroblast_caf": ("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "ACTA2", "PDGFRB", "COL6A1"),
    "endothelial": ("PECAM1", "VWF", "KDR", "FLT1", "ENG", "CDH5", "ESAM", "RAMP2"),
    "melanoma_tumor": ("PMEL", "MLANA", "TYR", "MITF", "SOX10", "DCT", "MAGEA3", "S100B"),
    "ifng_axis": ("IFNG", "CXCL9", "CXCL10", "IDO1", "STAT1", "IRF1", "GBP1", "HLA-DRA"),
    "proliferation": ("MKI67", "TOP2A", "UBE2C", "PCNA", "MCM2", "MCM6", "TYMS"),
}


@dataclass(frozen=True)
class ProcessedCohort:
    cohort: str
    expression_path: Path
    metadata_path: Path | None = None


def discover_processed_cohorts(input_dir: str | Path, include_demo: bool = False) -> list[ProcessedCohort]:
    root = Path(input_dir)
    cohorts: list[ProcessedCohort] = []
    for expression_path in sorted(root.glob("*.expr.tsv")):
        cohort = expression_path.name[: -len(".expr.tsv")]
        if (not include_demo) and cohort.lower().startswith("demo_"):
            continue
        metadata_path = expression_path.with_name(f"{cohort}.metadata.tsv")
        cohorts.append(
            ProcessedCohort(
                cohort=cohort,
                expression_path=expression_path,
                metadata_path=metadata_path if metadata_path.exists() else None,
            )
        )
    return cohorts


def read_expression_matrix(path: str | Path, sample_ids: set[str] | None = None) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", index_col=0)
    frame.index = frame.index.map(str)
    frame.columns = frame.columns.map(str)
    if sample_ids:
        row_overlap = len(set(frame.index) & sample_ids)
        col_overlap = len(set(frame.columns) & sample_ids)
        if col_overlap > row_overlap:
            frame = frame.T
    frame = frame.apply(pd.to_numeric, errors="coerce")
    frame.columns = [str(col).strip().upper() for col in frame.columns]
    if frame.columns.duplicated().any():
        frame = frame.T.groupby(level=0).mean().T
    return frame.dropna(axis=1, how="all")


def read_metadata(path: str | Path | None) -> pd.DataFrame:
    if path is None or not Path(path).exists():
        return pd.DataFrame()
    metadata = pd.read_csv(path, sep="\t", dtype=str)
    if "sample_id" not in metadata.columns:
        metadata = metadata.reset_index(names="sample_id")
    metadata["sample_id"] = metadata["sample_id"].astype(str)
    return metadata


def _standardize_expression(frame: pd.DataFrame) -> pd.DataFrame:
    numeric = frame.astype(float)
    finite = numeric.replace([np.inf, -np.inf], np.nan)
    finite_values = finite.to_numpy(dtype=float)
    finite_max = np.nanmax(finite_values) if np.isfinite(finite_values).any() else np.nan
    finite_min = np.nanmin(finite_values) if np.isfinite(finite_values).any() else np.nan
    if np.isfinite(finite_max) and np.isfinite(finite_min) and finite_min >= 0 and finite_max > 50:
        finite = np.log2(finite + 1.0)
    centered = finite - finite.mean(axis=0, skipna=True)
    scale = finite.std(axis=0, skipna=True).replace(0, np.nan)
    z = centered.divide(scale, axis=1)
    return z.fillna(0.0)


def score_marker_abundance(
    expression: pd.DataFrame,
    cohort: str,
    signatures: Mapping[str, tuple[str, ...]] = CELL_MARKER_SIGNATURES,
    min_markers: int = 2,
) -> pd.DataFrame:
    if expression.empty:
        return pd.DataFrame(
            columns=[
                "sample_id",
                "cohort",
                "cell_type",
                "abundance_score",
                "relative_abundance",
                "n_markers_available",
                "marker_genes_used",
                "status",
            ]
        )
    z = _standardize_expression(expression)
    rows: list[pd.DataFrame] = []
    for cell_type, markers in signatures.items():
        available = [marker for marker in markers if marker in z.columns]
        if len(available) < min_markers:
            scores = pd.Series(np.nan, index=z.index, dtype=float)
            status = "insufficient_marker_overlap"
            relative = pd.Series(np.nan, index=z.index, dtype=float)
        else:
            scores = z[available].mean(axis=1)
            status = "scored_marker_z_baseline"
            relative = scores.rank(pct=True, method="average")
        rows.append(
            pd.DataFrame(
                {
                    "sample_id": z.index.astype(str),
                    "cohort": cohort,
                    "cell_type": cell_type,
                    "abundance_score": scores.to_numpy(dtype=float),
                    "relative_abundance": relative.to_numpy(dtype=float),
                    "n_markers_available": len(available),
                    "marker_genes_used": ";".join(available),
                    "status": status,
                }
            )
        )
    return pd.concat(rows, ignore_index=True)


def score_processed_cohorts(
    input_dir: str | Path,
    include_demo: bool = False,
    min_markers: int = 2,
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    score_frames: list[pd.DataFrame] = []
    metadata_by_cohort: dict[str, pd.DataFrame] = {}
    for cohort in discover_processed_cohorts(input_dir, include_demo=include_demo):
        metadata = read_metadata(cohort.metadata_path)
        sample_ids = set(metadata["sample_id"].astype(str)) if "sample_id" in metadata.columns else None
        expression = read_expression_matrix(cohort.expression_path, sample_ids=sample_ids)
        scores = score_marker_abundance(expression, cohort=cohort.cohort, min_markers=min_markers)
        scores["expression_file"] = str(cohort.expression_path)
        if cohort.metadata_path is not None:
            scores["metadata_file"] = str(cohort.metadata_path)
            metadata_by_cohort[cohort.cohort] = metadata
        else:
            scores["metadata_file"] = ""
        if not scores.empty:
            score_frames.append(scores)
    if not score_frames:
        return pd.DataFrame(), metadata_by_cohort
    return pd.concat(score_frames, ignore_index=True), metadata_by_cohort


def summarize_abundance(scores: pd.DataFrame) -> pd.DataFrame:
    if scores.empty:
        return pd.DataFrame(
            columns=[
                "cohort",
                "cell_type",
                "n_samples",
                "n_scored_samples",
                "mean_abundance_score",
                "median_abundance_score",
                "sd_abundance_score",
                "n_markers_available",
                "status",
            ]
        )
    summary = (
        scores.groupby(["cohort", "cell_type"], dropna=False)
        .agg(
            n_samples=("sample_id", "nunique"),
            n_scored_samples=("abundance_score", lambda value: int(value.notna().sum())),
            mean_abundance_score=("abundance_score", "mean"),
            median_abundance_score=("abundance_score", "median"),
            sd_abundance_score=("abundance_score", "std"),
            n_markers_available=("n_markers_available", "max"),
            status=("status", lambda value: "PASS" if (value == "scored_marker_z_baseline").any() else "RESULT_PENDING"),
        )
        .reset_index()
    )
    return summary


def _metadata_labels(metadata_by_cohort: Mapping[str, pd.DataFrame]) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for cohort, metadata in metadata_by_cohort.items():
        if metadata.empty or "sample_id" not in metadata.columns or "label" not in metadata.columns:
            continue
        labels = metadata[["sample_id", "label"]].copy()
        labels["cohort"] = cohort
        labels["label"] = pd.to_numeric(labels["label"], errors="coerce")
        labels = labels.dropna(subset=["label"])
        frames.append(labels)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["sample_id", "cohort", "label"])


def evaluate_abundance_baselines(scores: pd.DataFrame, metadata_by_cohort: Mapping[str, pd.DataFrame]) -> pd.DataFrame:
    labels = _metadata_labels(metadata_by_cohort)
    if scores.empty or labels.empty:
        return pd.DataFrame(
            columns=[
                "cohort",
                "cell_type",
                "n_labeled",
                "n_label_1",
                "n_label_0",
                "auroc_label1_as_high_score",
                "auprc_label1_as_high_score",
                "status",
                "interpretation",
            ]
        )
    merged = scores.merge(labels, on=["sample_id", "cohort"], how="inner")
    merged = merged.dropna(subset=["abundance_score", "label"])
    rows: list[dict[str, object]] = []
    groupers = [("all_public_processed", ["cell_type"]), ("by_cohort", ["cohort", "cell_type"])]
    for scope, columns in groupers:
        for key, group in merged.groupby(columns, dropna=False):
            labels_series = group["label"].astype(int)
            score_series = group["abundance_score"].astype(float)
            n_label_1 = int((labels_series == 1).sum())
            n_label_0 = int((labels_series == 0).sum())
            if len(group) >= 4 and n_label_1 > 0 and n_label_0 > 0:
                auroc = float(roc_auc_score(labels_series, score_series))
                auprc = float(average_precision_score(labels_series, score_series))
                status = "PASS"
            else:
                auroc = np.nan
                auprc = np.nan
                status = "RESULT_PENDING"
            if isinstance(key, tuple):
                values = key
            else:
                values = (key,)
            row = {
                "cohort": "all_public_processed" if scope == "all_public_processed" else values[0],
                "cell_type": values[-1],
                "n_labeled": int(len(group)),
                "n_label_1": n_label_1,
                "n_label_0": n_label_0,
                "auroc_label1_as_high_score": auroc,
                "auprc_label1_as_high_score": auprc,
                "status": status,
                "interpretation": "descriptive marker-abundance baseline; higher score is tested against label=1 without direction tuning",
            }
            rows.append(row)
    return pd.DataFrame(rows).sort_values(["status", "cohort", "cell_type"], ascending=[True, True, True])
