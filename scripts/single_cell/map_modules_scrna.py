from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def main() -> None:
    module_path = ROOT / "results/real/econiche_module.tsv"
    expr_path = ROOT / "data/raw/GSE115978/suppl/GSE115978_tpm.csv.gz"
    anno_path = ROOT / "data/raw/GSE115978/suppl/GSE115978_cell.annotations.csv.gz"
    out_dir = ROOT / "results/scrna"
    out_dir.mkdir(parents=True, exist_ok=True)
    if not module_path.exists() or not expr_path.exists() or not anno_path.exists():
        (out_dir / "module_scores_by_cell.tsv").write_text(
            "status\treason\nRESULT_PENDING\tMissing module table or GSE115978 files\n",
            encoding="utf-8",
        )
        return
    module = pd.read_csv(module_path, sep="\t")
    genes = set(module["gene"].dropna().astype(str))
    selected = []
    for chunk in pd.read_csv(expr_path, compression="gzip", chunksize=2000):
        gene_col = chunk.columns[0]
        chunk[gene_col] = chunk[gene_col].astype(str)
        keep = chunk[chunk[gene_col].isin(genes)]
        if not keep.empty:
            selected.append(keep)
    if not selected:
        (out_dir / "module_scores_by_cell.tsv").write_text(
            "status\treason\nRESULT_PENDING\tNo EcoNiche module genes found in GSE115978 TPM\n",
            encoding="utf-8",
        )
        return
    expr = pd.concat(selected, ignore_index=True)
    gene_col = expr.columns[0]
    expr = expr.set_index(gene_col)
    expr = expr.apply(pd.to_numeric, errors="coerce")
    annotations = pd.read_csv(anno_path, compression="gzip")
    annotations = annotations.rename(columns={"cells": "cell_id", "samples": "patient_id", "cell.types": "cell_type", "treatment.group": "treatment_group"})
    annotations = annotations.set_index("cell_id", drop=False)
    score_frames = []
    for state, state_module in module.groupby("state"):
        state_genes = [gene for gene in state_module["gene"].astype(str) if gene in expr.index]
        if not state_genes:
            continue
        values = expr.loc[state_genes].mean(axis=0).rename("module_score").reset_index()
        values = values.rename(columns={"index": "cell_id"})
        values["state"] = state
        score_frames.append(values)
    scores = pd.concat(score_frames, ignore_index=True) if score_frames else pd.DataFrame()
    scores = scores.merge(annotations.reset_index(drop=True), on="cell_id", how="left")
    scores.to_csv(out_dir / "module_scores_by_cell.tsv", sep="\t", index=False)
    patient = scores.groupby(["patient_id", "state"], dropna=False)["module_score"].mean().reset_index()
    patient.to_csv(out_dir / "module_scores_by_patient.tsv", sep="\t", index=False)
    enrichment = scores.groupby(["cell_type", "state"], dropna=False)["module_score"].agg(["mean", "median", "count"]).reset_index()
    enrichment.to_csv(out_dir / "cell_type_enrichment.tsv", sep="\t", index=False)
    print(f"Wrote scRNA module scores to {out_dir}")


if __name__ == "__main__":
    main()
