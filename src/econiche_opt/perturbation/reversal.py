from __future__ import annotations

import math
import re
from collections import defaultdict
from typing import Iterable

import pandas as pd
import requests


ENRICHR_ADDLIST_URL = "https://maayanlab.cloud/Enrichr/addList"
ENRICHR_ENRICH_URL = "https://maayanlab.cloud/Enrichr/enrich"
REVERSAL_COLUMNS = [
    "query_direction",
    "library",
    "rank",
    "perturbation_id",
    "perturbation_name",
    "cell_line",
    "timepoint",
    "dose",
    "pvalue",
    "zscore",
    "combined_score",
    "adjusted_pvalue",
    "overlapping_genes",
    "signature_genes_submitted",
    "reversal_score",
    "status",
    "interpretation",
]


def score_reversal(module_table: pd.DataFrame, perturbation_table: pd.DataFrame) -> pd.DataFrame:
    if module_table.empty or perturbation_table.empty:
        return pd.DataFrame(columns=["drug", "reversal_score", "status"])
    genes = set(module_table["gene"].astype(str))
    rows = []
    for drug, group in perturbation_table.groupby("drug"):
        overlap = group[group["gene"].astype(str).isin(genes)]
        score = -float(overlap.get("effect", pd.Series(dtype=float)).mean()) if not overlap.empty else 0.0
        rows.append({"drug": drug, "reversal_score": score, "status": "hypothesis_only"})
    return pd.DataFrame(rows).sort_values("reversal_score", ascending=False)


def _clean_gene(gene: object) -> str:
    return str(gene).strip().upper()


def _numeric_column(module: pd.DataFrame, name: str, default: float) -> pd.Series:
    if name in module.columns:
        return pd.to_numeric(module[name], errors="coerce").fillna(default)
    return pd.Series(default, index=module.index, dtype=float)


def build_reversal_signature(module_table: pd.DataFrame, top_genes: int = 150) -> dict[str, list[str]]:
    if module_table.empty or "gene" not in module_table.columns:
        return {"resistance_up": [], "resistance_down": []}
    module = module_table.copy()
    module["gene"] = module["gene"].map(_clean_gene)
    module = module[module["gene"].ne("")]
    if module.empty:
        return {"resistance_up": [], "resistance_down": []}

    direction = _numeric_column(module, "direction", 1.0)
    coefficient = _numeric_column(module, "coefficient", 1.0)
    selection_frequency = _numeric_column(module, "selection_frequency", 1.0).abs()
    module["resistance_weight"] = direction * coefficient
    if (module["resistance_weight"].abs() == 0).all():
        module["resistance_weight"] = direction
    module["rank_weight"] = module["resistance_weight"].abs() * selection_frequency.clip(lower=0.01)
    module = module.sort_values("rank_weight", ascending=False).drop_duplicates("gene")
    positive = module[module["resistance_weight"] > 0].head(top_genes)["gene"].tolist()
    negative = module[module["resistance_weight"] < 0].head(top_genes)["gene"].tolist()
    if not positive and not negative:
        genes = module.head(top_genes)["gene"].tolist()
        positive = genes
    return {"resistance_up": positive, "resistance_down": negative}


def parse_lincs_term(term: str) -> dict[str, str]:
    parts = str(term).split()
    perturbation_id = parts[0] if parts else str(term)
    cell_line = parts[1] if len(parts) > 1 else ""
    payload = " ".join(parts[2:]) if len(parts) > 2 else str(term)
    match = re.match(r"(?P<timepoint>\d+H)-(?P<drug>.+)-(?P<dose>[0-9.]+)$", payload)
    if match:
        drug = match.group("drug")
        timepoint = match.group("timepoint")
        dose = match.group("dose")
    else:
        drug = payload or str(term)
        timepoint = ""
        dose = ""
    return {
        "perturbation_id": perturbation_id,
        "perturbation_name": drug,
        "cell_line": cell_line,
        "timepoint": timepoint,
        "dose": dose,
    }


def _submit_enrichr_list(genes: Iterable[str], description: str, timeout: int = 30) -> int:
    gene_text = "\n".join(dict.fromkeys(_clean_gene(gene) for gene in genes if _clean_gene(gene)))
    if not gene_text:
        raise ValueError("no genes available for Enrichr query")
    response = requests.post(
        ENRICHR_ADDLIST_URL,
        files={"list": (None, gene_text), "description": (None, description)},
        timeout=timeout,
    )
    response.raise_for_status()
    return int(response.json()["userListId"])


def query_enrichr_library(
    genes: Iterable[str],
    library: str,
    query_direction: str,
    max_results: int = 100,
    timeout: int = 30,
) -> pd.DataFrame:
    unique_genes = list(dict.fromkeys(_clean_gene(gene) for gene in genes if _clean_gene(gene)))
    if not unique_genes:
        return pd.DataFrame(columns=REVERSAL_COLUMNS)
    user_list_id = _submit_enrichr_list(unique_genes, f"EcoNiche-Opt {query_direction}", timeout=timeout)
    response = requests.get(
        ENRICHR_ENRICH_URL,
        params={"userListId": user_list_id, "backgroundType": library},
        timeout=timeout,
    )
    response.raise_for_status()
    rows = []
    for entry in response.json().get(library, [])[:max_results]:
        parsed = parse_lincs_term(entry[1])
        adjusted_pvalue = float(entry[6]) if len(entry) > 6 and entry[6] is not None else math.nan
        combined_score = float(entry[4]) if len(entry) > 4 and entry[4] is not None else math.nan
        if math.isfinite(adjusted_pvalue) and adjusted_pvalue > 0:
            p_component = -math.log10(adjusted_pvalue)
        else:
            p_component = 0.0
        score = (combined_score if math.isfinite(combined_score) else 0.0) + p_component
        overlap = entry[5] if len(entry) > 5 else []
        rows.append(
            {
                "query_direction": query_direction,
                "library": library,
                "rank": int(entry[0]),
                **parsed,
                "pvalue": float(entry[2]) if len(entry) > 2 and entry[2] is not None else math.nan,
                "zscore": float(entry[3]) if len(entry) > 3 and entry[3] is not None else math.nan,
                "combined_score": combined_score,
                "adjusted_pvalue": adjusted_pvalue,
                "overlapping_genes": ";".join(overlap) if isinstance(overlap, list) else str(overlap),
                "signature_genes_submitted": len(unique_genes),
                "reversal_score": score,
                "status": "LINCS_Enrichr_hypothesis",
                "interpretation": "signature-reversal hypothesis only; not a clinical recommendation",
            }
        )
    return pd.DataFrame(rows, columns=REVERSAL_COLUMNS)


def query_lincs_reversal(
    module_table: pd.DataFrame,
    top_genes: int = 150,
    max_results: int = 100,
    timeout: int = 30,
) -> pd.DataFrame:
    signature = build_reversal_signature(module_table, top_genes=top_genes)
    frames = []
    if signature["resistance_up"]:
        frames.append(
            query_enrichr_library(
                signature["resistance_up"],
                library="LINCS_L1000_Chem_Pert_down",
                query_direction="downregulate_resistance_up_genes",
                max_results=max_results,
                timeout=timeout,
            )
        )
    if signature["resistance_down"]:
        frames.append(
            query_enrichr_library(
                signature["resistance_down"],
                library="LINCS_L1000_Chem_Pert_up",
                query_direction="upregulate_resistance_down_genes",
                max_results=max_results,
                timeout=timeout,
            )
        )
    frames = [frame for frame in frames if not frame.empty]
    if not frames:
        return pd.DataFrame(columns=REVERSAL_COLUMNS)
    return pd.concat(frames, ignore_index=True).sort_values("reversal_score", ascending=False)


def _gene_union(values: Iterable[str]) -> str:
    genes: set[str] = set()
    for value in values:
        for gene in str(value).split(";"):
            gene = _clean_gene(gene)
            if gene:
                genes.add(gene)
    return ";".join(sorted(genes))


def aggregate_reversal_candidates(reversal_results: pd.DataFrame, top_n: int = 100) -> pd.DataFrame:
    columns = [
        "perturbation_id",
        "perturbation_name",
        "target_gene",
        "mechanism",
        "reversal_score",
        "depmap_score",
        "dgidb_evidence",
        "target_state",
        "priority_score",
        "status",
        "interpretation",
    ]
    if reversal_results.empty:
        return pd.DataFrame(columns=columns)
    rows = []
    for name, group in reversal_results.groupby("perturbation_name", dropna=False):
        score = float(group["reversal_score"].sum())
        libraries = ";".join(sorted(group["library"].dropna().astype(str).unique()))
        ids = ";".join(sorted(group["perturbation_id"].dropna().astype(str).unique()))
        directions = ";".join(sorted(group["query_direction"].dropna().astype(str).unique()))
        rows.append(
            {
                "perturbation_id": ids,
                "perturbation_name": name,
                "target_gene": _gene_union(group["overlapping_genes"]),
                "mechanism": f"LINCS Enrichr signature reversal ({directions}; {libraries})",
                "reversal_score": score,
                "depmap_score": "",
                "dgidb_evidence": "",
                "target_state": "EcoNiche module-derived resistance signature",
                "priority_score": score,
                "status": "hypothesis_only",
                "interpretation": "candidate perturbation hypothesis only; requires experimental validation",
            }
        )
    return pd.DataFrame(rows, columns=columns).sort_values("priority_score", ascending=False).head(top_n)


def combine_with_dgidb(priority: pd.DataFrame, dgidb_hits: pd.DataFrame) -> pd.DataFrame:
    if priority.empty or dgidb_hits.empty or "drug_name" not in dgidb_hits.columns:
        return priority
    evidence: dict[str, list[str]] = defaultdict(list)
    for _, row in dgidb_hits.iterrows():
        drug = str(row.get("drug_name", "")).strip().upper()
        if not drug:
            continue
        gene = str(row.get("target_gene", "")).strip().upper()
        interaction_value = row.get("interaction_types", "")
        sources_value = row.get("sources", "")
        interaction = "" if pd.isna(interaction_value) else str(interaction_value).strip()
        sources = "" if pd.isna(sources_value) else str(sources_value).strip()
        evidence[drug].append(f"{gene}:{interaction or 'interaction'}:{sources}")
    updated = priority.copy()
    merged_evidence = []
    for name in updated["perturbation_name"]:
        key = str(name).strip().upper()
        matches = list(evidence.get(key, []))
        for drug, items in evidence.items():
            if key and drug and (key in drug or drug in key):
                matches.extend(items)
        merged_evidence.append(";".join(sorted(set(matches))))
    updated["dgidb_evidence"] = merged_evidence
    updated["priority_score"] = updated["priority_score"].astype(float) + updated["dgidb_evidence"].ne("").astype(int) * 10.0
    return updated.sort_values("priority_score", ascending=False)
