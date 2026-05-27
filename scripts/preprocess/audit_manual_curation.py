from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.registry import load_registry, normalize_access_status


NEWLY_CURATED = {"GSE67501", "GSE93157", "GSE140901", "GSE145996", "GSE165252"}


def _first_present(row: pd.Series, columns: list[str], default: object = pd.NA) -> object:
    for col in columns:
        if col in row and pd.notna(row[col]) and str(row[col]).strip():
            return row[col]
    return default


def _counts(series: pd.Series) -> str:
    counts = series.value_counts(dropna=False).to_dict()
    normalized = {str(k): int(v) for k, v in counts.items()}
    return json.dumps(normalized, ensure_ascii=True, sort_keys=True)


def evidence_field(accession: str) -> str:
    evidence = {
        "GSE67501": "description (RCC sample ID); characteristics_ch1 response to anti-PD-1 nivolumab; GEO series design states pre-treatment tumors",
        "GSE93157": "characteristics_ch1 best.resp/drug/biopsy; expression columns matched by GEO sample order",
        "GSE140901": "characteristics_ch1 best_response/clinical_benefit_response; title Sxx matched to expression column Sxx",
        "GSE145996": "MDPI Cancers supplementary workbook: SuppTable1 Patient Number/Best Response plus SuppTable2 Tumor ID number/RNASEQ; expression columns matched by MB tumor ID",
        "GSE165252": "characteristics_ch1 response; title sample_N baseline/on_treatment/resection; expression columns matched by GEO sample order",
        "GSE91061": "characteristics_ch1 visit_pre_or_on_treatment and response",
        "GSE78220": "clinical/FPKM workbook metadata response fields",
        "GSE115821": "series/sample metadata response fields",
        "GSE136961": "GEO sample title D/N durable clinical benefit encoding",
        "GSE168204": "series/sample metadata response fields",
        "GSE176307": "series/sample metadata response fields",
    }
    return evidence.get(accession, "response_raw harmonized from available sample metadata")


def decision_for(accession: str, cohort: dict, qc_row: pd.Series | None, meta: pd.DataFrame, proc: pd.DataFrame) -> tuple[str, str, str]:
    n_labels = int(proc["label"].notna().sum()) if "label" in proc else 0
    n_classes = int(proc["label"].nunique(dropna=True)) if "label" in proc else 0
    status = str(qc_row.get("status", "")) if qc_row is not None else ""
    role = str(cohort.get("role", ""))
    notes = str(cohort.get("notes", ""))
    uses = set(cohort.get("uses", []) or [])

    if n_labels >= 4 and n_classes >= 2:
        if accession == "GSE145996":
            return (
                "usable_melanoma_response",
                "Paper supplement maps RNA-seq tumor IDs to patient-level RECIST best response; expression columns match MB tumor IDs.",
                "Include in melanoma endpoint-stratified analysis with explicit supplement citation/evidence.",
            )
        if accession == "GSE165252":
            return (
                "usable_secondary_dynamic_confounded",
                "Binary pathologic response labels and baseline samples recovered, but therapy is atezolizumab plus chemoradiotherapy; keep as secondary/context-specific validation.",
                "Use in pan-cancer/secondary analyses; report confounding explicitly.",
            )
        if accession in {"GSE67501", "GSE93157", "GSE140901"}:
            return (
                "usable_pan_cancer_response",
                "Sample-level response labels, therapy, and expression sample matching are traceable from GEO metadata/supplementary files.",
                "Include in real processed data and stratified pan-cancer LODO.",
            )
        return (
            "usable_response_benchmark",
            "Sample-level response labels and expression matrix are available after harmonization.",
            "Keep in primary response benchmark if endpoint/setting matches registry role.",
        )

    if accession == "GSE145996":
        return (
            "hold_missing_sample_level_response_map",
            "GEO/paper summary states RECIST responder/non-responder groups, but no local supplement-derived per-sample response map was available at audit time.",
            "Download the MDPI Cancers supplementary workbook and map Sample_XX_MB_YYYY columns to Tumor ID number and Best Response.",
        )
    if accession == "GSE121810":
        return (
            "exclude_primary_response",
            "GEO metadata describes neoadjuvant/adjuvant pembrolizumab treatment-arm and survival setting in recurrent glioblastoma, not per-sample RECIST response.",
            "Use only as context-specific survival/mechanism cohort unless a traceable binary response table is added.",
        )
    if accession == "GSE244982" or "mechanism" in uses or "mechanism" in role:
        return (
            "exclude_primary_response",
            "Registry and sample metadata describe acquired resistance/progression mechanism groups rather than primary binary ICB response.",
            "Keep for mechanism analyses; exclude from primary response LODO.",
        )
    if accession == "GSE183924" or "survival" in uses or "survival" in role or "Survival" in notes:
        return (
            "survival_only_not_response",
            "Clinical workbook contains relapse/RFS/OS fields but no CR/PR/SD/PD or DCB/NDB response endpoint.",
            "Use in survival analyses, not response classification.",
        )
    if "unusable_expression_shape" in status:
        return (
            "expression_parser_required",
            "Expression file did not parse into a usable sample-by-gene matrix under the current reader.",
            "Write a cohort-specific parser or platform annotation before use.",
        )
    if "missing_expression" in status:
        return (
            "missing_expression_file",
            "No public expression matrix is present locally for this accession.",
            "Download/verify expression data before curation.",
        )
    return (
        "hold_needs_manual_curation",
        "Current local files do not provide enough matched binary response labels for training.",
        "Search paper supplements/clinical workbooks and add sample-level evidence before use.",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", default="config/data_registry.yml")
    parser.add_argument("--metadata", default="data/metadata/metadata_harmonized.tsv")
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--qc-report", default="tables/expression_qc_report_manual_curation.tsv")
    parser.add_argument("--sample-out", default="data/metadata/manual_response_curation.tsv")
    parser.add_argument("--audit-out", default="results/curation/manual_curation_audit.tsv")
    parser.add_argument("--markdown-out", default="results/curation/MANUAL_CURATION_AUDIT.md")
    args = parser.parse_args()

    registry = load_registry(ROOT / args.registry)
    cohorts = [
        c
        for c in registry.get("cohorts", [])
        if str(c.get("accession", "")).startswith("GSE")
        and normalize_access_status(c.get("access")) == "public"
        and "scRNA" not in str(c.get("platform", ""))
    ]
    cohort_by_accession = {c["accession"]: c for c in cohorts}

    metadata = pd.read_csv(ROOT / args.metadata, sep="\t", low_memory=False)
    qc = pd.read_csv(ROOT / args.qc_report, sep="\t") if (ROOT / args.qc_report).exists() else pd.DataFrame()
    qc_by_accession = {str(r["cohort"]): r for _, r in qc.iterrows()} if not qc.empty else {}

    sample_rows: list[dict[str, object]] = []
    audit_rows: list[dict[str, object]] = []
    processed_dir = ROOT / args.processed_dir
    for accession, cohort in cohort_by_accession.items():
        meta = metadata[metadata["accession"].astype(str) == accession].copy()
        proc_path = processed_dir / f"{accession}.metadata.tsv"
        proc = pd.read_csv(proc_path, sep="\t") if proc_path.exists() else pd.DataFrame()
        qc_row = qc_by_accession.get(accession)
        decision, reason, next_action = decision_for(accession, cohort, qc_row, meta, proc)
        label_counts = _counts(proc["label"]) if "label" in proc else "{}"
        response_counts = _counts(meta["response_raw"]) if "response_raw" in meta else "{}"
        audit_rows.append(
            {
                "accession": accession,
                "name": cohort.get("name", ""),
                "cancer_type": cohort.get("cancer_type", ""),
                "therapy_registry": cohort.get("therapy", ""),
                "role": cohort.get("role", ""),
                "source_file": qc_row.get("source_file", "") if qc_row is not None else "",
                "qc_status": qc_row.get("status", "") if qc_row is not None else "",
                "n_metadata_rows": len(meta),
                "n_processed_samples": len(proc),
                "n_genes": int(qc_row.get("n_genes", 0)) if qc_row is not None and pd.notna(qc_row.get("n_genes", pd.NA)) else 0,
                "response_raw_counts": response_counts,
                "processed_label_counts": label_counts,
                "evidence_field": evidence_field(accession),
                "decision": decision,
                "reason": reason,
                "next_action": next_action,
                "geo_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
            }
        )
        if proc.empty or "label" not in proc:
            continue
        labeled = proc[proc["label"].notna()].copy()
        for _, row in labeled.iterrows():
            sample_rows.append(
                {
                    "accession": accession,
                    "cohort_name": cohort.get("name", ""),
                    "curation_tier": "manual_curated_2026_05_05" if accession in NEWLY_CURATED else "previously_parsed_or_pipeline_curated",
                    "sample_id": row.get("sample_id", ""),
                    "expression_sample_id": row.get("expression_sample_id", row.get("sample_id", "")),
                    "geo_accession": row.get("geo_accession", ""),
                    "patient_id": row.get("patient_id", row.get("patient_id_raw", "")),
                    "cancer_type": cohort.get("cancer_type", ""),
                    "therapy": _first_present(row, ["therapy", "treatment", "drug", "io_therapy"], cohort.get("therapy", "")),
                    "timepoint": row.get("timepoint", ""),
                    "response_raw": row.get("response_raw", ""),
                    "label": int(float(row.get("label"))),
                    "label_definition": "0=responder/clinical benefit; 1=non-responder/no benefit; SD treated as non-response for primary_recist",
                    "match_method": row.get("match_method", ""),
                    "evidence_source_file": qc_row.get("source_file", "") if qc_row is not None else "",
                    "evidence_field": evidence_field(accession),
                    "source_characteristics": row.get("characteristics_ch1", ""),
                    "geo_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
                }
            )

    sample_out = ROOT / args.sample_out
    audit_out = ROOT / args.audit_out
    markdown_out = ROOT / args.markdown_out
    sample_out.parent.mkdir(parents=True, exist_ok=True)
    audit_out.parent.mkdir(parents=True, exist_ok=True)
    markdown_out.parent.mkdir(parents=True, exist_ok=True)
    sample_df = pd.DataFrame(sample_rows)
    audit_df = pd.DataFrame(audit_rows)
    sample_df.to_csv(sample_out, sep="\t", index=False)
    audit_df.to_csv(audit_out, sep="\t", index=False)

    recovered = audit_df[audit_df["accession"].isin(sorted(NEWLY_CURATED))]
    usable = audit_df[audit_df["decision"].str.startswith("usable", na=False)]
    held = audit_df[~audit_df["decision"].str.startswith("usable", na=False)]
    lines = [
        "# Manual ICB Response Curation Audit",
        "",
        f"- Sample-level evidence table: `{sample_out.relative_to(ROOT)}` ({len(sample_df)} labeled samples)",
        f"- Cohort audit table: `{audit_out.relative_to(ROOT)}` ({len(audit_df)} cohorts)",
        f"- Newly recovered cohorts: {', '.join(sorted(NEWLY_CURATED))}",
        f"- Usable response cohorts after curation: {len(usable)}",
        "",
        "## Newly recovered cohorts",
        "",
        "| Accession | Samples | Genes | Labels | Decision | Evidence |",
        "|---|---:|---:|---|---|---|",
    ]
    for _, row in recovered.iterrows():
        lines.append(
            f"| {row['accession']} | {row['n_processed_samples']} | {row['n_genes']} | "
            f"`{row['processed_label_counts']}` | {row['decision']} | {row['evidence_field']} |"
        )
    lines.extend(["", "## Held out / excluded cohorts", "", "| Accession | Decision | Reason |", "|---|---|---|"])
    for _, row in held.iterrows():
        if row["n_processed_samples"] > 0:
            continue
        lines.append(f"| {row['accession']} | {row['decision']} | {row['reason']} |")
    lines.extend(
        [
            "",
            "## Label rule",
            "",
            "The response endpoint is harmonized as `0=responder/clinical benefit` and `1=non-responder/no benefit`. For the primary RECIST-style endpoint, `SD` is conservatively grouped with non-response unless an explicit durable-clinical-benefit endpoint is used.",
        ]
    )
    markdown_out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote sample evidence to {sample_out}")
    print(f"Wrote cohort audit to {audit_out}")
    print(f"Wrote markdown audit to {markdown_out}")


if __name__ == "__main__":
    main()
