from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]


def exists(path: str) -> bool:
    return (ROOT / path).exists()


def file_size(path: str) -> int:
    p = ROOT / path
    return p.stat().st_size if p.exists() and p.is_file() else 0


def main() -> None:
    rows = [
        {
            "accession": "GSE145996",
            "priority": "high",
            "status": "COMPLETED_INCLUDED",
            "evidence_strength": "strong_sample_level",
            "local_expression": "data/processed/bulk/GSE145996.expr.tsv",
            "local_response_map": "data/processed/bulk/GSE145996.metadata.tsv",
            "local_external_evidence": "data/external/GSE145996/cancers-12-01943-s001/cancers-836356-supplemenraty-final/cancers-836356 - supplementary - revised.xlsx",
            "source_url": "https://res.mdpi.com/d_attachment/cancers/cancers-12-01943/article_deploy/cancers-12-01943-s001.zip",
            "completed_action": "Downloaded MDPI supplement; mapped expression Sample_XX_MB_YYYY to SuppTable2 Tumor ID number and SuppTable1 Best Response.",
            "remaining_gap": "Small RNA-seq subset only; SD handling must remain endpoint-stratified.",
            "next_action": "Use in melanoma/endpoint-stratified analysis; cite exact supplement table evidence.",
        },
        {
            "accession": "PRJEB23709",
            "priority": "high",
            "status": "CLINICAL_MAP_COMPLETED_EXPRESSION_PENDING",
            "evidence_strength": "strong_clinical_map_no_expression_matrix",
            "local_expression": "",
            "local_response_map": "data/external/PRJEB23709/gide_prjeb23709_response_map.tsv",
            "local_external_evidence": "data/external/PRJEB23709/mmc2.xlsx;data/external/PRJEB23709/mmc3.xlsx;data/raw/PRJEB23709/ena_run_metadata.tsv",
            "source_url": "https://ars.els-cdn.com/content/image/1-s2.0-S1535610819300376-mmc2.xlsx;https://ars.els-cdn.com/content/image/1-s2.0-S1535610819300376-mmc3.xlsx",
            "completed_action": "Downloaded Gide supplementary clinical tables and joined them to ENA sample_title/run metadata for 91 RNA-seq runs.",
            "remaining_gap": "No processed expression matrix is local; FASTQ quantification or a validated processed matrix is required before modeling.",
            "next_action": "Either run Salmon/STAR quantification for PRE samples or locate an author/repository processed count matrix; then restrict primary analysis to PRE samples and split anti-PD1 monotherapy vs combo.",
        },
        {
            "accession": "IMvigor210",
            "priority": "high",
            "status": "EXPRESSION_AND_CLINICAL_PENDING_ENVIRONMENT_BLOCKED",
            "evidence_strength": "known_public_package_not_local",
            "local_expression": "",
            "local_response_map": "",
            "local_external_evidence": "data/raw/IMvigor210/status.tsv",
            "source_url": "https://bioconductor.org/packages/IMvigor210CoreBiologies/",
            "completed_action": "Recorded access blocker; Rscript/Bioconductor package execution is unavailable in this environment.",
            "remaining_gap": "Need expression matrix plus response/therapy/sample annotations from Bioconductor/ExperimentHub or a verified mirror.",
            "next_action": "Install/enable Rscript or fetch package RData/RDS assets via a direct verified source; only then add as urothelial anti-PDL1 validation.",
        },
        {
            "accession": "GSE123728",
            "priority": "medium",
            "status": "RAW_AND_RESPONSE_PENDING",
            "evidence_strength": "not_local",
            "local_expression": "",
            "local_response_map": "",
            "local_external_evidence": "",
            "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123728",
            "completed_action": "Registry candidate identified; no local raw directory or expression matrix found.",
            "remaining_gap": "Need verify whether tumor bulk expression and binary ICB response exist at sample level.",
            "next_action": "Download GEO supplementary files and sample metadata; include only if expression columns can be matched to response evidence.",
        },
        {
            "accession": "GSE165745",
            "priority": "medium",
            "status": "RAW_AND_RESPONSE_PENDING",
            "evidence_strength": "not_local",
            "local_expression": "",
            "local_response_map": "",
            "local_external_evidence": "",
            "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165745",
            "completed_action": "Registry candidate identified; no local raw directory or expression matrix found.",
            "remaining_gap": "Need verify expression modality, baseline status, therapy, and binary response.",
            "next_action": "Download GEO supplementary files and clinical/sample annotations; reject if only single-cell/mechanism without patient response.",
        },
        {
            "accession": "GSE122220",
            "priority": "medium",
            "status": "RAW_AND_RESPONSE_PENDING",
            "evidence_strength": "not_local",
            "local_expression": "",
            "local_response_map": "",
            "local_external_evidence": "",
            "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122220",
            "completed_action": "Registry candidate identified; no local raw directory or expression matrix found.",
            "remaining_gap": "Need verify sample-level expression and response labels.",
            "next_action": "Download GEO supplementary files and audit sample-to-patient/response mapping before use.",
        },
        {
            "accession": "GSE93157",
            "priority": "medium",
            "status": "INCLUDED_ORDER_MATCH_NEEDS_EXTRA_SUPPLEMENT_CONFIRMATION",
            "evidence_strength": "usable_but_order_match",
            "local_expression": "data/processed/bulk/GSE93157.expr.tsv",
            "local_response_map": "data/processed/bulk/GSE93157.metadata.tsv",
            "local_external_evidence": "data/metadata/manual_response_curation.tsv",
            "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93157",
            "completed_action": "Included through GEO sample-order expression/metadata match with best response fields.",
            "remaining_gap": "A supplementary sample map would further harden the expression-column to metadata match.",
            "next_action": "Search article supplement for a NanoString sample sheet; retain current cohort but flag order-match evidence in Methods.",
        },
        {
            "accession": "GSE165252",
            "priority": "medium",
            "status": "INCLUDED_SECONDARY_CONFOUNDED_ORDER_MATCH",
            "evidence_strength": "usable_secondary_but_order_match",
            "local_expression": "data/processed/bulk/GSE165252.expr.tsv",
            "local_response_map": "data/processed/bulk/GSE165252.metadata.tsv",
            "local_external_evidence": "data/metadata/manual_response_curation.tsv",
            "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165252",
            "completed_action": "Included as secondary transfer cohort with baseline samples and response labels.",
            "remaining_gap": "Treatment combines chemoradiotherapy and atezolizumab; expression columns are GEO-order matched.",
            "next_action": "Keep out of primary melanoma/pan-cancer headline; use as secondary/confounded transfer validation.",
        },
    ]
    frame = pd.DataFrame(rows)
    for column in ["local_expression", "local_response_map", "local_external_evidence"]:
        frame[f"{column}_exists"] = frame[column].map(lambda value: all(exists(part) for part in str(value).split(";") if part) if str(value) else False)
    for column in ["local_expression", "local_response_map"]:
        frame[f"{column}_bytes"] = frame[column].map(lambda value: file_size(str(value)) if str(value) else 0)

    out_dir = ROOT / "results/curation"
    out_dir.mkdir(parents=True, exist_ok=True)
    tsv = out_dir / "external_cohort_gap_audit.tsv"
    md = out_dir / "EXTERNAL_COHORT_GAP_AUDIT.md"
    frame.to_csv(tsv, sep="\t", index=False)
    lines = [
        "# External Cohort Gap Audit",
        "",
        f"- Audit table: `{tsv.relative_to(ROOT)}`",
        "- Goal: prioritize cohorts that can improve melanoma primary and pan-cancer transfer without weakening traceability.",
        "",
        "| Accession | Status | Evidence | Remaining gap | Next action |",
        "|---|---|---|---|---|",
    ]
    for _, row in frame.iterrows():
        lines.append(
            f"| {row['accession']} | {row['status']} | {row['evidence_strength']} | {row['remaining_gap']} | {row['next_action']} |"
        )
    md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {tsv}")
    print(f"Wrote {md}")


if __name__ == "__main__":
    main()
