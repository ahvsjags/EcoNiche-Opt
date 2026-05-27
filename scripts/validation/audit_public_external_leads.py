from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


LEADS: list[dict[str, str]] = [
    {
        "lead_id": "LIU_MGSP_PHS000452",
        "source_database": "dbGaP/TIGER/cBioPortal",
        "candidate_accession": "phs000452; mel_dfci_2019",
        "title": "Liu/MGSP pretreatment melanoma anti-PD-1 tumor RNA-seq",
        "cancer_type": "melanoma",
        "therapy": "anti-PD-1",
        "sample_type": "tumor tissue",
        "platform": "bulk RNA-seq",
        "timepoint": "pretreatment",
        "response_label": "CR/PR/MR/R versus SD/PD/NR/NDB; strict RECIST excludes MR/SD downstream",
        "independence_role": "strict_external_source_crosscheck",
        "eligibility_status": "eligible_processed_duplicate_sensitive",
        "evidence_url": "https://www.cbioportal.org/study/summary?id=mel_dfci_2019",
        "evidence_note": "Public cBioPortal expression and clinical attributes are processed locally; TIGER phs000452 is a related processed source. Do not count duplicate Liu-derived versions as independent validations.",
        "next_action": "Use only one Liu-derived source per locked external analysis; keep duplicate-source boundaries explicit.",
    },
    {
        "lead_id": "GIDE_PRJEB23709",
        "source_database": "ENA/TIGER/cBioPortal-iAtlas",
        "candidate_accession": "PRJEB23709",
        "title": "Gide melanoma anti-PD-1 and anti-PD-1 plus anti-CTLA-4 pretreatment RNA-seq",
        "cancer_type": "melanoma",
        "therapy": "anti-PD-1; combination anti-PD-1 plus anti-CTLA-4",
        "sample_type": "tumor tissue",
        "platform": "bulk RNA-seq",
        "timepoint": "pretreatment",
        "response_label": "RECIST-like response harmonized from supplementary clinical tables",
        "independence_role": "primary_discovery_or_external_if_frozen",
        "eligibility_status": "eligible_processed",
        "evidence_url": "https://www.ebi.ac.uk/ena/browser/view/PRJEB23709",
        "evidence_note": "Processed monotherapy and combination pretreatment cohorts exist locally. iAtlas/cBioPortal versions should be treated as duplicate source cross-checks.",
        "next_action": "Keep monotherapy as a primary melanoma component or freeze before any external use; do not mix duplicate iAtlas data as independent evidence.",
    },
    {
        "lead_id": "RIAZ_GSE91061",
        "source_database": "GEO",
        "candidate_accession": "GSE91061",
        "title": "Riaz nivolumab metastatic melanoma RNA-seq",
        "cancer_type": "melanoma",
        "therapy": "anti-PD-1 nivolumab",
        "sample_type": "tumor tissue",
        "platform": "bulk RNA-seq",
        "timepoint": "pretreatment and on-treatment; current primary uses pretreatment only",
        "response_label": "RECIST-like response harmonized at patient level",
        "independence_role": "primary_discovery_internal_lodo",
        "eligibility_status": "eligible_processed",
        "evidence_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061",
        "evidence_note": "Processed locally and used under patient-level LODO. cBioPortal/iAtlas Riaz is a duplicate source cross-check.",
        "next_action": "Maintain pretreatment-only filtering and patient-level fold boundaries.",
    },
    {
        "lead_id": "HUGO_GSE78220",
        "source_database": "GEO",
        "candidate_accession": "GSE78220",
        "title": "Hugo pretreatment melanoma anti-PD-1 RNA-seq",
        "cancer_type": "melanoma",
        "therapy": "anti-PD-1 pembrolizumab/nivolumab",
        "sample_type": "tumor tissue",
        "platform": "bulk RNA-seq",
        "timepoint": "pretreatment",
        "response_label": "irRECIST/response harmonized from publication and GEO records",
        "independence_role": "primary_discovery_internal_lodo",
        "eligibility_status": "eligible_processed",
        "evidence_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78220",
        "evidence_note": "Processed locally; small classic anti-PD-1 melanoma cohort.",
        "next_action": "Keep as part of high-evidence primary melanoma LODO or freeze before external scoring.",
    },
    {
        "lead_id": "MGH_GSE115821",
        "source_database": "GEO",
        "candidate_accession": "GSE115821",
        "title": "MGH melanoma anti-PD-1 pre/on-treatment RNA-seq",
        "cancer_type": "melanoma",
        "therapy": "anti-PD-1",
        "sample_type": "tumor tissue",
        "platform": "bulk RNA-seq",
        "timepoint": "pretreatment and on-treatment",
        "response_label": "RECIST response after manual filtering",
        "independence_role": "secondary_external_or_sensitivity",
        "eligibility_status": "eligible_processed_small",
        "evidence_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115821",
        "evidence_note": "Processed locally; small sample size and mixed timing limit its strength as a sole strict external set.",
        "next_action": "Use in endpoint/timing sensitivity analyses with baseline-only filtering.",
    },
    {
        "lead_id": "MGH_GSE168204",
        "source_database": "GEO",
        "candidate_accession": "GSE168204",
        "title": "New MGH melanoma anti-PD-1 pre/on-treatment RNA-seq",
        "cancer_type": "melanoma",
        "therapy": "anti-PD-1/PD-L1",
        "sample_type": "tumor tissue",
        "platform": "bulk RNA-seq",
        "timepoint": "pretreatment and on-treatment",
        "response_label": "RECIST response after manual filtering",
        "independence_role": "secondary_external_or_sensitivity",
        "eligibility_status": "eligible_processed_small",
        "evidence_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168204",
        "evidence_note": "Processed locally; related literature describes 27 pre/on-treatment tumor specimens with exclusions for missing RECIST/RNA-seq/duplicates.",
        "next_action": "Use as a small melanoma sensitivity layer, not as the only top-tier external validation.",
    },
    {
        "lead_id": "GSE145996",
        "source_database": "GEO",
        "candidate_accession": "GSE145996",
        "title": "Melanoma tumor RNA-seq cohort with limited response-label evidence",
        "cancer_type": "melanoma",
        "therapy": "ICB context",
        "sample_type": "tumor tissue",
        "platform": "bulk RNA-seq",
        "timepoint": "patient tumor specimen; pretreatment status requires curated evidence",
        "response_label": "curated binary response label used in current strict external stress test",
        "independence_role": "strict_external_stress_test",
        "eligibility_status": "eligible_processed_evidence_sensitive",
        "evidence_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145996",
        "evidence_note": "Processed locally but GEO sample annotations emphasize tumor specimens and need stronger source/timepoint/RECIST mapping.",
        "next_action": "Keep as strict stress-test support; harden sample-level response and baseline evidence before headline use.",
    },
    {
        "lead_id": "LEE_RIZOS_EGAS00001001552",
        "source_database": "EGA",
        "candidate_accession": "EGAS00001001552; EGAD00001005738",
        "title": "Lee/Rizos melanoma immune checkpoint blockade RNA-seq",
        "cancer_type": "melanoma",
        "therapy": "immune checkpoint blockade",
        "sample_type": "tumor tissue",
        "platform": "bulk RNA-seq",
        "timepoint": "pre/on-treatment in controlled dataset",
        "response_label": "publication-level RECIST response; raw data controlled",
        "independence_role": "highest_value_controlled_external_candidate",
        "eligibility_status": "controlled_access_required",
        "evidence_url": "https://ega-archive.org/studies/EGAS00001001552",
        "evidence_note": "EGA lists EGAD00001005738 as 79 RNA-seq samples from 56 melanoma patients who underwent immune checkpoint blockade immunotherapy.",
        "next_action": "Request EGA access and run registered preprocessing; do not substitute or fabricate expression matrices.",
    },
    {
        "lead_id": "ABRIL_RODRIGUEZ_PHS001919",
        "source_database": "dbGaP",
        "candidate_accession": "phs001919.v1.p1",
        "title": "Abril-Rodriguez melanoma anti-PD-1 RNA-seq immune-exclusion cohort",
        "cancer_type": "melanoma",
        "therapy": "anti-PD-1/checkpoint blockade",
        "sample_type": "tumor biopsy",
        "platform": "bulk RNA-seq",
        "timepoint": "treated melanoma biopsy; baseline and response timing require controlled clinical files",
        "response_label": "response and immune-infiltration variables require dbGaP phenotype access",
        "independence_role": "high_value_controlled_external_candidate",
        "eligibility_status": "controlled_access_required",
        "evidence_url": "https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001919.v1.p1",
        "evidence_note": "dbGaP describes melanoma patients treated with human anti-PD-1, RNA sequencing on biopsies, and immune-exclusion/non-response context; expression and phenotype files require authorized access.",
        "next_action": "Request dbGaP access, verify pretreatment tumor timing and RECIST response variables, then run the registered locked external scorer.",
    },
    {
        "lead_id": "MGH_HACOHEN_PHS002683",
        "source_database": "dbGaP",
        "candidate_accession": "phs002683.v1.p1",
        "title": "Combined tumor and immune signals in melanoma checkpoint blockade",
        "cancer_type": "melanoma",
        "therapy": "anti-CTLA4, anti-PD1, anti-PDL1, or combination checkpoint blockade",
        "sample_type": "tumor biopsy",
        "platform": "bulk RNA-seq plus WES",
        "timepoint": "biopsy timepoint available after controlled phenotype access",
        "response_label": "response and survival outcomes require dbGaP phenotype access",
        "independence_role": "high_value_controlled_external_candidate",
        "eligibility_status": "controlled_access_required",
        "evidence_url": "https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002683.v1.p1",
        "evidence_note": "dbGaP reports 67 checkpoint-blockade melanoma patients and aggregated bulk RNA-seq data for 178 patients, with transcriptomic factors associated with response and survival.",
        "next_action": "Request dbGaP access, isolate pretreatment anti-PD-1 or PD-1-like tumor RNA-seq samples, and keep this cohort locked outside model development.",
    },
    {
        "lead_id": "GSE123728",
        "source_database": "GEO",
        "candidate_accession": "GSE123728",
        "title": "Neoadjuvant single-dose PD-1 blockade melanoma NanoString panel",
        "cancer_type": "melanoma",
        "therapy": "neoadjuvant anti-PD-1",
        "sample_type": "tumor tissue",
        "platform": "NanoString custom panel",
        "timepoint": "pre and post samples",
        "response_label": "clinical/response characteristics; requires endpoint harmonization",
        "independence_role": "panel_transfer_candidate",
        "eligibility_status": "panel_transfer_metadata_pending",
        "evidence_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123728",
        "evidence_note": "GEO lists 24 samples, a custom NanoString panel, pre/post labels, and neoadjuvant PD-1 context. This is not bulk RNA-seq and is best suited to panel-transfer sensitivity.",
        "next_action": "Curate pre-only samples and response source if used; do not count as strict bulk RNA-seq external validation.",
    },
    {
        "lead_id": "GSE165745",
        "source_database": "GEO",
        "candidate_accession": "GSE165745",
        "title": "Metastatic melanoma prior to anti-PD-1 NanoString Wnt pathway panel",
        "cancer_type": "melanoma",
        "therapy": "anti-PD-1 pembrolizumab/nivolumab",
        "sample_type": "FFPE tumor tissue",
        "platform": "NanoString nCounter Vantage 3D Wnt panel plus custom immune markers",
        "timepoint": "prior to anti-PD-1",
        "response_label": "sample labels encode responder/nonresponder",
        "independence_role": "panel_transfer_candidate",
        "eligibility_status": "panel_transfer_public",
        "evidence_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165745",
        "evidence_note": "GEO lists 24 pretreatment FFPE tumor samples with responder/nonresponder labels and a 204-gene NanoString panel, not bulk RNA-seq.",
        "next_action": "Implement a targeted-panel scorer only if enough locked panel genes overlap; keep separate from strict bulk RNA-seq validation.",
    },
    {
        "lead_id": "GSE122220",
        "source_database": "GEO",
        "candidate_accession": "GSE122220",
        "title": "Melanoma tumor samples before checkpoint inhibitors",
        "cancer_type": "melanoma",
        "therapy": "checkpoint inhibitors",
        "sample_type": "tumor tissue",
        "platform": "Illumina HumanHT-12 V4.0 expression beadchip",
        "timepoint": "before treatment",
        "response_label": "GEO sample characteristics list CR/PR/SD/PD best-response categories",
        "independence_role": "small_array_external_candidate",
        "eligibility_status": "low_n_array_public_processed",
        "evidence_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122220",
        "evidence_note": "GEO lists 10 pretreatment melanoma tumor biopsy array samples; the series matrix sample characteristics include treatment and CR/PR/SD/PD response fields.",
        "next_action": "Run the registered array preprocessing and score only as low-n platform sensitivity, not as strict bulk RNA-seq validation.",
    },
    {
        "lead_id": "GSE93157",
        "source_database": "GEO",
        "candidate_accession": "GSE93157",
        "title": "Mixed tumor anti-PD-1 NanoString/pan-cancer ICB panel",
        "cancer_type": "mixed; includes melanoma subset",
        "therapy": "anti-PD-1",
        "sample_type": "tumor tissue",
        "platform": "NanoString/panel expression",
        "timepoint": "pretreatment/mixed depending sample",
        "response_label": "response harmonized in current processed metadata",
        "independence_role": "panel_transfer_or_pan_cancer_support",
        "eligibility_status": "panel_transfer_processed",
        "evidence_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93157",
        "evidence_note": "Processed locally; better for panel transfer and pan-cancer portability than strict melanoma-only headline validation.",
        "next_action": "Use for panel compatibility and transfer analyses; separate from pure melanoma anti-PD-1 tumor RNA-seq claims.",
    },
    {
        "lead_id": "IMVIGOR210",
        "source_database": "IMvigor210 package/public processed",
        "candidate_accession": "IMvigor210",
        "title": "Urothelial atezolizumab trial transcriptomic cohort",
        "cancer_type": "urothelial carcinoma",
        "therapy": "anti-PD-L1 atezolizumab",
        "sample_type": "tumor tissue",
        "platform": "bulk RNA-seq",
        "timepoint": "trial baseline",
        "response_label": "RECIST response available in package resources",
        "independence_role": "pan_cancer_external_only",
        "eligibility_status": "ineligible_for_melanoma_primary",
        "evidence_url": "http://research-pub.gene.com/IMvigor210CoreBiologies/",
        "evidence_note": "Biologically valuable pan-cancer/urothelial external set but not melanoma anti-PD-1; should not repair melanoma-specific external AUROC.",
        "next_action": "Keep for pan-cancer transfer if processed; do not use as melanoma primary validation.",
    },
    {
        "lead_id": "ICBATLAS_TIGER_AGGREGATES",
        "source_database": "ICBatlas/TIGER",
        "candidate_accession": "ICBatlas; TIGER immunotherapy tables",
        "title": "Aggregated immunotherapy transcriptomic resources",
        "cancer_type": "multi-cancer including melanoma",
        "therapy": "PD-1/PD-L1/CTLA-4 inhibitors",
        "sample_type": "tumor tissue or processed expression depending dataset",
        "platform": "RNA-seq and microarray aggregates",
        "timepoint": "mixed",
        "response_label": "harmonized processed labels",
        "independence_role": "discovery_resource_and_duplicate_locator",
        "eligibility_status": "aggregate_duplicate_screen",
        "evidence_url": "http://tiger.canceromics.org/tiger/downloadtable.php?table=immunotherapy",
        "evidence_note": "Useful to identify processed expression and response labels, but aggregate resources can duplicate GEO/ENA/cBioPortal source cohorts.",
        "next_action": "Use as a locator and provenance cross-check; trace each cohort back to original accession before claim wording.",
    },
]


ALIASES: dict[str, list[tuple[str, str]]] = {
    "LIU_MGSP_PHS000452": [
        ("data/processed/bulk", "PHS000452_LIU_LIKE_PRE"),
        ("data/processed/cbioportal_melanoma", "CBIO_LIU_DFCI_2019_PRE"),
    ],
    "GIDE_PRJEB23709": [
        ("data/processed/bulk", "PRJEB23709_PD1_PRE"),
        ("data/processed/bulk", "PRJEB23709_COMBO_PRE"),
    ],
    "RIAZ_GSE91061": [("data/processed/bulk", "GSE91061")],
    "HUGO_GSE78220": [("data/processed/bulk", "GSE78220")],
    "MGH_GSE115821": [("data/processed/bulk", "GSE115821")],
    "MGH_GSE168204": [("data/processed/bulk", "GSE168204")],
    "GSE145996": [("data/processed/bulk", "GSE145996")],
    "GSE123728": [("data/processed/bulk", "GSE123728")],
    "GSE165745": [("data/processed/bulk", "GSE165745")],
    "GSE122220": [("data/processed/bulk", "GSE122220")],
    "GSE93157": [("data/processed/bulk", "GSE93157")],
    "IMVIGOR210": [("data/processed/bulk", "IMvigor210")],
}


ACCESSION_KEYS: dict[str, list[str]] = {
    "LIU_MGSP_PHS000452": ["Liu_DFCI_melanoma", "CBIO_LIU_DFCI_2019_PRE"],
    "GIDE_PRJEB23709": ["PRJEB23709", "CBIO_IATLAS_GIDE_2019_PRE"],
    "RIAZ_GSE91061": ["GSE91061", "CBIO_IATLAS_RIAZ_2017_PRE"],
    "HUGO_GSE78220": ["GSE78220", "CBIO_IATLAS_HUGO_2016_PRE"],
    "MGH_GSE115821": ["GSE115821"],
    "MGH_GSE168204": ["GSE168204"],
    "GSE145996": ["GSE145996"],
    "LEE_RIZOS_EGAS00001001552": ["EGAS00001001552"],
    "ABRIL_RODRIGUEZ_PHS001919": ["phs001919"],
    "MGH_HACOHEN_PHS002683": ["phs002683"],
    "GSE123728": ["GSE123728"],
    "GSE165745": ["GSE165745"],
    "GSE122220": ["GSE122220"],
    "GSE93157": ["GSE93157"],
    "IMVIGOR210": ["IMvigor210"],
}


def _load_registry(registry_path: Path) -> dict[str, dict[str, object]]:
    if not registry_path.exists():
        return {}
    registry = yaml.safe_load(registry_path.read_text(encoding="utf-8")) or {}
    return {str(row.get("accession")): row for row in registry.get("cohorts", [])}


def _processed_status(root: Path, lead_id: str) -> tuple[str, str]:
    aliases = ALIASES.get(lead_id, [])
    if not aliases:
        return "not_applicable", ""
    present: list[str] = []
    missing: list[str] = []
    for directory, cohort in aliases:
        expr = root / directory / f"{cohort}.expr.tsv"
        meta = root / directory / f"{cohort}.metadata.tsv"
        if expr.exists() and meta.exists():
            present.append(cohort)
        else:
            missing.append(cohort)
    if present and not missing:
        return "processed", ",".join(present)
    if present:
        return "partially_processed", f"present={','.join(present)}; missing={','.join(missing)}"
    return "not_processed", ",".join(missing)


def build_lead_triage(root: Path = ROOT, registry_path: Path | None = None) -> pd.DataFrame:
    registry_path = registry_path if registry_path is not None else root / "config" / "data_registry.yml"
    registry = _load_registry(registry_path)
    rows: list[dict[str, object]] = []
    for lead in LEADS:
        row = dict(lead)
        status, cohorts = _processed_status(root, lead["lead_id"])
        keys = ACCESSION_KEYS.get(lead["lead_id"], [lead["candidate_accession"]])
        registry_hits = [key for key in keys if key in registry]
        registry_missing = [key for key in keys if key not in registry]
        row["current_registry_status"] = (
            f"registered={','.join(registry_hits)}; missing={','.join(registry_missing)}"
            if registry_hits or registry_missing
            else "not_registered"
        )
        row["processed_status"] = status
        row["processed_cohorts"] = cohorts
        row["strict_melanoma_primary_suitability"] = _strict_suitability(row)
        rows.append(row)
    return pd.DataFrame(rows)


def _strict_suitability(row: dict[str, object]) -> str:
    status = str(row["eligibility_status"])
    if status.startswith("eligible_processed") and str(row["platform"]) == "bulk RNA-seq":
        return "usable_with_predeclared_train_external_boundary"
    if status == "controlled_access_required":
        return "potentially_high_value_after_access"
    if "panel" in status:
        return "panel_transfer_not_bulk_primary"
    if "array" in status:
        return "low_n_platform_sensitivity_only"
    if status.startswith("ineligible"):
        return "not_melanoma_primary_validation"
    return "curation_or_duplicate_screen_only"


def write_markdown(triage: pd.DataFrame, out_md: Path) -> None:
    counts = triage["eligibility_status"].value_counts().to_dict()
    lines = [
        "# Public External Melanoma ICB Lead Triage",
        "",
        "This file records public, aggregate, and controlled-access leads for the melanoma anti-PD-1 tumor-tissue validation target. It is a no-fabrication control: a row is not validation evidence until expression, response labels, sample timing, and independence are verified by a registered pipeline.",
        "",
        "## Status Counts",
        "",
    ]
    for status, count in sorted(counts.items()):
        lines.append(f"- {status}: {count}")
    lines.extend(["", "## Leads", ""])
    for _, row in triage.iterrows():
        lines.append(
            "- **{lead_id}**: {title}; status={status}; processed={processed}; suitability={suitability}; next={next_action}".format(
                lead_id=row["lead_id"],
                title=row["title"],
                status=row["eligibility_status"],
                processed=row["processed_status"],
                suitability=row["strict_melanoma_primary_suitability"],
                next_action=row["next_action"],
            )
        )
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "NanoString, microarray, pan-cancer, aggregate-resource, and duplicate-source leads can support panel transfer, provenance checks, or sensitivity analyses, but they do not by themselves satisfy the strict independent melanoma bulk RNA-seq external AUROC target.",
            "",
        ]
    )
    out_md.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", default="config/data_registry.yml")
    parser.add_argument("--out", default="deliverables/public_external_lead_triage_20260527.tsv")
    parser.add_argument("--out-md", default="deliverables/public_external_lead_triage_20260527.md")
    args = parser.parse_args()

    triage = build_lead_triage(ROOT, ROOT / args.registry)
    out = ROOT / args.out
    out_md = ROOT / args.out_md
    out.parent.mkdir(parents=True, exist_ok=True)
    triage.to_csv(out, sep="\t", index=False)
    write_markdown(triage, out_md)
    print(
        json.dumps(
            {
                "n_leads": int(len(triage)),
                "processed": int(triage["processed_status"].eq("processed").sum()),
                "controlled": int(triage["eligibility_status"].eq("controlled_access_required").sum()),
                "panel_transfer": int(triage["eligibility_status"].str.contains("panel").sum()),
            },
            ensure_ascii=False,
        )
    )
    print(f"Wrote {out}")
    print(f"Wrote {out_md}")


if __name__ == "__main__":
    main()
