from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche_opt.data.registry import load_registry
from econiche.registry import normalize_access_status


PRIORITY_ACCESSIONS = [
    "Liu_DFCI_melanoma",
    "CBIO_LIU_DFCI_2019_PRE",
    "CBIO_IATLAS_LIU_2019_PRE",
    "CBIO_IATLAS_GIDE_2019_PRE",
    "CBIO_IATLAS_RIAZ_2017_PRE",
    "CBIO_IATLAS_HUGO_2016_PRE",
    "EGAS00001001552",
    "PRJEB23709",
    "GSE91061",
    "GSE78220",
    "GSE115821",
    "GSE168204",
    "GSE145996",
    "GSE93157",
    "IMvigor210",
]


def _processed_status(accession: str) -> tuple[str, str]:
    cohort_aliases = {
        "Liu_DFCI_melanoma": ["PHS000452_LIU_LIKE_PRE"],
        "PRJEB23709": ["PRJEB23709_PD1_PRE", "PRJEB23709_COMBO_PRE"],
        "EGAS00001001552": ["EGAS00001001552"],
    }
    cohorts = cohort_aliases.get(accession, [accession])
    present = []
    missing = []
    for cohort in cohorts:
        processed_root = ROOT / "data" / "processed" / ("cbioportal_melanoma" if cohort.startswith("CBIO_") else "bulk")
        expr = processed_root / f"{cohort}.expr.tsv"
        meta = processed_root / f"{cohort}.metadata.tsv"
        if expr.exists() and meta.exists():
            present.append(cohort)
        else:
            missing.append(cohort)
    if present and not missing:
        return "processed", ",".join(present)
    if present:
        return "partially_processed", f"present={','.join(present)}; missing={','.join(missing)}"
    return "not_processed", ",".join(missing)


def build_candidate_audit(registry_path: Path) -> pd.DataFrame:
    registry = load_registry(registry_path)
    by_accession = {str(row.get("accession")): row for row in registry.get("cohorts", [])}
    rows: list[dict[str, object]] = []
    for accession in PRIORITY_ACCESSIONS:
        cohort = by_accession.get(accession)
        if cohort is None:
            rows.append(
                {
                    "accession": accession,
                    "name": "",
                    "priority": "",
                    "access": "MISSING_FROM_REGISTRY",
                    "normalized_access_status": "missing",
                    "processed_status": "not_processed",
                    "processed_cohorts": "",
                    "primary_use": "candidate_or_requested_priority",
                    "evidence_basis": "",
                    "next_action": "Add a no-fabrication registry entry before use.",
                }
            )
            continue
        access = str(cohort.get("access", ""))
        normalized = normalize_access_status(access)
        processed_status, processed_cohorts = _processed_status(accession)
        if normalized == "controlled":
            next_action = "Request controlled access; do not create substitute expression data."
        elif processed_status == "processed":
            next_action = "Keep as available evidence and score only under predeclared train/test boundaries."
        else:
            next_action = "Download and curate sample-level response metadata before modeling."
        rows.append(
            {
                "accession": accession,
                "name": cohort.get("name", accession),
                "priority": cohort.get("priority", ""),
                "access": access,
                "normalized_access_status": normalized,
                "processed_status": processed_status,
                "processed_cohorts": processed_cohorts,
                "primary_use": ",".join(cohort.get("uses", [])) if isinstance(cohort.get("uses"), list) else cohort.get("uses", ""),
                "evidence_basis": cohort.get("notes", ""),
                "next_action": next_action,
            }
        )
    return pd.DataFrame(rows)


def write_markdown(audit: pd.DataFrame, out_md: Path) -> None:
    lines = [
        "# Melanoma External Data Candidate Audit",
        "",
        "This audit tracks the public and controlled melanoma/ICB datasets relevant to the top-tier strict external validation target. Controlled datasets are evidence opportunities, not current results.",
        "",
    ]
    for _, row in audit.iterrows():
        lines.append(
            f"- **{row['accession']}** ({row['name']}): access={row['normalized_access_status']}; "
            f"processed={row['processed_status']} ({row['processed_cohorts']}); next={row['next_action']}"
        )
    lines.extend(
        [
            "",
            "## No-fabrication boundary",
            "",
            "ACCESS_RESTRICTED EGA/dbGaP datasets can strengthen the external validation only after approved access and registered preprocessing. Until then, they remain candidate evidence and must not be represented as completed validation.",
            "",
        ]
    )
    out_md.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", default="config/data_registry.yml")
    parser.add_argument("--out", default="deliverables/melanoma_external_data_candidates_20260527.tsv")
    parser.add_argument("--out-md", default="deliverables/melanoma_external_data_candidates_20260527.md")
    args = parser.parse_args()

    audit = build_candidate_audit(ROOT / args.registry)
    out = ROOT / args.out
    out_md = ROOT / args.out_md
    out.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(out, sep="\t", index=False)
    write_markdown(audit, out_md)
    controlled = int(audit["normalized_access_status"].eq("controlled").sum())
    processed = int(audit["processed_status"].eq("processed").sum())
    print(f"Melanoma external candidate audit: processed={processed}, controlled_or_restricted={controlled}")
    print(f"Wrote {out}")
    print(f"Wrote {out_md}")


if __name__ == "__main__":
    main()
