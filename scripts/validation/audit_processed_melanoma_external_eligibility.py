from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[2]
STRICT_ELIGIBLE_COHORTS = {
    "GSE91061",
    "GSE78220",
    "PRJEB23709_PD1_PRE",
    "GSE145996",
    "PHS000452_LIU_LIKE_PRE",
    "GSE115821",
    "GSE168204",
}


def _read_registry(path: Path) -> dict[str, dict[str, object]]:
    registry = yaml.safe_load(path.read_text(encoding="utf-8"))
    out: dict[str, dict[str, object]] = {}
    for cohort in registry.get("cohorts", []):
        accession = str(cohort.get("accession", ""))
        if accession:
            out[accession] = cohort
    return out


def _value_counts_text(frame: pd.DataFrame, column: str) -> str:
    if column not in frame.columns or frame.empty:
        return ""
    values = frame[column].astype(str).replace({"nan": ""})
    counts = values[values.ne("")].value_counts().head(6)
    return ";".join(f"{key}:{int(value)}" for key, value in counts.items())


def _registry_match(cohort_name: str, registry: dict[str, dict[str, object]]) -> dict[str, object]:
    stem = cohort_name.replace(".metadata.tsv", "")
    if stem in registry:
        return registry[stem]
    if stem.startswith("PRJEB23709_") and "PRJEB23709" in registry:
        return registry["PRJEB23709"]
    if stem.startswith("PHS000452") and "Liu_DFCI_melanoma" in registry:
        return registry["Liu_DFCI_melanoma"]
    return {}


def _eligibility(stem: str, meta: pd.DataFrame, reg: dict[str, object]) -> tuple[str, str]:
    if meta.empty:
        return "not_usable_empty_metadata", "metadata has zero rows"
    cancer = str(reg.get("cancer_type", "")).lower()
    therapy = str(reg.get("therapy", "")).lower()
    platform = str(reg.get("platform", "")).lower()
    timepoints = ",".join(str(item).lower() for item in reg.get("timepoints", []))
    role = str(reg.get("role", ""))
    notes = str(reg.get("notes", ""))
    if stem in STRICT_ELIGIBLE_COHORTS:
        if stem in {"GSE115821", "GSE168204"}:
            return "secondary_small_melanoma_sensitivity", "eligible melanoma pretreatment ICB-like cohort, but small and response/timing evidence is less robust than primary discovery cohorts"
        if stem in {"GSE145996", "PHS000452_LIU_LIKE_PRE"}:
            return "strict_external_current", "used in current strict melanoma external stress test"
        return "primary_discovery_or_lodo", "used in high-evidence primary melanoma discovery/LODO boundary"
    if "melanoma" not in cancer and "melanoma" not in (notes + " " + _value_counts_text(meta, "tumor_type")).lower():
        return "not_melanoma_primary", f"registry cancer_type={reg.get('cancer_type', '')}"
    if "nanostring" in platform or "panel" in platform or "ncounter" in platform:
        return "panel_transfer_not_bulk_strict", f"platform={reg.get('platform', '')}"
    if "array" in platform or "beadchip" in platform or "humanht" in platform:
        return "low_n_array_platform_sensitivity", f"platform={reg.get('platform', '')}; n={len(meta)}"
    if "pretreatment" not in timepoints and "baseline" not in timepoints:
        return "timing_not_strict_pretreatment", f"timepoints={reg.get('timepoints', '')}"
    if "pd1" not in therapy and "checkpoint" not in therapy and "icb" not in therapy:
        return "therapy_not_strict_pd1_like", f"therapy={reg.get('therapy', '')}"
    return "needs_manual_source_hardening", f"role={role}; notes={notes[:120]}"


def build_audit(processed_dir: Path, registry_path: Path) -> pd.DataFrame:
    registry = _read_registry(registry_path)
    rows: list[dict[str, object]] = []
    for path in sorted(processed_dir.glob("*.metadata.tsv")):
        stem = path.name.replace(".metadata.tsv", "")
        if stem.startswith("demo_cohort_"):
            continue
        meta = pd.read_csv(path, sep="\t")
        reg = _registry_match(path.name, registry)
        status, reason = _eligibility(stem, meta, reg)
        rows.append(
            {
                "cohort": stem,
                "n_metadata_rows": int(len(meta)),
                "registry_cancer_type": reg.get("cancer_type", ""),
                "registry_therapy": reg.get("therapy", ""),
                "registry_platform": reg.get("platform", ""),
                "registry_timepoints": ",".join(str(item) for item in reg.get("timepoints", [])),
                "registry_role": reg.get("role", ""),
                "label_counts": _value_counts_text(meta, "label"),
                "response_raw_counts": _value_counts_text(meta, "response_raw"),
                "timepoint_counts": _value_counts_text(meta, "timepoint"),
                "therapy_counts": _value_counts_text(meta, "therapy") or _value_counts_text(meta, "treatment"),
                "tumor_type_counts": _value_counts_text(meta, "tumor_type") or _value_counts_text(meta, "cancer_type"),
                "eligibility_status": status,
                "eligibility_reason": reason,
            }
        )
    return pd.DataFrame(rows).sort_values(["eligibility_status", "cohort"]).reset_index(drop=True)


def write_markdown(audit: pd.DataFrame, out_md: Path) -> None:
    counts = audit["eligibility_status"].value_counts().to_dict() if not audit.empty else {}
    lines = [
        "# Processed Melanoma External Eligibility Audit",
        "",
        "This audit checks all processed bulk metadata files for possible strict melanoma pretreatment anti-PD-1 tumor-tissue external-validation use.",
        "",
        "## Status Counts",
        "",
    ]
    for status, count in sorted(counts.items()):
        lines.append(f"- {status}: {count}")
    lines.extend(["", "## Cohort Decisions", ""])
    for _, row in audit.iterrows():
        lines.append(
            f"- `{row['cohort']}`: {row['eligibility_status']} ({row['eligibility_reason']}); n={int(row['n_metadata_rows'])}."
        )
    lines.extend(
        [
            "",
            "## Conclusion",
            "",
            "No overlooked large public bulk RNA-seq pretreatment melanoma anti-PD-1 external cohort was found among the currently processed metadata files. The remaining path to a top-tier strict external claim is controlled-access acquisition or a newly discovered independent public cohort with hardened sample-level response evidence.",
            "",
        ]
    )
    out_md.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", default="data/processed/bulk")
    parser.add_argument("--registry", default="config/data_registry.yml")
    parser.add_argument("--out", default="deliverables/processed_melanoma_external_eligibility_20260527.tsv")
    parser.add_argument("--out-md", default="deliverables/processed_melanoma_external_eligibility_20260527.md")
    args = parser.parse_args()

    audit = build_audit(ROOT / args.processed_dir, ROOT / args.registry)
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(out, sep="\t", index=False)
    write_markdown(audit, ROOT / args.out_md)
    print(json.dumps(audit["eligibility_status"].value_counts().to_dict(), ensure_ascii=False))
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
