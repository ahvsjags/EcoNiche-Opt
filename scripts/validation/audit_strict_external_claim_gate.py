from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TARGET_MODEL = "EcoNiche-Opt-HeuristicEcology-LockedPanel"
TRANSFER_HEAD = "EcoNiche-Opt-PD1LikeTransferHead"
STRICT_FAMILY = "strict_pd1_like_external"
STRICT_COHORTS = ["GSE145996", "PHS000452_LIU_LIKE_PRE"]
PRIMARY_ENDPOINT = "strict_recist"


def _fmt(value: object, digits: int = 3) -> str:
    numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(numeric):
        return "NA"
    return f"{float(numeric):.{digits}f}"


def _claim_status(
    target_auroc: float,
    family_delta: float | None,
    fdr_q: float | None,
    min_primary_auroc: float,
    alpha: float,
) -> str:
    if np.isnan(target_auroc):
        return "RESULT_PENDING"
    if (
        target_auroc >= min_primary_auroc
        and family_delta is not None
        and fdr_q is not None
        and family_delta > 0
        and fdr_q <= alpha
    ):
        return "primary_external_claim_supported"
    if family_delta is not None and family_delta > 0:
        return "modest_point_estimate_only"
    return "not_supported_for_external_superiority"


def build_strict_external_claim_gate(
    metrics: pd.DataFrame,
    family: pd.DataFrame,
    rescue: pd.DataFrame,
    min_primary_auroc: float = 0.70,
    alpha: float = 0.05,
) -> tuple[pd.DataFrame, str]:
    rows: list[dict[str, object]] = []

    strict_family = family[family["validation_family"].astype(str).eq(STRICT_FAMILY)].copy()
    for _, row in strict_family.iterrows():
        endpoint = str(row["endpoint"])
        target_auroc = float(row["target_AUROC"])
        mean_delta = float(row["mean_delta_vs_signature_family"])
        fdr_q = float(row["two_sided_fdr_q"])
        status = _claim_status(target_auroc, mean_delta, fdr_q, min_primary_auroc, alpha)
        rows.append(
            {
                "gate_id": f"strict_family_{endpoint}",
                "evidence_type": "strict_pd1_like_family_omnibus",
                "endpoint": endpoint,
                "cohort_set": "+".join(STRICT_COHORTS),
                "n_samples": int(row["n_samples"]),
                "target_AUROC": target_auroc,
                "family_mean_AUROC": float(row["mean_signature_AUROC"]),
                "best_signature_AUROC": float(row["best_signature_AUROC"]),
                "delta_vs_family_mean": mean_delta,
                "two_sided_fdr_q": fdr_q,
                "claim_status": status,
                "allowed_claim": (
                    "strict melanoma PD1-like external scoring may be reported as a locked, refitting-free, "
                    "modest point-estimate external support result"
                    if status == "modest_point_estimate_only"
                    else "strict melanoma PD1-like external superiority claim is supported"
                    if status == "primary_external_claim_supported"
                    else "report cohort-specific metrics only"
                ),
                "blocked_claim": (
                    "Do not claim AUROC >=0.70, FDR-supported strict melanoma external validation, "
                    "or clinical-level validation from this stratum."
                    if status != "primary_external_claim_supported"
                    else ""
                ),
                "required_next_evidence": (
                    "Freeze the score and validate it in an independent pretreatment melanoma tumor-tissue cohort "
                    "with anti-PD1/anti-PD1-based therapy and RECIST labels, ideally n>=50-100."
                    if status != "primary_external_claim_supported"
                    else "Maintain locked scoring and cite the FDR-supported strict external evidence."
                ),
            }
        )

    if {"model_name", "cohort"}.issubset(metrics.columns):
        target_metrics = metrics[
            metrics["model_name"].astype(str).eq(TARGET_MODEL) & metrics["cohort"].astype(str).isin(STRICT_COHORTS)
        ].copy()
    else:
        target_metrics = pd.DataFrame()
    for _, row in target_metrics.iterrows():
        auroc = float(row["AUROC"])
        endpoint = str(row["endpoint"])
        rows.append(
            {
                "gate_id": f"strict_cohort_{endpoint}_{row['cohort']}",
                "evidence_type": "strict_pd1_like_single_cohort_metric",
                "endpoint": endpoint,
                "cohort_set": str(row["cohort"]),
                "n_samples": int(row["n_samples"]),
                "target_AUROC": auroc,
                "family_mean_AUROC": np.nan,
                "best_signature_AUROC": np.nan,
                "delta_vs_family_mean": np.nan,
                "two_sided_fdr_q": np.nan,
                "claim_status": "cohort_support_moderate" if auroc >= 0.60 else "cohort_support_weak_or_modest",
                "allowed_claim": "report as a cohort-specific locked external metric with n, endpoint, and calibration metrics",
                "blocked_claim": "Do not generalize this single cohort into clinical validation or pan-cohort superiority.",
                "required_next_evidence": "Use only as one component of the strict external evidence hierarchy.",
            }
        )

    if {"model_name", "cohort", "endpoint"}.issubset(rescue.columns):
        rescue_pooled = rescue[
            rescue["model_name"].astype(str).eq(TRANSFER_HEAD)
            & rescue["cohort"].astype(str).eq("+".join(STRICT_COHORTS))
            & rescue["endpoint"].astype(str).eq(PRIMARY_ENDPOINT)
        ].copy()
    else:
        rescue_pooled = pd.DataFrame()
    for _, row in rescue_pooled.iterrows():
        rows.append(
            {
                "gate_id": "pd1_like_transfer_head_pooled_strict_recist",
                "evidence_type": "secondary_pd1_like_rescue",
                "endpoint": str(row["endpoint"]),
                "cohort_set": str(row["cohort"]),
                "n_samples": int(row["n_samples"]),
                "target_AUROC": float(row["AUROC"]),
                "family_mean_AUROC": np.nan,
                "best_signature_AUROC": np.nan,
                "delta_vs_family_mean": np.nan,
                "two_sided_fdr_q": np.nan,
                "claim_status": "secondary_model_development_rescue_not_primary_claim",
                "allowed_claim": "describe as a transparent rescue analysis generated after identifying weak strict external performance",
                "blocked_claim": "Do not use the transfer head as the primary locked external validation claim without a fresh cohort.",
                "required_next_evidence": "Freeze the transfer head first, then test it on a newly independent melanoma tumor-tissue cohort.",
            }
        )

    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["evidence_type", "endpoint", "cohort_set"]).reset_index(drop=True)

    primary = out[out["gate_id"].eq(f"strict_family_{PRIMARY_ENDPOINT}")]
    if primary.empty:
        headline = "Strict melanoma PD1-like external gate: RESULT_PENDING because the strict family row is missing."
    else:
        row = primary.iloc[0]
        headline = (
            "Strict melanoma PD1-like external gate: "
            f"{row['claim_status']} for {PRIMARY_ENDPOINT}; "
            f"n={int(row['n_samples'])}, AUROC={_fmt(row['target_AUROC'])}, "
            f"family mean AUROC={_fmt(row['family_mean_AUROC'])}, "
            f"delta={_fmt(row['delta_vs_family_mean'])}, q={_fmt(row['two_sided_fdr_q'])}. "
            "This supports locked refitting-free external scoring, not a high-strength clinical validation claim."
        )
    return out, headline


def write_markdown(report: pd.DataFrame, headline: str, out_md: Path) -> None:
    lines = [
        "# Strict Melanoma PD1-like External Claim Gate",
        "",
        headline,
        "",
        "## What this gate checks",
        "",
        "- Strict external evidence is restricted to GSE145996 and PHS000452_LIU_LIKE_PRE.",
        "- The primary endpoint for this gate is strict RECIST.",
        "- The gate reads locked external outputs only; it does not refit, recalibrate, select thresholds, or change labels.",
        "- A primary external superiority claim requires AUROC >= 0.70 and FDR-supported improvement over the predeclared signature family.",
        "",
        "## Evidence summary",
        "",
    ]
    if report.empty:
        lines.append("No gate rows were produced.")
    else:
        for _, row in report.iterrows():
            lines.append(
                f"- {row['gate_id']}: {row['claim_status']}; endpoint={row['endpoint']}; "
                f"cohort_set={row['cohort_set']}; n={int(row['n_samples'])}; "
                f"AUROC={_fmt(row['target_AUROC'])}; q={_fmt(row['two_sided_fdr_q'])}."
            )
    lines.extend(
        [
            "",
            "## Allowed wording",
            "",
            "The manuscript may state that strict melanoma PD1-like external cohorts were scored with a locked, refitting-free rule and showed modest point-estimate support. It should include cohort set, endpoint, n, AUROC, family mean AUROC, FDR q value, and calibration metrics.",
            "",
            "## Blocked wording",
            "",
            "Do not state that the strict melanoma external validation is clinically strong, AUROC >=0.70, FDR-supported, or prospectively validated unless a newly independent tumor-tissue cohort proves those claims after the score is frozen.",
            "",
            "## Required next evidence",
            "",
            "The next decisive dataset is an independent pretreatment melanoma tumor-tissue cohort treated with anti-PD1/anti-PD1-based therapy, RECIST CR/PR/SD/PD labels, and RNA-seq/NanoString/qPCR measurement of the locked panel, ideally n>=50-100.",
            "",
        ]
    )
    out_md.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--metrics",
        default="results/locked_external_panel_validation_calibrated_20260519/locked_external_metrics.tsv",
    )
    parser.add_argument(
        "--family",
        default="results/locked_external_panel_validation_calibrated_20260519/locked_external_signature_family_omnibus.tsv",
    )
    parser.add_argument("--rescue", default="results/pd1_like_external_rescue/pd1_like_rescue_metrics.tsv")
    parser.add_argument("--out", default="deliverables/strict_melanoma_external_claim_gate_20260527.tsv")
    parser.add_argument("--out-md", default="deliverables/strict_melanoma_external_claim_gate_20260527.md")
    parser.add_argument("--min-primary-auroc", type=float, default=0.70)
    parser.add_argument("--alpha", type=float, default=0.05)
    args = parser.parse_args()

    metrics = pd.read_csv(ROOT / args.metrics, sep="\t")
    family = pd.read_csv(ROOT / args.family, sep="\t")
    rescue = pd.read_csv(ROOT / args.rescue, sep="\t")
    report, headline = build_strict_external_claim_gate(
        metrics,
        family,
        rescue,
        min_primary_auroc=args.min_primary_auroc,
        alpha=args.alpha,
    )
    out = ROOT / args.out
    out_md = ROOT / args.out_md
    out.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(out, sep="\t", index=False)
    write_markdown(report, headline, out_md)
    print(headline)
    print(f"Wrote {out}")
    print(f"Wrote {out_md}")


if __name__ == "__main__":
    main()
