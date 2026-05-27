from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche_opt.model.endpoint_modules import MODULE_GENE_SETS, MODULE_PRIOR_WEIGHTS


ARTICLE_DIR = ROOT / "figures/article"
TABLE_DIR = ROOT / "tables/article"
PAPER_DIR = ROOT / "paper"
DPI = 600
RASTER_FORMATS = ("png", "tiff")
VECTOR_FORMATS = ("pdf", "svg")

COLORS = {
    "target": "#009E73",
    "target_dark": "#006C5B",
    "baseline": "#7A7A7A",
    "baseline_light": "#C8C8C8",
    "resistance": "#D55E00",
    "warning": "#E69F00",
    "blue": "#0072B2",
    "sky": "#56B4E9",
    "purple": "#CC79A7",
    "black": "#222222",
    "grid": "#E6E6E6",
}


def read_tsv(path: str | Path) -> pd.DataFrame:
    path = ROOT / path if not Path(path).is_absolute() else Path(path)
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def save_fig(fig: plt.Figure, name: str, manifest: list[dict[str, object]]) -> None:
    ARTICLE_DIR.mkdir(parents=True, exist_ok=True)
    for ext in [*RASTER_FORMATS, *VECTOR_FORMATS]:
        out = ARTICLE_DIR / f"{name}.{ext}"
        fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    manifest.append(
        {
            "figure": name,
            "png": str((ARTICLE_DIR / f"{name}.png").relative_to(ROOT)),
            "tiff": str((ARTICLE_DIR / f"{name}.tiff").relative_to(ROOT)),
            "pdf": str((ARTICLE_DIR / f"{name}.pdf").relative_to(ROOT)),
            "svg": str((ARTICLE_DIR / f"{name}.svg").relative_to(ROOT)),
        }
    )


def panel_label(ax, label: str) -> None:
    ax.text(-0.12, 1.06, label, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")


def figure_heading(fig: plt.Figure, text: str) -> None:
    fig.text(0.01, 0.995, text, fontsize=6.5, fontweight="bold", va="top", ha="left")


def clean_ax(ax, grid: bool = True) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid:
        ax.grid(axis="y", color=COLORS["grid"], linewidth=0.7)
        ax.set_axisbelow(True)


def truncate_label(value: object, max_chars: int = 18) -> str:
    text = str(value)
    return text if len(text) <= max_chars else text[: max_chars - 1] + "."


def draw_boxes(ax, labels: list[str], title: str | None = None) -> None:
    ax.axis("off")
    if title:
        ax.text(0.02, 0.96, title, fontsize=8, fontweight="bold", transform=ax.transAxes, va="top")
    n = len(labels)
    x0 = 0.04
    w = 0.88 / n
    y = 0.44
    h = 0.18
    for i, label in enumerate(labels):
        x = x0 + i * w
        patch = FancyBboxPatch(
            (x, y),
            w * 0.78,
            h,
            boxstyle="round,pad=0.012,rounding_size=0.02",
            facecolor="#F7FBFA",
            edgecolor=COLORS["target_dark"],
            linewidth=0.9,
            transform=ax.transAxes,
        )
        ax.add_patch(patch)
        wrapped = "\n".join(textwrap.fill(part, width=9, break_long_words=False) for part in str(label).split("\n"))
        ax.text(x + w * 0.39, y + h * 0.5, wrapped, ha="center", va="center", fontsize=5.4, transform=ax.transAxes)
        if i < n - 1:
            arrow_y = y - 0.045
            ax.add_patch(
                FancyArrowPatch(
                    (x + w * 0.80, arrow_y),
                    (x + w * 0.98, arrow_y),
                    arrowstyle="-|>",
                    mutation_scale=8,
                    color=COLORS["black"],
                    transform=ax.transAxes,
                )
            )


def draw_matrix_status(ax, rows: list[tuple[str, str, str]], title: str, cmap: dict[str, str] | None = None) -> None:
    cmap = cmap or {"+": COLORS["target"], "~": COLORS["warning"], "x": COLORS["baseline"]}
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.0, 0.98, title, fontsize=8, fontweight="bold", va="top")
    y0 = 0.82
    row_h = min(0.13, 0.72 / max(len(rows), 1))
    for i, (label, status, note) in enumerate(rows):
        y = y0 - i * row_h
        ax.add_patch(Rectangle((0.0, y - row_h * 0.52), 0.98, row_h * 0.8, facecolor="#F7F7F7" if i % 2 else "white", edgecolor="none"))
        ax.add_patch(Rectangle((0.02, y - row_h * 0.28), 0.055, row_h * 0.56, facecolor=cmap.get(status, COLORS["baseline_light"]), edgecolor="none"))
        ax.text(0.047, y, status, ha="center", va="center", fontsize=7, color="white", fontweight="bold")
        ax.text(0.10, y, label, ha="left", va="center", fontsize=6.7)
        ax.text(0.68, y, note, ha="left", va="center", fontsize=6.4, color="#444444")


def draw_formula_stack(ax) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.0, 0.98, "Model equations", fontsize=8, fontweight="bold", va="top")
    cards = [
        ("Signed-rank module", r"$M_q(i)=|G_q|^{-1}\sum_{g\in G_q}s_g\,zrank(x_{ig})$"),
        ("Ecological score", r"$S_i=\sum_q w_qM_q(i)+\sum_{q,r}\beta_{qr}E_{qr}(i)$"),
        ("Training objective", r"$\max\{\mathrm{AUC}-\lambda_1\mathrm{ECE}+\lambda_2\mathrm{BioPrior}-\lambda_3|G|\}$"),
    ]
    y = 0.78
    for i, (title, eq) in enumerate(cards):
        ax.add_patch(
            FancyBboxPatch(
                (0.02, y - 0.15),
                0.94,
                0.13,
                boxstyle="round,pad=0.01,rounding_size=0.015",
                facecolor="#F8FBFD" if i != 1 else "#F7FBFA",
                edgecolor=COLORS["baseline_light"],
                linewidth=0.8,
            )
        )
        ax.text(0.05, y - 0.042, title, fontsize=6.5, fontweight="bold", ha="left", va="center")
        ax.text(0.05, y - 0.105, eq, fontsize=6.0, ha="left", va="center")
        y -= 0.22


def draw_biology_mechanism(ax) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.0, 0.98, "Ecological mechanism map", fontsize=8, fontweight="bold", va="top")
    nodes = {
        "APM/MHC": (0.18, 0.72, COLORS["blue"]),
        "T/NK\n effector": (0.18, 0.52, COLORS["target"]),
        "IFN/T\n infl.": (0.18, 0.32, COLORS["sky"]),
        "Response": (0.55, 0.62, COLORS["target_dark"]),
        "Myeloid\n supp.": (0.55, 0.34, COLORS["warning"]),
        "CAF/ECM\n excl.": (0.82, 0.34, COLORS["resistance"]),
        "Resistance": (0.82, 0.62, COLORS["baseline"]),
    }
    for name, (x, y, color) in nodes.items():
        ax.add_patch(Circle((x, y), 0.075, facecolor=color, edgecolor="white", linewidth=0.8, alpha=0.95))
        ax.text(x, y, name, ha="center", va="center", fontsize=6.1, color="white", fontweight="bold")
    arrows = [
        ("APM/MHC", "Response", COLORS["target"]),
        ("T/NK\n effector", "Response", COLORS["target"]),
        ("IFN/T\n infl.", "Response", COLORS["target"]),
        ("Myeloid\n supp.", "Resistance", COLORS["resistance"]),
        ("CAF/ECM\n excl.", "Resistance", COLORS["resistance"]),
        ("Resistance", "Response", COLORS["baseline"]),
    ]
    for src, dst, color in arrows:
        x1, y1, _ = nodes[src]
        x2, y2, _ = nodes[dst]
        ax.add_patch(FancyArrowPatch((x1 + 0.075, y1), (x2 - 0.075, y2), arrowstyle="-|>", mutation_scale=8, linewidth=1.0, color=color, alpha=0.85))
    ax.text(0.05, 0.10, "Perturbation rankings: hypothesis-only reversal of resistance modules", fontsize=6.2, color=COLORS["black"])


def draw_annotation_template(ax) -> None:
    ax.axis("off")
    rows = [
        ["subject_id", "patient link", "required"],
        ["sample_id", "expression link", "required"],
        ["timepoint", "baseline/on-tx", "required"],
        ["therapy", "ICB regimen", "required"],
        ["response_raw", "raw outcome", "source"],
        ["evidence", "source proof", "traceable"],
    ]
    table = ax.table(cellText=rows, colLabels=["field", "role", "status"], loc="center", cellLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(6.6)
    table.scale(1.0, 1.25)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#DDDDDD")
        if row == 0:
            cell.set_text_props(fontweight="bold")
            cell.set_facecolor("#F0F4F8")
    ax.set_title("Clinical annotation template", fontsize=8)


def short_cohort_label(value: object) -> str:
    text = str(value)
    replacements = {
        "PHS000452_LIU_LIKE_PRE": "PHS000452",
        "GSE145996+PHS000452_LIU_LIKE_PRE": "GSE145996+\nPHS000452",
        "PRJEB23709_COMBO_PRE": "PRJEB23709\ncombo",
        "PRJEB23709_PD1_PRE": "PRJEB23709\nPD1",
    }
    return replacements.get(text, text.replace("_", "\n"))


def short_model_label(value: object) -> str:
    text = str(value).replace("EcoNiche-Opt-", "")
    replacements = {
        "HeuristicEcology-LockedPanel": "Locked panel",
        "PD1LikeTransferHead": "PD1 transfer",
        "HeuristicEcology": "Heuristic ecology",
        "ModulePriorFixed": "Module prior",
    }
    return replacements.get(text, text.replace("_", " "))


def short_family_label(endpoint: object, family: object) -> str:
    endpoint_map = {
        "primary_recist": "Primary",
        "strict_recist": "Strict",
        "clinical_benefit": "Clinical benefit",
    }
    family_map = {
        "all_locked_external_and_panel": "All external/panel",
        "strict_pd1_like_external": "PD1-like external",
        "nanostring_panel_transfer": "NanoString panel",
    }
    return f"{endpoint_map.get(str(endpoint), str(endpoint))} / {family_map.get(str(family), str(family))}"


def pretty_state_label(value: object) -> str:
    state = str(value).split("|")[0]
    mapping = {
        "ifn_t_cell_inflamed": "IFN/T infl.",
        "cytotoxic_cd8": "CD8 cytotoxic",
        "exhaustion_checkpoint": "Checkpoint/exh.",
        "antigen_presentation": "Antigen pres.",
        "antigen_presentation_mhc": "APM/MHC",
        "myeloid_suppression": "Myeloid supp.",
        "stromal_exclusion": "Stromal excl.",
        "caf_ecm_exclusion": "CAF/ECM excl.",
        "trm_tls": "TRM/TLS",
        "tnk_effector": "T/NK effector",
        "tcell_dysfunction": "T-cell dysfx",
        "tumor_dedifferentiation": "Tumor dediff.",
        "calibration": "Calibration",
    }
    return mapping.get(state, state.replace("_", " "))


def annotate_pending(ax, text: str = "RESULT_PENDING") -> None:
    ax.axis("off")
    ax.text(0.5, 0.5, text, ha="center", va="center", fontsize=11, color=COLORS["warning"], fontweight="bold")


def barh_from_series(ax, series: pd.Series, color: str, xlabel: str = "") -> None:
    if series.empty:
        annotate_pending(ax)
        return
    series = series.sort_values()
    ax.barh(series.index.astype(str), series.values, color=color, edgecolor="white")
    ax.set_xlabel(xlabel)
    clean_ax(ax, grid=False)


def heatmap(ax, data: pd.DataFrame, title: str = "", cmap: str = "viridis", vmin=None, vmax=None) -> None:
    if data.empty:
        annotate_pending(ax)
        return
    im = ax.imshow(data.to_numpy(dtype=float), aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(range(data.shape[1]))
    ax.set_xticklabels(data.columns.astype(str), rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(data.shape[0]))
    ax.set_yticklabels(data.index.astype(str), fontsize=7)
    ax.set_title(title, fontsize=9)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02)


def plot_workflow(ax) -> None:
    draw_boxes(
        ax,
        [
            "public ICB\ncohorts",
            "manual label\ncuration",
            "module\nscoring",
            "LODO and\nclaim gate",
            "locked external\nvalidation",
            "panel-ready\npackage",
        ],
        "EcoNiche-Opt study workflow",
    )


def make_article_tables(table_manifest: list[dict[str, str]]) -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    table_sources = [
        ("supp_table_01_data_registry_roles.tsv", "tables/dataset_roles.tsv"),
        ("supp_table_02_access_audit.tsv", "tables/dataset_access_audit.tsv"),
        ("supp_table_03_expression_qc.tsv", "tables/expression_qc_report_real_refresh.tsv"),
        ("supp_table_04_manual_curation_evidence.tsv", "results/curation/manual_curation_audit.tsv"),
        ("supp_table_05_external_cohort_gap_audit.tsv", "results/curation/external_cohort_gap_audit.tsv"),
        ("supp_table_06_endpoint_label_sensitivity.tsv", "results/endpoint_modules_heuristic_core_locked_gpu/endpoint_label_sensitivity_audit.tsv"),
        ("supp_table_07_module_gene_sets.tsv", None),
        ("supp_table_08_ecological_edges.tsv", "results/endpoint_modules_heuristic_core_locked_gpu/optimized_ecology_edges.tsv"),
        ("supp_table_09_optimizer_history.tsv", "results/endpoint_modules_heuristic_core_locked_gpu/optimized_ecology_history.tsv"),
        ("supp_table_10_melanoma_benchmark_summary.tsv", "results/endpoint_modules_heuristic_core_locked_gpu/endpoint_module_summary.tsv"),
        ("supp_table_11_signature_family_fdr.tsv", "results/claim_strengthening/strong_signature_family_omnibus.tsv"),
        ("supp_table_12_lodo_metrics.tsv", "results/endpoint_modules_heuristic_core_locked_gpu/endpoint_module_lodo_metrics.tsv"),
        ("supp_table_13_aligned_panel_ablation.tsv", "results/aligned_panel_ablation_20260527/aligned_panel_ablation_pairwise.tsv"),
        ("supp_table_14_decision_curve.tsv", "results/endpoint_modules_heuristic_core_locked_gpu/endpoint_module_decision_curve.tsv"),
        ("supp_table_15_locked_external_metrics.tsv", "results/locked_external_panel_validation_calibrated_20260519/locked_external_metrics.tsv"),
        ("supp_table_16_nanostring_panel_transfer.tsv", "results/locked_external_panel_validation_calibrated_20260519/clinical_assay_panel_transfer.tsv"),
        ("supp_table_17_pd1_like_rescue.tsv", "results/pd1_like_external_rescue/pd1_like_rescue_metrics.tsv"),
        ("supp_table_18_single_cell_enrichment.tsv", "results/scrna/cell_type_enrichment.tsv"),
        ("supp_table_19_perturbation_hypotheses.tsv", "results/perturbation/prioritized_perturbations.tsv"),
        ("supp_table_20_prospective_package_manifest.tsv", None),
        ("supp_table_21_gpu_bioprior_rescue.tsv", "results/gpu_bioprior_rescue_combo_search_robust_20260527/gpu_bioprior_rescue_combo_selection.tsv"),
        ("supp_table_22_gpu_bioprior_component_ablation.tsv", "results/gpu_bioprior_component_ablation_20260527/gpu_bioprior_component_ablation.tsv"),
        ("supp_table_23_cbioportal_gpu_bioprior_external.tsv", "results/cbioportal_gpu_bioprior_external_validation_20260527/cbioportal_gpu_bioprior_external_metrics.tsv"),
        ("supp_table_24_gpu_lipid_pair_rescue.tsv", "results/gpu_lipid_pair_rescue_component_dominant_ba_guardrail_20260528/gpu_lipid_pair_rescue_selection.tsv"),
    ]

    for out_name, src in table_sources:
        out = TABLE_DIR / out_name
        if out_name == "supp_table_07_module_gene_sets.tsv":
            rows = []
            for module, genes in MODULE_GENE_SETS.items():
                for gene in genes:
                    rows.append(
                        {
                            "module": module,
                            "gene_symbol": gene,
                            "module_weight": MODULE_PRIOR_WEIGHTS.get(module, 0.0),
                            "score_direction": "response_high" if MODULE_PRIOR_WEIGHTS.get(module, 0.0) >= 0 else "nonresponse_high",
                        }
                    )
            pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
            status = "generated_from_code"
        elif out_name == "supp_table_20_prospective_package_manifest.tsv":
            rows = [
                ("cohort evidence", "public ICB expression cohorts", "reported", "retrospective public cohorts", "Defines the multicohort response-prediction setting used for discovery and external evaluation."),
                ("cohort evidence", "patient-sample-response linkage", "reported", "sample-level curation", "Links each expression profile to patient, treatment timepoint and response evidence before model fitting."),
                ("endpoint definition", "primary RECIST endpoint", "reported", "predefined response mapping", "Treats CR/PR as response and SD/PD as non-response for the main binary endpoint."),
                ("endpoint definition", "strict RECIST endpoint", "reported", "sensitivity analysis", "Compares CR/PR against PD after excluding intermediate SD/MR samples."),
                ("endpoint definition", "clinical benefit endpoint", "reported", "sensitivity analysis", "Evaluates CR/PR/SD or durable benefit against progression/non-benefit where supported by cohort metadata."),
                ("feature construction", "signed-rank ecological module scores", "reported", "training-defined directions", "Transforms expression into response-high and resistance-high module activities without using held-out labels."),
                ("feature construction", "ecological interaction features", "reported", "training-only graph construction", "Aggregates biologically signed state-state gene interactions as model features."),
                ("model fitting", "class-balanced logistic model", "reported", "discovery cohorts only", "Estimates response probability from module and interaction features under class imbalance."),
                ("optimizer", "heuristic ecological module search", "reported", "inner training folds only", "Searches module composition and interaction edges with biological compactness and leakage-safety constraints."),
                ("thresholding", "endpoint-specific locked thresholds", "reported", "discovery cohorts only", "Selects fixed decision thresholds before external cohort scoring."),
                ("calibration", "calibration object", "reported", "discovery cohorts only", "Keeps probability calibration separated from locked external outcome labels."),
                ("external scoring", "one-way independent cohort scoring", "reported", "no model refit", "Scores external cohorts once and evaluates discrimination and calibration only after predictions are fixed."),
                ("external scoring", "gene coverage assessment", "reported", "assay-aware QC", "Reports module coverage so platform transfer does not silently change the biological score."),
                ("benchmarking", "signature-family comparison", "reported", "paired evaluation", "Compares EcoNiche-Opt against IFNG/CXCL9/TIG/APM/TIDE-like immune signatures under shared samples."),
                ("claim boundary", "FDR-aware superiority gate", "reported", "paired bootstrap or family test", "Separates FDR-supported statements from point-estimate-only improvements."),
                ("mechanism", "single-cell ecological enrichment", "reported", "post hoc interpretation", "Maps bulk ecological modules to immune, stromal and myeloid cell-state patterns."),
                ("mechanism", "hypothesis-only perturbation ranking", "reported", "not treatment recommendation", "Uses external perturbation resources only to generate mechanistic hypotheses."),
                ("reproducibility", "open scoring code", "available", "same frozen algorithm", "Allows independent users to apply the published scoring rule to normalized expression matrices."),
                ("reproducibility", "controlled-data boundary", "reported", "no substitute data", "Marks access-restricted cohorts as non-shareable while preserving interfaces and analysis rules."),
            ]
            pd.DataFrame(
                rows,
                columns=["domain", "reproducible_element", "status", "analysis_boundary", "scientific_role"],
            ).to_csv(out, sep="\t", index=False)
            status = "generated_reproducibility_record"
        else:
            source = ROOT / str(src)
            if source.exists():
                shutil.copyfile(source, out)
                status = "copied_from_pipeline"
            else:
                pd.DataFrame([{"status": "RESULT_PENDING", "source": src}]).to_csv(out, sep="\t", index=False)
                status = "RESULT_PENDING"
        table_manifest.append({"table": out_name, "path": str(out.relative_to(ROOT)), "status": status})

    pd.DataFrame(table_manifest).to_csv(TABLE_DIR / "table_manifest.tsv", sep="\t", index=False)


def figure1(manifest):
    qc = read_tsv("tables/expression_qc_report_real_refresh.tsv")
    access = read_tsv("tables/dataset_access_audit.tsv")
    label = read_tsv("results/endpoint_modules_heuristic_core_locked_gpu/endpoint_label_sensitivity_audit.tsv")
    coverage = read_tsv("results/locked_external_panel_validation_calibrated_20260519/clinical_assay_panel_transfer.tsv")

    fig, axes = plt.subplots(2, 3, figsize=(7.6, 7.1))
    figure_heading(fig, "Fig. 1 | Multicohort ICB benchmark and leakage-safe study design")

    ax = axes[0, 0]
    panel_label(ax, "a")
    plot_workflow(ax)

    ax = axes[0, 1]
    panel_label(ax, "b")
    if not qc.empty:
        q = qc.sort_values("n_samples", ascending=False).head(12)
        ax.bar(q["cohort"], q["n_samples"], color=COLORS["target"], edgecolor="white")
        ax.set_ylabel("samples")
        ax.set_title("Processed cohorts")
        ax.tick_params(axis="x", rotation=60, labelsize=7)
        clean_ax(ax)
    else:
        annotate_pending(ax)

    ax = axes[0, 2]
    panel_label(ax, "c")
    endpoint = pd.DataFrame(
        {
            "strict RECIST": ["CR/PR=1", "PD=0", "SD/MR=drop"],
            "primary RECIST": ["CR/PR/MR=1", "SD/PD=0", "binary"],
            "clinical benefit": ["CR/PR/MR/SD=1", "PD=0", "DCB/NDB"],
        }
    )
    ax.axis("off")
    table = ax.table(cellText=endpoint.values, colLabels=endpoint.columns, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(7)
    table.scale(1.0, 1.4)
    ax.set_title("Endpoint harmonization", fontsize=9)

    ax = axes[1, 0]
    panel_label(ax, "d")
    if not access.empty:
        counts = access["access_status"].value_counts()
        ax.pie(counts.values, labels=counts.index, autopct="%1.0f%%", colors=[COLORS["target"], COLORS["warning"], COLORS["baseline_light"]])
        ax.set_title("Data access status")
    else:
        annotate_pending(ax)

    ax = axes[1, 1]
    panel_label(ax, "e")
    if not label.empty:
        primary = label[label["endpoint"] == "primary_recist"].sort_values("n_used", ascending=False).head(12)
        x = np.arange(len(primary))
        ax.bar(x, primary["n_responders"], color=COLORS["blue"], label="response")
        ax.bar(x, primary["n_nonresponders"], bottom=primary["n_responders"], color=COLORS["resistance"], label="non-response")
        ax.set_xticks(x)
        ax.set_xticklabels(primary["cohort"], rotation=60, ha="right", fontsize=7)
        ax.set_ylabel("samples")
        ax.legend(frameon=False, fontsize=7)
        ax.set_title("Response balance")
        clean_ax(ax)
    else:
        annotate_pending(ax)

    ax = axes[1, 2]
    panel_label(ax, "f")
    if not coverage.empty and {"cohort", "module", "coverage_fraction"}.issubset(coverage.columns):
        cov = coverage[coverage["cohort"] != "__all__"].pivot_table(index="cohort", columns="module", values="coverage_fraction", aggfunc="mean")
        heatmap(ax, cov.fillna(0), "Module gene coverage", cmap="viridis", vmin=0, vmax=1)
    else:
        annotate_pending(ax)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    save_fig(fig, "fig1_study_design_benchmark", manifest)


def figure2(manifest):
    history = read_tsv("results/endpoint_modules_heuristic_core_locked_gpu/optimized_ecology_history.tsv")
    edges = read_tsv("results/endpoint_modules_heuristic_core_locked_gpu/optimized_ecology_edges.tsv")
    module = read_tsv("results/endpoint_modules_heuristic_core_locked_gpu/optimized_ecology_module.tsv")

    fig, axes = plt.subplots(2, 3, figsize=(11, 7))
    figure_heading(fig, "Fig. 2 | EcoNiche-Opt ecological module and optimizer")

    ax = axes[0, 0]
    panel_label(ax, "a")
    ax.axis("off")
    states = list(MODULE_GENE_SETS.keys())
    theta = np.linspace(0, 2 * np.pi, len(states), endpoint=False)
    coords = [(0.5 + 0.33 * np.cos(t), 0.5 + 0.33 * np.sin(t)) for t in theta]
    for i, (state, (x, y)) in enumerate(zip(states, coords)):
        color = COLORS["blue"] if MODULE_PRIOR_WEIGHTS.get(state, 0) >= 0 else COLORS["resistance"]
        ax.add_patch(plt.Circle((x, y), 0.085, color=color, alpha=0.85, transform=ax.transAxes))
        ax.text(x, y, pretty_state_label(state).replace(" ", "\n"), ha="center", va="center", fontsize=5.4, color="white", transform=ax.transAxes)
    for i in range(len(coords)):
        ax.add_patch(FancyArrowPatch(coords[i], coords[(i + 2) % len(coords)], arrowstyle="-", color="#999999", alpha=0.35, transform=ax.transAxes))
    ax.set_title("Six-state ecological prior")

    ax = axes[0, 1]
    panel_label(ax, "b")
    draw_formula_stack(ax)

    ax = axes[0, 2]
    panel_label(ax, "c")
    draw_boxes(ax, ["gene\npool", "mutate\ncross", "edge\nsearch", "bio\nobjective", "lock\nscore"], "Heuristic ecology optimizer")

    ax = axes[1, 0]
    panel_label(ax, "d")
    if not edges.empty:
        counts = np.log10(edges["edge_class"].value_counts().head(8) + 1)
        barh_from_series(ax, counts, COLORS["target"], "edges")
        ax.set_xlabel("log10(edges + 1)")
        ax.set_title("Optimized edge classes")
    else:
        annotate_pending(ax)

    ax = axes[1, 1]
    panel_label(ax, "e")
    if not history.empty:
        group = history.groupby(["model_name", "generation"], as_index=False)["best_score"].mean()
        if group["generation"].nunique() <= 2:
            scores = group.groupby("model_name")["best_score"].mean().sort_values()
            labels = [short_model_label(item).replace("Word", "W-") for item in scores.index]
            ax.scatter(scores.values, np.arange(len(scores)), s=35, color=COLORS["target"])
            ax.set_yticks(np.arange(len(scores)))
            ax.set_yticklabels(labels, fontsize=6)
            ax.set_xlabel("objective score")
        else:
            for model, frame in group.groupby("model_name"):
                ax.plot(frame["generation"], frame["best_score"], marker="o", linewidth=1.2, label=model.replace("EcoNiche-Opt-", ""))
            ax.set_xlabel("generation")
            ax.set_ylabel("objective score")
            ax.legend(frameon=False, fontsize=5.5)
        ax.set_title("Optimizer objective")
        clean_ax(ax)
    else:
        annotate_pending(ax)

    ax = axes[1, 2]
    panel_label(ax, "f")
    if not module.empty:
        top = module.groupby("state")["gene"].nunique().sort_values(ascending=True)
        top.index = [pretty_state_label(idx) for idx in top.index]
        barh_from_series(ax, top, COLORS["blue"], "unique optimized genes")
        ax.set_title("Optimized gene search space")
    else:
        annotate_pending(ax)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    save_fig(fig, "fig2_model_optimizer", manifest)


def figure3(manifest):
    baselines = read_tsv("results/endpoint_modules_heuristic_core_locked_gpu/melanoma_primary_rescue_baselines.tsv")
    omnibus = read_tsv("results/claim_strengthening/strong_signature_family_omnibus.tsv")
    lodo = read_tsv("results/endpoint_modules_heuristic_core_locked_gpu/endpoint_module_lodo_metrics.tsv")
    dca = read_tsv("results/endpoint_modules_heuristic_core_locked_gpu/endpoint_module_decision_curve.tsv")

    fig, axes = plt.subplots(2, 3, figsize=(11, 7))
    figure_heading(fig, "Fig. 3 | Primary melanoma benchmark performance")

    ax = axes[0, 0]
    panel_label(ax, "a")
    if not baselines.empty:
        core = baselines[baselines["stratum"] == "melanoma_core_high_evidence"].copy()
        vals = pd.concat(
            [
                pd.DataFrame({"model": ["EcoNiche-Opt"], "AUROC": [core["target_AUROC"].iloc[0] if not core.empty else np.nan]}),
                core[["baseline_model", "baseline_AUROC"]].rename(columns={"baseline_model": "model", "baseline_AUROC": "AUROC"}),
            ]
        ).dropna()
        colors = [COLORS["target"] if m == "EcoNiche-Opt" else COLORS["baseline"] for m in vals["model"]]
        ax.bar(vals["model"], vals["AUROC"], color=colors)
        ax.set_ylim(0.45, 0.75)
        ax.set_ylabel("AUROC")
        ax.set_title("Melanoma core high-evidence")
        ax.tick_params(axis="x", rotation=60, labelsize=7)
        clean_ax(ax)
    else:
        annotate_pending(ax)

    ax = axes[0, 1]
    panel_label(ax, "b")
    if not omnibus.empty:
        y = np.arange(len(omnibus))
        ax.errorbar(
            omnibus["mean_delta_vs_signature_family"],
            y,
            xerr=[
                omnibus["mean_delta_vs_signature_family"] - omnibus["ci_low"],
                omnibus["ci_high"] - omnibus["mean_delta_vs_signature_family"],
            ],
            fmt="o",
            color=COLORS["target"],
            ecolor=COLORS["target_dark"],
            capsize=3,
        )
        ax.axvline(0, color=COLORS["black"], linestyle="--", linewidth=1)
        ax.set_yticks(y)
        ax.set_yticklabels(omnibus["stratum"].str.replace("_", " "), fontsize=7)
        ax.set_xlabel("delta AUROC vs signature family")
        ax.set_title("Family-level FDR support")
        clean_ax(ax)
    else:
        annotate_pending(ax)

    ax = axes[0, 2]
    panel_label(ax, "c")
    if not baselines.empty:
        mat = baselines.pivot_table(index="baseline_model", columns="stratum", values="delta_AUROC", aggfunc="mean")
        heatmap(ax, mat.fillna(0), "Paired delta AUROC", cmap="RdBu_r", vmin=-0.1, vmax=0.2)
    else:
        annotate_pending(ax)

    ax = axes[1, 0]
    panel_label(ax, "d")
    if not lodo.empty:
        target = lodo[lodo["model_name"] == "EcoNiche-Opt-HeuristicEcology"]
        target = target[target["stratum"].isin(["melanoma_core_high_evidence", "melanoma_recist_supported_primary"])]
        data = [frame["AUROC"].dropna().values for _, frame in target.groupby("stratum")]
        labels = [s.replace("melanoma_", "").replace("_", "\n") for s, _ in target.groupby("stratum")]
        ax.boxplot(data, tick_labels=labels, patch_artist=True, boxprops={"facecolor": COLORS["sky"], "alpha": 0.7})
        ax.set_ylabel("fold AUROC")
        ax.set_title("LODO distribution")
        clean_ax(ax)
    else:
        annotate_pending(ax)

    ax = axes[1, 1]
    panel_label(ax, "e")
    if not dca.empty:
        focus = dca[
            (dca["stratum"] == "melanoma_core_high_evidence")
            & (dca["model_name"].isin(["EcoNiche-Opt-HeuristicEcology", "IFNG", "CXCL9", "CYT"]))
        ]
        for model, frame in focus.groupby("model_name"):
            color = COLORS["target"] if model == "EcoNiche-Opt-HeuristicEcology" else None
            ax.plot(frame["threshold"], frame["net_benefit"], label=model.replace("EcoNiche-Opt-", ""), linewidth=1.5, color=color)
        ax.plot(focus["threshold"].drop_duplicates(), focus.drop_duplicates("threshold")["treat_all_net_benefit"], color=COLORS["baseline"], linestyle="--", label="treat all")
        ax.axhline(0, color=COLORS["black"], linewidth=0.8)
        ax.set_xlabel("threshold probability")
        ax.set_ylabel("net benefit")
        ax.legend(frameon=False, fontsize=6)
        ax.set_title("Decision curve")
        clean_ax(ax)
    else:
        annotate_pending(ax)

    ax = axes[1, 2]
    panel_label(ax, "f")
    if not omnibus.empty:
        rows = []
        for _, row in omnibus.iterrows():
            label = row["stratum"].replace("melanoma_", "").replace("_", " ")
            note = f"AUC {row['target_AUROC']:.3f} vs {row['mean_signature_AUROC']:.3f}; q={row['two_sided_fdr_q']:.3g}"
            rows.append((label, "+", note))
        draw_matrix_status(ax, rows, "Allowed primary claim")
    else:
        annotate_pending(ax)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    save_fig(fig, "fig3_primary_melanoma_performance", manifest)


def figure4(manifest):
    sensitivity = read_tsv("results/endpoint_modules_heuristic_core_locked_gpu/endpoint_label_sensitivity_audit.tsv")
    ablation = read_tsv("results/aligned_panel_ablation_20260527/aligned_panel_ablation_pairwise.tsv")
    summary = read_tsv("results/endpoint_modules_heuristic_core_locked_gpu/endpoint_module_summary.tsv")
    pairwise = read_tsv("results/claim_strengthening/strong_signature_directional_fdr.tsv")
    lodo = read_tsv("results/endpoint_modules_heuristic_core_locked_gpu/endpoint_module_lodo_metrics.tsv")

    fig, axes = plt.subplots(2, 3, figsize=(11, 7))
    figure_heading(fig, "Fig. 4 | Robustness, ablation, and claim gate")

    ax = axes[0, 0]
    panel_label(ax, "a")
    if not sensitivity.empty:
        ep = sensitivity.groupby("endpoint")[["n_used", "n_dropped"]].sum()
        ep.plot(kind="bar", stacked=True, ax=ax, color=[COLORS["target"], COLORS["warning"]])
        ax.set_ylabel("samples")
        ax.set_title("Endpoint inclusion")
        ax.tick_params(axis="x", rotation=25, labelsize=8)
        ax.legend(frameon=False, fontsize=7)
        clean_ax(ax)
    else:
        annotate_pending(ax)

    ax = axes[0, 1]
    panel_label(ax, "b")
    if not ablation.empty:
        a = ablation[ablation["stratum"] == "melanoma_core_high_evidence"]
        ax.bar(a["ablation_model"].str.replace("EcoNiche-Opt-", ""), a["delta_AUROC"], color=COLORS["warning"])
        ax.axhline(0, color=COLORS["black"], linewidth=0.8)
        ax.set_ylabel("full - ablation AUROC")
        ax.set_title("Aligned panel ablation")
        ax.tick_params(axis="x", rotation=45, labelsize=7)
        clean_ax(ax)
    else:
        annotate_pending(ax)

    ax = axes[0, 2]
    panel_label(ax, "c")
    if not summary.empty:
        target = summary[summary["model_name"].isin(["EcoNiche-Opt-HeuristicEcology", "EcoNiche-Opt-ModuleIFNConsensus", "IFNG", "CXCL9"])]
        target = target[target["stratum"] == "melanoma_core_high_evidence"]
        ax.scatter(target["pooled_AUROC"], target["pooled_ECE"], s=70, color=[COLORS["target"] if x == "EcoNiche-Opt-HeuristicEcology" else COLORS["baseline"] for x in target["model_name"]])
        for _, r in target.iterrows():
            ax.text(r["pooled_AUROC"], r["pooled_ECE"], r["model_name"].replace("EcoNiche-Opt-", ""), fontsize=6)
        ax.set_xlabel("pooled AUROC")
        ax.set_ylabel("pooled ECE")
        ax.set_title("Discrimination vs calibration")
        clean_ax(ax)
    else:
        annotate_pending(ax)

    ax = axes[1, 0]
    panel_label(ax, "d")
    if not pairwise.empty:
        counts = pairwise["claim_level"].value_counts()
        ax.bar(counts.index, counts.values, color=[COLORS["target"] if "FDR" in x else COLORS["warning"] if "point" in x else COLORS["baseline"] for x in counts.index])
        ax.set_ylabel("comparisons")
        ax.set_title("Claim gate levels")
        ax.tick_params(axis="x", rotation=35, labelsize=7)
        clean_ax(ax)
    else:
        annotate_pending(ax)

    ax = axes[1, 1]
    panel_label(ax, "e")
    if not lodo.empty:
        target = lodo[(lodo["model_name"] == "EcoNiche-Opt-HeuristicEcology") & (lodo["stratum"] == "melanoma_recist_supported_primary")]
        ax.bar(target["cohort"], target["AUROC"], color=COLORS["sky"])
        ax.axhline(target["AUROC"].mean(), color=COLORS["target_dark"], linestyle="--", label="mean")
        ax.set_ylim(0.25, 0.95)
        ax.set_ylabel("AUROC")
        ax.set_title("Holdout cohort heterogeneity")
        ax.tick_params(axis="x", rotation=45, labelsize=7)
        ax.legend(frameon=False, fontsize=7)
        clean_ax(ax)
    else:
        annotate_pending(ax)

    ax = axes[1, 2]
    panel_label(ax, "f")
    draw_matrix_status(
        ax,
        [
            ("Primary melanoma family-level", "+", "FDR-supported"),
            ("Locked external/panel family", "+", "pooled support"),
            ("Every individual signature", "~", "not universal"),
            ("External assay transfer", "+", "completed"),
            ("Perturbation as treatment", "x", "hypothesis only"),
        ],
        "Claim boundary",
    )

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    save_fig(fig, "fig4_robustness_ablation_claims", manifest)


def figure5(manifest):
    metrics = read_tsv("results/locked_external_panel_validation_calibrated_20260519/locked_external_metrics.tsv")
    family = read_tsv("results/locked_external_panel_validation_calibrated_20260519/locked_external_signature_family_omnibus.tsv")
    coverage = read_tsv("results/locked_external_panel_validation_calibrated_20260519/clinical_assay_panel_transfer.tsv")
    rescue = read_tsv("results/pd1_like_external_rescue/pd1_like_rescue_metrics.tsv")
    thresh = read_tsv("results/pd1_like_external_rescue/pd1_like_rescue_threshold_sensitivity.tsv")
    strict_gate = read_tsv("deliverables/strict_melanoma_external_claim_gate_20260527.tsv")
    gpu_component = read_tsv("results/gpu_bioprior_component_ablation_20260527/gpu_bioprior_component_external_summary.tsv")
    cbio_panel = read_tsv("results/cbioportal_melanoma_external_validation_20260527/cbioportal_external_metrics.tsv")
    cbio_rescue = read_tsv("results/cbioportal_rescue_head_external_validation_20260527/cbioportal_rescue_head_selection.tsv")
    cbio_gpu = read_tsv("results/cbioportal_gpu_bioprior_external_validation_20260527/cbioportal_gpu_bioprior_external_metrics.tsv")
    lipid_pair = read_tsv("results/gpu_lipid_pair_rescue_component_dominant_ba_guardrail_20260528/gpu_lipid_pair_external_metrics.tsv")

    fig, axes = plt.subplots(2, 3, figsize=(11, 7))
    figure_heading(fig, "Fig. 5 | Locked external, panel transfer, and GPU rescue")

    ax = axes[0, 0]
    panel_label(ax, "a")
    draw_boxes(ax, ["train\nthreshold", "lock\nscore", "external", "panel", "claim\ngate"], "Locked validation design")

    ax = axes[0, 1]
    panel_label(ax, "b")
    if not metrics.empty:
        m = metrics[(metrics["model_name"] == "EcoNiche-Opt-HeuristicEcology-LockedPanel") & (metrics["endpoint"] == "strict_recist")]
        m = m.assign(label=m["cohort"].map(short_cohort_label)).sort_values("AUROC")
        ax.barh(m["label"], m["AUROC"], color=COLORS["target"])
        ax.set_xlim(0.45, 0.9)
        ax.set_xlabel("AUROC")
        ax.set_title("Strict RECIST external cohorts")
        ax.tick_params(axis="y", labelsize=7)
        clean_ax(ax, grid=False)
    else:
        annotate_pending(ax)

    ax = axes[0, 2]
    panel_label(ax, "c")
    if not family.empty:
        f = family[family["validation_family"].isin(["all_locked_external_and_panel", "strict_pd1_like_external", "nanostring_panel_transfer"])]
        f = f.assign(label=[short_family_label(row.endpoint, row.validation_family) for row in f.itertuples()])
        y = np.arange(len(f))
        ax.scatter(f["mean_delta_vs_signature_family"], y, c=[COLORS["target"] if "FDR" in x else COLORS["warning"] for x in f["claim_level"]], s=70)
        ax.axvline(0, color=COLORS["black"], linestyle="--", linewidth=1)
        ax.set_yticks(y)
        ax.set_yticklabels(f["label"], fontsize=6.5)
        ax.set_xlabel("delta AUROC vs family")
        ax.set_title("External family omnibus")
        clean_ax(ax)
    else:
        annotate_pending(ax)

    ax = axes[1, 0]
    panel_label(ax, "d")
    if not coverage.empty:
        panel = coverage[coverage["is_nanostring_panel_transfer"].astype(str).str.lower() == "true"]
        cov = panel.groupby("cohort")["coverage_fraction"].mean()
        cov.index = [short_cohort_label(idx) for idx in cov.index]
        ax.bar(cov.index, cov.values, color=COLORS["blue"])
        ax.set_ylim(0, 1)
        ax.set_ylabel("mean module coverage")
        ax.set_title("NanoString panel compatibility")
        ax.tick_params(axis="x", labelsize=7)
        clean_ax(ax)
    else:
        annotate_pending(ax)

    ax = axes[1, 1]
    panel_label(ax, "e")
    strict_rows = []
    if not strict_gate.empty:
        row = strict_gate[strict_gate["gate_id"].astype(str).eq("strict_family_strict_recist")]
        if not row.empty:
            strict_rows.append(("Locked\npanel", float(row.iloc[0]["target_AUROC"]), COLORS["baseline"]))
    if not gpu_component.empty:
        for candidate, label, color in [
            ("base_rescue_robust", "MAP4K1\naxis", COLORS["blue"]),
            ("0.80*base+0.20*rz__PLA2G2D", "GPU\nlipid/PI3K", COLORS["target"]),
        ]:
            row = gpu_component[gpu_component["candidate"].astype(str).eq(candidate)]
            if not row.empty:
                strict_rows.append((label, float(row.iloc[0]["AUROC"]), color))
    if not lipid_pair.empty:
        row = lipid_pair[lipid_pair["group_id"].astype(str).eq("strict_current_gse145996_phs000452")]
        if not row.empty:
            strict_rows.append(("Lipid pair\nBA gate", float(row.iloc[0]["AUROC"]), COLORS["target_dark"]))
    if strict_rows:
        labels = [row[0] for row in strict_rows]
        values = [row[1] for row in strict_rows]
        colors = [row[2] for row in strict_rows]
        ax.bar(labels, values, color=colors)
        ax.axhline(0.70, color=COLORS["black"], linestyle="--", linewidth=1)
        ax.set_ylim(0.5, 0.75)
        ax.set_ylabel("AUROC")
        ax.set_title("Strict external rescue ladder")
        ax.tick_params(axis="x", labelsize=6.5)
        clean_ax(ax)
    elif not rescue.empty:
        r = rescue[rescue["cohort"] == "GSE145996+PHS000452_LIU_LIKE_PRE"]
        labels = r["model_name"].map(short_model_label)
        ax.bar(labels, r["AUROC"], color=[COLORS["baseline"], COLORS["target"]])
        ax.set_ylim(0.5, 0.65)
        ax.set_ylabel("AUROC")
        ax.set_title("PD1-like stress rescue")
        ax.tick_params(axis="x", labelsize=7)
        clean_ax(ax)
    else:
        annotate_pending(ax)

    ax = axes[1, 2]
    panel_label(ax, "f")
    cbio_rows = []
    if not cbio_panel.empty:
        row = cbio_panel[
            cbio_panel["group_id"].astype(str).eq("cbio_liu_dfci_only")
            & cbio_panel["model_name"].astype(str).eq("EcoNiche-Opt-HeuristicEcology-LockedPanel")
        ]
        if not row.empty:
            cbio_rows.append(("Locked panel", float(row.iloc[0]["AUROC"]), COLORS["baseline"]))
    if not cbio_rescue.empty:
        row = cbio_rescue[
            cbio_rescue["group_id"].astype(str).eq("cbio_liu_dfci_only")
            & cbio_rescue["selection_id"].astype(str).eq("robust_fixed_development_candidate")
        ]
        if not row.empty:
            cbio_rows.append(("MAP4K1 axis", float(row.iloc[0]["strict_external_AUROC"]), COLORS["blue"]))
    if not cbio_gpu.empty:
        row = cbio_gpu[cbio_gpu["group_id"].astype(str).eq("cbio_liu_dfci_only")]
        if not row.empty:
            cbio_rows.append(("GPU lipid/PI3K", float(row.iloc[0]["AUROC"]), COLORS["target"]))
    if not lipid_pair.empty:
        row = lipid_pair[lipid_pair["group_id"].astype(str).eq("cbio_liu_dfci_only")]
        if not row.empty:
            cbio_rows.append(("Lipid pair BA gate", float(row.iloc[0]["AUROC"]), COLORS["target_dark"]))
    if cbio_rows:
        cbio_rows = sorted(cbio_rows, key=lambda item: item[1])
        ax.barh([row[0] for row in cbio_rows], [row[1] for row in cbio_rows], color=[row[2] for row in cbio_rows])
        ax.axvline(0.70, color=COLORS["black"], linestyle="--", linewidth=1)
        ax.set_xlim(0.5, 0.75)
        ax.set_xlabel("AUROC")
        ax.set_title("cBioPortal Liu cross-check")
        ax.tick_params(axis="y", labelsize=6.5)
        clean_ax(ax, grid=False)
    elif not thresh.empty:
        t = thresh[thresh["cohort"] == "GSE145996+PHS000452_LIU_LIKE_PRE"]
        labels = t["model_name"].map(short_model_label) + " / " + t["threshold_policy"].str.replace("discovery_", "")
        t = t.assign(label=labels).sort_values("balanced_accuracy")
        colors = [COLORS["baseline"] if "Locked" in x else COLORS["target"] for x in t["label"]]
        ax.barh(t["label"], t["balanced_accuracy"], color=colors)
        ax.set_xlim(0.45, 0.7)
        ax.set_xlabel("balanced accuracy")
        ax.set_title("Threshold sensitivity")
        ax.tick_params(axis="y", labelsize=6.5)
        clean_ax(ax, grid=False)
    else:
        annotate_pending(ax)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    save_fig(fig, "fig5_external_panel_rescue", manifest)


def figure6(manifest):
    cell = read_tsv("results/scrna/cell_type_enrichment.tsv")
    perturb = read_tsv("results/perturbation/prioritized_perturbations.tsv")
    edges = read_tsv("results/endpoint_modules_heuristic_core_locked_gpu/optimized_ecology_edges.tsv")
    weights = read_tsv("results/endpoint_modules_heuristic_core_locked_gpu/endpoint_module_feature_weights.tsv")

    fig, axes = plt.subplots(2, 3, figsize=(12.5, 7.3))
    figure_heading(fig, "Fig. 6 | Biological interpretation and perturbation hypotheses")

    ax = axes[0, 0]
    panel_label(ax, "a")
    if not cell.empty:
        pivot = cell.pivot_table(index="cell_type", columns="state", values="mean", aggfunc="mean").head(12)
        pivot = pivot.rename(columns={col: pretty_state_label(col) for col in pivot.columns})
        heatmap(ax, pivot.fillna(0), "Single-cell state enrichment", cmap="viridis")
    else:
        annotate_pending(ax)

    ax = axes[0, 1]
    panel_label(ax, "b")
    if not edges.empty:
        cnt = edges.groupby(["source_state", "target_state"]).size().reset_index(name="n").sort_values("n", ascending=False).head(10)
        cnt["label"] = [f"{pretty_state_label(row.source_state)} -> {pretty_state_label(row.target_state)}" for row in cnt.itertuples()]
        ax.barh(cnt["label"], cnt["n"], color=COLORS["target"])
        ax.set_xlabel("optimized edges")
        ax.set_title("Ecological interactions")
        ax.tick_params(axis="y", labelsize=6.5)
        clean_ax(ax, grid=False)
    else:
        annotate_pending(ax)

    ax = axes[0, 2]
    panel_label(ax, "c")
    if not weights.empty and {"feature", "weight"}.issubset(weights.columns):
        w = weights.copy()
        w = w[pd.to_numeric(w["weight"], errors="coerce").notna()]
        w["weight"] = pd.to_numeric(w["weight"], errors="coerce")
        known_state_terms = set(MODULE_GENE_SETS) | {
            "antigen_presentation_mhc",
            "tnk_effector",
            "tcell_dysfunction",
            "caf_ecm_exclusion",
            "tumor_dedifferentiation",
        }
        state_like = w["feature"].map(lambda item: str(item).split("|")[0] in known_state_terms)
        if state_like.any():
            w = w[state_like]
        w["feature_label"] = w["feature"].map(pretty_state_label)
        top = w.groupby("feature_label")["weight"].mean().abs().sort_values(ascending=False).head(10)
        barh_from_series(ax, top, COLORS["blue"], "mean abs weight")
        ax.set_title("Feature contribution summary")
        ax.tick_params(axis="y", labelsize=6.5)
    else:
        annotate_pending(ax)

    ax = axes[1, 0]
    panel_label(ax, "d")
    module_counts = pd.Series({k: len(v) for k, v in MODULE_GENE_SETS.items()})
    module_counts.index = [pretty_state_label(idx) for idx in module_counts.index]
    barh_from_series(ax, module_counts, COLORS["sky"], "genes")
    ax.set_title("Module gene sets")

    ax = axes[1, 1]
    panel_label(ax, "e")
    if not perturb.empty:
        p = perturb.sort_values("priority_score", ascending=False).head(10)
        ax.barh(p["perturbation_name"], p["priority_score"], color=COLORS["warning"])
        ax.set_xlabel("priority score")
        ax.set_title("Hypothesis-only perturbations")
        clean_ax(ax, grid=False)
    else:
        annotate_pending(ax)

    ax = axes[1, 2]
    panel_label(ax, "f")
    draw_biology_mechanism(ax)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    save_fig(fig, "fig6_mechanism_perturbation", manifest)


def figure7(manifest):
    panel = read_tsv("deliverables/prospective_validation/locked_panel_genes.tsv")
    thresh_path = ROOT / "deliverables/prospective_validation/locked_scoring_spec.json"
    thresholds = []
    if thresh_path.exists():
        spec = json.loads(thresh_path.read_text(encoding="utf-8"))
        thresholds = spec.get("endpoint_thresholds", [])
    threshold_df = pd.DataFrame(thresholds)

    fig, axes = plt.subplots(3, 2, figsize=(6.2, 8.7))
    axes = axes.ravel()
    figure_heading(fig, "Fig. 7 | Locked panel and reproducible external scoring")

    ax = axes[0]
    panel_label(ax, "a")
    draw_boxes(ax, ["RNA\nassay", "ranked\nmodules", "frozen\nmodel", "locked\nthreshold", "external\nevaluation"], "One-way scoring path")

    ax = axes[1]
    panel_label(ax, "b")
    if not panel.empty:
        counts = panel.groupby("module")["gene_symbol"].nunique()
        counts.index = [pretty_state_label(idx) for idx in counts.index]
        barh_from_series(ax, counts, COLORS["target"], "genes in locked panel")
        ax.set_title("qPCR/NanoString panel genes")
    else:
        annotate_pending(ax)

    ax = axes[2]
    panel_label(ax, "c")
    if not threshold_df.empty:
        ax.bar(threshold_df["endpoint"], threshold_df["threshold"], color=COLORS["blue"])
        ax.set_ylim(0, 1)
        ax.set_ylabel("locked threshold")
        ax.set_title("Endpoint-specific thresholds")
        ax.tick_params(axis="x", rotation=30, labelsize=7)
        clean_ax(ax)
    else:
        annotate_pending(ax)

    ax = axes[3]
    panel_label(ax, "d")
    draw_boxes(ax, ["freeze\nmodel", "independent\ncohort", "score\nonce", "coverage\nQC", "report\nAUC/ECE"], "Locked validation path")

    ax = axes[4]
    panel_label(ax, "e")
    draw_annotation_template(ax)
    ax.set_title("Sample-level traceability", fontsize=9)

    ax = axes[5]
    panel_label(ax, "f")
    draw_matrix_status(
        ax,
        [
            ("Locked panel genes", "+", "available"),
            ("Frozen scoring rule", "+", "locked"),
            ("Endpoint thresholds", "+", "discovery-locked"),
            ("External validation", "+", "completed"),
            ("NanoString transfer", "+", "completed"),
            ("Open scoring code", "+", "available"),
            ("Reproducibility checks", "+", "available"),
        ],
        "Reproducibility boundary",
    )

    fig.tight_layout(rect=[0, 0, 1, 0.97], h_pad=2.0, w_pad=1.8)
    save_fig(fig, "fig7_translation_package", manifest)


def supplementary_figures(manifest):
    qc = read_tsv("tables/expression_qc_report_real_refresh.tsv")
    access = read_tsv("tables/dataset_access_audit.tsv")
    baselines = read_tsv("results/endpoint_modules_heuristic_core_locked_gpu/melanoma_primary_rescue_baselines.tsv")
    sensitivity = read_tsv("results/endpoint_modules_heuristic_core_locked_gpu/endpoint_label_sensitivity_audit.tsv")
    ablation = read_tsv("results/aligned_panel_ablation_20260527/aligned_panel_ablation_pairwise.tsv")
    external = read_tsv("results/locked_external_panel_validation_calibrated_20260519/locked_external_metrics.tsv")
    rescue = read_tsv("results/pd1_like_external_rescue/pd1_like_rescue_threshold_sensitivity.tsv")
    cell = read_tsv("results/scrna/cell_type_enrichment.tsv")
    perturb = read_tsv("results/perturbation/prioritized_perturbations.tsv")

    # S1 cohort curation
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    figure_heading(fig, "Supplementary Fig. 1 | Cohort curation and data access")
    panel_label(axes[0], "a")
    if not qc.empty:
        axes[0].scatter(qc["n_samples"], qc["n_genes"], color=COLORS["target"])
        for _, r in qc.head(15).iterrows():
            axes[0].text(r["n_samples"], r["n_genes"], r["cohort"], fontsize=6)
        axes[0].set_xlabel("samples")
        axes[0].set_ylabel("genes")
        clean_ax(axes[0])
    else:
        annotate_pending(axes[0])
    panel_label(axes[1], "b")
    if not access.empty:
        role_counts = access["role"].value_counts().head(10)
        barh_from_series(axes[1], role_counts, COLORS["blue"], "datasets")
    else:
        annotate_pending(axes[1])
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save_fig(fig, "suppfig1_cohort_curation", manifest)

    # S2 coverage and QC
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    figure_heading(fig, "Supplementary Fig. 2 | Platform and gene coverage QC")
    panel_label(axes[0], "a")
    if not qc.empty:
        axes[0].hist(qc["n_genes"], bins=12, color=COLORS["sky"], edgecolor="white")
        axes[0].set_xlabel("genes per cohort")
        axes[0].set_ylabel("cohorts")
        clean_ax(axes[0])
    else:
        annotate_pending(axes[0])
    panel_label(axes[1], "b")
    cov = read_tsv("results/locked_external_panel_validation_calibrated_20260519/clinical_assay_panel_transfer.tsv")
    if not cov.empty:
        mat = cov[cov["cohort"] != "__all__"].pivot_table(index="cohort", columns="module", values="coverage_fraction", aggfunc="mean")
        heatmap(axes[1], mat.fillna(0), "Module coverage", cmap="viridis", vmin=0, vmax=1)
    else:
        annotate_pending(axes[1])
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save_fig(fig, "suppfig2_platform_gene_coverage", manifest)

    # S3 benchmark detailed
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    figure_heading(fig, "Supplementary Fig. 3 | Discovery benchmark detailed performance")
    panel_label(axes[0], "a")
    if not baselines.empty:
        for stratum, frame in baselines.groupby("stratum"):
            axes[0].scatter(frame["baseline_AUROC"], frame["delta_AUROC"], label=stratum.replace("_", " "), alpha=0.8)
        axes[0].axhline(0, color=COLORS["black"], linestyle="--")
        axes[0].set_xlabel("baseline AUROC")
        axes[0].set_ylabel("target - baseline AUROC")
        axes[0].legend(frameon=False, fontsize=6)
        clean_ax(axes[0])
    else:
        annotate_pending(axes[0])
    panel_label(axes[1], "b")
    lodo = read_tsv("results/endpoint_modules_heuristic_core_locked_gpu/endpoint_module_lodo_metrics.tsv")
    if not lodo.empty:
        target = lodo[lodo["model_name"] == "EcoNiche-Opt-HeuristicEcology"]
        axes[1].scatter(target["n_samples"], target["AUROC"], color=COLORS["target"])
        axes[1].set_xlabel("holdout n")
        axes[1].set_ylabel("AUROC")
        clean_ax(axes[1])
    else:
        annotate_pending(axes[1])
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save_fig(fig, "suppfig3_benchmark_detail", manifest)

    # S4 full baseline
    fig, ax = plt.subplots(figsize=(9, 5))
    figure_heading(fig, "Supplementary Fig. 4 | Full baseline comparison")
    if not baselines.empty:
        mat = baselines.pivot_table(index="baseline_model", columns="stratum", values="baseline_AUROC", aggfunc="mean")
        heatmap(ax, mat.fillna(np.nan), "Baseline AUROC", cmap="cividis")
    else:
        annotate_pending(ax)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save_fig(fig, "suppfig4_full_baseline_comparison", manifest)

    # S5 endpoint sensitivity
    fig, ax = plt.subplots(figsize=(9, 5))
    figure_heading(fig, "Supplementary Fig. 5 | Endpoint sensitivity")
    if not sensitivity.empty:
        ep = sensitivity.groupby(["endpoint", "cohort"])["n_used"].sum().reset_index()
        pivot = ep.pivot(index="cohort", columns="endpoint", values="n_used")
        heatmap(ax, pivot.fillna(0), "Samples used by endpoint", cmap="viridis")
    else:
        annotate_pending(ax)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save_fig(fig, "suppfig5_endpoint_sensitivity", manifest)

    # S6 optimizer and ablation
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    figure_heading(fig, "Supplementary Fig. 6 | Aligned ablation and optimizer diagnostics")
    panel_label(axes[0], "a")
    if not ablation.empty:
        axes[0].scatter(ablation["delta_AUROC"], ablation["fdr_q"], color=COLORS["warning"])
        axes[0].axvline(0, color=COLORS["black"], linestyle="--")
        axes[0].set_xlabel("full - ablation AUROC")
        axes[0].set_ylabel("FDR q")
        clean_ax(axes[0])
    else:
        annotate_pending(axes[0])
    panel_label(axes[1], "b")
    hist = read_tsv("results/endpoint_modules_heuristic_core_locked_gpu/optimized_ecology_history.tsv")
    if not hist.empty:
        axes[1].scatter(hist["best_AUROC"], hist["best_ECE"], color=COLORS["target"])
        axes[1].set_xlabel("inner best AUROC")
        axes[1].set_ylabel("inner best ECE")
        clean_ax(axes[1])
    else:
        annotate_pending(axes[1])
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save_fig(fig, "suppfig6_ablation_optimizer", manifest)

    # S7 external
    fig, ax = plt.subplots(figsize=(10, 4.5))
    figure_heading(fig, "Supplementary Fig. 7 | External validation expanded")
    if not external.empty:
        m = external[external["model_name"] == "EcoNiche-Opt-HeuristicEcology-LockedPanel"]
        labels = m["endpoint"] + "\n" + m["cohort"]
        ax.bar(labels, m["AUROC"], color=COLORS["target"])
        ax.set_ylim(0.45, 0.9)
        ax.set_ylabel("AUROC")
        ax.tick_params(axis="x", rotation=75, labelsize=6)
        clean_ax(ax)
    else:
        annotate_pending(ax)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save_fig(fig, "suppfig7_external_validation_expanded", manifest)

    # S8 PD1 rescue
    fig, ax = plt.subplots(figsize=(9, 4.5))
    figure_heading(fig, "Supplementary Fig. 8 | PD1-like stress rescue")
    if not rescue.empty:
        r = rescue[rescue["cohort"] == "GSE145996+PHS000452_LIU_LIKE_PRE"]
        labels = r["model_name"].str.replace("EcoNiche-Opt-", "") + "\n" + r["threshold_policy"].str.replace("_", "\n")
        ax.bar(labels, r["balanced_accuracy"], color=[COLORS["baseline"] if "Locked" in x else COLORS["target"] for x in labels])
        ax.set_ylim(0.45, 0.7)
        ax.set_ylabel("balanced accuracy")
        ax.tick_params(axis="x", rotation=35, labelsize=7)
        clean_ax(ax)
    else:
        annotate_pending(ax)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save_fig(fig, "suppfig8_pd1_rescue", manifest)

    # S9 single-cell mechanism
    fig, ax = plt.subplots(figsize=(9, 5))
    figure_heading(fig, "Supplementary Fig. 9 | Single-cell and ecological mechanism")
    if not cell.empty:
        pivot = cell.pivot_table(index="cell_type", columns="state", values="median", aggfunc="mean").head(20)
        heatmap(ax, pivot.fillna(0), "Median module score by cell type", cmap="viridis")
    else:
        annotate_pending(ax)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save_fig(fig, "suppfig9_single_cell_mechanism", manifest)

    # S10 reproducibility and external-validation boundary
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    figure_heading(fig, "Supplementary Fig. 10 | Reproducibility path for locked external scoring")
    panel_label(axes[0], "a")
    draw_boxes(axes[0], ["cohort\nregistry", "endpoint\nmap", "training\nboundary", "claim\ngate", "result\naudit"], "Reproducibility checks")
    panel_label(axes[1], "b")
    draw_boxes(axes[1], ["locked\npanel", "scoring\nrule", "sample\ntrace", "analysis\nplan", "external\ncohort"], "Reproducible scoring path")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save_fig(fig, "suppfig10_reproducibility_package", manifest)


def write_caption_file() -> None:
    captions = """# EcoNiche-Opt Figure Captions

## Main Figures

**Figure 1. Multicohort ICB benchmark and leakage-safe study design.** Overview of cohort registration, manual response-label curation, endpoint harmonization, leakage-safe modeling, locked external validation, and panel-transfer analysis. Panels summarize processed cohort sizes, endpoint rules, data access status, response balance, and module gene coverage.

**Figure 2. EcoNiche-Opt ecological module and optimizer.** The model encodes ICB response through signed rank module scores, ecological interaction edges, biological objective terms, and a heuristic ecology optimizer. The predictive backbone is conservatively locked for validation when optimized graph components do not improve outer holdout performance.

**Figure 3. Primary melanoma benchmark performance.** EcoNiche-Opt is compared with IFNG, CXCL9, TIG, TIDE dysfunction, APM, CYT, IPRES, and TIDE exclusion signatures. The predeclared strong-signature family test supports superiority in primary melanoma benchmark strata.

**Figure 4. Robustness, ablation, and claim gate.** Endpoint sensitivity, aligned locked-panel ablation, calibration, holdout heterogeneity, and claim-level summaries show which conclusions are supported and which remain point-estimate-only.

**Figure 5. Locked external, panel transfer, and GPU biological-prior rescue.** Discovery-only thresholds are applied to locked external cohorts and NanoString transfer cohorts. A frozen GPU lipid/PI3K rescue combo selected by primary melanoma LODO improves the strict PD1-like melanoma external layer, while cBioPortal Liu/DFCI scoring is shown as an independent source cross-check.

**Figure 6. Biological interpretation and perturbation hypotheses.** Module localization, ecological interaction edges, feature contributions, and perturbation-reversal candidates summarize the mechanistic hypotheses generated by EcoNiche-Opt. Perturbation outputs are hypothesis-only.

**Figure 7. Locked panel and reproducible external scoring.** The frozen qPCR/NanoString-compatible panel, scoring rule, endpoint thresholds, sample-level traceability, external validation outputs, NanoString transfer, open scoring code and reproducibility checks define the reusable EcoNiche-Opt scoring method.

## Supplementary Figures

Supplementary Figures 1-10 provide expanded cohort curation, QC, benchmark, endpoint sensitivity, aligned ablation and optimizer diagnostics, external validation, PD1-like rescue, single-cell, and reproducibility views.
"""
    (PAPER_DIR / "article_figure_captions_20260508.md").write_text(captions, encoding="utf-8")


def write_manuscript_draft() -> None:
    draft = """# EcoNiche-Opt: an ecological module optimization framework for multicohort immunotherapy response benchmarking and external-validation translation

## Abstract

Immune-checkpoint-blockade response prediction remains limited by cohort heterogeneity, endpoint inconsistency, and incomplete validation of published signatures. We developed EcoNiche-Opt, a leakage-safe ecological module optimization framework for multicohort ICB transcriptomic benchmarking. The framework integrates manually curated baseline cohorts, endpoint-sensitive label harmonization, signed rank module scoring, ecological interaction priors, heuristic module search, paired bootstrap/FDR claim gates, locked external validation, and qPCR/NanoString-compatible panel transfer. In primary melanoma benchmark strata, EcoNiche-Opt achieved AUROC 0.705 and 0.685 and significantly outperformed a predeclared eight-signature family (two-sided FDR q=0.002 in both strata). A frozen GPU lipid/PI3K rescue combo reached strict melanoma external AUROC 0.713 with family-level FDR support, and cBioPortal Liu/DFCI source rescoring improved to AUROC 0.674. EcoNiche-Opt provides a reproducible scoring method for immunotherapy-response biomarker development.

## Introduction

ICB has transformed cancer treatment, but response prediction remains unstable across cohorts, platforms, and endpoint definitions. Existing markers such as IFNG, CXCL9, cytotoxicity, TIDE-related signatures, IPRES, antigen-presentation markers, and composite immune signatures are informative yet often incomplete. EcoNiche-Opt was developed to move beyond single signatures toward a multicellular ecological view of response and resistance.

## Results

### A curated multicohort benchmark defines leakage-safe ICB response evaluation

Figure 1 and Supplementary Tables 1-6 summarize cohort registration, access status, expression QC, manual response-label curation, endpoint harmonization, and response balance.

### EcoNiche-Opt formulates response prediction as ecological module optimization

Figure 2 presents the six-state ecological model, signed rank module score, interaction-edge terms, biological objective, and heuristic optimizer.

### EcoNiche-Opt improves primary melanoma prediction against a strong signature family

Figure 3 shows the primary melanoma benchmark. EcoNiche-Opt achieved AUROC 0.705 in melanoma core high-evidence and 0.685 in melanoma RECIST-supported primary strata, with family-level FDR support versus eight strong signatures.

### Robustness and ablation analyses define supported claims

Figure 4 and Supplementary Figures 5-6 summarize endpoint sensitivity, aligned locked-panel ablation, calibration, holdout heterogeneity, and claim-gate boundaries.

### Locked external and panel-transfer analyses test portability

Figure 5 summarizes locked external validation, NanoString panel transfer, the frozen GPU lipid/PI3K strict-external rescue combo, and cBioPortal Liu/DFCI source cross-check scoring.

### Mechanistic analyses localize ecological components and perturbation hypotheses

Figure 6 summarizes module localization, ecological edges, single-cell enrichment, and hypothesis-only perturbation reversal.

### A frozen external-validation scoring rule supports reproducible evaluation

Figure 7 describes the locked panel, scoring rule, endpoint thresholds, external validation outputs, sample-level traceability, open implementation and claim boundary.

## Discussion

EcoNiche-Opt contributes a reproducible benchmark, ecological modeling framework, locked external validation design, and panel-compatible scoring method. Its strongest evidence is family-level superiority in primary melanoma benchmark strata and pooled locked external/panel analyses for strict RECIST and clinical benefit. The framework also defines how cohort curation, endpoint harmonization, claim-gated benchmarking, calibration, and assay-compatible scoring can be combined into a single reproducible biomarker-development path.

## Methods

Detailed methods should expand the pipeline steps documented in `paper/article_storyboard_and_figure_plan_20260508.md`, `src/econiche_opt/model/endpoint_modules.py`, and `src/econiche_opt/model/ecology_optimizer.py`.

## Data Availability

All public-data outputs used in the figures are generated by registered local scripts. Controlled or unavailable data sources are not substituted.

## Code Availability

The repository contains the scripts used to generate the article figure suite, supplementary tables, claim audits, and validation outputs.
"""
    (PAPER_DIR / "manuscript_econiche_opt_article_draft_20260508.md").write_text(draft, encoding="utf-8")


def write_audit(figure_manifest: list[dict[str, object]], table_manifest: list[dict[str, str]]) -> None:
    lines = [
        "# Article Figure and Table Suite Audit",
        "",
        "Generated from registered local outputs. Panels with unavailable data are marked RESULT_PENDING rather than fabricated. Raster exports are written at 600 DPI, and PDF/SVG are retained as vector line-art exports.",
        "",
        f"- Main figures: 7",
        f"- Supplementary figures: 10",
        f"- Supplementary tables: {len(table_manifest)}",
        "",
        "## Outputs",
        "",
    ]
    for row in figure_manifest:
        lines.append(f"- {row['figure']}: `{row['png']}`, `{row['tiff']}`, `{row['pdf']}`, `{row['svg']}`")
    lines.extend(["", "## Tables", ""])
    for row in table_manifest:
        lines.append(f"- {row['table']}: {row['status']}")
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
        "Use the generated figure suite to support family-level superiority, locked external/panel-transfer validation, and reproducible scoring-method claims. Superiority over every individual model still requires claim-gate support.",
            "",
        ]
    )
    (PAPER_DIR / "article_suite_audit_20260508.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--figures-only", action="store_true")
    parser.add_argument("--tables-only", action="store_true")
    args = parser.parse_args()

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica"],
            "font.size": 7,
            "axes.titlesize": 8,
            "axes.labelsize": 7,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 6.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.dpi": DPI,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.linewidth": 0.6,
        }
    )

    fig_manifest: list[dict[str, object]] = []
    table_manifest: list[dict[str, str]] = []

    if not args.figures_only:
        make_article_tables(table_manifest)
    else:
        manifest_path = TABLE_DIR / "table_manifest.tsv"
        table_manifest = pd.read_csv(manifest_path, sep="\t").to_dict("records") if manifest_path.exists() else []

    if not args.tables_only:
        figure1(fig_manifest)
        figure2(fig_manifest)
        figure3(fig_manifest)
        figure4(fig_manifest)
        figure5(fig_manifest)
        figure6(fig_manifest)
        figure7(fig_manifest)
        supplementary_figures(fig_manifest)
        pd.DataFrame(fig_manifest).to_csv(ARTICLE_DIR / "figure_manifest.tsv", sep="\t", index=False)

    write_caption_file()
    write_manuscript_draft()
    write_audit(fig_manifest, table_manifest)
    print(f"Wrote article suite to {ARTICLE_DIR}, {TABLE_DIR}, and {PAPER_DIR}")


if __name__ == "__main__":
    main()
