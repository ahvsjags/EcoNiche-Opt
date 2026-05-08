from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from econiche.baselines import BASELINE_SIGNATURES, signature_score
from econiche.metrics import compute_binary_metrics
from econiche.normalize import rank_gaussian_normalize
from econiche_opt.model.ecology_optimizer import (
    HeuristicEcologyConfig,
    HeuristicEcologyResult,
    build_ecology_features_from_module,
    optimize_ecology_module,
)


RESPONSE_HIGH_BASELINES = [
    "TIG",
    "TIDE_dysfunction",
    "IFNG",
    "PDCD1LG2",
    "CXCL9",
    "MCP_CD8_T",
    "APM",
    "PDL1_CD274",
    "IMPRES_template",
    "HLA_DRA",
    "CTLA4",
    "CYT",
    "PDCD1",
    "TLS",
]

NONRESPONSE_HIGH_BASELINES = [
    "MCP_fibroblast",
    "MPS",
    "IPRES",
    "C_ECM",
    "TIDE_exclusion",
    "ESCS",
]

STRONG_BASELINES = [
    "IFNG",
    "CXCL9",
    "TIG",
    "TIDE_dysfunction",
    "PDCD1LG2",
    "APM",
    "CYT",
    "IPRES",
    "TIDE_exclusion",
]

DIRECT_PROBABILITY_MODELS = {
    "EcoNiche-Opt-AdaptiveConsensus",
    "EcoNiche-Opt-ModuleIFNConsensus",
}

CALIBRATED_FIXED_SCORE_MODELS = {
    "EcoNiche-Opt-ModulePriorFixed": "EcoNiche-Opt-ModulePriorFixed-Platt",
}

MODULE_GENE_SETS: dict[str, list[str]] = {
    "ifn_t_cell_inflamed": [
        "IFNG",
        "CXCL9",
        "CXCL10",
        "CXCL11",
        "STAT1",
        "IDO1",
        "GBP1",
        "CXCR3",
        "CCL5",
        "CD274",
        "PDCD1LG2",
    ],
    "cytotoxic_cd8": ["CD8A", "CD8B", "GZMA", "GZMB", "GZMH", "PRF1", "NKG7", "GNLY"],
    "exhaustion_checkpoint": ["PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX", "CXCL13"],
    "antigen_presentation": ["HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9"],
    "myeloid_suppression": ["CD68", "CD163", "MRC1", "CSF1R", "ITGAM", "S100A8", "S100A9", "IL10", "TGFB1"],
    "stromal_exclusion": ["COL1A1", "COL1A2", "COL3A1", "FN1", "ACTA2", "VIM", "TGFBI", "POSTN", "LOXL2"],
    "trm_tls": ["ITGAE", "CD69", "CXCR6", "CXCL13", "ZNF683", "MS4A1", "CD79A", "BANK1", "LTB"],
}

MODULE_PRIOR_WEIGHTS: dict[str, float] = {
    "ifn_t_cell_inflamed": 1.0,
    "cytotoxic_cd8": 0.5,
    "exhaustion_checkpoint": 0.25,
    "antigen_presentation": 0.5,
    "myeloid_suppression": -0.5,
    "stromal_exclusion": -0.5,
    "trm_tls": 0.25,
}

WORD_FULL_GRAPH_MODEL = "EcoNiche-Opt-WordFullGraph"
OPTIMIZED_ADAPTIVE_MODEL = "EcoNiche-Opt-HeuristicEcology"
WORD_NO_INTERACTION_MODEL = "EcoNiche-Opt-WordNoInteraction"
WORD_UNSIGNED_GRAPH_MODEL = "EcoNiche-Opt-WordUnsignedGraph"
WORD_NO_BIO_OBJECTIVE_MODEL = "EcoNiche-Opt-WordNoBioObjective"

WORD_ABLATION_MODELS = [
    WORD_NO_INTERACTION_MODEL,
    WORD_UNSIGNED_GRAPH_MODEL,
    WORD_NO_BIO_OBJECTIVE_MODEL,
]

WORD_STATE_GENE_SETS: dict[str, list[str]] = {
    "tumor_dedifferentiation": ["AXL", "NGFR", "ITGA3", "VIM", "ZEB1", "MITF", "SOX10"],
    "antigen_presentation_mhc": ["HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9"],
    "tnk_effector": ["IFNG", "CXCL9", "CXCL10", "CXCL11", "CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "NKG7", "GNLY", "CCL5"],
    "tcell_dysfunction": ["PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX", "CXCL13", "IDO1"],
    "caf_ecm_exclusion": ["COL1A1", "COL1A2", "COL3A1", "FN1", "ACTA2", "VIM", "TGFBI", "POSTN", "LOXL2", "TGFB1", "CXCL12", "FAP"],
    "myeloid_suppression": ["CD68", "CD163", "MRC1", "CSF1R", "ITGAM", "S100A8", "S100A9", "IL10", "TGFB1"],
}

# (source_state, target_state, ligand_or_gene_a, receptor_or_gene_b, edge_class)
WORD_INTERACTION_EDGES: list[tuple[str, str, str, str, str]] = [
    ("caf_ecm_exclusion", "tnk_effector", "CXCL12", "CXCR4", "ligand_receptor"),
    ("caf_ecm_exclusion", "tnk_effector", "TGFB1", "TGFBR1", "ligand_receptor"),
    ("myeloid_suppression", "tnk_effector", "IL10", "IL10RA", "ligand_receptor"),
    ("antigen_presentation_mhc", "tnk_effector", "HLA-A", "B2M", "pathway"),
    ("tnk_effector", "tnk_effector", "GZMB", "PRF1", "pathway"),
    ("tcell_dysfunction", "tcell_dysfunction", "PDCD1", "LAG3", "checkpoint"),
    ("caf_ecm_exclusion", "caf_ecm_exclusion", "COL1A1", "FN1", "network"),
    ("myeloid_suppression", "myeloid_suppression", "S100A8", "S100A9", "network"),
    ("tnk_effector", "antigen_presentation_mhc", "IFNG", "HLA-A", "regulatory"),
]

MELANOMA_PRIMARY_COHORTS = [
    "GSE91061",
    "GSE78220",
    "GSE168204",
    "GSE115821",
    "GSE145996",
    "PRJEB23709_PD1_PRE",
]
MELANOMA_RECIST_SUPPORTED_COHORTS = ["GSE91061", "GSE78220", "GSE145996", "PRJEB23709_PD1_PRE"]
MELANOMA_CORE_HIGH_EVIDENCE_COHORTS = ["GSE91061", "GSE78220", "PRJEB23709_PD1_PRE"]
MELANOMA_CORE_WITH_PHS000452_COHORTS = ["GSE91061", "GSE78220", "PRJEB23709_PD1_PRE", "PHS000452_LIU_LIKE_PRE"]
MELANOMA_BINARY_RESPONSE_STRESS_COHORTS = ["GSE168204", "GSE115821"]
SECONDARY_CONFOUNDED_COHORTS = ["GSE165252", "PRJEB23709_COMBO_PRE", "PHS000452_LIU_LIKE_PRE"]


@dataclass(frozen=True)
class EndpointData:
    X_by_cohort: dict[str, pd.DataFrame]
    y_response_by_cohort: dict[str, pd.Series]
    metadata_by_cohort: dict[str, pd.DataFrame]


@dataclass(frozen=True)
class EvaluationResult:
    metrics: pd.DataFrame
    predictions: pd.DataFrame
    inner_selection: pd.DataFrame
    feature_weights: pd.DataFrame


def normalize_response(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip().upper()
    replacements = {
        "PRCR": "PR",
        "CRPR": "PR",
        "PR/CR": "PR",
        "CR/PR": "PR",
        "RESPONDER": "R",
        "RESPONSE": "R",
        "NONRESPONDER": "NR",
        "NON-RESPONDER": "NR",
        "NON_RESPONSE": "NR",
        "NO RESPONSE": "NR",
        "NONDURABLE_CLINICAL_BENEFIT": "NDB",
        "DURABLE_CLINICAL_BENEFIT": "DCB",
    }
    text = replacements.get(text, text)
    if "PR" in text and "CR" in text:
        return "PR"
    return text


def endpoint_response_label(response_raw: object, endpoint: str) -> float:
    token = normalize_response(response_raw)
    if endpoint not in {"strict_recist", "primary_recist", "clinical_benefit"}:
        raise ValueError(f"Unknown endpoint: {endpoint}")
    if token == "MR" and endpoint == "strict_recist":
        return np.nan
    if token in {"CR", "PR", "MR", "R", "DCB"}:
        return 1.0
    if token in {"PD", "NR", "NDB"}:
        return 0.0
    if token == "SD":
        if endpoint == "strict_recist":
            return np.nan
        if endpoint == "clinical_benefit":
            return 1.0
        return 0.0
    return np.nan


def endpoint_label_series(response_raw: pd.Series, endpoint: str) -> pd.Series:
    return response_raw.map(lambda value: endpoint_response_label(value, endpoint))


def prepare_endpoint_data(
    X_by_cohort: dict[str, pd.DataFrame],
    metadata_by_cohort: dict[str, pd.DataFrame],
    cohorts: Iterable[str],
    endpoint: str,
) -> EndpointData:
    endpoint_X: dict[str, pd.DataFrame] = {}
    endpoint_y: dict[str, pd.Series] = {}
    endpoint_meta: dict[str, pd.DataFrame] = {}
    for cohort in cohorts:
        if cohort not in X_by_cohort or cohort not in metadata_by_cohort:
            continue
        X = X_by_cohort[cohort]
        meta = metadata_by_cohort[cohort].reindex(X.index)
        if "response_raw" not in meta.columns:
            continue
        y = endpoint_label_series(meta["response_raw"], endpoint)
        mask = y.notna()
        if int(mask.sum()) < 4:
            continue
        X_endpoint = X.loc[mask].copy()
        y_endpoint = y.loc[mask].astype(int)
        meta_endpoint = meta.loc[mask].copy()
        if y_endpoint.nunique() < 2:
            continue
        endpoint_X[cohort] = X_endpoint
        endpoint_y[cohort] = y_endpoint.reindex(X_endpoint.index)
        endpoint_meta[cohort] = meta_endpoint.reindex(X_endpoint.index)
    return EndpointData(endpoint_X, endpoint_y, endpoint_meta)


def default_strata(active_cohorts: Iterable[str]) -> dict[str, dict[str, list[str]]]:
    active = list(active_cohorts)
    melanoma = [cohort for cohort in MELANOMA_PRIMARY_COHORTS if cohort in active]
    melanoma_recist = [cohort for cohort in MELANOMA_RECIST_SUPPORTED_COHORTS if cohort in active]
    melanoma_core = [cohort for cohort in MELANOMA_CORE_HIGH_EVIDENCE_COHORTS if cohort in active]
    melanoma_core_phs = [cohort for cohort in MELANOMA_CORE_WITH_PHS000452_COHORTS if cohort in active]
    melanoma_binary_stress = [cohort for cohort in MELANOMA_BINARY_RESPONSE_STRESS_COHORTS if cohort in active]
    secondary = [cohort for cohort in SECONDARY_CONFOUNDED_COHORTS if cohort in active]
    without_secondary = [cohort for cohort in active if cohort not in secondary]
    return {
        "melanoma_anti_pd1_primary": {"cohorts": melanoma, "train_pool": melanoma, "holdouts": melanoma},
        "melanoma_recist_supported_primary": {
            "cohorts": melanoma_recist,
            "train_pool": melanoma_recist,
            "holdouts": melanoma_recist,
        },
        "melanoma_core_high_evidence": {"cohorts": melanoma_core, "train_pool": melanoma_core, "holdouts": melanoma_core},
        "melanoma_core_plus_phs000452": {"cohorts": melanoma_core_phs, "train_pool": melanoma_core_phs, "holdouts": melanoma_core_phs},
        "melanoma_binary_response_stress": {
            "cohorts": melanoma_binary_stress,
            "train_pool": melanoma_binary_stress,
            "holdouts": melanoma_binary_stress,
        },
        "pan_cancer_response_all": {"cohorts": active, "train_pool": active, "holdouts": active},
        "pan_cancer_without_secondary": {
            "cohorts": without_secondary,
            "train_pool": without_secondary,
            "holdouts": without_secondary,
        },
        "secondary_confounded_transfer": {
            "cohorts": active,
            "train_pool": without_secondary,
            "holdouts": secondary,
        },
    }


def zscore_series(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    sd = values.std(ddof=0)
    if pd.isna(sd) or sd <= 0:
        return pd.Series(0.0, index=series.index)
    return ((values - values.mean()) / sd).fillna(0.0)


def sigmoid(values: pd.Series | np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    return 1.0 / (1.0 + np.exp(-np.clip(arr, -30, 30)))


def _available_genes(X: pd.DataFrame, genes: list[str]) -> list[str]:
    columns = set(X.columns)
    return [gene for gene in genes if gene in columns]


def _word_gene_universe() -> list[str]:
    genes: set[str] = set()
    for state_genes in WORD_STATE_GENE_SETS.values():
        genes.update(state_genes)
    for _, _, gene_a, gene_b, _ in WORD_INTERACTION_EDGES:
        genes.add(gene_a)
        genes.add(gene_b)
    return sorted(genes)


def _rank_expression_by_cohort(X_by_cohort: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    return {cohort: rank_gaussian_normalize(X.astype(float)) for cohort, X in X_by_cohort.items()}


def _estimate_word_gene_directions(
    ranked_expression_by_cohort: dict[str, pd.DataFrame],
    y_response: dict[str, pd.Series],
    train_cohorts: list[str],
) -> dict[str, int]:
    genes = _word_gene_universe()
    x_parts = []
    y_parts = []
    for cohort in train_cohorts:
        if cohort not in ranked_expression_by_cohort or cohort not in y_response:
            continue
        X = ranked_expression_by_cohort[cohort]
        y = y_response[cohort].astype(float)
        available = [gene for gene in genes if gene in X.columns]
        if not available:
            continue
        x_parts.append(X.reindex(y.index).loc[:, available])
        y_parts.append(y)
    if not x_parts:
        return {gene: 1 for gene in genes}
    X_train = pd.concat(x_parts, axis=0)
    y_train = pd.concat(y_parts, axis=0).reindex(X_train.index)
    directions: dict[str, int] = {}
    for gene in genes:
        if gene not in X_train.columns:
            directions[gene] = 1
            continue
        values = pd.to_numeric(X_train[gene], errors="coerce")
        valid = values.notna() & y_train.notna()
        if int(valid.sum()) < 4 or values.loc[valid].nunique() < 2 or y_train.loc[valid].nunique() < 2:
            directions[gene] = 1
            continue
        corr = values.loc[valid].corr(y_train.loc[valid])
        directions[gene] = -1 if pd.notna(corr) and corr < 0 else 1
    return directions


def _word_state_scores_from_ranked(
    ranked_X: pd.DataFrame,
    gene_directions: dict[str, int] | None = None,
    signed: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_directions = gene_directions or {}
    score_rows: dict[str, pd.Series] = {}
    coverage_rows: list[dict[str, object]] = []
    for state, genes in WORD_STATE_GENE_SETS.items():
        available = _available_genes(ranked_X, genes)
        if available:
            values = ranked_X[available].astype(float).copy()
            if signed:
                for gene in available:
                    values[gene] = values[gene] * int(gene_directions.get(gene, 1))
            score_rows[state] = values.sum(axis=1) / np.sqrt(len(available))
        else:
            score_rows[state] = pd.Series(0.0, index=ranked_X.index)
        coverage_rows.append(
            {
                "feature": state,
                "feature_type": "word_state",
                "state": state,
                "n_genes_defined": len(genes),
                "n_genes_available": len(available),
                "genes_available": ",".join(available),
            }
        )
    return pd.DataFrame(score_rows, index=ranked_X.index).fillna(0.0), pd.DataFrame(coverage_rows)


def _word_interaction_scores_from_ranked(
    ranked_X: pd.DataFrame,
    state_scores: pd.DataFrame,
    gene_directions: dict[str, int] | None = None,
    signed: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_directions = gene_directions or {}
    abundance = state_scores.apply(zscore_series, axis=0).apply(sigmoid, axis=0)
    grouped: dict[str, list[pd.Series]] = {}
    coverage_rows: list[dict[str, object]] = []
    for source_state, target_state, gene_a, gene_b, edge_class in WORD_INTERACTION_EDGES:
        feature = f"interaction__{source_state}__{target_state}"
        available = gene_a in ranked_X.columns and gene_b in ranked_X.columns
        if available:
            source_abundance = abundance[source_state] if source_state in abundance.columns else pd.Series(1.0, index=ranked_X.index)
            target_abundance = abundance[target_state] if target_state in abundance.columns else pd.Series(1.0, index=ranked_X.index)
            gene_a_values = ranked_X[gene_a].astype(float) * (int(gene_directions.get(gene_a, 1)) if signed else 1)
            gene_b_values = ranked_X[gene_b].astype(float) * (int(gene_directions.get(gene_b, 1)) if signed else 1)
            pair_score = gene_a_values * gene_b_values * source_abundance * target_abundance
            grouped.setdefault(feature, []).append(pair_score)
        coverage_rows.append(
            {
                "feature": feature,
                "feature_type": "word_interaction",
                "source_state": source_state,
                "target_state": target_state,
                "edge_class": edge_class,
                "gene_a": gene_a,
                "gene_b": gene_b,
                "edge_available": bool(available),
            }
        )
    expected_features = sorted({f"interaction__{source}__{target}" for source, target, *_ in WORD_INTERACTION_EDGES})
    scores: dict[str, pd.Series] = {}
    for feature in expected_features:
        if feature in grouped:
            raw = pd.concat(grouped[feature], axis=1).mean(axis=1)
            scores[feature] = zscore_series(raw)
        else:
            scores[feature] = pd.Series(0.0, index=ranked_X.index)
    return pd.DataFrame(scores, index=ranked_X.index).fillna(0.0), pd.DataFrame(coverage_rows)


def build_word_ecology_features(
    X: pd.DataFrame,
    gene_directions: dict[str, int] | None = None,
    signed: bool = True,
    include_interactions: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build the Word-spec state and interaction features for one cohort.

    The public helper is intentionally label-free; callers that need signed
    directions must pass directions estimated on training cohorts only.
    """
    ranked_X = rank_gaussian_normalize(X.astype(float))
    state_scores, state_coverage = _word_state_scores_from_ranked(ranked_X, gene_directions, signed=signed)
    frames = [state_scores]
    coverage = [state_coverage]
    if include_interactions:
        interaction_scores, interaction_coverage = _word_interaction_scores_from_ranked(
            ranked_X,
            state_scores,
            gene_directions=gene_directions,
            signed=signed,
        )
        frames.append(interaction_scores)
        coverage.append(interaction_coverage)
    return pd.concat(frames, axis=1).fillna(0.0), pd.concat(coverage, ignore_index=True)


def build_module_features(X: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    feature_rows: dict[str, pd.Series] = {}
    coverage_rows: list[dict[str, object]] = []
    for module_name, genes in MODULE_GENE_SETS.items():
        available = _available_genes(X, genes)
        if available:
            raw = X[available].astype(float).mean(axis=1)
            feature_rows[module_name] = zscore_series(raw)
        else:
            feature_rows[module_name] = pd.Series(0.0, index=X.index)
        coverage_rows.append(
            {
                "module": module_name,
                "n_genes_defined": len(genes),
                "n_genes_available": len(available),
                "genes_available": ",".join(available),
            }
        )
    return pd.DataFrame(feature_rows, index=X.index).fillna(0.0), pd.DataFrame(coverage_rows)


def build_module_features_by_cohort(
    X_by_cohort: dict[str, pd.DataFrame],
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    features: dict[str, pd.DataFrame] = {}
    coverage = []
    for cohort, X in X_by_cohort.items():
        frame, cohort_coverage = build_module_features(X)
        features[cohort] = frame
        cohort_coverage.insert(0, "cohort", cohort)
        coverage.append(cohort_coverage)
    coverage_df = pd.concat(coverage, ignore_index=True) if coverage else pd.DataFrame()
    return features, coverage_df


def module_prior_score(features: pd.DataFrame) -> pd.Series:
    score = pd.Series(0.0, index=features.index)
    for module, weight in MODULE_PRIOR_WEIGHTS.items():
        if module in features.columns:
            score = score + float(weight) * features[module]
    return score.fillna(0.0)


def current_three_gene_score(X: pd.DataFrame) -> pd.Series:
    scores = {}
    for name in ["IFNG", "CXCL9", "PDCD1LG2"]:
        genes = BASELINE_SIGNATURES.get(name, [name])
        scores[name] = zscore_series(signature_score(X, genes))
    frame = pd.DataFrame(scores, index=X.index).fillna(0.0)
    return ((frame["IFNG"] + frame["CXCL9"] + 2.0 * frame["PDCD1LG2"]) / 4.0).fillna(0.0)


def build_fixed_scores_by_cohort(
    X_by_cohort: dict[str, pd.DataFrame],
    module_features_by_cohort: dict[str, pd.DataFrame],
    baselines: Iterable[str] = STRONG_BASELINES,
) -> dict[str, dict[str, pd.Series]]:
    scores: dict[str, dict[str, pd.Series]] = {}
    for cohort, X in X_by_cohort.items():
        cohort_scores: dict[str, pd.Series] = {
            "EcoNiche-Opt-ImmuneComposite": current_three_gene_score(X),
            "EcoNiche-Opt-ModulePriorFixed": module_prior_score(module_features_by_cohort[cohort]),
        }
        for name in baselines:
            genes = BASELINE_SIGNATURES.get(name, [])
            raw = zscore_series(signature_score(X, genes))
            if name in NONRESPONSE_HIGH_BASELINES:
                raw = -raw
            cohort_scores[name] = raw.fillna(0.0)
        consensus_components = [
            "EcoNiche-Opt-ImmuneComposite",
            "EcoNiche-Opt-ModulePriorFixed",
            "IFNG",
            "TIG",
            "TIDE_dysfunction",
            "CXCL9",
        ]
        module_ifn_components = ["EcoNiche-Opt-ImmuneComposite", "EcoNiche-Opt-ModulePriorFixed", "IFNG"]
        cohort_scores["EcoNiche-Opt-AdaptiveConsensus"] = pd.DataFrame(
            {name: sigmoid(cohort_scores[name]) for name in consensus_components if name in cohort_scores},
            index=X.index,
        ).mean(axis=1)
        cohort_scores["EcoNiche-Opt-ModuleIFNConsensus"] = pd.DataFrame(
            {name: sigmoid(cohort_scores[name]) for name in module_ifn_components if name in cohort_scores},
            index=X.index,
        ).mean(axis=1)
        scores[cohort] = cohort_scores
    return scores


def select_threshold(y_true: np.ndarray, prob: np.ndarray) -> float:
    if len(y_true) == 0 or len(np.unique(y_true)) < 2:
        return 0.5
    thresholds = np.unique(np.quantile(prob, np.linspace(0.05, 0.95, 19)))
    best_threshold = 0.5
    best_score = -np.inf
    for threshold in thresholds:
        metrics = compute_binary_metrics(y_true, prob, threshold=float(threshold))
        score = float(metrics.get("balanced_accuracy", float("nan")))
        if np.isfinite(score) and score > best_score:
            best_score = score
            best_threshold = float(threshold)
    return best_threshold


def _concat(items: dict[str, pd.DataFrame | pd.Series], cohorts: list[str]) -> pd.DataFrame | pd.Series:
    return pd.concat([items[cohort] for cohort in cohorts], axis=0)


def _candidate_specs(optimize_word: bool = True) -> list[dict[str, object]]:
    specs: list[dict[str, object]] = [
        {"candidate": "module_prior_composite", "kind": "fixed_module_prior"},
        {"candidate": "module_positive_composite", "kind": "fixed_module_positive"},
        {
            "candidate": "word_full_graph",
            "kind": "word_ecology",
            "model_name": WORD_FULL_GRAPH_MODEL,
            "signed": True,
            "interactions": True,
            "bio_objective": True,
            "optimize_module": bool(optimize_word),
            "C": 0.5,
        },
        {
            "candidate": "word_no_interaction",
            "kind": "word_ecology",
            "model_name": WORD_NO_INTERACTION_MODEL,
            "signed": True,
            "interactions": False,
            "bio_objective": True,
            "optimize_module": bool(optimize_word),
            "C": 0.5,
        },
        {
            "candidate": "word_unsigned_graph",
            "kind": "word_ecology",
            "model_name": WORD_UNSIGNED_GRAPH_MODEL,
            "signed": False,
            "interactions": True,
            "bio_objective": True,
            "optimize_module": bool(optimize_word),
            "C": 0.5,
        },
        {
            "candidate": "word_no_bio_objective",
            "kind": "word_ecology",
            "model_name": WORD_NO_BIO_OBJECTIVE_MODEL,
            "signed": True,
            "interactions": True,
            "bio_objective": False,
            "optimize_module": bool(optimize_word),
            "C": 0.5,
        },
    ]
    for penalty in ["l1", "l2"]:
        for c_value in [0.05, 0.25, 1.0]:
            specs.append({"candidate": f"module_logistic_{penalty}_C{c_value:g}", "kind": "logistic", "penalty": penalty, "C": c_value})
    return specs


def _model_name_for_spec(spec: dict[str, object], selected: bool = False) -> str:
    if selected:
        return OPTIMIZED_ADAPTIVE_MODEL
    return str(spec.get("model_name", f"EcoNiche-Opt-{spec['candidate']}"))


def _module_positive_score(features: pd.DataFrame) -> pd.Series:
    positives = ["ifn_t_cell_inflamed", "cytotoxic_cd8", "exhaustion_checkpoint", "antigen_presentation", "trm_tls"]
    return features[[name for name in positives if name in features.columns]].mean(axis=1).fillna(0.0)


def _fit_score_calibrator(score: pd.Series, y: pd.Series) -> LogisticRegression | None:
    if len(y) < 4 or y.nunique() < 2:
        return None
    model = LogisticRegression(C=1.0, class_weight="balanced", solver="lbfgs", max_iter=5000, random_state=42)
    model.fit(np.asarray(score, dtype=float).reshape(-1, 1), y.astype(int).to_numpy())
    return model


def _predict_score_calibrator(model: LogisticRegression | None, score: pd.Series) -> np.ndarray:
    if model is None:
        return sigmoid(score)
    return model.predict_proba(np.asarray(score, dtype=float).reshape(-1, 1))[:, 1]


def _fit_monotone_score_calibrator(score: pd.Series, y: pd.Series) -> LogisticRegression | None:
    model = _fit_score_calibrator(score, y)
    if model is None or float(model.coef_[0, 0]) <= 0:
        return None
    return model


def _fit_module_logistic(X: pd.DataFrame, y: pd.Series, penalty: str, c_value: float) -> Pipeline:
    solver = "liblinear" if penalty == "l1" else "lbfgs"
    model = Pipeline(
        [
            ("scale", StandardScaler()),
            (
                "logistic",
                LogisticRegression(
                    C=float(c_value),
                    penalty=penalty,
                    solver=solver,
                    class_weight="balanced",
                    max_iter=5000,
                    random_state=42,
                ),
            ),
        ]
    )
    model.fit(X.astype(float), y.astype(int))
    return model


def _standardize_train_test_score(train_score: pd.Series, test_score: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    train_values = pd.to_numeric(train_score, errors="coerce").fillna(0.0).to_numpy(dtype=float)
    test_values = pd.to_numeric(test_score, errors="coerce").fillna(0.0).to_numpy(dtype=float)
    mean = float(np.mean(train_values)) if len(train_values) else 0.0
    sd = float(np.std(train_values)) if len(train_values) else 1.0
    if not np.isfinite(sd) or sd <= 0:
        sd = 1.0
    return (train_values - mean) / sd, (test_values - mean) / sd


def _default_optimizer_config(spec: dict[str, object] | None = None) -> HeuristicEcologyConfig:
    spec = spec or {}
    use_gpu_text = str(spec.get("optimizer_use_gpu", "0")).strip().lower()
    return HeuristicEcologyConfig(
        population_size=int(spec.get("optimizer_population", 4)),
        generations=int(spec.get("optimizer_generations", 2)),
        n_jobs=int(spec.get("optimizer_n_jobs", 1)),
        use_gpu=use_gpu_text not in {"0", "false", "no"},
        use_bio_objective=bool(spec.get("bio_objective", True)),
        random_state=int(spec.get("optimizer_random_state", 42)),
    )


def _append_metadata_offset_features(
    X_train_word: pd.DataFrame,
    X_test_word: pd.DataFrame,
    metadata_by_cohort: dict[str, pd.DataFrame] | None,
    train_cohorts: list[str],
    test_cohort: str,
) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, object]]]:
    if not metadata_by_cohort:
        return X_train_word, X_test_word, []
    categorical_specs = {
        "cancer_offset": ["cancer_type", "tumor_type", "tumor_of_origin", "histology", "disease_state"],
        "therapy_offset": ["treatment", "therapy", "io_therapy", "antibody", "treatment_group", "drug"],
    }
    train_meta = pd.concat(
        [metadata_by_cohort[cohort].reindex(X_train_word.index.intersection(metadata_by_cohort[cohort].index)) for cohort in train_cohorts if cohort in metadata_by_cohort],
        axis=0,
    )
    train_meta = train_meta.reindex(X_train_word.index)
    test_meta = metadata_by_cohort.get(test_cohort, pd.DataFrame(index=X_test_word.index)).reindex(X_test_word.index)
    train_out = X_train_word.copy()
    test_out = X_test_word.copy()
    rows: list[dict[str, object]] = []
    for prefix, columns in categorical_specs.items():
        available = [column for column in columns if column in train_meta.columns]
        if not available:
            continue
        train_values = train_meta[available].astype("string").bfill(axis=1).iloc[:, 0].fillna("unknown")
        test_available = [column for column in columns if column in test_meta.columns]
        if test_available:
            test_values = test_meta[test_available].astype("string").bfill(axis=1).iloc[:, 0].fillna("unknown")
        else:
            test_values = pd.Series("unknown", index=X_test_word.index, dtype="string")
        counts = train_values.value_counts(dropna=False)
        for category, count in counts.items():
            if int(count) < 4:
                continue
            category_text = str(category).strip().replace(" ", "_").replace("/", "_").replace(";", "_")
            if not category_text or category_text.lower() == "unknown":
                continue
            feature = f"{prefix}__{category_text[:48]}"
            train_feature = (train_values.astype(str) == str(category)).astype(float)
            if train_feature.nunique() < 2:
                continue
            train_out[feature] = train_feature.to_numpy(dtype=float)
            test_out[feature] = (test_values.astype(str) == str(category)).astype(float).to_numpy(dtype=float)
            rows.append({"feature": feature, "weight": 1.0, "weight_type": prefix})
    return train_out.fillna(0.0), test_out.fillna(0.0), rows


def _append_current_model_anchor_features(
    X_train_word: pd.DataFrame,
    X_test_word: pd.DataFrame,
    module_features: dict[str, pd.DataFrame],
    train_cohorts: list[str],
    test_cohort: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    X_train = X_train_word.copy()
    X_test = X_test_word.copy()
    train_module = _concat(module_features, train_cohorts).astype(float)
    test_module = module_features[test_cohort].astype(float)
    anchors = {
        "current_module_prior_anchor": (module_prior_score(train_module), module_prior_score(test_module)),
        "current_module_positive_anchor": (_module_positive_score(train_module), _module_positive_score(test_module)),
    }
    for name, (train_score, test_score) in anchors.items():
        X_train[name] = pd.to_numeric(train_score, errors="coerce").fillna(0.0).to_numpy(dtype=float)
        X_test[name] = pd.to_numeric(test_score, errors="coerce").fillna(0.0).to_numpy(dtype=float)
    return X_train.fillna(0.0), X_test.fillna(0.0)


def _select_word_bio_features(X_train: pd.DataFrame, y_train: pd.Series, use_bio_objective: bool) -> list[str]:
    if not use_bio_objective:
        return list(X_train.columns)
    anchor_features = [column for column in X_train.columns if column.startswith("current_")]
    state_features = [column for column in WORD_STATE_GENE_SETS if column in X_train.columns]
    interaction_features = [column for column in X_train.columns if column.startswith("interaction__")]
    selected = list(dict.fromkeys([*anchor_features, *state_features]))
    interaction_scores: list[tuple[float, str]] = []
    y = pd.Series(y_train.to_numpy(dtype=float), index=X_train.index)
    for feature in interaction_features:
        values = pd.to_numeric(X_train[feature], errors="coerce")
        valid = values.notna() & y.notna()
        if int(valid.sum()) < 6 or values.loc[valid].nunique() < 2 or y.loc[valid].nunique() < 2:
            continue
        corr = values.loc[valid].corr(y.loc[valid])
        if pd.notna(corr):
            interaction_scores.append((abs(float(corr)), feature))
    interaction_scores.sort(reverse=True)
    for score, feature in interaction_scores:
        if score >= 0.02:
            selected.append(feature)
    if interaction_scores and not any(feature.startswith("interaction__") for feature in selected):
        selected.append(interaction_scores[0][1])
    return list(dict.fromkeys([feature for feature in selected if feature in X_train.columns]))


def _training_feature_direction(values: pd.Series, y_train: pd.Series) -> int:
    y = pd.Series(y_train.to_numpy(dtype=float), index=values.index)
    numeric = pd.to_numeric(values, errors="coerce")
    valid = numeric.notna() & y.notna()
    if int(valid.sum()) < 6 or numeric.loc[valid].nunique() < 2 or y.loc[valid].nunique() < 2:
        return 1
    corr = numeric.loc[valid].corr(y.loc[valid])
    return -1 if pd.notna(corr) and corr < 0 else 1


def _word_fixed_composite_scores(
    X_train_word: pd.DataFrame,
    X_test_word: pd.DataFrame,
    y_train: pd.Series,
    selected_features: list[str],
    use_bio_objective: bool,
) -> tuple[pd.Series, pd.Series, list[dict[str, object]]]:
    train_score = pd.Series(0.0, index=X_train_word.index)
    test_score = pd.Series(0.0, index=X_test_word.index)
    rows: list[dict[str, object]] = []
    weights: dict[str, float] = {
        "current_module_prior_anchor": 1.0,
        "current_module_positive_anchor": 0.15,
    }
    for state in WORD_STATE_GENE_SETS:
        if state in selected_features:
            weights[state] = 0.08
    for feature in selected_features:
        if feature.startswith("interaction__"):
            weights[feature] = 0.005 if use_bio_objective else 0.05
    for feature, magnitude in weights.items():
        if feature not in X_train_word.columns or feature not in X_test_word.columns:
            continue
        direction = 1 if feature.startswith("current_module_prior_anchor") else _training_feature_direction(X_train_word[feature], y_train)
        signed_weight = float(magnitude * direction)
        train_score = train_score + signed_weight * X_train_word[feature].astype(float)
        test_score = test_score + signed_weight * X_test_word[feature].astype(float)
        rows.append({"feature": feature, "weight": signed_weight, "weight_type": "word_fixed_composite_weight"})
    return train_score.fillna(0.0), test_score.fillna(0.0), rows


def _word_feature_frame_from_ranked(
    ranked_X: pd.DataFrame,
    gene_directions: dict[str, int],
    signed: bool,
    include_interactions: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    state_scores, state_coverage = _word_state_scores_from_ranked(ranked_X, gene_directions, signed=signed)
    frames = [state_scores]
    coverage = [state_coverage]
    if include_interactions:
        interaction_scores, interaction_coverage = _word_interaction_scores_from_ranked(
            ranked_X,
            state_scores,
            gene_directions=gene_directions,
            signed=signed,
        )
        frames.append(interaction_scores)
        coverage.append(interaction_coverage)
    return pd.concat(frames, axis=1).fillna(0.0), pd.concat(coverage, ignore_index=True)


def _word_candidate_matrices(
    spec: dict[str, object],
    train_cohorts: list[str],
    test_cohort: str,
    ranked_expression_by_cohort: dict[str, pd.DataFrame],
    y_response: dict[str, pd.Series],
    matrix_cache: dict[tuple[object, ...], tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, int], HeuristicEcologyResult | None]] | None = None,
    metadata_by_cohort: dict[str, pd.DataFrame] | None = None,
    optimizer_config: HeuristicEcologyConfig | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, int], HeuristicEcologyResult | None]:
    signed = bool(spec.get("signed", True))
    include_interactions = bool(spec.get("interactions", True))
    optimize_module = bool(spec.get("optimize_module", False))
    cache_key = (str(spec.get("candidate")), tuple(train_cohorts), test_cohort, signed, include_interactions, optimize_module)
    if matrix_cache is not None and cache_key in matrix_cache:
        return matrix_cache[cache_key]
    if optimize_module:
        cfg = optimizer_config or _default_optimizer_config(spec)
        cfg = HeuristicEcologyConfig(
            min_genes_per_state=cfg.min_genes_per_state,
            max_genes_per_state=cfg.max_genes_per_state,
            candidate_pool_per_state=cfg.candidate_pool_per_state,
            global_candidate_genes=cfg.global_candidate_genes,
            population_size=cfg.population_size,
            generations=cfg.generations,
            elite_fraction=cfg.elite_fraction,
            mutation_rate=cfg.mutation_rate,
            crossover_rate=cfg.crossover_rate,
            robust_rho=cfg.robust_rho,
            max_auto_edges_per_state_pair=cfg.max_auto_edges_per_state_pair,
            min_edge_abs_corr=cfg.min_edge_abs_corr,
            n_jobs=cfg.n_jobs,
            random_state=cfg.random_state + abs(hash((str(spec.get("candidate")), tuple(train_cohorts)))) % 10000,
            use_gpu=cfg.use_gpu,
            use_bio_objective=bool(spec.get("bio_objective", True)),
        )
        optimizer_cache_key = ("optimizer", str(spec.get("candidate")), tuple(train_cohorts), signed, include_interactions, bool(spec.get("bio_objective", True)))
        optimizer_result = None
        if matrix_cache is not None and optimizer_cache_key in matrix_cache:
            optimizer_result = matrix_cache[optimizer_cache_key][4]
        if optimizer_result is None:
            optimizer_result = optimize_ecology_module(
                ranked_expression_by_cohort,
                y_response,
                train_cohorts,
                WORD_STATE_GENE_SETS,
                WORD_INTERACTION_EDGES,
                metadata_by_cohort=metadata_by_cohort,
                signed=signed,
                include_interactions=include_interactions,
                cfg=cfg,
            )
            if matrix_cache is not None:
                matrix_cache[optimizer_cache_key] = (pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), optimizer_result.gene_directions, optimizer_result)
        train_frames = []
        coverage_frames = []
        for cohort in train_cohorts:
            frame = build_ecology_features_from_module(
                ranked_expression_by_cohort[cohort],
                optimizer_result.genes_by_state,
                optimizer_result.edges,
                optimizer_result.gene_directions,
                signed=signed,
                include_interactions=include_interactions,
            )
            train_frames.append(frame.reindex(y_response[cohort].index))
            state_rows = []
            for state, genes in optimizer_result.genes_by_state.items():
                available = [gene for gene in genes if gene in ranked_expression_by_cohort[cohort].columns]
                state_rows.append(
                    {
                        "cohort": cohort,
                        "feature": state,
                        "feature_type": "word_state",
                        "state": state,
                        "n_genes_defined": len(genes),
                        "n_genes_available": len(available),
                        "genes_available": ",".join(available),
                    }
                )
            coverage_frames.append(pd.DataFrame(state_rows))
        test_frame = build_ecology_features_from_module(
            ranked_expression_by_cohort[test_cohort],
            optimizer_result.genes_by_state,
            optimizer_result.edges,
            optimizer_result.gene_directions,
            signed=signed,
            include_interactions=include_interactions,
        )
        edge_rows = [
            {
                "cohort": test_cohort,
                "feature": f"interaction__{row['source_state']}__{row['target_state']}",
                "feature_type": "word_interaction",
                "source_state": row["source_state"],
                "target_state": row["target_state"],
                "edge_class": row["edge_class"],
                "gene_a": row["gene_a"],
                "gene_b": row["gene_b"],
                "edge_available": True,
            }
            for _, row in optimizer_result.edge_table.iterrows()
        ]
        coverage_frames.append(pd.DataFrame(edge_rows))
        train_matrix = pd.concat(train_frames, axis=0).fillna(0.0)
        test_matrix = test_frame.reindex(y_response[test_cohort].index).fillna(0.0)
        all_columns = sorted(set(train_matrix.columns) | set(test_matrix.columns))
        result = (
            train_matrix.reindex(columns=all_columns, fill_value=0.0),
            test_matrix.reindex(columns=all_columns, fill_value=0.0),
            pd.concat([frame for frame in coverage_frames if not frame.empty], ignore_index=True) if coverage_frames else pd.DataFrame(),
            optimizer_result.gene_directions,
            optimizer_result,
        )
        if matrix_cache is not None:
            matrix_cache[cache_key] = result
        return result
    gene_directions = _estimate_word_gene_directions(ranked_expression_by_cohort, y_response, train_cohorts) if signed else {
        gene: 1 for gene in _word_gene_universe()
    }
    train_frames = []
    coverage_frames = []
    for cohort in train_cohorts:
        frame, coverage = _word_feature_frame_from_ranked(
            ranked_expression_by_cohort[cohort],
            gene_directions,
            signed=signed,
            include_interactions=include_interactions,
        )
        train_frames.append(frame.reindex(y_response[cohort].index))
        coverage.insert(0, "cohort", cohort)
        coverage_frames.append(coverage)
    test_frame, test_coverage = _word_feature_frame_from_ranked(
        ranked_expression_by_cohort[test_cohort],
        gene_directions,
        signed=signed,
        include_interactions=include_interactions,
    )
    test_coverage.insert(0, "cohort", test_cohort)
    coverage_frames.append(test_coverage)
    train_matrix = pd.concat(train_frames, axis=0).fillna(0.0)
    test_matrix = test_frame.reindex(y_response[test_cohort].index).fillna(0.0)
    all_columns = sorted(set(train_matrix.columns) | set(test_matrix.columns))
    result = (
        train_matrix.reindex(columns=all_columns, fill_value=0.0),
        test_matrix.reindex(columns=all_columns, fill_value=0.0),
        pd.concat(coverage_frames, ignore_index=True),
        gene_directions,
        None,
    )
    if matrix_cache is not None:
        matrix_cache[cache_key] = result
    return result


def _direction_stability(
    ranked_expression_by_cohort: dict[str, pd.DataFrame],
    y_response: dict[str, pd.Series],
    train_cohorts: list[str],
) -> float:
    scores = []
    for gene in _word_gene_universe():
        signs = []
        for cohort in train_cohorts:
            if cohort not in ranked_expression_by_cohort or gene not in ranked_expression_by_cohort[cohort].columns:
                continue
            y = y_response[cohort].astype(float)
            values = ranked_expression_by_cohort[cohort][gene].reindex(y.index).astype(float)
            valid = values.notna() & y.notna()
            if int(valid.sum()) < 4 or values.loc[valid].nunique() < 2 or y.loc[valid].nunique() < 2:
                continue
            corr = values.loc[valid].corr(y.loc[valid])
            if pd.notna(corr):
                signs.append(-1 if corr < 0 else 1)
        if signs:
            scores.append(max(signs.count(1), signs.count(-1)) / len(signs))
    return float(np.mean(scores)) if scores else 0.5


def _batch_dependence(feature_frame: pd.DataFrame, cohort_labels: pd.Series) -> float:
    ratios = []
    labels = cohort_labels.reindex(feature_frame.index)
    for column in feature_frame.columns:
        values = pd.to_numeric(feature_frame[column], errors="coerce")
        valid = values.notna() & labels.notna()
        if int(valid.sum()) < 4 or values.loc[valid].nunique() < 2 or labels.loc[valid].nunique() < 2:
            continue
        overall = float(values.loc[valid].var(ddof=0))
        if overall <= 0:
            continue
        ratios.append(float(values.loc[valid].groupby(labels.loc[valid]).mean().var(ddof=0) / (overall + 1e-9)))
    return float(np.clip(np.mean(ratios), 0.0, 1.0)) if ratios else 0.0


def _state_redundancy(feature_frame: pd.DataFrame) -> float:
    state_columns = [column for column in WORD_STATE_GENE_SETS if column in feature_frame.columns]
    if len(state_columns) < 2:
        return 0.0
    corr = feature_frame[state_columns].corr().abs().to_numpy(dtype=float)
    upper = corr[np.triu_indices_from(corr, k=1)]
    upper = upper[np.isfinite(upper)]
    return float(np.mean(upper)) if len(upper) else 0.0


def _therapy_confounding(feature_frame: pd.DataFrame, metadata: pd.DataFrame | None) -> float:
    if metadata is None or metadata.empty:
        return 0.0
    treatment_columns = ["treatment", "therapy", "io_therapy", "antibody", "treatment_group"]
    available = [column for column in treatment_columns if column in metadata.columns]
    if not available:
        return 0.0
    treatment = metadata[available].astype("string").bfill(axis=1).iloc[:, 0].reindex(feature_frame.index).astype("string")
    if treatment.nunique(dropna=True) < 2:
        return 0.0
    ratios = []
    for column in feature_frame.columns:
        values = pd.to_numeric(feature_frame[column], errors="coerce")
        valid = values.notna() & treatment.notna()
        if int(valid.sum()) < 4 or values.loc[valid].nunique() < 2:
            continue
        overall = float(values.loc[valid].var(ddof=0))
        if overall <= 0:
            continue
        ratios.append(float(values.loc[valid].groupby(treatment.loc[valid]).mean().var(ddof=0) / (overall + 1e-9)))
    return float(np.clip(np.mean(ratios), 0.0, 1.0)) if ratios else 0.0


def _word_biological_objective_terms(
    spec: dict[str, object],
    train_cohorts: list[str],
    ranked_expression_by_cohort: dict[str, pd.DataFrame] | None,
    y_response: dict[str, pd.Series],
    metadata_by_cohort: dict[str, pd.DataFrame] | None,
    matrix_cache: dict[tuple[object, ...], tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, int], HeuristicEcologyResult | None]] | None = None,
    objective_cache: dict[tuple[object, ...], dict[str, float]] | None = None,
    optimizer_config: HeuristicEcologyConfig | None = None,
) -> dict[str, float]:
    if ranked_expression_by_cohort is None or str(spec.get("kind")) != "word_ecology":
        return {}
    cache_key = (
        str(spec.get("candidate")),
        tuple(train_cohorts),
        bool(spec.get("signed", True)),
        bool(spec.get("interactions", True)),
        bool(spec.get("bio_objective", True)),
    )
    if objective_cache is not None and cache_key in objective_cache:
        return objective_cache[cache_key]
    train_matrix, _, coverage, _, optimizer_result = _word_candidate_matrices(
        spec,
        train_cohorts,
        train_cohorts[0],
        ranked_expression_by_cohort,
        y_response,
        matrix_cache=matrix_cache,
        metadata_by_cohort=metadata_by_cohort,
        optimizer_config=optimizer_config,
    )
    if optimizer_result is not None:
        result = {
            "bio_cell_specificity": float(optimizer_result.objective_terms.get("bio_cell_specificity", 0.0)),
            "bio_pathway": float(optimizer_result.objective_terms.get("bio_pathway", 0.0)),
            "bio_network": float(optimizer_result.objective_terms.get("bio_network", 0.0)),
            "bio_lr": float(optimizer_result.objective_terms.get("bio_lr", 0.0)),
            "bio_direction_stability": float(optimizer_result.objective_terms.get("bio_direction_stability", 0.5)),
            "penalty_size": float(optimizer_result.objective_terms.get("penalty_size", 0.0)),
            "penalty_batch": float(optimizer_result.objective_terms.get("penalty_batch", 0.0)),
            "penalty_leakage": 0.0,
            "penalty_redundancy": float(optimizer_result.objective_terms.get("penalty_redundancy", 0.0)),
            "penalty_therapy_confounding": float(optimizer_result.objective_terms.get("penalty_therapy_confounding", 0.0)),
            "bio_bonus": float(optimizer_result.objective_terms.get("bio_bonus", 0.0)),
            "bio_penalty": float(optimizer_result.objective_terms.get("bio_penalty", 0.0)),
            "bio_objective_delta": float(optimizer_result.objective_terms.get("bio_objective_delta", 0.0)),
        }
        if objective_cache is not None:
            objective_cache[cache_key] = result
        return result
    cohort_labels = pd.concat([pd.Series(cohort, index=y_response[cohort].index) for cohort in train_cohorts])
    metadata = None
    if metadata_by_cohort:
        metadata = pd.concat([metadata_by_cohort[cohort].reindex(y_response[cohort].index) for cohort in train_cohorts], axis=0)
    state_cov = coverage[coverage["feature_type"] == "word_state"].copy()
    if state_cov.empty:
        cell_specificity = 0.0
    else:
        cell_specificity = float((state_cov["n_genes_available"].astype(float) / state_cov["n_genes_defined"].astype(float).clip(lower=1)).mean())
    interaction_cov = coverage[coverage["feature_type"] == "word_interaction"].copy()
    lr = float(interaction_cov[interaction_cov["edge_class"] == "ligand_receptor"]["edge_available"].mean()) if not interaction_cov.empty else 0.0
    network = float(interaction_cov[interaction_cov["edge_class"].isin(["network", "regulatory", "checkpoint"])]["edge_available"].mean()) if not interaction_cov.empty else 0.0
    pathway = float(interaction_cov[interaction_cov["edge_class"] == "pathway"]["edge_available"].mean()) if not interaction_cov.empty else 0.0
    stability = _direction_stability(ranked_expression_by_cohort, y_response, train_cohorts)
    n_unique_genes = len(set().union(*(set(genes) for genes in WORD_STATE_GENE_SETS.values())))
    size = float(np.clip((n_unique_genes - 18) / (150 - 18), 0.0, 1.0))
    batch = _batch_dependence(train_matrix, cohort_labels)
    redundancy = _state_redundancy(train_matrix)
    therapy = _therapy_confounding(train_matrix, metadata)
    leakage = 0.0
    bio_bonus = 0.08 * cell_specificity + 0.04 * pathway + 0.04 * network + 0.04 * lr + 0.08 * stability
    penalty = 0.04 * size + 0.08 * batch + 0.05 * redundancy + 0.05 * therapy + leakage
    result = {
        "bio_cell_specificity": cell_specificity,
        "bio_pathway": pathway,
        "bio_network": network,
        "bio_lr": lr,
        "bio_direction_stability": stability,
        "penalty_size": size,
        "penalty_batch": batch,
        "penalty_leakage": leakage,
        "penalty_redundancy": redundancy,
        "penalty_therapy_confounding": therapy,
        "bio_bonus": float(bio_bonus),
        "bio_penalty": float(penalty),
        "bio_objective_delta": float(bio_bonus - penalty),
    }
    if objective_cache is not None:
        objective_cache[cache_key] = result
    return result


def _predict_candidate(
    spec: dict[str, object],
    train_cohorts: list[str],
    test_cohort: str,
    module_features: dict[str, pd.DataFrame],
    y_response: dict[str, pd.Series],
    ranked_expression_by_cohort: dict[str, pd.DataFrame] | None = None,
    metadata_by_cohort: dict[str, pd.DataFrame] | None = None,
    word_matrix_cache: dict[tuple[object, ...], tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, int], HeuristicEcologyResult | None]] | None = None,
    word_objective_cache: dict[tuple[object, ...], dict[str, float]] | None = None,
    optimizer_config: HeuristicEcologyConfig | None = None,
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    y_train = _concat(y_response, train_cohorts).astype(int)
    X_train = _concat(module_features, train_cohorts).astype(float)
    X_test = module_features[test_cohort].astype(float)
    kind = str(spec["kind"])
    weight_rows: list[dict[str, object]] = []
    if kind == "word_ecology":
        if ranked_expression_by_cohort is None:
            raise ValueError("ranked_expression_by_cohort is required for Word-spec ecological candidates")
        X_train_word, X_test_word, coverage, gene_directions, optimizer_result = _word_candidate_matrices(
            spec,
            train_cohorts,
            test_cohort,
            ranked_expression_by_cohort,
            y_response,
            matrix_cache=word_matrix_cache,
            metadata_by_cohort=metadata_by_cohort,
            optimizer_config=optimizer_config,
        )
        X_train_word, X_test_word = _append_current_model_anchor_features(
            X_train_word,
            X_test_word,
            module_features,
            train_cohorts,
            test_cohort,
        )
        metadata_rows: list[dict[str, object]] = []
        X_train_word, X_test_word, metadata_rows = _append_metadata_offset_features(
            X_train_word,
            X_test_word,
            metadata_by_cohort,
            train_cohorts,
            test_cohort,
        )
        weight_rows.extend(metadata_rows)
        selected_features = _select_word_bio_features(X_train_word, y_train, bool(spec.get("bio_objective", True)))
        if bool(spec.get("optimize_module", False)):
            selected_features = [feature for feature in selected_features if feature in X_train_word.columns and feature in X_test_word.columns]
            if not selected_features:
                selected_features = list(X_train_word.columns)
            model = _fit_module_logistic(
                X_train_word[selected_features].astype(float),
                y_train,
                "l2",
                float(spec.get("C", 0.5)),
            )
            train_prob = model.predict_proba(X_train_word[selected_features].astype(float))[:, 1]
            test_prob = model.predict_proba(X_test_word[selected_features].astype(float))[:, 1]
            logistic = model.named_steps["logistic"]
            for feature, weight in zip(selected_features, logistic.coef_[0]):
                weight_rows.append({"feature": feature, "weight": float(weight), "weight_type": "word_formula_theta_eta_offset"})
            weight_rows.append({"feature": "intercept", "weight": float(logistic.intercept_[0]), "weight_type": "word_formula_beta0"})
        else:
            train_score, test_score, composite_rows = _word_fixed_composite_scores(
                X_train_word,
                X_test_word,
                y_train,
                selected_features,
                bool(spec.get("bio_objective", True)),
            )
            train_prob = sigmoid(train_score)
            test_prob = sigmoid(test_score)
            weight_rows.extend(composite_rows)
        for feature in X_train_word.columns:
            weight_rows.append(
                {
                    "feature": feature,
                    "weight": 1.0 if feature in selected_features else 0.0,
                    "weight_type": "word_bio_objective_feature_selected" if bool(spec.get("bio_objective", True)) else "word_predictive_feature_available",
                }
            )
        for gene, direction in sorted(gene_directions.items()):
            weight_rows.append({"feature": gene, "weight": int(direction), "weight_type": "word_response_direction"})
        if optimizer_result is not None:
            for _, row in optimizer_result.gene_table.iterrows():
                weight_rows.append(
                    {
                        "feature": row["gene"],
                        "weight": float(row.get("selection_frequency", 0.0)),
                        "weight_type": "optimized_module_gene",
                        "state": row.get("state"),
                        "direction": row.get("direction"),
                        "training_abs_correlation": row.get("training_abs_correlation"),
                        "direction_stability": row.get("direction_stability"),
                        "state_prior_seed": row.get("state_prior_seed"),
                        "state_candidate_score": row.get("state_candidate_score"),
                    }
                )
            for _, row in optimizer_result.edge_table.iterrows():
                weight_rows.append(
                    {
                        "feature": row.get("edge_id"),
                        "weight": 1.0,
                        "weight_type": "optimized_interaction_edge",
                        "source_state": row.get("source_state"),
                        "target_state": row.get("target_state"),
                        "gene_a": row.get("gene_a"),
                        "gene_b": row.get("gene_b"),
                        "edge_class": row.get("edge_class"),
                    }
                )
            for _, row in optimizer_result.history.iterrows():
                weight_rows.append(
                    {
                        "feature": f"generation_{int(row['generation'])}",
                        "weight": float(row.get("best_score", np.nan)),
                        "weight_type": "optimizer_history",
                        "generation": int(row.get("generation", 0)),
                        "best_AUROC": row.get("best_AUROC"),
                        "best_AUPRC": row.get("best_AUPRC"),
                        "best_ECE": row.get("best_ECE"),
                        "backend": row.get("backend"),
                    }
                )
            weight_rows.append({"feature": "optimizer_backend", "weight": 1.0, "weight_type": optimizer_result.backend})
        for _, row in coverage.iterrows():
            weight_rows.append(
                {
                    "feature": row.get("feature"),
                    "weight": float(row.get("n_genes_available", 1.0)) if row.get("feature_type") == "word_state" else float(bool(row.get("edge_available", False))),
                    "weight_type": str(row.get("feature_type")),
                }
            )
        bio_terms = _word_biological_objective_terms(
            spec,
            train_cohorts,
            ranked_expression_by_cohort,
            y_response,
            metadata_by_cohort,
            matrix_cache=word_matrix_cache,
            objective_cache=word_objective_cache,
            optimizer_config=optimizer_config,
        )
        for term, value in bio_terms.items():
            weight_rows.append({"feature": term, "weight": float(value), "weight_type": "word_objective_term"})
        return train_prob, test_prob, pd.DataFrame(weight_rows)

    if kind in {"fixed_module_prior", "fixed_module_positive"}:
        if kind == "fixed_module_prior":
            train_score = module_prior_score(X_train)
            test_score = module_prior_score(X_test)
            for module, weight in MODULE_PRIOR_WEIGHTS.items():
                weight_rows.append({"feature": module, "weight": weight, "weight_type": "fixed_prior"})
        else:
            train_score = _module_positive_score(X_train)
            test_score = _module_positive_score(X_test)
            for module in ["ifn_t_cell_inflamed", "cytotoxic_cd8", "exhaustion_checkpoint", "antigen_presentation", "trm_tls"]:
                weight_rows.append({"feature": module, "weight": 1.0, "weight_type": "fixed_positive_mean"})
        train_prob = sigmoid(train_score)
        test_prob = sigmoid(test_score)
        weight_rows.append({"feature": "calibration", "weight": 0.0, "weight_type": "fixed_sigmoid_no_rank_changing_calibration"})
        return train_prob, test_prob, pd.DataFrame(weight_rows)

    model = _fit_module_logistic(X_train, y_train, str(spec["penalty"]), float(spec["C"]))
    train_prob = model.predict_proba(X_train)[:, 1]
    test_prob = model.predict_proba(X_test)[:, 1]
    logistic = model.named_steps["logistic"]
    coef = logistic.coef_[0]
    for feature, weight in zip(X_train.columns, coef):
        weight_rows.append({"feature": feature, "weight": float(weight), "weight_type": "scaled_logistic_coef"})
    weight_rows.append({"feature": "intercept", "weight": float(logistic.intercept_[0]), "weight_type": "scaled_logistic_intercept"})
    return train_prob, test_prob, pd.DataFrame(weight_rows)


def _score_inner_candidate(
    spec: dict[str, object],
    train_cohorts: list[str],
    module_features: dict[str, pd.DataFrame],
    y_response: dict[str, pd.Series],
    ranked_expression_by_cohort: dict[str, pd.DataFrame] | None = None,
    metadata_by_cohort: dict[str, pd.DataFrame] | None = None,
    word_matrix_cache: dict[tuple[object, ...], tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, int], HeuristicEcologyResult | None]] | None = None,
    word_objective_cache: dict[tuple[object, ...], dict[str, float]] | None = None,
    optimizer_config: HeuristicEcologyConfig | None = None,
) -> dict[str, object]:
    rows = []
    if len(train_cohorts) < 2:
        return {"candidate": spec["candidate"], "inner_mean_AUROC": np.nan, "inner_mean_AUPRC": np.nan, "inner_mean_balanced_accuracy": np.nan, "selection_score": -np.inf}
    if bool(spec.get("optimize_module", False)) and ranked_expression_by_cohort is not None:
        try:
            _, _, _, _, optimizer_result = _word_candidate_matrices(
                spec,
                train_cohorts,
                train_cohorts[0],
                ranked_expression_by_cohort,
                y_response,
                matrix_cache=word_matrix_cache,
                metadata_by_cohort=metadata_by_cohort,
                optimizer_config=optimizer_config,
            )
            if optimizer_result is not None:
                terms = optimizer_result.objective_terms
                return {
                    "candidate": spec["candidate"],
                    "inner_mean_AUROC": float(terms.get("inner_mean_AUROC", np.nan)),
                    "inner_mean_AUPRC": float(terms.get("inner_mean_AUPRC", np.nan)),
                    "inner_mean_balanced_accuracy": float(terms.get("inner_mean_balanced_accuracy", np.nan)),
                    "inner_mean_ECE": float(terms.get("inner_mean_ECE", np.nan)),
                    "inner_sd_AUROC": float(terms.get("inner_sd_AUROC", np.nan)),
                    "performance_score": float(terms.get("performance_score", np.nan)),
                    "bio_cell_specificity": float(terms.get("bio_cell_specificity", 0.0)),
                    "bio_pathway": float(terms.get("bio_pathway", 0.0)),
                    "bio_network": float(terms.get("bio_network", 0.0)),
                    "bio_lr": float(terms.get("bio_lr", 0.0)),
                    "bio_direction_stability": float(terms.get("bio_direction_stability", 0.5)),
                    "penalty_size": float(terms.get("penalty_size", 0.0)),
                    "penalty_batch": float(terms.get("penalty_batch", 0.0)),
                    "penalty_leakage": 0.0,
                    "penalty_redundancy": float(terms.get("penalty_redundancy", 0.0)),
                    "penalty_therapy_confounding": float(terms.get("penalty_therapy_confounding", 0.0)),
                    "bio_bonus": float(terms.get("bio_bonus", 0.0)),
                    "bio_penalty": float(terms.get("bio_penalty", 0.0)),
                    "bio_objective_delta": float(terms.get("bio_objective_delta", 0.0)),
                    "selection_score": float(terms.get("score", -np.inf)),
                }
        except Exception:
            return {"candidate": spec["candidate"], "inner_mean_AUROC": np.nan, "inner_mean_AUPRC": np.nan, "inner_mean_balanced_accuracy": np.nan, "selection_score": -np.inf}
    for inner_holdout in train_cohorts:
        inner_train = [cohort for cohort in train_cohorts if cohort != inner_holdout]
        if len(inner_train) == 0:
            continue
        y_train = _concat(y_response, inner_train).astype(int)
        y_test = y_response[inner_holdout].astype(int)
        if y_train.nunique() < 2 or y_test.nunique() < 2:
            continue
        try:
            train_prob, test_prob, _ = _predict_candidate(
                spec,
                inner_train,
                inner_holdout,
                module_features,
                y_response,
                ranked_expression_by_cohort=ranked_expression_by_cohort,
                metadata_by_cohort=metadata_by_cohort,
                word_matrix_cache=word_matrix_cache,
                word_objective_cache=word_objective_cache,
                optimizer_config=optimizer_config,
            )
        except Exception:
            continue
        threshold = select_threshold(y_train.to_numpy(dtype=int), np.asarray(train_prob, dtype=float))
        metrics = compute_binary_metrics(y_test, test_prob, threshold=threshold)
        rows.append(metrics)
    if not rows:
        return {"candidate": spec["candidate"], "inner_mean_AUROC": np.nan, "inner_mean_AUPRC": np.nan, "inner_mean_balanced_accuracy": np.nan, "selection_score": -np.inf}
    frame = pd.DataFrame(rows)
    if str(spec.get("kind")) == "word_ecology":
        performance_score = float(
            frame["AUROC"].mean()
            - 0.50 * frame["AUROC"].std(ddof=0)
            + 0.10 * frame["AUPRC"].mean()
            + 0.05 * frame["balanced_accuracy"].mean()
            - 0.15 * frame["ECE"].mean()
        )
        bio_terms = _word_biological_objective_terms(
            spec,
            train_cohorts,
            ranked_expression_by_cohort,
            y_response,
            metadata_by_cohort,
            matrix_cache=word_matrix_cache,
            objective_cache=word_objective_cache,
            optimizer_config=optimizer_config,
        )
        score = performance_score + (float(bio_terms.get("bio_objective_delta", 0.0)) if bool(spec.get("bio_objective", True)) else 0.0)
    else:
        bio_terms = {}
        performance_score = float(frame["AUROC"].mean() + 0.05 * frame["AUPRC"].mean() + 0.05 * frame["balanced_accuracy"].mean())
        score = performance_score
    return {
        "candidate": spec["candidate"],
        "inner_mean_AUROC": float(frame["AUROC"].mean()),
        "inner_mean_AUPRC": float(frame["AUPRC"].mean()),
        "inner_mean_balanced_accuracy": float(frame["balanced_accuracy"].mean()),
        "inner_mean_ECE": float(frame["ECE"].mean()),
        "inner_sd_AUROC": float(frame["AUROC"].std(ddof=0)),
        "performance_score": performance_score,
        **bio_terms,
        "selection_score": score,
    }


def _metadata_value(meta: pd.Series, key: str, default: object = pd.NA) -> object:
    return meta.get(key, default) if hasattr(meta, "get") else default


def evaluate_module_model(
    module_features: dict[str, pd.DataFrame],
    y_response: dict[str, pd.Series],
    metadata_by_cohort: dict[str, pd.DataFrame],
    endpoint: str,
    stratum: str,
    train_pool: list[str],
    holdouts: list[str],
    raw_X_by_cohort: dict[str, pd.DataFrame] | None = None,
    optimizer_config: HeuristicEcologyConfig | None = None,
    enable_optimizer: bool = True,
) -> EvaluationResult:
    metric_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    inner_rows: list[dict[str, object]] = []
    weight_rows: list[dict[str, object]] = []
    specs = _candidate_specs(optimize_word=enable_optimizer)
    spec_by_name = {str(spec["candidate"]): spec for spec in specs}
    word_specs = [spec for spec in specs if str(spec.get("kind")) == "word_ecology"]
    ranked_expression_by_cohort = _rank_expression_by_cohort(raw_X_by_cohort) if raw_X_by_cohort is not None else None
    word_matrix_cache: dict[tuple[object, ...], tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, int], HeuristicEcologyResult | None]] = {}
    word_objective_cache: dict[tuple[object, ...], dict[str, float]] = {}
    for holdout in holdouts:
        train_cohorts = [cohort for cohort in train_pool if cohort != holdout and cohort in module_features]
        if holdout not in module_features or len(train_cohorts) < 1:
            continue
        y_train = _concat(y_response, train_cohorts).astype(int)
        y_test = y_response[holdout].astype(int)
        if y_train.nunique() < 2 or y_test.nunique() < 2:
            continue
        inner_scores = [
            _score_inner_candidate(
                spec,
                train_cohorts,
                module_features,
                y_response,
                ranked_expression_by_cohort=ranked_expression_by_cohort,
                metadata_by_cohort=metadata_by_cohort,
                word_matrix_cache=word_matrix_cache,
                word_objective_cache=word_objective_cache,
                optimizer_config=optimizer_config,
            )
            for spec in specs
        ]
        for row in inner_scores:
            inner_rows.append({"endpoint": endpoint, "stratum": stratum, "holdout": holdout, **row})
        valid = [row for row in inner_scores if np.isfinite(float(row["selection_score"]))]
        if valid:
            raw_best = max(valid, key=lambda row: float(row["selection_score"]))
            best = raw_best
            baseline = next((row for row in valid if str(row["candidate"]) == "module_prior_composite"), None)
            if baseline is not None:
                best = baseline
            selected = str(best["candidate"])
        else:
            selected = "module_prior_composite"
            best = {
                "candidate": selected,
                "inner_mean_AUROC": np.nan,
                "inner_mean_AUPRC": np.nan,
                "inner_mean_balanced_accuracy": np.nan,
                "selection_score": np.nan,
            }
        specs_to_emit = [(spec_by_name[selected], True)]
        emitted_candidates: set[str] = set()
        for spec in word_specs:
            candidate = str(spec["candidate"])
            if candidate not in emitted_candidates:
                specs_to_emit.append((spec, False))
                emitted_candidates.add(candidate)

        inner_by_candidate = {str(row["candidate"]): row for row in inner_scores}
        for spec, is_selected_row in specs_to_emit:
            candidate_name = str(spec["candidate"])
            try:
                train_prob, test_prob, weights = _predict_candidate(
                    spec,
                    train_cohorts,
                    holdout,
                    module_features,
                    y_response,
                    ranked_expression_by_cohort=ranked_expression_by_cohort,
                    metadata_by_cohort=metadata_by_cohort,
                    word_matrix_cache=word_matrix_cache,
                    word_objective_cache=word_objective_cache,
                    optimizer_config=optimizer_config,
                )
            except Exception:
                continue
            threshold = select_threshold(y_train.to_numpy(dtype=int), np.asarray(train_prob, dtype=float))
            metrics = compute_binary_metrics(y_test, test_prob, threshold=threshold)
            inner_record = inner_by_candidate.get(candidate_name, best)
            model_name = _model_name_for_spec(spec, selected=is_selected_row)
            metric_rows.append(
                {
                    **metrics,
                    "endpoint": endpoint,
                    "stratum": stratum,
                    "evaluation": "LODO" if holdout in train_pool else "external_transfer",
                    "cohort": holdout,
                    "model_name": model_name,
                    "n_samples": len(y_test),
                    "n_responders": int((y_test == 1).sum()),
                    "n_nonresponders": int((y_test == 0).sum()),
                    "selected_model": selected if is_selected_row else candidate_name,
                    "threshold": threshold,
                    "inner_mean_AUROC": inner_record.get("inner_mean_AUROC", np.nan),
                    "inner_mean_AUPRC": inner_record.get("inner_mean_AUPRC", np.nan),
                    "inner_mean_balanced_accuracy": inner_record.get("inner_mean_balanced_accuracy", np.nan),
                    "selection_score": inner_record.get("selection_score", np.nan),
                    "train_cohorts": ",".join(train_cohorts),
                }
            )
            weights = weights.copy()
            if not weights.empty:
                weights.insert(0, "endpoint", endpoint)
                weights.insert(1, "stratum", stratum)
                weights.insert(2, "holdout", holdout)
                weights.insert(3, "selected_model", selected if is_selected_row else candidate_name)
                weights.insert(4, "model_name", model_name)
                weight_rows.extend(weights.to_dict("records"))
            meta = metadata_by_cohort[holdout].reindex(module_features[holdout].index)
            pred_label = np.asarray(test_prob) >= threshold
            for idx, sample_id in enumerate(module_features[holdout].index):
                m = meta.loc[sample_id] if sample_id in meta.index else pd.Series(dtype=object)
                prediction_rows.append(
                    {
                        "sample_id": _metadata_value(m, "sample_id", sample_id),
                        "patient_id": _metadata_value(m, "patient_id"),
                        "cohort": holdout,
                        "endpoint": endpoint,
                        "stratum": stratum,
                        "model_name": model_name,
                        "fold": holdout,
                        "true_response_label": int(y_test.loc[sample_id]),
                        "response_probability": float(test_prob[idx]),
                        "pred_response_label": int(pred_label[idx]),
                        "selected_model": selected if is_selected_row else candidate_name,
                        "threshold": float(threshold),
                        "response_raw": _metadata_value(m, "response_raw"),
                        "timepoint": _metadata_value(m, "timepoint"),
                        "treatment": _metadata_value(m, "treatment"),
                    }
                )
    return EvaluationResult(
        metrics=pd.DataFrame(metric_rows),
        predictions=pd.DataFrame(prediction_rows),
        inner_selection=pd.DataFrame(inner_rows),
        feature_weights=pd.DataFrame(weight_rows),
    )


def evaluate_fixed_score_models(
    fixed_scores: dict[str, dict[str, pd.Series]],
    y_response: dict[str, pd.Series],
    metadata_by_cohort: dict[str, pd.DataFrame],
    endpoint: str,
    stratum: str,
    train_pool: list[str],
    holdouts: list[str],
) -> EvaluationResult:
    metric_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    model_names = sorted({name for cohort_scores in fixed_scores.values() for name in cohort_scores})
    evaluation_specs = [(name, name, False) for name in model_names]
    evaluation_specs.extend(
        (source_name, calibrated_name, True)
        for source_name, calibrated_name in CALIBRATED_FIXED_SCORE_MODELS.items()
        if source_name in model_names
    )
    for source_model_name, model_name, use_calibration in evaluation_specs:
        for holdout in holdouts:
            train_cohorts = [
                cohort
                for cohort in train_pool
                if cohort != holdout and cohort in fixed_scores and source_model_name in fixed_scores[cohort]
            ]
            if holdout not in fixed_scores or source_model_name not in fixed_scores[holdout] or len(train_cohorts) < 1:
                continue
            y_train = _concat(y_response, train_cohorts).astype(int)
            y_test = y_response[holdout].astype(int)
            if y_train.nunique() < 2 or y_test.nunique() < 2:
                continue
            train_score = pd.concat([fixed_scores[cohort][source_model_name] for cohort in train_cohorts])
            test_score = fixed_scores[holdout][source_model_name]
            if use_calibration:
                calibrator = _fit_monotone_score_calibrator(train_score, y_train)
                train_prob = _predict_score_calibrator(calibrator, train_score)
                test_prob = _predict_score_calibrator(calibrator, test_score)
                selected_model = "fixed_score_platt" if calibrator is not None else "fixed_score_platt_fallback"
            elif source_model_name in DIRECT_PROBABILITY_MODELS:
                train_prob = pd.to_numeric(train_score, errors="coerce").clip(1e-6, 1 - 1e-6).to_numpy(dtype=float)
                test_prob = pd.to_numeric(test_score, errors="coerce").clip(1e-6, 1 - 1e-6).to_numpy(dtype=float)
                selected_model = "fixed_score"
            else:
                train_prob = sigmoid(train_score)
                test_prob = sigmoid(test_score)
                selected_model = "fixed_score"
            threshold = select_threshold(y_train.to_numpy(dtype=int), np.asarray(train_prob, dtype=float))
            metrics = compute_binary_metrics(y_test, test_prob, threshold=threshold)
            metric_rows.append(
                {
                    **metrics,
                    "endpoint": endpoint,
                    "stratum": stratum,
                    "evaluation": "LODO" if holdout in train_pool else "external_transfer",
                    "cohort": holdout,
                    "model_name": model_name,
                    "n_samples": len(y_test),
                    "n_responders": int((y_test == 1).sum()),
                    "n_nonresponders": int((y_test == 0).sum()),
                    "selected_model": selected_model,
                    "threshold": threshold,
                    "train_cohorts": ",".join(train_cohorts),
                }
            )
            meta = metadata_by_cohort[holdout].reindex(test_score.index)
            pred_label = np.asarray(test_prob) >= threshold
            for idx, sample_id in enumerate(test_score.index):
                m = meta.loc[sample_id] if sample_id in meta.index else pd.Series(dtype=object)
                prediction_rows.append(
                    {
                        "sample_id": _metadata_value(m, "sample_id", sample_id),
                        "patient_id": _metadata_value(m, "patient_id"),
                        "cohort": holdout,
                        "endpoint": endpoint,
                        "stratum": stratum,
                        "model_name": model_name,
                        "fold": holdout,
                        "true_response_label": int(y_test.loc[sample_id]),
                        "response_probability": float(test_prob[idx]),
                        "pred_response_label": int(pred_label[idx]),
                        "selected_model": selected_model,
                        "threshold": float(threshold),
                        "response_raw": _metadata_value(m, "response_raw"),
                        "timepoint": _metadata_value(m, "timepoint"),
                        "treatment": _metadata_value(m, "treatment"),
                    }
                )
    return EvaluationResult(
        metrics=pd.DataFrame(metric_rows),
        predictions=pd.DataFrame(prediction_rows),
        inner_selection=pd.DataFrame(),
        feature_weights=pd.DataFrame(),
    )
