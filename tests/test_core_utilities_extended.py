import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from econiche.baselines import score_baselines, signature_score
from econiche.expression import (
    choose_expression_file,
    clean_expression_matrix,
    collapse_duplicate_genes,
    load_entrez_symbol_map,
    normalize_gene_index,
    orient_expression,
    read_table_matrix,
)
from econiche.geo import (
    harmonized_geo_metadata,
    infer_patient_id,
    infer_response,
    infer_timepoint,
    parse_series_matrix_metadata,
)
from econiche.module import EcoNicheModule
from econiche.networks import lr_coherence, network_coherence, pathway_coherence
from econiche.plotting import plot_metric_bar, write_pending_figure
from econiche.statistics import benjamini_hochberg, make_lodo_folds, paired_bootstrap_delta
from econiche.survival import cox_placeholder
from econiche.utils import write_json


def test_baseline_scores_cover_available_and_missing_genes():
    X = pd.DataFrame({"CD274": [1.0, 3.0], "GZMA": [2.0, 4.0]}, index=["S1", "S2"])
    meta = pd.DataFrame({"patient_id": ["P1", "P2"], "cohort": ["C1", "C1"], "label": [0, 1]}, index=X.index)
    assert signature_score(X, ["CD274"]).notna().all()
    assert signature_score(X, ["MISSING"]).isna().all()
    scored = score_baselines(X, meta, {"available": ["CD274"], "missing": ["NOPE"]})
    assert set(scored["status"]) == {"available", "unavailable_with_reason"}
    assert scored.loc[scored["model_name"] == "available", "n_genes_available"].iloc[0] == 1


def test_expression_helpers_read_clean_and_choose_files(tmp_path: Path):
    gene_info = tmp_path / "gene_info.gz"
    with gzip.open(gene_info, "wt", encoding="utf-8") as handle:
        handle.write("#tax\tGeneID\tSymbol\n9606\t1\tA1BG\n")
    assert load_entrez_symbol_map(gene_info) == {"1": "A1BG"}
    assert normalize_gene_index(['"1.5"', "GENE_12"], {"1": "A1BG"}) == ["A1BG", "GENE"]

    duplicated = pd.DataFrame({"S1": ["1", "3"], "S2": ["2", "4"]}, index=["A", "A"])
    collapsed = collapse_duplicate_genes(duplicated)
    assert collapsed.loc["A", "S1"] == 2.0

    matrix_path = tmp_path / "matrix.tsv"
    matrix_path.write_text("gene\tS1\tS2\nA\t1\t2\nB\t3\t4\n", encoding="utf-8")
    raw = read_table_matrix(matrix_path)
    assert raw.shape == (2, 2)
    cleaned = clean_expression_matrix(raw)
    assert list(cleaned.index) == ["S1", "S2"]

    sample_rows = pd.DataFrame({"A": [1, 2], "B": [3, 4], "C": [5, 6]}, index=["Pt1", "Pt2"])
    cleaned_sample_rows = clean_expression_matrix(sample_rows)
    assert "A" in cleaned_sample_rows.columns

    oriented = orient_expression(pd.DataFrame([[1, 2, 3]], index=["S1"], columns=["A", "B", "C"]), {"A", "B", "C"})
    assert oriented.shape == (3, 1)

    suppl = tmp_path / "GSE91061" / "suppl"
    suppl.mkdir(parents=True)
    preferred = suppl / "sample_fpkm.tsv"
    fallback = suppl / "other.tsv"
    fallback.write_text("x", encoding="utf-8")
    preferred.write_text("x", encoding="utf-8")
    unknown_suppl = tmp_path / "UNKNOWN" / "suppl"
    unknown_suppl.mkdir(parents=True)
    unknown_fallback = unknown_suppl / "other.tsv"
    unknown_fallback.write_text("x", encoding="utf-8")
    assert choose_expression_file("GSE91061", tmp_path) == preferred
    assert choose_expression_file("UNKNOWN", tmp_path) == unknown_fallback


def test_geo_metadata_inference_and_series_matrix(tmp_path: Path):
    matrix = tmp_path / "GSE_test_series_matrix.txt"
    matrix.write_text(
        "\n".join(
            [
                '!Sample_geo_accession\t"GSM1"\t"GSM2"',
                '!Sample_title\t"Pt1 complete response pre"\t"Pt2 progressive disease on treatment"',
                '!Sample_characteristics_ch1\t"patient: P1"\t"patient: P2"',
                '!Sample_characteristics_ch1\t"response: CR"\t"response: PD"',
                "!series_matrix_table_begin",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    parsed = parse_series_matrix_metadata(matrix)
    assert parsed.shape[0] == 2
    assert infer_patient_id(parsed.iloc[0]) == "P1"
    assert infer_response(parsed.iloc[0]) == "CR"
    assert infer_timepoint(parsed.iloc[0]) == "pretreatment"
    harmonized = harmonized_geo_metadata(parsed, "GSETEST", platform="RNA-seq")
    assert set(["sample_id", "patient_id", "response_raw", "timepoint"]).issubset(harmonized.columns)

    assert infer_patient_id(pd.Series({"title": "Patient 7 sample"})) == "Patient7"
    assert infer_response(pd.Series({"title": "no response"})) == "NR"
    assert infer_timepoint(pd.Series({"title": "progression sample"})) == "progression"


def test_network_statistics_survival_plotting_and_utils(tmp_path: Path):
    module = EcoNicheModule({"s1": {"A", "B"}, "s2": {"C"}})
    assert pathway_coherence(module, {"p": {"A", "B"}}) > 0
    assert network_coherence(module, {("A", "B")}) == 1.0
    assert lr_coherence(module, {("A", "C")}, [("s1", "s2")]) == 0.5

    meta = pd.DataFrame(
        {
            "cohort": ["A", "A", "B", "B"],
            "label": [0, 1, 0, 1],
            "patient_id": ["P1", "P2", "P3", "P4"],
        }
    )
    folds = make_lodo_folds(meta)
    assert set(folds) == {"A", "B"}
    assert np.all(benjamini_hochberg([0.01, 0.2]) <= 1)
    delta = paired_bootstrap_delta([0, 0, 1, 1], [0.1, 0.2, 0.8, 0.9], [0.4, 0.3, 0.6, 0.7], n_bootstrap=20)
    assert {"mean_delta", "p_value"}.issubset(delta)
    empty_delta = paired_bootstrap_delta([1, 1, 1], [0.1, 0.2, 0.3], [0.2, 0.2, 0.2], n_bootstrap=5)
    assert np.isnan(empty_delta["mean_delta"])

    pending = cox_placeholder(pd.DataFrame({"EcoNicheScore": [1.0]}))
    assert pending.loc[0, "status"] == "RESULT_PENDING"

    out_json = tmp_path / "nested" / "data.json"
    write_json({"b": 1}, out_json)
    assert out_json.exists()
    write_pending_figure(tmp_path / "pending.svg", "Title")
    plot_metric_bar(pd.DataFrame({"cohort": ["A"], "AUROC": [0.75]}), tmp_path / "metric.svg")
    assert (tmp_path / "pending.svg").exists()
    assert (tmp_path / "metric.svg").exists()
