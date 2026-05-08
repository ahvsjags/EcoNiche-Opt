from pathlib import Path

import pandas as pd
import yaml

from econiche_opt.cli import main
from econiche_opt.data.download_plan import build_download_plan, write_download_plan
from econiche_opt.data.registry import audit_accession, load_registry, validate_registry
from econiche_opt.panel.compress import compress_module_panel
from econiche_opt.perturbation.reversal import aggregate_reversal_candidates, build_reversal_signature, parse_lincs_term, score_reversal
from econiche_opt.reporting.citations import validate_source_registry
from econiche_opt.reporting.reproducibility import generate_reproducibility_report
from econiche_opt.reporting.safety_checks import validate_manuscript_safety
from econiche_opt.utils.diagnostics import EcoNicheDiagnosticError, require_file
from econiche_opt.utils.random import set_global_seed
from econiche_opt.utils.resources import assert_matrix_within_memory, estimate_matrix_megabytes
from econiche_opt.validation.project import validate_project
from econiche_opt.validation.schema import validate_result_schema


def test_registry_download_plan_and_source_validation(tmp_path: Path):
    registry = {
        "cohorts": [
            {
                "accession": "GSE1",
                "layer": "A",
                "cancer_type": "melanoma",
                "therapy": "anti-PD1",
                "platform": "RNA-seq",
                "timepoints": ["baseline"],
                "endpoint": ["response"],
                "role": "demo",
                "access": "public",
                "download_script": "download.py",
                "preprocessing_script": "pre.py",
                "uses": ["demo"],
            },
            {
                "accession": "CTRL",
                "layer": "A",
                "cancer_type": "melanoma",
                "therapy": "anti-PD1",
                "platform": "RNA-seq",
                "timepoints": ["baseline"],
                "endpoint": ["response"],
                "role": "controlled",
                "access": "controlled",
                "download_script": "ACCESS_RESTRICTED",
                "preprocessing_script": "pre.py",
                "uses": ["external"],
            },
        ]
    }
    path = tmp_path / "registry.yml"
    path.write_text(yaml.safe_dump(registry), encoding="utf-8")
    loaded = load_registry(path)
    assert validate_registry(loaded, require_minimum=False)["is_valid"].all()
    assert set(audit_accession(loaded)["download_action"]) == {"download_or_metadata_extract", "emit_access_instructions"}
    plan = build_download_plan(path)
    assert set(plan["planned_status"]) == {"READY", "ACCESS_RESTRICTED"}
    out = write_download_plan(path, tmp_path / "plan.tsv")
    assert len(out) == 2

    source_path = tmp_path / "sources.yml"
    source_path.write_text(
        yaml.safe_dump({"sources": [{"id": "x", "type": "db", "title": "T", "source": "S", "access_status": "public"}]}),
        encoding="utf-8",
    )
    assert validate_source_registry(source_path)["is_valid"].all()


def test_reporting_validation_cli_and_utility_interfaces(tmp_path: Path):
    manuscript = tmp_path / "paper.md"
    manuscript.write_text("EcoNiche-Opt reports descriptive benchmark results.\n", encoding="utf-8")
    assert validate_manuscript_safety(manuscript).empty
    report = generate_reproducibility_report(tmp_path / "repro.md")
    assert report.exists()

    results = tmp_path / "results"
    results.mkdir()
    pd.DataFrame({"sample_id": ["S1"], "patient_id": ["P1"], "cohort": ["C"], "pred_prob": [0.2]}).to_csv(
        results / "lodo_predictions.tsv", sep="\t", index=False
    )
    pd.DataFrame({"cohort": ["C"], "model_name": ["M"]}).to_csv(results / "lodo_metrics.tsv", sep="\t", index=False)
    pd.DataFrame({"state": ["s"], "gene": ["G"], "direction": [1]}).to_csv(
        results / "econiche_module.tsv", sep="\t", index=False
    )
    assert validate_result_schema(results)["is_valid"].all()

    root = tmp_path / "project"
    (root / "config").mkdir(parents=True)
    (root / "docs" / "reproducibility").mkdir(parents=True)
    (root / "results" / "demo").mkdir(parents=True)
    for name in ["README.md", "AGENTS.md"]:
        (root / name).write_text("x", encoding="utf-8")
    (root / "docs" / "reproducibility" / "no_fabrication_policy.md").write_text("x", encoding="utf-8")
    (root / "docs" / "goal_status.yml").write_text(yaml.safe_dump({"goals": {}}), encoding="utf-8")
    (root / "config" / "data_registry.yml").write_text(yaml.safe_dump({"cohorts": []}), encoding="utf-8")
    project_report = validate_project(root)
    assert "goal_status" in set(project_report["check"])

    assert main(["validate-results", "--results-dir", str(results)]) == 0

    set_global_seed(123)
    assert estimate_matrix_megabytes(10, 10) > 0
    assert_matrix_within_memory(10, 10, max_megabytes=1)
    try:
        require_file(tmp_path / "missing.txt")
    except EcoNicheDiagnosticError:
        pass
    else:
        raise AssertionError("require_file should raise for missing files")


def test_panel_and_perturbation_interfaces():
    module = pd.DataFrame(
        {
            "state": ["s1", "s1", "s2"],
            "gene": ["A", "B", "C"],
            "selection_frequency": [0.2, 0.9, 0.5],
        }
    )
    panel = compress_module_panel(module, max_genes=2)
    assert panel["gene"].tolist() == ["B", "C"]
    perturbation = pd.DataFrame({"drug": ["D1", "D1", "D2"], "gene": ["A", "C", "X"], "effect": [-1.0, 0.5, -3.0]})
    scored = score_reversal(module, perturbation)
    assert "reversal_score" in scored.columns

    signature = build_reversal_signature(
        module.assign(direction=[1, -1, 1], coefficient=[2.0, 1.0, -1.0]),
        top_genes=2,
    )
    assert signature["resistance_up"] == ["A"]
    assert set(signature["resistance_down"]) == {"B", "C"}
    parsed = parse_lincs_term("LJP008 MCF7 24H-UNC0638-0.12")
    assert parsed["perturbation_name"] == "UNC0638"
    summary = aggregate_reversal_candidates(
        pd.DataFrame(
            {
                "perturbation_id": ["LJP008"],
                "perturbation_name": ["UNC0638"],
                "library": ["LINCS_L1000_Chem_Pert_down"],
                "query_direction": ["downregulate_resistance_up_genes"],
                "overlapping_genes": ["A;B"],
                "reversal_score": [10.0],
            }
        )
    )
    assert summary.iloc[0]["perturbation_name"] == "UNC0638"
