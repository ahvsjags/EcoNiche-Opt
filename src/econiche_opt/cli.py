from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from econiche_opt.analysis.precomputed_scores import import_precomputed_scores
from econiche_opt.data.download_plan import write_download_plan
from econiche_opt.data.registry import audit_accession, load_registry, validate_registry
from econiche_opt.deploy import fit_package_model, score_package_model
from econiche_opt.pipeline.real import RealPipelineConfig, assert_real_pipeline_ok, run_real_pipeline
from econiche_opt.reporting.citations import validate_source_registry
from econiche_opt.reporting.reproducibility import generate_reproducibility_report
from econiche_opt.reporting.safety_checks import validate_manuscript_safety
from econiche_opt.validation.goals import assert_goal_status_valid, validate_goal_status
from econiche_opt.validation.project import assert_project_valid, validate_project
from econiche_opt.validation.schema import assert_result_schema_valid, validate_result_schema

ROOT = Path(__file__).resolve().parents[2]


def _run_script(script: str, *args: str) -> int:
    command = [sys.executable, str(ROOT / script), *args]
    return subprocess.run(command, cwd=ROOT, check=True).returncode


def _write_report(frame, out: str | None, default: str) -> None:
    if out:
        out_path = ROOT / out if not Path(out).is_absolute() else Path(out)
    else:
        out_path = ROOT / default
    out_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out_path, sep="\t", index=False)
    print(frame.to_string(index=False))
    print(f"Wrote {out_path}")


def cmd_validate_registry(args: argparse.Namespace) -> int:
    report = validate_registry(load_registry(args.registry), require_minimum=not args.skip_minimum)
    _write_report(report, args.out, "results/audit/registry_validation.tsv")
    if not report["is_valid"].all():
        raise SystemExit(1)
    return 0


def cmd_audit_dataset(args: argparse.Namespace) -> int:
    report = audit_accession(load_registry(args.registry))
    _write_report(report, args.out, "results/audit/data_access_audit.tsv")
    return 0


def cmd_download_data(args: argparse.Namespace) -> int:
    out = args.out or "results/audit/download_plan.tsv"
    plan = write_download_plan(args.registry, ROOT / out)
    print(plan.to_string(index=False))
    if args.dry_run:
        print("Dry run only; no data downloaded.")
        return 0
    return _run_script("scripts/download/download_geo.py", "--registry", args.registry, "--metadata-only")


def cmd_preprocess_data(args: argparse.Namespace) -> int:
    if args.demo:
        return _run_script("scripts/make_demo_data.py")
    for script in [
        "scripts/preprocess/build_metadata.py",
        "scripts/preprocess/harmonize_labels.py",
        "scripts/preprocess/deduplicate_patients.py",
        "scripts/preprocess/preprocess_bulk.py",
    ]:
        _run_script(script)
    return 0


def cmd_make_demo(args: argparse.Namespace) -> int:
    _run_script("scripts/make_demo_data.py")
    _run_script("scripts/model/run_econiche.py", "--config", "config/model_config.yml", "--demo")
    return 0


def cmd_train_econiche(args: argparse.Namespace) -> int:
    script_args = ["--config", args.config]
    if args.demo:
        script_args.append("--demo")
    if args.out:
        script_args.extend(["--out", args.out])
    return _run_script("scripts/model/run_econiche.py", *script_args)


def cmd_run_response_composite(args: argparse.Namespace) -> int:
    script_args = [
        "--processed-dir",
        args.processed_dir,
        "--baseline-predictions",
        args.baseline_predictions,
        "--ml-baseline-predictions",
        args.ml_baseline_predictions,
        "--current-predictions",
        args.current_predictions,
        "--out",
        args.out,
        "--preferred-candidate",
        args.preferred_candidate,
        "--preferred-tolerance",
        str(args.preferred_tolerance),
    ]
    if args.include_demo:
        script_args.append("--include-demo")
    return _run_script("scripts/model/run_response_composite.py", *script_args)


def cmd_fit_package_model(args: argparse.Namespace) -> int:
    optimizer_kwargs = {}
    if args.mode == "heuristic_ecology":
        optimizer_kwargs = {
            "population_size": args.population_size,
            "generations": args.generations,
            "n_jobs": args.n_jobs,
            "use_gpu": args.use_gpu,
            "random_state": args.random_state,
        }
    model = fit_package_model(
        expression_path=args.expression,
        labels_path=args.labels,
        model_path=args.model_out,
        out_dir=args.out_dir,
        mode=args.mode,
        transpose=args.transpose,
        label_column=args.label_column,
        cohort_column=args.cohort_column,
        calibration=args.calibration,
        random_state=args.random_state,
        optimizer_kwargs=optimizer_kwargs,
    )
    print(f"Wrote model: {args.model_out}")
    print(f"Wrote module/edge/metadata artifacts: {args.out_dir}")
    print(model.package_metadata())
    return 0


def cmd_score_package_model(args: argparse.Namespace) -> int:
    scores = score_package_model(
        model_path=args.model,
        expression_path=args.expression,
        out_path=args.out,
        transpose=args.transpose,
        coverage_path=args.coverage_out,
    )
    print(scores.head().to_string(index=False))
    print(f"Wrote scores: {args.out}")
    if args.coverage_out:
        print(f"Wrote feature coverage: {args.coverage_out}")
    return 0


def cmd_score_locked_validation(args: argparse.Namespace) -> int:
    script_args = [
        "--package-dir",
        args.package_dir,
        "--expression",
        args.expression,
        "--sample-manifest",
        args.sample_manifest,
        "--out-dir",
        args.out_dir,
    ]
    if args.clinical_annotation:
        script_args.extend(["--clinical-annotation", args.clinical_annotation])
    if args.transpose:
        script_args.append("--transpose")
    return _run_script("scripts/validation/score_locked_validation_cohort.py", *script_args)


def cmd_run_benchmark(args: argparse.Namespace) -> int:
    _run_script("scripts/benchmark/run_lodo.py")
    _run_script("scripts/benchmark/bootstrap_compare.py")
    _run_script("scripts/benchmark/calibration.py")
    _run_script("scripts/benchmark/decision_curve.py")
    return 0


def cmd_run_pancancer(args: argparse.Namespace) -> int:
    return _run_script("scripts/benchmark/run_pan_cancer_transfer.py")


def cmd_run_survival(args: argparse.Namespace) -> int:
    return _run_script("scripts/benchmark/survival_analysis.py")


def cmd_run_single_cell(args: argparse.Namespace) -> int:
    return _run_script("scripts/single_cell/map_modules_scrna.py")


def cmd_run_deconvolution(args: argparse.Namespace) -> int:
    script_args = ["--out-dir", args.out_dir, "--min-markers", str(args.min_markers)]
    if args.input_dir:
        script_args.extend(["--input-dir", args.input_dir])
    if args.include_demo:
        script_args.append("--include-demo")
    return _run_script("scripts/preprocess/run_deconvolution.py", *script_args)


def cmd_download_xena(args: argparse.Namespace) -> int:
    script_args = ["--out-dir", args.out_dir, "--manifest", args.manifest, "--max-download-mb", str(args.max_download_mb)]
    if args.download:
        script_args.append("--download")
    if args.skip_gdc_expression:
        script_args.append("--skip-gdc-expression")
    if args.strict:
        script_args.append("--strict")
    return _run_script("scripts/download/download_xena.py", *script_args)


def cmd_run_perturbation(args: argparse.Namespace) -> int:
    _run_script("scripts/perturbation/lincs_reversal.py")
    _run_script("scripts/perturbation/depmap_prioritize.py")
    _run_script("scripts/perturbation/dgidb_lookup.py")
    if args.demo:
        _run_script("scripts/analysis/run_perturbation_prioritization.py", "--demo")
    else:
        _run_script("scripts/analysis/run_perturbation_prioritization.py")
    return 0


def cmd_make_figures(args: argparse.Namespace) -> int:
    scripts = [
        "scripts/figures/make_fig1_overview.py",
        "scripts/figures/make_fig2_benchmark.py",
        "scripts/figures/make_fig3_module_network.py",
        "scripts/figures/make_fig4_single_cell.py",
        "scripts/figures/make_fig5_survival.py",
        "scripts/figures/make_fig6_perturbation.py",
        "scripts/figures/figure1_overview.py",
        "scripts/figures/figure2_model.py",
        "scripts/figures/figure3_benchmark.py",
        "scripts/figures/figure4_pancancer.py",
        "scripts/figures/figure5_single_cell.py",
        "scripts/figures/figure6_perturbation.py",
    ]
    for script in scripts:
        if args.demo and Path(script).name.startswith("figure"):
            _run_script(script, "--demo")
        else:
            _run_script(script)
    return 0


def cmd_validate_goals(args: argparse.Namespace) -> int:
    report = validate_goal_status(args.goal_file, demo_mode=not args.real_mode)
    print(report.to_string(index=False))
    assert_goal_status_valid(args.goal_file, demo_mode=not args.real_mode)
    return 0


def cmd_validate_results(args: argparse.Namespace) -> int:
    results_dir = args.results_dir or ("results/demo" if args.demo else "results/real")
    report = validate_result_schema(ROOT / results_dir, demo=args.demo)
    print(report.to_string(index=False))
    assert_result_schema_valid(ROOT / results_dir, demo=args.demo)
    return 0


def cmd_validate_sources(args: argparse.Namespace) -> int:
    report = validate_source_registry(args.source_registry)
    _write_report(report, args.out, "results/audit/source_registry_validation.tsv")
    if not report["is_valid"].all():
        raise SystemExit(1)
    return 0


def cmd_import_precomputed_scores(args: argparse.Namespace) -> int:
    frame = import_precomputed_scores(args.input, args.out)
    print(frame.head().to_string(index=False))
    return 0


def cmd_check_claims(args: argparse.Namespace) -> int:
    report = validate_manuscript_safety(args.manuscript, args.evidence)
    _write_report(report, args.out, "results/audit/manuscript_claim_safety.tsv")
    if not report.empty:
        raise SystemExit(1)
    return 0


def cmd_make_reproducibility_report(args: argparse.Namespace) -> int:
    out = args.out or "paper/reproducibility_report.md"
    path = generate_reproducibility_report(ROOT / out, mode=args.mode)
    print(f"Wrote {path}")
    return 0


def cmd_validate_project(args: argparse.Namespace) -> int:
    report = validate_project(ROOT, mode=args.mode)
    print(report.to_string(index=False))
    assert_project_valid(ROOT, mode=args.mode)
    return 0


def cmd_run_real_pipeline(args: argparse.Namespace) -> int:
    cfg = RealPipelineConfig(
        root=ROOT,
        registry=ROOT / args.registry,
        results_dir=ROOT / args.results_dir,
        out_dir=ROOT / args.out_dir,
        execute_download=args.execute_download,
        execute_preprocess=args.execute_preprocess,
        execute_training=args.execute_training,
        execute_secondary=not args.skip_secondary,
        strict_existing_results=not args.allow_missing_results,
    )
    manifest = run_real_pipeline(cfg)
    print(manifest[["stage", "status", "reason"]].to_string(index=False))
    assert_real_pipeline_ok(manifest)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="econiche-opt")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("validate-registry")
    p.add_argument("--registry", default="config/data_registry.yml")
    p.add_argument("--out")
    p.add_argument("--skip-minimum", action="store_true")
    p.set_defaults(func=cmd_validate_registry)

    p = sub.add_parser("audit-dataset")
    p.add_argument("--registry", default="config/data_registry.yml")
    p.add_argument("--out")
    p.set_defaults(func=cmd_audit_dataset)

    p = sub.add_parser("download-data")
    p.add_argument("--registry", default="config/data_registry.yml")
    p.add_argument("--out")
    p.add_argument("--dry-run", action="store_true")
    p.set_defaults(func=cmd_download_data)

    p = sub.add_parser("preprocess-data")
    p.add_argument("--demo", action="store_true")
    p.set_defaults(func=cmd_preprocess_data)

    p = sub.add_parser("check-leakage")
    p.add_argument("--metadata", default="")
    p.set_defaults(func=lambda args: _run_script("scripts/preprocess/deduplicate_patients.py"))

    p = sub.add_parser("make-demo")
    p.set_defaults(func=cmd_make_demo)

    p = sub.add_parser("train-econiche")
    p.add_argument("--config", default="config/model_config.yml")
    p.add_argument("--demo", action="store_true")
    p.add_argument("--out")
    p.set_defaults(func=cmd_train_econiche)

    p = sub.add_parser("run-response-composite")
    p.add_argument("--processed-dir", default="data/processed/bulk")
    p.add_argument("--baseline-predictions", default="results/real/baseline_predictions.tsv")
    p.add_argument("--ml-baseline-predictions", default="results/real/ml_baseline_predictions.tsv")
    p.add_argument("--current-predictions", default="results/real/lodo_predictions.tsv")
    p.add_argument("--out", default="results/real_optimized")
    p.add_argument("--include-demo", action="store_true")
    p.add_argument("--preferred-candidate", default="ifn_core_pdcd1lg2_weighted")
    p.add_argument("--preferred-tolerance", type=float, default=0.02)
    p.set_defaults(func=cmd_run_response_composite)

    p = sub.add_parser("fit-package-model")
    p.add_argument("--expression", required=True, help="TSV/CSV expression matrix; samples rows, genes columns by default")
    p.add_argument("--labels", required=True, help="TSV/CSV with sample_id, response_label, and optional cohort")
    p.add_argument("--model-out", required=True, help="Output .joblib model path")
    p.add_argument("--out-dir", default="results/package_model", help="Output directory for module, edge, and metadata artifacts")
    p.add_argument("--mode", choices=["word_full_graph", "heuristic_ecology"], default="word_full_graph")
    p.add_argument("--transpose", action="store_true", help="Use when input expression has genes as rows and samples as columns")
    p.add_argument("--label-column", default="response_label")
    p.add_argument("--cohort-column", default="cohort")
    p.add_argument("--calibration", choices=["isotonic"], default=None)
    p.add_argument("--random-state", type=int, default=42)
    p.add_argument("--population-size", type=int, default=14)
    p.add_argument("--generations", type=int, default=8)
    p.add_argument("--n-jobs", type=int, default=1)
    p.add_argument("--use-gpu", action="store_true")
    p.set_defaults(func=cmd_fit_package_model)

    p = sub.add_parser("score-package-model")
    p.add_argument("--model", required=True, help="Model .joblib written by fit-package-model or EcoNicheOptClassifier.save")
    p.add_argument("--expression", required=True, help="TSV/CSV expression matrix to score")
    p.add_argument("--out", required=True, help="Output score TSV")
    p.add_argument("--transpose", action="store_true")
    p.add_argument("--coverage-out", help="Optional feature coverage TSV")
    p.set_defaults(func=cmd_score_package_model)

    p = sub.add_parser("score-locked-validation")
    p.add_argument("--package-dir", default="deliverables/prospective_validation", help="Frozen validation package directory")
    p.add_argument("--expression", required=True, help="TSV/CSV expression matrix; samples rows, genes columns by default")
    p.add_argument("--sample-manifest", required=True, help="Assay/sample manifest matching the prospective validation template")
    p.add_argument("--clinical-annotation", help="Optional clinical annotation table with endpoint labels for metric reporting")
    p.add_argument("--out-dir", required=True, help="Output directory for locked validation scores and audits")
    p.add_argument("--transpose", action="store_true", help="Use when expression has genes as rows and samples as columns")
    p.set_defaults(func=cmd_score_locked_validation)

    p = sub.add_parser("run-benchmark")
    p.add_argument("--demo", action="store_true")
    p.set_defaults(func=cmd_run_benchmark)

    p = sub.add_parser("run-real-pipeline")
    p.add_argument("--registry", default="config/data_registry.yml")
    p.add_argument("--results-dir", default="results/real")
    p.add_argument("--out-dir", default="results/real_pipeline")
    p.add_argument("--execute-download", action="store_true")
    p.add_argument("--execute-preprocess", action="store_true")
    p.add_argument("--execute-training", action="store_true")
    p.add_argument("--skip-secondary", action="store_true")
    p.add_argument("--allow-missing-results", action="store_true")
    p.set_defaults(func=cmd_run_real_pipeline)

    for name, script in [
        ("run-ablation", "scripts/analysis/run_ablation.py"),
        ("run-sensitivity", "scripts/analysis/run_sensitivity.py"),
        ("run-stratification", "scripts/analysis/run_stratification.py"),
        ("run-public-smoke-tests", "scripts/analysis/run_public_smoke_tests.py"),
        ("make-tables", "scripts/reporting/make_tables.py"),
        ("make-manuscript", "scripts/reporting/make_manuscript.py"),
    ]:
        p = sub.add_parser(name)
        p.add_argument("--demo", action="store_true")
        p.add_argument("--out")
        p.set_defaults(func=lambda args, script=script: _run_script(script, *(["--demo"] if args.demo else []), *(["--out", args.out] if args.out else [])))

    p = sub.add_parser("run-pancancer")
    p.set_defaults(func=cmd_run_pancancer)

    p = sub.add_parser("run-survival")
    p.set_defaults(func=cmd_run_survival)

    p = sub.add_parser("run-single-cell")
    p.set_defaults(func=cmd_run_single_cell)

    p = sub.add_parser("run-deconvolution")
    p.add_argument("--input-dir", default="data/processed/bulk")
    p.add_argument("--out-dir", default="results/real")
    p.add_argument("--include-demo", action="store_true")
    p.add_argument("--min-markers", type=int, default=2)
    p.set_defaults(func=cmd_run_deconvolution)

    p = sub.add_parser("download-xena")
    p.add_argument("--out-dir", default="data/raw/TCGA_SKCM_Xena")
    p.add_argument("--manifest", default="results/real/xena_tcga_manifest.tsv")
    p.add_argument("--download", action="store_true")
    p.add_argument("--skip-gdc-expression", action="store_true")
    p.add_argument("--max-download-mb", type=float, default=200.0)
    p.add_argument("--strict", action="store_true")
    p.set_defaults(func=cmd_download_xena)

    p = sub.add_parser("run-perturbation")
    p.add_argument("--demo", action="store_true")
    p.set_defaults(func=cmd_run_perturbation)

    p = sub.add_parser("make-figures")
    p.add_argument("--demo", action="store_true")
    p.set_defaults(func=cmd_make_figures)

    p = sub.add_parser("make-reproducibility-report")
    p.add_argument("--mode", default="demo")
    p.add_argument("--out")
    p.set_defaults(func=cmd_make_reproducibility_report)

    p = sub.add_parser("validate-goals")
    p.add_argument("--goal-file", default="docs/goal_status.yml")
    p.add_argument("--real-mode", action="store_true")
    p.set_defaults(func=cmd_validate_goals)

    p = sub.add_parser("validate-results")
    p.add_argument("--results-dir")
    p.add_argument("--demo", action="store_true")
    p.set_defaults(func=cmd_validate_results)

    p = sub.add_parser("validate-sources")
    p.add_argument("--source-registry", default="config/source_registry.yml")
    p.add_argument("--out")
    p.set_defaults(func=cmd_validate_sources)

    p = sub.add_parser("import-precomputed-scores")
    p.add_argument("--input", required=True)
    p.add_argument("--out")
    p.set_defaults(func=cmd_import_precomputed_scores)

    p = sub.add_parser("check-claims")
    p.add_argument("--manuscript", default="paper/manuscript.md")
    p.add_argument("--evidence")
    p.add_argument("--out")
    p.set_defaults(func=cmd_check_claims)

    p = sub.add_parser("validate-project")
    p.add_argument("--mode", choices=["demo", "real"], default="demo")
    p.set_defaults(func=cmd_validate_project)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
