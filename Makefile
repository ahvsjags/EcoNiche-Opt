PYTHON ?= python

.PHONY: install setup demo test lint registry-audit registry_audit download-geo download_geo download-xena download_xena download-xena-real download_xena_real download-dry-run download_dry_run preprocess priors deconvolution train response-composite response_composite baselines benchmark real-pipeline real_pipeline benchmark-demo benchmark_demo pancancer survival single-cell single_cell perturbation figures paper report-demo report_demo clean-demo clean_demo validate all

install: setup

setup:
	$(PYTHON) -m pip install -r requirements.txt

demo:
	$(PYTHON) -m econiche_opt.cli make-demo

test:
	$(PYTHON) -m pytest -q

lint:
	$(PYTHON) -m compileall -q src scripts tests

registry-audit:
	$(PYTHON) -m econiche_opt.cli audit-dataset --registry config/data_registry.yml --out results/audit/data_access_audit.tsv

registry_audit: registry-audit

download-geo:
	$(PYTHON) scripts/download/download_geo.py --registry config/data_registry.yml --metadata-only

download_geo: download-geo

download-xena:
	$(PYTHON) scripts/download/download_xena.py

download_xena: download-xena

download-xena-real:
	$(PYTHON) scripts/download/download_xena.py --download --strict --max-download-mb 200

download_xena_real: download-xena-real

download-dry-run:
	$(PYTHON) scripts/download/download_geo.py --registry config/data_registry.yml --dry-run

download_dry_run: download-dry-run

preprocess:
	$(PYTHON) scripts/preprocess/build_metadata.py
	$(PYTHON) scripts/preprocess/harmonize_labels.py
	$(PYTHON) scripts/preprocess/deduplicate_patients.py
	$(PYTHON) scripts/preprocess/preprocess_bulk.py

priors:
	$(PYTHON) scripts/preprocess/build_gene_universe.py
	$(PYTHON) scripts/preprocess/build_priors.py
	$(PYTHON) scripts/preprocess/build_network.py
	$(PYTHON) scripts/preprocess/build_lr_edges.py

deconvolution:
	$(PYTHON) scripts/preprocess/run_deconvolution.py --out-dir results/real

train:
	$(PYTHON) scripts/model/run_econiche.py --config config/model_config.yml

response-composite:
	$(PYTHON) scripts/model/run_response_composite.py

response_composite: response-composite

baselines:
	$(PYTHON) scripts/baselines/run_baselines.py --config config/baselines.yml
	$(PYTHON) scripts/baselines/train_ml_baselines.py

benchmark:
	$(PYTHON) scripts/benchmark/run_lodo.py
	$(PYTHON) scripts/benchmark/bootstrap_compare.py
	$(PYTHON) scripts/benchmark/calibration.py
	$(PYTHON) scripts/benchmark/decision_curve.py

real-pipeline:
	$(PYTHON) -m econiche_opt.cli run-real-pipeline

real_pipeline: real-pipeline

benchmark-demo:
	$(PYTHON) scripts/benchmark/run_lodo.py --demo --out results/demo
	$(PYTHON) scripts/analysis/run_ablation.py --demo
	$(PYTHON) scripts/analysis/run_sensitivity.py --demo
	$(PYTHON) scripts/analysis/run_panel_compression.py --demo

benchmark_demo: benchmark-demo

pancancer:
	$(PYTHON) scripts/benchmark/run_pan_cancer_transfer.py

survival:
	$(PYTHON) scripts/benchmark/survival_analysis.py

single-cell:
	Rscript scripts/single_cell/preprocess_scrna.R
	Rscript scripts/single_cell/map_modules_scrna.R

single_cell: single-cell

perturbation:
	$(PYTHON) scripts/perturbation/lincs_reversal.py
	$(PYTHON) scripts/perturbation/depmap_prioritize.py
	$(PYTHON) scripts/perturbation/dgidb_lookup.py
	$(PYTHON) scripts/analysis/run_perturbation_prioritization.py

figures:
	$(PYTHON) scripts/figures/make_fig1_overview.py
	$(PYTHON) scripts/figures/make_fig2_benchmark.py
	$(PYTHON) scripts/figures/make_fig3_module_network.py
	$(PYTHON) scripts/figures/make_fig4_single_cell.py
	$(PYTHON) scripts/figures/make_fig5_survival.py
	$(PYTHON) scripts/figures/make_fig6_perturbation.py
	$(PYTHON) scripts/figures/figure1_overview.py --demo
	$(PYTHON) scripts/figures/figure2_model.py --demo
	$(PYTHON) scripts/figures/figure3_benchmark.py --demo
	$(PYTHON) scripts/figures/figure4_pancancer.py --demo
	$(PYTHON) scripts/figures/figure5_single_cell.py --demo
	$(PYTHON) scripts/figures/figure6_perturbation.py --demo

paper:
	$(PYTHON) scripts/paper/generate_methods_text.py
	$(PYTHON) scripts/paper/generate_result_summaries.py

report-demo: figures
	$(PYTHON) scripts/reporting/make_tables.py --demo
	$(PYTHON) scripts/reporting/make_qc_reports.py --demo
	$(PYTHON) scripts/reporting/make_module_report.py --demo
	$(PYTHON) scripts/reporting/make_failure_analysis.py --demo
	$(PYTHON) scripts/reporting/make_reproducibility_report.py --demo
	$(PYTHON) scripts/reporting/make_manuscript.py --demo

report_demo: report-demo

clean-demo:
	$(PYTHON) -c "import shutil, pathlib; [shutil.rmtree(path, ignore_errors=True) for path in [pathlib.Path('results/demo'), pathlib.Path('results/demo_ablation'), pathlib.Path('results/demo_sensitivity'), pathlib.Path('results/demo_panel'), pathlib.Path('tables/demo')]]"

clean_demo: clean-demo

validate:
	$(PYTHON) -m econiche_opt.cli validate-goals --goal-file docs/goal_status.yml
	$(PYTHON) -m econiche_opt.cli validate-results --demo
	$(PYTHON) -m econiche_opt.cli validate-sources --source-registry config/source_registry.yml
	$(PYTHON) -m econiche_opt.cli validate-project --mode demo

all: install demo test registry-audit preprocess priors train baselines benchmark figures paper validate
