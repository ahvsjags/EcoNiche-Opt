# EcoNiche-Opt v2 全流程落地实施手册

> 用途：把 **“生态位约束多组学优化模型 EcoNiche-Opt”** 从论文想法落地为可复现工程、数据基准、模型实现、基线对照、统计验证、图表和论文初稿。  
> 目标执行者：Claude Code、OpenAI Codex、Cursor Agent、Devin 或任何能读写仓库的 coding agent。  
> 重要原则：**不得伪造数据、不得伪造结果、不得声称超越所有模型，除非 benchmark 结果和统计检验支持。**

---

## 0. 给 Claude Code / Codex 的总提示词

把下面这段直接复制给 Claude Code 或 Codex，作为项目总指令。

```text
你是一个资深计算肿瘤学、生物信息学和机器学习工程 agent。请在当前仓库中实现 EcoNiche-Opt v2 项目。

项目目标：构建一个严格去重、多队列、可复现的免疫检查点抑制剂 ICB response benchmark，并实现 EcoNiche-Opt：一种生态位约束多组学优化模型，用于识别 melanoma anti-PD-1 response / resistance 的多细胞生态位模块。模型需要对比 TIDE、IMPRES、IPRES、IRIS/RDI、EaSIeR、PredictIO_100、经典免疫 signatures、deconvolution signatures 和常规 ML baselines，在 leave-one-dataset-out、locked external melanoma validation、pan-cancer transfer、single-cell mechanism validation、survival validation 和 perturbation prioritization 中完成验证。

执行要求：
1. 不要伪造真实数据结果。真实数据不可下载或需申请时，标记为 ACCESS_RESTRICTED，并提供可插拔接口。
2. 先实现可在 synthetic/demo 数据上跑通的完整 pipeline，再接入真实数据。
3. 所有数据集必须进入 data_registry.yml，记录 accession、cancer type、therapy、platform、sample count、patient count、timepoint、endpoint、access status、download script、preprocessing script、是否用于 discovery/internal/external/pan-cancer/scRNA/survival/perturbation。
4. 任何训练/测试划分必须 patient-level split，不允许同一患者不同时间点或重复样本跨 train/test 泄漏。
5. primary benchmark 必须使用 leave-one-dataset-out；调参只能在 training datasets 内部 nested CV 进行。
6. 所有指标必须输出 AUROC、AUPRC、balanced accuracy、MCC、F1、sensitivity、specificity、ECE、Brier score、calibration curve、decision curve、bootstrap CI。
7. 与 baselines 的比较必须做 paired bootstrap 或 DeLong 检验，并做 FDR correction。
8. 输出清晰的 figures、tables、supplementary tables、methods text 和 reproducibility report。
9. 所有核心函数都要有单元测试，所有 pipeline 要有 smoke test。
10. 必须逐项执行 GOAL-000 到 GOAL-080；每完成一个 GOAL，更新 docs/goal_status.yml。
11. 最终交付一个可复现仓库，包括 README、environment、Dockerfile、Snakemake/Nextflow 或 Makefile、scripts、notebooks、results schema、paper draft skeleton。
```

---

## 0A. Codex 全目标执行清单（必须逐项实现）

> 本节是给 Claude Code / Codex 的 **Goal-driven implementation contract**。后续 Sprint 是时间顺序，本节是验收顺序。Codex 必须把每个 GOAL 转化为代码、脚本、测试、文档或明确的 `ACCESS_RESTRICTED` / `UNAVAILABLE_WITH_REASON` 记录。不得用虚构数据填充真实结果。

### 0A.1 Codex 执行总原则

```text
1. 逐项实现 GOAL-000 至 GOAL-080，不要只实现看起来容易的部分。
2. 每完成一个 GOAL，更新 docs/goal_status.yml，记录 status、files_changed、tests_run、notes、blocking_issues。
3. 对无法访问的数据，不要跳过目标；实现接口、registry 条目、下载占位脚本、错误处理和说明。
4. 所有真实结果必须来自 pipeline 输出，不允许手工填假数。
5. 所有训练、特征选择、方向符号估计、阈值选择、校准和模型选择都不得使用 locked external / holdout 数据。
6. 所有结果表必须包含 cohort、sample_id、patient_id、timepoint、therapy、endpoint、model、metric、estimate、ci_low、ci_high、p_value、q_value。
7. 所有图表必须能由脚本重建，不能只放手工图片。
8. 所有 benchmark 中的缺失 baseline 必须写明 unavailable_with_reason。
9. 任何“outperforms / superior / best”类结论必须由 paired bootstrap / DeLong / FDR 支持，否则输出 RESULT_PENDING 或 NOT_SIGNIFICANT。
10. 项目完成条件不是“代码能跑”，而是 make demo、make test、make benchmark_demo、make report_demo 全部通过。
```

### 0A.2 必须创建的 Goal 状态文件

Codex 必须创建：

```text
docs/goal_status.yml
```

模板：

```yaml
project: EcoNiche-Opt
version: 0.2.0
goal_status_schema: 1
last_updated: null
summary:
  total_goals: 81
  completed: 0
  blocked: 0
  pending: 81
  failed: 0
goals:
  GOAL-000:
    title: Repository governance and no-fabrication policy
    priority: P0
    status: pending
    files_changed: []
    tests_run: []
    blocking_issues: []
    notes: ""
```

同时实现一个验证器：

```bash
python -m econiche_opt.cli validate-goals --goal-file docs/goal_status.yml
```

验证器要求：

```text
1. 所有 GOAL ID 存在。
2. P0 / P1 goal 不能缺 status。
3. completed goal 必须有 files_changed 或 tests_run。
4. blocked goal 必须有 blocking_issues。
5. demo 模式下允许真实数据 goal 标记为 interface_completed，但不允许 silent skip。
```

---

# 0B. 全部 Implementation Goals

## GOAL-000 — 项目治理与不得伪造政策

**Priority:** P0  
**Goal:** 创建项目治理文件，明确真实数据、受限数据、缺失结果和模型比较的边界。

**Codex must implement:**

```text
1. README.md 中加入 no-fabrication policy。
2. docs/reproducibility/no_fabrication_policy.md。
3. results 生成器中凡无真实输出，一律写 RESULT_PENDING。
4. manuscript 自动填充器中禁止写虚构性能数字。
5. 所有 claims 通过 claim_gate 检查。
```

**Deliverables:**

```text
docs/reproducibility/no_fabrication_policy.md
src/econiche_opt/reporting/claim_gate.py
tests/test_claim_gate.py
```

**Acceptance:**

```bash
pytest -q tests/test_claim_gate.py
```

---

## GOAL-001 — 仓库结构初始化

**Priority:** P0  
**Goal:** 创建可复现工程目录结构。

**Codex must implement:**

```text
src/econiche_opt/
scripts/
workflow/
config/
data/raw/
data/processed/
results/
figures/
paper/
docs/
tests/
notebooks/
```

**Deliverables:**

```text
pyproject.toml
requirements.txt
environment.yml
Makefile
Dockerfile
README.md
AGENTS.md
```

**Acceptance:**

```bash
python -m pip install -e .
python -c "import econiche_opt; print(econiche_opt.__version__)"
```

---

## GOAL-002 — AGENTS.md 给 Codex 的仓库级指令

**Priority:** P0  
**Goal:** 创建 Codex / Claude Code 每次进入仓库都会读取的工程说明。

**Codex must implement:**

```text
AGENTS.md 必须包含：
- 项目目标；
- 不得伪造数据；
- patient-level split 要求；
- locked external 禁止泄漏；
- 如何运行测试；
- 如何更新 docs/goal_status.yml；
- 如何处理 ACCESS_RESTRICTED 数据；
- 如何生成 figures / tables / manuscript。
```

**Deliverables:**

```text
AGENTS.md
```

**Acceptance:**

```bash
test -s AGENTS.md
grep -i "no fabrication" AGENTS.md
grep -i "patient-level" AGENTS.md
```

---

## GOAL-003 — data_registry.yml 完整数据登记

**Priority:** P0  
**Goal:** 所有 melanoma、pan-cancer、single-cell、survival、perturbation 数据进入统一 registry。

**Codex must implement:**

```text
1. config/data_registry.yml。
2. 每个数据集记录 accession、PMID/DOI、cancer_type、therapy、platform、timepoint、endpoint、sample_count_reported、patient_count_reported、access_status、download_method、role。
3. 不允许下载脚本里硬编码 accession；必须从 registry 读取。
4. 实现 registry validator。
```

**Minimum registry entries:**

```text
Melanoma bulk ICB:
GSE91061, GSE78220, GSE145996, GSE168204, GSE115821, GSE93157, GSE244982, GSE123728, GSE165745, GSE122220, Liu_DFCI_s41591_2019, Gide_PRJEB23709

Pan-cancer ICB:
GSE136961, GSE176307, IMvigor210_EGAS00001002556, GSE67501, GSE121810, GSE140901, GSE165252, GSE183924

Single-cell / mechanism:
GSE115978, GSE123139

Survival / reference:
TCGA_SKCM_Xena, GDC_TCGA_SKCM

Perturbation / dependency:
LINCS_L1000, CMap, DepMap, DGIdb, DrugBank_optional
```

**Deliverables:**

```text
config/data_registry.yml
src/econiche_opt/data/registry.py
tests/test_registry.py
```

**Acceptance:**

```bash
python -m econiche_opt.cli validate-registry --registry config/data_registry.yml
pytest -q tests/test_registry.py
```

---

## GOAL-004 — 数据访问审计器

**Priority:** P0  
**Goal:** 自动区分 public、controlled、manual、unavailable 数据，保证 pipeline 不因为受限数据崩溃。

**Codex must implement:**

```text
1. check_access_status.py。
2. 对 ACCESS_RESTRICTED 数据生成占位 metadata。
3. 受限数据输出 instructions，而不是假数据。
4. 生成 reports/data_access_audit.tsv。
```

**Deliverables:**

```text
scripts/check_access_status.py
results/audit/data_access_audit.tsv
```

**Acceptance:**

```bash
python scripts/check_access_status.py --registry config/data_registry.yml --out results/audit/data_access_audit.tsv
```

---

## GOAL-005 — GEO 下载器与 metadata 抽取

**Priority:** P0  
**Goal:** 从 GEO 下载 public series metadata、supplementary files、platform annotation。

**Codex must implement:**

```text
1. R 脚本：scripts/data_download/download_geo.R。
2. Python wrapper：src/econiche_opt/data/download_geo.py。
3. 能下载 public GEO metadata。
4. 能处理 GEO supplement 缺失、压缩文件、多平台。
5. 输出 standardized raw metadata。
```

**Deliverables:**

```text
scripts/data_download/download_geo.R
src/econiche_opt/data/download_geo.py
```

**Acceptance:**

```bash
Rscript scripts/data_download/download_geo.R --accession GSE78220 --out data/raw/GSE78220
python -m econiche_opt.cli audit-dataset --accession GSE78220
```

---

## GOAL-006 — 表达矩阵读取与统一格式

**Priority:** P0  
**Goal:** 将 RNA-seq、NanoString、microarray 表达矩阵转换为统一 long/wide schema。

**Codex must implement:**

```text
1. 读取 count/TPM/FPKM/intensity/nCounter。
2. 统一 gene symbol。
3. 处理 probe-to-gene collapse。
4. 输出 data/processed/{dataset}/expression.tsv。
5. 输出 data/processed/{dataset}/expression_qc.tsv。
```

**Deliverables:**

```text
src/econiche_opt/preprocess/expression.py
scripts/preprocess/preprocess_expression.py
tests/test_expression_preprocess.py
```

**Acceptance:**

```bash
pytest -q tests/test_expression_preprocess.py
```

---

## GOAL-007 — metadata harmonization

**Priority:** P0  
**Goal:** 统一 sample、patient、cohort、therapy、timepoint、response、survival 字段。

**Standard columns:**

```text
sample_id
patient_id
cohort
accession
cancer_type
therapy_type
therapy_line
biopsy_timepoint
response_raw
response_binary
response_strict
RECIST_category
PFS_time
PFS_event
OS_time
OS_event
platform
source_file
```

**Codex must implement:**

```text
1. MetadataHarmonizer 类。
2. response label 映射。
3. timepoint 映射。
4. therapy type 映射。
5. 缺失字段报告。
```

**Deliverables:**

```text
src/econiche_opt/preprocess/metadata.py
config/label_maps.yml
tests/test_metadata_harmonization.py
```

**Acceptance:**

```bash
pytest -q tests/test_metadata_harmonization.py
```

---

## GOAL-008 — response label 统一与敏感性标签

**Priority:** P0  
**Goal:** 同时生成主标签和敏感性标签。

**Codex must implement:**

```text
response_binary:
  responder = CR + PR
  non_responder = SD + PD

response_strict:
  responder = CR + PR
  non_responder = PD
  exclude = SD

durable_benefit:
  DCB = CR/PR/SD with PFS >= 6 months
  NDB = PD or PFS < 6 months
```

**Deliverables:**

```text
src/econiche_opt/preprocess/labels.py
tests/test_label_mapping.py
```

**Acceptance:**

```bash
pytest -q tests/test_label_mapping.py
```

---

## GOAL-009 — patient-level 去重与泄漏检测

**Priority:** P0  
**Goal:** 防止同一患者、同一活检、同源样本跨训练/测试泄漏。

**Codex must implement:**

```text
1. patient_id 标准化。
2. sample duplication check。
3. 同一 patient 多 timepoint 只能进入同一个 split。
4. PMID/accession/source overlap audit。
5. 生成 leakage_report.tsv。
```

**Deliverables:**

```text
src/econiche_opt/validation/leakage.py
tests/test_leakage_guard.py
```

**Acceptance:**

```bash
pytest -q tests/test_leakage_guard.py
python -m econiche_opt.cli check-leakage --metadata data/processed/all_metadata.tsv
```

---

## GOAL-010 — rank normalization 与跨平台稳健转换

**Priority:** P0  
**Goal:** 实现 sample-wise rank Gaussian、gene-pair、z-score 等多种转换。

**Codex must implement:**

```text
1. rank_gaussian_normalize。
2. within-cohort zscore。
3. gene-pair feature builder。
4. missing gene imputation strategy。
5. 输出 normalization audit。
```

**Deliverables:**

```text
src/econiche_opt/preprocess/normalization.py
tests/test_normalization.py
```

**Acceptance:**

```bash
pytest -q tests/test_normalization.py
```

---

## GOAL-011 — synthetic demo 数据生成器

**Priority:** P0  
**Goal:** 创建能恢复已知生态位模块的模拟多队列数据，用于 smoke test。

**Codex must implement:**

```text
1. 生成 K 个 cohort。
2. 嵌入 tumor、MHC、effector、dysfunction、CAF、myeloid 六个模块。
3. 加入 batch effect 和平台噪声。
4. 加入 patient/timepoint metadata。
5. 生成 scRNA-like cell state priors、pathways、network、ligand-receptor pairs。
```

**Deliverables:**

```text
scripts/demo/make_synthetic_benchmark.py
data/demo/
tests/test_synthetic_recovery.py
```

**Acceptance:**

```bash
python scripts/demo/make_synthetic_benchmark.py --out data/demo --seed 42
pytest -q tests/test_synthetic_recovery.py
```

---

## GOAL-012 — pathway / network / prior 构建

**Priority:** P0  
**Goal:** 构建 Reactome/MSigDB/STRING/ImmPort/ligand-receptor/cell-state prior 的可插拔接口。

**Codex must implement:**

```text
1. GMT reader。
2. Network edge reader。
3. ligand-receptor pair reader。
4. immune gene prior reader。
5. scRNA cell specificity prior reader。
6. demo prior generator。
```

**Deliverables:**

```text
src/econiche_opt/priors/pathways.py
src/econiche_opt/priors/network.py
src/econiche_opt/priors/cell_state.py
src/econiche_opt/priors/ligand_receptor.py
tests/test_priors.py
```

**Acceptance:**

```bash
pytest -q tests/test_priors.py
```

---

## GOAL-013 — 六状态 EcoNiche 模块定义

**Priority:** P0  
**Goal:** 将候选解定义为六个生态状态模块和跨状态 interaction edges。

**States:**

```text
tumor_dedifferentiation
antigen_presentation_mhc
t_nk_effector
t_cell_dysfunction
caf_ecm_exclusion
myeloid_suppression
```

**Codex must implement:**

```text
1. EcoNicheModule dataclass。
2. 每个 state 的 gene set。
3. interaction edges。
4. module validation。
5. serialization to/from JSON/TSV。
```

**Deliverables:**

```text
src/econiche_opt/model/module.py
tests/test_module.py
```

**Acceptance:**

```bash
pytest -q tests/test_module.py
```

---

## GOAL-014 — 基因方向符号估计

**Priority:** P0  
**Goal:** 只用 training cohorts 估计 gene direction，避免 holdout 泄漏。

**Codex must implement:**

```text
1. estimate_gene_direction(X_train_by_cohort, y_train_by_cohort)。
2. 支持 correlation、logFC、AUC direction。
3. 输出 direction_stability。
4. 在 LODO 中每折重新估计。
5. 测试确保 holdout 不参与方向估计。
```

**Deliverables:**

```text
src/econiche_opt/model/direction.py
tests/test_direction_no_leakage.py
```

**Acceptance:**

```bash
pytest -q tests/test_direction_no_leakage.py
```

---

## GOAL-015 — 模块活性与 EcoNicheScore

**Priority:** P0  
**Goal:** 实现 state activity、interaction score、patient-level EcoNiche resistance score。

**Codex must implement:**

```text
A_iq = sum(s_g * X_ig) / sqrt(|G_q|)
I_iqr = ligand-receptor or pathway interaction score
R_i = beta0 + state terms + interaction terms + optional cancer/therapy intercepts
p_i = sigmoid(R_i)
```

**Deliverables:**

```text
src/econiche_opt/model/scoring.py
tests/test_scoring.py
```

**Acceptance:**

```bash
pytest -q tests/test_scoring.py
```

---

## GOAL-016 — 目标函数全部项实现

**Priority:** P0  
**Goal:** 实现 robust objective，不只优化 AUC。

**Objective terms:**

```text
mean_AUROC_LODO
sd_AUROC_LODO
mean_AUPRC_LODO
Cindex_PFS_OS
ECE
CellSpec
PathwayCoh
NetworkCoh
LRcoh
DirectionStability
ModuleSizePenalty
BatchDependencePenalty
DatasetLeakagePenalty
StateRedundancyPenalty
TherapyConfoundingPenalty
```

**Codex must implement:**

```text
1. 每个 objective term 独立函数。
2. objective_terms.tsv 记录每项数值。
3. 权重来自 config/model_config.yml。
4. ablation 可关闭每一项。
```

**Deliverables:**

```text
src/econiche_opt/model/objective.py
config/model_config.yml
tests/test_objective_terms.py
```

**Acceptance:**

```bash
pytest -q tests/test_objective_terms.py
```

---

## GOAL-017 — 生态位优化算子

**Priority:** P0  
**Goal:** 实现 Niche Initialization、Migration、Predation、Symbiosis、Extinction、Niche Jump。

**Codex must implement:**

```text
1. initialization from priors。
2. mutation using network neighbors。
3. crossover/symbiosis across modules。
4. predation removing overfit genes。
5. extinction/elite selection。
6. migration under leave-one-dataset stress。
```

**Deliverables:**

```text
src/econiche_opt/optimization/operators.py
src/econiche_opt/optimization/evolution.py
tests/test_operators.py
```

**Acceptance:**

```bash
pytest -q tests/test_operators.py
```

---

## GOAL-018 — EcoNicheOpt 主模型类

**Priority:** P0  
**Goal:** 提供统一 fit / predict / score / save / load API。

**Codex must implement:**

```python
model = EcoNicheOpt(config)
model.fit(train_data, priors)
model.predict(test_expression, metadata)
model.score(test_expression, y)
model.save(path)
EcoNicheOpt.load(path)
```

**Deliverables:**

```text
src/econiche_opt/model/econiche.py
tests/test_econiche_api.py
```

**Acceptance:**

```bash
pytest -q tests/test_econiche_api.py
```

---

## GOAL-019 — LODO / nested CV / locked external validation

**Priority:** P0  
**Goal:** 严格实现 leave-one-dataset-out 和 nested CV。

**Codex must implement:**

```text
1. Outer split by dataset/cohort。
2. Inner split only inside training datasets。
3. Locked external datasets never used for model selection。
4. patient-level grouping。
5. therapy/timepoint filters。
```

**Deliverables:**

```text
src/econiche_opt/validation/splits.py
src/econiche_opt/validation/runner.py
tests/test_lodo_splits.py
```

**Acceptance:**

```bash
pytest -q tests/test_lodo_splits.py
```

---

## GOAL-020 — 指标计算完整实现

**Priority:** P0  
**Goal:** 所有模型统一输出同一指标 schema。

**Metrics:**

```text
AUROC
AUPRC
balanced_accuracy
MCC
F1
sensitivity
specificity
precision
NPV
Brier_score
ECE
calibration_slope
calibration_intercept
decision_curve_net_benefit
C_index
HR
logrank_p
```

**Deliverables:**

```text
src/econiche_opt/evaluation/metrics.py
tests/test_metrics.py
```

**Acceptance:**

```bash
pytest -q tests/test_metrics.py
```

---

## GOAL-021 — Bootstrap、DeLong、FDR 统计比较

**Priority:** P0  
**Goal:** 实现模型之间显著性比较。

**Codex must implement:**

```text
1. paired bootstrap CI for AUROC/AUPRC/delta metrics。
2. DeLong test where applicable。
3. permutation test option。
4. Benjamini-Hochberg FDR。
5. superiority table。
```

**Deliverables:**

```text
src/econiche_opt/evaluation/statistics.py
tests/test_statistics.py
```

**Acceptance:**

```bash
pytest -q tests/test_statistics.py
```

---

## GOAL-022 — Baseline signature 框架

**Priority:** P0  
**Goal:** 所有基线统一为可注册、可运行、可失败解释的 scoring interface。

**Codex must implement:**

```python
class BaselineModel:
    name: str
    required_genes: list[str]
    fit(...)
    predict(...)
    score(...)
    unavailable_with_reason: str | None
```

**Deliverables:**

```text
src/econiche_opt/baselines/base.py
src/econiche_opt/baselines/registry.py
tests/test_baseline_interface.py
```

**Acceptance:**

```bash
pytest -q tests/test_baseline_interface.py
```

---

## GOAL-023 — 已发表 ICB biomarker baselines

**Priority:** P0  
**Goal:** 实现或封装主要 published biomarkers。

**Minimum baselines:**

```text
TIDE-compatible score placeholder/interface
IMPRES
IPRES
IFN-gamma score
CYT score
T cell-inflamed GEP / TIG
APM score
TLS score
MPS
C-ECM
ESCS
TIRP
PredictIO_100 placeholder/interface
EaSIeR placeholder/interface
IRIS/RDI placeholder/interface
```

**Rule:**

```text
如果外部工具没有 public implementation 或无法自动运行，必须实现 unavailable_with_reason，并允许用户提供预计算 score。
```

**Deliverables:**

```text
src/econiche_opt/baselines/signatures.py
config/baseline_signatures.yml
tests/test_signature_baselines.py
```

**Acceptance:**

```bash
pytest -q tests/test_signature_baselines.py
```

---

## GOAL-024 — 常规 ML baselines

**Priority:** P0  
**Goal:** 实现同等验证设计下的 ML 对照。

**Models:**

```text
LASSO logistic
Elastic Net logistic
Random Forest
XGBoost if installed, otherwise sklearn GradientBoosting fallback
SVM
MLP
WGCNA + LASSO placeholder/interface
ssGSEA + logistic
```

**Deliverables:**

```text
src/econiche_opt/baselines/ml.py
tests/test_ml_baselines.py
```

**Acceptance:**

```bash
pytest -q tests/test_ml_baselines.py
```

---

## GOAL-025 — Deconvolution / cell abundance baselines

**Priority:** P1  
**Goal:** 实现 xCell、MCP-counter、CIBERSORTx、EPIC 类接口。

**Codex must implement:**

```text
1. 允许用户提供预计算 cell abundance TSV。
2. demo 模式生成 mock cell abundance。
3. sklearn logistic wrapper。
4. 不得伪造真实 deconvolution 输出。
```

**Deliverables:**

```text
src/econiche_opt/baselines/deconvolution.py
config/deconvolution_sources.yml
```

**Acceptance:**

```bash
python -m econiche_opt.cli score-deconvolution --help
```

---

## GOAL-026 — Benchmark 运行器

**Priority:** P0  
**Goal:** 一条命令运行 EcoNiche-Opt 与所有可用 baseline。

**Codex must implement:**

```bash
python -m econiche_opt.cli run-benchmark \
  --registry config/data_registry.yml \
  --config config/model_config.yml \
  --out results/benchmark
```

**Outputs:**

```text
results/benchmark/predictions.tsv
results/benchmark/metrics.tsv
results/benchmark/bootstrap.tsv
results/benchmark/model_rankings.tsv
results/benchmark/unavailable_baselines.tsv
```

**Acceptance:**

```bash
python -m econiche_opt.cli run-benchmark --demo --out results/demo_benchmark
```

---

## GOAL-027 — Superiority gate

**Priority:** P0  
**Goal:** 自动判断是否允许写“优于现有模型”。

**Codex must implement:**

```text
1. 与 strongest baseline 的 delta AUROC / delta AUPRC。
2. paired bootstrap 95% CI。
3. FDR q-value。
4. 至少 70% holdout cohorts 胜出。
5. calibration 不显著劣化。
6. 不满足时 claim = NOT_SUPPORTED。
```

**Deliverables:**

```text
src/econiche_opt/reporting/superiority_gate.py
tests/test_superiority_gate.py
```

**Acceptance:**

```bash
pytest -q tests/test_superiority_gate.py
```

---

## GOAL-028 — Pan-cancer transfer validation

**Priority:** P1  
**Goal:** melanoma-trained module 在 pan-cancer ICB 队列上外推测试。

**Codex must implement:**

```text
1. no-retraining transfer。
2. cancer-type intercept transfer。
3. therapy-specific sensitivity analysis。
4. per-cancer and pooled metrics。
5. failures must be reported, not hidden。
```

**Deliverables:**

```text
src/econiche_opt/validation/pancancer.py
scripts/analysis/run_pancancer_transfer.py
```

**Acceptance:**

```bash
python scripts/analysis/run_pancancer_transfer.py --demo --out results/demo_pancancer
```

---

## GOAL-029 — Survival validation

**Priority:** P1  
**Goal:** 对 TCGA-SKCM / ICB cohorts 的 PFS/OS 做生存验证。

**Codex must implement:**

```text
1. Cox model。
2. Kaplan-Meier high/low risk。
3. log-rank test。
4. C-index。
5. time-dependent AUC placeholder/interface。
6. 多变量协变量接口：age、stage、TMB、PD-L1、therapy。
```

**Deliverables:**

```text
src/econiche_opt/survival/survival.py
scripts/analysis/run_survival_validation.py
tests/test_survival.py
```

**Acceptance:**

```bash
pytest -q tests/test_survival.py
python scripts/analysis/run_survival_validation.py --demo --out results/demo_survival
```

---

## GOAL-030 — TCGA/GDC/Xena 接口

**Priority:** P1  
**Goal:** 实现 TCGA-SKCM expression 与 survival metadata 的下载/读取接口。

**Codex must implement:**

```text
1. Xena public hub download template。
2. GDC API manifest template。
3. 用户本地文件导入。
4. TCGA 仅作为 prognosis / immune microenvironment reference，不作为 ICB response。
```

**Deliverables:**

```text
scripts/data_download/download_tcga_xena.R
src/econiche_opt/data/tcga.py
```

**Acceptance:**

```bash
python -m econiche_opt.cli import-tcga --help
```

---

## GOAL-031 — Single-cell mechanism validation

**Priority:** P1  
**Goal:** 将模块映射到 scRNA cell types 和 cell states。

**Codex must implement:**

```text
1. 读取 scRNA expression + metadata。
2. module score per cell。
3. cell type enrichment。
4. malignant/T/NK/CAF/myeloid 分层。
5. patient-level aggregation。
6. 输出 UMAP-ready table，而非强制重跑所有 scRNA。
```

**Deliverables:**

```text
src/econiche_opt/single_cell/module_mapping.py
scripts/analysis/run_single_cell_mapping.py
```

**Acceptance:**

```bash
python scripts/analysis/run_single_cell_mapping.py --demo --out results/demo_single_cell
```

---

## GOAL-032 — GSE115978 / GSE123139 专用 scRNA 配置

**Priority:** P1  
**Goal:** 为两个 melanoma scRNA 机制数据提供配置和解析模板。

**Codex must implement:**

```text
1. config/single_cell_registry.yml。
2. GSE115978 parser placeholder or downloader。
3. GSE123139 parser placeholder or downloader。
4. cell type label harmonization。
5. 不下载失败时仍能通过 demo。
```

**Deliverables:**

```text
config/single_cell_registry.yml
src/econiche_opt/single_cell/registry.py
```

**Acceptance:**

```bash
python -m econiche_opt.cli validate-single-cell-registry --registry config/single_cell_registry.yml
```

---

## GOAL-033 — Ligand-receptor ecological interaction scoring

**Priority:** P1  
**Goal:** 用 ligand-receptor / pathway edges 构建跨细胞生态位互作分数。

**Codex must implement:**

```text
1. LR pair table reader。
2. state-state LR score。
3. optional cell abundance weighting。
4. interaction ablation。
```

**Deliverables:**

```text
src/econiche_opt/model/interactions.py
tests/test_interactions.py
```

**Acceptance:**

```bash
pytest -q tests/test_interactions.py
```

---

## GOAL-034 — Perturbation reversal prioritization

**Priority:** P1  
**Goal:** 使用 LINCS/CMap/DepMap/DGIdb 生成可干预假设。

**Codex must implement:**

```text
1. resistance vector v_M。
2. perturbation expression delta reader。
3. reversal score = -corr(v_M, delta_d)。
4. DepMap dependency merge。
5. DGIdb drug-gene merge。
6. 输出 candidate perturbations，不写临床治疗结论。
```

**Deliverables:**

```text
src/econiche_opt/perturbation/reversal.py
scripts/analysis/run_perturbation_prioritization.py
tests/test_perturbation.py
```

**Acceptance:**

```bash
pytest -q tests/test_perturbation.py
python scripts/analysis/run_perturbation_prioritization.py --demo --out results/demo_perturbation
```

---

## GOAL-035 — Compressed clinical panel

**Priority:** P1  
**Goal:** 将完整模块压缩为 20–40 gene panel，并测试性能保持。

**Codex must implement:**

```text
1. panel selection by stability + importance + cell specificity。
2. panel size sweep：12, 20, 30, 40, 60 genes。
3. panel performance table。
4. NanoString/qPCR-ready gene list。
```

**Deliverables:**

```text
src/econiche_opt/panel/compress.py
scripts/analysis/run_panel_compression.py
```

**Acceptance:**

```bash
python scripts/analysis/run_panel_compression.py --demo --out results/demo_panel
```

---

## GOAL-036 — Ablation studies

**Priority:** P0  
**Goal:** 证明 EcoNiche-Opt 不是简单换皮模型。

**Required ablations:**

```text
full model
no_cell_specificity
no_pathway_coherence
no_network_coherence
no_LR_interactions
no_batch_penalty
no_direction_stability
no_state_redundancy
AUC_only_objective
random_priors
four_state_model
six_state_model
```

**Deliverables:**

```text
scripts/analysis/run_ablation.py
results/ablation/ablation_metrics.tsv
```

**Acceptance:**

```bash
python scripts/analysis/run_ablation.py --demo --out results/demo_ablation
```

---

## GOAL-037 — Robustness and sensitivity analyses

**Priority:** P1  
**Goal:** 对标签定义、timepoint、therapy、normalization、module size 做敏感性分析。

**Codex must implement:**

```text
1. CR/PR vs SD/PD。
2. CR/PR vs PD excluding SD。
3. DCB vs NDB。
4. pretreatment only。
5. on-treatment separate。
6. anti-PD1 only vs combo。
7. rank-normalized vs z-score vs gene-pair。
8. module size range sensitivity。
```

**Deliverables:**

```text
scripts/analysis/run_sensitivity.py
results/sensitivity/sensitivity_metrics.tsv
```

**Acceptance:**

```bash
python scripts/analysis/run_sensitivity.py --demo --out results/demo_sensitivity
```

---

## GOAL-038 — Calibration and decision curve analysis

**Priority:** P0  
**Goal:** 除 AUC 外评估实际风险分层价值。

**Codex must implement:**

```text
1. calibration curve data。
2. ECE。
3. Brier score。
4. calibration slope/intercept。
5. decision curve net benefit across thresholds。
```

**Deliverables:**

```text
src/econiche_opt/evaluation/calibration.py
src/econiche_opt/evaluation/decision_curve.py
tests/test_calibration_decision_curve.py
```

**Acceptance:**

```bash
pytest -q tests/test_calibration_decision_curve.py
```

---

## GOAL-039 — Figure 1：研究设计和数据流

**Priority:** P1  
**Goal:** 自动生成研究总览图数据和草图。

**Codex must implement:**

```text
1. 数据层级图：melanoma / pan-cancer / scRNA / perturbation。
2. pipeline flow。
3. training / validation / locked test separation。
4. 输出 SVG/PDF/PNG。
```

**Deliverables:**

```text
scripts/figures/figure1_overview.py
figures/figure1_overview.svg
```

**Acceptance:**

```bash
python scripts/figures/figure1_overview.py --demo --out figures/figure1_overview.svg
```

---

## GOAL-040 — Figure 2：EcoNiche-Opt 模型机制图

**Priority:** P1  
**Goal:** 生成六状态生态位模块、目标函数和优化算子图。

**Codex must implement:**

```text
1. six-state niche schematic。
2. objective function block。
3. operators block。
4. module output block。
```

**Deliverables:**

```text
scripts/figures/figure2_model.py
figures/figure2_model.svg
```

**Acceptance:**

```bash
python scripts/figures/figure2_model.py --demo --out figures/figure2_model.svg
```

---

## GOAL-041 — Figure 3：Melanoma benchmark 表现

**Priority:** P1  
**Goal:** 生成 AUROC/AUPRC/MCC/ranking/CI 图。

**Codex must implement:**

```text
1. per-cohort performance。
2. mean ranking。
3. delta vs best baseline。
4. bootstrap CI。
5. significance annotation only if supported。
```

**Deliverables:**

```text
scripts/figures/figure3_benchmark.py
figures/figure3_benchmark.svg
```

**Acceptance:**

```bash
python scripts/figures/figure3_benchmark.py --demo --out figures/figure3_benchmark.svg
```

---

## GOAL-042 — Figure 4：Pan-cancer transfer

**Priority:** P1  
**Goal:** 生成跨癌种迁移图。

**Codex must implement:**

```text
1. per-cancer AUROC/AUPRC。
2. no-retraining vs cancer-intercept。
3. component-level consistency。
4. failure cases visible。
```

**Deliverables:**

```text
scripts/figures/figure4_pancancer.py
figures/figure4_pancancer.svg
```

**Acceptance:**

```bash
python scripts/figures/figure4_pancancer.py --demo --out figures/figure4_pancancer.svg
```

---

## GOAL-043 — Figure 5：Single-cell mechanism mapping

**Priority:** P1  
**Goal:** 生成模块在细胞类型和状态上的定位图。

**Codex must implement:**

```text
1. module score by cell type。
2. malignant/T/NK/CAF/myeloid 分布。
3. UMAP-ready scatter input。
4. state enrichment heatmap。
```

**Deliverables:**

```text
scripts/figures/figure5_single_cell.py
figures/figure5_single_cell.svg
```

**Acceptance:**

```bash
python scripts/figures/figure5_single_cell.py --demo --out figures/figure5_single_cell.svg
```

---

## GOAL-044 — Figure 6：Perturbation and actionability

**Priority:** P1  
**Goal:** 生成扰动反转和候选机制图。

**Codex must implement:**

```text
1. reversal score ranked barplot。
2. DepMap dependency overlay。
3. drug-gene network table。
4. 标注 testable hypothesis，不写 clinical recommendation。
```

**Deliverables:**

```text
scripts/figures/figure6_perturbation.py
figures/figure6_perturbation.svg
```

**Acceptance:**

```bash
python scripts/figures/figure6_perturbation.py --demo --out figures/figure6_perturbation.svg
```

---

## GOAL-045 — Tables and supplementary tables

**Priority:** P0  
**Goal:** 自动生成主表和补充表。

**Required tables:**

```text
Table 1 dataset summary
Table 2 model performance summary
Table 3 superiority tests
Supplementary Table S1 registry
S2 label harmonization
S3 sample deduplication
S4 baseline definitions
S5 per-cohort metrics
S6 ablation
S7 sensitivity
S8 module genes
S9 scRNA mapping
S10 perturbation candidates
S11 unavailable datasets/baselines with reasons
```

**Deliverables:**

```text
scripts/reporting/make_tables.py
results/tables/*.tsv
```

**Acceptance:**

```bash
python scripts/reporting/make_tables.py --demo --out results/demo_tables
```

---

## GOAL-046 — Manuscript skeleton generator

**Priority:** P1  
**Goal:** 自动生成论文初稿骨架，真实结果缺失时写 RESULT_PENDING。

**Codex must implement:**

```text
1. Abstract。
2. Introduction。
3. Results sections 1-6。
4. Methods。
5. Data availability。
6. Code availability。
7. Limitations。
8. Supplementary note。
9. claim_gate 集成。
```

**Deliverables:**

```text
src/econiche_opt/reporting/manuscript.py
scripts/reporting/make_manuscript.py
paper/manuscript.md
```

**Acceptance:**

```bash
python scripts/reporting/make_manuscript.py --demo --out paper/demo_manuscript.md
grep -i "RESULT_PENDING\|demo" paper/demo_manuscript.md
```

---

## GOAL-047 — Reproducibility report

**Priority:** P0  
**Goal:** 生成完整复现报告。

**Codex must implement:**

```text
1. git commit hash。
2. package versions。
3. data registry snapshot。
4. command history。
5. random seeds。
6. unavailable data/baseline reasons。
7. train/test split manifest。
8. leakage report。
```

**Deliverables:**

```text
src/econiche_opt/reporting/reproducibility.py
scripts/reporting/make_reproducibility_report.py
results/reproducibility_report.md
```

**Acceptance:**

```bash
python scripts/reporting/make_reproducibility_report.py --demo --out results/demo_reproducibility_report.md
```

---

## GOAL-048 — CLI 完整命令集合

**Priority:** P0  
**Goal:** 提供统一命令行入口。

**Required commands:**

```text
validate-registry
audit-dataset
download-data
preprocess-data
check-leakage
make-demo
train-econiche
run-benchmark
run-ablation
run-sensitivity
run-pancancer
run-survival
run-single-cell
run-perturbation
make-figures
make-tables
make-manuscript
make-reproducibility-report
validate-goals
```

**Deliverables:**

```text
src/econiche_opt/cli.py
tests/test_cli.py
```

**Acceptance:**

```bash
python -m econiche_opt.cli --help
pytest -q tests/test_cli.py
```

---

## GOAL-049 — Makefile 一键复现

**Priority:** P0  
**Goal:** 提供一键运行 demo、测试、报告的 make targets。

**Required targets:**

```text
make install
make demo
make test
make lint
make benchmark-demo
make report-demo
make clean-demo
make validate
```

**Deliverables:**

```text
Makefile
```

**Acceptance:**

```bash
make demo
make test
make benchmark-demo
make report-demo
```

---

## GOAL-050 — Workflow manager

**Priority:** P1  
**Goal:** 提供 Snakemake 或 Nextflow 工作流。

**Codex must implement:**

```text
1. workflow/Snakefile 或 nextflow.config。
2. rule graph。
3. demo profile。
4. public-data profile。
5. controlled-data skip logic。
```

**Deliverables:**

```text
workflow/Snakefile
workflow/config/demo.yml
```

**Acceptance:**

```bash
snakemake -n -s workflow/Snakefile --configfile workflow/config/demo.yml
```

---

## GOAL-051 — Docker / conda reproducibility

**Priority:** P0  
**Goal:** 固化运行环境。

**Codex must implement:**

```text
1. Dockerfile。
2. environment.yml。
3. requirements.txt。
4. optional R packages list。
5. README 中写安装步骤。
```

**Deliverables:**

```text
Dockerfile
environment.yml
requirements.txt
requirements-r.txt
```

**Acceptance:**

```bash
python -m pip install -r requirements.txt
python -m pip install -e .
```

---

## GOAL-052 — Continuous integration tests

**Priority:** P1  
**Goal:** 提供 GitHub Actions 或本地 CI 配置。

**Codex must implement:**

```text
1. pytest。
2. smoke demo。
3. registry validation。
4. no large data download in CI。
5. artifact-free test mode。
```

**Deliverables:**

```text
.github/workflows/test.yml
```

**Acceptance:**

```bash
pytest -q
```

---

## GOAL-053 — Data manifest and split manifest

**Priority:** P0  
**Goal:** 所有样本和 split 都有可审计 manifest。

**Codex must implement:**

```text
1. data_manifest.tsv。
2. split_manifest.tsv。
3. cohort role: discovery/internal/external/locked/pancancer/scrna/survival。
4. patient-level grouping。
5. timepoint filters。
```

**Deliverables:**

```text
src/econiche_opt/data/manifest.py
results/manifests/data_manifest.tsv
results/manifests/split_manifest.tsv
```

**Acceptance:**

```bash
python -m econiche_opt.cli make-demo --out data/demo
python -m econiche_opt.cli check-leakage --metadata data/demo/metadata.tsv
```

---

## GOAL-054 — Dataset-level QC reports

**Priority:** P1  
**Goal:** 每个数据集输出 QC 报告。

**QC contents:**

```text
n_samples
n_patients
response distribution
timepoint distribution
therapy distribution
missing labels
missing survival
expression gene count
missing gene rate
platform
normalization used
```

**Deliverables:**

```text
src/econiche_opt/qc/dataset_qc.py
scripts/reporting/make_qc_reports.py
```

**Acceptance:**

```bash
python scripts/reporting/make_qc_reports.py --demo --out results/demo_qc
```

---

## GOAL-055 — Module interpretability reports

**Priority:** P1  
**Goal:** 对最终模块生成可解释报告。

**Codex must implement:**

```text
1. state-wise gene lists。
2. gene direction。
3. pathway enrichment。
4. network coherence。
5. cell specificity。
6. cohort stability。
7. feature importance。
```

**Deliverables:**

```text
src/econiche_opt/reporting/module_report.py
scripts/reporting/make_module_report.py
```

**Acceptance:**

```bash
python scripts/reporting/make_module_report.py --demo --out results/demo_module_report.md
```

---

## GOAL-056 — Failure analysis

**Priority:** P1  
**Goal:** 显示模型在哪些 cohort / cancer type / therapy 设置下失败。

**Codex must implement:**

```text
1. low-performance cohort report。
2. label imbalance diagnosis。
3. missing gene diagnosis。
4. therapy/timepoint mismatch diagnosis。
5. domain shift summary。
```

**Deliverables:**

```text
src/econiche_opt/evaluation/failure_analysis.py
scripts/reporting/make_failure_analysis.py
```

**Acceptance:**

```bash
python scripts/reporting/make_failure_analysis.py --demo --out results/demo_failure_analysis.md
```

---

## GOAL-057 — Therapy and timepoint stratification

**Priority:** P1  
**Goal:** 分层评估 anti-PD1、anti-PDL1、combo、baseline、on-treatment。

**Codex must implement:**

```text
1. filter API。
2. stratified metrics。
3. interaction terms for therapy/timepoint。
4. report insufficient sample size。
```

**Deliverables:**

```text
src/econiche_opt/analysis/stratification.py
```

**Acceptance:**

```bash
python -m econiche_opt.cli run-benchmark --demo --stratify therapy_type,biopsy_timepoint --out results/demo_stratified
```

---

## GOAL-058 — Missing gene and platform compatibility

**Priority:** P0  
**Goal:** Baseline 和 EcoNiche 在缺基因平台上仍可审计运行。

**Codex must implement:**

```text
1. required gene coverage。
2. missing gene rate。
3. minimum coverage threshold。
4. unavailable_with_reason if too many missing genes。
5. optional imputation only if declared。
```

**Deliverables:**

```text
src/econiche_opt/preprocess/gene_coverage.py
tests/test_gene_coverage.py
```

**Acceptance:**

```bash
pytest -q tests/test_gene_coverage.py
```

---

## GOAL-059 — Configuration system

**Priority:** P0  
**Goal:** 所有参数从 YAML 读取，可复现实验。

**Codex must implement:**

```text
1. config/model_config.yml。
2. config/benchmark_config.yml。
3. config/figure_config.yml。
4. config/demo_config.yml。
5. config hash 写入结果。
```

**Deliverables:**

```text
src/econiche_opt/config.py
tests/test_config.py
```

**Acceptance:**

```bash
pytest -q tests/test_config.py
```

---

## GOAL-060 — Random seed and determinism

**Priority:** P0  
**Goal:** 优化和 benchmark 可复现。

**Codex must implement:**

```text
1. global seed setter。
2. seed 写入 results metadata。
3. demo repeated run stable within tolerance。
4. stochastic optimizer history saved。
```

**Deliverables:**

```text
src/econiche_opt/utils/random.py
tests/test_determinism.py
```

**Acceptance:**

```bash
pytest -q tests/test_determinism.py
```

---

## GOAL-061 — Optimizer history and diagnostics

**Priority:** P1  
**Goal:** 保存优化过程，便于证明不是一次性碰运气。

**Codex must implement:**

```text
1. generation。
2. best_score。
3. objective terms。
4. module size。
5. mutation/crossover statistics。
6. convergence plot input。
```

**Deliverables:**

```text
src/econiche_opt/optimization/history.py
results/optimization_history.tsv
```

**Acceptance:**

```bash
python -m econiche_opt.cli train-econiche --demo --out results/demo_train
test -s results/demo_train/optimization_history.tsv
```

---

## GOAL-062 — Result schema validator

**Priority:** P0  
**Goal:** 所有输出表的字段和类型可验证。

**Codex must implement:**

```text
1. schemas/*.yml。
2. validate-results CLI。
3. metrics/predictions/modules/tables schema。
4. 缺字段直接报错。
```

**Deliverables:**

```text
schemas/predictions.yml
schemas/metrics.yml
schemas/modules.yml
src/econiche_opt/validation/schema.py
```

**Acceptance:**

```bash
python -m econiche_opt.cli validate-results --demo
```

---

## GOAL-063 — Paper methods text from config

**Priority:** P1  
**Goal:** Methods 中的参数、队列和模型描述由配置和结果自动生成。

**Codex must implement:**

```text
1. data curation methods。
2. preprocessing methods。
3. model formulation methods。
4. validation methods。
5. statistical analysis methods。
6. limitations text。
```

**Deliverables:**

```text
src/econiche_opt/reporting/methods_text.py
paper/methods_generated.md
```

**Acceptance:**

```bash
python -m econiche_opt.cli make-manuscript --demo --out paper/demo_manuscript.md
```

---

## GOAL-064 — Source citation registry

**Priority:** P1  
**Goal:** 论文和 README 中引用的数据集、工具、signature 来源。

**Codex must implement:**

```text
1. config/source_registry.yml。
2. 每个 dataset / baseline / database 有 citation_key、title、source、url、pmid/doi。
3. manuscript generator 使用 citation_key。
4. 缺 citation 的 registry entry 报警。
```

**Deliverables:**

```text
config/source_registry.yml
src/econiche_opt/reporting/citations.py
```

**Acceptance:**

```bash
python -m econiche_opt.cli validate-sources --source-registry config/source_registry.yml
```

---

## GOAL-065 — Dataset download dry-run

**Priority:** P0  
**Goal:** 在不下载大文件的情况下审计下载计划。

**Codex must implement:**

```text
1. --dry-run。
2. expected files。
3. access status。
4. disk space estimate if available。
5. skip controlled access。
```

**Deliverables:**

```text
src/econiche_opt/data/download_plan.py
```

**Acceptance:**

```bash
python -m econiche_opt.cli download-data --registry config/data_registry.yml --dry-run
```

---

## GOAL-066 — Public cohort integration smoke tests

**Priority:** P1  
**Goal:** 对至少一个 public GEO 队列完成 metadata 下载和预处理 smoke test。

**Codex must implement:**

```text
1. GSE78220 smoke test。
2. GSE91061 smoke test if files accessible。
3. 不要求 CI 下载大文件。
4. 本地运行说明。
```

**Deliverables:**

```text
tests/integration/test_geo_smoke.py
```

**Acceptance:**

```bash
pytest -q tests/integration/test_geo_smoke.py -m integration
```

---

## GOAL-067 — Model serialization and reload reproducibility

**Priority:** P0  
**Goal:** 模型保存后可复现预测。

**Codex must implement:**

```text
1. save model JSON/pickle。
2. save selected genes and directions。
3. save scaler/normalization info。
4. reload and predict same probabilities。
```

**Deliverables:**

```text
src/econiche_opt/model/io.py
tests/test_model_serialization.py
```

**Acceptance:**

```bash
pytest -q tests/test_model_serialization.py
```

---

## GOAL-068 — Threshold selection without leakage

**Priority:** P0  
**Goal:** 分类阈值只能在 training/inner CV 中选择。

**Codex must implement:**

```text
1. threshold selection strategy。
2. threshold stored per fold。
3. holdout never used。
4. test with sentinel labels。
```

**Deliverables:**

```text
src/econiche_opt/model/thresholds.py
tests/test_threshold_no_leakage.py
```

**Acceptance:**

```bash
pytest -q tests/test_threshold_no_leakage.py
```

---

## GOAL-069 — Calibration without leakage

**Priority:** P0  
**Goal:** 概率校准只在 training 内部完成。

**Codex must implement:**

```text
1. Platt / isotonic optional。
2. calibration fitted in inner CV only。
3. holdout only evaluated。
4. calibration report。
```

**Deliverables:**

```text
src/econiche_opt/model/calibration.py
tests/test_calibration_no_leakage.py
```

**Acceptance:**

```bash
pytest -q tests/test_calibration_no_leakage.py
```

---

## GOAL-070 — External precomputed scores import

**Priority:** P1  
**Goal:** 允许导入 TIDE、IRIS、EaSIeR、PredictIO 等外部工具预计算分数。

**Codex must implement:**

```text
1. precomputed_scores.tsv schema。
2. sample_id mapping。
3. missing sample audit。
4. benchmark runner 合并预计算 baseline。
```

**Deliverables:**

```text
src/econiche_opt/baselines/precomputed.py
schemas/precomputed_scores.yml
```

**Acceptance:**

```bash
python -m econiche_opt.cli import-precomputed-scores --help
```

---

## GOAL-071 — Secure handling of controlled data paths

**Priority:** P0  
**Goal:** 对 controlled data 只读取用户本地路径，不上传、不硬编码、不泄露。

**Codex must implement:**

```text
1. local_path environment variable support。
2. .gitignore data controlled paths。
3. no controlled raw data committed。
4. audit report redacts sensitive paths。
```

**Deliverables:**

```text
.gitignore
src/econiche_opt/data/controlled_access.py
tests/test_controlled_access_safety.py
```

**Acceptance:**

```bash
pytest -q tests/test_controlled_access_safety.py
```

---

## GOAL-072 — Documentation for real-data execution

**Priority:** P1  
**Goal:** 用户能从 README 按步骤跑 demo 和 public 数据。

**Codex must implement:**

```text
1. Quickstart demo。
2. Public GEO run。
3. Controlled data instructions。
4. Baseline external score instructions。
5. Troubleshooting。
6. Expected outputs。
```

**Deliverables:**

```text
README.md
docs/run_public_data.md
docs/troubleshooting.md
```

**Acceptance:**

```bash
test -s docs/run_public_data.md
test -s docs/troubleshooting.md
```

---

## GOAL-073 — Demo end-to-end execution

**Priority:** P0  
**Goal:** 在没有任何真实数据的情况下完整演示项目。

**Codex must implement:**

```text
make demo 完成：
1. synthetic data。
2. train EcoNiche。
3. run baselines。
4. run benchmark。
5. run ablation。
6. run figures。
7. run tables。
8. run manuscript。
9. validate schemas。
10. update goal status。
```

**Deliverables:**

```text
results/demo_benchmark/
figures/demo/
paper/demo_manuscript.md
```

**Acceptance:**

```bash
make demo
```

---

## GOAL-074 — Real-data pipeline entrypoint

**Priority:** P1  
**Goal:** 用户能用 registry 中的 public 数据运行完整流程。

**Codex must implement:**

```text
make public-data
```

该命令应：

```text
1. 下载 public datasets。
2. 跳过 controlled datasets。
3. 预处理。
4. 生成 manifests。
5. 运行 benchmark。
6. 生成报告。
```

**Deliverables:**

```text
Makefile target public-data
workflow config public_data.yml
```

**Acceptance:**

```bash
make public-data-dry-run
```

---

## GOAL-075 — Claims and manuscript safety checks

**Priority:** P0  
**Goal:** 自动防止论文夸大。

**Codex must implement:**

```text
1. ban phrases unless superiority_gate passed。
2. replace unsupported statements with RESULT_PENDING or NOT_SUPPORTED。
3. distinguish prediction, mechanism, perturbation hypothesis, clinical utility。
4. TCGA survival 不写成 ICB response。
```

**Deliverables:**

```text
src/econiche_opt/reporting/safety_checks.py
tests/test_manuscript_safety.py
```

**Acceptance:**

```bash
pytest -q tests/test_manuscript_safety.py
```

---

## GOAL-076 — License, citation and data-use notes

**Priority:** P1  
**Goal:** 添加开源和数据使用说明。

**Codex must implement:**

```text
1. LICENSE。
2. CITATION.cff。
3. data use disclaimer。
4. controlled data access disclaimer。
5. external tool license notes placeholder。
```

**Deliverables:**

```text
LICENSE
CITATION.cff
docs/data_use.md
```

**Acceptance:**

```bash
test -s CITATION.cff
test -s docs/data_use.md
```

---

## GOAL-077 — Unit test coverage floor

**Priority:** P0  
**Goal:** 核心代码有最低测试覆盖。

**Codex must implement:**

```text
1. pytest。
2. coverage config。
3. P0 modules coverage >= 70%。
4. smoke tests stable。
```

**Deliverables:**

```text
.coveragerc
pytest.ini
```

**Acceptance:**

```bash
pytest --cov=econiche_opt --cov-report=term-missing
```

---

## GOAL-078 — Performance and memory safeguards

**Priority:** P1  
**Goal:** 防止大数据运行无提示失败。

**Codex must implement:**

```text
1. chunked readers where possible。
2. memory usage warnings。
3. max genes / max samples demo options。
4. optimizer timeout / max generations。
5. progress logging。
```

**Deliverables:**

```text
src/econiche_opt/utils/logging.py
src/econiche_opt/utils/resources.py
```

**Acceptance:**

```bash
python -m econiche_opt.cli train-econiche --demo --generations 2 --population-size 5 --verbose
```

---

## GOAL-079 — Error handling and user-friendly diagnostics

**Priority:** P0  
**Goal:** 常见失败给出明确原因。

**Codex must implement:**

```text
1. missing genes。
2. missing labels。
3. one-class cohort。
4. controlled data absent。
5. invalid registry。
6. baseline unavailable。
7. insufficient samples。
```

**Deliverables:**

```text
src/econiche_opt/exceptions.py
tests/test_error_handling.py
```

**Acceptance:**

```bash
pytest -q tests/test_error_handling.py
```

---

## GOAL-080 — Final project validation gate

**Priority:** P0  
**Goal:** 项目最终验收入口。

**Codex must implement:**

```text
python -m econiche_opt.cli validate-project --mode demo
```

该命令检查：

```text
1. registry valid。
2. goals status valid。
3. tests pass marker exists。
4. demo benchmark outputs exist。
5. schema validation pass。
6. leakage report pass。
7. figures/tables/manuscript exist。
8. unsupported claims absent。
9. reproducibility report exists。
```

**Deliverables:**

```text
src/econiche_opt/validation/project_validation.py
```

**Acceptance:**

```bash
python -m econiche_opt.cli validate-project --mode demo
```

---

## 0C. Codex 可直接执行的总任务提示词

将下面这段作为单独任务发给 Codex，要求它在仓库里实现所有目标。

```text
请完整实现 EcoNiche-Opt v2 仓库，不要只写文档。你必须按 Markdown 实施手册中的 GOAL-000 到 GOAL-080 逐项完成。

执行顺序：
1. 创建仓库结构、AGENTS.md、pyproject、Makefile、Dockerfile、README。
2. 创建 docs/goal_status.yml，并为每个 GOAL 写入 pending。
3. 实现 demo synthetic benchmark，让 make demo 能无真实数据跑通。
4. 实现 data_registry.yml、registry validator、access audit、download dry-run。
5. 实现 metadata harmonization、label mapping、patient-level leakage guard、expression preprocessing、normalization。
6. 实现六状态 EcoNicheModule、gene direction、module activity、interaction score、objective terms、生态位优化算子、EcoNicheOpt API。
7. 实现 LODO/nested CV/locked external validation，确保任何 holdout 数据不进入 feature selection、direction estimation、thresholding、calibration 或 model selection。
8. 实现 baseline registry、published signatures、ML baselines、precomputed external score import 和 unavailable_with_reason 机制。
9. 实现 metrics、bootstrap、DeLong/permutation、FDR、superiority gate。
10. 实现 pan-cancer transfer、survival validation、single-cell mapping、perturbation reversal、panel compression、ablation 和 sensitivity analysis 的 demo 与真实数据接口。
11. 实现 figures、tables、manuscript skeleton、methods text、reproducibility report、failure analysis、QC report。
12. 实现 schema validators、claim safety checks、project validation gate。
13. 为所有 P0 目标写 pytest；P1 目标至少提供 smoke test 或 CLI --help 测试。
14. 每完成一个目标，更新 docs/goal_status.yml。
15. 最后运行：make demo && make test && make benchmark-demo && make report-demo && python -m econiche_opt.cli validate-project --mode demo。

硬性禁止：
- 不得伪造真实 cohort 样本数、AUC、P 值或显著性。
- 不得把 controlled access 数据写成 public。
- 不得把 TCGA 生存验证写成 ICB response 验证。
- 不得用 holdout 数据做 feature selection、gene direction、threshold、calibration 或调参。
- 不得因为某个 baseline 难实现就静默跳过；必须记录 unavailable_with_reason。
- 不得在 superiority_gate 未通过时写 outperforms all models。

完成后输出：
- 修改文件清单。
- 已完成 GOAL 清单。
- 未完成/阻塞 GOAL 及原因。
- 运行过的测试命令和结果。
- 下一步需要人工提供的数据或授权。
```

---

## 0D. Goal 分组完成判据

```text
P0 最低可交付：GOAL-000,001,002,003,004,005,006,007,008,009,010,011,012,013,014,015,016,017,018,019,020,021,022,023,024,026,027,036,038,045,047,048,049,051,053,058,059,060,062,065,067,068,069,071,073,075,077,079,080 全部 completed。

P1 投稿级可交付：GOAL-025,028,029,030,031,032,033,034,035,037,039,040,041,042,043,044,046,050,052,054,055,056,057,063,064,066,070,072,074,076,078 全部 completed 或 interface_completed。

最终项目不能只停在 P0。若目标是 Nature 系列投稿前工程底稿，P0 + P1 都必须完成；controlled access 数据只能标记 blocked_by_access，但接口和审计必须完成。
```

---


---

## 1. 项目科学主线

### 1.1 论文题目候选

主标题建议：

```text
EcoNiche-Opt identifies a multicellular resistance niche predicting melanoma anti-PD-1 response across a multi-cohort immunotherapy benchmark
```

更 Nature 风格的备选：

```text
A multicellular ecological resistance niche predicts immune checkpoint blockade response across melanoma cohorts
```

如果 pan-cancer 结果也强，可以升级为：

```text
Ecological niche-constrained optimization reveals multicellular resistance programs to immune checkpoint blockade across melanoma and pan-cancer cohorts
```

### 1.2 中文主线

我们提出 **EcoNiche-Opt**，一种生态位约束多组学优化模型。该模型把免疫治疗反应预测从单基因特征选择提升为跨细胞类型、跨通路、跨队列稳定的肿瘤生态位模块搜索。通过 melanoma ICB 多队列 benchmark、pan-cancer ICB 外部验证、单细胞机制验证、TCGA/GDC 或 UCSC Xena 生存验证和 perturbation databases，我们识别一个由 tumor dedifferentiation、antigen-presentation alteration、T/NK effector insufficiency、T-cell dysfunction、CAF/ECM exclusion 和 myeloid suppression 共同组成的 anti-PD-1 resistance niche。该模块不仅用于 response prediction，也用于 resistance mechanism stratification 和 perturbation-reversal hypothesis generation。

### 1.3 必须保守表述

不要写：

```text
EcoNiche-Opt 已经超越所有模型。
```

应该写：

```text
We benchmarked EcoNiche-Opt against published transcriptomic biomarkers and machine-learning baselines under leave-one-dataset-out and locked external validation. Superiority was claimed only when paired bootstrap confidence intervals and FDR-adjusted tests supported a significant improvement over the strongest baseline.
```

中文：

```text
我们只在严格 paired bootstrap / DeLong / FDR 统计检验支持时，才声称 EcoNiche-Opt 优于最强基线模型。
```

---

## 2. 框架从哪里来：论文引言逻辑

### 2.1 现有方法解决了什么

1. **TIDE**：将 ICB 不响应建模为 T-cell dysfunction 和 T-cell exclusion 两类免疫逃逸机制。  
2. **IMPRES**：使用 immune checkpoint genes 的 pairwise expression relations 预测 melanoma ICB response。  
3. **IPRES**：anti-PD-1 innate resistance 与 mesenchymal transition、cell adhesion、ECM organization、wound healing、angiogenesis 等肿瘤内源性和微环境过程相关。  
4. **EaSIeR**：强调 tumor immune microenvironment 是复杂动态系统，使用 prior knowledge 和 interpretable ML 构建 immune response biomarkers。  
5. **IRIS / RDI**：从 cell-type-specific ligand-receptor interactions 角度识别 ICB resistance interactions。  
6. **PredictIO**：在 pan-cancer ICB 大规模 meta-analysis 中比较 genomic / transcriptomic biomarkers，并提出 de novo expression signature。  
7. **Transcriptomic biomarker benchmarks**：已有系统 benchmark 包含多种 published transcriptomic signatures，提示单一 signature 的可复现性和跨队列泛化仍是挑战。

### 2.2 共同缺口

现有模型大多只覆盖以下部分能力：

```text
response score
single signature
immune activity
T-cell dysfunction / exclusion
pairwise gene relation
ligand-receptor interaction
pan-cancer expression signature
```

但仍缺少一个统一框架，同时满足：

```text
1. 跨队列泛化；
2. 跨细胞类型解释；
3. pathway / network coherence；
4. ligand-receptor ecological interaction；
5. batch / platform 稳健；
6. 模块足够小，能转成 qPCR/NanoString panel；
7. 输出 resistance subtype，而不是只有 response probability；
8. 可生成 perturbation reversal hypothesis。
```

### 2.3 EcoNiche-Opt 的提出

EcoNiche-Opt 将 ICB resistance prediction 从：

```text
flat gene signature scoring
```

升级为：

```text
ecological niche-constrained multicellular module optimization
```

核心创新是搜索一个多细胞 resistance niche，而不是一组随机筛出来的 biomarker genes。

---

## 3. 论文贡献 Contribution

最终论文可以写成 5 个贡献。

### Contribution 1：Benchmark 数据贡献

```text
We establish a rigorously deduplicated multi-cohort melanoma and pan-cancer ICB transcriptomic benchmark with harmonized response labels, treatment annotation, biopsy timing and patient-level splits.
```

中文：

```text
我们构建一个严格去重、多队列、跨癌种的 ICB transcriptomic benchmark，统一 response labels、therapy annotation、biopsy timing 和 patient-level split。
```

### Contribution 2：模型贡献

```text
We formulate ICB resistance prediction as ecological niche-constrained module optimization rather than single-marker or flat gene-signature scoring.
```

### Contribution 3：机制贡献

```text
EcoNiche-Opt identifies a compact multicellular resistance niche linking tumor dedifferentiation, antigen-presentation alteration, CAF/ECM exclusion, T-cell dysfunction and myeloid suppression.
```

### Contribution 4：泛化贡献

```text
EcoNiche-Opt is evaluated under leave-one-dataset-out, locked external melanoma validation and pan-cancer transfer validation, reducing dataset leakage and batch-driven overfitting.
```

### Contribution 5：转化贡献

```text
The final module can be reduced to a clinically measurable panel and used to generate perturbation-reversal hypotheses for combination immunotherapy research.
```

---

## 4. 数据集总体设计

将数据分成四层。

```text
Layer A: melanoma anti-PD-1 / anti-PD-1 + CTLA-4 核心训练与锁定外测
Layer B: pan-cancer ICB 外部泛化验证
Layer C: single-cell / spatial / cell-state 机制验证
Layer D: perturbation / dependency / drug-gene 数据库做可干预假设
```

重要原则：

```text
1. 所有 accession 必须先进入 registry，不要在代码中硬编码。
2. 所有 sample count / patient count 必须由下载后的 metadata 自动统计，不要手填为最终数字。
3. 如果数据需要 dbGaP / EGA / controlled access，标记 ACCESS_RESTRICTED。
4. 如果数据集有多个时间点，primary analysis 只用 pretreatment / baseline。
5. on-treatment 和 progression 样本只能做 secondary / exploratory analysis。
6. 同一患者多个样本必须 patient-level split，不能跨训练测试泄漏。
7. 所有队列都要记录 therapy type：anti-PD1、anti-PDL1、anti-CTLA4、combo、chemo-ICI、neoadjuvant/adjuvant。
```

---

## 5. 数据集清单 data_registry.yml

在仓库中创建：

```text
config/data_registry.yml
```

模板如下。

```yaml
version: 0.2.0
last_updated: "2026-05-04"
project: EcoNiche-Opt

resources:
  geo:
    description: "NCBI Gene Expression Omnibus; public functional genomics repository."
    url: "https://www.ncbi.nlm.nih.gov/geo/"
  gdc:
    description: "NCI Genomic Data Commons for TCGA harmonized data."
    url: "https://gdc.cancer.gov/"
  ucsc_xena:
    description: "UCSC Xena browser and hubs for TCGA / public cohorts."
    url: "https://xenabrowser.net/"
  lincs:
    description: "LINCS L1000 / Connectivity Map perturbation resource."
    url: "https://clue.io/"
  depmap:
    description: "Cancer Dependency Map."
    url: "https://depmap.org/"
  dgidb:
    description: "Drug Gene Interaction Database."
    url: "https://www.dgidb.org/"

states:
  - tumor_dedifferentiation
  - antigen_presentation_mhc
  - tnk_effector
  - tcell_dysfunction
  - caf_ecm_exclusion
  - myeloid_suppression

cohorts:
  - accession: GSE91061
    name: Riaz_CA209038_melanoma
    layer: A
    cancer_type: melanoma
    therapy: nivolumab_or_ICB
    platform: RNA-seq
    timepoints: [pretreatment, on_treatment]
    endpoint: [RECIST_response, PFS, OS]
    role: discovery_or_lodo
    access: public
    priority: high
    notes: "Advanced melanoma ICB cohort; use patient-level split; primary uses pretreatment only."

  - accession: GSE78220
    name: Hugo_antiPD1_melanoma
    layer: A
    cancer_type: melanoma
    therapy: anti-PD1
    platform: RNA-seq
    timepoints: [pretreatment]
    endpoint: [RECIST_response]
    role: classic_external_benchmark
    access: public
    priority: high
    notes: "Classic IPRES-associated melanoma anti-PD1 cohort. Small sample size; do not overclaim."

  - accession: GSE145996
    name: JerbyArnon_or_related_melanoma_antiPD1
    layer: A
    cancer_type: melanoma
    therapy: anti-PD1
    platform: RNA-seq
    timepoints: [pretreatment]
    endpoint: [RECIST_response]
    role: external_small
    access: public
    priority: medium
    notes: "GEO summary reports 52 patients monitored by RECIST 1.1; RNA-seq subset may be smaller. Verify from metadata."

  - accession: GSE168204
    name: Du_pathway_signature_melanoma_antiPD1
    layer: A
    cancer_type: melanoma
    therapy: anti-PD1_based
    platform: RNA-seq
    timepoints: [pretreatment, on_treatment]
    endpoint: [RECIST_response]
    role: pathway_signature_comparator
    access: public
    priority: high
    notes: "Pathway-based signatures for metastatic melanoma anti-PD1 response; use as strong comparator."

  - accession: GSE115821
    name: Auslander_IMPRES_melanoma
    layer: A
    cancer_type: melanoma
    therapy: ICB
    platform: expression
    timepoints: [pretreatment]
    endpoint: [response]
    role: IMPRES_comparator
    access: public
    priority: high
    notes: "Used for IMPRES/pairwise immune checkpoint relation comparisons. Verify overlap with other studies."

  - accession: GSE93157
    name: Prat_mixed_antiPD1_panel
    layer: A_B
    cancer_type: melanoma_NSCLC_HNSC
    therapy: pembrolizumab_or_nivolumab
    platform: NanoString_nCounter_immune_panel
    timepoints: [pretreatment]
    endpoint: [response]
    role: mixed_tumor_external
    access: public
    priority: medium
    notes: "Use melanoma subset for melanoma external; full set for pan-cancer panel transfer."

  - accession: PRJEB23709
    name: Gide_melanoma_ICB
    layer: A
    cancer_type: melanoma
    therapy: anti-PD1_or_combo
    platform: RNA-seq
    timepoints: [pretreatment]
    endpoint: [RECIST_response, PFS, OS]
    role: locked_external_if_accessible
    access: public_or_EBI_verify
    priority: high
    notes: "Verify availability and metadata; split anti-PD1 monotherapy and anti-PD1+anti-CTLA4 combo."

  - accession: Liu_DFCI_melanoma
    name: Liu_NatMed_2019_melanoma
    layer: A
    cancer_type: melanoma
    therapy: anti-PD1
    platform: RNA-seq_WES
    timepoints: [pretreatment]
    endpoint: [RECIST_response, PFS, OS]
    role: locked_external_if_accessible
    access: dbGaP_or_controlled_verify
    priority: high
    notes: "Major external cohort; do not hard-code. If unavailable, document access restriction."

  - accession: GSE244982
    name: acquired_resistance_melanoma_progression
    layer: A_C
    cancer_type: melanoma
    therapy: ICB_progression
    platform: RNA-seq
    timepoints: [progression]
    endpoint: [acquired_resistance]
    role: mechanism_only
    access: public
    priority: medium
    notes: "Do not use for primary response prediction; use for acquired resistance mechanism validation."

  - accession: GSE136961
    name: NSCLC_antiPD1_DCB_NDB
    layer: B
    cancer_type: NSCLC
    therapy: anti-PD1
    platform: expression
    timepoints: [pretreatment]
    endpoint: [durable_clinical_benefit]
    role: pan_cancer_external
    access: public
    priority: high
    notes: "Use for cross-cancer transfer validation."

  - accession: GSE176307
    name: urothelial_ICB_real_world
    layer: B
    cancer_type: urothelial_cancer
    therapy: ICB
    platform: RNA-seq
    timepoints: [pretreatment]
    endpoint: [response, survival]
    role: pan_cancer_external
    access: public
    priority: high
    notes: "Verify response endpoint and treated sample subset."

  - accession: IMvigor210
    name: IMvigor210_urothelial_atezolizumab
    layer: B
    cancer_type: urothelial_cancer
    therapy: anti-PDL1_atezolizumab
    platform: RNA-seq
    timepoints: [pretreatment]
    endpoint: [response, PFS, OS]
    role: pan_cancer_external_large
    access: public_or_package_verify
    priority: high
    notes: "Large anti-PDL1 cohort. Implement optional downloader if data source available."

  - accession: GSE67501
    name: RCC_nivolumab_small
    layer: B
    cancer_type: renal_cell_carcinoma
    therapy: nivolumab
    platform: expression
    timepoints: [pretreatment]
    endpoint: [response]
    role: pan_cancer_external_small
    access: public
    priority: medium
    notes: "Small sample size; use exploratory only."

  - accession: GSE121810
    name: glioblastoma_pembrolizumab
    layer: B
    cancer_type: glioblastoma
    therapy: pembrolizumab
    platform: RNA-seq
    timepoints: [baseline, on_treatment]
    endpoint: [response, OS]
    role: pan_cancer_external_context_specific
    access: public
    priority: medium
    notes: "Different disease biology and neoadjuvant setting; do not pool blindly."

  - accession: GSE140901
    name: HCC_ICI_nCounter
    layer: B
    cancer_type: hepatocellular_carcinoma
    therapy: anti-PD1_or_anti-PDL1_based
    platform: NanoString_nCounter
    timepoints: [pretreatment]
    endpoint: [response]
    role: pan_cancer_external
    access: public
    priority: medium
    notes: "Panel data; intersect genes only."

  - accession: GSE165252
    name: PERFECT_esophageal_atezolizumab
    layer: B
    cancer_type: esophageal_cancer
    therapy: atezolizumab_plus_chemoradiotherapy
    platform: RNA-seq
    timepoints: [baseline, on_treatment, resection]
    endpoint: [pathologic_response, survival]
    role: pan_cancer_dynamic
    access: public
    priority: medium
    notes: "Chemo-radiation confounded; use as secondary dynamic validation."

  - accession: GSE183924
    name: esophageal_GEj_durvalumab
    layer: B
    cancer_type: esophageal_or_GEj_adenocarcinoma
    therapy: durvalumab
    platform: RNA-seq
    timepoints: [baseline]
    endpoint: [relapse_free_survival]
    role: survival_external
    access: public
    priority: medium
    notes: "Survival endpoint rather than RECIST response."

  - accession: GSE115978
    name: melanoma_ecosystem_scRNA
    layer: C
    cancer_type: melanoma
    therapy: mixed_or_ICB_context
    platform: scRNA-seq
    endpoint: [cell_state, immune_evasion]
    role: single_cell_mechanism
    access: public
    priority: high
    notes: "Use for malignant-cell resistance program, T cell exclusion, immune evasion mapping."

  - accession: GSE123139
    name: melanoma_Tcell_states_scRNA
    layer: C
    cancer_type: melanoma
    therapy: ICB_context
    platform: scRNA-seq
    endpoint: [T_cell_state]
    role: single_cell_Tcell_validation
    access: public
    priority: high
    notes: "Use for dysfunctional CD8 T cell and effector/transitional states."
```

---

## 6. 仓库结构

Claude Code / Codex 应创建如下目录。

```text
econiche-opt/
  README.md
  LICENSE
  CITATION.cff
  pyproject.toml
  requirements.txt
  environment.yml
  Dockerfile
  Makefile
  config/
    data_registry.yml
    model_config.yml
    baselines.yml
    figure_config.yml
  data/
    raw/
    interim/
    processed/
    metadata/
    priors/
    external/
    README.md
  notebooks/
    00_registry_audit.ipynb
    01_data_qc.ipynb
    02_model_exploration.ipynb
    03_single_cell_mapping.ipynb
    04_figures.ipynb
  scripts/
    download/
      download_geo.R
      download_geo.py
      download_xena.R
      download_lincs_template.py
      download_depmap_template.py
    preprocess/
      build_metadata.py
      preprocess_bulk.py
      preprocess_nanostring.py
      harmonize_labels.py
      deduplicate_patients.py
      build_gene_universe.py
      build_priors.py
      build_network.py
      build_lr_edges.py
      deconvolution_template.R
    model/
      run_econiche.py
      tune_econiche.py
      score_econiche.py
    baselines/
      run_baselines.py
      score_tide_template.py
      score_impres.py
      score_ipres.py
      score_irirs_rdi_template.py
      score_predictio_template.py
      score_gene_sets.py
      train_ml_baselines.py
    benchmark/
      run_lodo.py
      run_locked_external.py
      run_pan_cancer_transfer.py
      bootstrap_compare.py
      delong_compare.py
      calibration.py
      decision_curve.py
      survival_analysis.py
    single_cell/
      preprocess_scrna.R
      map_modules_scrna.R
      compute_cell_state_priors.R
      ligand_receptor_activity.R
    perturbation/
      lincs_reversal.py
      depmap_prioritize.py
      dgidb_lookup.py
    figures/
      make_fig1_overview.py
      make_fig2_benchmark.py
      make_fig3_module_network.py
      make_fig4_single_cell.py
      make_fig5_survival.py
      make_fig6_perturbation.py
    paper/
      generate_methods_text.py
      generate_result_summaries.py
  src/
    econiche/
      __init__.py
      config.py
      registry.py
      io.py
      qc.py
      labels.py
      normalize.py
      priors.py
      networks.py
      module.py
      scoring.py
      objective.py
      optim.py
      model.py
      baselines.py
      metrics.py
      statistics.py
      survival.py
      plotting.py
      utils.py
  tests/
    test_registry.py
    test_labels.py
    test_normalize.py
    test_module.py
    test_objective.py
    test_optim_smoke.py
    test_lodo_split.py
    test_metrics.py
    test_no_leakage.py
  results/
    demo/
    real/
  figures/
  tables/
  paper/
    manuscript.md
    methods.md
    supplement.md
    references.bib
```

---

## 7. 环境配置

### 7.1 Python requirements.txt

```text
numpy>=1.26
pandas>=2.2
scipy>=1.12
scikit-learn>=1.4
statsmodels>=0.14
lifelines>=0.28
matplotlib>=3.8
networkx>=3.2
pyyaml>=6.0
tqdm>=4.66
joblib>=1.3
pyarrow>=15.0
openpyxl>=3.1
requests>=2.31
beautifulsoup4>=4.12
GEOparse>=2.0.4
```

### 7.2 R packages

```r
install.packages(c(
  "GEOquery",
  "limma",
  "edgeR",
  "sva",
  "GSVA",
  "msigdbr",
  "survival",
  "survminer",
  "UCSCXenaTools",
  "ComplexHeatmap",
  "ggplot2",
  "data.table",
  "dplyr",
  "tidyr",
  "Seurat",
  "SingleCellExperiment",
  "scater",
  "scran"
))
```

如需 Bioconductor：

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c(
  "GEOquery", "limma", "edgeR", "sva", "GSVA", "DESeq2",
  "SingleCellExperiment", "scater", "scran", "celldex"
))
```

---

## 8. 模型数学定义

### 8.1 六个生态状态

EcoNiche-Opt v2 使用六个状态：


a. tumor dedifferentiation  
b. antigen-presentation / MHC  
c. T/NK effector  
d. T-cell dysfunction  
e. CAF/ECM exclusion  
f. myeloid suppression

记为：

\[
\mathcal{Q}=\{q_1,q_2,q_3,q_4,q_5,q_6\}
\]

### 8.2 输入数据

有 \(K\) 个 bulk 队列：

\[
D_k=(X^{(k)}, y^{(k)}, m^{(k)})
\]

其中：

\[
X^{(k)}\in \mathbb{R}^{n_k\times p}
\]

是表达矩阵，\(y^{(k)}\) 是 response label，\(m^{(k)}\) 是 metadata，包括 patient_id、sample_id、cohort、therapy、timepoint、cancer_type、platform。

标签默认：

\[
y_i=1 \quad \text{non-responder / resistant}
\]

\[
y_i=0 \quad \text{responder / sensitive}
\]

主分析：

```text
Responder = CR + PR
Non-responder = SD + PD
```

敏感性分析：

```text
Responder = CR + PR
Non-responder = PD
Exclude = SD
```

DCB / NDB 分析：

```text
DCB = CR/PR/SD with PFS >= 6 months
NDB = PD or PFS < 6 months
```

### 8.3 候选生态位模块

一个候选解：

\[
M=\left(\{G_q\}_{q\in\mathcal{Q}}, \{E_{qr}\}_{q<r}\right)
\]

其中：

```text
G_q   = 状态 q 的基因集合
E_qr  = 状态 q 和 r 之间的 ligand-receptor / pathway / network interaction edges
```

示例：

```python
M = {
    "tumor_dedifferentiation": {"AXL", "NGFR", "ITGA3"},
    "antigen_presentation_mhc": {"HLA-A", "B2M", "TAP1"},
    "tnk_effector": {"GZMB", "PRF1", "NKG7"},
    "tcell_dysfunction": {"PDCD1", "LAG3", "HAVCR2"},
    "caf_ecm_exclusion": {"COL1A1", "FN1", "TGFB1"},
    "myeloid_suppression": {"S100A8", "S100A9", "LILRB2"},
}
```

### 8.4 Rank normalization

跨平台时优先使用 within-sample rank-Gaussian normalization：

\[
\widetilde X_{ig}^{(k)}
=
\Phi^{-1}
\left(
\frac{\mathrm{rank}_i(X_{ig}^{(k)})}{p+1}
\right)
\]

实现要求：

```python
def rank_gaussian_normalize(X: pd.DataFrame) -> pd.DataFrame:
    ranks = X.rank(axis=1, method="average", pct=True)
    ranks = ranks.clip(lower=1e-4, upper=1 - 1e-4)
    return pd.DataFrame(stats.norm.ppf(ranks), index=X.index, columns=X.columns)
```

### 8.5 方向符号

每个基因朝向 non-response 的方向：

\[
s_g=\mathrm{sign}\left(
\frac{1}{K_{train}}\sum_{k\in train}\mathrm{cor}(X_g^{(k)},y^{(k)})
\right)
\]

重要：\(s_g\) 只能在 training cohorts 中估计，不能使用 holdout cohort。

### 8.6 状态模块活性

\[
A_{iq}(M)=
\frac{1}{\sqrt{|G_q|}}
\sum_{g\in G_q}
s_g\widetilde X_{ig}
\]

### 8.7 生态位互作分数

如果有 ligand-receptor 或 pathway interaction edges：

\[
I_{iqr}(M)=
\frac{1}{|E_{qr}|}
\sum_{(l,r)\in E_{qr}}
\widetilde X_{il}\widetilde X_{ir}\widehat C_{iq}\widehat C_{ir}
\]

其中 \(\widehat C_{iq}\) 是样本 \(i\) 中状态 \(q\) 的估计丰度，可来自 MCP-counter、xCell、CIBERSORTx、EPIC 或自定义 marker score。

### 8.8 EcoNiche Resistance Score

\[
R_i=
\beta_0+
\alpha_{c_i}+
\delta_{t_i}+
\sum_{q\in\mathcal{Q}}\theta_q A_{iq}
+
\sum_{q<r}\eta_{qr}I_{iqr}
\]

其中：

```text
alpha_c = cancer-type intercept
Delta_t = therapy-type intercept
A_iq    = state module activity
I_iqr   = ecological interaction activity
```

预测 non-response：

\[
\hat p_i=\sigma(R_i)=\frac{1}{1+e^{-R_i}}
\]

### 8.9 Robust optimization objective

不要只最大化 mean AUC。使用 lower-confidence-bound robust objective：

\[
\mathcal{J}(M)=
\bar{\mathrm{AUROC}}_{LODO}
-
\rho\,\mathrm{SD}(\mathrm{AUROC}_{LODO})
+
\lambda_1\bar{\mathrm{AUPRC}}_{LODO}
+
\lambda_2\mathrm{Cindex}_{PFS/OS}
-
\lambda_3\mathrm{ECE}
+
\Omega_{bio}(M)
-
\Omega_{penalty}(M)
\]

生物学奖励项：

\[
\Omega_{bio}(M)=
a_1\mathrm{CellSpec}(M)
+a_2\mathrm{PathwayCoh}(M)
+a_3\mathrm{NetworkCoh}(M)
+a_4\mathrm{LRcoh}(M)
+a_5\mathrm{DirectionStability}(M)
\]

惩罚项：

\[
\Omega_{penalty}(M)=
b_1\mathrm{ModuleSize}(M)
+b_2\mathrm{BatchDependence}(M)
+b_3\mathrm{DatasetLeakage}(M)
+b_4\mathrm{StateRedundancy}(M)
+b_5\mathrm{TherapyConfounding}(M)
\]

### 8.10 具体目标函数实现

在 `src/econiche/objective.py` 中实现：

```python
score = (
    cfg.w_auc * (auc_mean - cfg.robust_rho * auc_sd)
    + cfg.w_auprc * auprc_mean
    + cfg.w_cindex * cindex_mean
    - cfg.w_ece * ece_mean
    + cfg.w_cell_specificity * cell_spec
    + cfg.w_pathway * pathway_coh
    + cfg.w_network * network_coh
    + cfg.w_lr * lr_coh
    + cfg.w_stability * direction_stability
    - cfg.w_size * size_penalty
    - cfg.w_batch * batch_dependence
    - cfg.w_leakage * leakage_penalty
    - cfg.w_redundancy * state_redundancy
    - cfg.w_therapy_confounding * therapy_confounding
)
```

---

## 9. 生态位优化算子

这些算子是算法创新核心。不要写成动物优化算法，要写成疾病机制约束优化。

### 9.1 Niche Initialization

初始候选模块来自：

```text
DE genes
WGCNA modules
single-cell cell-state markers
ImmPort immune genes
Reactome / MSigDB pathways
STRING network neighbors
known ICB signatures: TIDE, IPRES, IMPRES, IFN-gamma, CYT, TIG, TLS, MPS, C-ECM, ESCS
```

采样概率：

\[
P(g\in G_q)\propto
\alpha_1\pi_{gq}
+
\alpha_2I(g\in P_{immune})
+
\alpha_3I(g\in P_{ICB})
+
\alpha_4|cor(X_g,y)|
+
\alpha_5\mathrm{NetworkDegree}(g)
\]

### 9.2 Migration

跨队列迁移算子：

```text
保留在多个训练队列中方向一致、效果稳定的基因；降低只在单一队列有效的基因权重。
```

### 9.3 Predation

过拟合捕食算子：

\[
\mathrm{Overfit}(g)=
\mathrm{AUC}_{train}(g)-\mathrm{AUC}_{external}(g)+\mathrm{BatchDep}(g)
\]

如果超过阈值，删除或降权。

### 9.4 Symbiosis

共生组合算子：

```text
将两个候选生态位模块的稳定核心保留，互补状态模块重组，形成 child candidate。
```

### 9.5 Extinction

灭绝算子：

```text
淘汰低 robust objective 的候选模块，保留 elite candidates，并按 softmax fitness 采样。
```

### 9.6 Niche Jump

网络邻域跳跃：

\[
P(g_{new}|G_q)\propto
\sum_{h\in G_q}A_{gh}+\pi_{gq}
\]

用 STRING / Reactome / ligand-receptor edges 做 guided mutation。

---

## 10. 编码任务分解

下面每一节都是给 Claude Code / Codex 的独立 sprint。

---

# Sprint 0：初始化仓库

## 任务

1. 创建目录结构。
2. 创建 `pyproject.toml`、`requirements.txt`、`environment.yml`、`Dockerfile`、`Makefile`。
3. 创建 demo synthetic data generator。
4. 创建 smoke test：pipeline 能在 demo 数据上跑完。

## 交付

```text
scripts/make_demo_data.py
scripts/model/run_econiche.py
results/demo/econiche_module.tsv
results/demo/lodo_metrics.tsv
results/demo/objective_terms.json
tests/test_optim_smoke.py
```

## 验收命令

```bash
make setup
make demo
make test
```

## 验收标准

```text
1. make demo 能在 5 分钟内完成。
2. 输出 module table、metrics table、optimization history。
3. pytest 全部通过。
```

---

# Sprint 1：数据 registry 和访问审计

## 任务

实现 `src/econiche/registry.py`。

功能：

```python
load_registry(path) -> dict
validate_registry(registry) -> pd.DataFrame
audit_accession(registry) -> pd.DataFrame
write_registry_report(registry, out)
```

要求：

```text
1. 检查每个 cohort 必须有 accession、layer、cancer_type、therapy、platform、timepoints、endpoint、role、access。
2. 标记 public / controlled / unknown。
3. 生成 tables/dataset_access_audit.tsv。
4. 不要自动跳过 unknown，要报警告。
```

## 输出

```text
tables/dataset_access_audit.tsv
tables/dataset_roles.tsv
```

---

# Sprint 2：GEO 下载和 metadata 抽取

## 任务

实现：

```text
scripts/download/download_geo.R
scripts/download/download_geo.py
scripts/preprocess/build_metadata.py
```

R 版优先用 GEOquery。

## 功能

1. 根据 `config/data_registry.yml` 中 public GEO accession 下载 Series Matrix 或 supplementary files。
2. 提取 sample metadata：

```text
sample_id
patient_id_raw
cohort
accession
platform
title
source_name
characteristics_ch1
therapy
timepoint
response_raw
response_harmonized
pfs
os
pfs_event
os_event
```

3. 所有无法自动解析的字段写入 `needs_manual_curation.tsv`。

## 输出

```text
data/raw/GSE*/
data/metadata/all_geo_samples_raw.tsv
data/metadata/needs_manual_curation.tsv
```

## 验收标准

```text
1. 至少能下载 GSE91061、GSE78220、GSE145996、GSE168204、GSE115821、GSE93157 的 metadata。
2. 不能解析 response 的样本必须进入 needs_manual_curation.tsv。
3. 不允许在代码中硬编码每个样本标签；可以用 accession-specific parser，但必须写在 parser registry 中。
```

---

# Sprint 3：标签 harmonization

## 任务

实现 `src/econiche/labels.py` 和 `scripts/preprocess/harmonize_labels.py`。

## 标签规则

```text
CR, complete response -> R
PR, partial response  -> R
SD, stable disease    -> NR in primary; excluded in sensitivity analysis
PD, progressive disease -> NR
DCB -> R-like endpoint for DCB/NDB datasets
NDB -> NR-like endpoint for DCB/NDB datasets
```

## 输出

```text
data/metadata/metadata_harmonized.tsv
```

## 必须支持三种 endpoint

```text
primary_recist: CR/PR vs SD/PD
strict_recist: CR/PR vs PD, exclude SD
clinical_benefit: DCB vs NDB
```

## 单元测试

```text
tests/test_labels.py
```

---

# Sprint 4：patient-level 去重和泄漏检查

## 任务

实现：

```text
src/econiche/qc.py
scripts/preprocess/deduplicate_patients.py
tests/test_no_leakage.py
```

## 规则

1. 同一 patient_id 多个 sample 不得跨 train/test。
2. 同一 accession 中 baseline 优先于 on-treatment。
3. primary analysis 只保留 pretreatment。
4. on-treatment 单独进入 secondary analysis。
5. progression / acquired resistance 单独进入 mechanism_only。
6. 如果无法解析 patient_id，构造 conservative patient_id：

```text
accession + sample_source + title + possible_patient_token
```

并标记 `patient_id_confidence=low`。

## 输出

```text
data/metadata/metadata_dedup_primary.tsv
data/metadata/metadata_secondary_on_treatment.tsv
data/metadata/metadata_progression.tsv
tables/deduplication_report.tsv
```

## 验收标准

```text
No patient_id appears in both train and test folds.
```

---

# Sprint 5：表达矩阵预处理

## 任务

实现：

```text
src/econiche/normalize.py
scripts/preprocess/preprocess_bulk.py
scripts/preprocess/preprocess_nanostring.py
```

## 输入

```text
GEO series matrix
supplementary processed expression files
raw counts if available
NanoString panel files
```

## 处理规则

RNA-seq：

```text
raw counts -> CPM/TPM-like normalization + log2(x+1)
TPM/FPKM -> log2(x+1)
已 log expression -> detect and keep
```

Microarray：

```text
probe -> gene symbol
multiple probes -> median or max variance
quantile normalization if needed
```

NanoString：

```text
panel gene intersection
log2 transform if needed
rank normalization recommended
```

跨队列模型输入：

```text
primary feature space = intersection of genes across selected cohorts
secondary feature space = union + missingness-aware scoring for gene sets
recommended expression scale = within-sample rank Gaussian
```

## 输出

```text
data/processed/bulk/{cohort}.expr.tsv
data/processed/bulk/{cohort}.metadata.tsv
data/processed/bulk/common_genes_primary.txt
tables/expression_qc_report.tsv
```

---

# Sprint 6：构建基因宇宙和先验

## 任务

实现：

```text
scripts/preprocess/build_gene_universe.py
scripts/preprocess/build_priors.py
src/econiche/priors.py
```

## 基因宇宙来源

```text
1. all common genes in bulk cohorts
2. ICB-related signatures
3. ImmPort immune genes
4. Reactome immune pathways
5. MSigDB Hallmark / C2 / C7
6. single-cell cell-state marker genes
7. STRING network high-confidence genes
8. ligand-receptor genes
```

## 先验矩阵

创建：

```text
data/priors/cell_state_priors.tsv
```

格式：

```text
gene    tumor_dedifferentiation    antigen_presentation_mhc    tnk_effector    tcell_dysfunction    caf_ecm_exclusion    myeloid_suppression
AXL     0.9                         0.0                          0.0             0.0                  0.2                0.0
HLA-A   0.0                         0.9                          0.1             0.0                  0.0                0.0
GZMB    0.0                         0.0                          0.9             0.2                  0.0                0.0
```

## 注意

cell-state prior 可以来自 public scRNA 数据或 marker database，但不能用 holdout cohort 的 response label。

---

# Sprint 7：pathway / network / ligand-receptor 构建

## 任务

实现：

```text
src/econiche/networks.py
scripts/preprocess/build_network.py
scripts/preprocess/build_lr_edges.py
```

## 输入来源

```text
STRING edges
Reactome pathways
MSigDB GMT
CellChatDB / CellPhoneDB ligand-receptor pairs
OmniPath optional
```

## 输出

```text
data/priors/pathways.gmt
data/priors/string_edges.tsv
data/priors/ligand_receptor_edges.tsv
data/priors/pathway_membership.parquet
data/priors/network_adjacency.parquet
```

## 函数

```python
pathway_coherence(module, pathways) -> float
network_coherence(module, edges) -> float
lr_coherence(module, lr_edges, state_pairs) -> float
```

---

# Sprint 8：EcoNiche-Opt 核心模型实现

## 任务

实现：

```text
src/econiche/module.py
src/econiche/scoring.py
src/econiche/objective.py
src/econiche/optim.py
src/econiche/model.py
scripts/model/run_econiche.py
```

## 类设计

```python
@dataclass
class EcoNicheConfig:
    states: list[str]
    min_genes_per_state: int = 3
    max_genes_per_state: int = 25
    population_size: int = 120
    generations: int = 120
    elite_fraction: float = 0.10
    mutation_rate: float = 0.20
    crossover_rate: float = 0.50
    robust_rho: float = 0.50
    random_state: int = 42

    w_auc: float = 1.0
    w_auprc: float = 0.20
    w_cindex: float = 0.10
    w_ece: float = 0.15
    w_cell_specificity: float = 0.25
    w_pathway: float = 0.12
    w_network: float = 0.12
    w_lr: float = 0.12
    w_stability: float = 0.25
    w_size: float = 0.10
    w_batch: float = 0.20
    w_leakage: float = 1.00
    w_redundancy: float = 0.08
    w_therapy_confounding: float = 0.10
```

```python
@dataclass
class EcoNicheModule:
    genes_by_state: dict[str, set[str]]
    edges_by_state_pair: dict[tuple[str, str], set[tuple[str, str]]]
```

```python
@dataclass
class EcoNicheResult:
    best_module: EcoNicheModule
    objective_terms: dict
    lodo_metrics: pd.DataFrame
    predictions: pd.DataFrame
    history: pd.DataFrame
    coefficients: pd.DataFrame

    def module_table(self) -> pd.DataFrame: ...
```

```python
class EcoNicheOpt:
    def __init__(self, config, priors, pathways, network_edges, lr_edges): ...
    def fit(self, X_by_cohort, y_by_cohort, metadata_by_cohort) -> EcoNicheResult: ...
    def score_samples(self, X, metadata=None) -> pd.DataFrame: ...
    def predict_proba(self, X, metadata=None) -> np.ndarray: ...
```

## 验收标准

```text
1. 可在 demo 数据上收敛并找回 planted genes。
2. objective_terms 包含每个项。
3. history 包含 generation、best_score、mean_score、best_auc、best_auprc、size 等。
4. 预测输出包含 sample_id、cohort、true_label、pred_prob、EcoNicheScore、六个 state subscore。
```

---

# Sprint 9：LODO / nested CV / locked external validation

## 任务

实现：

```text
scripts/benchmark/run_lodo.py
scripts/benchmark/run_locked_external.py
src/econiche/statistics.py
```

## LODO 规则

```text
for holdout_cohort in cohorts:
    train_cohorts = all other cohorts
    inner_cv = leave-one-train-cohort-out for tuning
    estimate gene directions only in train_cohorts
    optimize module only in train_cohorts
    train final logistic / calibrated model in train_cohorts
    predict holdout_cohort
```

## 输出

```text
results/real/lodo_predictions.tsv
results/real/lodo_metrics_by_cohort.tsv
results/real/lodo_summary.tsv
results/real/lodo_objective_history.tsv
```

## locked external

```text
1. Train on predefined discovery cohorts.
2. Freeze module and parameters.
3. Evaluate locked external cohorts.
4. No retuning on locked external.
```

---

# Sprint 10：评价指标

## 任务

实现 `src/econiche/metrics.py`。

## 指标

```text
AUROC
AUPRC
balanced_accuracy
accuracy
MCC
F1
sensitivity
specificity
PPV
NPV
Brier score
ECE
calibration slope
calibration intercept
C-index for PFS/OS
integrated Brier score optional
decision curve net benefit
```

## 输出

```text
results/real/metrics_all_models.tsv
results/real/calibration_bins.tsv
results/real/decision_curve.tsv
```

## 注意

阈值选择必须只在 training data 中完成，不能用 holdout data 调阈值。

---

# Sprint 11：统计比较

## 任务

实现：

```text
scripts/benchmark/bootstrap_compare.py
scripts/benchmark/delong_compare.py
src/econiche/statistics.py
```

## 方法

```text
1. paired bootstrap for AUROC / AUPRC / MCC / ECE difference
2. DeLong test for AUROC where applicable
3. Benjamini-Hochberg FDR correction
4. Cohort-stratified bootstrap
5. Patient-level bootstrap
```

## superiority criterion

EcoNiche-Opt 可声称优于 best baseline，必须满足：

```text
1. mean AUROC 或 AUPRC 高于 best baseline；
2. paired bootstrap 95% CI of delta 不跨 0；
3. FDR-adjusted P < 0.05；
4. 在至少 70% locked melanoma external cohorts 上优于 best baseline；
5. calibration / Brier score 不显著劣化。
```

输出：

```text
results/real/model_comparison_bootstrap.tsv
results/real/model_comparison_fdr.tsv
```

---

# Sprint 12：baseline 模型实现

## 任务

实现：

```text
scripts/baselines/run_baselines.py
scripts/baselines/score_gene_sets.py
scripts/baselines/score_impres.py
scripts/baselines/score_ipres.py
scripts/baselines/score_tide_template.py
scripts/baselines/score_irirs_rdi_template.py
scripts/baselines/score_predictio_template.py
scripts/baselines/train_ml_baselines.py
src/econiche/baselines.py
```

## Baseline 清单

### 12.1 Classical biomarkers

```text
CD274 / PD-L1
PDCD1
PDCD1LG2
CXCL9
HLA-DRA
CTLA4
```

### 12.2 Immune activity signatures

```text
CYT score = mean(GZMA, PRF1) or literature version
IFN-gamma score
T cell-inflamed GEP / TIG
EIGS
TLS score
APM score
```

### 12.3 Resistance signatures

```text
IPRES
TIDE dysfunction score
TIDE exclusion score
TIRP
MPS
C-ECM
ESCS
```

### 12.4 Pairwise / rank methods

```text
IMPRES
gene-pair score
within-sample relative expression ordering
```

### 12.5 Systems / interaction models

```text
EaSIeR if code available
IRIS / RDI if code/data available
PredictIO_100 if signature available
EcoTyper ecotype score if feasible
```

### 12.6 Deconvolution

```text
MCP-counter CD8 T cells
MCP-counter cytotoxic lymphocytes
MCP-counter fibroblasts
MCP-counter monocytic lineage
xCell CD8 T
xCell CAF
EPIC immune/stroma estimates
CIBERSORTx optional
```

### 12.7 ML baselines

```text
LASSO logistic
Elastic Net
Random Forest
XGBoost optional
SVM
MLP optional
WGCNA + LASSO
ssGSEA + logistic
```

## Baseline 输出标准

每个 baseline 必须输出统一 prediction table：

```text
sample_id
patient_id
cohort
model_name
score
pred_prob
pred_label
true_label
analysis_endpoint
```

---

# Sprint 13：pan-cancer transfer validation

## 任务

实现：

```text
scripts/benchmark/run_pan_cancer_transfer.py
```

## 三种模式

### 模式 A：No-retraining transfer

```text
Train EcoNiche-Opt on melanoma discovery cohorts.
Freeze module and coefficients.
Evaluate on pan-cancer ICB cohorts.
```

### 模式 B：Cancer-type intercept only

```text
Freeze module genes and state coefficients.
Fit only cancer-type intercept on training subset or via calibration cohort.
```

### 模式 C：State-level recalibration

```text
Freeze module genes.
Refit logistic coefficients for six state scores in each cancer type.
```

## 输出

```text
results/real/pancancer_transfer_metrics.tsv
results/real/pancancer_state_direction.tsv
figures/pancancer_transfer_auc.pdf
```

## 表述要求

不要要求所有癌种都强。正确表述是：

```text
EcoNiche components showed partial pan-cancer transferability, with cancer-type and therapy-specific variation.
```

---

# Sprint 14：生存验证

## 任务

实现：

```text
scripts/benchmark/survival_analysis.py
src/econiche/survival.py
```

## 数据来源

```text
TCGA-SKCM via GDC or UCSC Xena
Liu / DFCI if accessible
Gide / PRJEB23709 if accessible
IMvigor210 if accessible
other cohorts with PFS/OS
```

## 分析

```text
1. Cox proportional hazards model
2. Kaplan-Meier high vs low EcoNiche risk
3. Time-dependent AUC optional
4. Multivariable Cox if covariates available:
   age, sex, stage, cancer type, therapy, TMB, PD-L1, purity
```

## 输出

```text
results/real/survival_cox.tsv
results/real/survival_km.tsv
figures/survival_km.pdf
figures/forest_plot_hr.pdf
```

## 注意

TCGA-SKCM 不是 ICB-treated cohort。它只能证明 prognosis / immune microenvironment association，不能证明 ICB response prediction。

---

# Sprint 15：single-cell mechanism validation

## 任务

实现：

```text
scripts/single_cell/preprocess_scrna.R
scripts/single_cell/compute_cell_state_priors.R
scripts/single_cell/map_modules_scrna.R
scripts/single_cell/ligand_receptor_activity.R
```

## 数据

```text
GSE115978 melanoma ecosystem scRNA-seq
GSE123139 melanoma T cell states scRNA-seq
optional: TISCH2 annotations if using database-level validation
```

## 分析

```text
1. QC cells and genes.
2. Annotate malignant, T/NK, CAF/stromal, myeloid, B/plasma, endothelial.
3. Compute module scores for six EcoNiche states.
4. Compare module scores across cell types.
5. For GSE123139, map tcell_dysfunction to dysfunctional CD8 T states.
6. For GSE115978, map tumor_dedifferentiation and CAF/ECM to immune exclusion / cold niche programs.
7. If patient response available, compare high-risk vs low-risk tumors without leaking labels into model training.
```

## 输出

```text
results/scrna/module_scores_by_cell.tsv
results/scrna/module_scores_by_patient.tsv
results/scrna/cell_type_enrichment.tsv
figures/scrna_umap_module_scores.pdf
figures/scrna_violin_by_cell_type.pdf
figures/scrna_dotplot_state_markers.pdf
```

---

# Sprint 16：perturbation prioritization

## 任务

实现：

```text
scripts/perturbation/lincs_reversal.py
scripts/perturbation/depmap_prioritize.py
scripts/perturbation/dgidb_lookup.py
```

## 目标

不是发现临床治疗方案，而是生成可实验验证的假设：

```text
EcoNiche-Opt prioritizes perturbations predicted to reverse the resistance niche.
```

## Reversal score

\[
\mathrm{Reversal}(d,M)=
-\mathrm{cor}(v_M,\Delta_d)
\]

其中：

```text
v_M(g)=+1 if g is resistance-up gene
v_M(g)=-1 if g is response-up gene
v_M(g)=0 otherwise
Delta_d(g)=perturbation-induced expression change
```

## 整合优先级

```text
priority_score =
  z(reversal_score)
+ z(depmap_dependency_if_targeted)
+ z(drug_gene_evidence)
+ z(target_expressed_in_relevant_state)
- z(toxicity_or_panessential_penalty)
```

## 输出

```text
results/perturbation/lincs_reversal.tsv
results/perturbation/depmap_targets.tsv
results/perturbation/dgidb_hits.tsv
results/perturbation/prioritized_perturbations.tsv
figures/perturbation_heatmap.pdf
```

## 表述限制

只写：

```text
candidate perturbations
putative actionable axes
testable combination hypotheses
```

不要写：

```text
recommended treatment
clinical therapy
validated drug
```

---

# Sprint 17：图表生成

## Figure 1：Study design and EcoNiche-Opt overview

输入：

```text
dataset registry
model schematic
six ecological states
optimization workflow
```

输出：

```text
figures/fig1_overview.pdf
```

内容：

```text
multi-cohort bulk RNA-seq -> cell-state priors -> pathway/network/LR priors -> EcoNiche-Opt -> response score + resistance niche + perturbation hypotheses
```

## Figure 2：Benchmark performance

```text
AUROC / AUPRC per cohort
mean rank across models
paired delta vs best baseline
calibration
decision curve
```

## Figure 3：Identified multicellular resistance niche

```text
six state modules
network graph
pathway enrichment
state coefficients
representative genes
```

## Figure 4：Single-cell mechanism mapping

```text
UMAP module score
cell-type enrichment
T-cell dysfunctional state mapping
malignant-cell resistance state mapping
```

## Figure 5：Survival and clinical association

```text
Kaplan-Meier
Cox forest plot
PFS / OS C-index
```

## Figure 6：Perturbation reversal

```text
LINCS reversal heatmap
DepMap target priority
drug-gene network
```

---

# Sprint 18：论文文本自动生成

## 任务

实现：

```text
scripts/paper/generate_methods_text.py
scripts/paper/generate_result_summaries.py
paper/manuscript.md
paper/methods.md
paper/supplement.md
```

## Methods 必须包含

```text
1. Dataset curation and access status
2. Response label harmonization
3. Patient-level deduplication
4. Expression preprocessing
5. Construction of biological priors
6. EcoNiche-Opt mathematical formulation
7. Optimization algorithm
8. Baseline model scoring
9. Benchmark protocol
10. Statistical comparison
11. Single-cell validation
12. Survival analysis
13. Perturbation prioritization
14. Code and data availability
```

## Result sections

```text
Result 1: A rigorously curated melanoma and pan-cancer ICB transcriptomic benchmark
Result 2: EcoNiche-Opt formulates ICB resistance as multicellular ecological niche optimization
Result 3: EcoNiche-Opt identifies a compact multicellular resistance niche in melanoma
Result 4: EcoNiche-Opt improves response prediction across locked melanoma cohorts
Result 5: Resistance niche components show pan-cancer transferability with context-specific variation
Result 6: Single-cell mapping localizes EcoNiche states to malignant, T/NK, CAF/ECM and myeloid compartments
Result 7: Perturbation-reversal analysis prioritizes candidate axes for experimental validation
```

## 自动生成时的限制

如果结果文件不存在，不要生成“positive result”。生成：

```text
[RESULT_PENDING: run benchmark first]
```

---

## 19. Makefile

Claude Code / Codex 应实现如下命令。

```makefile
setup:
	pip install -r requirements.txt

demo:
	python scripts/make_demo_data.py
	python scripts/model/run_econiche.py --config config/model_config.yml --demo

test:
	pytest -q

registry-audit:
	python scripts/preprocess/audit_registry.py --registry config/data_registry.yml --out tables/dataset_access_audit.tsv

download-geo:
	Rscript scripts/download/download_geo.R config/data_registry.yml

preprocess:
	python scripts/preprocess/build_metadata.py
	python scripts/preprocess/harmonize_labels.py
	python scripts/preprocess/deduplicate_patients.py
	python scripts/preprocess/preprocess_bulk.py

priors:
	python scripts/preprocess/build_gene_universe.py
	python scripts/preprocess/build_priors.py
	python scripts/preprocess/build_network.py
	python scripts/preprocess/build_lr_edges.py

train:
	python scripts/model/run_econiche.py --config config/model_config.yml

baselines:
	python scripts/baselines/run_baselines.py --config config/baselines.yml

benchmark:
	python scripts/benchmark/run_lodo.py
	python scripts/benchmark/bootstrap_compare.py
	python scripts/benchmark/calibration.py
	python scripts/benchmark/decision_curve.py

pancancer:
	python scripts/benchmark/run_pan_cancer_transfer.py

survival:
	python scripts/benchmark/survival_analysis.py

single-cell:
	Rscript scripts/single_cell/preprocess_scrna.R
	Rscript scripts/single_cell/map_modules_scrna.R

perturbation:
	python scripts/perturbation/lincs_reversal.py
	python scripts/perturbation/depmap_prioritize.py
	python scripts/perturbation/dgidb_lookup.py

figures:
	python scripts/figures/make_fig1_overview.py
	python scripts/figures/make_fig2_benchmark.py
	python scripts/figures/make_fig3_module_network.py
	python scripts/figures/make_fig4_single_cell.py
	python scripts/figures/make_fig5_survival.py
	python scripts/figures/make_fig6_perturbation.py

paper:
	python scripts/paper/generate_methods_text.py
	python scripts/paper/generate_result_summaries.py

all: setup demo test registry-audit preprocess priors train baselines benchmark figures paper
```

---

## 20. 关键测试

### test_no_leakage.py

必须检查：

```text
1. 同一 patient_id 不跨 train/test。
2. holdout cohort 不参与 gene direction estimation。
3. holdout cohort 不参与 module optimization。
4. holdout cohort 不参与 threshold selection。
5. locked external 不参与调参。
```

### test_lodo_split.py

检查每个 LODO fold：

```text
train cohorts != holdout cohort
no patient overlap
train labels available
holdout labels available
```

### test_objective.py

检查 objective：

```text
1. 大模块会被 size penalty 惩罚。
2. batch-dependent fake module 会被 batch penalty 惩罚。
3. pathway-coherent module 得分高于随机基因。
4. cell-state-specific module 得分高于错配状态模块。
```

### test_optim_smoke.py

在 synthetic data 上：

```text
1. planted resistance module 被部分找回。
2. AUROC > 0.75。
3. module size 在设定范围内。
4. 运行不报错。
```

---

## 21. 真实数据运行顺序

当 demo 跑通后，按下面顺序接真实数据。

```bash
make registry-audit
make download-geo
make preprocess
make priors
make train
make baselines
make benchmark
make pancancer
make survival
make single-cell
make perturbation
make figures
make paper
```

推荐先跑小范围：

```bash
python scripts/benchmark/run_lodo.py \
  --cohorts GSE91061,GSE78220,GSE145996,GSE168204 \
  --endpoint primary_recist \
  --out results/real/melanoma_core_lodo
```

再扩展：

```bash
python scripts/benchmark/run_lodo.py \
  --cohorts config/cohorts_melanoma_all.txt \
  --endpoint primary_recist \
  --out results/real/melanoma_all_lodo
```

最后 pan-cancer：

```bash
python scripts/benchmark/run_pan_cancer_transfer.py \
  --train-cohorts config/cohorts_melanoma_discovery.txt \
  --test-cohorts config/cohorts_pancancer_external.txt \
  --out results/real/pancancer_transfer
```

---

## 22. 结果报告模板

### 22.1 Benchmark 总表

```text
model_name
cohort
cancer_type
therapy
endpoint
n_patients
n_responders
n_nonresponders
AUROC
AUROC_CI_low
AUROC_CI_high
AUPRC
AUPRC_CI_low
AUPRC_CI_high
MCC
balanced_accuracy
ECE
Brier
calibration_slope
decision_curve_net_benefit
```

### 22.2 Superiority 表

```text
comparison
metric
mean_delta
ci_low
ci_high
p_value
fdr_q
n_cohorts_better
fraction_cohorts_better
claim_supported
```

### 22.3 Module 表

```text
state
gene
direction
selection_frequency
mean_train_correlation
cross_cohort_stability
cell_state_prior
pathway_annotation
network_degree
coefficient
```

### 22.4 Perturbation 表

```text
perturbation_id
perturbation_name
target_gene
mechanism
reversal_score
depmap_score
dgidb_evidence
target_state
priority_score
interpretation
```

---

## 23. Nature 系列投稿前必须补齐

纯公开数据挖掘很难冲高档 Nature 系列。若目标是 Nature Communications / npj Precision Oncology / Nature Computational Science，建议补：

```text
1. 本地独立 melanoma anti-PD1 队列 qPCR / NanoString validation。
2. 盲法预测验证。
3. IHC 或 multiplex immunofluorescence 验证 CAF/ECM、T-cell exclusion、myeloid suppression。
4. 至少一个 perturbation hypothesis 的体外实验：例如 tumor-immune co-culture、CAF-conditioned media、myeloid suppression assay。
5. 代码和数据处理流程完整开源。
6. 预注册 benchmark protocol 或至少在 supplement 中固定 superiority criteria。
```

如果只能做公开数据，目标更现实：

```text
Briefings in Bioinformatics
Bioinformatics
npj Precision Oncology if validation strong
Nature Communications if benchmark + mechanism + external validation very strong
iScience
BMC Medicine / Genome Medicine depending on novelty
Scientific Reports as保底 Nature Portfolio
```

---

## 24. 风险控制清单

### 高风险 1：数据重叠

处理：

```text
建立 sample_id / patient_id / PMID / accession / title / biopsy_timepoint 去重表。
```

### 高风险 2：平台差异导致假信号

处理：

```text
within-sample rank features
batch dependence penalty
leave-one-dataset-out
NanoString panel 单独分析
```

### 高风险 3：小样本 AUC 偏高

处理：

```text
bootstrap CI
cohort-level meta-analysis
report sample size
avoid overclaim
```

### 高风险 4：baseline 不公平

处理：

```text
每个 baseline 只在 training data 调参
同样的 folds
同样的 endpoint
同样的 missing gene handling
```

### 高风险 5：单细胞 label leakage

处理：

```text
single-cell 只做 cell-state prior 和 mechanism validation
不要用 test cohort response label 训练 prior
```

### 高风险 6：模型太复杂审稿人不信

处理：

```text
ablation study
small panel reduction
transparent module genes
open-source code
synthetic recovery test
```

---

## 25. Ablation studies

必须运行：

```text
EcoNiche-Opt full
without cell-state prior
without pathway coherence
without network coherence
without ligand-receptor edges
without batch penalty
without direction stability
without robust lower-bound objective
four-state version
six-state version
flat gene signature version
random modules matched by size
```

输出：

```text
results/real/ablation_metrics.tsv
figures/ablation_barplot.pdf
```

---

## 26. 小 panel 压缩

临床转化要从完整模块压缩到 20–40 genes。

实现：

```text
1. 基于 selection_frequency 排序。
2. 每个 state 至少保留 3 个 genes。
3. 优先保留跨队列稳定、cell-state specific、pathway coherent、测量可靠的 genes。
4. 重新计算 compressed EcoNicheScore。
5. 验证 compressed score 性能。
```

输出：

```text
results/real/compressed_panel_genes.tsv
results/real/compressed_panel_metrics.tsv
```

论文表述：

```text
The compressed panel preserved most of the predictive and mechanistic signal of the full model and may facilitate prospective NanoString/qPCR validation.
```

---

## 27. Manuscript skeleton

创建 `paper/manuscript.md`：

```markdown
# EcoNiche-Opt identifies a multicellular resistance niche predicting melanoma anti-PD-1 response across a multi-cohort immunotherapy benchmark

## Abstract
[RESULT_PENDING]

## Introduction
- ICB improves outcomes but response prediction remains difficult.
- Existing biomarkers include PD-L1, TMB, TIDE, IMPRES, IPRES, pathway signatures, interaction models.
- Existing signatures often lack cross-cohort robustness and multicellular mechanism resolution.
- We propose EcoNiche-Opt.

## Results

### A rigorously curated melanoma and pan-cancer ICB transcriptomic benchmark
[generated dataset summary]

### EcoNiche-Opt formulates ICB resistance as multicellular ecological niche optimization
[model overview]

### EcoNiche-Opt identifies a compact multicellular resistance niche in melanoma
[module results]

### EcoNiche-Opt improves prediction across melanoma external cohorts
[benchmark results]

### EcoNiche components show pan-cancer transferability with context-specific variation
[pan-cancer results]

### Single-cell mapping localizes resistance components to tumor and microenvironment compartments
[single-cell results]

### Perturbation reversal analysis prioritizes candidate axes for experimental validation
[perturbation results]

## Discussion
- Main findings.
- Mechanistic interpretation.
- Clinical translation potential.
- Limitations.
- Need prospective validation.

## Methods
See paper/methods.md

## Data availability
[generated]

## Code availability
[generated]
```

---

## 28. 参考资料和需要核查的来源

Coding agent 不需要在代码里引用这些，但需要用它们理解数据和 baseline。正式论文要在 references.bib 中添加准确文献。

```text
NCBI GEO overview:
https://www.ncbi.nlm.nih.gov/geo/
https://www.ncbi.nlm.nih.gov/geo/info/overview.html

TIDE:
Jiang et al. Signatures of T cell dysfunction and exclusion predict cancer immunotherapy response.
https://pubmed.ncbi.nlm.nih.gov/30127393/
https://pmc.ncbi.nlm.nih.gov/articles/PMC6487502/

IRIS / RDI:
Sahni et al. A machine learning model reveals expansive downregulation of ligand-receptor interactions that enhance lymphocyte infiltration in melanoma with developed resistance to immune checkpoint blockade.
https://www.nature.com/articles/s41467-024-52555-4

PredictIO:
https://github.com/bhklab/PredictIO
Bareche et al. Leveraging big data of immune checkpoint blockade response identifies novel potential therapeutic targets.

Transcriptomic ICI datasets systematic review:
Kovács & Győrffy. Transcriptomic datasets of cancer patients treated with immune-checkpoint inhibitors: a systematic review.
https://pmc.ncbi.nlm.nih.gov/articles/PMC9153191/

Transcriptomic biomarker benchmark:
Kang et al. A Comprehensive Benchmark of Transcriptomic Biomarkers for Immune Checkpoint Inhibitor Response.
https://pmc.ncbi.nlm.nih.gov/articles/PMC10452274/

GEO datasets:
GSE91061: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061
GSE78220: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78220
GSE145996: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145996
GSE168204: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168204
GSE115821: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115821
GSE93157: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93157
GSE136961: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136961
GSE176307: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176307
GSE115978: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115978
GSE123139: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123139
```

---

## 29. 最终验收标准

项目最终必须满足：

```text
1. Demo pipeline 完整跑通。
2. Public GEO cohorts 可下载、可预处理、可进入 benchmark。
3. Controlled access cohorts 被正确标记，且不会导致 pipeline 崩溃。
4. EcoNiche-Opt 模型和 ablation 可运行。
5. 至少 20 个 baseline / signature score 可运行或被标记为 unavailable_with_reason。
6. LODO / locked external / pan-cancer transfer 结果统一输出。
7. single-cell module mapping 可运行。
8. perturbation prioritization 可运行。
9. 所有指标有 bootstrap CI。
10. 所有模型比较有 FDR。
11. figures 和 tables 可自动生成。
12. manuscript skeleton 可自动填入真实结果；结果缺失时显示 RESULT_PENDING。
13. pytest 通过。
14. README 说明如何从 demo 到真实数据复现。
```

---

## 30. 第一周最小可行任务

如果时间有限，第一周只做这些：

```text
Day 1:
  - 创建仓库结构
  - synthetic demo data
  - EcoNiche-Opt core class skeleton

Day 2:
  - rank normalization
  - module scoring
  - LODO split
  - metrics

Day 3:
  - objective function
  - evolutionary optimizer
  - smoke test synthetic recovery

Day 4:
  - GEO metadata downloader
  - label harmonization
  - GSE91061 / GSE78220 / GSE145996 / GSE168204 metadata audit

Day 5:
  - expression preprocessing
  - first real-data run on 2–4 public cohorts
  - LASSO / Elastic Net / gene-set baselines

Day 6:
  - benchmark tables
  - bootstrap comparison
  - calibration / decision curve

Day 7:
  - Figure 1/2 prototype
  - methods text generator
  - README update
```

---

## 31. 第二阶段扩展任务

```text
Week 2:
  - add more melanoma cohorts
  - implement IMPRES / IPRES / TIDE-compatible scoring
  - run ablations
  - implement compressed panel

Week 3:
  - pan-cancer cohorts
  - survival analysis
  - scRNA mapping GSE115978 / GSE123139

Week 4:
  - perturbation analysis
  - figure polish
  - manuscript draft
  - reproducibility report
```

---

## 32. 最后提醒给 coding agent

```text
不要为了让结果好看而改变 validation design。
不要把 holdout data 用于 feature selection、gene direction、threshold selection、calibration 或 module optimization。
不要把 TCGA prognosis 说成 ICB response。
不要把 perturbation database 分析说成已验证药物。
不要在没有真实 benchmark 结果时写“outperformed all models”。
```

最终目标不是做一个普通 GEO signature，而是建立：

```text
public multi-cohort ICB benchmark
+ ecological niche-constrained model
+ strong baseline comparison
+ single-cell mechanism localization
+ perturbation hypothesis generation
+ reproducible codebase
```

这才是 EcoNiche-Opt 可以向 Nature 系列靠近的工程化路径。
