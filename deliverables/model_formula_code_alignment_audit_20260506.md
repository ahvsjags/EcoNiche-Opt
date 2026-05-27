# EcoNiche-Opt 模型公式与代码一致性核查

核查对象：`D:\EcoNiche-Opt\EcoNiche-Opt_模型公式原理创新点贡献.docx`

结论：该 Word 文档与仓库代码不是 100% 一致。它更接近“完整生态位优化框架”的方法学蓝图；仓库中早期 `src/econiche` 实现了其中一部分公式，但当前主结果主要由 `src/econiche_opt/model/endpoint_modules.py` 的 endpoint module pipeline 产生，属于更实证、更简化的“固定模块先验分数 + LODO/校准/基线比较”版本。

## 一致部分

| Word 中的设计 | 代码位置 | 核查结论 |
|---|---|---|
| 多队列输入 `D_k=(X^(k), y^(k))`，跨 cohort 验证 | `src/econiche/model.py`; `src/econiche_opt/model/endpoint_modules.py`; `scripts/model/run_endpoint_module_analysis.py` | 一致。代码有 cohort 字典输入、LODO 留队列验证、train/holdout 分离。 |
| 样本内 rank-based normalization | `src/econiche/normalize.py`; `src/econiche/model.py`; `src/econiche/scoring.py` | 早期 `econiche` 主类一致，使用 sample-wise rank Gaussian normalization。 |
| 训练集中估计 gene direction | `src/econiche/scoring.py` | 早期 `econiche` 实现一致，根据表达与标签相关方向给基因加正负号。 |
| 六个生态状态 | `src/econiche/module.py`; `config/model_config.yml` | 早期 `econiche` 实现一致：tumor dedifferentiation、antigen presentation、T/NK effector、T-cell dysfunction、CAF/ECM exclusion、myeloid suppression。 |
| 状态模块分数按 signed genes 汇总并用 `sqrt(|G_q|)` 归一化 | `src/econiche/scoring.py` | 早期 `econiche` 实现一致。 |
| 稳健目标函数包含 AUROC、AUPRC、ECE、生物学奖励、惩罚项 | `src/econiche/objective.py`; `src/econiche/optim.py` | 部分一致。公式框架存在，但很多项目前是占位或简化值。 |
| LODO、calibration、decision curve、FDR claim gate | `scripts/model/run_endpoint_module_analysis.py`; `src/econiche_opt/reporting/claim_gate.py`; `scripts/benchmark/calibration.py`; `scripts/benchmark/decision_curve.py` | 当前项目实现较充分，符合审稿稳健性要求。 |

## 主要不一致

| 公式/文档内容 | 当前代码实现 | 判断 |
|---|---|---|
| 文档定义 `y=1` 为 non-responder/resistant，`y=0` 为 responder/sensitive | `endpoint_response_label()` 中 CR/PR/R/DCB 返回 `1`，PD/NR/NDB 返回 `0`；当前预测列叫 `response_probability` | 标签方向相反。当前 endpoint 主分析预测的是 response probability，不是 non-response probability。 |
| 文档主分数 `R_i` 表示 non-response/high-risk probability | 当前 endpoint module 里 `EcoNiche-Opt-ModulePriorFixed` 的正向分数经 sigmoid/Platt 后用于 response probability | 语义不一致。论文不能直接写成“耐药概率”，除非改代码或重新定义分数方向。 |
| 文档为六个生态状态 `Q={q1,...,q6}` | 当前 endpoint module 使用 7 个模块：IFN/T-cell inflamed、cytotoxic CD8、exhaustion/checkpoint、antigen presentation、myeloid suppression、stromal exclusion、TRM/TLS | 当前主结果是七模块版本，不是 Word 的六状态版本。 |
| 文档包含候选模块互作边 `{E_qr}` 和 ligand-receptor/pathway/regulatory interaction terms | `EcoNicheModule` 数据结构有 `edges_by_state_pair` 字段，但 endpoint 主分析没有真正使用互作边；`objective_for_metrics()` 中 LR/network/pathway 多为 0 或简化 | 互作边目前主要是概念/接口，不是当前性能结果的真实核心贡献。 |
| 文档公式 (6) 要求样本内 rank-based normalization | 当前 endpoint module 的 `build_module_features()` 是 gene set 均值后做 cohort 内 z-score | 早期 `econiche` 一致；当前主结果不一致。 |
| 文档公式 (8) 使用 signed gene expression 和 `1/sqrt(|G_q|)` | 当前 endpoint module 对模块基因取均值再 z-score，不估计每个 gene 的训练方向 | 当前主结果不一致。 |
| 文档公式 (10) 包含 cancer-type/treatment correction、state activity、interaction score | 当前 endpoint module 没有显式 `alpha_c(i)`、`delta_t(i)` 参数，而是通过 strata、train_pool/holdout 分层处理异质性 | 思想相关，但数学实现不同。 |
| 文档公式 (24)(25) 写了初始化概率和网络邻域跳跃优化算子 | `select_module_from_priors()` 是 prior + correlation 排序选择；`make_optimizer_history()` 是历史记录模拟，不是真实进化优化过程 | 不能把当前代码描述成完整自然启发式/生态位进化优化器。 |

## 最关键风险

1. 如果论文按这个 Word 写，当前代码证据不足以支持“完整 ecological niche-constrained module optimization with interaction edges”的强方法学 claim。
2. 当前 endpoint 主分析的标签方向和 Word 里的耐药概率方向相反，这是必须修正的文字/公式问题。
3. 当前最好性能结果来自七模块固定先验/consensus 思路，而不是 Word 里完整的六状态图互作优化器。
4. 互作边、pathway coherence、network coherence、LR coherence、batch-dependence penalty、therapy confounding penalty 在早期代码中有接口或目标函数字段，但不是充分实装的强证据模块。

## 建议

推荐按当前代码重写方法说明，而不是让代码去硬追 Word 的旧公式。当前更稳妥的模型定义应写为：

`S_i = w_IFN Z_i,IFN + w_CD8 Z_i,CD8 + w_exh Z_i,exh + w_APM Z_i,APM - w_myeloid Z_i,myeloid - w_stroma Z_i,stroma + w_TRM/TLS Z_i,TRM/TLS`

其中每个 `Z_i,m` 是模块基因平均表达的 cohort-level z-score，`S_i` 通过训练 cohort 上选择的阈值或 Platt calibration 映射到 response probability。论文贡献应表述为“cell-state module-prior response model with strict multicohort validation”，不要写成已经完整实现 ligand-receptor graph optimization。

如果坚持让代码严格符合 Word 公式，则需要下一步真正实现：non-response 标签方向、six-state module、rank normalization、signed gene module score、interaction edge score、cancer/treatment offsets、真实优化算子，以及所有生物学奖励/惩罚项的非占位计算。
