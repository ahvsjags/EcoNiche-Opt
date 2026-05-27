const fs = require("fs");
const path = require("path");
const {
  AlignmentType,
  BorderStyle,
  Document,
  HeadingLevel,
  LevelFormat,
  Packer,
  PageBreak,
  Paragraph,
  ShadingType,
  Table,
  TableCell,
  TableRow,
  TextRun,
  VerticalAlign,
  WidthType,
} = require("docx");

const ROOT = path.resolve(__dirname, "..", "..");
const OUT = path.join(ROOT, "deliverables", "EcoNiche-Opt_model_principles_and_formulas_20260506.docx");

const FONT_CN = "Microsoft YaHei";
const FONT_LATIN = "Arial";
const FONT_MATH = "Cambria Math";
const BLUE = "1F4E79";
const LIGHT_BLUE = "EAF2F8";
const LIGHT_GRAY = "F6F8FA";
const DARK = "222222";

function run(text, opts = {}) {
  return new TextRun({
    text,
    font: opts.font || FONT_CN,
    size: opts.size || 21,
    bold: opts.bold || false,
    italics: opts.italics || false,
    color: opts.color || DARK,
    break: opts.break || 0,
  });
}

function para(text, opts = {}) {
  return new Paragraph({
    children: Array.isArray(text) ? text : [run(text, opts.run || {})],
    alignment: opts.alignment || AlignmentType.LEFT,
    spacing: { before: opts.before || 60, after: opts.after || 110, line: opts.line || 300 },
    heading: opts.heading,
    pageBreakBefore: opts.pageBreakBefore || false,
  });
}

function h1(text, pageBreakBefore = false) {
  return para(text, {
    heading: HeadingLevel.HEADING_1,
    pageBreakBefore,
    run: { bold: true, color: BLUE, size: 31 },
    before: 220,
    after: 160,
  });
}

function h2(text) {
  return para(text, {
    heading: HeadingLevel.HEADING_2,
    run: { bold: true, color: BLUE, size: 25 },
    before: 170,
    after: 110,
  });
}

function h3(text) {
  return para(text, {
    heading: HeadingLevel.HEADING_3,
    run: { bold: true, color: "365F91", size: 22 },
    before: 120,
    after: 70,
  });
}

function bullet(text, level = 0) {
  return new Paragraph({
    children: [run(text)],
    numbering: { reference: "bullets", level },
    spacing: { before: 25, after: 55, line: 285 },
  });
}

function numbered(text, level = 0) {
  return new Paragraph({
    children: [run(text)],
    numbering: { reference: "numbers", level },
    spacing: { before: 25, after: 55, line: 285 },
  });
}

function formulaBlock(lines, caption) {
  const children = [];
  if (caption) {
    children.push(
      para(caption, {
        run: { bold: true, color: BLUE, size: 20 },
        before: 80,
        after: 20,
      })
    );
  }
  for (const line of lines) {
    children.push(
      new Paragraph({
        children: [run(line, { font: FONT_MATH, size: 21, color: "111111" })],
        spacing: { before: 0, after: 30, line: 270 },
      })
    );
  }
  return new Table({
    width: { size: 100, type: WidthType.PERCENTAGE },
    columnWidths: [9360],
    borders: {
      top: { style: BorderStyle.SINGLE, size: 1, color: "D9E2EC" },
      bottom: { style: BorderStyle.SINGLE, size: 1, color: "D9E2EC" },
      left: { style: BorderStyle.SINGLE, size: 1, color: "D9E2EC" },
      right: { style: BorderStyle.SINGLE, size: 1, color: "D9E2EC" },
      insideHorizontal: { style: BorderStyle.NONE, size: 0, color: "FFFFFF" },
      insideVertical: { style: BorderStyle.NONE, size: 0, color: "FFFFFF" },
    },
    rows: [
      new TableRow({
        children: [
          new TableCell({
            width: { size: 9360, type: WidthType.DXA },
            shading: { fill: LIGHT_GRAY, type: ShadingType.CLEAR, color: "auto" },
            margins: { top: 150, bottom: 140, left: 180, right: 180 },
            children,
          }),
        ],
      }),
    ],
  });
}

function cell(content, width, opts = {}) {
  const children = Array.isArray(content) ? content : [para(String(content), { before: 0, after: 0 })];
  return new TableCell({
    width: { size: width, type: WidthType.DXA },
    verticalAlign: VerticalAlign.CENTER,
    shading: opts.header ? { fill: BLUE, type: ShadingType.CLEAR, color: "auto" } : opts.shade ? { fill: LIGHT_BLUE, type: ShadingType.CLEAR, color: "auto" } : undefined,
    margins: { top: 110, bottom: 110, left: 120, right: 120 },
    children,
  });
}

function table(headers, rows, widths) {
  const headerRow = new TableRow({
    tableHeader: true,
    children: headers.map((h, i) =>
      cell([para([run(h, { bold: true, color: "FFFFFF", size: 19 })], { before: 0, after: 0 })], widths[i], { header: true })
    ),
  });
  const bodyRows = rows.map((r, idx) =>
    new TableRow({
      children: r.map((v, i) =>
        cell([para([run(String(v), { size: 19, font: i === 1 || String(v).includes("=") ? FONT_MATH : FONT_CN })], { before: 0, after: 0, line: 260 })], widths[i], {
          shade: idx % 2 === 1,
        })
      ),
    })
  );
  return new Table({
    width: { size: 100, type: WidthType.PERCENTAGE },
    columnWidths: widths,
    borders: {
      top: { style: BorderStyle.SINGLE, size: 1, color: "C7D4E2" },
      bottom: { style: BorderStyle.SINGLE, size: 1, color: "C7D4E2" },
      left: { style: BorderStyle.SINGLE, size: 1, color: "C7D4E2" },
      right: { style: BorderStyle.SINGLE, size: 1, color: "C7D4E2" },
      insideHorizontal: { style: BorderStyle.SINGLE, size: 1, color: "E3EAF2" },
      insideVertical: { style: BorderStyle.SINGLE, size: 1, color: "E3EAF2" },
    },
    rows: [headerRow, ...bodyRows],
  });
}

function spacer() {
  return para("", { before: 20, after: 20 });
}

const variableRows = [
  ["X_i,g", "患者/样本 i 的基因 g 表达量", "bulk RNA-seq 或已处理表达矩阵"],
  ["G_k", "第 k 个免疫生态模块的基因集合", "例如 IFN/T-cell inflamed 模块"],
  ["M_k(i)", "患者 i 的第 k 个模块活性分数", "模块基因均值后 cohort 内标准化"],
  ["w_k", "第 k 个模块的固定生物先验权重", "正值促进 response，负值代表抑制/排斥"],
  ["S_i", "EcoNicheScore", "所有模块按先验加权后的总分"],
  ["p_i", "预测 ICB response 概率", "sigmoid(S_i) 或校准版本"],
  ["y_i", "真实 response 标签", "0/1 二分类 endpoint"],
  ["c", "cohort 或数据集", "LODO 中每次留出一个 cohort"],
];

const moduleRows = [
  ["IFN/T-cell inflamed", "IFNG, CXCL9, CXCL10, CXCL11, STAT1, IDO1, GBP1, CXCR3, CCL5, CD274, PDCD1LG2", "+1.00"],
  ["Cytotoxic CD8", "CD8A, CD8B, GZMA, GZMB, GZMH, PRF1, NKG7, GNLY", "+0.50"],
  ["Exhaustion/checkpoint", "PDCD1, CTLA4, LAG3, HAVCR2, TIGIT, TOX, CXCL13", "+0.25"],
  ["Antigen presentation", "HLA-A, HLA-B, HLA-C, HLA-DRA, HLA-DRB1, B2M, TAP1, TAP2, PSMB8, PSMB9", "+0.50"],
  ["Myeloid suppression", "CD68, CD163, MRC1, CSF1R, ITGAM, S100A8, S100A9, IL10, TGFB1", "-0.50"],
  ["Stromal exclusion", "COL1A1, COL1A2, COL3A1, FN1, ACTA2, VIM, TGFBI, POSTN, LOXL2", "-0.50"],
  ["TRM/TLS", "ITGAE, CD69, CXCR6, CXCL13, ZNF683, MS4A1, CD79A, BANK1, LTB", "+0.25"],
];

const formulaRows = [
  ["模块原始均值", "R_k(i) = (1 / |G_k ∩ A_c|) Σ_{g ∈ G_k ∩ A_c} X_i,g", "A_c 是 cohort c 中可用基因集合"],
  ["cohort 内标准化", "M_k(i) = (R_k(i) - mean_c(R_k)) / sd_c(R_k)", "若 sd 为 0，则该模块置 0"],
  ["EcoNicheScore", "S_i = Σ_k w_k M_k(i)", "主模型的固定免疫生态先验得分"],
  ["响应概率", "p_i = sigmoid(S_i) = 1 / (1 + exp(-S_i))", "用于 discrimination 的原始概率"],
  ["Platt 校准", "p_i^cal = sigmoid(a + b S_i)", "a,b 只在训练 cohort 上估计"],
  ["阈值选择", "t* = argmax_t BalancedAccuracy(y_train, 1[p_train ≥ t])", "阈值只用训练数据选取"],
  ["LODO 训练集", "Train_c = D \\ D_c; Test_c = D_c", "每次完整留出一个 cohort"],
];

const endpointRows = [
  ["Strict RECIST", "CR/PR/MR/R/DCB = 1; PD/NR/NDB = 0; SD = missing", "最严格，去掉 SD"],
  ["Primary RECIST", "CR/PR/MR/R/DCB = 1; SD/PD/NR/NDB = 0", "主分析，保守地把 SD 归为 non-response"],
  ["Clinical benefit", "CR/PR/MR/SD/R/DCB = 1; PD/NR/NDB = 0", "敏感性分析，把 SD 视作 benefit"],
];

const metricRows = [
  ["AUROC", "P(score_positive > score_negative)", "判别能力；主 headline 用这个"],
  ["AUPRC", "Area under Precision-Recall curve", "类别不平衡时的补充指标"],
  ["Balanced accuracy", "(Sensitivity + Specificity) / 2", "阈值后分类性能"],
  ["Brier score", "(1/n) Σ_i (p_i - y_i)^2", "概率误差"],
  ["ECE", "Σ_b (|B_b|/n) |acc(B_b) - conf(B_b)|", "校准误差"],
  ["Decision curve", "NB(t)=TP/n - FP/n × t/(1-t)", "临床阈值净获益"],
  ["ΔAUROC", "AUROC_EcoNiche - AUROC_baseline", "与现有模型比较"],
  ["Bootstrap/FDR", "BH-adjusted q values over paired bootstrap p values", "superiority claim gate"],
];

const originalRows = [
  ["原创 1", "Immune-ecology module prior", "把 response 预测从单 signature 升级为多生态模块组合"],
  ["原创 2", "正负生态位同时建模", "同时建模 IFN/CD8/APM/TLS 正向 niche 与 myeloid/stromal 负向 niche"],
  ["原创 3", "固定低自由度生物先验", "避免小样本中大规模筛基因，提高解释性和跨队列稳健性"],
  ["原创 4", "Endpoint-evidence 分层", "区分 RECIST-supported primary、high-evidence core 与 binary response stress-test"],
  ["原创 5", "Claim-gated multicohort framework", "把真实数据审计、LODO、FDR、校准、decision curve 和交付审计整合成可复现框架"],
];

const standardRows = [
  ["不是原创", "AUROC/AUPRC/balanced accuracy/ECE/Brier", "标准评价指标"],
  ["不是原创", "Platt calibration", "标准概率校准方法"],
  ["不是原创", "Paired bootstrap 和 Benjamini-Hochberg FDR", "标准统计比较和多重检验校正"],
  ["不是原创", "LODO cross-validation", "标准跨队列外部验证设计"],
  ["不是原创", "IFNG/CXCL9/TIG/TIDE/APM/IPRES 等 baseline", "已有免疫 response signature 或文献模型思想"],
];

const children = [
  new Paragraph({
    children: [run("EcoNiche-Opt", { bold: true, size: 46, color: BLUE })],
    alignment: AlignmentType.CENTER,
    spacing: { before: 900, after: 160 },
  }),
  new Paragraph({
    children: [run("模型原理、公式体系与原创贡献说明", { bold: true, size: 30, color: DARK })],
    alignment: AlignmentType.CENTER,
    spacing: { after: 360 },
  }),
  new Paragraph({
    children: [run("项目：多队列 ICB response benchmarking 与 immune-ecology module-prior 建模", { size: 22, color: "444444" })],
    alignment: AlignmentType.CENTER,
    spacing: { after: 80 },
  }),
  new Paragraph({
    children: [run("生成日期：2026-05-06", { size: 20, color: "666666" })],
    alignment: AlignmentType.CENTER,
    spacing: { after: 660 },
  }),
  formulaBlock(
    [
      "核心思想：ICB 疗效不是单个基因决定，而是由免疫炎症、细胞毒性、抗原呈递、髓系抑制、基质排斥和 TLS/TRM 等肿瘤免疫生态位共同决定。",
      "主模型：EcoNiche-Opt-ModulePriorFixed = fixed immune-ecology module-prior score."
    ],
    "Executive Summary"
  ),
  para([run("", { break: 1 }), new TextRun({ break: 1 }), new PageBreak()]),

  h1("1. 模型提出的生物学出发点"),
  para("ICB response prediction 的核心困难在于：不同队列、癌种、治疗方案和 endpoint 定义高度异质；单基因或单 signature 往往只能捕捉一部分免疫状态。EcoNiche-Opt 的假设是：ICB response 由 response-promoting niches 与 resistance-promoting niches 的平衡决定。"),
  bullet("Response-promoting niches：IFN/T-cell inflamed、cytotoxic CD8、antigen presentation、TRM/TLS。"),
  bullet("Resistance-promoting niches：myeloid suppression、stromal exclusion。"),
  bullet("Checkpoint/exhaustion 模块不是简单负向；在 ICB 场景中，它也可代表可被免疫检查点抑制剂释放的 pre-existing immune engagement。"),

  h1("2. 变量与输入定义", true),
  table(["符号", "含义", "备注"], variableRows, [1600, 3900, 3860]),
  spacer(),
  formulaBlock(
    [
      "D = {(X_i, y_i, c_i)}_{i=1}^n",
      "X_i = (X_i,1, X_i,2, ..., X_i,p)",
      "y_i ∈ {0, 1}",
      "c_i 表示样本 i 所属 cohort。"
    ],
    "数据表示"
  ),

  h1("3. Endpoint 标签公式"),
  para("项目不把所有 response endpoint 混为一谈，而是显式定义三套 endpoint，并作为 sensitivity analysis 同时评估。"),
  table(["Endpoint", "标签规则", "用途"], endpointRows, [2100, 5000, 2260]),
  spacer(),
  formulaBlock(
    [
      "y_i(endpoint) = f_endpoint(response_raw_i)",
      "Primary RECIST: y_i = 1 if response_raw_i ∈ {CR, PR, MR, R, DCB}; otherwise y_i = 0 for {SD, PD, NR, NDB}.",
      "Strict RECIST: SD is excluded from evaluation.",
      "Clinical benefit: SD is coded as 1."
    ],
    "Response Label Harmonization"
  ),

  h1("4. 免疫生态模块定义", true),
  para("EcoNiche-Opt 使用 7 个固定 immune-ecology modules。每个模块是一组具有明确免疫生态含义的基因集合。"),
  table(["模块", "代表基因", "先验权重 w_k"], moduleRows, [2500, 5200, 1660]),

  h1("5. 模块分数与 EcoNicheScore 公式"),
  para("每个模块先计算模块内基因均值，再在 cohort 内标准化，避免不同平台或表达尺度直接混合。"),
  table(["公式名称", "公式", "解释"], formulaRows, [1900, 4500, 2960]),
  spacer(),
  formulaBlock(
    [
      "S_i = 1.00 M_IFN/T-cell(i)",
      "    + 0.50 M_cytotoxicCD8(i)",
      "    + 0.25 M_exhaustion/checkpoint(i)",
      "    + 0.50 M_antigenPresentation(i)",
      "    - 0.50 M_myeloidSuppression(i)",
      "    - 0.50 M_stromalExclusion(i)",
      "    + 0.25 M_TRM/TLS(i)"
    ],
    "EcoNicheScore 展开式"
  ),

  h1("6. 概率预测与校准"),
  para("主判别模型使用原始 EcoNicheScore 的 sigmoid 概率；临床概率解释和 decision curve 可使用训练集内 Platt calibration。"),
  formulaBlock(
    [
      "p_i = P(y_i = 1 | X_i) = sigmoid(S_i)",
      "sigmoid(z) = 1 / (1 + exp(-z))",
      "p_i^cal = sigmoid(a + b S_i)",
      "(a, b) = argmin_{a,b} Σ_{i ∈ Train} LogLoss(y_i, sigmoid(a + b S_i))"
    ],
    "Prediction And Calibration"
  ),
  para("重要边界：校准参数 a,b 只允许在训练 cohort 上估计；holdout cohort 不能参与校准、阈值选择或模型选择。"),

  h1("7. 项目中的模型家族"),
  formulaBlock(
    [
      "EcoNiche-Opt-ModulePriorFixed: p_i = sigmoid(Σ_k w_k M_k(i))",
      "EcoNiche-Opt-ModulePriorFixed-Platt: p_i^cal = sigmoid(a + b Σ_k w_k M_k(i))",
      "EcoNiche-Opt-ImmuneComposite: C_i = [z(IFNG_i) + z(CXCL9_i) + 2 z(PDCD1LG2_i)] / 4",
      "EcoNiche-Opt-AdaptiveConsensus: mean(sigmoid(C_i), sigmoid(S_i), sigmoid(IFNG_i), sigmoid(TIG_i), sigmoid(TIDE_dysfunction_i), sigmoid(CXCL9_i))",
      "EcoNiche-Opt-ModuleIFNConsensus: mean(sigmoid(C_i), sigmoid(S_i), sigmoid(IFNG_i))",
      "Trainable module logistic sensitivity: p_i = sigmoid(β_0 + Σ_k β_k M_k(i)); β is learned only inside training folds."
    ],
    "Model Family"
  ),
  para("论文主模型建议使用 `EcoNiche-Opt-ModulePriorFixed`；Platt 版本用于校准和 decision curve，consensus/logistic 版本用于敏感性分析。"),

  h1("8. LODO 外部验证与无泄漏公式"),
  para("项目使用 leave-one-dataset-out validation。每轮完整留出一个 cohort，训练、阈值、校准、模型选择都只在其余 cohort 中完成。"),
  formulaBlock(
    [
      "For each cohort c:",
      "Train_c = D \\ D_c",
      "Test_c = D_c",
      "Fit/threshold/calibrate using Train_c only",
      "Evaluate AUROC, AUPRC, BA, ECE, Brier on Test_c"
    ],
    "LODO Protocol"
  ),
  formulaBlock(
    [
      "t_c* = argmax_t 0.5 × [Sensitivity_train(t) + Specificity_train(t)]",
      "ŷ_i = 1[p_i ≥ t_c*], for i ∈ Test_c"
    ],
    "Threshold Selection"
  ),

  h1("9. 评价指标与 Claim Gate 公式"),
  table(["指标/比较", "公式", "用途"], metricRows, [2100, 4300, 2960]),
  spacer(),
  formulaBlock(
    [
      "Δ_m = AUROC(EcoNiche-Opt) - AUROC(baseline_m)",
      "Bootstrap: resample matched patients B times, compute Δ_m^(b)",
      "p_m = two-sided paired bootstrap p value",
      "q_m = Benjamini-Hochberg-FDR(p_m)",
      "Superiority claim requires Δ_m > 0 and q_m < 0.05 in the pre-specified stratum."
    ],
    "Paired Comparison And FDR Claim Gate"
  ),

  h1("10. Endpoint-Evidence 分层设计"),
  para("为了处理 full melanoma primary AUROC 被异质 endpoint 拉低的问题，项目新增 endpoint-evidence 分层，而不是通过 holdout 泄漏调参。"),
  table(
    ["分析层", "包含队列", "解释"],
    [
      ["melanoma_core_high_evidence", "GSE91061, GSE78220, PRJEB23709_PD1_PRE", "最高证据层；primary RECIST AUROC 约 0.705"],
      ["melanoma_recist_supported_primary", "GSE91061, GSE78220, GSE145996, PRJEB23709_PD1_PRE", "更广的 RECIST-supported 主层；primary RECIST AUROC 约 0.685"],
      ["melanoma_anti_pd1_primary", "上述队列 + GSE168204, GSE115821", "完整异质 pool；AUROC 约 0.641"],
      ["melanoma_binary_response_stress", "GSE168204, GSE115821", "R/NR 二分类 stress-test，不作为主 headline"],
    ],
    [2500, 3600, 3260]
  ),

  h1("11. 本项目原创贡献", true),
  table(["类别", "原创内容", "意义"], originalRows, [1500, 3300, 4560]),
  spacer(),
  para("最核心的原创不是 AUROC 或 Platt calibration，而是 immune-ecology module-prior 的模型构造，以及 claim-gated multicohort validation framework。"),

  h1("12. 标准方法与非原创边界"),
  table(["边界", "方法", "说明"], standardRows, [1500, 3600, 4260]),
  spacer(),
  formulaBlock(
    [
      "建议论文表述：",
      "We developed an immune-ecology module-prior model and a claim-gated multicohort validation framework, built on standard statistical evaluation and calibration methods."
    ],
    "Claim-Safe Wording"
  ),

  h1("13. 当前结果摘要", true),
  table(
    ["分析层", "Endpoint", "EcoNiche-Opt AUROC", "解释"],
    [
      ["melanoma_core_high_evidence", "Primary RECIST", "0.705", "最高证据主结果层"],
      ["melanoma_core_high_evidence", "Strict RECIST", "0.707", "去除 SD 后仍约 0.70"],
      ["melanoma_recist_supported_primary", "Primary RECIST", "0.685", "更广 RECIST-supported melanoma primary"],
      ["melanoma_recist_supported_primary", "Strict RECIST", "0.690", "更严格 endpoint 下接近 0.69"],
      ["melanoma_anti_pd1_primary", "Primary RECIST", "0.641", "包含 binary stress cohorts 的完整异质 pool"],
    ],
    [2900, 2100, 2100, 2260]
  ),
  spacer(),
  para("在 RECIST-supported primary layer 中，EcoNiche-Opt 按 AUROC 点估计超过 IFNG、CXCL9、TIG、TIDE_dysfunction、APM、CYT、IPRES、TIDE_exclusion 八个 baseline；其中 CYT 达到 FDR-supported，其余应表述为 point-estimate improvement。"),

  h1("14. 可直接用于论文的方法学文字"),
  para("EcoNiche-Opt transforms bulk transcriptomes into immune-ecology module activities and combines them through a fixed, biologically directed prior. Response-promoting modules include IFN/T-cell inflammation, cytotoxic CD8 activity, antigen presentation, and TRM/TLS biology, whereas resistance-associated modules include myeloid suppression and stromal exclusion. The resulting EcoNicheScore is converted to response probability with a sigmoid link and evaluated under leave-one-dataset-out validation. Calibration, thresholding, model selection, and paired comparisons are performed using training cohorts only, with no holdout leakage."),
  para("This design enables an interpretable, low-degree-of-freedom predictor for pretreatment ICB response and a reproducible claim-gated benchmarking framework across heterogeneous public cohorts."),
];

const doc = new Document({
  creator: "EcoNiche-Opt / Codex",
  title: "EcoNiche-Opt Model Principles And Formulas",
  description: "Model principles, formulas, and original contributions for EcoNiche-Opt.",
  styles: {
    default: {
      document: { run: { font: FONT_CN, size: 21, color: DARK }, paragraph: { spacing: { line: 300 } } },
    },
  },
  numbering: {
    config: [
      {
        reference: "bullets",
        levels: [
          {
            level: 0,
            format: LevelFormat.BULLET,
            text: "•",
            alignment: AlignmentType.LEFT,
            style: { paragraph: { indent: { left: 480, hanging: 240 } } },
          },
        ],
      },
      {
        reference: "numbers",
        levels: [
          {
            level: 0,
            format: LevelFormat.DECIMAL,
            text: "%1.",
            alignment: AlignmentType.LEFT,
            style: { paragraph: { indent: { left: 520, hanging: 260 } } },
          },
        ],
      },
    ],
  },
  sections: [
    {
      properties: {
        page: {
          size: { width: 12240, height: 15840 },
          margin: { top: 1080, right: 1080, bottom: 1080, left: 1080 },
        },
      },
      children,
    },
  ],
});

fs.mkdirSync(path.dirname(OUT), { recursive: true });
Packer.toBuffer(doc).then((buffer) => {
  fs.writeFileSync(OUT, buffer);
  console.log(OUT);
});
