/* EcoNiche-Opt static portal: browser-side audit and locked-score preview. */
const state = {
  manifest: null,
  orientation: "samples",
  matrixText: "",
  matrix: null,
  scores: null,
  audit: null,
  moduleFilter: "all",
};

const $ = (selector, root = document) => root.querySelector(selector);
const $$ = (selector, root = document) => [...root.querySelectorAll(selector)];
const escapeHtml = (value) => String(value ?? "").replace(/[&<>"']/g, (char) => ({"&":"&amp;","<":"&lt;",">":"&gt;","\"":"&quot;","'":"&#039;"}[char]));
const pretty = (value) => String(value || "").replaceAll("_", " ").replace(/\b\w/g, (letter) => letter.toUpperCase());
const number = (value, digits = 3) => Number.isFinite(Number(value)) ? Number(value).toFixed(digits) : "--";

function showToast(message) {
  const toast = $("#toast");
  toast.textContent = message;
  toast.classList.add("is-visible");
  clearTimeout(showToast.timer);
  showToast.timer = setTimeout(() => toast.classList.remove("is-visible"), 3000);
}

function refreshIcons() {
  if (window.lucide) window.lucide.createIcons();
}

function setView(view) {
  $$("[data-view-panel]").forEach((panel) => panel.classList.toggle("is-visible", panel.dataset.viewPanel === view));
  $$(".nav-item").forEach((item) => item.classList.toggle("is-active", item.dataset.view === view));
  history.replaceState(null, "", `#${view}`);
  window.scrollTo({ top: 0, behavior: "smooth" });
  if (view === "evidence") renderEvidence();
  if (view === "panel") renderPanel();
  refreshIcons();
}

function renderOverview() {
  const { model, modules, endpoints, review } = state.manifest;
  const metrics = $("#overview-metrics").children;
  const values = [model.panel_unique_genes, modules.length, endpoints.length, `${review.resolved}/${review.total}`];
  [...metrics].forEach((card, index) => { $(".metric-value", card).textContent = values[index]; });
  $("#release-tag").textContent = model.release_tag;
  $("#build-stamp").textContent = model.release_tag.split("-")[0];
  $("#footer-release").textContent = model.release_tag;
  $("#release-digest").innerHTML = [
    ["Model", model.name],
    ["Intended context", model.context],
    ["Discovery source", model.threshold_source],
    ["Input contract", "samples-by-genes · sample_id row key"],
  ].map(([key, value]) => `<div><dt>${escapeHtml(key)}</dt><dd>${escapeHtml(value)}</dd></div>`).join("");
  $("#unsupported-context").textContent = model.unsupported_context;
  $("#score-context").textContent = model.name;
  $("#intended-context").textContent = model.context;
  $("#unsupported-context-repro").textContent = model.unsupported_context;
  $("#allowed-claim").textContent = model.allowed_claim;
  $("#forbidden-claim").textContent = model.forbidden_claim;
  $("#panel-count").textContent = model.panel_unique_genes;
  $("#audit-score").childNodes[0].textContent = `${state.manifest.audit.passed} / ${state.manifest.audit.total}`;
  refreshIcons();
}

function parseDelimited(text) {
  const lines = text.replace(/^\uFEFF/, "").split(/\r?\n/).filter((line) => line.trim() !== "");
  if (!lines.length) throw new Error("The file is empty.");
  const delimiter = lines[0].includes("\t") ? "\t" : ",";
  const parseLine = (line) => {
    const fields = [];
    let current = "";
    let quoted = false;
    for (let i = 0; i < line.length; i += 1) {
      const char = line[i];
      if (char === '"' && line[i + 1] === '"' && quoted) { current += '"'; i += 1; continue; }
      if (char === '"') { quoted = !quoted; continue; }
      if (char === delimiter && !quoted) { fields.push(current.trim()); current = ""; continue; }
      current += char;
    }
    fields.push(current.trim());
    return fields;
  };
  const headers = parseLine(lines[0]);
  const rows = lines.slice(1).map(parseLine);
  return { headers, rows, delimiter };
}

function parseExpression(text, orientation = state.orientation) {
  const parsed = parseDelimited(text);
  if (parsed.headers.length < 3) throw new Error("At least two numeric gene/sample columns are required.");
  if (orientation === "genes") {
    const geneNames = parsed.rows.map((row) => row[0]).filter(Boolean);
    const sampleNames = parsed.headers.slice(1);
    const rows = sampleNames.map((sample, sampleIndex) => {
      const values = {};
      parsed.rows.forEach((row, geneIndex) => { if (geneNames[geneIndex]) values[geneNames[geneIndex]] = Number(row[sampleIndex + 1]); });
      return { id: sample, values };
    });
    return { rows, genes: geneNames, orientation, invalid: rows.flatMap((row) => Object.values(row.values).filter((value) => !Number.isFinite(value))).length };
  }
  const genes = parsed.headers.slice(1).filter(Boolean);
  const rows = parsed.rows.map((row, index) => {
    const id = row[0] || `sample_${index + 1}`;
    const values = {};
    genes.forEach((gene, geneIndex) => { values[gene] = Number(row[geneIndex + 1]); });
    return { id, values };
  });
  return { rows, genes, orientation, invalid: rows.flatMap((row) => Object.values(row.values).filter((value) => !Number.isFinite(value))).length };
}

function uniquePanelGenes() { return [...new Set(state.manifest.modules.flatMap((module) => module.genes))]; }

function auditMatrix(matrix) {
  const panelGenes = uniquePanelGenes();
  const available = new Set(matrix.genes);
  const modules = state.manifest.modules.map((module) => {
    const covered = module.genes.filter((gene) => available.has(gene));
    return { ...module, covered, fraction: covered.length / module.genes.length };
  });
  const checks = [
    ["Expression matrix parsed", matrix.rows.length >= 2, `${matrix.rows.length} samples detected`],
    ["Numeric values complete", matrix.invalid === 0, matrix.invalid ? `${matrix.invalid} non-numeric cells` : "no non-numeric cells"],
    ["Sample identifiers unique", new Set(matrix.rows.map((row) => row.id)).size === matrix.rows.length, "sample IDs are row keys"],
    ["Locked panel represented", panelGenes.some((gene) => available.has(gene)), `${panelGenes.filter((gene) => available.has(gene)).length}/${panelGenes.length} genes present`],
    ["Pretreatment tumor context", $("#attest-tumor").checked, $("#attest-tumor").checked ? "user attested" : "attestation not selected"],
    ["Anti-PD-1 therapy context", $("#attest-therapy").checked, $("#attest-therapy").checked ? "user attested" : "attestation not selected"],
  ];
  return { checks, modules, panelGenes, available: panelGenes.filter((gene) => available.has(gene)), sampleCount: matrix.rows.length, ready: checks.every((check) => check[1]) };
}

function zscore(values) {
  const numeric = values.map((value) => Number.isFinite(value) ? value : 0);
  const mean = numeric.reduce((sum, value) => sum + value, 0) / numeric.length;
  const variance = numeric.reduce((sum, value) => sum + ((value - mean) ** 2), 0) / numeric.length;
  const sd = Math.sqrt(variance);
  return sd > 0 ? numeric.map((value) => (value - mean) / sd) : numeric.map(() => 0);
}

function sigmoid(value) { return 1 / (1 + Math.exp(-Math.max(-30, Math.min(30, value)))); }

function scoreMatrix(matrix) {
  const audit = auditMatrix(matrix);
  const moduleScores = Object.fromEntries(state.manifest.modules.map((module) => [module.id, []]));
  state.manifest.modules.forEach((module) => {
    const raw = matrix.rows.map((row) => {
      const values = module.genes.map((gene) => row.values[gene]).filter((value) => Number.isFinite(value));
      return values.length ? values.reduce((sum, value) => sum + value, 0) / values.length : 0;
    });
    moduleScores[module.id] = zscore(raw);
  });
  const rawScores = matrix.rows.map((_, rowIndex) => state.manifest.modules.reduce((sum, module) => sum + (module.weight * moduleScores[module.id][rowIndex]), 0));
  const rows = [];
  matrix.rows.forEach((row, index) => state.manifest.endpoints.forEach((endpoint) => {
    const probability = sigmoid((endpoint.calibration_coef * rawScores[index]) + endpoint.calibration_intercept);
    rows.push({ sample_id: row.id, endpoint: endpoint.endpoint, raw_module_prior_score: rawScores[index], response_probability: probability, locked_threshold: endpoint.threshold, predicted_label: probability >= endpoint.threshold ? 1 : 0, modules: Object.fromEntries(state.manifest.modules.map((module) => [module.id, moduleScores[module.id][index]])) });
  }));
  return { audit, rows, moduleScores, rawScores, matrix };
}

function renderAudit(audit) {
  $("#audit-meter-fill").style.width = `${(audit.checks.filter((check) => check[1]).length / audit.checks.length) * 100}%`;
  $("#audit-list").innerHTML = audit.checks.map(([name, valid, detail]) => `<div class="audit-row ${valid ? "" : "is-fail"}"><i data-lucide="${valid ? "check-circle-2" : "circle-x"}"></i><div><strong>${escapeHtml(name)}</strong><span>${escapeHtml(detail)}</span></div></div>`).join("");
  refreshIcons();
}

function renderResults(result) {
  const selectedEndpoint = $("#endpoint-select").value;
  const selected = result.rows.filter((row) => row.endpoint === selectedEndpoint);
  const positive = selected.filter((row) => row.predicted_label === 1).length;
  const selectedSpec = state.manifest.endpoints.find((endpoint) => endpoint.endpoint === selectedEndpoint);
  $("#output-status").textContent = result.audit.ready ? "audit ready" : "audit warnings";
  $("#empty-output").classList.add("is-hidden");
  $("#score-results").classList.remove("is-hidden");
  $("#score-summary").innerHTML = `<div class="summary-box primary"><span class="summary-label">${escapeHtml(pretty(selectedEndpoint))}</span><strong class="summary-value">${positive}/${selected.length}</strong><span class="summary-sub">samples called response · threshold ${number(selectedSpec.threshold, 3)}</span></div><div class="summary-box"><span class="summary-label">Mean probability</span><strong class="summary-value">${number(selected.reduce((sum, row) => sum + row.response_probability, 0) / selected.length, 3)}</strong><span class="summary-sub">label-free browser preview</span></div>`;
  const covered = result.audit.available.length;
  $("#coverage-strip").innerHTML = `<div class="coverage-chip ${covered < result.audit.panelGenes.length ? "warning" : ""}"><strong>${covered}/${result.audit.panelGenes.length}</strong><span>panel genes available</span></div><div class="coverage-chip"><strong>${result.matrix.rows.length}</strong><span>samples loaded</span></div><div class="coverage-chip"><strong>${result.audit.modules.filter((module) => module.fraction === 1).length}/${result.audit.modules.length}</strong><span>modules complete</span></div>`;
  const moduleMeans = state.manifest.modules.map((module) => ({ module, value: result.moduleScores[module.id].reduce((sum, value) => sum + value, 0) / result.matrix.rows.length }));
  const maxAbs = Math.max(...moduleMeans.map((item) => Math.abs(item.value)), 1);
  $("#module-chart").innerHTML = moduleMeans.map(({ module, value }) => `<div class="module-row"><span class="module-name" title="${escapeHtml(module.label)}">${escapeHtml(module.label)}</span><span class="bar-track"><span class="bar-fill ${value < 0 ? "negative" : ""}" style="width:${Math.max(3, Math.abs(value) / maxAbs * 100)}%"></span></span><span class="module-value">${value >= 0 ? "+" : ""}${number(value, 2)}</span></div>`).join("");
  $("#score-table thead").innerHTML = "<tr><th>Sample</th><th>Probability</th><th>Call</th><th>Raw score</th><th>Threshold</th></tr>";
  $("#score-table tbody").innerHTML = selected.map((row) => `<tr><td>${escapeHtml(row.sample_id)}</td><td class="probability">${number(row.response_probability, 4)}</td><td class="${row.predicted_label ? "call-response" : "call-nonresponse"}">${row.predicted_label ? "response" : "nonresponse"}</td><td>${number(row.raw_module_prior_score, 3)}</td><td>${number(row.locked_threshold, 3)}</td></tr>`).join("");
  renderAudit(result.audit);
  refreshIcons();
}

function downloadBlob(content, filename, type) { const blob = new Blob([content], { type }); const url = URL.createObjectURL(blob); const link = document.createElement("a"); link.href = url; link.download = filename; link.click(); setTimeout(() => URL.revokeObjectURL(url), 500); }

function downloadScores() {
  if (!state.scores) return;
  const header = ["sample_id", "endpoint", "raw_module_prior_score", "response_probability", "locked_threshold", "predicted_label", ...state.manifest.modules.map((module) => `module_score__${module.id}`)];
  const lines = [header.join("\t")];
  state.scores.rows.forEach((row) => lines.push([row.sample_id, row.endpoint, row.raw_module_prior_score, row.response_probability, row.locked_threshold, row.predicted_label, ...state.manifest.modules.map((module) => row.modules[module.id])].join("\t")));
  downloadBlob(`${lines.join("\n")}\n`, "econiche_opt_browser_scores.tsv", "text/tab-separated-values");
}

function downloadAudit() {
  if (!state.scores) return;
  downloadBlob(JSON.stringify({ release_tag: state.manifest.model.release_tag, model: state.manifest.model.name, audit: state.scores.audit, scores: state.scores.rows }, null, 2), "econiche_opt_browser_audit.json", "application/json");
}

function makeDemoText() {
  const genes = uniquePanelGenes();
  const header = ["sample_id", ...genes];
  const rows = ["DEMO_ONLY_01", "DEMO_ONLY_02", "DEMO_ONLY_03", "DEMO_ONLY_04", "DEMO_ONLY_05", "DEMO_ONLY_06"].map((sample, rowIndex) => [sample, ...genes.map((_, geneIndex) => (((rowIndex + 3) * (geneIndex + 5)) % 29) / 10)]);
  return [header, ...rows].map((row) => row.join("\t")).join("\n");
}

function loadMatrix(text, filename = "inline matrix") {
  try {
    const matrix = parseExpression(text, state.orientation);
    state.matrixText = text; state.matrix = matrix; state.scores = null;
    $("#selected-file").textContent = `${filename} · ${matrix.rows.length} samples · ${matrix.genes.length} columns`;
    $("#output-status").textContent = "matrix loaded";
    $("#empty-output").classList.remove("is-hidden"); $("#score-results").classList.add("is-hidden");
    showToast("Matrix loaded. Select the context attestations and run the audit.");
  } catch (error) { showToast(error.message); }
}

function readFile(file) { const reader = new FileReader(); reader.onload = () => loadMatrix(String(reader.result), file.name); reader.onerror = () => showToast("The file could not be read."); reader.readAsText(file); }

function renderEvidence() {
  if (!state.manifest) return;
  const { audit, review, benchmark } = state.manifest;
  $("#audit-score").childNodes[0].textContent = `${audit.passed} / ${audit.total}`;
  $("#audit-subtitle").textContent = `${audit.failed} failed · source audit`;
  $("#reviewer-subtitle").textContent = `${review.total} comments · resolution matrix`;
  $("#reviewer-badge").textContent = `${review.resolved}/${review.total} resolved`;
  $("#audit-meter-fill").style.width = `${audit.passed / audit.total * 100}%`;
  $("#audit-list").innerHTML = audit.rows.map((row) => `<div class="audit-row ${row.is_valid === "True" ? "" : "is-fail"}"><i data-lucide="${row.is_valid === "True" ? "check-circle-2" : "circle-x"}"></i><div><strong>${escapeHtml(row.check)}</strong><span>${escapeHtml(row.detail || "registered check")}</span></div></div>`).join("");
  renderReviewerTable();
  const preferred = benchmark.filter((row) => row.stratum === "melanoma_core_high_evidence").slice(0, 3);
  $("#benchmark-grid").innerHTML = preferred.length ? preferred.map((row) => `<div class="benchmark-item"><span>${escapeHtml(pretty(row.endpoint))}</span><strong>${number(row.pooled_AUROC, 3)}</strong><small>pooled AUROC · ${row.n_samples} samples · ${row.n_cohorts} cohorts</small><b>${escapeHtml(row.evaluation_modes || "LODO")}</b></div>`).join("") : `<div class="benchmark-item"><span>Registered output</span><strong>RESULT_PENDING</strong><small>No benchmark summary available in this bundle.</small></div>`;
  refreshIcons();
}

function renderReviewerTable() {
  const query = $("#review-search").value.trim().toLowerCase();
  const rows = state.manifest.review.rows.filter((row) => !query || Object.values(row).join(" ").toLowerCase().includes(query));
  $("#review-table tbody").innerHTML = rows.map((row) => `<tr><td>${escapeHtml(row.comment_id)}</td><td>${escapeHtml(row.issue)}</td><td>${escapeHtml(row.evidence)}</td><td>${escapeHtml(row.status)}</td></tr>`).join("");
}

function renderPanel() {
  if (!state.manifest) return;
  const query = $("#gene-search").value.trim().toLowerCase();
  const modules = state.manifest.modules.filter((module) => state.moduleFilter === "all" || module.id === state.moduleFilter);
  const rows = modules.flatMap((module) => module.genes.map((gene) => ({ gene, module }))).filter(({ gene, module }) => !query || `${gene} ${module.label}`.toLowerCase().includes(query));
  $("#gene-table tbody").innerHTML = rows.map(({ gene, module }) => `<tr><td>${escapeHtml(gene)}</td><td>${escapeHtml(module.label)}</td><td class="${module.weight < 0 ? "weight-negative" : "weight-positive"}">${module.weight > 0 ? "+" : ""}${module.weight.toFixed(2)}</td><td><span class="direction-pill">${escapeHtml(module.direction)}</span></td><td class="role-text">${module.weight < 0 ? "resistance / exclusion" : "response-supporting"}</td></tr>`).join("");
  $$(".filter-chip").forEach((chip) => chip.classList.toggle("is-active", chip.dataset.module === state.moduleFilter));
}

function initPanelFilters() {
  $("#module-filters").innerHTML = [`<button class="filter-chip is-active" type="button" data-module="all">All modules</button>`, ...state.manifest.modules.map((module) => `<button class="filter-chip" type="button" data-module="${escapeHtml(module.id)}">${escapeHtml(module.label)}</button>`)].join("");
  $$(".filter-chip").forEach((chip) => chip.addEventListener("click", () => { state.moduleFilter = chip.dataset.module; renderPanel(); }));
}

function init() {
  const storedTheme = localStorage.getItem("econiche-theme");
  if (storedTheme) document.documentElement.dataset.theme = storedTheme;
  $("#theme-toggle").addEventListener("click", () => { const next = document.documentElement.dataset.theme === "dark" ? "light" : "dark"; document.documentElement.dataset.theme = next; localStorage.setItem("econiche-theme", next); });
  $$(".nav-item").forEach((item) => item.addEventListener("click", () => setView(item.dataset.view)));
  $$("[data-go-view]").forEach((item) => item.addEventListener("click", () => setView(item.dataset.goView)));
  $$("[data-orientation]").forEach((button) => button.addEventListener("click", () => { state.orientation = button.dataset.orientation; $$("[data-orientation]").forEach((item) => item.classList.toggle("is-active", item === button)); if (state.matrixText) loadMatrix(state.matrixText, $("#selected-file").textContent.split(" · ")[0]); }));
  $("#choose-file").addEventListener("click", (event) => { event.stopPropagation(); $("#expression-file").click(); });
  $("#expression-file").addEventListener("change", (event) => { if (event.target.files[0]) readFile(event.target.files[0]); });
  $("#drop-zone").addEventListener("click", () => $("#expression-file").click());
  $("#drop-zone").addEventListener("keydown", (event) => { if (event.key === "Enter" || event.key === " ") { event.preventDefault(); $("#expression-file").click(); } });
  $("#drop-zone").addEventListener("dragover", (event) => { event.preventDefault(); $("#drop-zone").classList.add("is-dragging"); });
  $("#drop-zone").addEventListener("dragleave", () => $("#drop-zone").classList.remove("is-dragging"));
  $("#drop-zone").addEventListener("drop", (event) => { event.preventDefault(); $("#drop-zone").classList.remove("is-dragging"); if (event.dataTransfer.files[0]) readFile(event.dataTransfer.files[0]); });
  $("#load-demo").addEventListener("click", () => { loadMatrix(makeDemoText(), "deterministic demo fixture"); $("#attest-tumor").checked = true; $("#attest-therapy").checked = true; });
  $("#clear-input").addEventListener("click", () => { state.matrixText = ""; state.matrix = null; state.scores = null; $("#selected-file").textContent = "No file loaded"; $("#output-status").textContent = "waiting for input"; $("#empty-output").classList.remove("is-hidden"); $("#score-results").classList.add("is-hidden"); $("#expression-file").value = ""; });
  $("#run-score").addEventListener("click", () => { if (!state.matrix) { showToast("Load an expression matrix first."); return; } state.scores = scoreMatrix(state.matrix); renderResults(state.scores); if (!state.scores.audit.ready) showToast("Score generated with audit warnings. Review coverage and context before using the export."); });
  $("#download-csv").addEventListener("click", downloadScores); $("#download-audit").addEventListener("click", downloadAudit);
  $("#review-search").addEventListener("input", renderReviewerTable); $("#gene-search").addEventListener("input", renderPanel); $("#endpoint-select").addEventListener("change", () => { if (state.scores) renderResults(state.scores); });
  fetch("data/portal_manifest.json").then((response) => { if (!response.ok) throw new Error("Portal manifest unavailable"); return response.json(); }).then((manifest) => { state.manifest = manifest; renderOverview(); initPanelFilters(); renderPanel(); const hash = location.hash.replace("#", ""); setView(["overview", "score", "evidence", "panel", "reproducibility"].includes(hash) ? hash : "overview"); }).catch((error) => showToast(error.message));
  refreshIcons();
}

document.addEventListener("DOMContentLoaded", init);
