// ===================================================================
// LOOM CRISPR Search — Application Logic
// Loads pre-computed CRISPR targets and optionally a per-pathogen
// FM-index (via Web Worker + WASM) for live exact-match search.
// ===================================================================

// Index hosting: set to Azure Blob Storage URL after deployment.
// Run scripts/deploy-azure.sh to deploy, then paste the URL here.
// Local dev: use "indexes" to load from a local directory.
const INDEX_BASE_URL = "https://loomcrisprdata.blob.core.windows.net/indexes";

// Indexes available for browser live search.
// With Azure Blob Storage, all indexes are available (no size limit).
const BROWSER_INDEXES = new Set([
  "zika", "mers", "ebola", "hepatitis-b", "dengue", "rsv",
  "refseq-viral", "cholera", "tuberculosis",
  "hiv-1", "mpox", "influenza-a", "human-grch38"
]);

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------
const state = {
  data: null,              // Full JSON data keyed by pathogen
  selected: null,          // Currently selected pathogen key
  filteredTargets: [],     // Current filtered view
  filterText: "",
  filterPam: "all",
  filterGene: "all",       // Gene filter: "all", "protein-coding", or a symbol
  filterNovelty: "all",    // Novelty filter: "all", "novel", "published"
  page: 0,
  pageSize: 50,
  sortCol: 'o',           // Default sort column: coverage
  sortAsc: false,         // Default: descending
  // Live search
  worker: null,
  indexLoaded: false,
  indexPathogen: null,
  // PubMed scan
  pubmedScan: null,        // Pre-computed PubMed scan results
  // Diagnostic priority
  diagnosticPriority: null, // Diagnostic priority scores per pathogen
  // Guide scores
  guideScores: null,       // Pre-computed guide quality scores
  // Cross-reactivity
  crossReactivity: null,   // Cross-reactivity matrix (pathogen guides vs animal hosts)
  // Resistance overlay
  resistanceData: null,    // Drug-resistance variant overlay
  // Gene symbol maps (ontology bridge: product name → canonical symbol)
  _geneSymbolMaps: {},     // Cached per pathogen key
};

// ---------------------------------------------------------------------------
// DOM refs
// ---------------------------------------------------------------------------
const $ = (sel) => document.querySelector(sel);
const $$ = (sel) => document.querySelectorAll(sel);

const dom = {
  pathogenGrid:   $("#pathogenGrid"),
  filterInput:    $("#filterInput"),
  pamFilter:      $("#pamFilter"),
  resultsStatus:  $("#resultsStatus"),
  targetsBody:    $("#targetsBody"),
  emptyState:     $("#emptyState"),
  pagination:     $("#pagination"),
  exportCsv:      $("#exportCsv"),
  exportJson:     $("#exportJson"),
  totalTargets:   $("#totalTargets"),
  totalGenomes:   $("#totalGenomes"),
  // Tabs
  tabBar:         $("#tabBar"),
  databasePanel:  $("#databasePanel"),
  searchPanel:    $("#searchPanel"),
  // Live search
  indexBadge:     $("#indexBadge"),
  indexStatusText: $("#indexStatusText"),
  loadIndexBtn:   $("#loadIndexBtn"),
  loadProgress:   $("#loadProgress"),
  progressFill:   $("#progressFill"),
  progressText:   $("#progressText"),
  sequenceInput:  $("#sequenceInput"),
  sequenceInfo:   $("#sequenceInfo"),
  liveResults:    $("#liveResults"),
  // Literature scanner
  literaturePanel:  $("#literaturePanel"),
  litSequenceInput: $("#litSequenceInput"),
  litSearchBtn:     $("#litSearchBtn"),
  litScanSelected:  $("#litScanSelected"),
  litScanNGG:       $("#litScanNGG"),
  litScanPathogen:  $("#litScanPathogen"),
  litStatus:        $("#litStatus"),
  litProgressFill:  $("#litProgressFill"),
  litProgressText:  $("#litProgressText"),
  litResults:       $("#litResults"),
  // Disease context
  contextPanel:     $("#contextPanel"),
  contextContent:   $("#contextContent"),
  // Data sources
  sourcesPanel:     $("#sourcesPanel"),
  copyCitation:     $("#copyCitation"),
  citeBlock:        $("#citeBlock"),
  // Gene filter & top opportunities
  geneFilter:       $("#geneFilter"),
  noveltyFilter:    $("#noveltyFilter"),
  topOpportunities: $("#topOpportunities"),
};

// ---------------------------------------------------------------------------
// Init
// ---------------------------------------------------------------------------
async function init() {
  bindEvents();
  await loadData();
  initPanelDesigner();
  await Promise.all([
    loadPubmedScan(),
    loadOntologyEnrichment(),
    loadDiagnosticPriority(),
    loadGuideScores(),
    loadCrossReactivity(),
    loadResistanceData(),
  ]);
}

async function loadData() {
  const isAnimal = document.body.dataset.pageType === "animal";
  const jsonFile = isAnimal ? "data/animal-targets.json" : "data/crispr-targets.json";
  try {
    const resp = await fetch(jsonFile);
    if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
    state.data = await resp.json();
    renderPathogenGrid();
    updateHeaderStats();
  } catch (err) {
    dom.resultsStatus.innerHTML =
      `<span class="status-text" style="color:var(--red)">Failed to load data: ${escapeHtml(err.message)}</span>`;
  }
}

async function loadPubmedScan() {
  try {
    const resp = await fetch("data/pubmed-scan-results.json");
    if (!resp.ok) return; // optional file, silently skip
    state.pubmedScan = await resp.json();
    renderPubmedBanner();
  } catch { /* scan results not available */ }
}

async function loadGuideScores() {
  try {
    const resp = await fetch("data/guide-scores.json");
    if (!resp.ok) return;
    state.guideScores = await resp.json();
  } catch { /* optional */ }
}

async function loadCrossReactivity() {
  try {
    const resp = await fetch("data/cross-reactivity.json");
    if (!resp.ok) return;
    state.crossReactivity = await resp.json();
  } catch { /* optional */ }
}

async function loadResistanceData() {
  try {
    const resp = await fetch("data/clinvar-resistance.json");
    if (!resp.ok) return;
    state.resistanceData = await resp.json();
  } catch { /* optional */ }
}

function getTopTargets(pathogenKey, count = 5) {
  if (!state.data || !state.data[pathogenKey]) return [];
  const targets = state.data[pathogenKey].targets || [];
  // Sort: NGG first, then by genome coverage descending
  return [...targets]
    .sort((a, b) => {
      const aNGG = a.m === 'NGG' ? 1 : 0;
      const bNGG = b.m === 'NGG' ? 1 : 0;
      if (bNGG !== aNGG) return bNGG - aNGG;
      return b.o - a.o;
    })
    .slice(0, count);
}

function renderSeqRows(targets, expanded, totalGenomes) {
  const show = expanded ? targets : targets.slice(0, 5);
  let html = '';
  for (const t of show) {
    const covPct = totalGenomes > 0 ? (Math.min(t.o, totalGenomes) / totalGenomes * 100).toFixed(1) : '—';
    const covStr = totalGenomes > 0
      ? `${covPct}% (${t.o.toLocaleString()}/${totalGenomes.toLocaleString()})`
      : `${t.o.toLocaleString()} genomes`;
    html += `<div class="lit-seq-row">
      <span class="dna-seq">${colorDNA(t.s)}</span>
      ${pamBadge(t.m)}
      <span class="lit-seq-meta">${covStr}</span>
      <span class="lit-seq-meta">pos ${t.p.toLocaleString()}</span>
      <button class="copy-seq-btn" onclick="copySeq('${t.s}')" title="Copy sequence">Copy</button>
    </div>`;
  }
  return html;
}

function renderPubmedBanner() {
  if (!state.pubmedScan || !dom.litResults) return;
  const scan = state.pubmedScan;
  const s = scan.summary;
  const pats = scan.pathogens || {};
  const opps = scan.opportunities || [];

  let html = `
    <div class="lit-summary">
      <h4>PubMed Literature Landscape</h4>
      <p style="font-size:0.85rem;color:var(--text-2);margin-bottom:var(--sp-3)">
        Generated ${new Date(scan.generated).toLocaleDateString()} — region-level scan across ${s.total_pathogens} pathogens.
        Identifies under-studied gene targets and application areas where our ${s.total_targets_in_db?.toLocaleString()} pre-computed guides could enable new research.
      </p>
      <div class="lit-stats-row">
        <div class="lit-stat"><span class="lit-stat-val">${(s.total_targets_in_db || 0).toLocaleString()}</span><span class="lit-stat-lbl">CRISPR Targets</span></div>
        <div class="lit-stat"><span class="lit-stat-val">${(s.ngg_targets || 0).toLocaleString()}</span><span class="lit-stat-lbl">NGG/SpCas9</span></div>
        <div class="lit-stat novel"><span class="lit-stat-val">${s.gene_gaps || 0}</span><span class="lit-stat-lbl">Gene Gaps</span></div>
        <div class="lit-stat novel"><span class="lit-stat-val">${s.application_gaps || 0}</span><span class="lit-stat-lbl">Application Gaps</span></div>
      </div>
    </div>`;

  // Research Opportunities — sequences-first
  if (opps.length > 0) {
    // Group opportunities by pathogen
    const byPathogen = {};
    for (const opp of opps) {
      if (!byPathogen[opp.pathogen_key]) byPathogen[opp.pathogen_key] = [];
      byPathogen[opp.pathogen_key].push(opp);
    }

    html += `
    <div style="margin-top:var(--sp-4)">
      <h4 style="font-family:var(--font-display)">Research Opportunities <span class="novelty-badge novel">${opps.length} gaps found</span></h4>
      <p style="font-size:0.82rem;color:var(--text-3);margin-bottom:var(--sp-3)">
        Under-studied pathogen + gene/application combinations with zero PubMed results. Top guide sequences shown for each — click <strong>Copy</strong> to grab a sequence.
      </p>`;

    for (const [key, pathogenOpps] of Object.entries(byPathogen)) {
      const pdata = pats[key] || {};
      const topTargets = getTopTargets(key, 20);
      const nggTargets = topTargets.filter(t => t.m === 'NGG');
      const showTargets = nggTargets.length >= 3 ? nggTargets.slice(0, 10) : topTargets.slice(0, 10);

      html += `<div class="lit-opp-card" id="opp-${key}">
        <div class="lit-opp-header">
          <div>
            <div class="paper-title">${escapeHtml(pdata.name || key)}</div>
            <div class="paper-meta">
              <span><strong>${(pdata.total_papers || 0).toLocaleString()}</strong> CRISPR papers</span>
              <span><strong>${(pdata.targets_in_db || 0).toLocaleString()}</strong> targets in DB</span>
            </div>
          </div>
          <span class="novelty-badge novel">${pathogenOpps.length} gap${pathogenOpps.length > 1 ? 's' : ''}</span>
        </div>
        <div class="lit-opp-gaps">`;

      for (const opp of pathogenOpps) {
        const icon = opp.type === 'gene_gap' ? '🧬' : '🔬';
        html += `<span class="lit-gap-tag">${icon} ${escapeHtml(opp.type === 'gene_gap' ? opp.gene : opp.area)}</span>`;
      }

      html += `</div>
        <div class="lit-seq-label">Top guide sequences (NGG/SpCas9 prioritized, by genome coverage):</div>
        <div class="lit-seq-list">
          ${renderSeqRows(showTargets, false, state.data[key]?.genomes || 0)}
        </div>`;

      if (topTargets.length > 5) {
        html += `<button class="lit-show-more" onclick="toggleMoreSeqs(this, '${key}')">Show ${Math.min(topTargets.length, 20) - 5} more sequences</button>`;
      }

      html += `</div>`;
    }

    html += `</div>`;
  }

  // Well-studied pathogens (no gaps) — compact summary + sequences
  const wellStudied = Object.entries(pats).filter(([k, d]) =>
    !(d.gap_applications?.length || d.unstudied_genes?.length)
  );
  if (wellStudied.length > 0) {
    html += `
    <div style="margin-top:var(--sp-5)">
      <h4 style="font-family:var(--font-display)">Well-Studied Pathogens</h4>
      <p style="font-size:0.82rem;color:var(--text-3);margin-bottom:var(--sp-3)">All gene regions and application areas have existing CRISPR literature.</p>`;

    for (const [key, data] of wellStudied) {
      const topTargets = getTopTargets(key, 5);
      html += `<div class="lit-opp-card">
        <div class="lit-opp-header">
          <div>
            <div class="paper-title">${escapeHtml(data.name)}</div>
            <div class="paper-meta">
              <span><strong>${(data.total_papers || 0).toLocaleString()}</strong> CRISPR papers</span>
              <span class="novelty-badge published">well-studied</span>
            </div>
          </div>
        </div>
        ${data.studied_genes?.length ? `<div style="font-size:0.78rem;color:var(--text-3)">Studied genes: ${data.studied_genes.join(', ')}</div>` : ''}
        <div class="lit-seq-label" style="margin-top:var(--sp-2)">Top guide sequences:</div>
        <div class="lit-seq-list">${renderSeqRows(topTargets, true, state.data[key]?.genomes || 0)}</div>
      </div>`;
    }
    html += `</div>`;
  }

  // Application coverage heatmap — collapsed
  const apps = scan.applications || {};
  if (Object.keys(apps).length > 0) {
    html += `
    <details style="margin-top:var(--sp-5)">
      <summary style="font-family:var(--font-display);font-weight:600;cursor:pointer;padding:var(--sp-2) 0">Application Coverage Details</summary>
      <p style="font-size:0.82rem;color:var(--text-3);margin-bottom:var(--sp-3)">Papers per pathogen + application combination.</p>
      <div class="paper-list">`;

    for (const [key, appData] of Object.entries(apps)) {
      const pname = pats[key]?.name || key;
      const sorted = Object.entries(appData).sort((a, b) => b[1].count - a[1].count);
      html += `<div class="paper-card"><div class="paper-title">${escapeHtml(pname)}</div><div style="display:flex;flex-wrap:wrap;gap:var(--sp-2);margin-top:var(--sp-2)">`;
      for (const [app, r] of sorted) {
        const cls = r.count > 0 ? 'published' : 'novel';
        html += `<span class="novelty-badge ${cls}" style="font-size:0.72rem">${escapeHtml(app)}: ${r.count}</span>`;
      }
      html += `</div></div>`;
    }
    html += `</div></details>`;
  }

  dom.litResults.innerHTML = html;
}

function copySeq(seq) {
  navigator.clipboard.writeText(seq).then(() => {
    // Brief visual feedback
    const btns = document.querySelectorAll('.copy-seq-btn');
    for (const b of btns) {
      if (b.closest('.lit-seq-row')?.querySelector('.dna-seq')?.textContent === seq) {
        b.textContent = 'Copied!';
        b.classList.add('copied');
        setTimeout(() => { b.textContent = 'Copy'; b.classList.remove('copied'); }, 1500);
        break;
      }
    }
  });
}
window.copySeq = copySeq;

function toggleMoreSeqs(btn, pathogenKey) {
  const card = btn.closest('.lit-opp-card');
  const list = card.querySelector('.lit-seq-list');
  const topTargets = getTopTargets(pathogenKey, 20);
  const nggTargets = topTargets.filter(t => t.m === 'NGG');
  const showTargets = nggTargets.length >= 3 ? nggTargets.slice(0, 20) : topTargets.slice(0, 20);

  if (btn.dataset.expanded === 'true') {
    list.innerHTML = renderSeqRows(showTargets.slice(0, 5), false, state.data[pathogenKey]?.genomes || 0);
    btn.textContent = `Show ${showTargets.length - 5} more sequences`;
    btn.dataset.expanded = 'false';
  } else {
    list.innerHTML = renderSeqRows(showTargets, true, state.data[pathogenKey]?.genomes || 0);
    btn.textContent = 'Show fewer';
    btn.dataset.expanded = 'true';
  }
}
window.toggleMoreSeqs = toggleMoreSeqs;

// ---------------------------------------------------------------------------
// Tab Switching
// ---------------------------------------------------------------------------
const panelMap = { database: "databasePanel", search: "searchPanel", literature: "literaturePanel", context: "contextPanel", sources: "sourcesPanel", panel: "panelPanel" };

function switchTab(tab) {
  dom.tabBar.querySelectorAll(".tab-btn").forEach((b) => b.classList.toggle("active", b.dataset.tab === tab));
  document.querySelectorAll(".tab-panel").forEach((p) => p.classList.remove("active"));
  document.getElementById(panelMap[tab] || "databasePanel").classList.add("active");
  history.replaceState(null, '', '#' + tab);
}

function viewTargetsFor(pathogenKey) {
  if (!state.data || !state.data[pathogenKey]) return;
  selectPathogen(pathogenKey);
  switchTab('database');
  // Scroll to results
  dom.pathogenGrid.querySelector(`.pathogen-card[data-key="${pathogenKey}"]`)?.scrollIntoView({ behavior: 'smooth', block: 'center' });
}
// Expose for inline onclick in lit scanner (module scope)
window.viewTargetsFor = viewTargetsFor;

// ---------------------------------------------------------------------------
// Header Stats
// ---------------------------------------------------------------------------
function updateHeaderStats() {
  if (!state.data) return;
  let targets = 0, genomes = 0;
  for (const p of Object.values(state.data)) {
    targets += p.target_count;
    genomes = Math.max(genomes, p.genomes); // biggest set for headline
  }
  dom.totalTargets.textContent = targets.toLocaleString();
  // Sum unique genomes across all pathogens
  const totalG = Object.values(state.data).reduce((s, p) => s + p.genomes, 0);
  dom.totalGenomes.textContent = totalG >= 1000 ? `${(totalG / 1000).toFixed(0)}K+` : totalG.toLocaleString();
}

// ---------------------------------------------------------------------------
// Pathogen Grid
// ---------------------------------------------------------------------------
function renderPathogenGrid() {
  if (!state.data) return;
  const keys = Object.keys(state.data);
  const prio = state.diagnosticPriority?.pathogens;
  dom.pathogenGrid.innerHTML = keys.map((key) => {
    const p = state.data[key];
    let badgeHtml = '';
    if (prio && prio[key]) {
      const score = prio[key].total_score;
      const cls = score >= 70 ? 'priority-high' : score >= 40 ? 'priority-med' : 'priority-low';
      const label = score >= 70 ? 'High' : score >= 40 ? 'Medium' : 'Low';
      badgeHtml = `<span class="priority-badge ${cls}" title="Diagnostic Priority: ${score}/100 — composite of CRISPR dx gap, disease burden, case fatality, WHO class, outbreak potential">Dx ${label}</span>`;
    }
    return `
      <div class="pathogen-card" data-key="${key}">
        ${badgeHtml}
        <div class="pathogen-name">${escapeHtml(p.name)}</div>
        <div class="pathogen-family">${escapeHtml(p.family)}</div>
        <div class="pathogen-meta">
          <span><strong>${p.target_count.toLocaleString()}</strong> targets</span>
          <span><strong>${p.genomes.toLocaleString()}</strong> genomes</span>
        </div>
      </div>
    `;
  }).join("");

  // Click handlers
  dom.pathogenGrid.querySelectorAll(".pathogen-card").forEach((card) => {
    card.addEventListener("click", () => selectPathogen(card.dataset.key));
  });
}

function selectPathogen(key) {
  state.selected = key;
  state.page = 0;
  state.filterText = "";
  state.filterPam = "all";
  state.filterGene = "all";
  state.filterNovelty = "all";
  dom.filterInput.value = "";
  if (dom.noveltyFilter) dom.noveltyFilter.value = "all";

  // Populate gene filter dropdown
  if (dom.geneFilter) {
    const targets = state.data[key]?.targets || [];
    const geneSet = new Set();
    for (const t of targets) {
      if (t.g) geneSet.add(t.g);
    }
    const genes = [...geneSet].sort();
    dom.geneFilter.innerHTML = '<option value="all">All genes</option><option value="protein-coding">Protein-coding only</option>';
    for (const g of genes) {
      dom.geneFilter.innerHTML += `<option value="${escapeHtml(g)}">${escapeHtml(g)}</option>`;
    }
    dom.geneFilter.value = "all";
  }

  // Highlight card
  dom.pathogenGrid.querySelectorAll(".pathogen-card").forEach((c) =>
    c.classList.toggle("selected", c.dataset.key === key)
  );

  // Reset PAM chips
  dom.pamFilter.querySelectorAll(".chip").forEach((c) =>
    c.classList.toggle("active", c.dataset.pam === "all")
  );

  // Enable export buttons
  dom.exportCsv.disabled = false;
  dom.exportJson.disabled = false;

  // Update coverage column header with total genomes
  const coverageH = document.getElementById('coverageHeader');
  if (coverageH) coverageH.textContent = `Coverage (${state.data[key].genomes.toLocaleString()} genomes)`;

  // Enable index loading (only for indexes that fit browser download)
  if (BROWSER_INDEXES.has(key)) {
    dom.loadIndexBtn.disabled = false;
    dom.loadIndexBtn.textContent = `Load ${state.data[key].name} Index (${state.data[key].index_mb} MB)`;
    dom.loadIndexBtn.title = "";
  } else {
    dom.loadIndexBtn.disabled = true;
    dom.loadIndexBtn.textContent = `Index too large for browser (${state.data[key].index_mb} MB)`;
    dom.loadIndexBtn.title = "This index exceeds 2 GB. Use the CLI for live search on this pathogen.";
  }

  // Enable literature quick-scan buttons
  dom.litScanSelected.disabled = false;
  dom.litScanNGG.disabled = false;
  dom.litScanPathogen.textContent = state.data[key].name;

  // Render disease context if ontology data available
  renderDiseaseContext(key);

  applyFilters();
}

// ---------------------------------------------------------------------------
// Filters
// ---------------------------------------------------------------------------
function applyFilters() {
  if (!state.selected || !state.data[state.selected]) return;
  const pathogen = state.data[state.selected];
  let targets = pathogen.targets;

  // PAM filter
  if (state.filterPam !== "all") {
    targets = targets.filter((t) => {
      if (state.filterPam === "none") return t.m === "none";
      return t.m === state.filterPam;
    });
  }

  // Gene filter
  if (state.filterGene === "protein-coding") {
    targets = targets.filter((t) => t.gt === "protein-coding");
  } else if (state.filterGene !== "all") {
    targets = targets.filter((t) => t.g === state.filterGene);
  }

  // Novelty filter (cross-references PubMed scan data via ontology bridge)
  if (state.filterNovelty !== "all" && state.pubmedScan) {
    const lookup = getNoveltyLookup(state.selected);
    if (lookup) {
      if (state.filterNovelty === "novel") {
        // "Confirmed gap" = gene IS in unstudied_genes (4-strategy v4 verified with ontology, ~8 genes)
        targets = targets.filter((t) => {
          const sym = resolveGeneSymbol(t.g, state.selected);
          return sym && lookup.unstudied.has(sym);
        });
      } else if (state.filterNovelty === "published") {
        // "Published" = gene IS in studied_genes
        targets = targets.filter((t) => {
          const sym = resolveGeneSymbol(t.g, state.selected);
          return sym && lookup.studied.has(sym);
        });
      }
    }
  }

  // Text filter (sequence match)
  if (state.filterText) {
    const q = state.filterText.toUpperCase();
    targets = targets.filter((t) => t.s.includes(q));
  }

  state.filteredTargets = targets;
  sortTargets();
  state.page = 0;
  renderResults();
}

function sortTargets() {
  const col = state.sortCol;
  const dir = state.sortAsc ? 1 : -1;
  const scores = state.guideScores?.pathogens?.[state.selected] || {};
  state.filteredTargets.sort((a, b) => {
    if (col === 's') return dir * a.s.localeCompare(b.s);
    if (col === 'g') return dir * (a.g || '').localeCompare(b.g || '');
    if (col === 'p') return dir * (a.p - b.p);
    if (col === 'm') return dir * a.m.localeCompare(b.m);
    if (col === 'o') return dir * (a.o - b.o);
    if (col === 'sc') return dir * ((scores[a.s]?.score || 0) - (scores[b.s]?.score || 0));
    return 0;
  });
}

// ---------------------------------------------------------------------------
// Render Targets Table
// ---------------------------------------------------------------------------
function renderResults() {
  const targets = state.filteredTargets;
  const pathogen = state.data[state.selected];
  const total = pathogen.targets.length;
  const filtered = targets.length;

  // Status
  const pamInfo = state.filterPam !== "all" ? ` | PAM: ${state.filterPam}` : "";
  const filterInfo = state.filterText ? ` | Filter: "${escapeHtml(state.filterText)}"` : "";
  const noveltyInfo = state.filterNovelty !== "all" ? ` | ${state.filterNovelty === "novel" ? "Novel only" : "Published only"}` : "";
  dom.resultsStatus.innerHTML =
    `<span class="status-text"><strong>${escapeHtml(pathogen.name)}</strong> — ` +
    `<span class="count">${filtered.toLocaleString()}</span> of ${total.toLocaleString()} targets${pamInfo}${noveltyInfo}${filterInfo}</span>`;

  // Paginate
  const start = state.page * state.pageSize;
  const pageTargets = targets.slice(start, start + state.pageSize);

  // Render rows
  if (pageTargets.length === 0) {
    dom.targetsBody.innerHTML = "";
    dom.emptyState.querySelector("p").textContent =
      state.selected ? "No targets match your filter" : "Select a pathogen above to view its CRISPR targets";
    dom.emptyState.style.display = "flex";
  } else {
    dom.emptyState.style.display = "none";
    const totalG = pathogen.genomes;
    // Build novelty lookup from PubMed scan data (case-insensitive via ontology bridge)
    const noveltyLookup = getNoveltyLookup(state.selected);
    // Guide scores and cross-reactivity lookups
    const scores = state.guideScores?.pathogens?.[state.selected] || {};
    const xr = state.crossReactivity?.pathogens?.[state.selected] || [];
    const xrMap = {};
    for (const g of xr) { xrMap[g.seq] = g; }
    // Resistance overlay
    const resMuts = state.resistanceData?.pathogens?.[state.selected]?.mutations || [];
    dom.targetsBody.innerHTML = pageTargets.map((t, i) => {
      const pct = totalG > 0 ? (Math.min(t.o, totalG) / totalG * 100) : 0;
      const barW = Math.max(1, Math.round(pct / 100 * 48));
      const geneLabel = t.g || 'intergenic';
      const geneTitle = t.gn ? escapeHtml(t.gn) : '';
      const geneCls = t.gt === 'protein-coding' ? 'gene-tag gene-pc' : t.g ? 'gene-tag gene-other' : 'gene-tag gene-inter';
      let noveltyBadge = '';
      if (noveltyLookup && t.g) {
        const sym = resolveGeneSymbol(t.g, state.selected);
        if (sym && noveltyLookup.unstudied.has(sym)) {
          noveltyBadge = ' <span class="novelty-badge novel">research gap</span>';
        } else if (sym && noveltyLookup.studied.has(sym)) {
          noveltyBadge = ' <span class="novelty-badge published">studied</span>';
        }
      }
      // Guide quality score
      const sc = scores[t.s];
      const scoreVal = sc ? sc.score : null;
      const scoreCls = scoreVal >= 70 ? 'score-high' : scoreVal >= 40 ? 'score-med' : 'score-low';
      const scoreHtml = scoreVal != null
        ? `<span class="guide-score ${scoreCls}" title="GC:${sc.gc}% Seed:${sc.seed_gc}% Homo:${sc.max_homopolymer}">${scoreVal}</span>`
        : '<span class="guide-score score-na">\u2014</span>';
      // Cross-reactivity badge
      const xrInfo = xrMap[t.s];
      let xrBadge = '';
      if (xrInfo) {
        xrBadge = xrInfo.specific
          ? ' <span class="xr-badge xr-specific" title="No hits in any host genome">specific</span>'
          : ` <span class="xr-badge xr-reactive" title="Hits: ${Object.entries(xrInfo.host_hits).filter(([,v])=>v>0).map(([k,v])=>`${k}:${v}`).join(', ')}">${xrInfo.total_host_hits} host hits</span>`;
      }
      // Resistance region overlap
      let resBadge = '';
      const tEnd = t.p + 23;
      for (const m of resMuts) {
        const [ms, me] = m.region.split('-').map(Number);
        if (t.p < me && tEnd > ms) {
          resBadge = ` <span class="res-badge" title="${escapeHtml(m.drug)}: ${escapeHtml(m.name)}">\u2622 ${escapeHtml(m.gene)}</span>`;
          break;
        }
      }
      // BLAST link
      const blastUrl = `https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&PROGRAM=blastn&DATABASE=nt&QUERY=${encodeURIComponent(t.s)}`;
      return `
      <tr>
        <td class="col-num">${(start + i + 1).toLocaleString()}</td>
        <td class="col-seq"><span class="dna-seq">${colorDNA(t.s)}</span>${xrBadge}${resBadge}</td>
        <td class="col-gene"><span class="${geneCls}" title="${geneTitle}" data-gene="${escapeHtml(t.g || '')}">${escapeHtml(geneLabel)}</span>${noveltyBadge}</td>
        <td class="col-pos">${t.p.toLocaleString()}</td>
        <td class="col-pam">${pamBadge(t.m)}</td>
        <td class="col-occ">${t.o.toLocaleString()}<span class="coverage-pct">${pct.toFixed(1)}%</span> <span class="coverage-bar" style="width:${barW}px"></span></td>
        <td class="col-score">${scoreHtml}</td>
        <td class="col-actions"><a href="${blastUrl}" target="_blank" rel="noopener" class="blast-link" title="BLAST this sequence against NCBI nt">BLAST</a></td>
      </tr>`;
    }).join("");

    // Gene tag click → switch to Disease Context + Gene Map
    dom.targetsBody.querySelectorAll('.gene-tag[data-gene]').forEach(tag => {
      if (!tag.dataset.gene) return;
      tag.style.cursor = 'pointer';
      tag.addEventListener('click', (e) => {
        e.stopPropagation();
        const gene = tag.dataset.gene;
        if (!gene) return;
        switchTab('context');
        // Open Gene Map details and try to highlight
        setTimeout(() => {
          const details = dom.contextContent?.querySelector('details:has(.gene-table)');
          if (details) details.open = true;
          const rows = dom.contextContent?.querySelectorAll('.gene-table tbody tr');
          if (rows) {
            for (const row of rows) {
              const sym = row.querySelector('.gene-symbol');
              if (sym && sym.textContent.trim() === gene) {
                row.classList.add('gene-highlight');
                row.scrollIntoView({ behavior: 'smooth', block: 'center' });
                break;
              }
            }
          }
        }, 100);
      });
    });
  }

  renderPagination(filtered);
}

function renderPagination(total) {
  const pages = Math.ceil(total / state.pageSize);
  if (pages <= 1) {
    dom.pagination.innerHTML = "";
    return;
  }

  const cur = state.page;
  let html = "";

  // Prev
  html += `<button class="page-btn" data-page="${cur - 1}" ${cur === 0 ? "disabled" : ""}>&larr;</button>`;

  // Page numbers (show up to 7 with ellipsis)
  const range = paginationRange(cur, pages, 7);
  for (const p of range) {
    if (p === "...") {
      html += `<span class="page-info">&hellip;</span>`;
    } else {
      html += `<button class="page-btn ${p === cur ? "active" : ""}" data-page="${p}">${p + 1}</button>`;
    }
  }

  // Next
  html += `<button class="page-btn" data-page="${cur + 1}" ${cur >= pages - 1 ? "disabled" : ""}>&rarr;</button>`;

  // Info
  html += `<span class="page-info">${(cur * state.pageSize + 1).toLocaleString()}–${Math.min((cur + 1) * state.pageSize, total).toLocaleString()} of ${total.toLocaleString()}</span>`;

  dom.pagination.innerHTML = html;

  dom.pagination.querySelectorAll(".page-btn").forEach((btn) => {
    btn.addEventListener("click", () => {
      if (btn.disabled) return;
      state.page = parseInt(btn.dataset.page, 10);
      renderResults();
      dom.targetsBody.scrollIntoView({ behavior: "smooth", block: "start" });
    });
  });
}

function paginationRange(current, total, maxVisible) {
  if (total <= maxVisible) return Array.from({ length: total }, (_, i) => i);
  const half = Math.floor(maxVisible / 2);
  let start = Math.max(0, current - half);
  let end = start + maxVisible - 1;
  if (end >= total) {
    end = total - 1;
    start = end - maxVisible + 1;
  }
  const range = [];
  for (let i = start; i <= end; i++) range.push(i);
  if (range[0] > 0) { range[0] = 0; range[1] = "..."; }
  if (range[range.length - 1] < total - 1) { range[range.length - 1] = total - 1; range[range.length - 2] = "..."; }
  return range;
}

// ---------------------------------------------------------------------------
// Web Worker / WASM Live Search
// ---------------------------------------------------------------------------
function ensureWorker() {
  if (state.worker) return;
  state.worker = new Worker("loom-worker.js", { type: "module" });
  state.worker.onmessage = handleWorkerMessage;
  state.worker.onerror = (e) => {
    console.error("Worker error:", e);
    setIndexStatus("error", `Worker error: ${e.message}`);
  };
}

function handleWorkerMessage(e) {
  const msg = e.data;
  if (msg.type === "loaded") {
    state.indexLoaded = true;
    state.indexPathogen = state.selected;
    setIndexStatus("ready",
      `${state.data[state.selected].name} — ${msg.files.toLocaleString()} files, ${formatBytes(msg.corpus)} corpus`);
    dom.loadProgress.hidden = true;
    dom.sequenceInput.disabled = false;
    dom.sequenceInput.focus();
  } else if (msg.type === "error") {
    setIndexStatus("error", `Failed to load: ${msg.error}`);
    dom.loadProgress.hidden = true;
  } else if (msg.type === "results") {
    renderLiveResults(msg);
  }
}

function loadIndex() {
  if (!state.selected) return;
  ensureWorker();
  state.indexLoaded = false;
  const key = state.selected;
  const url = `${INDEX_BASE_URL}/${key}.idx`;

  setIndexStatus("loading", `Loading ${state.data[key].name} index...`);
  dom.loadProgress.hidden = false;
  dom.sequenceInput.disabled = true;

  // Simulate progress (we can't get real fetch progress from worker)
  let fakeProgress = 0;
  const interval = setInterval(() => {
    fakeProgress = Math.min(fakeProgress + Math.random() * 8, 90);
    dom.progressFill.style.width = `${fakeProgress}%`;
    dom.progressText.textContent = `Loading ${state.data[key].name}... ${Math.round(fakeProgress)}%`;
  }, 300);

  // Store interval to clear on load complete
  state._progressInterval = interval;

  state.worker.postMessage({ type: "load", url, id: 0 });

  // Override message handler to clear interval
  const origHandler = state.worker.onmessage;
  state.worker.onmessage = (e) => {
    if (e.data.type === "loaded" || e.data.type === "error") {
      clearInterval(state._progressInterval);
      dom.progressFill.style.width = "100%";
      dom.progressText.textContent = e.data.type === "loaded" ? "Ready!" : "Error";
      setTimeout(() => { dom.loadProgress.hidden = true; }, 600);
    }
    handleWorkerMessage(e);
  };
}

function setIndexStatus(status, text) {
  const dot = dom.indexBadge.querySelector(".badge-dot");
  dot.className = "badge-dot " + (status === "ready" ? "ready" : status === "loading" ? "loading" : "");
  dom.indexStatusText.textContent = text;
}

function doLiveSearch(sequence) {
  if (!state.indexLoaded || !state.worker) return;
  const cleaned = sequence.replace(/[^ACGTacgt]/g, "").toUpperCase();
  if (cleaned.length < 4) {
    dom.liveResults.innerHTML = `<div class="empty-state"><p>Enter at least 4 bases to search</p></div>`;
    return;
  }
  dom.sequenceInfo.textContent = `Searching: ${cleaned.length} bases`;
  const rid = `s_${Date.now()}`;
  state._searchStart = performance.now();
  state._searchRid = rid;
  state.worker.postMessage({ type: "search", pattern: cleaned, limit: 100, ctxChars: 40, rid });
}

function renderLiveResults(msg) {
  if (msg.rid !== state._searchRid) return;
  const elapsed = (performance.now() - state._searchStart).toFixed(2);
  const hits = msg.hits || [];

  if (hits.length === 0 && !msg.error) {
    dom.liveResults.innerHTML = `
      <div class="live-results-count">
        <span><span class="count">0</span> matches</span>
        <span class="search-time">${elapsed} ms</span>
      </div>
      <div class="empty-state"><p>No exact matches found in this pathogen's genome set</p></div>`;
    return;
  }

  if (msg.error) {
    dom.liveResults.innerHTML = `<div class="empty-state"><p style="color:var(--red)">${escapeHtml(msg.error)}</p></div>`;
    return;
  }

  let html = `
    <div class="live-results-count">
      <span><span class="count">${hits.length}</span> matches (showing up to 100)</span>
      <span class="search-time">${elapsed} ms</span>
    </div>`;

  for (const hit of hits) {
    const ctx = typeof hit === "string" ? hit : (hit.context || hit.ctx || JSON.stringify(hit));
    const file = hit.file || hit.f || "";
    const line = hit.line || hit.l || "";
    const col = hit.column || hit.col || hit.c || "";
    html += `
      <div class="result-card">
        <div class="result-header">
          <span class="result-file">${escapeHtml(file)}</span>
          <span class="result-pos">${line ? `L${line}` : ""}${col ? `:${col}` : ""}</span>
        </div>
        <div class="result-context">${highlightContext(ctx, dom.sequenceInput.value.replace(/[^ACGTacgt]/g, "").toUpperCase())}</div>
      </div>`;
  }

  dom.liveResults.innerHTML = html;
}

function highlightContext(ctx, pattern) {
  if (!pattern || !ctx) return escapeHtml(ctx);
  const escaped = escapeHtml(ctx);
  const patEsc = escapeHtml(pattern);
  // Case-insensitive highlight
  const re = new RegExp(`(${escapeRegex(patEsc)})`, "gi");
  return escaped.replace(re, '<span class="match">$1</span>');
}

// ---------------------------------------------------------------------------
// Export
// ---------------------------------------------------------------------------
function exportCSV() {
  if (!state.filteredTargets.length) return;
  const pathogen = state.data[state.selected];
  const scores = state.guideScores?.pathogens?.[state.selected] || {};
  const xr = state.crossReactivity?.pathogens?.[state.selected] || [];
  const xrMap = {};
  for (const g of xr) { xrMap[g.seq] = g; }
  const animals = state.crossReactivity?.animals || [];
  const header = "pathogen,sequence_23mer,gene_symbol,gene_name,gene_type,position,pam_type,occurrences,total_genomes,conservation_pct,guide_score,gc_pct,seed_gc_pct,max_homopolymer,host_specific" +
    (animals.length ? "," + animals.join(",") : "") + "\n";
  const rows = state.filteredTargets.map((t) => {
    const pct = pathogen.genomes > 0 ? (Math.min(t.o, pathogen.genomes) / pathogen.genomes * 100).toFixed(2) : '0';
    const sc = scores[t.s];
    const xrInfo = xrMap[t.s];
    const hostCols = animals.map(a => xrInfo?.host_hits?.[a] ?? '').join(',');
    return `${pathogen.name},${t.s},${t.g || ''},${(t.gn || '').replace(/,/g, ';')},${t.gt || ''},${t.p},${t.m},${t.o},${pathogen.genomes},${pct},${sc?.score ?? ''},${sc?.gc ?? ''},${sc?.seed_gc ?? ''},${sc?.max_homopolymer ?? ''},${xrInfo ? (xrInfo.specific ? 'yes' : 'no') : ''}` +
      (animals.length ? ',' + hostCols : '');
  }).join("\n");
  downloadFile(`${state.selected}_crispr_targets.csv`, header + rows, "text/csv");
}

function exportJSON() {
  if (!state.filteredTargets.length) return;
  const pathogen = state.data[state.selected];
  const scores = state.guideScores?.pathogens?.[state.selected] || {};
  const xr = state.crossReactivity?.pathogens?.[state.selected] || [];
  const xrMap = {};
  for (const g of xr) { xrMap[g.seq] = g; }
  const payload = {
    pathogen: pathogen.name,
    family: pathogen.family,
    total_genomes: pathogen.genomes,
    target_count: state.filteredTargets.length,
    exported: new Date().toISOString(),
    targets: state.filteredTargets.map((t) => {
      const sc = scores[t.s];
      const xrInfo = xrMap[t.s];
      return {
        sequence: t.s,
        position: t.p,
        pam_type: t.m,
        occurrences: t.o,
        gene_symbol: t.g || null,
        gene_name: t.gn || null,
        gene_type: t.gt || null,
        conservation_pct: pathogen.genomes > 0 ? +(Math.min(t.o, pathogen.genomes) / pathogen.genomes * 100).toFixed(2) : 0,
        guide_score: sc?.score ?? null,
        gc_pct: sc?.gc ?? null,
        seed_gc_pct: sc?.seed_gc ?? null,
        host_specific: xrInfo?.specific ?? null,
        host_hits: xrInfo?.host_hits ?? null,
      };
    }),
  };
  downloadFile(`${state.selected}_crispr_targets.json`, JSON.stringify(payload, null, 2), "application/json");
}

function downloadFile(name, content, mime) {
  const blob = new Blob([content], { type: mime });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = name;
  a.click();
  URL.revokeObjectURL(url);
}

// ---------------------------------------------------------------------------
// Event Binding
// ---------------------------------------------------------------------------
function bindEvents() {
  // Tabs — persist via URL hash
  dom.tabBar.addEventListener("click", (e) => {
    const btn = e.target.closest(".tab-btn");
    if (!btn) return;
    switchTab(btn.dataset.tab);
  });

  // Restore tab from URL hash on load
  const hashTab = location.hash.replace('#', '');
  if (hashTab && ['database', 'search', 'literature', 'context', 'sources', 'panel'].includes(hashTab)) {
    switchTab(hashTab);
  }

  // Handle back/forward
  window.addEventListener('hashchange', () => {
    const t = location.hash.replace('#', '');
    if (t && ['database', 'search', 'literature', 'context', 'sources', 'panel'].includes(t)) switchTab(t);
  });

  // Filter input (debounced)
  let filterTimer;
  dom.filterInput.addEventListener("input", () => {
    clearTimeout(filterTimer);
    filterTimer = setTimeout(() => {
      state.filterText = dom.filterInput.value.trim();
      state.page = 0;
      applyFilters();
    }, 200);
  });

  // PAM chips
  dom.pamFilter.addEventListener("click", (e) => {
    const chip = e.target.closest(".chip");
    if (!chip) return;
    dom.pamFilter.querySelectorAll(".chip").forEach((c) => c.classList.remove("active"));
    chip.classList.add("active");
    state.filterPam = chip.dataset.pam;
    state.page = 0;
    applyFilters();
  });

  // Export
  dom.exportCsv.addEventListener("click", exportCSV);
  dom.exportJson.addEventListener("click", exportJSON);

  // Gene filter
  if (dom.geneFilter) {
    dom.geneFilter.addEventListener("change", () => {
      state.filterGene = dom.geneFilter.value;
      state.page = 0;
      applyFilters();
    });
  }

  // Novelty filter
  if (dom.noveltyFilter) {
    dom.noveltyFilter.addEventListener("change", () => {
      state.filterNovelty = dom.noveltyFilter.value;
      state.page = 0;
      applyFilters();
    });
  }

  // Load index
  dom.loadIndexBtn.addEventListener("click", loadIndex);

  // Live search input (debounced)
  let searchTimer;
  dom.sequenceInput.addEventListener("input", () => {
    clearTimeout(searchTimer);
    searchTimer = setTimeout(() => {
      doLiveSearch(dom.sequenceInput.value);
    }, 250);
  });

  // Literature scanner
  if (dom.litSearchBtn) {
    dom.litSearchBtn.addEventListener("click", litSearchSingle);
    dom.litSequenceInput.addEventListener("keydown", (e) => {
      if (e.key === "Enter") litSearchSingle();
    });
    dom.litScanSelected.addEventListener("click", () => litBatchScan("All"));
    dom.litScanNGG.addEventListener("click", () => litBatchScan("NGG"));
  }

  // Citation copy
  if (dom.copyCitation) dom.copyCitation.addEventListener("click", copyCitationText);

  // Column sorting
  document.querySelectorAll('.data-table th.sortable').forEach(th => {
    th.addEventListener('click', () => {
      const col = th.dataset.col;
      if (state.sortCol === col) {
        state.sortAsc = !state.sortAsc;
      } else {
        state.sortCol = col;
        state.sortAsc = col === 's' || col === 'p' || col === 'g'; // alpha/pos/gene default asc, coverage desc
      }
      updateSortIndicators();
      if (state.filteredTargets.length) {
        sortTargets();
        state.page = 0;
        renderResults();
      }
    });
  });
  updateSortIndicators();
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
function updateSortIndicators() {
  document.querySelectorAll('.data-table th.sortable').forEach(th => {
    th.classList.remove('sort-asc', 'sort-desc');
    if (th.dataset.col === state.sortCol) {
      th.classList.add(state.sortAsc ? 'sort-asc' : 'sort-desc');
    }
  });
}

function colorDNA(seq) {
  return seq.split("").map((c) => {
    const cls = { A: "nuc-a", T: "nuc-t", C: "nuc-c", G: "nuc-g" }[c];
    return cls ? `<span class="${cls}">${c}</span>` : c;
  }).join("");
}

function pamBadge(pam) {
  if (pam === "NGG") return '<span class="pam-badge pam-ngg">NGG</span>';
  if (pam === "TTTN") return '<span class="pam-badge pam-tttn">TTTN</span>';
  return '<span class="pam-badge pam-none">—</span>';
}

function escapeHtml(str) {
  const div = document.createElement("div");
  div.textContent = str;
  return div.innerHTML;
}

function escapeRegex(str) {
  return str.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
}

function formatBytes(bytes) {
  if (bytes < 1024) return `${bytes} B`;
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`;
  if (bytes < 1024 * 1024 * 1024) return `${(bytes / (1024 * 1024)).toFixed(1)} MB`;
  return `${(bytes / (1024 * 1024 * 1024)).toFixed(2)} GB`;
}

// ---------------------------------------------------------------------------
// PubMed Literature Scanner
// ---------------------------------------------------------------------------
const EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils";
const PUBMED_DELAY_MS = 350; // stay under 3 req/s

async function pubmedSearch(query, maxResults = 10) {
  const url = `${EUTILS_BASE}/esearch.fcgi?db=pubmed&term=${encodeURIComponent(query)}&retmax=${maxResults}&retmode=json`;
  const resp = await fetch(url);
  if (!resp.ok) throw new Error(`PubMed search failed: HTTP ${resp.status}`);
  const data = await resp.json();
  const ids = data.esearchresult?.idlist || [];
  const count = parseInt(data.esearchresult?.count || "0", 10);
  if (ids.length === 0) return { count, papers: [] };

  await sleep(PUBMED_DELAY_MS);

  const sumUrl = `${EUTILS_BASE}/esummary.fcgi?db=pubmed&id=${ids.join(",")}&retmode=json`;
  const sumResp = await fetch(sumUrl);
  if (!sumResp.ok) throw new Error(`PubMed summary failed: HTTP ${sumResp.status}`);
  const sumData = await sumResp.json();

  const papers = ids.map((id) => {
    const r = sumData.result?.[id];
    if (!r) return null;
    return {
      pmid: id,
      title: r.title || "Untitled",
      authors: (r.authors || []).map((a) => a.name).slice(0, 5).join(", "),
      journal: r.source || "",
      year: r.pubdate?.split(" ")[0] || "",
      doi: (r.elocationid || "").replace("doi: ", ""),
    };
  }).filter(Boolean);

  return { count, papers };
}

function sleep(ms) { return new Promise((r) => setTimeout(r, ms)); }

async function litSearchSingle() {
  const query = dom.litSequenceInput.value.trim();
  if (!query) return;

  dom.litResults.innerHTML = '<div class="empty-state"><p>Searching PubMed...</p></div>';

  try {
    // Build a CRISPR-relevant query
    const pubQuery = buildCrisprQuery(query);
    const result = await pubmedSearch(pubQuery, 20);
    renderLitSingleResult(query, pubQuery, result);
  } catch (err) {
    dom.litResults.innerHTML = `<div class="empty-state"><p style="color:var(--red)">Error: ${escapeHtml(err.message)}</p></div>`;
  }
}

function buildCrisprQuery(input) {
  // If it looks like a DNA sequence (only ACGT), search as sequence + CRISPR context
  const cleaned = input.replace(/\s/g, "").toUpperCase();
  if (/^[ACGTU]+$/.test(cleaned) && cleaned.length >= 15) {
    return `("${cleaned}"[All Fields] OR "${cleaned}"[Title/Abstract]) AND (CRISPR[Title/Abstract] OR guide RNA[Title/Abstract] OR Cas9[Title/Abstract] OR Cas12[Title/Abstract])`;
  }
  // Otherwise treat as gene/pathogen keyword
  return `(${input}) AND (CRISPR[Title/Abstract] OR guide RNA[Title/Abstract] OR diagnostic[Title/Abstract])`;
}

function renderLitSingleResult(original, pubQuery, result) {
  if (result.papers.length === 0) {
    dom.litResults.innerHTML = `
      <div class="lit-summary">
        <h4>Results for: ${escapeHtml(original)}</h4>
        <div class="lit-stats-row">
          <div class="lit-stat novel"><span class="lit-stat-val">0</span><span class="lit-stat-lbl">Publications</span></div>
        </div>
        <p>No published CRISPR literature found for this query. <span class="novelty-badge novel">Potentially Novel</span></p>
        <div class="paper-query">PubMed query: ${escapeHtml(pubQuery)}</div>
      </div>`;
    return;
  }

  let html = `
    <div class="lit-summary">
      <h4>Results for: ${escapeHtml(original)}</h4>
      <div class="lit-stats-row">
        <div class="lit-stat published"><span class="lit-stat-val">${result.count.toLocaleString()}</span><span class="lit-stat-lbl">Total Publications</span></div>
        <div class="lit-stat"><span class="lit-stat-val">${result.papers.length}</span><span class="lit-stat-lbl">Shown</span></div>
      </div>
      <div class="paper-query">PubMed query: ${escapeHtml(pubQuery)}</div>
    </div>
    <div class="paper-list">`;

  for (const p of result.papers) {
    html += renderPaperCard(p);
  }

  html += `</div>`;
  dom.litResults.innerHTML = html;
}

function renderPaperCard(p) {
  const doiLink = p.doi ? ` · <a href="https://doi.org/${encodeURIComponent(p.doi)}" target="_blank" rel="noopener">DOI</a>` : "";
  return `
    <div class="paper-card">
      <div class="paper-title"><a href="https://pubmed.ncbi.nlm.nih.gov/${encodeURIComponent(p.pmid)}/" target="_blank" rel="noopener">${escapeHtml(p.title)}</a></div>
      <div class="paper-meta">
        <span class="authors">${escapeHtml(p.authors)}</span>
        <span class="journal">${escapeHtml(p.journal)}</span>
        <span class="year">${escapeHtml(p.year)}</span>
        <span>PMID: ${escapeHtml(p.pmid)}${doiLink}</span>
      </div>
    </div>`;
}

// ---------------------------------------------------------------------------
// Batch Scan — scan multiple targets against PubMed
// ---------------------------------------------------------------------------
let litScanAbort = false;

async function litBatchScan(pamFilter) {
  if (!state.selected || !state.data[state.selected]) return;
  litScanAbort = false;

  const pathogen = state.data[state.selected];
  let targets = pathogen.targets;

  if (pamFilter === "NGG") targets = targets.filter((t) => t.m === "NGG");

  // Cap to avoid excessive API calls (batch max 50)
  const maxTargets = 50;
  const scanTargets = targets.slice(0, maxTargets);
  const total = scanTargets.length;

  dom.litStatus.hidden = false;
  dom.litProgressFill.style.width = "0%";
  dom.litProgressText.textContent = `Scanning 0 / ${total} targets...`;
  dom.litResults.innerHTML = '<div class="empty-state"><p>Scanning in progress...</p></div>';

  // Disable scan buttons during scan
  dom.litScanSelected.disabled = true;
  dom.litScanNGG.disabled = true;

  const novelTargets = [];
  const publishedTargets = [];
  let scanned = 0;

  for (const t of scanTargets) {
    if (litScanAbort) break;

    try {
      const query = buildCrisprQuery(t.s);
      const result = await pubmedSearch(query, 3);

      if (result.count === 0) {
        novelTargets.push({ target: t, papers: [] });
      } else {
        publishedTargets.push({ target: t, papers: result.papers, totalPapers: result.count });
      }
    } catch (err) {
      // On error (likely rate limit), pause and continue
      await sleep(2000);
    }

    scanned++;
    const pct = Math.round((scanned / total) * 100);
    dom.litProgressFill.style.width = `${pct}%`;
    dom.litProgressText.textContent = `Scanned ${scanned} / ${total}... (${novelTargets.length} novel, ${publishedTargets.length} published)`;

    await sleep(PUBMED_DELAY_MS);
  }

  dom.litStatus.hidden = true;
  dom.litScanSelected.disabled = false;
  dom.litScanNGG.disabled = false;

  renderBatchResults(pathogen, scanTargets.length, novelTargets, publishedTargets, pamFilter);
}

function renderBatchResults(pathogen, totalScanned, novel, published, pamLabel) {
  let html = `
    <div class="lit-summary">
      <h4>${escapeHtml(pathogen.name)} — ${pamLabel || "All"} Targets Scan</h4>
      <div class="lit-stats-row">
        <div class="lit-stat"><span class="lit-stat-val">${totalScanned}</span><span class="lit-stat-lbl">Scanned</span></div>
        <div class="lit-stat novel"><span class="lit-stat-val">${novel.length}</span><span class="lit-stat-lbl">Novel (0 papers)</span></div>
        <div class="lit-stat published"><span class="lit-stat-val">${published.length}</span><span class="lit-stat-lbl">Published</span></div>
      </div>
      <p style="font-size:0.85rem;color:var(--text-2)">Scanned up to 50 targets. Targets with 0 PubMed results for CRISPR-related queries are marked as potentially novel.</p>
    </div>`;

  if (novel.length > 0) {
    html += `<div class="novel-targets-list"><h4>Potentially Novel Targets <span class="novelty-badge novel">${novel.length} targets</span></h4>`;
    for (const n of novel) {
      html += `<div class="novel-target-item"><span class="dna-seq">${colorDNA(n.target.s)}</span> ${pamBadge(n.target.m)} <span style="color:var(--text-3);font-size:0.78rem">pos ${n.target.p.toLocaleString()}</span></div>`;
    }
    html += `</div>`;
  }

  if (published.length > 0) {
    html += `<div style="margin-top:var(--sp-5)"><h4>Targets with Published Literature <span class="novelty-badge published">${published.length}</span></h4><div class="paper-list">`;
    for (const p of published) {
      html += `<div class="paper-card">
        <div class="paper-title">${colorDNA(p.target.s)} ${pamBadge(p.target.m)} — ${p.totalPapers} paper${p.totalPapers !== 1 ? "s" : ""}</div>`;
      for (const paper of p.papers) {
        html += `<div class="paper-meta" style="margin-top:4px">
          <a href="https://pubmed.ncbi.nlm.nih.gov/${encodeURIComponent(paper.pmid)}/" target="_blank" rel="noopener">${escapeHtml(paper.title)}</a>
          <span class="journal">${escapeHtml(paper.journal)}</span> <span class="year">${escapeHtml(paper.year)}</span>
        </div>`;
      }
      html += `</div>`;
    }
    html += `</div></div>`;
  }

  dom.litResults.innerHTML = html;
}

// ---------------------------------------------------------------------------
// Citation Copy
// ---------------------------------------------------------------------------
function copyCitationText() {
  if (!dom.citeBlock) return;
  const text = dom.citeBlock.textContent;
  navigator.clipboard.writeText(text).then(() => {
    dom.copyCitation.textContent = "Copied!";
    setTimeout(() => { dom.copyCitation.textContent = "Copy Citation"; }, 2000);
  });
}

// ---------------------------------------------------------------------------
// Ontology Enrichment
// ---------------------------------------------------------------------------
async function loadOntologyEnrichment() {
  try {
    const resp = await fetch("data/ontology-enrichment.json");
    if (!resp.ok) return;
    state.ontology = await resp.json();
  } catch { /* optional file, silently skip */ }
}

// Build a product-name → canonical-symbol map from ontology gene data.
// Falls back gracefully: if ontology is missing, direct matching is used.
function getGeneSymbolMap(pathogenKey) {
  if (state._geneSymbolMaps[pathogenKey]) return state._geneSymbolMaps[pathogenKey];
  const map = new Map(); // lowercase key → lowercase canonical symbol
  const ontGenes = state.ontology?.pathogens?.[pathogenKey]?.genes;
  if (ontGenes) {
    const geneList = Array.isArray(ontGenes)
      ? ontGenes
      : Object.values(ontGenes).flat();
    for (const g of geneList) {
      const sym = (g.symbol || "").trim();
      const name = (g.name || "").trim();
      if (sym) {
        map.set(sym.toLowerCase(), sym.toLowerCase());
        if (name) map.set(name.toLowerCase(), sym.toLowerCase());
      }
    }
  }
  state._geneSymbolMaps[pathogenKey] = map;
  return map;
}

// Resolve a target's gene value to a canonical symbol (lowercase).
// Tries: direct value, ontology name→symbol bridge.
function resolveGeneSymbol(geneValue, pathogenKey) {
  if (!geneValue) return "";
  const lc = geneValue.toLowerCase();
  const map = getGeneSymbolMap(pathogenKey);
  return map.get(lc) || lc;
}

// Build studied/unstudied Sets with lowercase symbols for a pathogen.
function getNoveltyLookup(pathogenKey) {
  const scanData = state.pubmedScan?.pathogens?.[pathogenKey];
  if (!scanData) return null;
  const studied = new Set((scanData.studied_genes || []).map(g => g.toLowerCase()));
  const unstudied = new Set((scanData.unstudied_genes || []).map(g => g.toLowerCase()));
  return { studied, unstudied };
}

async function loadDiagnosticPriority() {
  try {
    const resp = await fetch("data/diagnostic-priority.json");
    if (!resp.ok) return;
    state.diagnosticPriority = await resp.json();
    renderTopOpportunities();
    // Re-render pathogen grid to add priority badges
    renderPathogenGrid();
  } catch { /* optional file, silently skip */ }
}

function renderTopOpportunities() {
  if (!state.diagnosticPriority?.pathogens || !dom.topOpportunities) return;
  const entries = Object.values(state.diagnosticPriority.pathogens)
    .sort((a, b) => b.total_score - a.total_score)
    .slice(0, 3);
  if (entries.length === 0) return;

  dom.topOpportunities.style.display = '';
  dom.topOpportunities.innerHTML = `
    <div class="opps-header">
      <h3>Top Diagnostic Opportunities</h3>
      <span class="opps-sub">Pathogens ranked by unmet need for CRISPR-based diagnostics — score weights: Dx gap (30%), disease burden (25%), case fatality (15%), WHO class (15%), outbreak potential (15%)</span>
    </div>
    <div class="opps-cards">
      ${entries.map(e => {
        const factors = e.factors || {};
        const keyFact = factors.crispr_dx_gap?.detail || '';
        const burden = factors.disease_burden?.detail || '';
        const summary = [
          keyFact === 'none published' ? 'No CRISPR dx published' : (keyFact === 'limited' ? 'Limited CRISPR dx' : ''),
          burden ? `${burden} deaths/yr` : ''
        ].filter(Boolean).join(' · ');
        return `
        <div class="opp-card" data-key="${escapeHtml(e.pathogen)}" title="Click to view targets">
          <div class="opp-score">${e.total_score}<span class="opp-max">/100</span></div>
          <div class="opp-info">
            <div class="opp-name">${escapeHtml(e.species || e.pathogen)}</div>
            <div class="opp-detail">${escapeHtml(summary)}</div>
          </div>
        </div>`;
      }).join('')}
    </div>`;

  dom.topOpportunities.querySelectorAll('.opp-card').forEach(card => {
    card.addEventListener('click', () => viewTargetsFor(card.dataset.key));
  });
}

function renderDiseaseContext(pathogenKey) {
  if (!dom.contextContent) return;
  if (!state.ontology || !state.ontology.pathogens || !state.ontology.pathogens[pathogenKey]) {
    dom.contextContent.innerHTML = `<div class="empty-state"><p>No ontology data available for this pathogen. Run the ontology pipeline to generate data.</p></div>`;
    return;
  }

  const d = state.ontology.pathogens[pathogenKey];
  const pData = state.data?.[pathogenKey];
  let html = '';

  // Ontology ID chips — external links to authoritative databases
  html += `<div class="context-id-grid">`;
  if (d.ncbi_taxid) html += `<a class="context-id-chip" href="https://www.ncbi.nlm.nih.gov/taxonomy/${d.ncbi_taxid}" target="_blank" rel="noopener"><span class="id-label">NCBI TaxID</span><span class="id-val">${d.ncbi_taxid}</span></a>`;
  if (d.doid) html += `<a class="context-id-chip" href="https://disease-ontology.org/?id=${d.doid}" target="_blank" rel="noopener"><span class="id-label">Disease Ontology</span><span class="id-val">${d.doid}</span></a>`;
  if (d.mondo_id) html += `<a class="context-id-chip" href="https://mondo.monarchinitiative.org/disease/${d.mondo_id}" target="_blank" rel="noopener"><span class="id-label">MONDO</span><span class="id-val">${d.mondo_id}</span></a>`;
  if (d.icd11) html += `<a class="context-id-chip" href="https://icd.who.int/browse/2024-01/mms/en#${d.icd11}" target="_blank" rel="noopener"><span class="id-label">ICD-11</span><span class="id-val">${d.icd11}</span></a>`;
  html += `</div>`;

  // Organism info
  html += `<div class="context-meta-row">`;
  html += `<span class="context-meta-tag"><strong>Species:</strong> <em>${escapeHtml(d.species || '')}</em></span>`;
  if (d.genome_type) html += `<span class="context-meta-tag"><strong>Genome:</strong> ${escapeHtml(d.genome_type)}</span>`;
  if (d.transmission) html += `<span class="context-meta-tag"><strong>Transmission:</strong> ${escapeHtml(d.transmission)}</span>`;
  if (d.who_priority) html += `<span class="context-meta-tag context-who-priority">WHO Priority Pathogen</span>`;
  html += `</div>`;

  // Disease definition — strip trailing OBO URL citations
  const rawDef = d.disease_ontology?.definition || d.mondo?.definition || '';
  const cleanDef = rawDef.replace(/\s*\[url:http[^\]]*\]/g, '').replace(/\s*\[.*?\]$/, '').trim();
  if (cleanDef) {
    html += `<div class="context-definition"><h4>What is it?</h4><p>${escapeHtml(cleanDef)}</p></div>`;
  }

  // WHO Epi Context — the most actionable section
  const who = d.who_context;
  if (who) {
    html += `<div class="context-who-section"><h4>Why it matters — Epidemiological Context</h4><div class="context-epi-grid">`;
    if (who.classification) html += `<div class="epi-card"><span class="epi-label">WHO Classification</span><span class="epi-val">${escapeHtml(who.classification)}</span></div>`;
    if (who.case_fatality_rate) html += `<div class="epi-card"><span class="epi-label">Case Fatality Rate</span><span class="epi-val">${escapeHtml(who.case_fatality_rate)}</span></div>`;
    if (who.annual_cases) html += `<div class="epi-card"><span class="epi-label">Annual Cases (est.)</span><span class="epi-val">${escapeHtml(who.annual_cases)}</span></div>`;
    if (who.annual_deaths) html += `<div class="epi-card"><span class="epi-label">Annual Deaths (est.)</span><span class="epi-val">${escapeHtml(who.annual_deaths)}</span></div>`;
    html += `</div>`;
    if (who.geographic_spread) html += `<p class="epi-note"><strong>Geographic spread:</strong> ${escapeHtml(who.geographic_spread)}</p>`;
    if (who.diagnostic_need) html += `<p class="epi-note"><strong>Diagnostic need:</strong> ${escapeHtml(who.diagnostic_need)}</p>`;
    if (who.crispr_dx_status) html += `<p class="epi-note"><strong>CRISPR diagnostic status:</strong> ${escapeHtml(who.crispr_dx_status)}</p>`;
    if (who.existing_rapid_tests?.length) {
      html += `<p class="epi-note"><strong>Existing rapid tests:</strong> ${who.existing_rapid_tests.map(t => escapeHtml(t)).join(', ')}</p>`;
    } else {
      html += `<p class="epi-note epi-gap"><strong>Existing rapid tests:</strong> None — diagnostic gap</p>`;
    }
    html += `<p class="epi-note"><strong>Vaccine available:</strong> ${who.vaccine_available ? 'Yes' : 'No'}</p>`;
    html += `</div>`;
  }

  // Taxonomy lineage — structured chips only (skip raw lineage string)
  const tax = d.taxonomy;
  if (tax?.lineage_structured?.length) {
    html += `<details class="context-details"><summary><h4 style="display:inline">Taxonomy</h4></summary>`;
    html += `<div class="lineage-chips">`;
    for (const item of tax.lineage_structured) {
      html += `<span class="lineage-chip"><span class="lineage-rank">${escapeHtml(item.rank)}</span> ${escapeHtml(item.name)}</span>`;
    }
    html += `</div></details>`;
  }

  // Gene annotations — only show if genes have actual data
  const genes = d.genes;
  if (genes?.genes?.length > 0) {
    const validGenes = genes.genes.filter(g => g.symbol || g.name);
    if (validGenes.length > 0) {
      html += `<details class="context-details"><summary><h4 style="display:inline">Gene Map</h4> <span class="context-badge">${validGenes.length} genes</span></summary>`;
      html += `<p style="font-size:0.82rem;color:var(--text-3)">Reference genome: <a href="https://www.ncbi.nlm.nih.gov/datasets/genome/${escapeHtml(genes.reference_accession)}" target="_blank" rel="noopener">${escapeHtml(genes.reference_accession)}</a></p>`;
      html += `<div class="gene-table-container"><table class="data-table gene-table"><thead><tr><th>Symbol</th><th>Name</th><th>Start</th><th>End</th><th>Type</th></tr></thead><tbody>`;
      for (const g of validGenes.slice(0, 100)) {
        html += `<tr class="gene-map-row" data-gene-symbol="${escapeHtml(g.symbol || '')}">
          <td class="gene-symbol">${escapeHtml(g.symbol || '—')}</td>
          <td>${escapeHtml(g.name || '—')}</td>
          <td>${g.start != null ? Number(g.start).toLocaleString() : '—'}</td>
          <td>${g.end != null ? Number(g.end).toLocaleString() : '—'}</td>
          <td>${escapeHtml(g.type || '—')}</td>
        </tr>`;
      }
      html += `</tbody></table>`;
      if (validGenes.length > 100) html += `<p style="font-size:0.78rem;color:var(--text-3);margin-top:var(--sp-2)">Showing first 100 of ${validGenes.length} genes</p>`;
      html += `</div></details>`;
    }
  }

  // Synonyms — clean up OBO qualifiers at render time too (belt + suspenders)
  const rawSynonyms = d.disease_ontology?.synonyms || d.mondo?.synonyms || [];
  const synonyms = rawSynonyms
    .map(s => (typeof s === 'string' ? s : s.val || '').replace(/\s+(EXACT|RELATED|NARROW|BROAD)\s*(?:\S+\s*)?\[\]$/, '').trim())
    .filter(Boolean);
  if (synonyms.length > 0) {
    html += `<details class="context-details"><summary><h4 style="display:inline">Also known as</h4> <span class="context-badge">${synonyms.length}</span></summary>`;
    html += `<div class="synonym-chips">${synonyms.map(s => `<span class="synonym-chip">${escapeHtml(s)}</span>`).join('')}</div>`;
    html += `</details>`;
  }

  // Cross-references
  const xrefs = d.disease_ontology?.xrefs || [];
  if (xrefs.length > 0) {
    html += `<details class="context-details"><summary><h4 style="display:inline">Cross-References</h4> <span class="context-badge">${xrefs.length}</span></summary>`;
    html += `<div class="xref-chips">${xrefs.map(x => `<span class="xref-chip">${escapeHtml(typeof x === 'string' ? x : JSON.stringify(x))}</span>`).join('')}</div>`;
    html += `</details>`;
  }

  // Drug-resistance variant overlay
  const resData = state.resistanceData?.pathogens?.[pathogenKey];
  if (resData?.mutations?.length) {
    html += `<details class="context-details" open><summary><h4 style="display:inline">Drug-Resistance Regions</h4> <span class="context-badge">${resData.mutations.length} regions</span></summary>`;
    html += `<p style="font-size:0.82rem;color:var(--text-3);margin-bottom:var(--sp-2)">Source: ${escapeHtml(resData.source)}. Guides overlapping these regions can detect resistance mutations.</p>`;
    html += `<div class="gene-table-container"><table class="data-table"><thead><tr><th>Gene</th><th>Resistance</th><th>Drug</th><th>Region</th><th>Overlapping Guides</th><th>NGG Guides</th></tr></thead><tbody>`;
    for (const m of resData.mutations) {
      const sigCls = m.significance === 'high' ? 'style="color:#991B1B;font-weight:600"' : '';
      html += `<tr>
        <td><strong>${escapeHtml(m.gene)}</strong></td>
        <td ${sigCls}>${escapeHtml(m.name)}</td>
        <td>${escapeHtml(m.drug)}</td>
        <td style="font-family:var(--font-mono);font-size:0.78rem">${escapeHtml(m.region)}</td>
        <td style="text-align:right">${m.overlapping_guides}</td>
        <td style="text-align:right">${m.ngg_overlapping_guides}</td>
      </tr>`;
    }
    html += `</tbody></table></div></details>`;
  }

  // Ontology sources footer
  html += `<div class="context-sources"><p>Sources: NCBI Taxonomy · Disease Ontology (CC0) · MONDO (CC BY 4.0) · WHO · NCBI Gene · WHO Mutation Catalogs · Stanford HIVDB. All data bundled — works fully offline.</p></div>`;

  dom.contextContent.innerHTML = html;

  // Gene Map → Targets bidirectional link
  dom.contextContent.querySelectorAll('.gene-map-row[data-gene-symbol]').forEach(row => {
    const sym = row.dataset.geneSymbol;
    if (!sym || sym === '—') return;
    row.style.cursor = 'pointer';
    row.title = `Click to view targets in ${sym}`;
    row.addEventListener('click', () => {
      state.filterGene = sym;
      if (dom.geneFilter) {
        // Add gene to dropdown if not present (Gene Map may have genes without targets)
        if (![...dom.geneFilter.options].some(o => o.value === sym)) {
          const opt = document.createElement('option');
          opt.value = sym;
          opt.textContent = `${sym} (from Gene Map)`;
          dom.geneFilter.appendChild(opt);
        }
        dom.geneFilter.value = sym;
      }
      switchTab('database');
      applyFilters();
    });
  });
}

// ---------------------------------------------------------------------------
// Multiplexed Panel Designer
// ---------------------------------------------------------------------------
const PANEL_PRESETS = {
  respiratory: ['influenza-a', 'rsv', 'mers', 'tuberculosis'],
  hemorrhagic: ['ebola', 'dengue', 'zika'],
  all: null // filled dynamically
};

function initPanelDesigner() {
  const grid = document.getElementById('panelPathogenGrid');
  if (!grid || !state.data) return;
  const keys = Object.keys(state.data).sort();
  grid.innerHTML = keys.map(k => {
    const p = state.data[k];
    return `<label class="panel-pathogen-item"><input type="checkbox" value="${escapeHtml(k)}" /> ${escapeHtml(p.name)}</label>`;
  }).join('');
}

function selectPanelPreset(preset) {
  const grid = document.getElementById('panelPathogenGrid');
  if (!grid) return;
  const keys = preset === 'all'
    ? Object.keys(state.data)
    : (PANEL_PRESETS[preset] || []);
  grid.querySelectorAll('input[type="checkbox"]').forEach(cb => {
    cb.checked = keys.includes(cb.value);
  });
}
window.selectPanelPreset = selectPanelPreset;

function designPanel() {
  const grid = document.getElementById('panelPathogenGrid');
  const resultsDiv = document.getElementById('panelResults');
  if (!grid || !resultsDiv) return;

  const selected = [...grid.querySelectorAll('input[type="checkbox"]:checked')].map(cb => cb.value);
  if (selected.length < 2) {
    resultsDiv.style.display = 'block';
    resultsDiv.innerHTML = '<p class="panel-error">Select at least 2 pathogens to design a multiplexed panel.</p>';
    return;
  }

  const minConservation = parseInt(document.getElementById('panelMinConservation')?.value || '70');
  const minScore = parseInt(document.getElementById('panelMinScore')?.value || '70');
  const hostSpecificOnly = document.getElementById('panelHostSpecific')?.checked ?? true;

  const scores = state.guideScores?.pathogens || {};
  const xr = state.crossReactivity?.pathogens || {};

  // For each selected pathogen, gather candidate guides
  const candidates = {};
  for (const pk of selected) {
    const pathogen = state.data[pk];
    if (!pathogen) continue;
    const pScores = scores[pk] || {};
    const xrGuides = xr[pk] || [];
    const xrMap = {};
    for (const g of xrGuides) xrMap[g.seq] = g;
    const totalG = pathogen.genomes || 1;

    candidates[pk] = pathogen.targets
      .filter(t => {
        if (t.m !== 'NGG' && t.m !== 'ngg') return false;
        const pct = totalG > 0 ? (Math.min(t.o, totalG) / totalG * 100) : 0;
        if (pct < minConservation) return false;
        const sc = pScores[t.s];
        if (minScore > 0 && (!sc || sc.score < minScore)) return false;
        if (hostSpecificOnly) {
          const xrInfo = xrMap[t.s];
          if (xrInfo && !xrInfo.specific) return false;
        }
        return true;
      })
      .map(t => ({
        seq: t.s,
        gene: t.g || 'intergenic',
        pos: t.p,
        occ: t.o,
        pct: totalG > 0 ? (Math.min(t.o, totalG) / totalG * 100) : 0,
        score: pScores[t.s]?.score || 0,
        pathogen: pk
      }))
      .sort((a, b) => b.score - a.score || b.pct - a.pct);
  }

  // Greedy set-cover: pick best guide for each pathogen not yet covered
  // Additionally check that the guide is unique to its pathogen (not found as top candidate in others)
  const panel = [];
  const covered = new Set();

  // Build cross-lookup: which sequences appear as candidates in which pathogens
  const seqToPathogens = {};
  for (const [pk, cands] of Object.entries(candidates)) {
    for (const c of cands) {
      if (!seqToPathogens[c.seq]) seqToPathogens[c.seq] = new Set();
      seqToPathogens[c.seq].add(pk);
    }
  }

  // Round-robin greedy: pick best distinguishing guide for each uncovered pathogen
  let maxRounds = selected.length * 3;
  while (covered.size < selected.length && maxRounds-- > 0) {
    for (const pk of selected) {
      if (covered.has(pk)) continue;
      const cands = candidates[pk] || [];
      // Prefer guides unique to this pathogen
      const best = cands.find(c =>
        !panel.some(p => p.seq === c.seq) &&
        (seqToPathogens[c.seq]?.size === 1 || !cands.some(c2 => c2.seq !== c.seq && seqToPathogens[c2.seq]?.size === 1))
      ) || cands.find(c => !panel.some(p => p.seq === c.seq));

      if (best) {
        panel.push(best);
        covered.add(pk);
      }
    }
  }

  // Render results
  resultsDiv.style.display = 'block';
  if (panel.length === 0) {
    resultsDiv.innerHTML = '<p class="panel-error">No qualifying guides found. Try relaxing the filters.</p>';
    return;
  }

  const uncoveredList = selected.filter(k => !covered.has(k));
  let html = `<h3>Designed Panel: ${panel.length} guide${panel.length > 1 ? 's' : ''}</h3>`;
  if (uncoveredList.length > 0) {
    html += `<p class="panel-warning">⚠ Could not find qualifying guides for: ${uncoveredList.map(k => escapeHtml(state.data[k]?.name || k)).join(', ')}. Try relaxing filters.</p>`;
  }
  html += `<table class="panel-table"><thead><tr><th>Pathogen</th><th>Sequence (23-mer)</th><th>Gene</th><th>Position</th><th>Conservation</th><th>Score</th><th>Action</th></tr></thead><tbody>`;
  for (const g of panel) {
    const blastUrl = `https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&PROGRAM=blastn&DATABASE=nt&QUERY=${encodeURIComponent(g.seq)}`;
    html += `<tr>
      <td>${escapeHtml(state.data[g.pathogen]?.name || g.pathogen)}</td>
      <td class="dna-seq">${colorDNA(g.seq)}</td>
      <td>${escapeHtml(g.gene)}</td>
      <td>${g.pos.toLocaleString()}</td>
      <td>${g.pct.toFixed(1)}%</td>
      <td><span class="guide-score ${g.score >= 70 ? 'score-high' : g.score >= 40 ? 'score-med' : 'score-low'}">${g.score}</span></td>
      <td><a href="${blastUrl}" target="_blank" rel="noopener" class="blast-link">BLAST</a></td>
    </tr>`;
  }
  html += `</tbody></table>`;

  // Export button
  html += `<button class="btn-sm" onclick="exportPanelCSV()">Export Panel CSV</button>`;
  resultsDiv.innerHTML = html;

  // Stash for export
  state._lastPanel = panel;
}
window.designPanel = designPanel;

function exportPanelCSV() {
  if (!state._lastPanel) return;
  const rows = ['pathogen,sequence_23mer,gene,position,conservation_pct,guide_score'];
  for (const g of state._lastPanel) {
    rows.push(`${g.pathogen},${g.seq},${g.gene},${g.pos},${g.pct.toFixed(1)},${g.score}`);
  }
  const blob = new Blob([rows.join('\n')], { type: 'text/csv' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = `loom-panel-${new Date().toISOString().slice(0,10)}.csv`;
  a.click();
  URL.revokeObjectURL(a.href);
}
window.exportPanelCSV = exportPanelCSV;

// ---------------------------------------------------------------------------
// Boot
// ---------------------------------------------------------------------------
init();
