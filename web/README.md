# LOOM Edge Web UI

This is a browser UI for layperson-style Q&A over your indexed codebase.

- Natural language parsing + synthesis: local Ollama
- Exact retrieval: Rust FM-index compiled to WebAssembly
- Default model: `qwen2.5-coder:7b`

---

## ⛔ PRE-LAUNCH CHECKLIST — MUST DO BEFORE GOING PUBLIC

> **LEGAL HOLD:** Do NOT push to any public GitHub repo until lawyer approval is obtained.
> See also: the user memory note `legal-holds.md`.

### 1. Replace Azure SWA URLs with public GitHub repo URL
- `crispr-search.html` — **Methods block** and **How To Cite block** currently show
  `https://calm-mushroom-0185d800f.4.azurestaticapps.net` as the tool URL.
  Replace with the real public GitHub URL once the repo is created.
- Search for `calm-mushroom` across all web files to catch any other occurrences.

### 2. Audit all external links
- All links verified working as of 2026-03-05. Re-check before launch:
  - BLAST (NCBI) — `blast.ncbi.nlm.nih.gov`
  - NCBI Datasets, Taxonomy, Gene, Books
  - Disease Ontology — `disease-ontology.org`
  - MONDO — `mondo.monarchinitiative.org`
  - WHO Disease Outbreak News
  - ICD-11 — `icd.who.int`
  - PubMed / E-utilities
  - DOI links (dynamic, via `doi.org`)

### 3. Remove or update placeholder text
- Grep for `[repo]`, `[TODO]`, `[TBD]`, `placeholder` across `web/` and `docs/`.
- Check `crispr-search.html` citation block — make sure author/org name is correct.

### 4. Azure Blob Storage access
- FM-index files are served from `loomcrisprdata.blob.core.windows.net/indexes`.
  Confirm the container has public read access, or switch to a CDN/signed-URL approach.

### 5. Privacy & analytics
- No tracking/analytics are currently embedded. Decide if any are needed.
- The site makes no outbound calls except: Google Fonts, NCBI E-utilities (user-initiated),
  Azure Blob (FM-index downloads). Verify this is acceptable.

### 6. License
- Add a LICENSE file to the repo (currently absent from `web/`).
- Confirm data licensing: Disease Ontology (CC0), MONDO (CC BY 4.0), NCBI (public domain).

### 7. Security review
- WASM binary is loaded from `pkg/`. Ensure subresource integrity (SRI) hashes if hosting
  on a CDN.
- All user input is passed through `escapeHtml()` — verify no XSS vectors remain.
- PubMed queries use `encodeURIComponent()` — verify no injection paths.

---

## Prerequisites

- `wasm-pack` installed
- Ollama running locally with at least `qwen2.5-coder:7b`
- Existing index file at repository root: `rabbit.idx`

## Build WASM package

Run from repository root:

```bash
wasm-pack build --target web --out-dir pkg
```

## Run local server

From repository root:

```bash
python -m http.server 8080
```

Open:

- `http://localhost:8080/web/`

## Usage

1. Wait for status to show: `Ready · ... files`
2. Ask a plain English question
3. Click `Ask LOOM`
4. Inspect generated terms + evidence cards + final answer

---

## Runtime modes

The UI exposes three runtime modes selectable via the **Runtime Mode** dropdown:

| Mode | Description | Requires Ollama |
|---|---|---|
| ⚡ Keyword-only | Pure WASM FM-index search. Fastest, instant results, no LLM. | No |
| 🤏 Local model | Ollama expands terms; WASM retrieves evidence. No synthesis step. | Yes |
| 🧠 Full synthesis | Ollama expands terms + ontology enrichment + WASM retrieval + LLM answer. | Yes |

### End-to-end flow (Full synthesis)

```
User question
  → Route planner (query intent classification)
  → LLM term expansion (Ollama)
  → Ontology enrichment (embedded ONTO section from rabbit.idx via engine.enrichFromOntology())
  → Merged term list (LLM terms + ontology terms, deduped, top-18)
  → Function candidate discovery (WASM FM-index)
  → Evidence retrieval (WASM FM-index search per term)
  → Synthesis (Ollama LLM with retrieved code snippets as context)
  → Answer + telemetry displayed
```

The ontology section is embedded directly in the `.idx` binary file. No external database or
network call is needed for ontology lookup — WASM reads it from the same byte buffer loaded at boot.

### Degradation chain

When Ollama is unreachable (network error, `ECONNREFUSED`, HTTP 5xx), the UI automatically
degrades to **keyword-only** mode rather than displaying an error:

```
Full synthesis  ──(Ollama down)──▶  keyword-only (WASM FM-index only)
Local model     ──(Ollama down)──▶  keyword-only (WASM FM-index only)
Keyword-only    ── always works (no Ollama dependency)
```

A yellow **"Ollama unavailable — keyword mode"** status banner is shown when degradation occurs.
Retrieved hits are still displayed so the query is never a complete failure.

---

## Smoke-test commands

### 1. JavaScript syntax check

```bash
node --check web/app.js
```

Expected: no output, exit code 0.

### 2. WASM build

```bash
wasm-pack build --target web --out-dir pkg
```

Expected: `pkg/loom_bg.wasm` and `pkg/loom.js` created.

### 3. Full browser smoke test (manual)

```bash
python -m http.server 8080
# then open http://localhost:8080/web/ in a browser
# verify status shows: "Ready · N files · X.XX MB corpus"
# run a keyword-only query (no Ollama needed): select "Keyword-only" mode and ask "queue"
```

### 4. Node.js WASM performance benchmark (headless)

Requires a `pkg-node/` build (`wasm-pack build --target nodejs --out-dir pkg-node`):

```bash
node scripts/browser_wasm_bench.js rabbit.idx
```

---

## Notes

- This is a local-first edge prototype. No cloud backend is required.
- The keyword-only mode is fully offline — only `rabbit.idx` and the WASM binary are needed.
- If Ollama is not reachable, the UI degrades gracefully to keyword-only (see above).
- You can switch to smaller models in the Model dropdown.
