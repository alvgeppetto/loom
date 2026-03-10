# Project LOOM — Agent Guidelines

## Overview

LOOM is a BWT (Burrows-Wheeler Transform) powered exact-retrieval engine. A 195 KB WASM binary
performs sub-millisecond search across millions of documents. Primary application: open-science
CRISPR target discovery and an in-browser search tool for pathogen genomics.

## Hardware Profile

- **Machine:** Mac Studio M3 Ultra — 28 cores (20P + 8E), 96 GB unified RAM, ARM64/NEON
- **Primary target:** ARM64 + NEON SIMD. Never assume x86/AVX/CUDA.
- **Edge targets:** Raspberry Pi 5 (8 GB), iPads, phones, browsers (WASM)
- **BLAS:** Apple Accelerate (native). No external BLAS deps.
- **All agents MUST read** `.github/skills/performance-mandates/SKILL.md` before modifying Rust code.

## Language & Stack

- **Core (indexer, search):** Rust (crate: `brenda`)
- **WASM build:** `wasm-pack` → `pkg/` (195 KB)
- **LLM orchestration (CLI only):** Python (Ollama/llama.cpp bindings)

## Workflow

### Git Branches

1. Work on feature branches: `feature/<short-name>` or `fix/<issue>`
2. Commit frequently with conventional commits: `feat:`, `fix:`, `docs:`, `test:`, `refactor:`
3. Merge to `main` only after tests pass
4. Delete branch after merge

### Task Flow

1. Check `TODO-v3.md` for current priorities
2. Pick a task, create branch
3. Implement with tests (Cucumber specs where applicable)
4. Update `docs/decisions/` if architectural choices made
5. Update `docs/lessons/` if something was learned
6. Mark task complete in `TODO-v3.md`
7. Merge to main

### After Each Task

Always:
1. Show current `TODO-v3.md` state
2. Propose specific next task from TODO-v3

## Code Style

### Rust

- Use `rustfmt` defaults
- Prefer `Result<T, E>` over panics
- Document public APIs with `///` doc comments
- Use `thiserror` for error types

### Python

- Black formatter, 88 char lines
- Type hints on all public functions
- Docstrings in Google style

## Testing

- **Rust:** `cargo test`
- Write tests BEFORE implementation when possible (TDD)

## Documentation

- API docs → inline doc comments

## Key Directories

```
brenda/           # FM-index Rust crate (the core engine)
src/              # Rust source
  cli/            # Command-line interface
  dna/            # DNA corpus/benchmark utilities
  wasm.rs         # WASM bindings
web/              # CRISPR search website
scripts/          # Pipeline scripts (openscience)
publications/     # CRISPR paper
docs/
  runbooks/       # Pipeline & setup runbooks
data/
  crispr_guides/  # CRISPR reference data
```

## Reference Files

- [TODO-v3.md](TODO-v3.md) — Current task list

## Agent Skills

- [Performance Mandates Skill](.github/skills/performance-mandates/SKILL.md) — **MANDATORY** for any Rust/perf work. Green compute, commodity hardware, math-first optimization, facade doctrine.
- [Overnight Pipeline Skill](.github/skills/overnight-pipeline/SKILL.md) — Genome download, FM-index build, and CRISPR scan pipeline. Read before modifying `scripts/overnight_openscience.sh` or checking pipeline status.
- [Ontology-First Skill](.github/skills/ontology-first/SKILL.md) — **MANDATORY** for any cross-referencing of gene/protein/pathogen names across data sources. Enforces ontology-backed synonym expansion to prevent false negatives.
- [Citation Check Skill](.github/skills/citation-check/SKILL.md) — **MANDATORY** before building any paper PDF. Run `scripts/verify-citations.py --offline <paper>` after any edit to the paper markdown. Checks for placeholders, missing refs, missing DOIs, and (optionally live) DOI resolution.

## Active Architecture Boundary

- **Production path:** `brenda/` (FM-index crate) + `src/cli/` + `src/dna/` + `web/` + `pkg/`
- **Publications:** `publications/` (CRISPR paper only)
- **Archive:** `archive/` is gitignored — pre-pivot research code lives locally only.
