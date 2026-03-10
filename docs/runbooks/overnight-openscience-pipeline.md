# Overnight Open-Science Genome Pipeline

> **Script:** `scripts/overnight_openscience.sh`
> **Author:** Alvaro (with GitHub Copilot on Claude Opus 4.6)
> **Created:** 2026-03-04
> **Hardware:** Mac Studio M3 Ultra (28 cores, 96 GB), Samsung T9 2 TB SSD
> **Runtime:** 2–6 hours unattended

## What it does

Downloads **every public genome** for 12 priority pathogens, builds a BWT-based
FM-index for each, then runs a CRISPR target scan — all in one unattended
overnight run. Produces a combined CSV database of candidate diagnostic targets.

```
                ┌──────────────────────────────────────┐
                │   scripts/overnight_openscience.sh   │
                └──────────┬───────────────────────────┘
                           │ preflight
                           ▼
        ┌─────────────────────────────────────────┐
        │            tmux "openscience"            │
        │                                          │
        │  ┌─ 0:monitor ── health dashboard ────┐  │
        │  │  heartbeats / disk / progress       │  │
        │  └─────────────────────────────────────┘  │
        │                                          │
        │  ┌─ 1:download ──────────────────────┐   │
        │  │  aria2c (RefSeq, GRCh38)          │   │
        │  │  + NCBI datasets CLI (×12)        │   │
        │  │  max 4 parallel pathogen downloads │   │
        │  └──────────┬────────────────────────┘   │
        │             │ .fna lands on SSD          │
        │             ▼                            │
        │  ┌─ 2:index ────────────────────────┐   │
        │  │  loom index --build-mode auto     │   │
        │  │  3 parallel FM-index builds       │   │
        │  │  fires as downloads complete      │   │
        │  └──────────┬───────────────────────┘   │
        │             │ .idx ready                 │
        │             ▼                            │
        │  ┌─ 3:scan ─────────────────────────┐   │
        │  │  Python CorpusSearcher            │   │
        │  │  motif search + PAM extraction    │   │
        │  │  4 parallel scans                 │   │
        │  │  → per-pathogen CSVs              │   │
        │  │  → crispr_targets_combined.csv    │   │
        │  └──────────────────────────────────┘   │
        │                                          │
        │  ┌─ 4:logs ── tail -f all logs ──────┐  │
        │  └───────────────────────────────────┘  │
        └─────────────────────────────────────────┘
```

The three processing stages overlap: indexing starts as soon as the first FASTA
file lands, and scanning starts as soon as the first index is built. This
**pipeline parallelism** means a 2 GB Influenza-A download can be indexing while
Dengue is still downloading and SARS-CoV-2 is already being scanned.

## Quick start

```bash
# Dry run — see the plan, touch nothing
bash scripts/overnight_openscience.sh --dry-run

# Full run — launches tmux, go to bed
bash scripts/overnight_openscience.sh

# Attach later to monitor
tmux attach -t openscience     # window 0 = dashboard
```

## Prerequisites

| Tool | Purpose | Auto-installed? |
|------|---------|-----------------|
| `tmux` | Session persistence + multi-window | No (must exist) |
| `aria2c` | Multi-connection parallel downloads | Yes (via `brew`) |
| `datasets` | NCBI bulk genome download by taxon | Yes (from GitHub release) |
| `loom` | FM-index builder (`cargo build --release`) | Yes (auto-built if missing) |
| External SSD | Mounted at `/Volumes/T9` | No (must be mounted) |

The preflight step will install `aria2` and `datasets` via Homebrew if missing,
and build the loom binary from source if it doesn't exist.

## Directory layout on SSD

```
/Volumes/T9/loom-openscience/
├── raw/                          # Downloaded genomes (FASTA)
│   ├── sars-cov-2/
│   │   └── sars-cov-2.fna       # Extracted from NCBI zip
│   ├── influenza-a/
│   ├── dengue/
│   ├── rsv/
│   ├── hiv-1/
│   ├── ebola/
│   ├── zika/
│   ├── mers/
│   ├── mpox/
│   ├── hepatitis-b/
│   ├── tuberculosis/
│   ├── cholera/
│   ├── refseq-viral/             # All RefSeq viral (4 bulk files)
│   └── human-grch38/             # GRCh38 for off-target checking
├── indexes/                      # FM-indexes (.idx)
│   ├── sars-cov-2.idx
│   ├── influenza-a.idx
│   └── ...
├── crispr-db/                    # Scan results
│   ├── sars-cov-2_targets.csv
│   ├── influenza-a_targets.csv
│   ├── ...
│   └── crispr_targets_combined.csv  ← THE DELIVERABLE
├── logs/
│   ├── download.log
│   ├── index.log
│   ├── scan.log
│   └── aria2_*.log
├── .heartbeats/                  # Worker liveness files
│   ├── download
│   ├── index
│   └── scan
├── .done-download                # Stage completion markers
├── .done-index
├── .done-scan
├── monitor.sh                    # Auto-generated helper scripts
├── download.sh
├── index.sh
└── scan.sh
```

## Pathogens covered

| Pathogen | NCBI Taxon ID | Category |
|----------|---------------|----------|
| SARS-CoV-2 | 2697049 | virus |
| Influenza A | 11320 | virus |
| Dengue | 12637 | virus |
| RSV | 11250 | virus |
| HIV-1 | 11676 | virus |
| Ebola | 186538 | virus |
| Zika | 64320 | virus |
| MERS | 1335626 | virus |
| Mpox | 10244 | virus |
| Hepatitis B | 10407 | virus |
| Tuberculosis | 1773 | bacteria |
| Cholera | 666 | bacteria |

Plus: **RefSeq viral bulk** (all reference viral genomes, 4 files) and **GRCh38
human reference** (for off-target filtering).

## Health & safety mechanisms

### Heartbeat monitoring
Each worker (download, index, scan) touches a file in `.heartbeats/` periodically.
The monitor dashboard checks these every 30 seconds:
- **✅ active** — touched within 60s
- **⏳ idle** — 60–300s (normal during long operations)
- **⚠️ STALE** — >300s (something may be stuck)

### Disk space guard
Before every pathogen download, free space is checked. If the SSD drops below
**50 GB free**, downloads halt immediately. Index and scan stages continue on
already-downloaded data.

### Zero API cost
The script makes **zero LLM/Copilot/cloud API calls**. Every operation is local:
shell commands, `aria2c` for FTP downloads, `datasets` CLI for NCBI, the `loom`
Rust binary for indexing, and Python `CorpusSearcher` for scanning. You can
safely leave it running overnight with no risk of runaway API costs.

### Resume support
- `aria2c --continue=true` resumes interrupted downloads
- Index and scan stages skip already-completed work (checks for existing `.idx`
  and `_targets.csv` files)
- You can kill the session and re-run the script; it picks up where it left off

## Pipeline parallelism explained

Traditional approach (sequential):
```
Download ALL → Index ALL → Scan ALL
[====30 min====][====60 min====][====30 min====]  = 120 min total
```

This script (pipeline-overlapped):
```
Download: [==A==][==B==][==C==][==D==]...
Index:         [==A==][==B==][==C==]...
Scan:               [==A==][==B==]...
                                        = ~70 min total (42% faster)
```

The stages communicate via the filesystem:
1. **Download → Index**: Index stage polls `raw/*/` for new `.fna` files every 15s
2. **Index → Scan**: Scan stage polls `indexes/` for new `.idx` files every 15s
3. **Completion markers**: `.done-download`, `.done-index`, `.done-scan` signal stage completion

## CRISPR scanning strategy

The scan stage uses LOOM's `CorpusSearcher` (Python wrapper around the FM-index)
to find CRISPR-compatible target sequences:

1. **Motif search**: Searches for 9 diagnostic-relevant DNA motifs (regulatory
   signals, start codon contexts, restriction sites)
2. **Context extraction**: For each motif hit, extracts 23-mers from the surrounding
   context
3. **PAM detection**: Classifies each 23-mer by CRISPR PAM compatibility:
   - **NGG** (SpCas9): `[20bp target][N]GG`
   - **TTTN** (Cas12a): `TTT[N][23bp target]`
4. **Conservation markers**: Counts genome-wide ATG/TGA/TAA/TAG codons as
   conservation signals
5. **Deduplication**: Tracks seen targets to avoid duplicates across motifs

Output CSV columns: `pathogen, gene, position, sequence_23mer, pam_type, source_motif, occurrences`

## Adapting for your own use

### Change the target pathogens

Edit the `PATHOGENS` array at the top of the script. Format: `name|ncbi_taxon_id|estimated_mb|category`.
Find taxon IDs at [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy).

### Change the storage location

Set `SSD="/your/mount/point"` — any path with enough free space works. The script
creates its workspace directory automatically.

### Change parallelism limits

| Variable | Default | Where | What it controls |
|----------|---------|-------|------------------|
| `PARALLEL_LIMIT` | 4 | download.sh | Max concurrent pathogen downloads |
| `INDEX_PARALLEL` | 3 | index.sh | Max concurrent FM-index builds |
| `SCAN_PARALLEL` | 4 | scan.sh | Max concurrent CRISPR scans |
| `ARIA2_CONNECTIONS` | 6 | main script | aria2c connection count |
| `MIN_DISK_GB` | 50 | throughout | Disk space safety threshold |

### Use without external SSD

Change line 40 to point to a local directory:
```bash
SSD="$HOME"
WORKSPACE="$SSD/loom-openscience"
```

## Troubleshooting

| Symptom | Cause | Fix |
|---------|-------|-----|
| `FATAL: SSD not mounted` | T9 not plugged in / not mounted | Mount the drive or change `SSD` path |
| `DISK GUARD: Only X GB free` | SSD nearly full | Free space or raise `MIN_DISK_GB` |
| Worker heartbeat STALE | Long-running index build (normal for GRCh38) | Check `logs/index_*.log` — if progress is being made, it's fine |
| `datasets download` fails | NCBI rate limit or API change | Wait 5 min, re-run script (resumes) |
| Index build OOM | Genome too large for RAM | Edit `--build-mode streaming` in index.sh |

## Architecture notes for AI community

This script demonstrates several patterns useful for AI-assisted overnight automation:

1. **Tmux-first execution**: All long-running work runs inside tmux. The user can
   `tmux attach` at any time to inspect progress. This is strictly better than
   `nohup` or background processes.

2. **Pipeline parallelism via filesystem**: Instead of complex IPC, stages
   communicate through the filesystem — "a file appearing" is the signal. This is
   robust, inspectable, and trivially resumable.

3. **Heartbeat-based health monitoring**: Each worker touches a sentinel file.
   A dashboard polls these files. Simple, zero-dependency, works across languages.

4. **Zero-cost overnight**: The script is deliberately designed to make NO API
   calls. Every tool is local. A previous bug in an AI-assisted overnight script
   racked up $150 in idle Copilot requests — this design prevents that entirely.

5. **Self-installing prerequisites**: The preflight function installs missing
   tools automatically. The user runs one command and walks away.

6. **Idempotent re-runs**: Every operation checks for existing output before
   starting. Kill and restart at any point — nothing is lost or duplicated.

## Observed performance (2026-03-04)

From the initial run on Mac Studio M3 Ultra with 10 Gbps fiber:

| Metric | Value |
|--------|-------|
| GRCh38 download | 81 MB/s (10 seconds for 832 MB) |
| Influenza-A extract | 2.5 GB FASTA |
| Dengue extract | 248 MB FASTA |
| Pipeline start to first index build | ~30 seconds |
| Concurrent downloads observed | 6+ simultaneous |
