# Skill: Overnight Open-Science Pipeline

## When to use this skill

Use this skill when:
- The user asks to download pathogen genomes or run overnight genomics work
- The user asks to rebuild indexes or re-scan for CRISPR targets
- The user asks about the status of a running pipeline
- The user asks to add new pathogens or modify the scanning strategy
- You need to understand what data exists on the external SSD

## Script location

`scripts/overnight_openscience.sh` — the single entry point.

## Quick reference

```bash
# Dry run (safe — prints plan only)
bash scripts/overnight_openscience.sh --dry-run

# Full run (launches tmux session)
bash scripts/overnight_openscience.sh

# Attach to running session
tmux attach -t openscience

# Check if pipeline is running
tmux has-session -t openscience 2>/dev/null && echo "running" || echo "stopped"
```

## How the pipeline works

The script runs 3 stages with **pipeline parallelism** inside a tmux session called `openscience`:

| tmux window | Stage | What it does |
|-------------|-------|-------------|
| 0:monitor | Health | Auto-refreshing dashboard: disk, heartbeats, progress |
| 1:download | Download | Parallel genome downloads via aria2c + NCBI datasets CLI |
| 2:index | Index | FM-index builds via `loom index`, fires as downloads land |
| 3:scan | Scan | CRISPR target extraction via Python CorpusSearcher |
| 4:logs | Logs | `tail -f` on all log files |

Stages overlap — indexing begins as soon as the first FASTA file appears, scanning begins as soon as the first index is built. They communicate via filesystem:

- `raw/{pathogen}/{pathogen}.fna` — download output (triggers indexing)
- `indexes/{pathogen}.idx` — index output (triggers scanning)
- `crispr-db/{pathogen}_targets.csv` — scan output
- `.done-download`, `.done-index`, `.done-scan` — stage completion markers
- `.heartbeats/{download,index,scan}` — worker liveness files

## Data locations

All pipeline data lives on the external SSD:

```
/Volumes/T9/loom-openscience/
├── raw/              # FASTA files per pathogen
├── indexes/          # FM-index files (.idx)
├── crispr-db/        # Scan results (CSV)
│   └── crispr_targets_combined.csv  ← main deliverable
├── logs/             # All logs
└── .heartbeats/      # Worker health files
```

## Checking pipeline status

To determine the current state of the pipeline:

```bash
# Which stages have completed?
ls /Volumes/T9/loom-openscience/.done-* 2>/dev/null

# How many pathogens downloaded?
ls /Volumes/T9/loom-openscience/raw/

# How many indexed?
ls /Volumes/T9/loom-openscience/indexes/*.idx 2>/dev/null

# How many scanned?
ls /Volumes/T9/loom-openscience/crispr-db/*_targets.csv 2>/dev/null

# Combined results count
wc -l /Volumes/T9/loom-openscience/crispr-db/crispr_targets_combined.csv 2>/dev/null

# Recent log activity
tail -5 /Volumes/T9/loom-openscience/logs/{download,index,scan}.log 2>/dev/null
```

## Adding a new pathogen

1. Find the NCBI Taxonomy ID at https://www.ncbi.nlm.nih.gov/taxonomy
2. Add an entry to the `PATHOGENS` array in `scripts/overnight_openscience.sh`:
   ```bash
   "pathogen-name|TAXON_ID|estimated_mb|virus"
   ```
3. Also add to the `PATHOGENS` array inside the heredoc for `download.sh`:
   ```bash
   "pathogen-name|TAXON_ID"
   ```
4. Re-run the script — it will skip already-completed work and only process the new pathogen.

## Important constraints

1. **External SSD required**: The script expects `/Volumes/T9` mounted. Change `SSD=` variable to use a different path.
2. **macOS-specific**: Uses `stat -f %m` (macOS syntax), `df -g`, and Homebrew. For Linux, these need adaptation.
3. **Zero API cost**: Makes NO LLM/Copilot/cloud API calls. All tools are local (aria2c, datasets CLI, loom binary, Python).
4. **Idempotent**: Safe to re-run. Skips completed downloads, indexes, and scans.
5. **Disk guard**: Stops downloads if SSD free space drops below 50 GB.

## Dependencies

| Tool | Version tested | Install |
|------|---------------|---------|
| tmux | 3.6a | `brew install tmux` |
| aria2c | 1.37.0 | `brew install aria2` (auto-installed by script) |
| datasets | 18.19.0 | Downloaded from NCBI GitHub releases (auto-installed) |
| loom | 2.2 MB release | `cargo build --release --bin loom` (auto-built if missing) |
| Python 3 | 3.14+ | Required for scan stage (CorpusSearcher) |

## The CRISPR scan

The scan stage (`scan.sh`) uses LOOM's Python `CorpusSearcher` to search each FM-index:

1. Searches for 9 diagnostic DNA motifs (TATAAA, ATGAAA, GGATCC, etc.)
2. Extracts 23-mers from the context surrounding each hit
3. Classifies 23-mers by CRISPR PAM compatibility (NGG for SpCas9, TTTN for Cas12a)
4. Counts genome-wide codon markers (ATG, TGA, TAA, TAG) for conservation signal
5. Deduplicates and writes per-pathogen CSV + combined CSV

Python API used:
```python
from loom.search_executor import CorpusSearcher, SearchResponse, SearchResult

searcher = CorpusSearcher(idx_path)          # Load an .idx file
resp = searcher.search(term, context_lines=2, max_results=100)
# resp.count — total occurrences
# resp.results — list[SearchResult] with .file, .line, .column, .context
count = searcher.count(term)                  # Fast count-only
```

## Troubleshooting for agents

- **"SSD not mounted"**: The user needs to plug in the T9 drive. Don't try to work around this.
- **Pipeline running but no new progress**: Check heartbeats. If all stale >10 min, a worker may have crashed. Check `logs/*.log` for errors.
- **Want to restart a specific stage**: Kill the tmux window (`tmux send-keys -t openscience:download C-c`), fix the issue, then re-run the script.
- **Out of disk**: The pipeline won't corrupt data if it runs out of space — the disk guard stops downloads, and index/scan stages only process what's already downloaded.
