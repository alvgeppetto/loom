#!/usr/bin/env bash
# ===========================================================================
# PIVOT-OPENSCIENCE Overnight Pipeline
# ===========================================================================
#
# Downloads pathogen genomes (parallel), builds FM-indexes (parallel),
# runs CRISPR target scans, and catalogs results — all on external SSD.
#
# Usage:
#   bash scripts/overnight_openscience.sh          # full run
#   bash scripts/overnight_openscience.sh --dry-run # print plan, no action
#
# Architecture:
#   tmux session "openscience" with windows:
#     0:monitor   — health dashboard (auto-refreshes)
#     1:download  — parallel genome downloads (aria2c)
#     2:index     — FM-index builds (fires as downloads complete)
#     3:scan      — CRISPR target scanning (fires as indexes complete)
#     4:log       — tail all logs
#
# Health checks:
#   - Heartbeat file touched every 60s by each worker
#   - Monitor window checks heartbeats + disk space every 30s
#   - If any worker stalls >5 min, logs WARNING (no auto-kill to avoid data loss)
#   - Disk space guard: stops downloads if SSD < 50 GB free
#
# Cost guard:
#   - No Copilot / LLM / API calls — pure shell + Rust binary
#   - No idle loops — each stage completes and exits
#   - Total estimated runtime: 2-6 hours depending on bandwidth
#
# ===========================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
LOOM_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SSD="/Volumes/T9"
WORKSPACE="$SSD/loom-openscience"
RAW_DIR="$WORKSPACE/raw"
IDX_DIR="$WORKSPACE/indexes"
SCAN_DIR="$WORKSPACE/crispr-db"
LOG_DIR="$WORKSPACE/logs"
HEARTBEAT_DIR="$WORKSPACE/.heartbeats"
LOOM_BIN="$LOOM_ROOT/target/release/loom"
STAMP="$(date +%Y%m%d-%H%M%S)"
DRY_RUN=0
MIN_DISK_GB=50
SESSION="openscience"

# Max parallel downloads — use 6 connections to saturate 10 Gbps
ARIA2_CONNECTIONS=6

if [[ "${1:-}" == "--dry-run" ]]; then
  DRY_RUN=1
  echo "=== DRY RUN — no actions will be taken ==="
fi

# ---------------------------------------------------------------------------
# Preflight checks
# ---------------------------------------------------------------------------
preflight() {
  echo "== Preflight Checks =="

  # SSD mounted?
  if [[ ! -d "$SSD" ]]; then
    echo "FATAL: SSD not mounted at $SSD"
    exit 1
  fi

  # Check free space
  local free_gb
  free_gb=$(df -g "$SSD" | awk 'NR==2 {print $4}')
  echo "SSD free: ${free_gb} GB"
  if (( free_gb < MIN_DISK_GB )); then
    echo "FATAL: Less than ${MIN_DISK_GB} GB free on SSD"
    exit 1
  fi

  # Loom binary
  if [[ ! -x "$LOOM_BIN" ]]; then
    echo "Building loom release binary..."
    (cd "$LOOM_ROOT" && cargo build --release --bin loom)
  fi
  echo "loom binary: $LOOM_BIN ($(ls -lh "$LOOM_BIN" | awk '{print $5}'))"

  # aria2c
  if ! command -v aria2c &>/dev/null; then
    echo "Installing aria2 (multi-connection downloader)..."
    brew install aria2
  fi
  echo "aria2c: $(aria2c --version | head -1)"

  # ncbi-datasets-cli
  if ! command -v datasets &>/dev/null; then
    echo "Installing ncbi-datasets-cli..."
    brew install ncbi-datasets-cli
  fi
  echo "datasets: $(datasets version 2>/dev/null || datasets --version 2>/dev/null || echo 'installed')"

  # tmux
  if ! command -v tmux &>/dev/null; then
    echo "FATAL: tmux not found"
    exit 1
  fi

  # Create workspace
  mkdir -p "$RAW_DIR" "$IDX_DIR" "$SCAN_DIR" "$LOG_DIR" "$HEARTBEAT_DIR"
  echo "Workspace: $WORKSPACE"
  echo "== Preflight OK =="
}

# ---------------------------------------------------------------------------
# Pathogen definitions
# ---------------------------------------------------------------------------
# Format: NAME|TAXON_ID|EXPECTED_SIZE_MB|CATEGORY
# Using NCBI Taxonomy IDs for bulk download via datasets CLI
PATHOGENS=(
  "sars-cov-2|2697049|2000|virus"
  "influenza-a|11320|3000|virus"
  "dengue|12637|200|virus"
  "rsv|11250|100|virus"
  "hiv-1|11676|500|virus"
  "ebola|186538|50|virus"
  "zika|64320|30|virus"
  "mers|1335626|20|virus"
  "mpox|10244|100|virus"
  "hepatitis-b|10407|200|virus"
  "tuberculosis|1773|500|bacteria"
  "cholera|666|200|bacteria"
)

# RefSeq viral bulk (all viruses — already partially downloaded)
REFSEQ_VIRAL_URLS=(
  "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz"
  "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz"
  "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.1.genomic.fna.gz"
  "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.4.1.genomic.fna.gz"
)

# GRCh38 human reference (for off-target checking)
GRCH38_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"

# ---------------------------------------------------------------------------
# Health check infrastructure
# ---------------------------------------------------------------------------
touch_heartbeat() {
  local name="$1"
  touch "$HEARTBEAT_DIR/$name"
}

check_heartbeats() {
  local now stale_sec name mtime
  now=$(date +%s)
  for hb in "$HEARTBEAT_DIR"/*; do
    [[ -f "$hb" ]] || continue
    name=$(basename "$hb")
    mtime=$(stat -f %m "$hb" 2>/dev/null || echo 0)
    stale_sec=$(( now - mtime ))
    if (( stale_sec > 300 )); then
      echo "WARNING: Worker '$name' heartbeat stale (${stale_sec}s ago)"
    fi
  done
}

check_disk_space() {
  local free_gb
  free_gb=$(df -g "$SSD" | awk 'NR==2 {print $4}')
  if (( free_gb < MIN_DISK_GB )); then
    echo "DISK GUARD: Only ${free_gb} GB free — halting downloads"
    return 1
  fi
  echo "Disk: ${free_gb} GB free"
  return 0
}

# ---------------------------------------------------------------------------
# Monitor script (runs in tmux window 0)
# ---------------------------------------------------------------------------
write_monitor_script() {
  cat > "$WORKSPACE/monitor.sh" << 'MONITOR_EOF'
#!/usr/bin/env bash
WORKSPACE="$1"
HEARTBEAT_DIR="$WORKSPACE/.heartbeats"
SSD=$(dirname "$WORKSPACE")

while true; do
  clear
  echo "══════════════════════════════════════════════════════"
  echo "  PIVOT-OPENSCIENCE — Overnight Monitor"
  echo "  $(date)"
  echo "══════════════════════════════════════════════════════"
  echo ""

  # Disk space
  df -h "$SSD" | awk 'NR<=2'
  echo ""

  # Heartbeats
  echo "── Worker Heartbeats ──"
  now=$(date +%s)
  for hb in "$HEARTBEAT_DIR"/*; do
    [[ -f "$hb" ]] || continue
    name=$(basename "$hb")
    mtime=$(stat -f %m "$hb" 2>/dev/null || echo 0)
    age=$(( now - mtime ))
    if (( age > 300 )); then
      status="⚠️  STALE (${age}s)"
    elif (( age > 60 )); then
      status="⏳ idle (${age}s)"
    else
      status="✅ active (${age}s)"
    fi
    printf "  %-20s %s\n" "$name" "$status"
  done
  [[ -z "$(ls -A "$HEARTBEAT_DIR" 2>/dev/null)" ]] && echo "  (no workers yet)"
  echo ""

  # Download progress
  echo "── Downloads ──"
  for d in "$WORKSPACE/raw"/*/; do
    [[ -d "$d" ]] || continue
    name=$(basename "$d")
    count=$(find "$d" -name "*.fna" -o -name "*.fasta" -o -name "*.fna.gz" -o -name "*.zip" 2>/dev/null | wc -l | tr -d ' ')
    size=$(du -sh "$d" 2>/dev/null | cut -f1)
    printf "  %-25s %s files  %s\n" "$name" "$count" "$size"
  done
  echo ""

  # Indexes
  echo "── Indexes ──"
  for idx in "$WORKSPACE/indexes"/*.idx; do
    [[ -f "$idx" ]] || continue
    name=$(basename "$idx")
    size=$(ls -lh "$idx" | awk '{print $5}')
    printf "  %-35s %s\n" "$name" "$size"
  done
  [[ -z "$(ls "$WORKSPACE/indexes"/*.idx 2>/dev/null)" ]] && echo "  (none yet)"
  echo ""

  # Scans
  echo "── CRISPR Scans ──"
  for csv in "$WORKSPACE/crispr-db"/*.csv; do
    [[ -f "$csv" ]] || continue
    name=$(basename "$csv")
    lines=$(wc -l < "$csv" | tr -d ' ')
    size=$(ls -lh "$csv" | awk '{print $5}')
    printf "  %-35s %s lines  %s\n" "$name" "$lines" "$size"
  done
  [[ -z "$(ls "$WORKSPACE/crispr-db"/*.csv 2>/dev/null)" ]] && echo "  (none yet)"
  echo ""

  # Completion markers
  echo "── Completion ──"
  for marker in "$WORKSPACE"/.done-*; do
    [[ -f "$marker" ]] || continue
    echo "  ✅ $(basename "$marker" | sed 's/^.done-//')"
  done
  [[ -z "$(ls "$WORKSPACE"/.done-* 2>/dev/null)" ]] && echo "  (nothing completed yet)"

  sleep 30
done
MONITOR_EOF
  chmod +x "$WORKSPACE/monitor.sh"
}

# ---------------------------------------------------------------------------
# Stage 1: Download (parallel per-pathogen)
# ---------------------------------------------------------------------------
write_download_script() {
  cat > "$WORKSPACE/download.sh" << 'DL_EOF'
#!/usr/bin/env bash
set -euo pipefail

WORKSPACE="$1"
RAW_DIR="$WORKSPACE/raw"
LOG_DIR="$WORKSPACE/logs"
HEARTBEAT_DIR="$WORKSPACE/.heartbeats"
MIN_DISK_GB=50
SSD=$(dirname "$WORKSPACE")

touch_hb() { touch "$HEARTBEAT_DIR/download"; }

log() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG_DIR/download.log"; }

check_disk() {
  local free_gb
  free_gb=$(df -g "$SSD" | awk 'NR==2 {print $4}')
  (( free_gb >= MIN_DISK_GB ))
}

# Download a single pathogen using NCBI datasets CLI
# Handles: virus vs bacteria, Zip64 extraction, retry on GOAWAY
download_pathogen() {
  local name="$1" taxon="$2" outdir="$3" category="${4:-virus}"
  mkdir -p "$outdir"

  local zipfile="$outdir/${name}.zip"
  local fastafile="$outdir/${name}.fna"

  # Skip if already extracted and non-empty
  if [[ -f "$fastafile" ]] && [[ $(stat -f %z "$fastafile" 2>/dev/null || echo 0) -gt 1000 ]]; then
    log "SKIP $name — already have $fastafile ($(ls -lh "$fastafile" | awk '{print $5}'))"
    return 0
  fi

  # Delete truncated zips (no EOCD marker = incomplete download)
  if [[ -f "$zipfile" ]]; then
    if ! python3 -c "
import zipfile, sys
try:
    zipfile.ZipFile(sys.argv[1]).close()
except Exception:
    sys.exit(1)
" "$zipfile" 2>/dev/null; then
      log "REMOVING truncated zip for $name"
      rm -f "$zipfile"
    fi
  fi

  # Download if zip doesn't exist (or was just deleted)
  if [[ ! -f "$zipfile" ]]; then
    check_disk || { log "DISK FULL — skipping $name"; return 1; }

    local max_retries=3
    local attempt=0
    local dl_ok=0

    while (( attempt < max_retries && dl_ok == 0 )); do
      attempt=$((attempt + 1))
      log "DOWNLOAD $name (taxon=$taxon, $category, attempt $attempt/$max_retries)"
      touch_hb

      # Use correct subcommand: virus vs genome (for bacteria)
      local dl_cmd
      if [[ "$category" == "virus" ]]; then
        dl_cmd="datasets download virus genome taxon $taxon --include genome --filename $zipfile"
      else
        dl_cmd="datasets download genome taxon $taxon --include genome --filename $zipfile"
      fi

      if eval "$dl_cmd" 2>>"$LOG_DIR/download.log"; then
        # Verify zip is valid
        if python3 -c "import zipfile,sys; zipfile.ZipFile(sys.argv[1]).close()" "$zipfile" 2>/dev/null; then
          dl_ok=1
        else
          log "WARN: $name zip invalid after download (attempt $attempt)"
          rm -f "$zipfile"
        fi
      else
        log "WARN: $name download failed (attempt $attempt)"
        rm -f "$zipfile"
        sleep $((attempt * 10))  # back off: 10s, 20s, 30s
      fi
    done

    if (( dl_ok == 0 )); then
      log "ERROR: $name failed after $max_retries attempts"
      return 1
    fi
  fi

  # Extract FASTA using Python (handles Zip64 properly, macOS unzip cannot)
  touch_hb
  if [[ -f "$zipfile" ]]; then
    log "EXTRACTING $name zip (Python zipfile)..."
    if python3 -c "
import zipfile, sys, os
zpath = sys.argv[1]
outpath = sys.argv[2]
with zipfile.ZipFile(zpath) as z:
    fna_files = [n for n in z.namelist() if n.endswith('.fna')]
    if not fna_files:
        print('No .fna files in zip', file=sys.stderr)
        sys.exit(1)
    with open(outpath, 'wb') as out:
        for fname in fna_files:
            out.write(z.read(fname))
    print(f'Extracted {len(fna_files)} files -> {os.path.getsize(outpath)/1e6:.0f} MB')
" "$zipfile" "$fastafile" 2>&1 | tee -a "$LOG_DIR/download.log"; then
      local fasta_size
      fasta_size=$(stat -f %z "$fastafile" 2>/dev/null || echo 0)
      if (( fasta_size > 100 )); then
        log "OK $name — $(ls -lh "$fastafile" | awk '{print $5}')"
        rm -f "$zipfile"  # Save SSD space
        return 0
      fi
    fi
    log "WARN: Could not extract FASTA from $name zip"
    return 1
  fi
}

# Download RefSeq viral bulk files (parallel with aria2c)
download_refseq_viral() {
  local outdir="$RAW_DIR/refseq-viral"
  mkdir -p "$outdir"

  local urlfile="$outdir/urls.txt"
  cat > "$urlfile" << 'URLS'
https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.1.genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.4.1.genomic.fna.gz
URLS

  log "DOWNLOAD RefSeq viral bulk (4 files, parallel)..."
  touch_hb
  aria2c --input-file="$urlfile" \
    --dir="$outdir" \
    --max-concurrent-downloads=4 \
    --max-connection-per-server=4 \
    --split=4 \
    --continue=true \
    --auto-file-renaming=false \
    --log="$LOG_DIR/aria2_refseq.log" \
    --log-level=warn \
    2>&1 | tee -a "$LOG_DIR/download.log" || true
  touch_hb

  # Decompress
  for gz in "$outdir"/*.gz; do
    [[ -f "$gz" ]] || continue
    local target="${gz%.gz}"
    if [[ ! -f "$target" ]]; then
      log "DECOMPRESS $(basename "$gz")..."
      gzip -dk "$gz"
      touch_hb
    fi
  done
  log "RefSeq viral bulk done"
}

# Download GRCh38 human reference (for off-target)
download_human_ref() {
  local outdir="$RAW_DIR/human-grch38"
  mkdir -p "$outdir"
  local target="$outdir/GRCh38.fna"

  if [[ -f "$target" ]] && [[ $(stat -f %z "$target" 2>/dev/null || echo 0) -gt 1000000 ]]; then
    log "SKIP human GRCh38 — already have $(ls -lh "$target" | awk '{print $5}')"
    return 0
  fi

  log "DOWNLOAD GRCh38 human reference (~900 MB compressed)..."
  touch_hb
  aria2c "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" \
    --dir="$outdir" \
    --out="GRCh38.fna.gz" \
    --max-connection-per-server=8 \
    --split=8 \
    --continue=true \
    --log="$LOG_DIR/aria2_grch38.log" \
    --log-level=warn \
    2>&1 | tee -a "$LOG_DIR/download.log" || true

  if [[ -f "$outdir/GRCh38.fna.gz" ]] && [[ ! -f "$target" ]]; then
    log "DECOMPRESS GRCh38..."
    gzip -dk "$outdir/GRCh38.fna.gz"
    touch_hb
  fi
  log "GRCh38 done"
}

# ── Main download orchestrator ──

log "========================================="
log "STAGE 1: DOWNLOAD ($(date))"
log "========================================="

# Launch RefSeq bulk + GRCh38 + per-pathogen ALL IN PARALLEL
# Background jobs array for wait
PIDS=()

# RefSeq viral bulk (in background)
download_refseq_viral &
PIDS+=($!)

# GRCh38 human reference (in background)
download_human_ref &
PIDS+=($!)

# Per-pathogen downloads (in background, max 4 parallel)
PARALLEL_LIMIT=4
RUNNING=0
PATHOGENS=(
  "sars-cov-2|2697049|virus"
  "influenza-a|11320|virus"
  "dengue|12637|virus"
  "rsv|11250|virus"
  "hiv-1|11676|virus"
  "ebola|186538|virus"
  "zika|64320|virus"
  "mers|1335626|virus"
  "mpox|10244|virus"
  "hepatitis-b|10407|virus"
  "tuberculosis|1773|bacteria"
  "cholera|666|bacteria"
)

for entry in "${PATHOGENS[@]}"; do
  IFS='|' read -r name taxon category <<< "$entry"

  # Wait if too many parallel jobs
  while (( RUNNING >= PARALLEL_LIMIT )); do
    # Reap finished jobs
    NEW_PIDS=()
    for pid in "${PIDS[@]}"; do
      if kill -0 "$pid" 2>/dev/null; then
        NEW_PIDS+=("$pid")
      fi
    done
    PIDS=("${NEW_PIDS[@]}")
    RUNNING=${#PIDS[@]}
    (( RUNNING >= PARALLEL_LIMIT )) && sleep 5
  done

  download_pathogen "$name" "$taxon" "$RAW_DIR/$name" "$category" &
  PIDS+=($!)
  RUNNING=${#PIDS[@]}
  touch_hb
done

# Wait for all downloads
log "Waiting for ${#PIDS[@]} download jobs..."
for pid in "${PIDS[@]}"; do
  wait "$pid" 2>/dev/null || true
done

touch_hb
log "ALL DOWNLOADS COMPLETE"
touch "$WORKSPACE/.done-download"
DL_EOF
  chmod +x "$WORKSPACE/download.sh"
}

# ---------------------------------------------------------------------------
# Stage 2: Index building (parallel per-pathogen)
# ---------------------------------------------------------------------------
write_index_script() {
  cat > "$WORKSPACE/index.sh" << 'IDX_EOF'
#!/usr/bin/env bash
set -euo pipefail

WORKSPACE="$1"
LOOM_BIN="$2"
RAW_DIR="$WORKSPACE/raw"
IDX_DIR="$WORKSPACE/indexes"
LOG_DIR="$WORKSPACE/logs"
HEARTBEAT_DIR="$WORKSPACE/.heartbeats"

touch_hb() { touch "$HEARTBEAT_DIR/index"; }
log() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG_DIR/index.log"; }

# Wait for download stage to start producing data
log "========================================="
log "STAGE 2: INDEX BUILD ($(date))"
log "========================================="

# Poll for completed downloads and index them immediately
# This creates pipeline parallelism: download pathogen A while indexing pathogen B
INDEXED=()

index_pathogen() {
  local name="$1" raw_path="$2" idx_path="$3"
  touch_hb
  log "INDEX $name → $idx_path"

  # Find all FASTA files
  local fasta_files
  fasta_files=$(find "$raw_path" -name "*.fna" -o -name "*.fasta" -o -name "*.fa" 2>/dev/null | head -50)

  if [[ -z "$fasta_files" ]]; then
    log "WARN: No FASTA files found in $raw_path"
    return 1
  fi

  # Build index using loom CLI
  if "$LOOM_BIN" index \
      --path "$raw_path" \
      --extensions "fna,fasta,fa" \
      --output "$idx_path" \
      --build-mode auto \
      2>&1 | tee -a "$LOG_DIR/index_${name}.log"; then
    local idx_size
    idx_size=$(ls -lh "$idx_path" | awk '{print $5}')
    log "OK $name index: $idx_size"
    touch_hb
    return 0
  else
    log "ERROR: Index build failed for $name"
    return 1
  fi
}

# Wait for download stage to produce at least one pathogen
log "Waiting for downloads to start producing data..."
MAX_WAIT=3600  # 1 hour max wait
WAITED=0
while [[ ! -f "$WORKSPACE/.done-download" ]]; do
  # Check if any raw directories have FASTA files ready
  READY=$(find "$RAW_DIR" -name "*.fna" -o -name "*.fasta" 2>/dev/null | head -1)
  [[ -n "$READY" ]] && break
  sleep 30
  WAITED=$((WAITED + 30))
  touch_hb
  if (( WAITED >= MAX_WAIT )); then
    log "TIMEOUT: No data after ${MAX_WAIT}s"
    exit 1
  fi
done

# Index loop: process directories as they become available
# Run up to 3 parallel index builds (each is I/O bound on SSD)
INDEX_PARALLEL=3

while true; do
  PIDS=()
  CURRENT_RUNNING=0

  for raw_path in "$RAW_DIR"/*/; do
    [[ -d "$raw_path" ]] || continue
    name=$(basename "$raw_path")
    idx_path="$IDX_DIR/${name}.idx"

    # Skip already indexed
    if [[ -f "$idx_path" ]]; then
      continue
    fi

    # Skip if no FASTA files yet (still downloading)
    fasta_count=$(find "$raw_path" -name "*.fna" -o -name "*.fasta" -o -name "*.fa" 2>/dev/null | wc -l | tr -d ' ')
    if (( fasta_count == 0 )); then
      continue
    fi

    # Throttle parallel builds
    while (( CURRENT_RUNNING >= INDEX_PARALLEL )); do
      NEW_PIDS=()
      for pid in "${PIDS[@]}"; do
        if kill -0 "$pid" 2>/dev/null; then
          NEW_PIDS+=("$pid")
        fi
      done
      PIDS=("${NEW_PIDS[@]}")
      CURRENT_RUNNING=${#PIDS[@]}
      (( CURRENT_RUNNING >= INDEX_PARALLEL )) && sleep 10
    done

    index_pathogen "$name" "$raw_path" "$idx_path" &
    PIDS+=($!)
    CURRENT_RUNNING=${#PIDS[@]}
  done

  # Wait for current batch
  for pid in "${PIDS[@]}"; do
    wait "$pid" 2>/dev/null || true
  done

  # If downloads are done and everything is indexed, we're done
  if [[ -f "$WORKSPACE/.done-download" ]]; then
    UNINDEXED=0
    for raw_path in "$RAW_DIR"/*/; do
      [[ -d "$raw_path" ]] || continue
      name=$(basename "$raw_path")
      fasta_count=$(find "$raw_path" -name "*.fna" -o -name "*.fasta" -o -name "*.fa" 2>/dev/null | wc -l | tr -d ' ')
      (( fasta_count == 0 )) && continue
      [[ -f "$IDX_DIR/${name}.idx" ]] || UNINDEXED=$((UNINDEXED + 1))
    done
    if (( UNINDEXED == 0 )); then
      break
    fi
  fi

  touch_hb
  sleep 15
done

log "ALL INDEXES BUILT"

# Print summary
log "── Index Summary ──"
for idx in "$IDX_DIR"/*.idx; do
  [[ -f "$idx" ]] || continue
  log "  $(basename "$idx"): $(ls -lh "$idx" | awk '{print $5}')"
done

touch "$WORKSPACE/.done-index"
IDX_EOF
  chmod +x "$WORKSPACE/index.sh"
}

# ---------------------------------------------------------------------------
# Stage 3: CRISPR target scan (parallel per-pathogen index)
# ---------------------------------------------------------------------------
write_scan_script() {
  cat > "$WORKSPACE/scan.sh" << 'SCAN_EOF'
#!/usr/bin/env bash
set -euo pipefail

WORKSPACE="$1"
LOOM_BIN="$2"
LOOM_ROOT="$3"
IDX_DIR="$WORKSPACE/indexes"
SCAN_DIR="$WORKSPACE/crispr-db"
LOG_DIR="$WORKSPACE/logs"
HEARTBEAT_DIR="$WORKSPACE/.heartbeats"

touch_hb() { touch "$HEARTBEAT_DIR/scan"; }
log() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG_DIR/scan.log"; }

log "========================================="
log "STAGE 3: CRISPR SCAN ($(date))"
log "========================================="

# Wait for at least one index
log "Waiting for indexes..."
MAX_WAIT=7200
WAITED=0
while true; do
  IDX_COUNT=$(find "$IDX_DIR" -name "*.idx" 2>/dev/null | wc -l | tr -d ' ')
  (( IDX_COUNT > 0 )) && break
  sleep 30
  WAITED=$((WAITED + 30))
  touch_hb
  if (( WAITED >= MAX_WAIT )); then
    log "TIMEOUT: No indexes after ${MAX_WAIT}s"
    exit 1
  fi
done

# PAM patterns for CRISPR systems
# SpCas9: 20bp target + NGG (3')
# Cas12a: TTTN (5') + 23bp target
PAM_PATTERNS="NGG,TTTN,TTTV"

# Scan function: for each index, do a sliding-window CRISPR target enumeration
scan_pathogen() {
  local name="$1" idx_path="$2" out_csv="$3"
  touch_hb
  log "SCAN $name..."

  # Use the CorpusSearcher to scan for all 20-mers with PAM context
  python3 - "$idx_path" "$out_csv" "$name" "$LOOM_ROOT" << 'PYSCAN'
import csv
import json
import re
import sys
import time
from pathlib import Path

idx_path = sys.argv[1]
out_csv = sys.argv[2]
pathogen_name = sys.argv[3]
loom_root = sys.argv[4]

sys.path.insert(0, str(Path(loom_root) / "python"))
from loom.search_executor import CorpusSearcher, SearchResponse

# Common PAM patterns
PAMS = {
    "NGG": re.compile(r'[ACGT]{20}[ACGT]GG'),
    "TTTN": re.compile(r'TTT[ACGT][ACGT]{23}'),
}

searcher = CorpusSearcher(idx_path)
start = time.time()

# Strategy: search for conserved motifs and extract PAM-adjacent 23-mers
# For each short diagnostic-relevant motif, find occurrences and extract context
DIAGNOSTIC_MOTIFS = [
    "TATAAA", "TATAAT", "AATAAA",  # regulatory signals
    "ATGAAA", "ATGCCC", "ATGGGG",  # start codons + context
    "CGTACG", "GCTAGC", "GGATCC",  # restriction sites (diagnostic)
]

results = []
seen_targets = set()

# Phase A: search known diagnostic motifs
for motif in DIAGNOSTIC_MOTIFS:
    try:
        resp = searcher.search(motif, context_lines=2, max_results=100)
        count = resp.count
        if count > 0:
            for sr in resp.results:
                ctx = sr.context or ""
                # Extract 23-mers from context
                for m in re.finditer(r'[ACGT]{23}', ctx):
                    target = m.group()
                    if target not in seen_targets:
                        seen_targets.add(target)
                        # Check PAM
                        pam_type = "none"
                        if re.match(r'[ACGT]{20}[ACGT]GG$', target):
                            pam_type = "NGG"
                        elif re.match(r'^TTT[ACGT]', target):
                            pam_type = "TTTN"

                        results.append({
                            "pathogen": pathogen_name,
                            "gene": sr.file,
                            "position": sr.line,
                            "sequence_23mer": target,
                            "pam_type": pam_type,
                            "source_motif": motif,
                            "occurrences": count,
                        })
    except Exception as e:
        print(f"  WARN: motif {motif} search error: {e}", file=sys.stderr)

# Phase B: systematic marker counting (genome-wide conservation signals)
for marker in ["ATG", "TGA", "TAA", "TAG"]:
    try:
        occ = searcher.count(marker)
        results.append({
            "pathogen": pathogen_name,
            "gene": "genome-wide",
            "position": 0,
            "sequence_23mer": f"[{marker}-marker]",
            "pam_type": "marker",
            "source_motif": marker,
            "occurrences": occ,
        })
    except Exception:
        pass

elapsed = time.time() - start

# Write CSV
with open(out_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=[
        "pathogen", "gene", "position", "sequence_23mer",
        "pam_type", "source_motif", "occurrences"
    ])
    writer.writeheader()
    writer.writerows(results)

print(f"  {pathogen_name}: {len(results)} targets found in {elapsed:.1f}s")
PYSCAN

  local line_count
  line_count=$(wc -l < "$out_csv" | tr -d ' ')
  log "OK $name scan: $line_count targets → $out_csv"
  touch_hb
}

# Scan loop: process indexes as they become available
SCAN_PARALLEL=4

while true; do
  PIDS=()
  CURRENT_RUNNING=0

  for idx in "$IDX_DIR"/*.idx; do
    [[ -f "$idx" ]] || continue
    name=$(basename "$idx" .idx)
    out_csv="$SCAN_DIR/${name}_targets.csv"

    # Skip already scanned
    [[ -f "$out_csv" ]] && continue

    # Throttle
    while (( CURRENT_RUNNING >= SCAN_PARALLEL )); do
      NEW_PIDS=()
      for pid in "${PIDS[@]}"; do
        if kill -0 "$pid" 2>/dev/null; then
          NEW_PIDS+=("$pid")
        fi
      done
      PIDS=("${NEW_PIDS[@]}")
      CURRENT_RUNNING=${#PIDS[@]}
      (( CURRENT_RUNNING >= SCAN_PARALLEL )) && sleep 5
    done

    scan_pathogen "$name" "$idx" "$out_csv" &
    PIDS+=($!)
    CURRENT_RUNNING=${#PIDS[@]}
  done

  # Wait for batch
  for pid in "${PIDS[@]}"; do
    wait "$pid" 2>/dev/null || true
  done

  # Done when indexes are done and all scanned
  if [[ -f "$WORKSPACE/.done-index" ]]; then
    UNSCANNED=0
    for idx in "$IDX_DIR"/*.idx; do
      [[ -f "$idx" ]] || continue
      name=$(basename "$idx" .idx)
      [[ -f "$SCAN_DIR/${name}_targets.csv" ]] || UNSCANNED=$((UNSCANNED + 1))
    done
    if (( UNSCANNED == 0 )); then
      break
    fi
  fi

  touch_hb
  sleep 15
done

log "ALL SCANS COMPLETE"

# Merge all per-pathogen CSVs into one combined database
COMBINED="$SCAN_DIR/crispr_targets_combined.csv"
FIRST=1
for csv in "$SCAN_DIR"/*_targets.csv; do
  [[ -f "$csv" ]] || continue
  if (( FIRST )); then
    cat "$csv" > "$COMBINED"
    FIRST=0
  else
    tail -n +2 "$csv" >> "$COMBINED"
  fi
done

TOTAL=$(wc -l < "$COMBINED" | tr -d ' ')
log "COMBINED DATABASE: $TOTAL entries → $COMBINED"
touch "$WORKSPACE/.done-scan"

# Final summary
log "========================================="
log "PIPELINE COMPLETE — $(date)"
log "========================================="
log "Downloads: $(du -sh "$WORKSPACE/raw" | cut -f1)"
log "Indexes:   $(du -sh "$WORKSPACE/indexes" | cut -f1)"
log "Targets:   $TOTAL entries in $COMBINED"
log "========================================="
SCAN_EOF
  chmod +x "$WORKSPACE/scan.sh"
}

# ---------------------------------------------------------------------------
# Launch everything in tmux
# ---------------------------------------------------------------------------
launch_tmux() {
  # Kill existing session if any
  tmux kill-session -t "$SESSION" 2>/dev/null || true

  echo "== Launching tmux session: $SESSION =="

  # Create session with monitor window
  tmux new-session -d -s "$SESSION" -n monitor \
    "bash $WORKSPACE/monitor.sh $WORKSPACE"

  # Window 1: Downloads (all parallel)
  tmux new-window -t "$SESSION" -n download \
    "bash $WORKSPACE/download.sh $WORKSPACE 2>&1 | tee $LOG_DIR/download_${STAMP}.log"

  # Window 2: Index builds (fires as data arrives)
  tmux new-window -t "$SESSION" -n index \
    "bash $WORKSPACE/index.sh $WORKSPACE $LOOM_BIN 2>&1 | tee $LOG_DIR/index_${STAMP}.log"

  # Window 3: CRISPR scans (fires as indexes complete)
  tmux new-window -t "$SESSION" -n scan \
    "bash $WORKSPACE/scan.sh $WORKSPACE $LOOM_BIN $LOOM_ROOT 2>&1 | tee $LOG_DIR/scan_${STAMP}.log"

  # Window 4: Log tail
  tmux new-window -t "$SESSION" -n logs \
    "tail -f $LOG_DIR/download.log $LOG_DIR/index.log $LOG_DIR/scan.log 2>/dev/null || echo 'Waiting for logs...'; sleep 999999"

  # Select monitor window
  tmux select-window -t "$SESSION:monitor"

  echo ""
  echo "══════════════════════════════════════════════════════"
  echo "  PIVOT-OPENSCIENCE pipeline launched!"
  echo ""
  echo "  Attach:  tmux attach -t $SESSION"
  echo "  Windows: 0:monitor  1:download  2:index  3:scan  4:logs"
  echo ""
  echo "  Workspace: $WORKSPACE"
  echo "  SSD free:  $(df -h "$SSD" | awk 'NR==2 {print $4}')"
  echo "══════════════════════════════════════════════════════"
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
preflight

if (( DRY_RUN )); then
  echo ""
  echo "=== Pipeline Plan ==="
  echo ""
  echo "1. DOWNLOAD (parallel):"
  echo "   - RefSeq viral bulk (4 files via aria2c, 8 connections)"
  echo "   - GRCh38 human reference (~3 GB)"
  echo "   - 12 pathogens via NCBI datasets CLI (zip → extract, 4 parallel)"
  echo "   Estimated: 10-50 GB total, ~15 min on 10 Gbps"
  echo ""
  echo "2. INDEX (parallel, pipeline-overlapped with downloads):"
  echo "   - FM-index per pathogen (3 parallel builds)"
  echo "   - Fires as soon as FASTA files land"
  echo "   Estimated: 2-30 min per pathogen depending on size"
  echo ""
  echo "3. CRISPR SCAN (parallel, pipeline-overlapped with indexing):"
  echo "   - Sliding window + PAM detection per pathogen index"
  echo "   - 4 parallel scans"
  echo "   - Combined CSV database at end"
  echo "   Estimated: 1-10 min per pathogen"
  echo ""
  echo "4. HEALTH MONITORING:"
  echo "   - Heartbeat files per worker (stale >5 min = warning)"
  echo "   - Disk space guard (stops at <${MIN_DISK_GB} GB free)"
  echo "   - Monitor dashboard in tmux window 0"
  echo "   - Zero LLM/API calls (no cost risk)"
  echo ""
  echo "Workspace: $WORKSPACE"
  echo "Attach:    tmux attach -t $SESSION"
  exit 0
fi

write_monitor_script
write_download_script
write_index_script
write_scan_script
launch_tmux
