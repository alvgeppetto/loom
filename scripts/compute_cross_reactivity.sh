#!/usr/bin/env bash
# Cross-reactivity matrix: search top NGG pathogen guides against animal host indexes.
# Outputs web/data/cross-reactivity.json
#
# Usage: bash scripts/compute_cross_reactivity.sh
set -euo pipefail

# ── Memory-safety: cap RSS to 80% of physical RAM so we never OOM ──────────
# macOS: use sysctl; Linux: use /proc/meminfo
if command -v sysctl &>/dev/null; then
  TOTAL_RAM_KB=$(( $(sysctl -n hw.memsize) / 1024 ))
else
  TOTAL_RAM_KB=$(awk '/MemTotal/{print $2}' /proc/meminfo)
fi
RAM_LIMIT_KB=$(( TOTAL_RAM_KB * 80 / 100 ))
ulimit -v "$RAM_LIMIT_KB" 2>/dev/null || true   # soft RLIMIT_AS
echo "[cross-reactivity] RSS cap: ${RAM_LIMIT_KB} KB (80% of ${TOTAL_RAM_KB} KB)"

LOOM_BIN="./target/release/loom"
CRISPR_DB="/Volumes/T9/loom-openscience/crispr-db"
IDX_DIR="/Volumes/T9/loom-openscience/indexes"
OUT="web/data/cross-reactivity.json"
TMP_DIR="/tmp/cross-reactivity-$$"

PATHOGENS=(cholera dengue ebola hepatitis-b hiv-1 influenza-a mers mpox rsv sars-cov-2 tuberculosis zika)
ANIMALS=(human-grch38 pig bat-rousettus chicken cow camel mouse)

mkdir -p "$TMP_DIR"

echo "[cross-reactivity] Extracting top NGG guides per pathogen..."
for p in "${PATHOGENS[@]}"; do
  csv="$CRISPR_DB/${p}_targets.csv"
  if [[ ! -f "$csv" ]]; then
    echo "  WARN: $csv not found, skipping $p"
    continue
  fi
  # Extract top 100 NGG guides by occurrence count (most conserved first)
  # Use awk to filter + emit, then sort, avoiding SIGPIPE with pipefail
  awk -F, 'NR>1 && $5=="NGG" {print $7","$4}' "$csv" | sort -t, -k1 -rn | head -100 | cut -d, -f2 > "$TMP_DIR/${p}_ngg.txt" || true
  n=$(wc -l < "$TMP_DIR/${p}_ngg.txt" | tr -d ' ')
  echo "  $p: $n NGG guides"
done

# Combine all unique guides into one file for batch search
sort -u "$TMP_DIR"/*_ngg.txt > "$TMP_DIR/all_guides.txt"
TOTAL=$(wc -l < "$TMP_DIR/all_guides.txt" | tr -d ' ')
echo "[cross-reactivity] Total unique guides: $TOTAL"

# Search all guides against each animal index — ONE AT A TIME for memory safety.
# Each `batch-search` loads the full FM-index into RAM (1–3 GB per animal).
# Sequential execution ensures only one index is resident at a time.
for animal in "${ANIMALS[@]}"; do
  idx="$IDX_DIR/${animal}.idx"
  if [[ ! -f "$idx" ]]; then
    echo "  WARN: $idx not found, skipping $animal"
    continue
  fi
  idx_size=$(ls -lh "$idx" | awk '{print $5}')
  echo "[cross-reactivity] Searching $TOTAL guides against $animal ($idx_size)..."
  # nice -n 10: lower scheduler priority so other processes aren't starved
  nice -n 10 $LOOM_BIN batch-search --index "$idx" --patterns "$TMP_DIR/all_guides.txt" \
    > "$TMP_DIR/hits_${animal}.csv" 2>"$TMP_DIR/log_${animal}.txt"
  echo "  $(grep 'completed' "$TMP_DIR/log_${animal}.txt" || true)"
  echo "  done — releasing index memory before next animal"
done

# Assemble JSON with Python
echo "[cross-reactivity] Assembling JSON..."
python3 - "$TMP_DIR" "$OUT" "${PATHOGENS[*]}" "${ANIMALS[*]}" <<'PYEOF'
import sys, json, csv, os
from pathlib import Path

tmp_dir = Path(sys.argv[1])
out_path = sys.argv[2]
pathogens = sys.argv[3].split()
animals = sys.argv[4].split()

# Load per-animal hit counts: guide -> count
animal_hits = {}
for animal in animals:
    hits_file = tmp_dir / f"hits_{animal}.csv"
    if not hits_file.exists():
        continue
    counts = {}
    with open(hits_file) as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            if len(row) >= 2:
                counts[row[0]] = int(row[1])
    animal_hits[animal] = counts

# Load per-pathogen guide lists
pathogen_guides = {}
for p in pathogens:
    gfile = tmp_dir / f"{p}_ngg.txt"
    if not gfile.exists():
        continue
    with open(gfile) as f:
        pathogen_guides[p] = [line.strip() for line in f if line.strip()]

# Build result: for each pathogen, for each guide, show hits per animal
result = {
    "generated": __import__('datetime').datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
    "animals": animals,
    "pathogens": {}
}

for p, guides in pathogen_guides.items():
    p_data = []
    for g in guides:
        hits = {}
        for animal in animals:
            hits[animal] = animal_hits.get(animal, {}).get(g, 0)
        total_host_hits = sum(hits.values())
        p_data.append({
            "seq": g,
            "host_hits": hits,
            "total_host_hits": total_host_hits,
            "specific": total_host_hits == 0  # no hits in any host = specific
        })
    result["pathogens"][p] = p_data

# Summary stats
total_guides = sum(len(v) for v in pathogen_guides.values())
specific_guides = sum(1 for p in result["pathogens"].values() for g in p if g["specific"])
result["summary"] = {
    "total_guides_checked": total_guides,
    "specific_guides": specific_guides,
    "cross_reactive_guides": total_guides - specific_guides,
    "specificity_rate": round(specific_guides / total_guides * 100, 1) if total_guides else 0
}

with open(out_path, 'w') as f:
    json.dump(result, f, indent=2)

print(f"[cross-reactivity] Done: {total_guides} guides, {specific_guides} specific ({result['summary']['specificity_rate']}%)")
print(f"[cross-reactivity] Output: {out_path}")
PYEOF

rm -rf "$TMP_DIR"
echo "[cross-reactivity] Complete."
