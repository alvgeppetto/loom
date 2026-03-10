#!/usr/bin/env bash
# Cross-reactivity screen for Cas12a (TTTN) guides against animal host indexes.
# Mirrors compute_cross_reactivity.sh but for TTTN PAM instead of NGG.
# Outputs web/data/cross-reactivity-cas12a.json
set -euo pipefail

# ── Memory-safety: cap RSS to 80% of physical RAM so we never OOM ──────────
if command -v sysctl &>/dev/null; then
  TOTAL_RAM_KB=$(( $(sysctl -n hw.memsize) / 1024 ))
else
  TOTAL_RAM_KB=$(awk '/MemTotal/{print $2}' /proc/meminfo)
fi
RAM_LIMIT_KB=$(( TOTAL_RAM_KB * 80 / 100 ))
ulimit -v "$RAM_LIMIT_KB" 2>/dev/null || true
echo "[cas12a-cr] RSS cap: ${RAM_LIMIT_KB} KB (80% of ${TOTAL_RAM_KB} KB)"

LOOM_BIN="./target/release/loom"
CRISPR_DB="/Volumes/T9/loom-openscience/crispr-db"
IDX_DIR="/Volumes/T9/loom-openscience/indexes"
OUT="web/data/cross-reactivity-cas12a.json"
TMP_DIR="/tmp/cross-reactivity-cas12a-$$"

# 11 non-SARS-CoV-2 pathogens (same as NGG screen)
PATHOGENS=(cholera dengue ebola hepatitis-b hiv-1 influenza-a mers mpox rsv tuberculosis zika)
ANIMALS=(human-grch38 pig bat-rousettus chicken cow camel mouse)

mkdir -p "$TMP_DIR"

echo "[cas12a-cr] Extracting top TTTN guides per pathogen..."
for p in "${PATHOGENS[@]}"; do
  csv="$CRISPR_DB/${p}_targets.csv"
  if [[ ! -f "$csv" ]]; then
    echo "  WARN: $csv not found, skipping $p"
    continue
  fi
  # Extract top 100 TTTV (Cas12a) guides by occurrence count (most conserved first)
  # Filter: pam_type=="TTTN" AND 4th base of sequence is NOT T (i.e., TTTV where V=A/C/G)
  awk -F, 'NR>1 && $5=="TTTN" && substr($4,4,1)!="T" {print $7","$4}' "$csv" | sort -t, -k1 -rn | head -100 | cut -d, -f2 > "$TMP_DIR/${p}_tttn.txt" || true
  n=$(wc -l < "$TMP_DIR/${p}_tttn.txt" | tr -d ' ')
  echo "  $p: $n TTTV guides"
done

# Combine all unique guides into one file for batch search
sort -u "$TMP_DIR"/*_tttn.txt > "$TMP_DIR/all_guides.txt"
TOTAL=$(wc -l < "$TMP_DIR/all_guides.txt" | tr -d ' ')
echo "[cas12a-cr] Total unique TTTV guides: $TOTAL"

# Search against each animal index (sequential for memory safety)
for animal in "${ANIMALS[@]}"; do
  idx="$IDX_DIR/${animal}.idx"
  if [[ ! -f "$idx" ]]; then
    echo "  WARN: $idx not found, skipping $animal"
    continue
  fi
  idx_size=$(ls -lh "$idx" | awk '{print $5}')
  echo "[cas12a-cr] Searching $TOTAL guides against $animal ($idx_size)..."
  nice -n 10 $LOOM_BIN batch-search --index "$idx" --patterns "$TMP_DIR/all_guides.txt" \
    > "$TMP_DIR/hits_${animal}.csv" 2>"$TMP_DIR/log_${animal}.txt"
  echo "  $(grep 'completed' "$TMP_DIR/log_${animal}.txt" || true)"
  echo "  done"
done

# Assemble JSON
echo "[cas12a-cr] Assembling JSON..."
python3 - "$TMP_DIR" "$OUT" "${PATHOGENS[*]}" "${ANIMALS[*]}" <<'PYEOF'
import sys, json, csv
from pathlib import Path
from datetime import datetime, timezone

tmp_dir = Path(sys.argv[1])
out_path = sys.argv[2]
pathogens = sys.argv[3].split()
animals = sys.argv[4].split()

# Load per-animal hit counts
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
    gfile = tmp_dir / f"{p}_tttn.txt"
    if not gfile.exists():
        continue
    with open(gfile) as f:
        pathogen_guides[p] = [line.strip() for line in f if line.strip()]

# Build result
result = {
    "generated": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    "pam_type": "TTTV (Cas12a)",
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
            "specific": total_host_hits == 0
        })
    result["pathogens"][p] = p_data

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

print(f"[cas12a-cr] Done: {total_guides} guides, {specific_guides} specific ({result['summary']['specificity_rate']}%)")
print(f"[cas12a-cr] Output: {out_path}")
PYEOF

rm -rf "$TMP_DIR"
echo "[cas12a-cr] Complete."
