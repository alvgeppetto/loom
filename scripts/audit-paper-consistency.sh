#!/usr/bin/env bash
set -euo pipefail

# Consistency sweep for denominator-sensitive SARS-CoV-2 manuscript claims.
#
# Usage:
#   scripts/audit-paper-consistency.sh            # report mode
#   scripts/audit-paper-consistency.sh --strict   # fail on mixed-era claims
#
# Canonical same-corpus denominator (March 2026 rerun): 9,193,298.
# Legacy denominator token set: 9,419,528 / 9.4 million / 9419528.

STRICT=0
if [[ "${1:-}" == "--strict" ]]; then
  STRICT=1
fi

PAPER_DIR="publications/papers"
CANON_RE='9,193,298|9193298|9\.2 million'
LEGACY_RE='9,419,528|9419528|9\.4 million'

status=0

printf "\n[paper-consistency] Sweeping %s\n" "$PAPER_DIR"
printf "%-42s  %7s  %7s  %s\n" "file" "canon" "legacy" "status"
printf "%-42s  %7s  %7s  %s\n" "------------------------------------------" "-------" "-------" "------"

for f in "$PAPER_DIR"/*.md; do
  base=$(basename "$f")
  canon_count=$( (grep -Eo "$CANON_RE" "$f" 2>/dev/null || true) | wc -l | tr -d ' ')
  legacy_count=$( (grep -Eo "$LEGACY_RE" "$f" 2>/dev/null || true) | wc -l | tr -d ' ')

  file_status="OK"

  if [[ "$base" == "pangenomic-crispr-targets-merged.md" ]]; then
    if [[ "$legacy_count" -gt 0 ]]; then
      file_status="FAIL:merged-has-legacy"
      status=1
    elif [[ "$canon_count" -eq 0 ]]; then
      file_status="FAIL:merged-missing-canon"
      status=1
    fi
  else
    if [[ "$canon_count" -gt 0 && "$legacy_count" -gt 0 ]]; then
      file_status="WARN:mixed-era"
      if [[ "$STRICT" -eq 1 ]]; then
        file_status="FAIL:mixed-era"
        status=1
      fi
    elif [[ "$legacy_count" -gt 0 ]]; then
      file_status="WARN:legacy-only"
    fi
  fi

  printf "%-42s  %7s  %7s  %s\n" "$base" "$canon_count" "$legacy_count" "$file_status"
done

if [[ "$STRICT" -eq 1 && "$status" -ne 0 ]]; then
  echo "\n[paper-consistency] STRICT FAILED"
  exit 1
fi

echo "\n[paper-consistency] Done"
