#!/usr/bin/env bash
# ==========================================================================
# Build PDF preprints from markdown papers via pandoc + LaTeX.
#
# Usage:
#   ./scripts/build-preprints.sh            # Build all papers
#   ./scripts/build-preprints.sh <paper>    # Build one (e.g. pan-pathogen)
#
# Output: publications/preprints/*.pdf
# ==========================================================================
set -euo pipefail

PAPER_DIR="publications/papers"
OUT_DIR="publications/preprints"
mkdir -p "$OUT_DIR"

# ---------------------------------------------------------------------------
# Pre-build gate: citation integrity check (structural only — no network)
# Run with SKIP_CITE_CHECK=1 to bypass (e.g. in CI without internet)
# ---------------------------------------------------------------------------
PYTHON="${PYTHON:-python3}"

# Pandoc settings for academic preprint
PANDOC_OPTS=(
  --pdf-engine=xelatex
  --variable=geometry:margin=2.5cm
  --variable=fontsize:11pt
  --variable=documentclass:article
  --variable=linestretch:1.15
  --variable=linkcolor:blue
  --variable=urlcolor:blue
  --variable=header-includes:'\usepackage{booktabs}\usepackage{longtable}\usepackage{float}\usepackage{fancyhdr}\usepackage{adjustbox}\usepackage{newunicodechar}\usepackage{amsmath}\usepackage{amssymb}\usepackage{microtype}\usepackage{xurl}\newunicodechar{≥}{$\geq$}\newunicodechar{≤}{$\leq$}\newunicodechar{≠}{$\neq$}\newunicodechar{→}{$\rightarrow$}\pagestyle{fancy}\fancyhead[L]{\small\textit{PREPRINT --- Not peer reviewed}}\fancyhead[R]{\small March 2026}\fancyfoot[C]{\thepage}\setlength{\tabcolsep}{3pt}\renewcommand{\arraystretch}{1.2}\sloppy'
  --lua-filter=scripts/fit-tables.lua
  --columns=9999
  --number-sections
  --standalone
)

check_regressions() {
  local md_file="$1"
  local name
  name=$(basename "$md_file" .md)

  # Prevent known denominator regressions in the merged SARS-CoV-2 manuscript.
  # §2.5 intentionally mentions the discovery denominator (9,419,528) for provenance;
  # allow up to 2 occurrences but flag if it appears more broadly.
  if [[ "$name" == "pangenomic-crispr-targets-merged" ]]; then
    local stale_count
    stale_count=$(grep -cE '9,419,528|9419528' "$md_file" || true)
    if (( stale_count > 2 )); then
      echo ""
      echo "ERROR: Regression guard failed for $name."
      echo "Found $stale_count mentions of old denominator (max 2 provenance allowed)."
      exit 1
    fi
    if grep -Eq '9\.4 million genomes|stricter header-level genome deduplication' "$md_file"; then
      echo ""
      echo "ERROR: Regression guard failed for $name."
      echo "Found stale denominator/dedup phrasing in $md_file."
      echo "Expected same-corpus denominator language with 9,193,298 genomes."
      exit 1
    fi
    if ! grep -q '9,193,298' "$md_file"; then
      echo ""
      echo "ERROR: Regression guard failed for $name."
      echo "Expected denominator token '9,193,298' was not found in $md_file."
      exit 1
    fi
  fi
}

build_paper() {
  local md_file="$1"
  local name
  name=$(basename "$md_file" .md)
  local out_pdf="$OUT_DIR/${name}.pdf"

  # Per-paper citation check gate
  if [[ "${SKIP_CITE_CHECK:-0}" != "1" ]]; then
    echo "--- Citation check: $name ---"
    if ! "$PYTHON" scripts/verify-citations.py --offline "$md_file"; then
      echo ""
      echo "ERROR: Citation check FAILED for $name. Fix above issues before building."
      echo "       To bypass: SKIP_CITE_CHECK=1 ./scripts/build-preprints.sh"
      exit 1
    fi
    echo ""
  fi

  check_regressions "$md_file"

  # Artifact-backed claim verification gate
  if [[ "${SKIP_CLAIM_CHECK:-0}" != "1" ]]; then
    echo "--- Claim verification: $name ---"
    if ! "$PYTHON" scripts/verify-paper-claims.py "$md_file"; then
      echo ""
      echo "ERROR: Claim verification FAILED for $name. Fix above issues before building."
      echo "       To bypass: SKIP_CLAIM_CHECK=1 ./scripts/build-preprints.sh"
      exit 1
    fi
  fi

  echo "Building $name → $out_pdf"
  pandoc "$md_file" "${PANDOC_OPTS[@]}" -o "$out_pdf" 2>&1
  echo "  OK: $(du -h "$out_pdf" | cut -f1)"

  # Overfull hbox detection gate — prevents margin overflow regressions.
  # Generate TeX, compile with xelatex, and fail if any overfull warnings.
  if [[ "${SKIP_OVERFLOW_CHECK:-0}" != "1" ]] && command -v xelatex >/dev/null 2>&1; then
    echo "--- Overflow check: $name ---"
    local tmpdir
    tmpdir=$(mktemp -d)
    local tex_file="$tmpdir/paper.tex"
    pandoc "$md_file" "${PANDOC_OPTS[@]}" -t latex -o "$tex_file" 2>/dev/null
    local overflows
    overflows=$(cd "$tmpdir" && xelatex -interaction=nonstopmode paper.tex 2>&1 | grep -c 'Overfull \\hbox' || true)
    rm -rf "$tmpdir"
    if (( overflows > 0 )); then
      echo ""
      echo "ERROR: Overflow check FAILED for $name."
      echo "Found $overflows overfull hbox warning(s). Fix layout before publishing."
      echo "       To bypass: SKIP_OVERFLOW_CHECK=1 ./scripts/build-preprints.sh"
      exit 1
    fi
    echo "  No overfull hbox warnings."
  fi
}

if [[ $# -gt 0 ]]; then
  # Build specific paper
  target="$PAPER_DIR/$1.md"
  if [[ ! -f "$target" ]]; then
    # Try with fuzzy matching
    target=$(find "$PAPER_DIR" -name "*$1*.md" | head -1)
  fi
  if [[ -z "$target" ]] || [[ ! -f "$target" ]]; then
    echo "ERROR: Paper not found: $1"
    echo "Available papers:"
    ls "$PAPER_DIR"/*.md 2>/dev/null | sed 's/.*\//  /' | sed 's/\.md$//'
    exit 1
  fi
  build_paper "$target"
else
  # Build the active paper (merged manuscript is the single source of truth)
  build_paper "$PAPER_DIR/pangenomic-crispr-targets-merged.md"
  # Build the brief (no claim-check or regression gates — it's a condensed derivative)
  if [[ -f "$PAPER_DIR/pangenomic-crispr-targets-brief.md" ]]; then
    SKIP_CITE_CHECK=1 SKIP_CLAIM_CHECK=1 build_paper "$PAPER_DIR/pangenomic-crispr-targets-brief.md"
  fi
fi

echo ""
echo "=== Preprints ==="
ls -lh "$OUT_DIR"/*.pdf 2>/dev/null
