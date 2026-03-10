#!/usr/bin/env bash
# ===========================================================================
# Ontology Enrichment Overnight Pipeline
# ===========================================================================
#
# Fetches biomedical ontology data from 5 sources in parallel,
# assembles into web/data/ontology-enrichment.json for offline use.
#
# Usage:
#   bash scripts/ontology/overnight_ontology.sh
#
# Architecture:
#   tmux session "ontology" with windows:
#     0:monitor   — watches logs + shows progress
#     1:taxonomy   — NCBI Taxonomy fetcher
#     2:disease    — Disease Ontology fetcher
#     3:mondo      — MONDO fetcher
#     4:who        — WHO context fetcher
#     5:genes      — Gene annotation fetcher
#     6:assemble   — waits for all, then assembles final JSON
#
# All fetchers are independent → full parallelism.
# Assembly waits for all fetchers to complete.
# ===========================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOOM_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CACHE_DIR="$SCRIPT_DIR/ontology_cache"
LOG_DIR="$CACHE_DIR/logs"
SESSION="ontology"
STAMP="$(date +%Y%m%d-%H%M%S)"

mkdir -p "$CACHE_DIR" "$LOG_DIR"

echo "============================================"
echo " LOOM Ontology Enrichment Pipeline"
echo " Started: $(date)"
echo " Session: tmux attach -t $SESSION"
echo "============================================"

# Kill existing session if any
tmux kill-session -t "$SESSION" 2>/dev/null || true

# Create session with monitor window
tmux new-session -d -s "$SESSION" -n monitor

# ---------------------------------------------------------------------------
# Window 0: Monitor — watches all log files
# ---------------------------------------------------------------------------
tmux send-keys -t "$SESSION:monitor" "echo '=== Ontology Pipeline Monitor ===' && echo 'Waiting for fetchers to start...' && sleep 3 && tail -f $LOG_DIR/*.log 2>/dev/null || echo 'Logs will appear as fetchers start...'" Enter

# ---------------------------------------------------------------------------
# Window 1: NCBI Taxonomy (parallel)
# ---------------------------------------------------------------------------
tmux new-window -t "$SESSION" -n taxonomy
tmux send-keys -t "$SESSION:taxonomy" "cd $LOOM_ROOT && python3 scripts/ontology/fetch_ncbi_taxonomy.py 2>&1 | tee $LOG_DIR/taxonomy_${STAMP}.log && echo '=== TAXONOMY DONE ===' && touch $CACHE_DIR/.done_taxonomy" Enter

# ---------------------------------------------------------------------------
# Window 2: Disease Ontology (parallel)
# ---------------------------------------------------------------------------
tmux new-window -t "$SESSION" -n disease
tmux send-keys -t "$SESSION:disease" "cd $LOOM_ROOT && python3 scripts/ontology/fetch_disease_ontology.py 2>&1 | tee $LOG_DIR/disease_${STAMP}.log && echo '=== DISEASE ONTOLOGY DONE ===' && touch $CACHE_DIR/.done_disease" Enter

# ---------------------------------------------------------------------------
# Window 3: MONDO (parallel)
# ---------------------------------------------------------------------------
tmux new-window -t "$SESSION" -n mondo
tmux send-keys -t "$SESSION:mondo" "cd $LOOM_ROOT && python3 scripts/ontology/fetch_mondo.py 2>&1 | tee $LOG_DIR/mondo_${STAMP}.log && echo '=== MONDO DONE ===' && touch $CACHE_DIR/.done_mondo" Enter

# ---------------------------------------------------------------------------
# Window 4: WHO Context (parallel)
# ---------------------------------------------------------------------------
tmux new-window -t "$SESSION" -n who
tmux send-keys -t "$SESSION:who" "cd $LOOM_ROOT && python3 scripts/ontology/fetch_who_context.py 2>&1 | tee $LOG_DIR/who_${STAMP}.log && echo '=== WHO DONE ===' && touch $CACHE_DIR/.done_who" Enter

# ---------------------------------------------------------------------------
# Window 5: Gene annotations (parallel — this one takes longest)
# ---------------------------------------------------------------------------
tmux new-window -t "$SESSION" -n genes
tmux send-keys -t "$SESSION:genes" "cd $LOOM_ROOT && python3 scripts/ontology/fetch_gene_annotations.py 2>&1 | tee $LOG_DIR/genes_${STAMP}.log && echo '=== GENES DONE ===' && touch $CACHE_DIR/.done_genes" Enter

# ---------------------------------------------------------------------------
# Window 6: Assemble — waits for all fetchers, then builds final JSON
# ---------------------------------------------------------------------------
tmux new-window -t "$SESSION" -n assemble
tmux send-keys -t "$SESSION:assemble" "cd $LOOM_ROOT && echo 'Waiting for all fetchers to complete...' && while [ ! -f $CACHE_DIR/.done_taxonomy ] || [ ! -f $CACHE_DIR/.done_disease ] || [ ! -f $CACHE_DIR/.done_mondo ] || [ ! -f $CACHE_DIR/.done_who ] || [ ! -f $CACHE_DIR/.done_genes ]; do printf '.'; sleep 5; done && echo '' && echo 'All fetchers done! Assembling...' && python3 scripts/ontology/assemble_enrichment.py 2>&1 | tee $LOG_DIR/assemble_${STAMP}.log && echo '' && echo '============================================' && echo ' ONTOLOGY ENRICHMENT COMPLETE' && echo ' Output: web/data/ontology-enrichment.json' && echo ' $(date)' && echo '============================================'" Enter

# Clean up sentinel files from previous runs
rm -f "$CACHE_DIR"/.done_*

echo ""
echo "Pipeline launched in tmux session '$SESSION'"
echo ""
echo "  tmux attach -t $SESSION         # watch progress"
echo "  tmux attach -t $SESSION:monitor  # tail all logs"
echo "  tmux attach -t $SESSION:genes    # watch slowest fetcher"
echo "  tmux attach -t $SESSION:assemble # watch final assembly"
echo ""
echo "Fetchers running in parallel:"
echo "  1. NCBI Taxonomy (16 organisms)"
echo "  2. Disease Ontology (12 diseases)"
echo "  3. MONDO (12 diseases)"
echo "  4. WHO context (curated + GHO API)"
echo "  5. Gene annotations (12 genomes — slowest)"
echo ""
echo "Assembly will auto-run when all fetchers complete."
echo "Expected total time: 5-15 minutes."
echo ""
echo "Output: web/data/ontology-enrichment.json"
