#!/usr/bin/env bash
# ===========================================================================
# Overnight Ontology ↔ Targets Integration Pipeline
# ===========================================================================
#
# Three parallel worktrees, each tackling one tier of the integration:
#
#   Tier 1 (gene-annotate):  Positional join — map every CRISPR target to its
#                            gene using [start,end] coordinate overlap.
#
#   Tier 2 (dx-priority):   Diagnostic priority scoring — rank pathogens by
#                            unmet need (no CRISPR dx + high burden).
#
#   Tier 3 (ui-integration): Surface gene annotations + priority scores in
#                            the web UI (target table, filters, top-opps view).
#
# Usage:
#   bash scripts/overnight_ontology_integration.sh          # launch
#   bash scripts/overnight_ontology_integration.sh --dry-run # plan only
#   tmux attach -t onto-integrate                           # monitor
#
# ===========================================================================

set -euo pipefail

LOOM_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SESSION="onto-integrate"
STAMP="$(date +%Y%m%d-%H%M%S)"
LOG_DIR="$LOOM_ROOT/logs/onto-integrate-$STAMP"
DRY_RUN=0

BRANCH_T1="feature/gene-annotate-targets"
BRANCH_T2="feature/diagnostic-priority"
BRANCH_T3="feature/ui-gene-integration"

WORKTREE_T1="$LOOM_ROOT/../loom-gene-annotate"
WORKTREE_T2="$LOOM_ROOT/../loom-dx-priority"
WORKTREE_T3="$LOOM_ROOT/../loom-ui-gene"

if [[ "${1:-}" == "--dry-run" ]]; then
  DRY_RUN=1
  echo "=== DRY RUN — no actions will be taken ==="
fi

# ---------------------------------------------------------------------------
# Preflight
# ---------------------------------------------------------------------------
preflight() {
  echo "== Preflight =="
  command -v python3 >/dev/null || { echo "FATAL: python3 not found"; exit 1; }
  command -v git >/dev/null     || { echo "FATAL: git not found"; exit 1; }
  command -v tmux >/dev/null    || { echo "FATAL: tmux not found"; exit 1; }

  # Ensure main is clean
  cd "$LOOM_ROOT"
  if [[ -n "$(git status --porcelain)" ]]; then
    echo "WARNING: main worktree has uncommitted changes"
    echo "Proceeding anyway — worktrees branch from HEAD"
  fi

  # Verify data files exist
  [[ -f "$LOOM_ROOT/web/data/crispr-targets.json" ]] || { echo "FATAL: crispr-targets.json missing"; exit 1; }
  [[ -f "$LOOM_ROOT/web/data/ontology-enrichment.json" ]] || { echo "FATAL: ontology-enrichment.json missing"; exit 1; }

  echo "Data files OK"
  echo "Main branch: $(git rev-parse --short HEAD)"
}

# ---------------------------------------------------------------------------
# Worktree setup
# ---------------------------------------------------------------------------
setup_worktrees() {
  echo "== Setting up worktrees =="
  cd "$LOOM_ROOT"

  for pair in "$BRANCH_T1:$WORKTREE_T1" "$BRANCH_T2:$WORKTREE_T2" "$BRANCH_T3:$WORKTREE_T3"; do
    branch="${pair%%:*}"
    wt="${pair##*:}"
    if [[ -d "$wt" ]]; then
      echo "  Worktree $wt already exists — removing"
      git worktree remove --force "$wt" 2>/dev/null || true
    fi
    # Delete branch if it exists (stale from previous run)
    git branch -D "$branch" 2>/dev/null || true
    echo "  Creating $wt on branch $branch"
    git worktree add -b "$branch" "$wt" HEAD
  done
}

# ---------------------------------------------------------------------------
# Tier 1: Gene-annotate targets (positional join)
# ---------------------------------------------------------------------------
tier1_script() {
  cat << 'TIER1_EOF'
#!/usr/bin/env python3
"""
Tier 1: Annotate CRISPR targets with gene information.

For each target, binary-search its position against gene [start,end] intervals
to find which gene (if any) it falls in. Adds gene_symbol, gene_name, gene_type
to each target. Writes annotated crispr-targets.json.

Coordinate alignment note:
- Target positions (p) are in concatenated corpus space (multiple genomes).
- Gene positions are from a single reference genome.
- For pathogens where many genomes are concatenated, positions may exceed
  the reference genome length. We handle this via modular mapping: try
  p directly, then p % genome_length, annotating what we can.
"""
import json, bisect, sys, os, time

def build_gene_index(genes_data):
    """Build sorted interval list for binary search."""
    genes = genes_data.get("genes", [])
    if not genes:
        return [], 0
    # Sort by start position
    intervals = []
    for g in genes:
        start = g.get("start")
        end = g.get("end")
        if start is not None and end is not None and start > 0 and end > 0:
            intervals.append((min(start, end), max(start, end), g))
    intervals.sort(key=lambda x: x[0])
    
    # Estimate genome length from max gene end
    max_end = max(iv[1] for iv in intervals) if intervals else 0
    return intervals, max_end

def find_gene(intervals, starts, pos):
    """Binary search for gene containing position."""
    idx = bisect.bisect_right(starts, pos) - 1
    if idx < 0:
        return None
    start, end, gene = intervals[idx]
    if start <= pos <= end:
        return gene
    # Check next interval too (overlapping genes possible)
    if idx + 1 < len(intervals):
        start2, end2, gene2 = intervals[idx + 1]
        if start2 <= pos <= end2:
            return gene2
    return None

def annotate_pathogen(pathogen_key, targets_data, enrichment_data):
    """Annotate all targets for one pathogen."""
    enrich = enrichment_data.get("pathogens", {}).get(pathogen_key)
    if not enrich or "genes" not in enrich:
        return targets_data, {"annotated": 0, "total": len(targets_data.get("targets", [])), "reason": "no enrichment"}
    
    genes_data = enrich["genes"]
    intervals, genome_len = build_gene_index(genes_data)
    if not intervals:
        return targets_data, {"annotated": 0, "total": len(targets_data.get("targets", [])), "reason": "no valid gene intervals"}
    
    starts = [iv[0] for iv in intervals]
    annotated = 0
    total = 0
    
    for target in targets_data.get("targets", []):
        total += 1
        pos = target.get("p", 0)
        
        # Try direct position first
        gene = find_gene(intervals, starts, pos)
        
        # If no hit and position exceeds genome, try modular mapping
        if gene is None and genome_len > 0 and pos > genome_len:
            mod_pos = pos % genome_len
            gene = find_gene(intervals, starts, mod_pos)
        
        if gene:
            target["g"] = gene.get("symbol") or gene.get("name", "")
            target["gn"] = gene.get("name", "")
            target["gt"] = gene.get("type", "")
            annotated += 1
    
    stats = {"annotated": annotated, "total": total, "gene_intervals": len(intervals), "genome_len": genome_len}
    return targets_data, stats

def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    targets_path = os.path.join(root, "web", "data", "crispr-targets.json")
    enrichment_path = os.path.join(root, "web", "data", "ontology-enrichment.json")
    output_path = os.path.join(root, "web", "data", "crispr-targets.json")
    
    print("=== Tier 1: Gene-Annotate Targets ===")
    print(f"Targets: {targets_path}")
    print(f"Enrichment: {enrichment_path}")
    
    with open(targets_path) as f:
        targets = json.load(f)
    with open(enrichment_path) as f:
        enrichment = json.load(f)
    
    all_stats = {}
    total_annotated = 0
    total_targets = 0
    
    for pathogen_key, pathogen_data in targets.items():
        t0 = time.time()
        targets[pathogen_key], stats = annotate_pathogen(pathogen_key, pathogen_data, enrichment)
        elapsed = time.time() - t0
        all_stats[pathogen_key] = stats
        total_annotated += stats.get("annotated", 0)
        total_targets += stats.get("total", 0)
        rate = stats.get("annotated", 0)
        total = stats.get("total", 0)
        pct = (rate / total * 100) if total > 0 else 0
        print(f"  {pathogen_key:25s} {rate:5d}/{total:5d} annotated ({pct:5.1f}%) in {elapsed:.2f}s")
    
    print(f"\nTotal: {total_annotated}/{total_targets} targets annotated "
          f"({total_annotated/total_targets*100:.1f}%)" if total_targets else "")
    
    # Write annotated targets
    with open(output_path, "w") as f:
        json.dump(targets, f, separators=(",", ":"))
    
    print(f"Written: {output_path} ({os.path.getsize(output_path) / 1024:.0f} KB)")
    
    # Write stats report
    stats_path = os.path.join(root, "web", "data", "gene-annotation-stats.json")
    with open(stats_path, "w") as f:
        json.dump({"generated": time.strftime("%Y-%m-%dT%H:%M:%SZ"), "stats": all_stats,
                    "totals": {"annotated": total_annotated, "targets": total_targets}}, f, indent=2)
    print(f"Stats: {stats_path}")

if __name__ == "__main__":
    main()
TIER1_EOF
}

# ---------------------------------------------------------------------------
# Tier 2: Diagnostic priority scoring
# ---------------------------------------------------------------------------
tier2_script() {
  cat << 'TIER2_EOF'
#!/usr/bin/env python3
"""
Tier 2: Compute diagnostic priority scores for each pathogen.

Scoring factors:
- CRISPR diagnostic gap (no published dx = highest novelty)
- Disease burden (CFR, annual cases/deaths)
- WHO classification (epidemic-prone, pandemic potential)
- Existing rapid test availability (fewer = more need)
- Gene coverage (targets falling in known genes = more designable)

Output: diagnostic-priority.json consumed by the web UI.
"""
import json, os, re, time

def parse_number_range(s):
    """Parse strings like '1.3-4M' or '21,000-143,000' into a midpoint estimate."""
    if not s:
        return 0
    s = s.strip().lower().replace(",", "")
    # Handle M/K suffixes
    multiplier = 1
    if s.endswith("m"):
        s = s[:-1]
        multiplier = 1_000_000
    elif s.endswith("k"):
        s = s[:-1]
        multiplier = 1_000
    # Try range
    m = re.match(r"([\d.]+)\s*[-–]\s*([\d.]+)", s)
    if m:
        low, high = float(m.group(1)) * multiplier, float(m.group(2)) * multiplier
        return (low + high) / 2
    # Try single number
    m = re.match(r"([\d.]+)", s)
    if m:
        return float(m.group(1)) * multiplier
    return 0

def parse_cfr(s):
    """Parse CFR string into a numeric estimate (0-1 scale)."""
    if not s:
        return 0
    s = s.lower()
    # Look for percentages
    pcts = re.findall(r"([\d.]+)\s*%", s)
    if pcts:
        # Take the highest mentioned CFR (worst case)
        return max(float(p) / 100 for p in pcts)
    return 0

def score_pathogen(key, enrichment_data, annotation_stats):
    """Compute composite diagnostic priority score (0-100)."""
    enrich = enrichment_data.get("pathogens", {}).get(key, {})
    who = enrich.get("who_context", {})
    genes = enrich.get("genes", {})
    
    # Factor 1: CRISPR diagnostic gap (0 or 30 points)
    dx_status = who.get("crispr_dx_status", "")
    if not dx_status or "none" in dx_status.lower():
        dx_gap_score = 30
    elif "limited" in dx_status.lower() or "research" in dx_status.lower():
        dx_gap_score = 15
    else:
        dx_gap_score = 0
    
    # Factor 2: Disease burden — annual deaths (0-25 points)
    deaths_str = who.get("annual_deaths", "")
    deaths_mid = parse_number_range(deaths_str)
    if deaths_mid > 100_000:
        burden_score = 25
    elif deaths_mid > 10_000:
        burden_score = 20
    elif deaths_mid > 1_000:
        burden_score = 15
    elif deaths_mid > 100:
        burden_score = 10
    else:
        burden_score = 5
    
    # Factor 3: CFR severity (0-15 points)
    cfr = parse_cfr(who.get("case_fatality_rate", ""))
    if cfr >= 0.25:
        cfr_score = 15
    elif cfr >= 0.05:
        cfr_score = 10
    elif cfr >= 0.01:
        cfr_score = 5
    else:
        cfr_score = 2
    
    # Factor 4: WHO classification (0-15 points)
    classification = who.get("classification", "").lower()
    if "pandemic" in classification:
        who_score = 15
    elif "epidemic" in classification:
        who_score = 12
    elif "endemic" in classification:
        who_score = 8
    else:
        who_score = 5
    
    # Factor 5: Existing rapid tests (0-15 points, fewer = higher need)
    rapid = who.get("existing_rapid_tests", "")
    rapid_str = ", ".join(rapid) if isinstance(rapid, list) else (rapid or "")
    if not rapid_str or "none" in rapid_str.lower():
        rapid_score = 15
    elif "limited" in rapid_str.lower():
        rapid_score = 10
    else:
        rapid_score = 3
    
    total = dx_gap_score + burden_score + cfr_score + who_score + rapid_score
    
    # Gene annotation rate from Tier 1
    stats = annotation_stats.get("stats", {}).get(key, {})
    gene_annotation_rate = 0
    if stats.get("total", 0) > 0:
        gene_annotation_rate = stats.get("annotated", 0) / stats["total"]
    
    return {
        "pathogen": key,
        "species": enrich.get("species", key),
        "total_score": total,
        "max_score": 100,
        "factors": {
            "crispr_dx_gap": {"score": dx_gap_score, "max": 30, "detail": dx_status or "unknown"},
            "disease_burden": {"score": burden_score, "max": 25, "detail": deaths_str or "unknown"},
            "cfr_severity": {"score": cfr_score, "max": 15, "detail": who.get("case_fatality_rate", "unknown")},
            "who_classification": {"score": who_score, "max": 15, "detail": classification or "unknown"},
            "rapid_test_gap": {"score": rapid_score, "max": 15, "detail": rapid_str or "unknown"},
        },
        "gene_annotation_rate": round(gene_annotation_rate, 3),
        "gene_count": genes.get("gene_count", 0),
        "target_gene_coverage": stats.get("annotated", 0),
    }

def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    enrichment_path = os.path.join(root, "web", "data", "ontology-enrichment.json")
    stats_path = os.path.join(root, "web", "data", "gene-annotation-stats.json")
    output_path = os.path.join(root, "web", "data", "diagnostic-priority.json")
    
    print("=== Tier 2: Diagnostic Priority Scoring ===")
    
    with open(enrichment_path) as f:
        enrichment = json.load(f)
    
    # Load Tier 1 stats if available
    annotation_stats = {}
    if os.path.exists(stats_path):
        with open(stats_path) as f:
            annotation_stats = json.load(f)
        print(f"Loaded gene annotation stats from Tier 1")
    else:
        print("WARNING: No Tier 1 stats found — gene_annotation_rate will be 0")
    
    results = []
    for key in enrichment.get("pathogens", {}):
        score = score_pathogen(key, enrichment, annotation_stats)
        results.append(score)
        print(f"  {key:25s} score={score['total_score']:3d}/100  "
              f"dx_gap={score['factors']['crispr_dx_gap']['score']:2d}  "
              f"burden={score['factors']['disease_burden']['score']:2d}  "
              f"cfr={score['factors']['cfr_severity']['score']:2d}")
    
    # Sort by score descending
    results.sort(key=lambda x: x["total_score"], reverse=True)
    
    print(f"\n=== Top Opportunities ===")
    for i, r in enumerate(results[:5]):
        print(f"  {i+1}. {r['species']:30s} {r['total_score']}/100  "
              f"({r['factors']['crispr_dx_gap']['detail']})")
    
    output = {
        "generated": time.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "scoring_version": "1.0",
        "pathogens": {r["pathogen"]: r for r in results},
        "ranked": [r["pathogen"] for r in results],
    }
    
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nWritten: {output_path}")

if __name__ == "__main__":
    main()
TIER2_EOF
}

# ---------------------------------------------------------------------------
# Tier 3: UI integration
# ---------------------------------------------------------------------------
tier3_instructions() {
  cat << 'TIER3_EOF'
# Tier 3: UI Gene Integration — Instructions for Copilot

## Context
Tier 1 has annotated crispr-targets.json with gene fields: g (symbol), gn (name), gt (type).
Tier 2 has produced diagnostic-priority.json with scores for each pathogen.

## Tasks

### 3a. Target table gene column
In web/crispr-app.js, update the Targets Database tab:
- Add "Gene" column to the targets table header
- Render target.g (gene symbol) in each row; if empty show "intergenic"
- Add gene_name as tooltip on hover (target.gn)
- Color protein-coding targets differently from rRNA/tRNA/intergenic

### 3b. Gene filter dropdown
- Add a dropdown above the targets table: "Filter by gene"
- Populate from unique gene symbols in current pathogen's targets
- Options: "All genes", "Protein-coding only", then each gene symbol
- Selecting a gene filters the table instantly

### 3c. Diagnostic priority badges on pathogen cards
- Load diagnostic-priority.json in loadData()
- On each pathogen card, show a small priority score badge (e.g., "87/100")
- Color code: >=70 red (high need), 40-69 amber, <40 green
- Tooltip: "Diagnostic novelty score — higher = greater unmet need"

### 3d. Top Opportunities panel
- Add a "Top Opportunities" section above the pathogen grid (or as a banner)
- Show the top 3 pathogens by diagnostic priority score
- Each entry: pathogen name, score, key factors (e.g., "No CRISPR dx · 143K deaths/yr")
- Clicking navigates to that pathogen's targets

### 3e. Gene Map ↔ Targets bidirectional link
In Disease Context tab:
- Make gene rows in the Gene Map clickable
- Clicking a gene switches to Targets Database tab pre-filtered to that gene
In Targets Database tab:
- Clicking a gene name in the Gene column switches to Disease Context + scrolls to Gene Map

### Commit
After all changes, commit: "feat: integrate gene annotations + diagnostic priority into UI"
TIER3_EOF
}

# ---------------------------------------------------------------------------
# Main launch
# ---------------------------------------------------------------------------
main() {
  preflight

  echo ""
  echo "== Plan =="
  echo "  Tier 1: Gene-annotate targets (Python, ~30s)"
  echo "    Branch: $BRANCH_T1"
  echo "    Worktree: $WORKTREE_T1"
  echo ""
  echo "  Tier 2: Diagnostic priority scoring (Python, ~5s)"
  echo "    Branch: $BRANCH_T2"
  echo "    Worktree: $WORKTREE_T2"
  echo "    Depends: Tier 1 stats"
  echo ""
  echo "  Tier 3: UI integration (JS/CSS, Copilot-driven)"
  echo "    Branch: $BRANCH_T3"
  echo "    Worktree: $WORKTREE_T3"
  echo "    Depends: Tier 1 + Tier 2 outputs"
  echo ""

  if (( DRY_RUN )); then
    echo "=== DRY RUN complete ==="
    exit 0
  fi

  # Setup
  mkdir -p "$LOG_DIR"
  setup_worktrees

  # Copy data files to each worktree (they need the source data)
  for wt in "$WORKTREE_T1" "$WORKTREE_T2" "$WORKTREE_T3"; do
    mkdir -p "$wt/web/data"
    cp "$LOOM_ROOT/web/data/crispr-targets.json" "$wt/web/data/"
    cp "$LOOM_ROOT/web/data/ontology-enrichment.json" "$wt/web/data/"
    # Also copy pubmed scan if it exists
    [[ -f "$LOOM_ROOT/web/data/pubmed-scan-results.json" ]] && \
      cp "$LOOM_ROOT/web/data/pubmed-scan-results.json" "$wt/web/data/" || true
  done

  # Write Tier 1 script
  tier1_script > "$WORKTREE_T1/scripts/annotate_targets_with_genes.py"
  chmod +x "$WORKTREE_T1/scripts/annotate_targets_with_genes.py"

  # Write Tier 2 script
  tier2_script > "$WORKTREE_T2/scripts/compute_diagnostic_priority.py"
  chmod +x "$WORKTREE_T2/scripts/compute_diagnostic_priority.py"

  # Write Tier 3 instructions
  tier3_instructions > "$WORKTREE_T3/TIER3_INSTRUCTIONS.md"

  # ---------------------------------------------------------------------------
  # tmux session
  # ---------------------------------------------------------------------------
  tmux kill-session -t "$SESSION" 2>/dev/null || true
  tmux new-session -d -s "$SESSION" -n "monitor"

  # Window 0: Monitor
  tmux send-keys -t "$SESSION:monitor" "echo '=== Ontology Integration Pipeline ===' && echo 'Started: $(date)' && echo '' && echo 'Worktrees:' && cd $LOOM_ROOT && git worktree list | grep -E 'gene-annotate|dx-priority|ui-gene' && echo '' && echo 'Logs: $LOG_DIR' && echo '' && echo 'Watch progress: tail -f $LOG_DIR/*.log'" Enter

  # Window 1: Tier 1 — Gene annotation
  tmux new-window -t "$SESSION" -n "tier1-genes"
  tmux send-keys -t "$SESSION:tier1-genes" "cd $WORKTREE_T1 && source $LOOM_ROOT/.venv/bin/activate 2>/dev/null; echo '=== Tier 1: Gene Annotation ===' && python3 scripts/annotate_targets_with_genes.py 2>&1 | tee $LOG_DIR/tier1-genes.log && echo '' && echo 'Tier 1 DONE — committing...' && git add -f web/data/crispr-targets.json web/data/gene-annotation-stats.json scripts/annotate_targets_with_genes.py && git commit -m 'feat: annotate CRISPR targets with gene positions

- Positional join: map each target to its gene via [start,end] overlap
- Adds g (symbol), gn (name), gt (type) fields to each target
- Handles multi-genome coordinate spaces via modular mapping
- Generates gene-annotation-stats.json with coverage metrics' && echo 'Tier 1 COMMITTED' && echo '' && echo 'Copying outputs to Tier 2 and Tier 3 worktrees...' && cp web/data/gene-annotation-stats.json $WORKTREE_T2/web/data/ && cp web/data/crispr-targets.json $WORKTREE_T2/web/data/ && cp web/data/gene-annotation-stats.json $WORKTREE_T3/web/data/ && cp web/data/crispr-targets.json $WORKTREE_T3/web/data/ && echo 'Outputs copied. Tier 2 can proceed.'" Enter

  # Window 2: Tier 2 — Diagnostic priority (waits for Tier 1 stats)
  tmux new-window -t "$SESSION" -n "tier2-priority"
  tmux send-keys -t "$SESSION:tier2-priority" "cd $WORKTREE_T2 && source $LOOM_ROOT/.venv/bin/activate 2>/dev/null; echo '=== Tier 2: Diagnostic Priority ===' && echo 'Waiting for Tier 1 stats...' && while [[ ! -f web/data/gene-annotation-stats.json ]]; do sleep 2; done && echo 'Tier 1 stats received. Scoring...' && python3 scripts/compute_diagnostic_priority.py 2>&1 | tee $LOG_DIR/tier2-priority.log && echo '' && echo 'Tier 2 DONE — committing...' && git add -f web/data/diagnostic-priority.json scripts/compute_diagnostic_priority.py && git commit -m 'feat: compute diagnostic priority scores for pathogens

- Score 0-100 based on: CRISPR dx gap, disease burden, CFR, WHO class, rapid test gap
- Integrates Tier 1 gene annotation coverage metrics
- Produces ranked list of diagnostic opportunities
- Top pathogens by unmet need surfaced for the UI' && echo 'Tier 2 COMMITTED' && cp web/data/diagnostic-priority.json $WORKTREE_T3/web/data/ && echo 'Priority data copied to Tier 3.'" Enter

  # Window 3: Tier 3 — UI integration (waits for both Tier 1 and 2)
  tmux new-window -t "$SESSION" -n "tier3-ui"
  tmux send-keys -t "$SESSION:tier3-ui" "cd $WORKTREE_T3 && echo '=== Tier 3: UI Integration ===' && echo 'Waiting for Tier 1 + Tier 2 outputs...' && while [[ ! -f web/data/diagnostic-priority.json ]] || [[ ! -f web/data/gene-annotation-stats.json ]]; do sleep 3; done && echo 'All data ready. Tier 3 instructions in TIER3_INSTRUCTIONS.md' && echo '' && cat TIER3_INSTRUCTIONS.md && echo '' && echo '=== Tier 3 is ready for Copilot to implement ===' && echo 'Run: cd $WORKTREE_T3 && code .' 2>&1 | tee $LOG_DIR/tier3-ui.log" Enter

  # Window 4: Logs
  tmux new-window -t "$SESSION" -n "logs"
  tmux send-keys -t "$SESSION:logs" "echo 'Waiting for logs...' && sleep 5 && tail -f $LOG_DIR/*.log 2>/dev/null || echo 'No logs yet'" Enter

  echo ""
  echo "==========================================="
  echo " Pipeline launched in tmux session: $SESSION"
  echo "==========================================="
  echo ""
  echo "  tmux attach -t $SESSION         # monitor"
  echo "  tmux select-window -t $SESSION:1  # Tier 1"
  echo "  tmux select-window -t $SESSION:2  # Tier 2"
  echo "  tmux select-window -t $SESSION:3  # Tier 3"
  echo ""
  echo "  Tier 1+2 run automatically."
  echo "  Tier 3 prints instructions for Copilot to implement."
  echo ""
  echo "After all tiers complete, merge to main:"
  echo "  cd $LOOM_ROOT"
  echo "  git merge $BRANCH_T1 $BRANCH_T2 $BRANCH_T3"
  echo "  git worktree remove $WORKTREE_T1"
  echo "  git worktree remove $WORKTREE_T2"
  echo "  git worktree remove $WORKTREE_T3"
  echo ""
}

main
