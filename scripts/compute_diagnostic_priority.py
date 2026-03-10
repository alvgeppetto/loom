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
    rapid = who.get("existing_rapid_tests", [])
    if isinstance(rapid, list):
        rapid_str = ", ".join(rapid) if rapid else ""
    else:
        rapid_str = str(rapid)
    if not rapid_str or "none" in rapid_str.lower():
        rapid_score = 15
    elif "limited" in rapid_str.lower() or "pcr" in rapid_str.lower():
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
