#!/usr/bin/env python3
"""Compute guide RNA quality scores for all pathogen CRISPR targets.

Scores are based on empirically validated sequence features:
- GC content of spacer (optimal 40-70%)
- No poly-T runs >=4 in spacer (terminates Pol III transcription)
- No extreme GC at seed region (last 12nt of spacer)
- No homopolymer runs >=5 in spacer
- Self-complementarity check (hairpin potential)

Output: web/data/guide-scores.json
"""
import csv
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

CRISPR_DB = Path("/Volumes/T9/loom-openscience/crispr-db")
OUT = Path("web/data/guide-scores.json")

PATHOGENS = [
    "cholera", "dengue", "ebola", "hepatitis-b", "hiv-1",
    "influenza-a", "mers", "mpox", "rsv", "tuberculosis", "zika",
]


def gc_content(seq: str) -> float:
    """GC fraction of sequence."""
    gc = sum(1 for c in seq.upper() if c in "GC")
    return gc / len(seq) if seq else 0.0


def longest_homopolymer(seq: str) -> int:
    """Length of longest single-nucleotide run."""
    if not seq:
        return 0
    max_run = 1
    cur_run = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur_run += 1
            max_run = max(max_run, cur_run)
        else:
            cur_run = 1
    return max_run


def has_poly_t(seq: str, min_run: int = 4) -> bool:
    """Check for poly-T run (terminates Pol III)."""
    return "T" * min_run in seq.upper()


def self_complementarity_score(seq: str) -> float:
    """Simple self-complementarity check: fraction of self-pairing positions."""
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    seq = seq.upper()
    rev_comp = "".join(comp.get(c, "N") for c in reversed(seq))
    # Check for internal hairpin: slide reverse complement
    best = 0
    n = len(seq)
    for offset in range(4, n - 4):
        matches = sum(1 for i in range(min(n - offset, n)) if i + offset < n and seq[i] == rev_comp[i + offset])
        best = max(best, matches)
    return best / n if n else 0.0


def score_guide(seq_23mer: str, pam_type: str) -> dict:
    """Score a 23-mer guide (20bp spacer + 3bp PAM context).

    Returns dict with individual scores and composite 0-100 score.
    """
    # Spacer is first 20nt for NGG, last 20nt for TTTN
    if pam_type == "NGG":
        spacer = seq_23mer[:20]
    elif pam_type == "TTTN":
        spacer = seq_23mer[4:24] if len(seq_23mer) >= 24 else seq_23mer[3:]
    else:
        spacer = seq_23mer[:20]

    gc = gc_content(spacer)
    seed_gc = gc_content(spacer[-12:]) if len(spacer) >= 12 else gc
    poly_t = has_poly_t(spacer)
    max_homo = longest_homopolymer(spacer)
    self_comp = self_complementarity_score(spacer)

    # Component scores (each 0-1, higher = better)

    # GC content: optimal 40-70%, penalize outside
    if 0.40 <= gc <= 0.70:
        gc_score = 1.0
    elif 0.30 <= gc <= 0.80:
        gc_score = 0.6
    elif 0.20 <= gc <= 0.90:
        gc_score = 0.3
    else:
        gc_score = 0.1

    # Seed region GC: optimal 30-70%
    if 0.30 <= seed_gc <= 0.70:
        seed_score = 1.0
    elif 0.20 <= seed_gc <= 0.80:
        seed_score = 0.5
    else:
        seed_score = 0.2

    # Poly-T penalty
    poly_t_score = 0.2 if poly_t else 1.0

    # Homopolymer penalty
    if max_homo <= 3:
        homo_score = 1.0
    elif max_homo == 4:
        homo_score = 0.5
    else:
        homo_score = 0.1

    # Self-complementarity (lower is better)
    comp_score = max(0.0, 1.0 - self_comp * 3)

    # PAM bonus: NGG (SpCas9) is gold standard
    pam_score = 1.0 if pam_type == "NGG" else 0.7

    # Composite: weighted average
    composite = (
        gc_score * 0.25 +
        seed_score * 0.15 +
        poly_t_score * 0.20 +
        homo_score * 0.15 +
        comp_score * 0.10 +
        pam_score * 0.15
    )

    return {
        "gc": round(gc * 100, 1),
        "seed_gc": round(seed_gc * 100, 1),
        "poly_t": poly_t,
        "max_homopolymer": max_homo,
        "self_comp": round(self_comp, 3),
        "score": round(composite * 100),
    }


def main():
    result = {
        "generated": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "scoring_version": "1.0",
        "method": "Composite score: GC content (25%), poly-T penalty (20%), seed GC (15%), homopolymer (15%), PAM type (15%), self-complementarity (10%)",
        "pathogens": {},
    }

    total_scored = 0
    for pathogen in PATHOGENS:
        csv_path = CRISPR_DB / f"{pathogen}_targets.csv"
        if not csv_path.exists():
            print(f"  WARN: {csv_path} not found, skipping")
            continue

        scores = {}
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                seq = row["sequence_23mer"]
                pam = row.get("pam_type", "none")
                sc = score_guide(seq, pam)
                scores[seq] = sc

        result["pathogens"][pathogen] = scores
        total_scored += len(scores)

        # Stats
        vals = [s["score"] for s in scores.values()]
        avg = sum(vals) / len(vals) if vals else 0
        high = sum(1 for v in vals if v >= 70)
        print(f"  {pathogen}: {len(scores)} guides scored, avg={avg:.0f}, high-quality={high}")

    # Summary
    all_scores = [s["score"] for p in result["pathogens"].values() for s in p.values()]
    result["summary"] = {
        "total_guides_scored": total_scored,
        "mean_score": round(sum(all_scores) / len(all_scores)) if all_scores else 0,
        "high_quality_count": sum(1 for v in all_scores if v >= 70),
        "high_quality_pct": round(sum(1 for v in all_scores if v >= 70) / len(all_scores) * 100, 1) if all_scores else 0,
    }

    with open(OUT, "w") as f:
        json.dump(result, f)  # compact — this file is large

    print(f"\n  Output: {OUT} ({OUT.stat().st_size / 1024:.0f} KB)")
    print(f"  Total: {total_scored} guides, mean={result['summary']['mean_score']}, "
          f"high-quality={result['summary']['high_quality_pct']}%")


if __name__ == "__main__":
    main()
