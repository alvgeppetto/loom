#!/usr/bin/env python3
"""Standalone CRISPR target scanner for a single pathogen.

Usage:
    python scripts/scan_crispr_targets.py <name> <idx_path> <out_csv>

Scans an FM-index for all 23-mer CRISPR targets with PAM annotation.
Uses Aho-Corasick for efficient multi-pattern matching.
"""
import csv
import json
import os
import re
import struct
import sys
import time
from pathlib import Path

import ahocorasick

def load_corpus(idx_path: str) -> str:
    """Load corpus text from a .idx file."""
    payload = Path(idx_path).read_bytes()
    cursor = 0
    metadata_len = struct.unpack("<Q", payload[cursor:cursor + 8])[0]
    cursor += 8 + metadata_len
    corpus_len = struct.unpack("<Q", payload[cursor:cursor + 8])[0]
    cursor += 8
    corpus = payload[cursor:cursor + corpus_len].decode("utf-8", errors="ignore")
    return corpus


def extract_23mers_with_pam(corpus: str) -> list[dict]:
    """Sliding window extraction of all unique 23-mers with PAM classification."""
    # Step 1: Extract all unique 23-mers via sliding window
    t0 = time.time()
    target_counts = {}
    first_pos = {}
    clean = corpus.upper()
    dna_chars = set("ACGT")
    
    i = 0
    n = len(clean)
    while i <= n - 23:
        # Find next valid 23-mer (all ACGT)
        seq = clean[i:i+23]
        if all(c in dna_chars for c in seq):
            if seq in target_counts:
                target_counts[seq] += 1
            else:
                target_counts[seq] = 1
                first_pos[seq] = i
            i += 1
        else:
            # Skip to after the invalid character
            bad = next((j for j, c in enumerate(seq) if c not in dna_chars), 0)
            i += bad + 1
    
    elapsed = time.time() - t0
    print(f"  Extracted {len(target_counts)} unique 23-mers in {elapsed:.1f}s")
    
    # Step 2: Classify PAM
    targets = []
    for seq, count in target_counts.items():
        pam = "none"
        # SpCas9: last 3 bases = XGG
        if seq[-2:] == "GG":
            pam = "NGG"
        # Cas12a: first 4 bases = TTTX
        elif seq[:3] == "TTT":
            pam = "TTTN"
        
        targets.append({
            "sequence_23mer": seq,
            "position": first_pos[seq],
            "pam_type": pam,
            "occurrences": count,
        })
    
    # Sort by occurrences descending
    targets.sort(key=lambda t: t["occurrences"], reverse=True)
    
    # Cap at top 10000 targets to keep data manageable
    if len(targets) > 10000:
        print(f"  Capping from {len(targets)} to 10000 targets")
        targets = targets[:10000]
    
    return targets


def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <name> <idx_path> <out_csv>")
        sys.exit(1)
    
    name = sys.argv[1]
    idx_path = sys.argv[2]
    out_csv = sys.argv[3]
    
    print(f"Scanning {name} from {idx_path}...")
    t0 = time.time()
    
    corpus = load_corpus(idx_path)
    print(f"  Corpus: {len(corpus):,} characters")
    
    targets = extract_23mers_with_pam(corpus)
    
    # Write CSV
    with open(out_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            "pathogen", "gene", "position", "sequence_23mer",
            "pam_type", "source_motif", "occurrences"
        ])
        writer.writeheader()
        for t in targets:
            writer.writerow({
                "pathogen": name,
                "gene": f"{name}.fna",
                "position": t["position"],
                "sequence_23mer": t["sequence_23mer"],
                "pam_type": t["pam_type"],
                "source_motif": "sliding_window",
                "occurrences": t["occurrences"],
            })
    
    elapsed = time.time() - t0
    print(f"  Done: {len(targets)} targets in {elapsed:.1f}s → {out_csv}")


if __name__ == "__main__":
    main()
