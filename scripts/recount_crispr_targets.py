#!/usr/bin/env python3
"""Recount CRISPR target occurrences using actual FM-index corpus data.

The original scan script incorrectly assigned the *motif's* occurrence count
to each 23-mer extracted from context windows.  This script loads each
pathogen's FM-index, counts every 23-mer properly, computes real genomic
positions, filters junk entries, and rebuilds crispr-targets.json.

Usage:
    python scripts/recount_crispr_targets.py
"""
from __future__ import annotations

import csv
import json
import os
import struct
import sys
import time
from collections import OrderedDict
from pathlib import Path

import ahocorasick

ROOT = Path(__file__).resolve().parents[1]
IDX_DIR = Path("/Volumes/T9/loom-openscience/indexes")
CSV_DIR = Path("/Volumes/T9/loom-openscience/crispr-db")
OUT_JSON = ROOT / "web" / "data" / "crispr-targets.json"

META = OrderedDict([
    ("dengue",       {"name": "Dengue Virus",      "index_mb": 248, "genomes": 30158,  "family": "Flaviviridae"}),
    ("ebola",        {"name": "Ebola Virus",        "index_mb": 62,  "genomes": 24055,  "family": "Filoviridae"}),
    ("hepatitis-b",  {"name": "Hepatitis B",        "index_mb": 130, "genomes": 47636,  "family": "Hepadnaviridae"}),
    ("mers",         {"name": "MERS-CoV",           "index_mb": 26,  "genomes": 8986,   "family": "Coronaviridae"}),
    ("refseq-viral", {"name": "RefSeq Viral (all)", "index_mb": 553, "genomes": 288215, "family": "Mixed"}),
    ("rsv",          {"name": "RSV",                "index_mb": 368, "genomes": 376998, "family": "Pneumoviridae"}),
    ("zika",         {"name": "Zika Virus",         "index_mb": 19,  "genomes": 1689,   "family": "Flaviviridae"}),
])


def load_corpus(idx_path: Path) -> tuple[str, list[tuple[str, int, int]]]:
    """Load corpus text and file boundaries from a .idx file."""
    payload = idx_path.read_bytes()
    cursor = 0
    metadata_len = struct.unpack("<Q", payload[cursor:cursor + 8])[0]
    cursor += 8
    metadata = json.loads(payload[cursor:cursor + metadata_len].decode("utf-8"))
    cursor += metadata_len

    corpus_len = struct.unpack("<Q", payload[cursor:cursor + 8])[0]
    cursor += 8
    corpus = payload[cursor:cursor + corpus_len].decode("utf-8", errors="ignore")

    # Build file boundaries
    files = []
    for item in metadata.get("files", []):
        path = item["path"]
        start = int(item["start_offset"])
        end = int(item["end_offset"])
        marker = "\x01" + path + "\x01"
        content_start = start + len(marker)
        if content_start < end <= len(corpus):
            files.append((path, content_start, end))

    return corpus, files


def build_fasta_index(corpus: str) -> list[int]:
    """Build sorted list of FASTA record start positions (position of '>' char).

    Used to map absolute corpus positions back to within-genome coordinates.
    """
    headers = []
    pos = 0
    while True:
        pos = corpus.find(">", pos)
        if pos < 0:
            break
        headers.append(pos)
        pos += 1
    return headers


def genome_relative_position(corpus: str, fasta_starts: list[int], abs_pos: int) -> int:
    """Map absolute corpus position to base-pair position within its FASTA record.

    Counts only DNA bases (A/C/G/T/N) between the end of the header line and
    the target position, excluding newlines and other formatting.
    """
    import bisect
    # Find the FASTA record containing this position
    idx = bisect.bisect_right(fasta_starts, abs_pos) - 1
    if idx < 0:
        return abs_pos  # before any header (shouldn't happen)

    record_start = fasta_starts[idx]
    # Skip past the header line (find first \n after >)
    seq_start = corpus.find("\n", record_start)
    if seq_start < 0 or seq_start >= abs_pos:
        return 0
    seq_start += 1  # skip the newline

    # Count DNA bases from seq_start to abs_pos
    region = corpus[seq_start:abs_pos]
    bp = sum(1 for c in region if c in "ACGTNacgtn")
    return bp


def process_pathogen(key: str, meta: dict) -> dict:
    """Recount all targets for one pathogen using a single-pass sliding window."""
    csv_path = CSV_DIR / f"{key}_targets.csv"
    idx_path = IDX_DIR / f"{key}.idx"

    if not csv_path.exists():
        print(f"  SKIP {key}: no CSV")
        return None
    if not idx_path.exists():
        print(f"  SKIP {key}: no index")
        return None

    # Read original targets (23-mers only, skip markers)
    target_set = set()
    pam_map = {}  # seq -> pam_type
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            seq = row["sequence_23mer"]
            if not seq or "[" in seq or not all(c in "ACGT" for c in seq):
                continue
            target_set.add(seq)
            pam = row["pam_type"]
            if pam == "NGG" or (seq not in pam_map):
                pam_map[seq] = pam

    print(f"  {key}: {len(target_set)} unique 23-mers to count")

    # Load corpus
    t0 = time.time()
    print(f"  Loading {idx_path.name} ({meta['index_mb']} MB)...", end=" ", flush=True)
    corpus, files = load_corpus(idx_path)
    corpus_len = len(corpus)
    print(f"loaded in {time.time() - t0:.1f}s ({corpus_len:,} chars)")

    # Single-pass Aho-Corasick: count all target 23-mers and find first positions
    t0 = time.time()
    A = ahocorasick.Automaton()
    for seq in target_set:
        A.add_word(seq, seq)
    A.make_automaton()

    counts = {}      # seq -> count
    first_pos = {}   # seq -> absolute position of first occurrence

    for end_idx, seq in A.iter(corpus):
        start_idx = end_idx - 22  # 23-mer starts 22 chars before end
        if seq in counts:
            counts[seq] += 1
        else:
            counts[seq] = 1
            first_pos[seq] = start_idx

    elapsed = time.time() - t0
    print(f"  Single-pass scan: {elapsed:.1f}s ({corpus_len / elapsed / 1e6:.1f} Mchar/s)")

    # Build FASTA index for position mapping
    t0 = time.time()
    fasta_starts = build_fasta_index(corpus)
    print(f"  Built FASTA index: {len(fasta_starts)} records ({time.time() - t0:.1f}s)")

    # Build targets with real data
    targets = []
    for seq in sorted(counts.keys()):
        occ = counts[seq]
        abs_pos = first_pos[seq]
        genome_pos = genome_relative_position(corpus, fasta_starts, abs_pos)

        targets.append({
            "p": genome_pos,
            "s": seq,
            "m": pam_map.get(seq, "none"),
            "o": occ,
        })

    # Show some stats
    occ_values = sorted(set(t["o"] for t in targets))
    pos_values = sorted(set(t["p"] for t in targets))
    print(f"  {len(targets)} targets found, {len(pos_values)} unique positions, {len(occ_values)} unique occurrence counts")
    if occ_values:
        print(f"  Occurrences range: {min(occ_values):,} – {max(occ_values):,}")
    if pos_values:
        print(f"  Position range: {min(pos_values):,} – {max(pos_values):,}")

    # Sort by: NGG first, then by descending occurrences, then by position
    def sort_key(t):
        pam_rank = 0 if t["m"] == "NGG" else 1
        return (pam_rank, -t["o"], t["p"])

    targets.sort(key=sort_key)

    # Free corpus memory before next pathogen
    del corpus

    return {**meta, "target_count": len(targets), "targets": targets}


def main():
    print("CRISPR Target Recount")
    print("=" * 50)

    output = OrderedDict()
    for key, m in META.items():
        print(f"\nProcessing {m['name']}...")
        result = process_pathogen(key, m)
        if result:
            output[key] = result

    # Write output
    os.makedirs(OUT_JSON.parent, exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(output, f, separators=(",", ":"))

    sz = os.path.getsize(OUT_JSON)
    print(f"\n{'=' * 50}")
    for k, v in output.items():
        print(f"  {k}: {v['target_count']} targets")
    total = sum(v["target_count"] for v in output.values())
    print(f"  Total: {total} targets")
    print(f"  JSON size: {sz / 1024:.1f} KB")
    print(f"  Written to: {OUT_JSON}")


if __name__ == "__main__":
    main()
