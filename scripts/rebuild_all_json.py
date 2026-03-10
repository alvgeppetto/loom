#!/usr/bin/env python3
"""Rebuild crispr-targets.json and animal-targets.json from all CSV files.

Reads all *_targets.csv from CRISPR_DIR, categorizes them as human-pathogen
or animal-reference, and writes two JSON files for the web app.
"""
import csv
import json
import os
from collections import OrderedDict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IDX_DIR = Path("/Volumes/T9/loom-openscience/indexes")
CSV_DIR = Path("/Volumes/T9/loom-openscience/crispr-db")
RAW_DIR = Path("/Volumes/T9/loom-openscience/raw")
OUT_HUMAN = ROOT / "web" / "data" / "crispr-targets.json"
OUT_ANIMAL = ROOT / "web" / "data" / "animal-targets.json"

# Pathogen metadata: name, family, type (human-pathogen or animal-ref)
METADATA = {
    # Existing indexed
    "dengue":       {"name": "Dengue Virus",       "family": "Flaviviridae",     "type": "human"},
    "ebola":        {"name": "Ebola Virus",         "family": "Filoviridae",      "type": "human"},
    "hepatitis-b":  {"name": "Hepatitis B",         "family": "Hepadnaviridae",   "type": "human"},
    "mers":         {"name": "MERS-CoV",            "family": "Coronaviridae",    "type": "human"},
    "refseq-viral": {"name": "RefSeq Viral (all)",  "family": "Mixed",            "type": "human"},
    "rsv":          {"name": "RSV",                 "family": "Pneumoviridae",    "type": "human"},
    "zika":         {"name": "Zika Virus",          "family": "Flaviviridae",     "type": "human"},
    "influenza-a":  {"name": "Influenza A",         "family": "Orthomyxoviridae", "type": "human"},
    "human-grch38": {"name": "Human GRCh38",        "family": "Homo sapiens (off-target)", "type": "human"},
    # New pathogens
    "hiv-1":        {"name": "HIV-1",               "family": "Retroviridae",     "type": "human"},
    "mpox":         {"name": "Mpox",                "family": "Poxviridae",       "type": "human"},
    "sars-cov-2":   {"name": "SARS-CoV-2",          "family": "Coronaviridae",    "type": "human"},
    "tuberculosis": {"name": "M. tuberculosis",     "family": "Mycobacteriaceae", "type": "human"},
    "cholera":      {"name": "V. cholerae",         "family": "Vibrionaceae",     "type": "human"},
    "hepatitis-c":  {"name": "Hepatitis C",         "family": "Flaviviridae",     "type": "human"},
    "hpv":          {"name": "HPV",                 "family": "Papillomaviridae", "type": "human"},
    "norovirus":    {"name": "Norovirus",           "family": "Caliciviridae",    "type": "human"},
    "rabies":       {"name": "Rabies Virus",        "family": "Rhabdoviridae",    "type": "human"},
    "measles":      {"name": "Measles Virus",       "family": "Paramyxoviridae",  "type": "human"},
    "rotavirus-a":  {"name": "Rotavirus A",         "family": "Reoviridae",       "type": "human"},
    "chikungunya":  {"name": "Chikungunya",         "family": "Togaviridae",      "type": "human"},
    "west-nile":    {"name": "West Nile Virus",     "family": "Flaviviridae",     "type": "human"},
    "yellow-fever": {"name": "Yellow Fever",        "family": "Flaviviridae",     "type": "human"},
    "marburg":      {"name": "Marburg Virus",       "family": "Filoviridae",      "type": "human"},
    "nipah":        {"name": "Nipah Virus",         "family": "Paramyxoviridae",  "type": "human"},
    "lassa":        {"name": "Lassa Virus",         "family": "Arenaviridae",     "type": "human"},
    "plasmodium-falciparum": {"name": "P. falciparum", "family": "Plasmodiidae",  "type": "human"},
    # Animal reference genomes
    "bat-rousettus": {"name": "Egyptian Fruit Bat",    "family": "Pteropodidae (Ebola reservoir)",    "type": "animal"},
    "chicken":       {"name": "Chicken (Gallus gallus)","family": "Phasianidae (Avian flu host)",     "type": "animal"},
    "pig":           {"name": "Pig (Sus scrofa)",       "family": "Suidae (Swine flu, Nipah host)",   "type": "animal"},
    "cow":           {"name": "Cow (Bos taurus)",       "family": "Bovidae (Bovine TB, MERS-like)",   "type": "animal"},
    "camel":         {"name": "Dromedary Camel",        "family": "Camelidae (MERS host)",            "type": "animal"},
    "mouse":         {"name": "Mouse (Mus musculus)",    "family": "Muridae (Lab model)",             "type": "animal"},
}


def count_genomes(name: str) -> int:
    """Count FASTA records in raw genome file (searches recursively)."""
    raw = RAW_DIR / name
    for ext in ("fna", "fa", "fasta"):
        for f in sorted(raw.rglob(f"*.{ext}"), key=lambda p: p.stat().st_size, reverse=True):
            if f.stat().st_size > 0:
                count = 0
                with open(f) as fh:
                    for line in fh:
                        if line.startswith(">"):
                            count += 1
                return count
    return 0


def process_pathogen(key: str) -> dict | None:
    """Read CSV and build JSON entry for one pathogen/genome."""
    meta = METADATA.get(key)
    if not meta:
        return None
    
    csv_path = CSV_DIR / f"{key}_targets.csv"
    idx_path = IDX_DIR / f"{key}.idx"
    
    if not csv_path.exists():
        print(f"  SKIP {key}: no CSV")
        return None
    
    idx_mb = round(idx_path.stat().st_size / (1024*1024)) if idx_path.exists() else 0
    
    targets = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            pam = row["pam_type"]
            if pam == "marker":
                pam = "none"
            targets.append({
                "s": row["sequence_23mer"],
                "p": int(row["position"]),
                "m": pam,
                "o": int(row["occurrences"]),
            })
    
    genomes = count_genomes(key)
    if genomes == 0:
        genomes = 1  # At least 1 for reference genomes
    
    return {
        "name": meta["name"],
        "family": meta["family"],
        "genomes": genomes,
        "target_count": len(targets),
        "index_mb": idx_mb,
        "targets": targets,
    }


def main():
    human_data = {}
    animal_data = {}
    
    # Process all CSVs that exist
    for csv_file in sorted(CSV_DIR.glob("*_targets.csv")):
        key = csv_file.stem.replace("_targets", "")
        meta = METADATA.get(key)
        if not meta:
            print(f"  UNKNOWN: {key} (not in metadata, skipping)")
            continue
        
        entry = process_pathogen(key)
        if entry:
            if meta["type"] == "animal":
                animal_data[key] = entry
            else:
                human_data[key] = entry
            print(f"  OK {key}: {entry['target_count']} targets, {entry['genomes']} genomes")
    
    # Write human pathogens JSON
    with open(OUT_HUMAN, 'w') as f:
        json.dump(human_data, f, separators=(',', ':'))
    print(f"\nHuman JSON: {len(human_data)} pathogens → {OUT_HUMAN} ({OUT_HUMAN.stat().st_size / 1024:.0f} KB)")
    
    # Write animal references JSON
    with open(OUT_ANIMAL, 'w') as f:
        json.dump(animal_data, f, separators=(',', ':'))
    print(f"Animal JSON: {len(animal_data)} references → {OUT_ANIMAL} ({OUT_ANIMAL.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
