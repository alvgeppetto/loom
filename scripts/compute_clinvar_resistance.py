#!/usr/bin/env python3
"""Map ClinVar drug-resistance and diagnostic-relevant variants to CRISPR targets.

For each pathogen with available resistance data, this script:
1. Downloads or reads pre-curated resistance mutation coordinates
2. Maps them to overlapping CRISPR guide target positions
3. Produces web/data/clinvar-resistance.json

Note: ClinVar primarily covers human variants. For pathogen resistance mutations
we use curated databases (WHO catalogs for TB, Stanford HIVDB, etc.) rather than
ClinVar directly. We map variant positions onto our guide coordinate space.
"""
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

CRISPR_DB = Path("/Volumes/T9/loom-openscience/crispr-db")
OUT = Path("web/data/clinvar-resistance.json")

# Curated resistance mutation databases (position ranges in reference genome coordinates)
# Sources: WHO mutation catalog (TB), Stanford HIVDB (HIV), literature reviews
RESISTANCE_DATA = {
    "tuberculosis": {
        "source": "WHO Catalogue of Mutations (2021)",
        "source_url": "https://www.who.int/publications/i/item/9789240028173",
        "mutations": [
            {"gene": "rpoB", "name": "Rifampicin resistance (RRDR)", "start": 759807, "end": 760353, "drug": "Rifampicin", "significance": "high"},
            {"gene": "katG", "name": "Isoniazid resistance (S315T)", "start": 2155168, "end": 2156111, "drug": "Isoniazid", "significance": "high"},
            {"gene": "inhA", "name": "Isoniazid/Ethionamide resistance", "start": 1674202, "end": 1675011, "drug": "Isoniazid/Ethionamide", "significance": "high"},
            {"gene": "embB", "name": "Ethambutol resistance", "start": 4247429, "end": 4249810, "drug": "Ethambutol", "significance": "high"},
            {"gene": "pncA", "name": "Pyrazinamide resistance", "start": 2288681, "end": 2289241, "drug": "Pyrazinamide", "significance": "high"},
            {"gene": "gyrA", "name": "Fluoroquinolone resistance (QRDR)", "start": 7302, "end": 9818, "drug": "Fluoroquinolones", "significance": "high"},
            {"gene": "gyrB", "name": "Fluoroquinolone resistance", "start": 5240, "end": 7267, "drug": "Fluoroquinolones", "significance": "moderate"},
            {"gene": "rrs", "name": "Aminoglycoside resistance", "start": 1471846, "end": 1473382, "drug": "Streptomycin/Kanamycin", "significance": "high"},
            {"gene": "eis", "name": "Kanamycin resistance (promoter)", "start": 2714124, "end": 2715332, "drug": "Kanamycin", "significance": "moderate"},
        ]
    },
    "hiv-1": {
        "source": "Stanford HIV Drug Resistance Database (HIVDB)",
        "source_url": "https://hivdb.stanford.edu/",
        "mutations": [
            {"gene": "RT", "name": "Reverse transcriptase NRTI resistance", "start": 2550, "end": 3869, "drug": "NRTIs (AZT, 3TC, TDF)", "significance": "high"},
            {"gene": "RT", "name": "Reverse transcriptase NNRTI resistance", "start": 2550, "end": 3869, "drug": "NNRTIs (EFV, NVP)", "significance": "high"},
            {"gene": "PR", "name": "Protease inhibitor resistance", "start": 2253, "end": 2549, "drug": "PIs (LPV, ATV, DRV)", "significance": "high"},
            {"gene": "IN", "name": "Integrase inhibitor resistance", "start": 4230, "end": 5096, "drug": "INSTIs (DTG, RAL, BIC)", "significance": "high"},
        ]
    },
    "influenza-a": {
        "source": "WHO Antiviral Resistance Monitoring",
        "source_url": "https://www.who.int/teams/global-influenza-programme",
        "mutations": [
            {"gene": "NA", "name": "Neuraminidase inhibitor resistance (H275Y)", "start": 0, "end": 1410, "drug": "Oseltamivir/Zanamivir", "significance": "high"},
            {"gene": "M2", "name": "M2 ion channel blocker resistance (S31N)", "start": 0, "end": 294, "drug": "Amantadine/Rimantadine", "significance": "high"},
            {"gene": "PA", "name": "Baloxavir resistance (I38T)", "start": 0, "end": 2151, "drug": "Baloxavir", "significance": "moderate"},
        ]
    },
    "hepatitis-b": {
        "source": "EASL Clinical Practice Guidelines (2017)",
        "source_url": "https://www.journal-of-hepatology.eu/article/S0168-8278(17)30185-X",
        "mutations": [
            {"gene": "P (polymerase)", "name": "Nucleos(t)ide analogue resistance", "start": 2307, "end": 3215, "drug": "Entecavir/Tenofovir/Lamivudine", "significance": "high"},
            {"gene": "preS/S", "name": "Vaccine escape mutations", "start": 2848, "end": 835, "drug": "Vaccine/HBIg", "significance": "moderate"},
        ]
    },
}


def find_overlapping_targets(targets: list, mutation: dict) -> list:
    """Find CRISPR targets whose position overlaps a resistance mutation region."""
    overlapping = []
    m_start = mutation["start"]
    m_end = mutation["end"]

    # Handle circular genomes (end < start)
    if m_end < m_start:
        for t in targets:
            pos = int(t["position"])
            t_end = pos + 23
            if pos >= m_start or t_end <= m_end:
                overlapping.append(t)
    else:
        for t in targets:
            pos = int(t["position"])
            t_end = pos + 23
            # Overlap check: target [pos, pos+23) intersects [m_start, m_end)
            if pos < m_end and t_end > m_start:
                overlapping.append(t)

    return overlapping


def main():
    result = {
        "generated": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "version": "1.0",
        "description": "Drug-resistance and diagnostic-relevant variant regions mapped to CRISPR guide target coordinates",
        "pathogens": {},
    }

    total_mutations = 0
    total_overlapping_guides = 0

    for pathogen, rdata in RESISTANCE_DATA.items():
        csv_path = CRISPR_DB / f"{pathogen}_targets.csv"
        if not csv_path.exists():
            print(f"  WARN: {csv_path} not found, skipping")
            continue

        # Load all targets
        targets = []
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                targets.append(row)

        # Map mutations to overlapping targets
        mutations_mapped = []
        for mutation in rdata["mutations"]:
            overlaps = find_overlapping_targets(targets, mutation)
            ngg_overlaps = [t for t in overlaps if t.get("pam_type") == "NGG"]

            m_result = {
                "gene": mutation["gene"],
                "name": mutation["name"],
                "drug": mutation["drug"],
                "significance": mutation["significance"],
                "region": f"{mutation['start']}-{mutation['end']}",
                "overlapping_guides": len(overlaps),
                "ngg_overlapping_guides": len(ngg_overlaps),
                "top_guides": [
                    {"seq": t["sequence_23mer"], "pos": int(t["position"]),
                     "pam": t.get("pam_type", "none"), "occ": int(t.get("occurrences", 0))}
                    for t in sorted(ngg_overlaps if ngg_overlaps else overlaps,
                                    key=lambda x: -int(x.get("occurrences", 0)))[:5]
                ],
            }
            mutations_mapped.append(m_result)
            total_mutations += 1
            total_overlapping_guides += len(overlaps)

        result["pathogens"][pathogen] = {
            "source": rdata["source"],
            "source_url": rdata["source_url"],
            "mutations": mutations_mapped,
            "total_resistance_regions": len(mutations_mapped),
            "total_overlapping_guides": sum(m["overlapping_guides"] for m in mutations_mapped),
        }

        print(f"  {pathogen}: {len(mutations_mapped)} resistance regions, "
              f"{sum(m['overlapping_guides'] for m in mutations_mapped)} overlapping guides")

    result["summary"] = {
        "pathogens_with_resistance_data": len(result["pathogens"]),
        "total_resistance_regions": total_mutations,
        "total_overlapping_guides": total_overlapping_guides,
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(result, f, indent=2)

    print(f"\n  Output: {OUT} ({OUT.stat().st_size / 1024:.1f} KB)")
    print(f"  {len(result['pathogens'])} pathogens, {total_mutations} resistance regions, "
          f"{total_overlapping_guides} overlapping guides")


if __name__ == "__main__":
    main()
