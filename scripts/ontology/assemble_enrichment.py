#!/usr/bin/env python3
"""
Assemble ontology enrichment data from cached fetcher outputs.
Reads: ontology_cache/*.json + pathogen_manifest.json
Writes: web/data/ontology-enrichment.json

This file is loaded by the web UI as optional sidecar data.
Zero dependency on any external service at runtime — fully offline.
"""
import json
import os
import sys
import time

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.join(SCRIPT_DIR, "ontology_cache")
LOOM_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
MANIFEST = os.path.join(SCRIPT_DIR, "pathogen_manifest.json")
OUTPUT = os.path.join(LOOM_ROOT, "web", "data", "ontology-enrichment.json")


def load_cache(filename: str) -> dict:
    path = os.path.join(CACHE_DIR, filename)
    if not os.path.exists(path):
        print(f"  WARNING: {filename} not found, skipping")
        return {}
    with open(path) as f:
        return json.load(f)


def main():
    print("Assembling ontology enrichment data...")

    # Load manifest
    with open(MANIFEST) as f:
        manifest = json.load(f)

    # Load cached data
    taxonomy = load_cache("ncbi_taxonomy.json")
    disease_ontology = load_cache("disease_ontology.json")
    mondo = load_cache("mondo.json")
    who_context = load_cache("who_context.json")
    gene_annotations = load_cache("gene_annotations.json")

    # Build enrichment for each pathogen
    enrichment = {}
    for key, info in manifest["pathogens"].items():
        entry = {
            # Core identifiers from manifest
            "species": info["species"],
            "ncbi_taxid": info["ncbi_taxid"],
            "doid": info.get("doid"),
            "mondo_id": info.get("mondo"),
            "icd11": info.get("icd11"),
            "transmission": info.get("transmission"),
            "genome_type": info.get("genome_type"),
            "who_priority": info.get("who_priority", False),
        }

        # NCBI Taxonomy enrichment
        tax_data = taxonomy.get("organisms", {}).get(key, {})
        if tax_data and "error" not in tax_data:
            entry["taxonomy"] = {
                "scientific_name": tax_data.get("scientific_name", ""),
                "rank": tax_data.get("rank", ""),
                "division": tax_data.get("division", ""),
                "lineage": tax_data.get("lineage", ""),
                "lineage_structured": tax_data.get("lineage_structured", []),
            }

        # Disease Ontology enrichment
        do_data = disease_ontology.get("diseases", {}).get(key, {})
        if do_data and "error" not in do_data:
            entry["disease_ontology"] = {
                "doid": do_data.get("doid", ""),
                "name": do_data.get("name", ""),
                "definition": do_data.get("definition", ""),
                "synonyms": do_data.get("synonyms", []),
                "xrefs": do_data.get("xrefs", []),
            }

        # MONDO enrichment
        mondo_data = mondo.get("diseases", {}).get(key, {})
        if mondo_data and "error" not in mondo_data:
            entry["mondo"] = {
                "id": mondo_data.get("id", ""),
                "name": mondo_data.get("name", ""),
                "definition": mondo_data.get("definition", ""),
                "synonyms": mondo_data.get("synonyms", []),
                "ancestors": mondo_data.get("ancestors", []),
            }

        # WHO context
        who_data = who_context.get("pathogens", {}).get(key, {})
        if who_data:
            entry["who_context"] = {
                "classification": who_data.get("who_classification", ""),
                "case_fatality_rate": who_data.get("case_fatality_rate", ""),
                "annual_cases": who_data.get("annual_cases_est", ""),
                "annual_deaths": who_data.get("annual_deaths_est", ""),
                "geographic_spread": who_data.get("geographic_spread", ""),
                "diagnostic_need": who_data.get("diagnostic_need", ""),
                "existing_rapid_tests": who_data.get("existing_rapid_tests", []),
                "crispr_dx_status": who_data.get("crispr_diagnostic_status", ""),
                "vaccine_available": who_data.get("vaccine_available"),
            }
            if who_data.get("gho_supplement"):
                entry["who_context"]["gho_data"] = who_data["gho_supplement"]

        # Gene annotations
        gene_data = gene_annotations.get("pathogens", {}).get(key, {})
        if gene_data and gene_data.get("gene_count", 0) > 0:
            entry["genes"] = {
                "reference_accession": gene_data.get("accession", ""),
                "gene_count": gene_data.get("gene_count", 0),
                "genes": gene_data.get("genes", []),
            }

        enrichment[key] = entry

    # Animal enrichment (lighter)
    animals = {}
    for key, info in manifest.get("animals", {}).items():
        tax_data = taxonomy.get("organisms", {}).get(key, {})
        animals[key] = {
            "species": info["species"],
            "ncbi_taxid": info["ncbi_taxid"],
        }
        if tax_data and "error" not in tax_data:
            animals[key]["taxonomy"] = {
                "scientific_name": tax_data.get("scientific_name", ""),
                "lineage": tax_data.get("lineage", ""),
            }

    # Compute summary stats
    n_with_taxonomy = sum(1 for v in enrichment.values() if "taxonomy" in v)
    n_with_disease = sum(1 for v in enrichment.values() if "disease_ontology" in v)
    n_with_mondo = sum(1 for v in enrichment.values() if "mondo" in v)
    n_with_who = sum(1 for v in enrichment.values() if "who_context" in v)
    n_with_genes = sum(1 for v in enrichment.values() if "genes" in v)
    total_genes = sum(v.get("genes", {}).get("gene_count", 0) for v in enrichment.values())

    output = {
        "version": 1,
        "generated": time.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "summary": {
            "pathogens": len(enrichment),
            "animals": len(animals),
            "with_taxonomy": n_with_taxonomy,
            "with_disease_ontology": n_with_disease,
            "with_mondo": n_with_mondo,
            "with_who_context": n_with_who,
            "with_gene_annotations": n_with_genes,
            "total_genes_annotated": total_genes,
        },
        "ontology_sources": [
            {"name": "NCBI Taxonomy", "url": "https://www.ncbi.nlm.nih.gov/taxonomy", "license": "Public Domain"},
            {"name": "Disease Ontology", "url": "https://disease-ontology.org/", "license": "CC0 1.0"},
            {"name": "MONDO", "url": "https://mondo.monarchinitiative.org/", "license": "CC BY 4.0"},
            {"name": "WHO GHO", "url": "https://www.who.int/data/gho", "license": "CC BY-NC-SA 3.0 IGO"},
            {"name": "ICD-11", "url": "https://icd.who.int/", "license": "WHO"},
        ],
        "pathogens": enrichment,
        "animals": animals,
    }

    os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
    with open(OUTPUT, "w") as f:
        json.dump(output, f, indent=2)

    size_kb = os.path.getsize(OUTPUT) / 1024
    print(f"\nAssembly complete:")
    print(f"  Pathogens: {len(enrichment)} ({n_with_taxonomy} taxonomy, {n_with_disease} DO, {n_with_mondo} MONDO, {n_with_who} WHO, {n_with_genes} genes)")
    print(f"  Animals: {len(animals)}")
    print(f"  Total genes: {total_genes}")
    print(f"  Output: {OUTPUT} ({size_kb:.1f} KB)")


if __name__ == "__main__":
    main()
