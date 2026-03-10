#!/usr/bin/env python3
"""
Fetch Disease Ontology (DOID) data for each pathogen's disease.
Outputs: ontology_cache/disease_ontology.json

Uses the Disease Ontology API (free, no key needed).
Fetches: disease name, definition, synonyms, cross-references (xrefs to MONDO, ICD, MESH, etc.)
"""
import json
import os
import re
import sys
import time
import urllib.request
import urllib.error

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MANIFEST = os.path.join(SCRIPT_DIR, "pathogen_manifest.json")
CACHE_DIR = os.path.join(SCRIPT_DIR, "ontology_cache")
OUTPUT = os.path.join(CACHE_DIR, "disease_ontology.json")

DO_API = "https://www.disease-ontology.org/api"


def fetch_doid(doid: str) -> dict:
    """Fetch a single DOID term from the Disease Ontology API."""
    # DO API: /api/metadata/DOID:1498
    doid_num = doid.replace("DOID:", "")
    url = f"{DO_API}/metadata/DOID:{doid_num}"

    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "LOOM-ontology/1.0", "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except (urllib.error.URLError, TimeoutError) as e:
            print(f"  Retry {attempt+1}/3 for {doid}: {e}", file=sys.stderr)
            time.sleep(2 ** attempt)
        except json.JSONDecodeError:
            # Try OLS as fallback
            return fetch_doid_ols(doid)
    return {"id": doid, "error": "fetch_failed"}


def fetch_doid_ols(doid: str) -> dict:
    """Fallback: fetch from EBI OLS (Ontology Lookup Service)."""
    encoded = doid.replace(":", "_")
    url = f"https://www.ebi.ac.uk/ols4/api/ontologies/doid/terms/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F{encoded}"

    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "LOOM-ontology/1.0", "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read().decode("utf-8"))
            return {
                "id": doid,
                "name": data.get("label", ""),
                "definition": (data.get("description") or [""])[0],
                "synonyms": data.get("synonyms", []),
                "xrefs": data.get("annotation", {}).get("database_cross_reference", []),
                "source": "OLS4"
            }
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            print(f"  OLS retry {attempt+1}/3 for {doid}: {e}", file=sys.stderr)
            time.sleep(2 ** attempt)
    return {"id": doid, "error": "fetch_failed"}


def main():
    os.makedirs(CACHE_DIR, exist_ok=True)

    with open(MANIFEST) as f:
        manifest = json.load(f)

    results = {}
    pathogens = manifest["pathogens"]

    print(f"Fetching Disease Ontology for {len(pathogens)} pathogens...")
    for key, info in pathogens.items():
        doid = info.get("doid")
        if not doid:
            print(f"  {key}: no DOID, skipping")
            continue

        print(f"  {key} ({doid})...", end=" ", flush=True)
        data = fetch_doid(doid)

        # Normalize to consistent shape
        raw_synonyms = data.get("synonyms", data.get("synonym", [])) or []
        # Strip OBO qualifiers like "EXACT []", "RELATED OMO:... []", etc.
        clean_synonyms = [re.sub(r'\s+(EXACT|RELATED|NARROW|BROAD)\s*(?:\S+\s*)?\[\]$', '', s).strip()
                          for s in raw_synonyms if isinstance(s, str)]
        # Deduplicate while preserving order
        seen = set()
        deduped = []
        for s in clean_synonyms:
            if s and s.lower() not in seen:
                seen.add(s.lower())
                deduped.append(s)

        results[key] = {
            "doid": doid,
            "name": data.get("name", data.get("label", "")),
            "definition": data.get("definition", data.get("def", "")),
            "synonyms": deduped,
            "xrefs": data.get("xrefs", data.get("xref", [])),
        }
        print(f"✓ {results[key]['name'] or '(fetched)'}")
        time.sleep(0.5)

    with open(OUTPUT, "w") as f:
        json.dump({"source": "Disease Ontology + OLS4", "fetched": time.strftime("%Y-%m-%dT%H:%M:%SZ"), "diseases": results}, f, indent=2)

    print(f"\nDone. Written to {OUTPUT}")


if __name__ == "__main__":
    main()
