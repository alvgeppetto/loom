#!/usr/bin/env python3
"""
Fetch MONDO disease terms for each pathogen.
Outputs: ontology_cache/mondo.json

Uses EBI OLS4 API (free, no key needed).
"""
import json
import os
import sys
import time
import urllib.request
import urllib.error

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MANIFEST = os.path.join(SCRIPT_DIR, "pathogen_manifest.json")
CACHE_DIR = os.path.join(SCRIPT_DIR, "ontology_cache")
OUTPUT = os.path.join(CACHE_DIR, "mondo.json")

OLS4_BASE = "https://www.ebi.ac.uk/ols4/api/ontologies/mondo/terms"


def fetch_mondo(mondo_id: str) -> dict:
    """Fetch MONDO term from OLS4."""
    encoded = mondo_id.replace(":", "_")
    iri = f"http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F{encoded}"
    url = f"{OLS4_BASE}/{iri}"

    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "LOOM-ontology/1.0", "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read().decode("utf-8"))
            return {
                "id": mondo_id,
                "name": data.get("label", ""),
                "definition": (data.get("description") or [""])[0],
                "synonyms": data.get("synonyms", []),
                "xrefs": data.get("annotation", {}).get("database_cross_reference", []),
                "exact_mappings": data.get("annotation", {}).get("exactMatch", []),
            }
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            print(f"  Retry {attempt+1}/3 for {mondo_id}: {e}", file=sys.stderr)
            time.sleep(2 ** attempt)
    return {"id": mondo_id, "error": "fetch_failed"}


def fetch_mondo_ancestors(mondo_id: str) -> list:
    """Fetch ancestor terms (hierarchy) from OLS4."""
    encoded = mondo_id.replace(":", "_")
    iri = f"http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F{encoded}"
    url = f"{OLS4_BASE}/{iri}/ancestors?size=50"

    try:
        req = urllib.request.Request(url, headers={"User-Agent": "LOOM-ontology/1.0", "Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        terms = data.get("_embedded", {}).get("terms", [])
        return [{"label": t.get("label", ""), "obo_id": t.get("obo_id", "")} for t in terms if t.get("obo_id")]
    except Exception:
        return []


def main():
    os.makedirs(CACHE_DIR, exist_ok=True)

    with open(MANIFEST) as f:
        manifest = json.load(f)

    results = {}
    pathogens = manifest["pathogens"]

    print(f"Fetching MONDO terms for {len(pathogens)} pathogens...")
    for key, info in pathogens.items():
        mondo_id = info.get("mondo")
        if not mondo_id:
            print(f"  {key}: no MONDO ID, skipping")
            continue

        print(f"  {key} ({mondo_id})...", end=" ", flush=True)
        data = fetch_mondo(mondo_id)
        data["ancestors"] = fetch_mondo_ancestors(mondo_id)
        results[key] = data
        print(f"✓ {data.get('name', '?')}")
        time.sleep(0.5)

    with open(OUTPUT, "w") as f:
        json.dump({"source": "MONDO via OLS4", "fetched": time.strftime("%Y-%m-%dT%H:%M:%SZ"), "diseases": results}, f, indent=2)

    print(f"\nDone. Written to {OUTPUT}")


if __name__ == "__main__":
    main()
