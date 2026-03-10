#!/usr/bin/env python3
"""
Fetch NCBI Taxonomy lineage for each pathogen TaxID.
Outputs: ontology_cache/ncbi_taxonomy.json

Uses NCBI E-utilities (free, no API key needed at low volume).
Rate limit: 3 requests/sec without key, 10/sec with.
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
OUTPUT = os.path.join(CACHE_DIR, "ncbi_taxonomy.json")

EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def fetch_lineage(taxid: int) -> dict:
    """Fetch full taxonomy lineage for a given TaxID."""
    url = f"{EFETCH_URL}?db=taxonomy&id={taxid}&retmode=xml"
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "LOOM-ontology/1.0"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                xml = resp.read().decode("utf-8")
            return parse_taxon_xml(xml, taxid)
        except (urllib.error.URLError, TimeoutError) as e:
            print(f"  Retry {attempt+1}/3 for TaxID {taxid}: {e}", file=sys.stderr)
            time.sleep(2 ** attempt)
    return {"taxid": taxid, "error": "fetch_failed"}


def parse_taxon_xml(xml: str, taxid: int) -> dict:
    """Minimal XML parsing — extract key fields without lxml dependency."""
    def extract_tag(tag: str, text: str) -> str:
        start = text.find(f"<{tag}>")
        if start == -1:
            return ""
        start += len(tag) + 2
        end = text.find(f"</{tag}>", start)
        return text[start:end].strip() if end != -1 else ""

    result = {
        "taxid": taxid,
        "scientific_name": extract_tag("ScientificName", xml),
        "rank": extract_tag("Rank", xml),
        "division": extract_tag("Division", xml),
        "lineage": extract_tag("Lineage", xml),
    }

    # Extract lineage items with ranks
    lineage_items = []
    pos = 0
    while True:
        start = xml.find("<Taxon>", pos + 1 if pos else 0)
        if start == -1 or start == pos:
            break
        end = xml.find("</Taxon>", start)
        if end == -1:
            break
        chunk = xml[start:end]
        name = extract_tag("ScientificName", chunk)
        rank = extract_tag("Rank", chunk)
        tid = extract_tag("TaxId", chunk)
        if name and rank and rank != "no rank":
            lineage_items.append({"name": name, "rank": rank, "taxid": int(tid) if tid.isdigit() else None})
        pos = end

    result["lineage_structured"] = lineage_items
    return result


def main():
    os.makedirs(CACHE_DIR, exist_ok=True)

    with open(MANIFEST) as f:
        manifest = json.load(f)

    results = {}
    all_items = list(manifest["pathogens"].items()) + list(manifest["animals"].items())

    print(f"Fetching NCBI Taxonomy for {len(all_items)} organisms...")
    for key, info in all_items:
        taxid = info["ncbi_taxid"]
        print(f"  {key} (TaxID: {taxid})...", end=" ", flush=True)
        data = fetch_lineage(taxid)
        results[key] = data
        print(f"✓ {data.get('scientific_name', '?')}")
        time.sleep(0.4)  # respect rate limit

    with open(OUTPUT, "w") as f:
        json.dump({"source": "NCBI Taxonomy", "fetched": time.strftime("%Y-%m-%dT%H:%M:%SZ"), "organisms": results}, f, indent=2)

    print(f"\nDone. Written to {OUTPUT}")


if __name__ == "__main__":
    main()
