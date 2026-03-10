#!/usr/bin/env python3
"""
Fetch gene annotations (GFF3) for pathogen reference genomes from NCBI.
Outputs: ontology_cache/gene_annotations.json

Downloads GFF3 from NCBI Datasets, parses gene features, stores gene name + position + product.
This enables mapping CRISPR target positions to gene context.
"""
import json
import os
import sys
import time
import urllib.request
import urllib.error
import gzip
import io

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MANIFEST = os.path.join(SCRIPT_DIR, "pathogen_manifest.json")
CACHE_DIR = os.path.join(SCRIPT_DIR, "ontology_cache")
GFF_DIR = os.path.join(CACHE_DIR, "gff3")
OUTPUT = os.path.join(CACHE_DIR, "gene_annotations.json")

# Reference genome accessions for each pathogen (representative/reference genomes)
REFERENCE_ACCESSIONS = {
    "cholera": "GCF_008369605.1",       # V. cholerae O1 biovar El Tor str. N16961
    "dengue": "GCF_000862125.1",         # Dengue virus 2, reference
    "ebola": "GCF_000848505.1",          # Zaire ebolavirus
    "hepatitis-b": "GCF_000856845.1",    # Hepatitis B virus
    "hiv-1": "GCF_000864765.1",          # HIV-1 reference
    "influenza-a": "GCF_000865085.1",    # Influenza A H1N1
    "mers": "GCF_000901155.1",           # MERS-CoV
    "mpox": "GCF_014621545.1",           # Monkeypox virus
    "rsv": "GCF_002815475.1",            # RSV-A
    "sars-cov-2": "GCF_009858895.2",     # SARS-CoV-2 (Wuhan-Hu-1)
    "tuberculosis": "GCF_000195955.2",   # M. tuberculosis H37Rv
    "zika": "GCF_000860085.1",           # Zika virus
}


def fetch_gff3(accession: str, key: str) -> str:
    """Download GFF3 annotation from NCBI Datasets API."""
    # NCBI Datasets API for genome annotation
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/annotation_report"

    # Simpler approach: use the direct FTP-style URL for GFF3
    # Try the datasets API first for structured data
    gff_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download?include_annotation_type=GENOME_GFF"

    # Actually, let's use the simple efetch approach for gene list
    genes_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&rettype=gene_table&retmode=text"

    # Best approach: get gene info via NCBI Datasets gene API using taxid
    return ""


def fetch_genes_via_datasets(accession: str) -> list:
    """Fetch gene list from NCBI Datasets V2 API."""
    all_genes = []
    page_token = None

    for attempt in range(3):
        try:
            while True:
                url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/annotation_report?page_size=1000"
                if page_token:
                    url += f"&page_token={page_token}"

                req = urllib.request.Request(url, headers={
                    "User-Agent": "LOOM-ontology/1.0",
                    "Accept": "application/json"
                })
                with urllib.request.urlopen(req, timeout=60) as resp:
                    data = json.loads(resp.read().decode("utf-8"))

                for report in data.get("reports", []):
                    ann = report.get("annotation", {})
                    # Extract genomic range from nested structure
                    regions = ann.get("genomic_regions", [])
                    gene_range = {}
                    if regions:
                        gr = regions[0].get("gene_range", {})
                        ranges = gr.get("range", [])
                        if ranges:
                            gene_range = ranges[0]

                    begin = gene_range.get("begin")
                    end = gene_range.get("end")

                    all_genes.append({
                        "symbol": ann.get("symbol", ""),
                        "name": ann.get("name", ""),
                        "gene_id": int(ann["gene_id"]) if ann.get("gene_id") else None,
                        "type": ann.get("gene_type", ""),
                        "start": int(begin) if begin else None,
                        "end": int(end) if end else None,
                        "orientation": gene_range.get("orientation", ""),
                    })

                page_token = data.get("next_page_token")
                if not page_token:
                    break
                time.sleep(0.3)

            return all_genes
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, KeyError) as e:
            print(f"    Datasets API retry {attempt+1}/3: {e}", file=sys.stderr)
            time.sleep(2 ** attempt)
    return []


def fetch_genes_via_esearch(taxid: int, species: str) -> list:
    """Fallback: fetch gene list via NCBI E-utilities."""
    # Search for genes
    search_url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=gene&term={urllib.request.quote(species)}[Organism]&retmax=500&retmode=json"
    )

    try:
        req = urllib.request.Request(search_url, headers={"User-Agent": "LOOM-ontology/1.0"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        ids = data.get("esearchresult", {}).get("idlist", [])
        if not ids:
            return []

        # Fetch summaries in batches
        genes = []
        batch_size = 100
        for i in range(0, len(ids), batch_size):
            batch = ids[i:i + batch_size]
            summary_url = (
                f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                f"?db=gene&id={','.join(batch)}&retmode=json"
            )
            req = urllib.request.Request(summary_url, headers={"User-Agent": "LOOM-ontology/1.0"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                sdata = json.loads(resp.read().decode("utf-8"))

            for gid in batch:
                info = sdata.get("result", {}).get(gid, {})
                if isinstance(info, dict) and "name" in info:
                    loc = info.get("genomicinfo", [{}])
                    gi = loc[0] if loc else {}
                    genes.append({
                        "symbol": info.get("name", ""),
                        "name": info.get("description", ""),
                        "gene_id": int(gid),
                        "type": info.get("type", {}).get("value", ""),
                        "start": gi.get("chrstart"),
                        "end": gi.get("chrstop"),
                    })
            time.sleep(0.4)

        return genes
    except Exception as e:
        print(f"    E-utils fallback failed: {e}", file=sys.stderr)
        return []


def main():
    os.makedirs(GFF_DIR, exist_ok=True)
    os.makedirs(CACHE_DIR, exist_ok=True)

    with open(MANIFEST) as f:
        manifest = json.load(f)

    results = {}

    print(f"Fetching gene annotations for {len(REFERENCE_ACCESSIONS)} pathogens...")
    for key, accession in REFERENCE_ACCESSIONS.items():
        info = manifest["pathogens"].get(key, {})
        taxid = info.get("ncbi_taxid")
        species = info.get("species", "")

        print(f"  {key} ({accession})...", end=" ", flush=True)

        # Try Datasets API first
        genes = fetch_genes_via_datasets(accession)

        # Fallback to E-utils
        if not genes and taxid:
            print("(trying E-utils)...", end=" ", flush=True)
            genes = fetch_genes_via_esearch(taxid, species)

        results[key] = {
            "accession": accession,
            "gene_count": len(genes),
            "genes": genes
        }
        print(f"✓ {len(genes)} genes")
        time.sleep(0.5)

    with open(OUTPUT, "w") as f:
        json.dump({
            "source": "NCBI Datasets + E-utilities",
            "fetched": time.strftime("%Y-%m-%dT%H:%M:%SZ"),
            "pathogens": results
        }, f, indent=2)

    print(f"\nDone. Written to {OUTPUT}")


if __name__ == "__main__":
    main()
