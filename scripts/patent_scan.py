#!/usr/bin/env python3
"""
Patent Landscape Scan for CRISPR Targets.

Queries Google Patents via the fetch_webpage tool (JS rendering required)
or the Lens.org API (if LENS_TOKEN is set) for patent family counts.

NOTE: Python urllib cannot extract result counts from Google Patents because
the page is JavaScript-rendered. The scan was executed using the fetch_webpage
tool in a Copilot conversation on 2026-03-07. The compiled results are in:
  data/crispr_guides/patent_scan_report.json

This script documents the methodology and can be rerun if API access is available.

Generates a structured report:
  1. Patent landscape per pathogen (CRISPR + pathogen)
  2. Gene-level counts for SARS-CoV-2 target regions
  3. 20-mer guide sequence searches

Output: data/crispr_guides/patent_scan_report.json

Usage:
  python scripts/patent_scan.py                     # uses Lens API if LENS_TOKEN set
  LENS_TOKEN=xxx python scripts/patent_scan.py      # explicit token
"""

import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data" / "crispr_guides"
NOVELTY_FILE = DATA_DIR / "novelty_report.json"
OUTPUT_FILE = DATA_DIR / "patent_scan_report.json"

LENS_API = "https://api.lens.org/patent/search"
LENS_TOKEN = os.environ.get("LENS_TOKEN", "")

# Rate limit: be polite
DELAY = 1.5  # seconds between requests

# ── Pathogen registry (ontology-expanded per ontology-first SKILL) ──────
#
# Synonyms sourced from web/data/ontology-enrichment.json (NCBI Taxonomy,
# Disease Ontology, MONDO) plus manual additions.  Ambiguous abbreviations
# like bare "RSV" (Rous sarcoma virus) and bare "TB" (terabyte) are excluded
# to avoid false matches in full-text patent search.
#
# Re-scanned 2026-03-08 with ontology expansion via Google Patents.

ONTOLOGY_FILE = ROOT / "web" / "data" / "ontology-enrichment.json"

PATHOGENS = {
    "SARS-CoV-2": {
        "query": ('"SARS-CoV-2" OR "COVID-19" OR "2019-nCoV" OR "COVID19"'
                  ' OR "SARS-coronavirus 2"'
                  ' OR "severe acute respiratory syndrome coronavirus 2"'),
    },
    "Dengue": {
        "query": '"Dengue" OR "Dengue virus" OR "DENV" OR "dengue disease"',
    },
    "Ebola": {
        "query": ('"Ebola" OR "Zaire ebolavirus" OR "EHF"'
                  ' OR "Ebola hemorrhagic fever"'),
    },
    "Hepatitis B": {
        "query": ('"Hepatitis B" OR "HBV" OR "Hepatitis B virus"'
                  ' OR "serum hepatitis"'),
    },
    "MERS": {
        "query": ('"MERS" OR "MERS-CoV"'
                  ' OR "Middle East respiratory syndrome" OR "camel flu"'),
    },
    "RSV": {
        "query": ('"respiratory syncytial virus"'
                  ' OR "Human respiratory syncytial virus" OR "HRSV"'),
        # NOTE: bare "RSV" excluded — ambiguous (Rous sarcoma virus etc.)
    },
    "Zika": {
        "query": '"Zika" OR "Zika virus" OR "ZikV" OR "Zika fever"',
    },
    "HIV": {
        "query": ('"HIV" OR "human immunodeficiency virus"'
                  ' OR "Human immunodeficiency virus 1"'),
    },
    "Mpox": {
        "query": '"Mpox" OR "Monkeypox" OR "Monkeypox virus"',
    },
    "Tuberculosis": {
        "query": ('"tuberculosis" OR "Mycobacterium tuberculosis"'
                  ' OR "Kochs disease"'),
        # NOTE: bare "TB" excluded — ambiguous (terabyte etc.)
    },
    "Influenza A": {
        "query": '"Influenza A" OR "Influenza A virus"',
    },
    "Cholera": {
        "query": ('"Vibrio cholerae" OR "cholera"'
                  ' OR "Vibrio cholerae infection"'),
    },
}

# Gene targets for SARS-CoV-2 (our 52 guides target these)
SARS_GENES = [
    "ORF1a", "ORF1b", "RdRp", "nsp12", "nsp13", "nsp14", "nsp15",
    "spike", "nucleocapsid", "ORF3a", "ORF8", "membrane", "envelope",
]


# ── Lens.org API approach ───────────────────────────────────────────────

def lens_api_search(query_text: str) -> dict | None:
    """Search Lens.org patent API. Returns {total, sample_ids} or None."""
    if not LENS_TOKEN:
        return None

    payload = json.dumps({
        "query": {"match": {"full_text": query_text}},
        "size": 5,
        "include": ["lens_id", "title", "applicants.extracted_name.value",
                     "publication_type", "date_published"],
    }).encode()

    req = urllib.request.Request(LENS_API, data=payload, headers={
        "Authorization": f"Bearer {LENS_TOKEN}",
        "Content-Type": "application/json",
    })

    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read())
        total = data.get("total", 0)
        sample = []
        for r in data.get("data", [])[:5]:
            sample.append({
                "lens_id": r.get("lens_id"),
                "title": r.get("title", "")[:120],
                "applicants": [a.get("extracted_name", {}).get("value", "")
                               for a in r.get("applicants", [])[:2]],
            })
        return {"total": total, "sample": sample}
    except urllib.error.HTTPError as e:
        print(f"  Lens API error {e.code}: {e.read().decode()[:200]}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"  Lens API error: {e}", file=sys.stderr)
        return None


# ── Google Patents fallback ─────────────────────────────────────────────

def google_patents_count(query_terms: str) -> dict:
    """
    Fetch Google Patents search page and extract the result count.
    Returns {total, url, sample_patents}.
    """
    params = urllib.parse.urlencode({"q": query_terms, "type": "PATENT"})
    url = f"https://patents.google.com/?{params}"

    req = urllib.request.Request(url, headers={
        "User-Agent": ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                       "AppleWebKit/537.36 (KHTML, like Gecko) "
                       "Chrome/120.0.0.0 Safari/537.36"),
        "Accept": "text/html,application/xhtml+xml",
        "Accept-Language": "en-US,en;q=0.9",
    })

    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            html = resp.read().decode("utf-8", errors="replace")

        # Extract "About X,XXX results"
        m = re.search(r'About\s+([\d,]+)\s+results?', html)
        if m:
            total = int(m.group(1).replace(",", ""))
        else:
            # JS-rendered page — count won't be in static HTML
            total = None

        # Extract patent IDs from links
        patent_ids = re.findall(r'patents\.google\.com/patent/([A-Z]{2}\d+[A-Z0-9]*)', html)
        unique_ids = list(dict.fromkeys(patent_ids))[:5]

        return {"total": total, "url": url, "sample_patent_ids": unique_ids}
    except Exception as e:
        return {"total": None, "url": url, "error": str(e)}


# ── Lens.org public web search (no API key needed) ─────────────────────

def lens_web_count(query_text: str) -> dict:
    """
    Query Lens.org patent search (public web) and extract result count.
    Returns {total, url}.
    """
    params = urllib.parse.urlencode({
        "q": query_text,
        "p": 0,
        "n": 1,
    })
    url = f"https://www.lens.org/lens/search/patent/list?{params}"

    req = urllib.request.Request(url, headers={
        "User-Agent": ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                       "AppleWebKit/537.36"),
        "Accept": "text/html",
    })

    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            html = resp.read().decode("utf-8", errors="replace")

        # Lens uses "X results" or "X patent results" in the page
        m = re.search(r'([\d,]+)\s+(?:patent\s+)?results?', html, re.IGNORECASE)
        total = int(m.group(1).replace(",", "")) if m else None
        return {"total": total, "url": url}
    except Exception as e:
        return {"total": None, "url": url, "error": str(e)}


# ── Unified search function ─────────────────────────────────────────────

def search_patents(query_text: str) -> dict:
    """Try Lens API first, then Google Patents, then Lens web."""
    # Try Lens API (best structured data)
    if LENS_TOKEN:
        result = lens_api_search(query_text)
        if result is not None:
            return {**result, "source": "lens_api"}

    # Fall back to Google Patents
    result = google_patents_count(query_text)
    if result.get("total") is not None:
        return {**result, "source": "google_patents"}

    # Last resort: Lens web (no result count available from static HTML usually)
    result2 = lens_web_count(query_text)
    if result2.get("total") is not None:
        return {**result2, "source": "lens_web"}

    # Return whatever we got from Google Patents (may have sample IDs even without count)
    return {**result, "source": "google_patents_partial"}


# ── Main scan ────────────────────────────────────────────────────────────

def run_scan():
    report = {
        "scan_date": datetime.now(timezone.utc).isoformat(),
        "method": "lens_api" if LENS_TOKEN else "google_patents+lens_web",
        "landscape": {},
        "applications": {},
        "sars_cov2_genes": {},
        "guide_sequences": {},
    }

    # ─── 1. Landscape: CRISPR + pathogen ────────────────────────────────
    print("═══ Patent Landscape Scan ═══\n")
    print("Phase 1: CRISPR patent landscape per pathogen")
    print("─" * 55)

    for name, info in PATHOGENS.items():
        query = f"CRISPR AND ({info['query']})"
        print(f"  {name:20s} → ", end="", flush=True)
        result = search_patents(query)
        total = result.get("total", "?")
        print(f"{total:>8} patent families  [{result.get('source', '?')}]")
        report["landscape"][name] = {
            "query": query,
            **result,
        }
        time.sleep(DELAY)

    # ─── 2. Application breakdown for key pathogens ─────────────────────
    print("\nPhase 2: Application breakdown (diagnostic vs therapeutic)")
    print("─" * 55)

    for app_type in ["diagnostic", "therapeutic"]:
        report["applications"][app_type] = {}
        for name, info in PATHOGENS.items():
            query = f"CRISPR AND ({info['query']}) AND {app_type}"
            print(f"  {name:20s} [{app_type:12s}] → ", end="", flush=True)
            result = search_patents(query)
            total = result.get("total", "?")
            print(f"{total:>8}")
            report["applications"][app_type][name] = {
                "query": query,
                "total": result.get("total"),
                "source": result.get("source"),
            }
            time.sleep(DELAY)

    # ─── 3. Gene-level scan for SARS-CoV-2 ──────────────────────────────
    print("\nPhase 3: SARS-CoV-2 gene-level patent counts")
    print("─" * 55)

    for gene in SARS_GENES:
        query = f'CRISPR AND "SARS-CoV-2" AND "{gene}"'
        print(f"  {gene:20s} → ", end="", flush=True)
        result = search_patents(query)
        total = result.get("total", "?")
        print(f"{total:>8}")
        report["sars_cov2_genes"][gene] = {
            "query": query,
            "total": result.get("total"),
            "source": result.get("source"),
        }
        time.sleep(DELAY)

    # ─── 4. Guide sequence search (top 10 of our 52) ────────────────────
    print("\nPhase 4: Guide sequence search (top 10 targets)")
    print("─" * 55)

    if NOVELTY_FILE.exists():
        novelty = json.loads(NOVELTY_FILE.read_text())
        targets = novelty.get("targets", [])[:10]
        for i, t in enumerate(targets):
            seq = t.get("sequence", "")
            gene = t.get("gene", "?")
            if not seq or len(seq) < 20:
                continue
            query = f'"{seq}"'
            print(f"  #{i+1:2d} {seq[:20]:20s} ({gene[:15]}) → ", end="", flush=True)
            result = search_patents(query)
            total = result.get("total", "?")
            status = "CLEAR" if total == 0 or total is None else f"⚠ {total} match(es)"
            print(f"{status}")
            report["guide_sequences"][seq] = {
                "gene": gene,
                "total": result.get("total"),
                "source": result.get("source"),
                "status": "clear" if (total == 0 or total is None) else "review_needed",
            }
            time.sleep(DELAY)
    else:
        print("  ⚠ novelty_report.json not found, skipping sequence search")

    # ─── Write report ────────────────────────────────────────────────────
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_FILE.write_text(json.dumps(report, indent=2))
    print(f"\n✓ Report written to {OUTPUT_FILE.relative_to(ROOT)}")

    # ─── Summary ─────────────────────────────────────────────────────────
    print("\n═══ Summary ═══")
    for name, data in report["landscape"].items():
        total = data.get("total", "?")
        print(f"  {name:20s}: {total:>8} CRISPR patent families")

    return report


if __name__ == "__main__":
    run_scan()
