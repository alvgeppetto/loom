#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════╗
# ║  DEPRECATED — DO NOT USE                                    ║
# ║  Replaced by: `loom pubmed-scan` (Rust CLI, v4)             ║
# ║  This file is kept for reference only.                      ║
# ║  See: src/dna/pubmed_scan.rs                                ║
# ╚══════════════════════════════════════════════════════════════╝
"""
PubMed Literature Scanner for LOOM CRISPR Targets — v2.

Honest novelty assessment using gene-region and application-level queries,
not raw 23-mer sequences (which never appear in paper abstracts).

Strategy:
  1. Landscape scan — total CRISPR papers per pathogen
  2. Application scan — per pathogen, check CRISPR papers in:
     diagnostics, therapeutics, gene drive, antiviral
  3. Gene-region scan — per pathogen, check known gene targets
     (NS1, NS5, E protein, polymerase, etc.)
  4. Flag genuinely under-studied areas as research opportunities

Rate limit: 10 req/sec with API key, 3 req/sec without.
"""

import http.client
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API_KEY = os.environ.get("NCBI_API_KEY", "")
DELAY = 0.15 if API_KEY else 0.35

DATA_DIR = Path(__file__).resolve().parent.parent / "web" / "data"
INPUT_FILE = DATA_DIR / "crispr-targets.json"
OUTPUT_FILE = DATA_DIR / "pubmed-scan-results.json"
ONTOLOGY_FILE = DATA_DIR / "ontology-enrichment.json"

# Known gene regions and application areas per pathogen
PATHOGEN_GENES = {
    "dengue": {
        "search_name": "Dengue",
        "genes": ["NS1", "NS2A", "NS2B", "NS3", "NS4A", "NS4B", "NS5",
                  "capsid", "prM", "envelope"],
        "applications": ["diagnostic", "detection", "therapeutic", "antiviral",
                         "gene drive", "SHERLOCK", "DETECTR", "Cas13", "isothermal"],
    },
    "ebola": {
        "search_name": "Ebola",
        "genes": ["NP", "VP35", "VP40", "GP", "VP30", "VP24", "L protein"],
        "applications": ["diagnostic", "detection", "therapeutic",
                         "SHERLOCK", "DETECTR", "Cas13", "point-of-care"],
    },
    "hepatitis-b": {
        "search_name": "Hepatitis B",
        "genes": ["HBsAg", "HBcAg", "HBx", "polymerase", "precore",
                  "cccDNA", "pgRNA", "surface antigen", "core protein"],
        "applications": ["therapeutic", "antiviral", "gene editing", "knockdown",
                         "diagnostic", "detection", "cure", "cccDNA disruption"],
    },
    "mers": {
        "search_name": "MERS",
        "genes": ["spike", "RdRp", "ORF1a", "ORF1b",
                  "nucleocapsid", "envelope", "membrane"],
        "applications": ["diagnostic", "detection", "therapeutic",
                         "SHERLOCK", "DETECTR", "Cas13", "point-of-care"],
    },
    "rsv": {
        "search_name": "RSV",
        "genes": ["F protein", "G protein", "N protein", "L protein"],
        "applications": ["diagnostic", "detection", "therapeutic", "antiviral",
                         "SHERLOCK", "DETECTR", "Cas13"],
    },
    "zika": {
        "search_name": "Zika",
        "genes": ["NS1", "NS2A", "NS2B", "NS3", "NS4A", "NS4B", "NS5",
                  "capsid", "prM", "envelope"],
        "applications": ["diagnostic", "detection", "therapeutic",
                         "SHERLOCK", "DETECTR", "Cas13", "point-of-care",
                         "congenital", "microcephaly"],
    },
    "refseq-viral": {
        "search_name": "virus",
        "genes": ["polymerase", "capsid", "envelope", "protease",
                  "nucleocapsid", "spike"],
        "applications": ["diagnostic", "detection", "SHERLOCK", "DETECTR",
                         "multiplexed", "pan-viral"],
    },
    "hiv-1": {
        "search_name": "HIV",
        "genes": ["gag", "pol", "env", "tat", "rev", "nef", "vif", "vpr",
                  "vpu", "LTR", "integrase", "reverse transcriptase",
                  "protease", "gp120", "gp41", "capsid"],
        "applications": ["therapeutic", "cure", "latency", "reservoir",
                         "excision", "gene editing", "diagnostic", "detection",
                         "proviral", "SHERLOCK", "Cas13"],
    },
    "mpox": {
        "search_name": "Mpox OR Monkeypox",
        "genes": ["A33R", "B5R", "B6R", "E8L", "J2R", "A56R",
                  "thymidine kinase", "hemagglutinin", "envelope protein"],
        "applications": ["diagnostic", "detection", "therapeutic",
                         "SHERLOCK", "DETECTR", "Cas13", "point-of-care",
                         "isothermal"],
    },
    "tuberculosis": {
        "search_name": "Mycobacterium tuberculosis",
        "genes": ["rpoB", "katG", "inhA", "gyrA", "gyrB", "embB",
                  "pncA", "ethA", "rrs", "eis", "IS6110",
                  "Ag85B", "ESAT-6", "CFP-10"],
        "applications": ["diagnostic", "detection", "drug resistance",
                         "point-of-care", "SHERLOCK", "DETECTR", "Cas12",
                         "isothermal", "MDR-TB", "XDR-TB"],
    },
    "influenza-a": {
        "search_name": "Influenza A",
        "genes": ["HA", "NA", "hemagglutinin", "neuraminidase",
                  "PB1", "PB2", "PA", "NP", "M1", "M2", "NS1", "NEP"],
        "applications": ["diagnostic", "detection", "therapeutic", "antiviral",
                         "SHERLOCK", "DETECTR", "Cas13", "point-of-care",
                         "subtyping", "surveillance"],
    },
    "cholera": {
        "search_name": "Vibrio cholerae",
        "genes": ["ctxA", "ctxB", "tcpA", "toxR", "ompU", "ompT",
                  "hapA", "hlyA", "rtxA", "cholera toxin"],
        "applications": ["diagnostic", "detection", "point-of-care",
                         "SHERLOCK", "DETECTR", "Cas12", "isothermal",
                         "environmental monitoring", "serogroup typing"],
    },
    "sars-cov-2": {
        "search_name": "SARS-CoV-2",
        "genes": ["spike", "S protein", "RBD", "RdRp", "ORF1a", "ORF1b",
                  "nucleocapsid", "N protein", "envelope", "membrane",
                  "ORF3a", "ORF7a", "ORF8", "nsp12", "nsp14"],
        "applications": ["diagnostic", "detection", "therapeutic", "antiviral",
                         "SHERLOCK", "DETECTR", "Cas13", "Cas12",
                         "point-of-care", "variant detection", "isothermal"],
    },
}

# Gene synonyms: maps gene name → list of alternative names for expanded queries.
# If a gene is not here, only the primary name is searched.
GENE_SYNONYMS = {
    # TB genes
    "inhA":   ["enoyl-ACP reductase", "isoniazid resistance", "Rv1484"],
    "gyrB":   ["gyrase B", "fluoroquinolone resistance"],
    "pncA":   ["pyrazinamidase", "pyrazinamide resistance", "Rv2043c"],
    "rrs":    ["16S rRNA", "aminoglycoside resistance"],
    "eis":    ["enhanced intracellular survival", "kanamycin resistance", "Rv2416c"],
    "CFP-10": ["Rv3874", "ESAT-6/CFP-10", "culture filtrate protein 10"],
    "ethA":   ["ethionamide resistance", "Rv3854c"],
    # Cholera genes
    "ctxA":   ["cholera toxin A", "CT-A", "cholera toxin subunit A"],
    "ctxB":   ["cholera toxin B", "CT-B", "cholera toxin subunit B"],
    # Ebola genes
    "NP":     ["nucleoprotein"],
    "GP":     ["glycoprotein"],
    "L protein": ["RNA-dependent RNA polymerase"],
    # MERS genes
    "spike":  ["S protein", "spike protein"],
    "ORF1a":  ["nsp1", "nsp2", "nsp3", "papain-like protease"],
    "ORF1b":  ["nsp12", "nsp13", "nsp14", "nsp15", "nsp16"],
    # RSV genes
    "F protein":  ["fusion protein"],
    "G protein":  ["attachment protein"],
    "N protein":  ["nucleoprotein"],
}


def _load_ontology_synonyms() -> dict[str, list[str]]:
    """Auto-generate gene synonyms from NCBI ontology annotations.

    For each gene in PATHOGEN_GENES, match it to ontology data by:
      1. Direct match: gene name == ontology symbol (case-insensitive)
      2. Protein suffix match: "X protein" → symbol "X"
      3. Common name match: gene name appears in ontology gene name

    Returns dict of gene → [additional synonyms] to merge with GENE_SYNONYMS.
    Does NOT overwrite existing manual synonyms — only adds new ones.
    """
    if not ONTOLOGY_FILE.exists():
        return {}

    with open(ONTOLOGY_FILE) as f:
        data = json.load(f)

    ontology_pathogens = data.get("pathogens", {})
    extra_synonyms: dict[str, set[str]] = {}

    for pathogen_key, cfg in PATHOGEN_GENES.items():
        ont_data = ontology_pathogens.get(pathogen_key)
        if not ont_data:
            continue

        # Extract gene annotations from ontology
        genes_struct = ont_data.get("genes", {})
        if isinstance(genes_struct, dict):
            gene_list = genes_struct.get("genes", [])
        elif isinstance(genes_struct, list):
            gene_list = genes_struct
        else:
            continue

        # Build lookup: lowercase symbol → {symbol, name}
        ont_by_symbol = {}
        ont_by_name_lower = {}
        for g in gene_list:
            if not isinstance(g, dict):
                continue
            sym = g.get("symbol", "").strip()
            name = g.get("name", "").strip()
            if sym:
                ont_by_symbol[sym.lower()] = {"symbol": sym, "name": name}
            if name:
                ont_by_name_lower[name.lower()] = {"symbol": sym, "name": name}

        # Match each of our scanner genes to ontology
        for gene in cfg.get("genes", []):
            gene_lower = gene.lower()
            matched = None

            # Strategy 1: Direct symbol match (e.g., "NP" → "NP")
            if gene_lower in ont_by_symbol:
                matched = ont_by_symbol[gene_lower]

            # Strategy 2: Strip " protein" suffix (e.g., "F protein" → "F")
            if not matched and gene_lower.endswith(" protein"):
                stem = gene_lower[:-len(" protein")]
                if stem in ont_by_symbol:
                    matched = ont_by_symbol[stem]

            # Strategy 3: Check if gene name IS the ontology name
            if not matched and gene_lower in ont_by_name_lower:
                matched = ont_by_name_lower[gene_lower]

            if matched:
                sym = matched["symbol"]
                name = matched["name"]
                additions = set()
                # Add the NCBI gene name as synonym (it's the official annotation)
                if name and name.lower() != gene_lower:
                    additions.add(name)
                # Add the NCBI symbol if different from our gene name
                if sym and sym.lower() != gene_lower:
                    additions.add(sym)

                if additions:
                    if gene not in extra_synonyms:
                        extra_synonyms[gene] = set()
                    extra_synonyms[gene].update(additions)

    return {k: sorted(v) for k, v in extra_synonyms.items()}


# Merge ontology-derived synonyms into GENE_SYNONYMS (manual entries take priority)
_ontology_synonyms = _load_ontology_synonyms()
for _gene, _aliases in _ontology_synonyms.items():
    existing = set(GENE_SYNONYMS.get(_gene, []))
    new_aliases = [a for a in _aliases if a not in existing
                   and a.lower() != _gene.lower()]
    if new_aliases:
        GENE_SYNONYMS.setdefault(_gene, []).extend(new_aliases)


def eutils_get(endpoint: str, params: dict) -> dict:
    """Make a GET request to NCBI E-utilities."""
    params["retmode"] = "json"
    if API_KEY:
        params["api_key"] = API_KEY
    url = f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(params)}"
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=20) as resp:
                return json.loads(resp.read())
        except (urllib.error.URLError, TimeoutError, http.client.RemoteDisconnected,
                ConnectionResetError, OSError) as exc:
            wait = 2 * (attempt + 1)
            if attempt < 3:
                print(f"\n  [retry {attempt+1}/3] {exc.__class__.__name__}, "
                      f"waiting {wait}s...", end="", flush=True)
                time.sleep(wait)
                continue
            raise RuntimeError(f"E-utilities failed: {exc}") from exc
    return {}


def esearch(query: str, max_ids: int = 0) -> tuple[int, list[str]]:
    """Return (count, [pmid...]) for a PubMed query."""
    data = eutils_get("esearch.fcgi",
                      {"db": "pubmed", "term": query, "retmax": max_ids})
    result = data.get("esearchresult", {})
    return int(result.get("count", 0)), result.get("idlist", [])


def esummary(pmids: list[str]) -> list[dict]:
    """Fetch paper summaries."""
    if not pmids:
        return []
    data = eutils_get("esummary.fcgi",
                      {"db": "pubmed", "id": ",".join(pmids)})
    results = data.get("result", {})
    papers = []
    for pmid in pmids:
        r = results.get(pmid)
        if not r or "error" in r:
            continue
        papers.append({
            "pmid": pmid,
            "title": r.get("title", ""),
            "authors": ", ".join(a["name"] for a in r.get("authors", [])[:4]),
            "journal": r.get("source", ""),
            "year": r.get("pubdate", "").split(" ")[0],
        })
    return papers


def _name_clause(name: str) -> str:
    """Build the pathogen name clause, handling OR-separated synonyms."""
    if " OR " in name:
        parts = [p.strip() for p in name.split(" OR ")]
        clauses = " OR ".join(f'"{p}"[Title/Abstract]' for p in parts)
        return f"({clauses})"
    return f'("{name}"[Title/Abstract])'


def landscape_scan(pathogens: dict) -> dict:
    """Total CRISPR papers per pathogen."""
    print("\n=== Landscape Scan: Total CRISPR Papers ===")
    results = {}
    for key in pathogens:
        cfg = PATHOGEN_GENES.get(key, {})
        name = cfg.get("search_name", pathogens[key]["name"])
        name_part = _name_clause(name)
        query = (
            f'{name_part} AND '
            f'(CRISPR[Title/Abstract] OR "guide RNA"[Title/Abstract] '
            f'OR Cas9[Title/Abstract] OR Cas12[Title/Abstract] '
            f'OR Cas13[Title/Abstract])'
        )
        count, _ = esearch(query)
        results[key] = {"query": query, "total_papers": count}
        print(f"  {pathogens[key]['name']:30s} → {count:>5} papers")
        time.sleep(DELAY)
    return results


def application_scan(pathogens: dict) -> dict:
    """Per pathogen: papers per application area."""
    print("\n=== Application Scan: CRISPR by Use Case ===")
    results = {}

    for key in pathogens:
        cfg = PATHOGEN_GENES.get(key, {})
        name = cfg.get("search_name", pathogens[key]["name"])
        apps = cfg.get("applications",
                       ["diagnostic", "therapeutic", "detection"])

        app_results = {}
        name_part = _name_clause(name)
        for app in apps:
            query = (f'{name_part} AND CRISPR[Title/Abstract] '
                     f'AND ("{app}"[Title/Abstract])')
            count, ids = esearch(query, 3)
            time.sleep(DELAY)

            papers = []
            if ids:
                papers = esummary(ids)
                time.sleep(DELAY)

            app_results[app] = {"count": count, "papers": papers}

        covered = [a for a, r in app_results.items() if r["count"] > 0]
        gaps = [a for a, r in app_results.items() if r["count"] == 0]

        results[key] = {
            "pathogen": pathogens[key]["name"],
            "applications": app_results,
            "covered_areas": covered,
            "gap_areas": gaps,
        }

        print(f"\n  {pathogens[key]['name']}:")
        for app, r in sorted(app_results.items(), key=lambda x: -x[1]["count"]):
            marker = "●" if r["count"] > 0 else "○"
            print(f"    {marker} {app:25s} {r['count']:>4} papers")
        if gaps:
            print(f"    → GAPS: {', '.join(gaps)}")

    return results


def gene_scan(pathogens: dict) -> dict:
    """Per pathogen: CRISPR papers per gene target."""
    print("\n=== Gene Region Scan: CRISPR by Target Gene ===")
    results = {}

    for key in pathogens:
        cfg = PATHOGEN_GENES.get(key, {})
        name = cfg.get("search_name", pathogens[key]["name"])
        genes = cfg.get("genes", [])

        gene_results = {}
        name_part = _name_clause(name)
        for gene in genes:
            # Build expanded gene clause with synonyms
            all_names = [gene] + GENE_SYNONYMS.get(gene, [])
            gene_clause = " OR ".join(
                f'"{g}"[Title/Abstract]' for g in all_names
            )
            query = (f'{name_part} AND CRISPR[Title/Abstract] '
                     f'AND ({gene_clause})')
            count, ids = esearch(query, 3)
            time.sleep(DELAY)

            papers = []
            if ids:
                papers = esummary(ids)
                time.sleep(DELAY)

            gene_results[gene] = {"count": count, "papers": papers}

        studied = [g for g, r in gene_results.items() if r["count"] > 0]
        unstudied = [g for g, r in gene_results.items() if r["count"] == 0]

        results[key] = {
            "pathogen": pathogens[key]["name"],
            "genes": gene_results,
            "studied_genes": studied,
            "unstudied_genes": unstudied,
        }

        print(f"\n  {pathogens[key]['name']}:")
        for g, r in sorted(gene_results.items(), key=lambda x: -x[1]["count"]):
            marker = "●" if r["count"] > 0 else "○"
            print(f"    {marker} {g:25s} {r['count']:>4} papers")
        if unstudied:
            print(f"    → UNSTUDIED: {', '.join(unstudied)}")

    return results


def generate_report(landscape, app_results, gene_results, pathogens):
    """Compose the full report with opportunities."""
    opportunities = []
    for key in pathogens:
        pname = pathogens[key]["name"]
        for gap in app_results.get(key, {}).get("gap_areas", []):
            opportunities.append({
                "type": "application_gap",
                "pathogen": pname,
                "pathogen_key": key,
                "area": gap,
                "description": f"No CRISPR + {gap} papers for {pname}",
            })
        for gene in gene_results.get(key, {}).get("unstudied_genes", []):
            opportunities.append({
                "type": "gene_gap",
                "pathogen": pname,
                "pathogen_key": key,
                "gene": gene,
                "description": f"No CRISPR papers targeting {gene} in {pname}",
            })

    total_targets = sum(p["target_count"] for p in pathogens.values())
    ngg_targets = sum(
        sum(1 for t in p["targets"] if t["m"] == "NGG")
        for p in pathogens.values()
    )

    # Per-pathogen summary stats
    pathogen_summaries = {}
    for key in pathogens:
        land = landscape.get(key, {})
        apps = app_results.get(key, {})
        genes = gene_results.get(key, {})
        pathogen_summaries[key] = {
            "name": pathogens[key]["name"],
            "total_papers": land.get("total_papers", 0),
            "targets_in_db": pathogens[key]["target_count"],
            "covered_applications": apps.get("covered_areas", []),
            "gap_applications": apps.get("gap_areas", []),
            "studied_genes": genes.get("studied_genes", []),
            "unstudied_genes": genes.get("unstudied_genes", []),
        }

    return {
        "generated": datetime.now(timezone.utc).isoformat(),
        "version": 2,
        "methodology": (
            "Region-level PubMed scan with ontology-enhanced synonym expansion. "
            "Gene synonyms are auto-generated from NCBI gene annotations via "
            "ontology-enrichment.json (mapping official gene symbols to common "
            "research names) and merged with curated manual synonyms. "
            "Searches pathogen + application area and pathogen + gene region "
            "combinations. A 'gap' means zero PubMed results for that "
            "combination — indicating an under-studied area where our "
            "pre-computed CRISPR targets could enable new research."
        ),
        "summary": {
            "total_pathogens": len(pathogens),
            "total_targets_in_db": total_targets,
            "ngg_targets": ngg_targets,
            "application_gaps": len([o for o in opportunities
                                     if o["type"] == "application_gap"]),
            "gene_gaps": len([o for o in opportunities
                              if o["type"] == "gene_gap"]),
            "total_opportunities": len(opportunities),
        },
        "landscape": landscape,
        "applications": {k: v["applications"] for k, v in app_results.items()},
        "genes": {k: v["genes"] for k, v in gene_results.items()},
        "pathogens": pathogen_summaries,
        "opportunities": opportunities,
    }


def print_exec_summary(report, pathogens):
    """Print scan results summary."""
    s = report["summary"]
    print("\n" + "=" * 70)
    print("  PUBMED SCAN v2 — RESULTS SUMMARY")
    print("=" * 70)
    print(f"  Database:            {s['total_targets_in_db']:,} CRISPR targets "
          f"({s['ngg_targets']:,} NGG/SpCas9)")
    print(f"  Pathogens:           {s['total_pathogens']}")
    print(f"  Application gaps:    {s['application_gaps']}")
    print(f"  Gene region gaps:    {s['gene_gaps']}")
    print(f"  Total opportunities: {s['total_opportunities']}")
    print("-" * 70)

    print("\n  PER-PATHOGEN SUMMARY:")
    for key, ps in report["pathogens"].items():
        n_gaps = len(ps["gap_applications"]) + len(ps["unstudied_genes"])
        print(f"    {ps['name']:25s} | {ps['total_papers']:>5} papers | "
              f"{ps['targets_in_db']:>5} targets | {n_gaps} gaps")

    opps = report["opportunities"]
    if opps:
        print(f"\n  TOP RESEARCH OPPORTUNITIES ({len(opps)} total):")
        for i, opp in enumerate(opps[:25], 1):
            print(f"    {i:2d}. {opp['description']}")

    print("=" * 70)


def main():
    if not INPUT_FILE.exists():
        print(f"Error: {INPUT_FILE} not found", file=sys.stderr)
        sys.exit(1)

    with open(INPUT_FILE) as f:
        pathogens = json.load(f)

    # Exclude non-pathogen entries (e.g., Human GRCh38 is for off-target filtering)
    EXCLUDE_KEYS = {"human-grch38", "human"}
    pathogens = {k: v for k, v in pathogens.items() if k not in EXCLUDE_KEYS}

    print(f"Loaded {len(pathogens)} pathogens, "
          f"{sum(p['target_count'] for p in pathogens.values()):,} total targets")
    print(f"API key: {'configured (10 req/sec)' if API_KEY else 'not set (3 req/sec)'}")

    total_queries = (
        len(pathogens)
        + sum(len(PATHOGEN_GENES.get(k, {}).get("applications", [])) * 2
              for k in pathogens)
        + sum(len(PATHOGEN_GENES.get(k, {}).get("genes", [])) * 2
              for k in pathogens)
    )
    print(f"Estimated: ~{total_queries} API calls ({total_queries * DELAY / 60:.1f} min)")

    landscape = landscape_scan(pathogens)
    app_results = application_scan(pathogens)
    gene_results = gene_scan(pathogens)
    report = generate_report(landscape, app_results, gene_results, pathogens)

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_FILE, "w") as f:
        json.dump(report, f, indent=2)
    print(f"\nResults saved to {OUTPUT_FILE}")

    print_exec_summary(report, pathogens)


if __name__ == "__main__":
    main()
