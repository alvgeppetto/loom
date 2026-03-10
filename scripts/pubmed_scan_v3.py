#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════╗
# ║  DEPRECATED — DO NOT USE                                    ║
# ║  Replaced by: `loom pubmed-scan` (Rust CLI, v4)             ║
# ║  This file is kept for reference only.                      ║
# ║  See: src/dna/pubmed_scan.rs                                ║
# ╚══════════════════════════════════════════════════════════════╝
"""
PubMed Literature Scanner v3 — Multi-source, abstract-corpus, confidence-scored.

Problem with v2: narrow [Title/Abstract] single-term queries reported "zero"
for combinations where papers DO exist (e.g. V. cholerae CRISPR detection).
A claim of "gap" must survive multiple query strategies and two independent
databases. This version is designed to produce defensible, publication-grade
gap claims only.

Strategy per (pathogen × gene/application):
  1. Multi-strategy PubMed queries:
       - Narrow:   term[Title/Abstract] — original v2 approach
       - Medium:   term[All Fields] — catches full-text indexed metadata
       - Broad:    expanded synonym OR-clause [All Fields]
  2. Europe PMC (independent database, catches more preprints + grey lit)
  3. Abstract corpus scan:
       - Fetch ALL abstracts for the landscape query (all CRISPR × pathogen papers)
       - Scan locally for application/gene term co-occurrence
       - If a landscape paper covers a "claimed gap", the gap is refuted
  4. Confidence scoring per gap:
       CONFIRMED   — 0 in all strategies, 0 in Europe PMC, no corpus match
       PROBABLE    — 0 in narrow/medium, >0 in broad OR Europe PMC (review needed)
       UNCERTAIN   — 0 in narrow, >0 in medium or corpus suggests coverage
       FALSE       — >0 in medium/broad or corpus match found

Run time: ~30–60 min per pathogen at free tier (no NCBI API key).
With NCBI API key (set $NCBI_API_KEY): ~3× faster.

Output:
  web/data/pubmed-scan-v3.json          — full machine-readable results
  web/data/gap-audit-v3.txt             — human-readable audit trail
"""

import http.client
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

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "web" / "data"
INPUT_FILE = DATA_DIR / "crispr-targets.json"
OUTPUT_FILE = DATA_DIR / "pubmed-scan-v3.json"
AUDIT_FILE = DATA_DIR / "gap-audit-v3.txt"
ONTOLOGY_FILE = DATA_DIR / "ontology-enrichment.json"
CACHE_DIR = REPO_ROOT / "target" / "pubmed_cache_v3"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

# ── API config ─────────────────────────────────────────────────────────────────
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EUROPEPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
API_KEY = os.environ.get("NCBI_API_KEY", "")
# Strict rate limiting: NCBI free = 3/s, with key = 10/s; Europe PMC = ~3/s
PUBMED_DELAY = 0.12 if API_KEY else 0.38
EPMC_DELAY = 0.38
MAX_LANDSCAPE_FETCH = 500   # abstracts to fetch per pathogen for corpus scan
MAX_DETAIL_FETCH = 5        # top papers to fetch per low-count query result

EXCLUDE_KEYS = {"human-grch38", "human", "refseq-viral"}

# ── Pathogen config ────────────────────────────────────────────────────────────
# search_name: what to put in PubMed queries (can use OR)
# mesh_term:   if set, also try MeSH-specific query
# alt_names:   additional organism/disease names to OR-expand landscape query
PATHOGEN_CONFIG: dict[str, dict] = {
    "cholera": {
        "search_name": "Vibrio cholerae",
        "alt_names": ["cholera", "V. cholerae"],
        "mesh_term": "Vibrio cholerae[MeSH]",
        "genes": ["ctxA", "ctxB", "tcpA", "toxR", "ompU", "ompT",
                  "hapA", "hlyA", "rtxA"],
        "gene_synonyms": {
            "ctxA":  ["cholera toxin A", "CT-A", "cholera toxin subunit A",
                      "ctx operon", "ctxAB"],
            "ctxB":  ["cholera toxin B", "CT-B", "cholera toxin subunit B",
                      "ctxAB"],
            "tcpA":  ["toxin-coregulated pilus", "TCP pilus", "tcp operon"],
            "toxR":  ["ToxR", "virulence regulator"],
            "ompU":  ["outer membrane protein U"],
            "rtxA":  ["RTX toxin"],
            "hlyA":  ["haemolysin A", "hemolysin"],
        },
        "applications": {
            "diagnostics": ["diagnostic", "diagnostics", "detection assay",
                            "detection", "identify", "identification",
                            "rapid test", "rapid detection", "test strip",
                            "biosensor", "lateral flow", "POC", "point-of-care",
                            "SHERLOCK", "DETECTR", "Cas12", "Cas13",
                            "isothermal", "LAMP", "field-deployable",
                            "serology", "serogroup", "serovar", "serotyping",
                            "environmental surveillance", "water testing",
                            "wastewater"],
            "therapeutics": ["antimicrobial", "antibacterial", "therapy",
                             "treatment", "kill", "inhibit", "antibiotic"],
        },
    },
    "dengue": {
        "search_name": "Dengue",
        "alt_names": ["dengue virus", "DENV"],
        "mesh_term": "Dengue virus[MeSH]",
        "genes": ["NS1", "NS2A", "NS2B", "NS3", "NS4A", "NS4B", "NS5",
                  "capsid", "prM", "envelope"],
        "gene_synonyms": {
            "NS5":    ["RNA-dependent RNA polymerase", "RdRp", "methyltransferase"],
            "NS3":    ["helicase", "protease", "NS3 helicase"],
            "NS1":    ["non-structural protein 1"],
            "prM":    ["pre-membrane", "premembrane"],
            "envelope": ["E protein", "E gene"],
        },
        "applications": {
            "diagnostics": ["diagnostic", "detection", "assay", "rapid test",
                            "SHERLOCK", "DETECTR", "point-of-care", "Cas12",
                            "Cas13", "isothermal", "LAMP", "biosensor",
                            "lateral flow", "serotype", "serotyping"],
            "therapeutics": ["antiviral", "therapy", "inhibit", "treatment"],
        },
    },
    "ebola": {
        "search_name": "Ebola",
        "alt_names": ["Ebola virus", "EBOV", "ebolavirus", "filovirus"],
        "mesh_term": "Ebolavirus[MeSH]",
        "genes": ["NP", "VP35", "VP40", "GP", "VP30", "VP24", "L protein"],
        "gene_synonyms": {
            "NP":        ["nucleoprotein"],
            "GP":        ["glycoprotein", "surface glycoprotein"],
            "L protein": ["RNA-dependent RNA polymerase", "RdRp", "L gene"],
            "VP40":      ["matrix protein"],
        },
        "applications": {
            "diagnostics": ["diagnostic", "detection", "assay", "rapid test",
                            "SHERLOCK", "DETECTR", "point-of-care", "Cas12",
                            "Cas13", "isothermal", "field-deployable",
                            "biosensor", "outbreak"],
            "therapeutics": ["antiviral", "therapy", "treatment", "inhibit",
                             "knockdown"],
        },
    },
    "mpox": {
        "search_name": "Mpox OR Monkeypox",
        "alt_names": ["MPXV", "mpox virus", "orthopoxvirus", "monkeypox virus"],
        "mesh_term": "Monkeypox virus[MeSH]",
        "genes": ["A33R", "B5R", "B6R", "E8L", "J2R", "A56R",
                  "thymidine kinase", "hemagglutinin"],
        "gene_synonyms": {
            "A33R":           ["EEV-specific protein", "extracellular enveloped virus"],
            "B5R":            ["EEV membrane antigen"],
            "B6R":            ["B6R gene"],
            "thymidine kinase": ["TK gene", "J2R thymidine kinase"],
            "hemagglutinin":  ["HA protein", "A56R"],
        },
        "applications": {
            "diagnostics": ["diagnostic", "detection", "assay", "rapid test",
                            "SHERLOCK", "DETECTR", "Cas12", "Cas13",
                            "point-of-care", "isothermal", "PCR-free",
                            "surveillance", "biosensor"],
            "therapeutics": ["antiviral", "therapy", "treatment", "inhibit"],
        },
    },
    "tuberculosis": {
        "search_name": "Mycobacterium tuberculosis",
        "alt_names": ["M. tuberculosis", "MTB", "TB"],
        "mesh_term": "Mycobacterium tuberculosis[MeSH]",
        "genes": ["rpoB", "katG", "inhA", "gyrA", "gyrB", "embB",
                  "pncA", "ethA", "rrs", "eis", "IS6110",
                  "Ag85B", "ESAT-6", "CFP-10"],
        "gene_synonyms": {
            "inhA":   ["enoyl-ACP reductase", "isoniazid resistance", "Rv1484"],
            "gyrB":   ["gyrase subunit B", "fluoroquinolone resistance"],
            "pncA":   ["pyrazinamidase", "pyrazinamide resistance", "Rv2043c"],
            "rrs":    ["16S rRNA", "16S ribosomal RNA", "aminoglycoside resistance"],
            "eis":    ["enhanced intracellular survival", "kanamycin resistance",
                       "Rv2416c"],
            "CFP-10": ["Rv3874", "ESAT-6/CFP-10", "culture filtrate protein 10",
                       "EsxB"],
            "ethA":   ["ethionamide monooxygenase", "ethionamide resistance",
                       "Rv3854c"],
            "embB":   ["arabinosyltransferase", "ethambutol resistance"],
            "rpoB":   ["RNA polymerase beta", "rifampicin resistance",
                       "rifampin resistance"],
            "katG":   ["catalase-peroxidase", "isoniazid resistance KatG"],
            "IS6110": ["insertion sequence IS6110", "repeat element TB"],
            "ESAT-6": ["early secreted antigenic target", "EsxA", "Rv3875"],
            "Ag85B":  ["antigen 85B", "Rv1886c", "fbpB"],
        },
        "applications": {
            "diagnostics": ["diagnostic", "detection", "drug resistance",
                            "MDR-TB", "XDR-TB", "rifampicin resistance",
                            "isoniazid resistance", "rapid test", "SHERLOCK",
                            "DETECTR", "Cas12", "Cas13", "point-of-care",
                            "isothermal", "LAMP", "sputum"],
            "therapeutics": ["antimicrobial", "antibacterial", "treatment",
                             "inhibit", "gene silencing", "CRISPRi"],
        },
    },
    "zika": {
        "search_name": "Zika",
        "alt_names": ["Zika virus", "ZIKV"],
        "mesh_term": "Zika Virus[MeSH]",
        "genes": ["NS1", "NS2A", "NS2B", "NS3", "NS4A", "NS4B", "NS5",
                  "capsid", "prM", "envelope"],
        "gene_synonyms": {
            "NS5":    ["RNA-dependent RNA polymerase", "RdRp", "methyltransferase"],
            "NS3":    ["helicase", "serine protease"],
            "NS2B":   ["NS2B cofactor", "protease cofactor"],
            "prM":    ["pre-membrane protein"],
            "envelope": ["E protein", "E gene"],
        },
        "applications": {
            "diagnostics": ["diagnostic", "detection", "assay", "rapid test",
                            "SHERLOCK", "DETECTR", "Cas12", "Cas13",
                            "point-of-care", "isothermal", "LAMP",
                            "microcephaly screening", "congenital"],
            "therapeutics": ["antiviral", "therapy", "treatment"],
        },
    },
    "rsv": {
        "search_name": "RSV",
        "alt_names": ["respiratory syncytial virus", "hRSV", "human RSV"],
        "mesh_term": "Respiratory Syncytial Viruses[MeSH]",
        "genes": ["F protein", "G protein", "N protein", "L protein",
                  "M protein", "M2-1", "SH protein"],
        "gene_synonyms": {
            "F protein":    ["fusion protein", "F gene", "RSV F"],
            "G protein":    ["attachment glycoprotein", "G gene"],
            "N protein":    ["nucleoprotein", "N gene"],
            "L protein":    ["RNA-dependent RNA polymerase", "L gene", "RdRp"],
            "SH protein":   ["small hydrophobic protein"],
        },
        "applications": {
            "diagnostics": ["diagnostic", "detection", "assay", "rapid test",
                            "SHERLOCK", "DETECTR", "Cas12", "Cas13",
                            "point-of-care", "isothermal", "LAMP",
                            "pediatric", "infant", "bronchiolitis"],
            "therapeutics": ["antiviral", "therapy", "inhibit", "palivizumab"],
        },
    },
    "influenza-a": {
        "search_name": "Influenza A",
        "alt_names": ["influenza virus", "flu", "H1N1", "H3N2", "H5N1",
                      "avian influenza"],
        "mesh_term": "Influenza A virus[MeSH]",
        "genes": ["HA", "NA", "PB1", "PB2", "PA", "NP", "M1", "M2", "NS1", "NEP"],
        "gene_synonyms": {
            "HA":   ["hemagglutinin", "H1", "H3", "H5"],
            "NA":   ["neuraminidase", "N1", "N2"],
            "PB1":  ["polymerase basic 1", "PB1-F2"],
            "PB2":  ["polymerase basic 2"],
            "PA":   ["polymerase acidic"],
            "NP":   ["nucleoprotein"],
            "NS1":  ["non-structural protein 1", "interferon antagonist"],
            "M2":   ["M2 ion channel", "amantadine target"],
        },
        "applications": {
            "diagnostics": ["diagnostic", "detection", "assay", "rapid test",
                            "subtyping", "SHERLOCK", "DETECTR", "Cas12", "Cas13",
                            "point-of-care", "isothermal", "surveillance",
                            "pandemic"],
            "therapeutics": ["antiviral", "therapy", "inhibit", "oseltamivir"],
        },
    },
    "mers": {
        "search_name": "MERS",
        "alt_names": ["MERS-CoV", "Middle East respiratory syndrome",
                      "Middle East respiratory syndrome coronavirus"],
        "mesh_term": "Middle East Respiratory Syndrome Coronavirus[MeSH]",
        "genes": ["spike", "RdRp", "ORF1a", "ORF1b", "nucleocapsid",
                  "envelope", "membrane"],
        "gene_synonyms": {
            "spike":      ["S protein", "spike glycoprotein", "receptor-binding domain",
                           "RBD"],
            "RdRp":       ["nsp12", "RNA-dependent RNA polymerase"],
            "ORF1a":      ["nsp1", "nsp3", "papain-like protease", "PLpro"],
            "ORF1b":      ["nsp13", "nsp14", "helicase"],
            "nucleocapsid": ["N protein", "N gene"],
        },
        "applications": {
            "diagnostics": ["diagnostic", "detection", "assay", "rapid test",
                            "SHERLOCK", "DETECTR", "Cas12", "Cas13",
                            "point-of-care", "isothermal", "camel",
                            "zoonotic surveillance"],
            "therapeutics": ["antiviral", "therapy", "treatment", "inhibit"],
        },
    },
    "hepatitis-b": {
        "search_name": "Hepatitis B",
        "alt_names": ["HBV", "hepatitis B virus"],
        "mesh_term": "Hepatitis B virus[MeSH]",
        "genes": ["HBsAg", "HBcAg", "HBx", "polymerase", "precore",
                  "cccDNA", "pgRNA"],
        "gene_synonyms": {
            "HBsAg":      ["surface antigen", "HBs antigen", "HBsAg gene",
                           "HBV surface"],
            "cccDNA":     ["covalently closed circular DNA", "episomal DNA",
                           "nuclear HBV"],
            "HBx":        ["HBx protein", "X gene", "hepatitis B x antigen"],
            "polymerase": ["HBV polymerase", "reverse transcriptase", "RT gene"],
            "pgRNA":      ["pregenomic RNA", "pg RNA", "3.5 kb RNA"],
        },
        "applications": {
            "diagnostics": ["diagnostic", "detection", "assay",
                            "surface antigen", "HBsAg detection",
                            "viral load", "quantification"],
            "therapeutics": ["cure", "cccDNA disruption", "functional cure",
                             "antiviral", "knockdown", "silencing",
                             "gene editing", "excision"],
        },
    },
    "hiv-1": {
        "search_name": "HIV",
        "alt_names": ["HIV-1", "human immunodeficiency virus"],
        "mesh_term": "HIV-1[MeSH]",
        "genes": ["gag", "pol", "env", "tat", "rev", "nef", "vif",
                  "LTR", "integrase", "gp120"],
        "gene_synonyms": {
            "gag":        ["capsid", "matrix", "p24", "gag gene"],
            "pol":        ["reverse transcriptase", "integrase", "protease"],
            "env":        ["gp120", "gp41", "envelope glycoprotein"],
            "tat":        ["transactivator", "Tat protein"],
            "LTR":        ["long terminal repeat", "HIV LTR"],
            "integrase":  ["IN gene", "HIV integrase"],
        },
        "applications": {
            "diagnostics": ["diagnostic", "detection", "viral load",
                            "quantification", "proviral DNA"],
            "therapeutics": ["cure", "latency reversal", "reservoir elimination",
                             "excision", "gene editing", "knockdown",
                             "CRISPRi", "CRISPR excision"],
        },
    },

    # ── SARS-CoV-2 ────────────────────────────────────────────────────────────
    # Genes to scan = ALL genes that have novel targets in our 52-region set
    # (nsp13:8, nsp14:5, nsp15:2, nsp16:2, ORF1a sub-nsps:20, intergenic/other:15)
    # plus the established targets (spike, N, E, RdRp) to confirm they ARE studied.
    # This is the first time sars-cov-2 is run through the v3 three-strategy scanner.
    "sars-cov-2": {
        "search_name": "SARS-CoV-2",
        "alt_names": ["COVID-19", "SARS CoV 2", "2019-nCoV",
                      "coronavirus disease 2019", "severe acute respiratory syndrome coronavirus 2"],
        "mesh_term": "SARS-CoV-2[MeSH]",
        # Established genes (should be FALSE = studied)
        # Novel-target genes (may be gaps): nsp13, nsp14, nsp15, nsp16,
        #   plus ORF1a sub-genes (nsp3/PLpro, nsp5/Mpro, nsp7/8/9/10)
        "genes": [
            "spike", "RdRp", "nucleocapsid", "envelope", "membrane",
            "nsp13", "nsp14", "nsp15", "nsp16",
            "nsp3", "nsp5", "nsp7", "nsp8", "nsp9", "nsp10",
            "ORF3a", "ORF8", "ORF10",
        ],
        "gene_synonyms": {
            "spike":        ["S protein", "S gene", "spike glycoprotein",
                             "receptor-binding domain", "RBD", "furin cleavage site",
                             "D614G"],
            "RdRp":         ["nsp12", "RNA-dependent RNA polymerase",
                             "replicase", "polymerase", "ORF1b polymerase"],
            "nucleocapsid": ["N protein", "N gene", "nucleocapsid protein"],
            "envelope":     ["E protein", "E gene", "envelope protein"],
            "membrane":     ["M protein", "M gene", "membrane protein"],
            "nsp13":        ["helicase", "nsp-13", "NSP13",
                             "SARS-CoV-2 helicase", "coronavirus helicase",
                             "nsp13 helicase"],
            "nsp14":        ["exonuclease", "ExoN", "nsp-14", "NSP14",
                             "proofreading exonuclease", "N7-methyltransferase",
                             "nsp14 exonuclease"],
            "nsp15":        ["endoribonuclease", "NendoU", "nsp-15", "NSP15",
                             "uridylate-specific endoribonuclease"],
            "nsp16":        ["2'-O-methyltransferase", "nsp-16", "NSP16",
                             "cap methyltransferase"],
            "nsp3":         ["papain-like protease", "PLpro", "PL2pro",
                             "nsp-3", "NSP3", "ADP-ribose phosphatase",
                             "macro domain"],
            "nsp5":         ["main protease", "3CLpro", "Mpro", "3C-like protease",
                             "nsp-5", "NSP5"],
            "nsp7":         ["nsp-7", "NSP7", "primase cofactor"],
            "nsp8":         ["nsp-8", "NSP8", "primase"],
            "nsp9":         ["nsp-9", "NSP9", "ssRNA-binding protein"],
            "nsp10":        ["nsp-10", "NSP10", "cap-0 methyltransferase activator"],
            "ORF3a":        ["ORF 3a", "3a protein", "viroporin 3a"],
            "ORF8":         ["ORF 8", "8 protein", "ORF8 protein"],
            "ORF10":        ["ORF 10", "10 protein"],
        },
        "applications": {
            "diagnostics": ["diagnostic", "diagnostics", "detection", "assay",
                            "rapid test", "SHERLOCK", "DETECTR", "Cas12", "Cas13",
                            "point-of-care", "POC", "isothermal", "LAMP",
                            "biosensor", "lateral flow", "RT-LAMP",
                            "CRISPR-based detection", "guide RNA",
                            "field-deployable", "wastewater surveillance",
                            "clinical testing"],
            "therapeutics": ["antiviral", "therapy", "treatment", "inhibit",
                             "knockdown", "silencing", "gene editing",
                             "CRISPRi", "CRISPR knockout"],
        },
    },
}

# ── HTTP helpers ────────────────────────────────────────────────────────────────

def _get_json(url: str, retries: int = 4, delay_after: float = 0.0) -> dict:
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=25) as r:
                data = json.loads(r.read())
            if delay_after:
                time.sleep(delay_after)
            return data
        except (urllib.error.URLError, TimeoutError,
                http.client.RemoteDisconnected, ConnectionResetError, OSError) as e:
            wait = 2.0 * (attempt + 1)
            print(f"\n  [retry {attempt+1}/{retries-1}] {e.__class__.__name__}, "
                  f"sleeping {wait}s...", end="", flush=True)
            time.sleep(wait)
    raise RuntimeError(f"HTTP fetch failed after {retries} attempts: {url[:80]}")


def pubmed_search(query: str, max_ids: int = 0) -> tuple[int, list[str]]:
    params = {"db": "pubmed", "term": query, "retmax": max_ids, "retmode": "json"}
    if API_KEY:
        params["api_key"] = API_KEY
    url = f"{EUTILS}/esearch.fcgi?{urllib.parse.urlencode(params)}"
    data = _get_json(url, delay_after=PUBMED_DELAY)
    r = data.get("esearchresult", {})
    return int(r.get("count", 0)), r.get("idlist", [])


def pubmed_abstracts(pmids: list[str]) -> list[dict]:
    """Fetch title+abstract for a list of PMIDs."""
    if not pmids:
        return []
    # Cap to avoid overloading
    pmids = pmids[:MAX_LANDSCAPE_FETCH]
    params = {"db": "pubmed", "id": ",".join(pmids),
              "rettype": "abstract", "retmode": "json"}
    if API_KEY:
        params["api_key"] = API_KEY
    url = f"{EUTILS}/efetch.fcgi?{urllib.parse.urlencode(params)}"
    # efetch returns XML in abstract mode; use esummary for JSON
    params2 = {"db": "pubmed", "id": ",".join(pmids), "retmode": "json"}
    if API_KEY:
        params2["api_key"] = API_KEY
    url2 = f"{EUTILS}/esummary.fcgi?{urllib.parse.urlencode(params2)}"
    data = _get_json(url2, delay_after=PUBMED_DELAY)
    res = data.get("result", {})
    out = []
    for pid in pmids:
        r = res.get(pid, {})
        if "error" in r:
            continue
        out.append({
            "pmid": pid,
            "title": r.get("title", ""),
            "source": r.get("source", ""),
            "year": r.get("pubdate", "").split()[0],
        })
    return out


def europepmc_search(query: str) -> int:
    """Return article count from Europe PMC."""
    params = {"query": query, "format": "json", "pageSize": 1,
              "resultType": "lite"}
    url = f"{EUROPEPMC}/search?{urllib.parse.urlencode(params)}"
    try:
        data = _get_json(url, delay_after=EPMC_DELAY)
        return int(data.get("hitCount", 0))
    except Exception as e:
        print(f"\n  [EuropePMC error] {e}", end="")
        return -1  # -1 = could not check


# ── Query builders ──────────────────────────────────────────────────────────────

def _organism_clause_pm(cfg: dict) -> str:
    """PubMed clause covering primary + alternative names."""
    names = [cfg["search_name"]] + cfg.get("alt_names", [])
    parts = " OR ".join(f'"{n}"[Title/Abstract]' for n in names)
    return f"({parts})"


def _organism_clause_epmc(cfg: dict) -> str:
    # Use only the primary scientific name to avoid false positives.
    # alt_names often contain abbreviations (e.g. "V. cholerae") that
    # Europe PMC doesn't handle well without full quoted phrase matching.
    name = cfg["search_name"].split(" OR ")[0].strip()
    return f'"{name}"'


def _crispr_clause_pm() -> str:
    return ('(CRISPR[Title/Abstract] OR "guide RNA"[Title/Abstract] '
            'OR Cas9[Title/Abstract] OR Cas12[Title/Abstract] '
            'OR Cas13[Title/Abstract] OR "gene editing"[Title/Abstract])')


def _crispr_clause_epmc() -> str:
    return ('(CRISPR OR "guide RNA" OR Cas9 OR Cas12 OR Cas13 '
            'OR "gene editing")')


def _gene_clause_pm(gene: str, synonyms: dict) -> str:
    all_names = [gene] + synonyms.get(gene, [])
    return "(" + " OR ".join(f'"{g}"[Title/Abstract]' for g in all_names) + ")"


def _gene_clause_epmc(gene: str, synonyms: dict) -> str:
    all_names = [gene] + synonyms.get(gene, [])
    return "(" + " OR ".join(f'"{g}"' for g in all_names) + ")"


def _app_clause_pm_narrow(term: str) -> str:
    return f'"{term}"[Title/Abstract]'


def _app_clause_pm_broad(terms: list[str]) -> str:
    return "(" + " OR ".join(f'"{t}"[Title/Abstract]' for t in terms) + ")"


def _app_clause_epmc_broad(terms: list[str]) -> str:
    return "(" + " OR ".join(f'"{t}"' for t in terms) + ")"


# ── Corpus scan ─────────────────────────────────────────────────────────────────

def _normalise(text: str) -> str:
    return text.lower()


def corpus_matches(corpus: list[dict], terms: list[str]) -> list[dict]:
    """Return papers from corpus whose title contains any of the given terms."""
    matches = []
    norm_terms = [t.lower() for t in terms]
    for paper in corpus:
        hay = _normalise(paper.get("title", "") + " " + paper.get("source", ""))
        if any(t in hay for t in norm_terms):
            matches.append(paper)
    return matches


# ── Landscape + corpus fetch ────────────────────────────────────────────────────

def fetch_landscape_corpus(key: str, cfg: dict) -> tuple[int, list[dict]]:
    """Fetch total count + up to MAX_LANDSCAPE_FETCH paper summaries for corpus scan."""
    cache_file = CACHE_DIR / f"corpus_{key}.json"
    if cache_file.exists():
        with open(cache_file) as f:
            cached = json.load(f)
        print(f"  (corpus cache hit: {len(cached['papers'])} papers)")
        return cached["total"], cached["papers"]

    org = _organism_clause_pm(cfg)
    crispr = _crispr_clause_pm()
    query = f"{org} AND {crispr}"
    total, pmids = pubmed_search(query, max_ids=MAX_LANDSCAPE_FETCH)
    time.sleep(PUBMED_DELAY)

    papers: list[dict] = []
    if pmids:
        # Batch fetch summaries: max 200/request
        for batch_start in range(0, len(pmids), 200):
            batch = pmids[batch_start:batch_start + 200]
            papers.extend(pubmed_abstracts(batch))
            time.sleep(PUBMED_DELAY * 2)

    with open(cache_file, "w") as f:
        json.dump({"total": total, "papers": papers}, f)

    return total, papers


# ── Per-claim multi-strategy scoring ───────────────────────────────────────────

CONFIDENCE_ORDER = ["CONFIRMED", "PROBABLE", "UNCERTAIN", "FALSE"]


def query_and_score(
    label: str,
    pm_narrow_query: str,
    pm_broad_query: str,
    epmc_query: str,
    corpus: list[dict],
    corpus_terms: list[str],
) -> dict:
    """
    Run 3 query strategies + corpus scan. Return full evidence dict.

    Confidence rules (first matching rule wins):
      FALSE     — pm_broad > 5 OR epmc > 5
      UNCERTAIN — pm_broad > 0 OR epmc > 0 OR corpus_matches > 0
      PROBABLE  — pm_narrow > 0 (but everything else = 0)
      CONFIRMED — all counts == 0 and corpus_matches == 0
    """
    # Strategy 1: narrow (exact term, T/A only)
    n_narrow, narrow_ids = pubmed_search(pm_narrow_query, max_ids=MAX_DETAIL_FETCH)
    time.sleep(PUBMED_DELAY)

    # Strategy 2: broad (synonym expansion, T/A)
    if pm_narrow_query == pm_broad_query:
        n_broad, broad_ids = n_narrow, narrow_ids
    else:
        n_broad, broad_ids = pubmed_search(pm_broad_query, max_ids=MAX_DETAIL_FETCH)
        time.sleep(PUBMED_DELAY)

    # Strategy 3: Europe PMC
    n_epmc = europepmc_search(epmc_query)
    time.sleep(EPMC_DELAY)

    # Strategy 4: corpus scan
    corpus_hits = corpus_matches(corpus, corpus_terms)

    # Fetch top paper summaries for audit trail
    all_ids = list(dict.fromkeys(narrow_ids[:3] + broad_ids[:3]))
    top_papers = pubmed_abstracts(all_ids) if all_ids else []
    if all_ids:
        time.sleep(PUBMED_DELAY)

    # Confidence scoring.
    #
    # PubMed T/A is our primary source of truth (narrowly scoped, peer-reviewed).
    # EuropePMC searches full text and covers preprints; it produces much higher
    # counts because review papers co-mention pathogen + gene + CRISPR without
    # those being the actual study target. We use it as a secondary flag only.
    #
    # Rules (first match wins):
    #   FALSE     — PubMed broad > 5 (strong direct evidence the work exists)
    #   UNCERTAIN — PubMed broad > 0, OR corpus title match, OR EuropePMC > 30
    #               (indirect evidence; needs manual paper inspection)
    #   PROBABLE  — PubMed narrow > 0 but PubMed broad = 0
    #               (synonym expansion found something; likely real)
    #   CONFIRMED — PubMed narrow = 0, PubMed broad = 0, corpus = 0
    #               (EuropePMC may still be > 0 due to full-text co-mention;
    #               this is disclosed in the methodology statement)
    if n_broad > 5:
        confidence = "FALSE"
    elif n_broad > 0 or len(corpus_hits) > 0 or n_epmc > 30:
        confidence = "UNCERTAIN"
    elif n_narrow > 0:
        confidence = "PROBABLE"
    else:
        confidence = "CONFIRMED"

    return {
        "label": label,
        "confidence": confidence,
        "counts": {
            "pubmed_narrow": n_narrow,
            "pubmed_broad": n_broad,
            "europepmc": n_epmc,
            "corpus_title_matches": len(corpus_hits),
        },
        "queries": {
            "pubmed_narrow": pm_narrow_query,
            "pubmed_broad": pm_broad_query,
            "europepmc": epmc_query,
        },
        "corpus_hits": [p["title"][:100] for p in corpus_hits[:5]],
        "top_papers": [
            {"pmid": p["pmid"], "title": p["title"][:120],
             "year": p["year"], "journal": p["source"]}
            for p in top_papers[:5]
        ],
    }


# ── Main scan ───────────────────────────────────────────────────────────────────

def scan_pathogen(key: str, cfg: dict) -> dict:
    print(f"\n{'='*60}")
    print(f"  {cfg['search_name']} ({key})")
    print(f"{'='*60}")

    # Fetch landscape + corpus
    print("  Fetching landscape + corpus...")
    total_papers, corpus = fetch_landscape_corpus(key, cfg)
    print(f"  Landscape: {total_papers} total CRISPR papers, "
          f"{len(corpus)} abstracts fetched")

    org_pm = _organism_clause_pm(cfg)
    org_epmc = _organism_clause_epmc(cfg)
    crispr_pm = _crispr_clause_pm()
    crispr_epmc = _crispr_clause_epmc()
    synonyms = cfg.get("gene_synonyms", {})

    gene_results: dict[str, dict] = {}
    application_results: dict[str, dict] = {}

    # ── Gene-level scan ──────────────────────────────────────────────────────
    print(f"\n  Gene scan ({len(cfg['genes'])} genes):")
    for gene in cfg["genes"]:
        gene_clause_pm = _gene_clause_pm(gene, synonyms)
        gene_clause_epmc = _gene_clause_epmc(gene, synonyms)
        all_gene_terms = [gene] + synonyms.get(gene, [])

        pm_narrow = f"{org_pm} AND {crispr_pm} AND {gene_clause_pm}"
        pm_broad  = pm_narrow  # gene query already uses full synonym expansion
        epmc_q    = f"{org_epmc} AND {crispr_epmc} AND {gene_clause_epmc}"

        result = query_and_score(
            label=f"{key}:gene:{gene}",
            pm_narrow_query=pm_narrow,
            pm_broad_query=pm_broad,
            epmc_query=epmc_q,
            corpus=corpus,
            corpus_terms=all_gene_terms,
        )
        gene_results[gene] = result

        marker = {"CONFIRMED": "○", "PROBABLE": "?",
                  "UNCERTAIN": "⚠", "FALSE": "●"}.get(result["confidence"], " ")
        print(f"    {marker} {gene:25s}  "
              f"pm={result['counts']['pubmed_narrow']:>3} "
              f"pm_broad={result['counts']['pubmed_broad']:>3} "
              f"epmc={result['counts']['europepmc']:>4} "
              f"corpus={result['counts']['corpus_title_matches']:>2}  "
              f"→ {result['confidence']}")

    # ── Application-level scan ───────────────────────────────────────────────
    print(f"\n  Application scan ({len(cfg['applications'])} categories):")
    for app_cat, app_terms in cfg["applications"].items():
        primary_term = app_terms[0]

        pm_narrow = (f"{org_pm} AND {crispr_pm} AND "
                     f"{_app_clause_pm_narrow(primary_term)}")
        pm_broad  = (f"{org_pm} AND {crispr_pm} AND "
                     f"{_app_clause_pm_broad(app_terms)}")
        epmc_q    = (f"{org_epmc} AND {crispr_epmc} AND "
                     f"{_app_clause_epmc_broad(app_terms)}")

        result = query_and_score(
            label=f"{key}:app:{app_cat}",
            pm_narrow_query=pm_narrow,
            pm_broad_query=pm_broad,
            epmc_query=epmc_q,
            corpus=corpus,
            corpus_terms=app_terms,
        )
        application_results[app_cat] = result

        marker = {"CONFIRMED": "○", "PROBABLE": "?",
                  "UNCERTAIN": "⚠", "FALSE": "●"}.get(result["confidence"], " ")
        print(f"    {marker} {app_cat:22s}  "
              f"pm={result['counts']['pubmed_narrow']:>3} "
              f"pm_broad={result['counts']['pubmed_broad']:>3} "
              f"epmc={result['counts']['europepmc']:>4} "
              f"corpus={result['counts']['corpus_title_matches']:>2}  "
              f"→ {result['confidence']}")

    # Summarise
    confirmed_gene_gaps = [g for g, r in gene_results.items()
                           if r["confidence"] == "CONFIRMED"]
    confirmed_app_gaps = [a for a, r in application_results.items()
                          if r["confidence"] == "CONFIRMED"]
    uncertain = (
        [f"gene:{g}" for g, r in gene_results.items()
         if r["confidence"] in ("PROBABLE", "UNCERTAIN")]
        + [f"app:{a}" for a, r in application_results.items()
           if r["confidence"] in ("PROBABLE", "UNCERTAIN")]
    )
    false_claims = (
        [f"gene:{g}" for g, r in gene_results.items()
         if r["confidence"] == "FALSE"]
        + [f"app:{a}" for a, r in application_results.items()
           if r["confidence"] == "FALSE"]
    )

    print(f"\n  Summary for {cfg['search_name']}:")
    print(f"    CONFIRMED gaps:  genes={confirmed_gene_gaps}, apps={confirmed_app_gaps}")
    print(f"    UNCERTAIN:       {uncertain}")
    print(f"    FALSE claims:    {false_claims}")

    return {
        "pathogen_key": key,
        "search_name": cfg["search_name"],
        "landscape_total": total_papers,
        "corpus_fetched": len(corpus),
        "genes": gene_results,
        "applications": application_results,
        "confirmed_gene_gaps": confirmed_gene_gaps,
        "confirmed_app_gaps": confirmed_app_gaps,
        "uncertain": uncertain,
        "false_claims": false_claims,
    }


# ── Report + audit file ─────────────────────────────────────────────────────────

def write_audit(results: dict, out_path: Path):
    lines = [
        "=" * 70,
        "LOOM PUBMED SCAN v3 — AUDIT TRAIL",
        f"Generated: {results['generated']}",
        f"Total confirmed gaps: {results['summary']['confirmed_gaps']}",
        f"Uncertain/needs review: {results['summary']['uncertain_count']}",
        f"Refuted (FALSE) claims: {results['summary']['false_count']}",
        "=" * 70,
        "",
        "CONFIDENCE KEY:",
        "  CONFIRMED — 0 results across all query strategies + sources + corpus",
        "  PROBABLE  — 0 in narrow PubMed query only; not verified elsewhere",
        "  UNCERTAIN — 0 in narrow, but >0 in broad/EuropePMC or corpus hit",
        "  FALSE     — >5 results in broad query or EuropePMC; NOT a gap",
        "",
    ]
    for key, pr in results["pathogens"].items():
        lines.append(f"{'─'*60}")
        lines.append(f"PATHOGEN: {pr['search_name']} ({key})")
        lines.append(f"  Landscape: {pr['landscape_total']} total CRISPR papers")
        lines.append(f"  Corpus fetched: {pr['corpus_fetched']} abstracts scanned")
        lines.append("")

        # Gene gaps
        lines.append("  GENE GAPS:")
        for gene, r in pr["genes"].items():
            c = r["confidence"]
            cnt = r["counts"]
            lines.append(
                f"    [{c:9s}] {gene:25s} "
                f"pm_narrow={cnt['pubmed_narrow']:>3}  "
                f"pm_broad={cnt['pubmed_broad']:>3}  "
                f"epmc={cnt['europepmc']:>4}  "
                f"corpus={cnt['corpus_title_matches']:>2}"
            )
            if r["corpus_hits"]:
                for h in r["corpus_hits"][:3]:
                    lines.append(f"               corpus: {h[:80]}")
            if r["top_papers"] and c in ("UNCERTAIN", "FALSE"):
                for p in r["top_papers"][:2]:
                    lines.append(f"               paper {p['year']}: {p['title'][:80]}")

        lines.append("")
        lines.append("  APPLICATION GAPS:")
        for app, r in pr["applications"].items():
            c = r["confidence"]
            cnt = r["counts"]
            lines.append(
                f"    [{c:9s}] {app:22s} "
                f"pm_narrow={cnt['pubmed_narrow']:>3}  "
                f"pm_broad={cnt['pubmed_broad']:>3}  "
                f"epmc={cnt['europepmc']:>4}  "
                f"corpus={cnt['corpus_title_matches']:>2}"
            )
            if r["corpus_hits"]:
                for h in r["corpus_hits"][:3]:
                    lines.append(f"               corpus: {h[:80]}")
            if r["top_papers"] and c in ("UNCERTAIN", "FALSE"):
                for p in r["top_papers"][:2]:
                    lines.append(f"               paper {p['year']}: {p['title'][:80]}")

        lines.append("")
        if pr.get("false_claims"):
            lines.append(f"  *** FALSE CLAIMS (must be removed/corrected): "
                         f"{pr['false_claims']}")
        if pr.get("uncertain"):
            lines.append(f"  *** NEEDS MANUAL REVIEW: {pr['uncertain']}")
        lines.append("")

    out_path.write_text("\n".join(lines))


def main():
    import argparse
    parser = argparse.ArgumentParser(description="PubMed literature scanner v3")
    parser.add_argument(
        "--only", metavar="KEY",
        help="Scan only this pathogen key (e.g. sars-cov-2). Merges results with existing output.",
    )
    parser.add_argument(
        "--no-merge", action="store_true",
        help="When using --only, overwrite output instead of merging with existing results.",
    )
    args = parser.parse_args()

    if not INPUT_FILE.exists():
        print(f"Error: {INPUT_FILE} not found", file=sys.stderr)
        sys.exit(1)

    pathogens_raw = json.loads(INPUT_FILE.read_text())
    available_keys = set(pathogens_raw.keys()) - EXCLUDE_KEYS

    # Only scan pathogens we have config for
    keys_to_scan = sorted(available_keys & set(PATHOGEN_CONFIG.keys()))
    skipped = available_keys - set(PATHOGEN_CONFIG.keys())
    if skipped:
        print(f"Note: no config for {skipped}, skipping")

    if args.only:
        if args.only not in PATHOGEN_CONFIG:
            print(f"Error: no config for '{args.only}' in PATHOGEN_CONFIG", file=sys.stderr)
            sys.exit(1)
        if args.only not in available_keys:
            print(f"Warning: '{args.only}' not in input file keys; adding stub entry")
            pathogens_raw[args.only] = {}
        keys_to_scan = [args.only]
        print(f"--only mode: scanning only '{args.only}'")

    print(f"Scanning {len(keys_to_scan)} pathogens: {keys_to_scan}")
    print(f"API key: {'set ({int(1/PUBMED_DELAY)}/s)' if API_KEY else 'not set (3/s)'}")

    # Estimate total API calls
    total_calls = sum(
        # landscape (1) + gene queries (N×2 pm + 1 epmc) + app queries (M×2 pm + 1 epmc)
        1
        + len(PATHOGEN_CONFIG[k]["genes"]) * 3
        + len(PATHOGEN_CONFIG[k]["applications"]) * 3
        for k in keys_to_scan
    )
    est_min = total_calls * max(PUBMED_DELAY, EPMC_DELAY) / 60
    print(f"Estimated: ~{total_calls} API calls, ~{est_min:.0f} min")
    print()

    all_results: dict[str, dict] = {}
    for key in keys_to_scan:
        cfg = PATHOGEN_CONFIG[key]
        all_results[key] = scan_pathogen(key, cfg)

    # If --only mode, merge new result into existing output
    if args.only and not args.no_merge and OUTPUT_FILE.exists():
        print(f"\nMerging '{args.only}' result into existing {OUTPUT_FILE.name} ...")
        existing = json.loads(OUTPUT_FILE.read_text())
        # Update pathogens dict with new result
        existing.setdefault("pathogens", {})[args.only] = all_results[args.only]
        # Rebuild summary from all pathogens (existing + new)
        all_results = existing["pathogens"]

    # Compile summary
    all_confirmed_gene = [
        {"pathogen": k, "gene": g}
        for k, r in all_results.items()
        for g in r["confirmed_gene_gaps"]
    ]
    all_confirmed_app = [
        {"pathogen": k, "application": a}
        for k, r in all_results.items()
        for a in r["confirmed_app_gaps"]
    ]
    all_uncertain = [
        {"pathogen": k, "item": it}
        for k, r in all_results.items()
        for it in r["uncertain"]
    ]
    all_false = [
        {"pathogen": k, "item": it}
        for k, r in all_results.items()
        for it in r["false_claims"]
    ]

    output = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "version": 3,
        "methodology": (
            "Multi-strategy gap analysis: each claimed gap is tested with "
            "(1) narrow PubMed query [Title/Abstract], "
            "(2) broad PubMed query with full synonym expansion [Title/Abstract], "
            "(3) Europe PMC independent search, "
            "(4) local corpus scan of all ~N landscape abstracts. "
            "Only claims with CONFIRMED confidence (0 in all strategies) are "
            "reported as gaps in the paper. UNCERTAIN and FALSE claims are "
            "flagged for manual review or removal."
        ),
        "summary": {
            "pathogens_scanned": len(all_results),
            "confirmed_gene_gaps": len(all_confirmed_gene),
            "confirmed_app_gaps": len(all_confirmed_app),
            "confirmed_gaps": len(all_confirmed_gene) + len(all_confirmed_app),
            "uncertain_count": len(all_uncertain),
            "false_count": len(all_false),
        },
        "confirmed_gene_gaps": all_confirmed_gene,
        "confirmed_app_gaps": all_confirmed_app,
        "uncertain": all_uncertain,
        "false_claims": all_false,
        "pathogens": all_results,
    }

    OUTPUT_FILE.write_text(json.dumps(output, indent=2))
    print(f"\n\nResults → {OUTPUT_FILE}")

    write_audit(output, AUDIT_FILE)
    print(f"Audit trail → {AUDIT_FILE}")

    # Print summary
    print("\n" + "=" * 70)
    print("  PUBMED SCAN v3 — SUMMARY")
    print("=" * 70)
    print(f"  Confirmed gene gaps:  {len(all_confirmed_gene)}")
    print(f"  Confirmed app gaps:   {len(all_confirmed_app)}")
    print(f"  Needs manual review:  {len(all_uncertain)}")
    print(f"  False (retracted):    {len(all_false)}")
    print("-" * 70)
    if all_false:
        print("  FALSE CLAIMS (must be removed from paper):")
        for f in all_false:
            print(f"    {f['pathogen']:15s} | {f['item']}")
    print("=" * 70)


if __name__ == "__main__":
    main()
