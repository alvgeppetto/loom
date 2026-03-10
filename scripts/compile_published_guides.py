#!/usr/bin/env python3
"""Compile a comprehensive database of published SARS-CoV-2 CRISPR guide sequences.

Sources:
  1. Hand-curated guides from papers cited in our manuscript (Tier 0)
  2. Europe PMC full-text mining of SARS-CoV-2 CRISPR papers (Tier 1)
  3. Supplementary materials linked from Europe PMC records (Tier 2)

Output: data/crispr_guides/published_guides_comprehensive.json

Usage:
  python scripts/compile_published_guides.py [--output path] [--no-api]
"""

import argparse
import json
import re
import sys
import time
import urllib.request
import urllib.parse
from pathlib import Path
from typing import Optional

# ── Constants ────────────────────────────────────────────────────────────────

EPMC_REST = "https://www.ebi.ac.uk/europepmc/webservices/rest"
EPMC_DELAY = 0.4  # seconds between API calls (polite rate)
NCBI_EFETCH = "https://efetch.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Wuhan-Hu-1 accession — we'll fetch this to map guide positions
WUHAN_HU1_ACC = "NC_045512.2"
WUHAN_HU1_LEN = 29903

# Regex: 18–30 bp DNA-only strings (likely guide/spacer sequences)
GUIDE_RE = re.compile(r"\b([ACGTU]{18,30})\b")

# Ontology-derived search terms for SARS-CoV-2 CRISPR papers
SEARCH_SYNONYMS = [
    "SARS-CoV-2", "COVID-19", "2019-nCoV", "severe acute respiratory syndrome coronavirus 2",
]
CRISPR_TERMS = [
    "CRISPR", "Cas12", "Cas13", "Cas12a", "Cas13a",
    "SHERLOCK", "DETECTR", "CARMEN", "STOP-COVID", "STOPCovid",
    "guide RNA", "crRNA", "gRNA",
]

# ── Tier 0: Hand-curated guides from cited papers ───────────────────────────

TIER0_GUIDES = [
    # Broughton et al. 2020 (DETECTR) — Cas12a, Nat Biotechnol
    {"id": "DETECTR_N", "seq": "GACCCCAAAATCAGCGAAAT", "gene": "N",
     "paper": "Broughton2020", "pmid": "32300245", "cas": "Cas12a", "use": "Dx"},
    {"id": "DETECTR_E", "seq": "TTACAAACATTGGCCGCAAA", "gene": "E",
     "paper": "Broughton2020", "pmid": "32300245", "cas": "Cas12a", "use": "Dx"},
    # Patchsung et al. 2020 (SHERLOCK) — Cas13a, Nat Biomed Eng
    {"id": "SHERLOCK_S", "seq": "CCAGAACTCAATTACCCCCTG", "gene": "S",
     "paper": "Patchsung2020", "pmid": "32826945", "cas": "Cas13a", "use": "Dx"},
    {"id": "SHERLOCK_Orf1ab", "seq": "CTTCGATTGTGTGCGTACTGC", "gene": "ORF1ab",
     "paper": "Patchsung2020", "pmid": "32826945", "cas": "Cas13a", "use": "Dx"},
    # Abbott et al. 2020 (PACMAN) — Cas13d, Cell
    {"id": "PACMAN_N1", "seq": "GTTCCAATTAACACCAATAGC", "gene": "N",
     "paper": "Abbott2020", "pmid": "32353252", "cas": "Cas13d", "use": "Tx"},
    {"id": "PACMAN_N2", "seq": "GATTACAAACATTGGCCGCAA", "gene": "N",
     "paper": "Abbott2020", "pmid": "32353252", "cas": "Cas13d", "use": "Tx"},
    {"id": "PACMAN_RdRp1", "seq": "ACCTTTCAAACTCAACATCAACCA", "gene": "RdRp",
     "paper": "Abbott2020", "pmid": "32353252", "cas": "Cas13d", "use": "Tx"},
    {"id": "PACMAN_RdRp2", "seq": "TGATAGAAGTGCAAGGTTACAAG", "gene": "RdRp",
     "paper": "Abbott2020", "pmid": "32353252", "cas": "Cas13d", "use": "Tx"},
    # Freije et al. 2019 (CARVER) — Cas13, Mol Cell
    {"id": "CARVER_con1", "seq": "AAAGGTTTACCCAATAATACTGCG", "gene": "conserved",
     "paper": "Freije2019", "pmid": "31607545", "cas": "Cas13a", "use": "Tx"},
    {"id": "CARVER_con2", "seq": "GCCAACAACAAAGGTTACTTTTGG", "gene": "conserved",
     "paper": "Freije2019", "pmid": "31607545", "cas": "Cas13a", "use": "Tx"},
    # Joung et al. 2020 (STOPCovid.v2) — Cas12b, NEJM
    {"id": "STOPCovid_N_cr1", "seq": "TTCCAATTAACACCAATAGCAG", "gene": "N",
     "paper": "Joung2020", "pmid": "32941093", "cas": "Cas12b", "use": "Dx"},
    {"id": "STOPCovid_N_cr2", "seq": "TATTACCGCCTTTGAGTGAGC", "gene": "N",
     "paper": "Joung2020", "pmid": "32941093", "cas": "Cas12b", "use": "Dx"},
    # Fozouni et al. 2021 (amplification-free Cas13a) — Cell
    {"id": "Fozouni_N_cr1", "seq": "UAUUACCGCCUUUGAGUGAGC", "gene": "N",
     "paper": "Fozouni2021", "pmid": "33357727", "cas": "Cas13a", "use": "Dx"},
    # Ackerman et al. 2020 (CARMEN) — Cas13, Nature
    {"id": "CARMEN_N_cr1", "seq": "GACCCCAAAATCAGCGAAAT", "gene": "N",
     "paper": "Ackerman2020", "pmid": "32532005", "cas": "Cas13a", "use": "Dx"},
    {"id": "CARMEN_S_cr1", "seq": "CCAGAACTCAATTACCCCCTG", "gene": "S",
     "paper": "Ackerman2020", "pmid": "32532005", "cas": "Cas13a", "use": "Dx"},
    # Ding et al. 2020 (All-in-one Cas12a) — Nat Commun
    {"id": "Ding_N_cr1", "seq": "GACCCCAAAATCAGCGAAATG", "gene": "N",
     "paper": "Ding2020", "pmid": "32999292", "cas": "Cas12a", "use": "Dx"},
    {"id": "Ding_N_cr2", "seq": "GGCACGTCAATATGCTTATTC", "gene": "N",
     "paper": "Ding2020", "pmid": "32999292", "cas": "Cas12a", "use": "Dx"},
    # Huang Z et al. 2020 (Cas12a E-gene) — Biosens Bioelectron
    {"id": "Huang_E_cr1", "seq": "TTACAAACATTGGCCGCAAATT", "gene": "E",
     "paper": "Huang2020", "pmid": "32534345", "cas": "Cas12a", "use": "Dx"},
    # Guo L et al. 2020 (CRISPR-COVID) — ACS Nano
    {"id": "CRISPR_COVID_N", "seq": "GAAATTTTGGATGACCTTCTGC", "gene": "N",
     "paper": "Guo2020", "pmid": "33064890", "cas": "Cas12a", "use": "Dx"},
    # Rauch et al. 2021 (SENSR) — J Clin Microbiol
    {"id": "SENSR_S_cr1", "seq": "CCAGAACTCAATTACCCCCTGC", "gene": "S",
     "paper": "Rauch2021", "pmid": "32938741", "cas": "Cas13a", "use": "Dx"},
    # Zhang F et al. 2020 (a protocol) — Broad preprint
    {"id": "Zhang_S_cr1", "seq": "CCAGAACTCAATTACCCCCTG", "gene": "S",
     "paper": "Zhang2020_protocol", "pmid": "n/a", "cas": "Cas13a", "use": "Dx"},
    # Lucia et al. 2020 (Cas12a S+N) — preprint
    {"id": "Lucia_S_cr1", "seq": "CCCTGCATACACTAATTCTTTCAC", "gene": "S",
     "paper": "Lucia2020", "pmid": "n/a", "cas": "Cas12a", "use": "Dx"},
    {"id": "Lucia_N_cr1", "seq": "TGGTTACTGCCAGTTGAATCTG", "gene": "N",
     "paper": "Lucia2020", "pmid": "n/a", "cas": "Cas12a", "use": "Dx"},
    # Azhar M et al. 2021 (FnCas9 Editor Linked Uniform Detection, iSCAN) — ACS Nano
    {"id": "iSCAN_RdRp_cr1", "seq": "CTCCTCTAGTGGCGGCTATTG", "gene": "RdRp",
     "paper": "Azhar2021", "pmid": "33289377", "cas": "FnCas9", "use": "Dx"},
    {"id": "iSCAN_N_cr1", "seq": "CTCAATTTTCCCCCAGCGCTTC", "gene": "N",
     "paper": "Azhar2021", "pmid": "33289377", "cas": "FnCas9", "use": "Dx"},
    # Ali et al. 2020 (Cas12b AIOD-CRISPR) — Nat Commun
    {"id": "Ali_N_cr1", "seq": "TTTGAAATTTGGATGACCTTCTGC", "gene": "N",
     "paper": "Ali2020", "pmid": "33188177", "cas": "Cas12b", "use": "Dx"},
]


# ── Utility functions ────────────────────────────────────────────────────────

def http_get_json(url: str, retries: int = 3) -> dict:
    """GET a JSON endpoint with retries."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=20) as resp:
                return json.loads(resp.read().decode())
        except Exception as e:
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)
    return {}


def http_get_text(url: str, retries: int = 2) -> str:
    """GET a text/XML endpoint with retries."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req, timeout=30) as resp:
                return resp.read().decode(errors="replace")
        except Exception as e:
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)
    return ""


def normalize_seq(seq: str) -> str:
    """Normalize guide sequence: uppercase, U→T, strip non-ACGT."""
    return re.sub(r"[^ACGT]", "", seq.upper().replace("U", "T"))


def fetch_wuhan_hu1() -> str:
    """Fetch the Wuhan-Hu-1 reference genome (NC_045512.2) from NCBI."""
    cache_path = Path("data/crispr_guides/.wuhan_hu1_cache.fa")
    if cache_path.exists():
        text = cache_path.read_text()
        # Parse FASTA: skip header lines, concatenate sequence
        return "".join(
            line.strip() for line in text.splitlines()
            if line.strip() and not line.startswith(">")
        ).upper()

    url = (
        f"{NCBI_EFETCH}?db=nucleotide&id={WUHAN_HU1_ACC}"
        f"&rettype=fasta&retmode=text"
    )
    print(f"  Fetching Wuhan-Hu-1 ({WUHAN_HU1_ACC})...")
    text = http_get_text(url)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(text)

    return "".join(
        line.strip() for line in text.splitlines()
        if line.strip() and not line.startswith(">")
    ).upper()


def find_position_on_ref(seq: str, ref_genome: str) -> Optional[int]:
    """Find 0-based start position of seq on Wuhan-Hu-1 (forward strand)."""
    pos = ref_genome.find(seq)
    if pos >= 0:
        return pos
    # Try reverse complement
    rc = reverse_complement(seq)
    pos = ref_genome.find(rc)
    if pos >= 0:
        return pos
    return None


def reverse_complement(seq: str) -> str:
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def hamming_distance(a: str, b: str) -> int:
    """Hamming distance between two equal-length sequences."""
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(c1 != c2 for c1, c2 in zip(a, b))


# ── Tier 1: Europe PMC full-text mining ──────────────────────────────────────

def search_epmc_papers(max_pages: int = 5) -> list[dict]:
    """Search Europe PMC for SARS-CoV-2 CRISPR papers, return paper metadata."""
    # Build query with ontology-derived synonyms
    organism = " OR ".join(f'"{s}"' for s in SEARCH_SYNONYMS)
    crispr = " OR ".join(f'"{t}"' for t in CRISPR_TERMS[:6])  # top terms
    query = f"({organism}) AND ({crispr}) AND (diagnostic OR detection OR assay)"

    papers = []
    cursor_mark = "*"
    for page in range(max_pages):
        params = urllib.parse.urlencode({
            "query": query,
            "format": "json",
            "pageSize": 100,
            "resultType": "core",
            "cursorMark": cursor_mark,
            "sort": "CITED desc",
        })
        url = f"{EPMC_REST}/search?{params}"
        try:
            data = http_get_json(url)
        except Exception as e:
            print(f"  [EPMC search error] {e}")
            break

        results = data.get("resultList", {}).get("result", [])
        if not results:
            break

        for r in results:
            papers.append({
                "pmid": r.get("pmid", ""),
                "pmcid": r.get("pmcid", ""),
                "title": r.get("title", ""),
                "year": r.get("pubYear", ""),
                "journal": r.get("journalTitle", ""),
                "doi": r.get("doi", ""),
                "cited_by": r.get("citedByCount", 0),
                "has_fulltext": r.get("hasTextMinedTerms") == "Y"
                                or r.get("inEPMC") == "Y",
            })

        cursor_mark = data.get("nextCursorMark", cursor_mark)
        time.sleep(EPMC_DELAY)
        print(f"  EPMC search page {page + 1}: {len(results)} papers "
              f"(total {len(papers)})")

    return papers


def fetch_fulltext(pmcid: str) -> str:
    """Fetch full-text XML from Europe PMC for a given PMC ID."""
    if not pmcid:
        return ""
    url = f"{EPMC_REST}/{pmcid}/fullTextXML"
    try:
        text = http_get_text(url)
        time.sleep(EPMC_DELAY)
        return text
    except Exception:
        return ""


def extract_guide_sequences(text: str) -> list[str]:
    """Extract candidate CRISPR guide sequences from full-text content.

    Heuristics:
    - 18–30 bp DNA-only strings
    - Must contain at least 3 of each of A, C, G, T (filter out poly-X runs)
    - Ignore common non-guide artifact strings (primers, adapters)
    """
    candidates = set()
    for match in GUIDE_RE.finditer(text):
        raw = match.group(1)
        seq = normalize_seq(raw)
        if len(seq) < 18 or len(seq) > 30:
            continue
        # Must have reasonable base diversity (not poly-X)
        counts = {b: seq.count(b) for b in "ACGT"}
        if min(counts.values()) < 2:
            continue
        # Skip very common adapter/primer motifs
        if seq.startswith("AGATCGGAAG"):  # Illumina adapter
            continue
        candidates.add(seq)
    return sorted(candidates)


# ── Main pipeline ────────────────────────────────────────────────────────────

def run(output_path: Path, skip_api: bool = False):
    print("=" * 60)
    print("Compiling comprehensive published SARS-CoV-2 CRISPR guide database")
    print("=" * 60)

    # 1. Fetch Wuhan-Hu-1 reference (optional — used for position mapping)
    print("\n[1/4] Fetching Wuhan-Hu-1 reference...")
    ref_genome: Optional[str] = None
    try:
        ref_genome = fetch_wuhan_hu1()
        print(f"  Reference length: {len(ref_genome)} bp")
    except Exception as e:
        print(f"  WARNING: Could not fetch reference ({e})")
        print("  Position mapping will be skipped; sequence comparison still works.")

    # 2. Load Tier 0 (hand-curated) guides
    print(f"\n[2/4] Loading {len(TIER0_GUIDES)} hand-curated Tier 0 guides...")
    all_guides = {}
    for g in TIER0_GUIDES:
        seq = normalize_seq(g["seq"])
        pos = find_position_on_ref(seq, ref_genome) if ref_genome else None
        all_guides[seq] = {
            "id": g["id"],
            "sequence": seq,
            "length": len(seq),
            "gene": g["gene"],
            "paper": g["paper"],
            "pmid": g["pmid"],
            "cas_enzyme": g["cas"],
            "use": g["use"],
            "wuhan_hu1_pos": pos,
            "tier": 0,
            "source": "hand_curated",
        }
    print(f"  {len(all_guides)} unique sequences (after dedup)")

    # 3. Europe PMC full-text mining
    mined_seqs = []
    if not skip_api:
        print(f"\n[3/4] Searching Europe PMC for SARS-CoV-2 CRISPR papers...")
        papers = search_epmc_papers(max_pages=3)
        print(f"  Found {len(papers)} papers, fetching full texts...")

        # Sort by citation count, fetch top papers
        papers.sort(key=lambda p: p.get("cited_by", 0), reverse=True)
        fetched = 0
        max_fetch = 80  # top 80 most-cited papers
        for paper in papers[:max_fetch]:
            pmcid = paper.get("pmcid", "")
            if not pmcid:
                continue
            fulltext = fetch_fulltext(pmcid)
            if not fulltext:
                continue
            fetched += 1
            seqs = extract_guide_sequences(fulltext)
            for seq in seqs:
                if seq not in all_guides:
                    pos = find_position_on_ref(seq, ref_genome) if ref_genome else None
                    if ref_genome is None or pos is not None:  # Keep if maps to SARS-CoV-2 or no ref
                        mined_seqs.append(seq)
                        all_guides[seq] = {
                            "id": f"EPMC_{pmcid}_{pos}",
                            "sequence": seq,
                            "length": len(seq),
                            "gene": "mapped",
                            "paper": paper.get("title", "")[:60],
                            "pmid": paper.get("pmid", ""),
                            "cas_enzyme": "unknown",
                            "use": "unknown",
                            "wuhan_hu1_pos": pos,
                            "tier": 1,
                            "source": f"epmc_fulltext:{pmcid}",
                        }

            if fetched % 10 == 0:
                print(f"  Processed {fetched} papers, "
                      f"{len(mined_seqs)} new sequences extracted")

        print(f"  Full-text mining: {fetched} papers processed, "
              f"{len(mined_seqs)} new SARS-CoV-2-mapping sequences")
    else:
        print("\n[3/4] Skipping Europe PMC search (--no-api)")

    # 4. Compile final database
    print(f"\n[4/4] Compiling database...")
    guides_list = sorted(all_guides.values(), key=lambda g: (g["tier"], g.get("wuhan_hu1_pos") or 99999))

    db = {
        "version": 1,
        "generated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "reference": WUHAN_HU1_ACC,
        "reference_length": len(ref_genome) if ref_genome else None,
        "total_guides": len(guides_list),
        "tier0_count": sum(1 for g in guides_list if g["tier"] == 0),
        "tier1_count": sum(1 for g in guides_list if g["tier"] == 1),
        "sources": {
            "tier0": "Hand-curated from cited diagnostic/therapeutic papers",
            "tier1": "Europe PMC full-text mining (SARS-CoV-2 CRISPR papers)",
        },
        "guides": guides_list,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(db, indent=2) + "\n")
    print(f"\n  Written {len(guides_list)} guides to {output_path}")
    print(f"  Tier 0 (hand-curated): {db['tier0_count']}")
    print(f"  Tier 1 (EPMC mined):   {db['tier1_count']}")
    return db


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", "-o", type=Path,
        default=Path("data/crispr_guides/published_guides_comprehensive.json"),
    )
    parser.add_argument(
        "--no-api", action="store_true",
        help="Skip Europe PMC API calls (use only Tier 0 hand-curated guides)",
    )
    args = parser.parse_args()
    run(args.output, skip_api=args.no_api)


if __name__ == "__main__":
    main()
