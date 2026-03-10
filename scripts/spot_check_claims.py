#!/usr/bin/env python3
"""Quick spot-check of highest-risk claims against live PubMed."""
import json
import time
import urllib.parse
import urllib.request


def check(label, query):
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
        + urllib.parse.urlencode(
            {"db": "pubmed", "term": query, "retmax": 5, "retmode": "json"}
        )
    )
    resp = urllib.request.urlopen(url, timeout=15)
    data = json.loads(resp.read())
    result = data.get("esearchresult", {})
    count = result.get("count", "0")
    ids = result.get("idlist", [])
    status = "FALSE GAP" if int(count) > 0 else "confirmed gap"
    print(f"  [{status:12s}] {label:50s} → {count:>3} papers  PMIDs: {ids[:3]}")
    time.sleep(0.4)


print("=== MPOX (corrected query — was 0 due to bug) ===")
check(
    "Mpox/Monkeypox + CRISPR (landscape)",
    '(Mpox[Title/Abstract] OR Monkeypox[Title/Abstract]) AND '
    '(CRISPR[Title/Abstract] OR "guide RNA"[Title/Abstract] '
    'OR Cas9[Title/Abstract] OR Cas12[Title/Abstract] OR Cas13[Title/Abstract])',
)
check(
    "Mpox + CRISPR + diagnostic",
    '(Mpox[Title/Abstract] OR Monkeypox[Title/Abstract]) AND '
    'CRISPR[Title/Abstract] AND ("diagnostic"[Title/Abstract])',
)
check(
    "Mpox + CRISPR + detection",
    '(Mpox[Title/Abstract] OR Monkeypox[Title/Abstract]) AND '
    'CRISPR[Title/Abstract] AND ("detection"[Title/Abstract])',
)
check(
    "Mpox + CRISPR + SHERLOCK",
    '(Mpox[Title/Abstract] OR Monkeypox[Title/Abstract]) AND '
    'CRISPR[Title/Abstract] AND ("SHERLOCK"[Title/Abstract])',
)
check(
    "Mpox + CRISPR + point-of-care",
    '(Mpox[Title/Abstract] OR Monkeypox[Title/Abstract]) AND '
    'CRISPR[Title/Abstract] AND ("point-of-care"[Title/Abstract])',
)

print("\n=== TB GENE GAPS (claims: all zero) ===")
check(
    "TB + CRISPR + inhA (expanded)",
    '("Mycobacterium tuberculosis"[Title/Abstract]) AND CRISPR[Title/Abstract] '
    'AND (inhA[Title/Abstract] OR "enoyl-ACP reductase"[Title/Abstract] '
    'OR "isoniazid resistance"[Title/Abstract])',
)
check(
    "TB + CRISPR + gyrB (expanded)",
    '("Mycobacterium tuberculosis"[Title/Abstract]) AND CRISPR[Title/Abstract] '
    'AND (gyrB[Title/Abstract] OR "gyrase B"[Title/Abstract] '
    'OR "fluoroquinolone resistance"[Title/Abstract])',
)
check(
    "TB + CRISPR + pncA (expanded)",
    '("Mycobacterium tuberculosis"[Title/Abstract]) AND CRISPR[Title/Abstract] '
    'AND (pncA[Title/Abstract] OR pyrazinamidase[Title/Abstract] '
    'OR "pyrazinamide resistance"[Title/Abstract])',
)
check(
    "TB + CRISPR + CFP-10 (expanded)",
    '("Mycobacterium tuberculosis"[Title/Abstract]) AND CRISPR[Title/Abstract] '
    'AND ("CFP-10"[Title/Abstract] OR Rv3874[Title/Abstract] '
    'OR "ESAT-6/CFP-10"[Title/Abstract])',
)
check(
    "TB + CRISPR + rrs (expanded)",
    '("Mycobacterium tuberculosis"[Title/Abstract]) AND CRISPR[Title/Abstract] '
    'AND (rrs[Title/Abstract] OR "16S rRNA"[Title/Abstract] '
    'OR "aminoglycoside resistance"[Title/Abstract])',
)
check(
    "TB + CRISPR + eis (expanded)",
    '("Mycobacterium tuberculosis"[Title/Abstract]) AND CRISPR[Title/Abstract] '
    'AND (eis[Title/Abstract] OR "enhanced intracellular survival"[Title/Abstract] '
    'OR "kanamycin resistance"[Title/Abstract])',
)

print("\n=== CHOLERA GAPS ===")
check(
    "Cholera + CRISPR + diagnostic",
    '("Vibrio cholerae"[Title/Abstract]) AND CRISPR[Title/Abstract] '
    'AND (diagnostic[Title/Abstract])',
)
check(
    "Cholera + CRISPR + diagnostic (broadened)",
    '("Vibrio cholerae"[Title/Abstract] OR cholera[Title/Abstract]) AND '
    'CRISPR[Title/Abstract] AND (diagnostic[Title/Abstract] OR detection[Title/Abstract])',
)
check(
    "Cholera + CRISPR + ctxA",
    '("Vibrio cholerae"[Title/Abstract]) AND CRISPR[Title/Abstract] '
    'AND (ctxA[Title/Abstract] OR "cholera toxin"[Title/Abstract])',
)

print("\n=== OTHER HIGH-RISK GAPS ===")
check(
    "Ebola + CRISPR + NP/nucleoprotein",
    '("Ebola"[Title/Abstract]) AND CRISPR[Title/Abstract] '
    'AND (NP[Title/Abstract] OR nucleoprotein[Title/Abstract])',
)
check(
    "RSV + CRISPR + N protein",
    '("RSV"[Title/Abstract]) AND CRISPR[Title/Abstract] '
    'AND ("N protein"[Title/Abstract] OR nucleoprotein[Title/Abstract])',
)
check(
    "RSV + CRISPR + G protein",
    '("RSV"[Title/Abstract]) AND CRISPR[Title/Abstract] '
    'AND ("G protein"[Title/Abstract] OR "attachment protein"[Title/Abstract])',
)

print("\nDone. Any 'FALSE GAP' above must be corrected in the papers.")
