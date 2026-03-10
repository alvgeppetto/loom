#!/usr/bin/env python3
"""
verify-citations.py — Citation integrity checker for LOOM papers.

Checks:
  1. Every [N] in-text citation has a matching [N] reference entry.
  2. Every reference entry is cited at least once.
  3. No [?] placeholder citations remain.
  4. Every reference entry carries a DOI or PMID.
  5. DOIs resolve via doi.org HEAD request (skipped with --offline).

Usage:
  python scripts/verify-citations.py publications/papers/pangenomic-crispr-targets-merged.md
  python scripts/verify-citations.py --offline publications/papers/*.md
"""
import re
import sys
import argparse
import urllib.request
import urllib.error
from pathlib import Path


PLACEHOLDER_RE = re.compile(r'\[\?+\]')
CITE_RE = re.compile(r'\[(\d+(?:,\s*\d+)*)\]')   # [1] or [1,2] or [1, 2]
REF_ENTRY_RE = re.compile(r'^\[(\d+)\]\s+(.+)', re.MULTILINE)
DOI_IN_REF_RE = re.compile(r'DOI:\s*10\.\S+', re.IGNORECASE)
PMID_IN_REF_RE = re.compile(r'PMID:\s*\d+', re.IGNORECASE)
DOI_BARE_RE = re.compile(r'10\.\d{4,}/\S+')


def extract_cites(text: str) -> set[int]:
    """Return every cited reference number from text (excluding reference section)."""
    # Split off the References section to avoid counting ref entries as cites
    ref_section_start = re.search(r'^##\s+References\s*$', text, re.MULTILINE)
    body = text[:ref_section_start.start()] if ref_section_start else text
    nums: set[int] = set()
    for m in CITE_RE.finditer(body):
        for part in m.group(1).split(','):
            nums.add(int(part.strip()))
    return nums


def extract_references(text: str) -> dict[int, str]:
    """Return {num: full_ref_text} for all [N] entries in References section.
    Handles multi-line references where continuation lines start with whitespace.
    """
    ref_section_match = re.search(r'^##\s+References\s*$', text, re.MULTILINE)
    if not ref_section_match:
        return {}
    ref_text = text[ref_section_match.start():]
    # Split into lines; aggregate continuation lines (starting with spaces/tabs)
    refs: dict[int, str] = {}
    current_num: int | None = None
    current_lines: list[str] = []

    for line in ref_text.splitlines():
        m = REF_ENTRY_RE.match(line)
        if m:
            # Save previous entry
            if current_num is not None:
                refs[current_num] = ' '.join(current_lines)
            current_num = int(m.group(1))
            current_lines = [m.group(2).strip()]
        elif current_num is not None and line.startswith((' ', '\t')) and line.strip():
            current_lines.append(line.strip())
        elif current_num is not None and line.startswith('---'):
            # End of references section
            break

    if current_num is not None:
        refs[current_num] = ' '.join(current_lines)

    return refs


def extract_dois(ref_text: str) -> list[str]:
    """Extract DOI strings from a reference entry."""
    dois = []
    for m in DOI_IN_REF_RE.finditer(ref_text):
        raw = m.group(0).replace('DOI:', '').replace('doi:', '').strip()
        dois.append(raw)
    if not dois:
        for m in DOI_BARE_RE.finditer(ref_text):
            dois.append(m.group(0).rstrip('.,)'))
    return dois


def doi_resolves(doi: str, timeout: int = 8) -> tuple[bool, str]:
    """Check if a DOI resolves via doi.org HEAD request."""
    url = f"https://doi.org/{doi}"
    req = urllib.request.Request(url, method='HEAD',
                                  headers={'User-Agent': 'LOOM-citation-checker/1.0'})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return True, str(resp.status)
    except urllib.error.HTTPError as e:
        # 302/301 redirects from doi.org are normal success
        if e.code in (301, 302, 303):
            return True, f"redirect {e.code}"
        return False, f"HTTP {e.code}"
    except Exception as e:
        return False, str(e)


def check_paper(path: Path, offline: bool = False, strict: bool = False) -> list[str]:
    """Check a single paper. Returns list of error/warning strings (empty = pass)."""
    text = path.read_text(encoding='utf-8')
    issues = []

    # 1. Placeholder citations
    for m in PLACEHOLDER_RE.finditer(text):
        lno = text[:m.start()].count('\n') + 1
        issues.append(f"  FAIL [line {lno}] Placeholder citation: {m.group(0)!r}")

    cites = extract_cites(text)
    refs = extract_references(text)

    if not refs:
        issues.append("  WARN  No References section found — skipping structural checks")
        return issues

    # 2. Cited but missing reference
    for n in sorted(cites):
        if n not in refs:
            issues.append(f"  FAIL  [ref {n}] Cited in text but no reference entry")

    # 3. Reference defined but never cited
    for n in sorted(refs):
        if n not in cites:
            issues.append(f"  WARN  [ref {n}] Reference entry exists but never cited")

    # 4. References missing DOI or PMID
    for n, rtext in sorted(refs.items()):
        has_doi = bool(DOI_IN_REF_RE.search(rtext) or DOI_BARE_RE.search(rtext))
        has_pmid = bool(PMID_IN_REF_RE.search(rtext))
        if not has_doi and not has_pmid:
            if "In preparation" in rtext or "in preparation" in rtext or "companion" in rtext.lower():
                pass  # unpublished companion manuscript — no DOI expected
            elif "WHO" in rtext or "World Health Organization" in rtext or "CDC" in rtext:
                pass  # institutional report — no DOI assigned
            else:
                issues.append(f"  FAIL  [ref {n}] Missing DOI/PMID: {rtext[:80]}…")

    # 5. DOI resolution (live check)
    if not offline:
        for n, rtext in sorted(refs.items()):
            dois = extract_dois(rtext)
            for doi in dois:
                ok, code = doi_resolves(doi)
                if ok:
                    issues.append(f"  OK    [ref {n}] DOI resolves ({code}): {doi}")
                else:
                    issues.append(f"  FAIL  [ref {n}] DOI did NOT resolve ({code}): {doi}")

    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description="Citation integrity checker for LOOM papers")
    parser.add_argument('papers', nargs='+', type=Path, help='Markdown paper(s) to check')
    parser.add_argument('--offline', action='store_true',
                        help='Skip live DOI resolution (structural checks only)')
    parser.add_argument('--strict', action='store_true',
                        help='Treat WARNs as failures')
    args = parser.parse_args()

    total_failures = 0
    for path in args.papers:
        if not path.exists():
            print(f"\n[SKIP] {path} — file not found")
            continue

        print(f"\n{'='*60}")
        print(f"Checking: {path}")
        print('='*60)
        issues = check_paper(path, offline=args.offline, strict=args.strict)

        if not issues:
            print("  PASS  All citations verified.")
        else:
            for line in issues:
                print(line)

        fails = sum(1 for i in issues if i.strip().startswith('FAIL'))
        warns = sum(1 for i in issues if i.strip().startswith('WARN'))
        oks   = sum(1 for i in issues if i.strip().startswith('OK'))
        if args.strict:
            fails += warns
        total_failures += fails

        print(f"\n  Summary: {fails} failures, {warns} warnings, {oks} DOIs verified")

    print(f"\n{'='*60}")
    if total_failures > 0:
        print(f"CITATION CHECK FAILED: {total_failures} issue(s) require attention.")
        return 1
    else:
        print("CITATION CHECK PASSED.")
        return 0


if __name__ == '__main__':
    sys.exit(main())
