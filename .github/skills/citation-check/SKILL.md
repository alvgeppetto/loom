# Citation Integrity Check Skill

## When to Use

**MANDATORY** before building any paper PDF. Use whenever:
- A paper markdown file is edited (any change to text, references, or tables)
- New in-text citations `[N]` are added
- The references section is modified
- Running `scripts/build-preprints.sh`

## What It Checks

`scripts/verify-citations.py` enforces five rules:

| Rule | What it catches |
|------|----------------|
| **R1 — No placeholders** | `[?]`, `[??]` etc. anywhere in body text |
| **R2 — Cite→Ref** | Every `[N]` in text has a matching `[N]` entry in References |
| **R3 — Ref→Cite** | Every reference entry is cited at least once (warns; not fatal) |
| **R4 — DOI/PMID** | Every reference has `DOI: 10.xxx/yyy` or `PMID: NNNNNNN` (exception: "in preparation" companion manuscripts) |
| **R5 — DOI resolves** | Live HEAD request to `doi.org` confirms each DOI exists (skipped with `--offline`) |

## Paper Reference Format

References must follow this exact format for R4 to pass:

```
[N] Author A, Author B, et al. Title of paper. *Journal*. YEAR;VOL(ISSUE):PAGES.
    DOI: 10.XXXX/XXXXXX
```

For unpublished companion manuscripts:
```
[N] Author. Title. Companion manuscript, YEAR. (In preparation — available as preprint alongside this work.)
```
No DOI required for entries containing "In preparation".

## Running Manually

```bash
# Structural check (fast, no network)
python scripts/verify-citations.py --offline publications/papers/pangenomic-crispr-targets-merged.md

# Full check including DOI resolution (requires internet)
python scripts/verify-citations.py publications/papers/pangenomic-crispr-targets-merged.md

# Check all papers
python scripts/verify-citations.py publications/papers/*.md

# Strict mode (warnings become failures)
python scripts/verify-citations.py --strict publications/papers/*.md
```

## Build Integration

`scripts/build-preprints.sh` automatically runs the `--offline` structural check before any build. To bypass (not recommended):
```bash
SKIP_CITE_CHECK=1 ./scripts/build-preprints.sh
```

## Adding a New Citation

When adding a new in-text citation:
1. Add `[N]` in body text where appropriate
2. Add corresponding `[N] Author...` entry to the References section
3. Include `DOI: 10.xxxx/xxxxx` on a continuation line
4. Run `python scripts/verify-citations.py --offline <paper>` to confirm no errors
5. Optionally run live check to confirm DOI resolves

## Common Failures and Fixes

| Error | Fix |
|-------|-----|
| `Placeholder citation: '[?]'` | Replace with actual `[N]` and add reference entry |
| `Cited in text but no reference entry` | Add `[N]` entry to References section |
| `Missing DOI/PMID` | Look up DOI on doi.org or PMID on pubmed.ncbi.nlm.nih.gov and add to reference |
| `DOI did NOT resolve` | Check DOI for typos; verify on doi.org manually |
