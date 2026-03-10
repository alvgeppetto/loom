---
name: medical-claims-verification
description: "MANDATORY for any code or document that makes medical, scientific, or literature-based claims. Enforces verification discipline for PubMed queries, gene naming, gap claims, and cross-database limitations. Use when: writing/editing PubMed scan code, updating papers/briefs with quantitative claims, adding new pathogens or genes, or reviewing gap analysis output."
argument-hint: "Describe the claim or query being verified"
---

# Medical Claims Verification Skill

**This skill exists because sloppy query construction in March 2026 produced false
medical claims that nearly destroyed project credibility.** Every rule below was
earned through a real bug that generated a false claim in a medical context.

## MANDATORY: Read Before Any Medical/Literature Work

### The Three Bugs That Almost Killed the Project

1. **Mpox OR bug (2026-03-04):** `"Mpox OR Monkeypox"` was treated as a literal
   phrase in PubMed, not boolean OR. Result: 0 papers found when 69 existed.
   17 false gap claims published. **Root cause:** `_name_clause()` didn't exist;
   raw string was wrapped in quotes.

2. **Gene synonym bug (2026-03-04):** Gene queries used only primary names
   (e.g., "inhA"). PubMed papers use expanded names ("enoyl-ACP reductase").
   Result: TB showed 6 gene gaps when only 2 were real. **Root cause:** No
   synonym expansion; narrow queries = false negatives = false gap claims.

3. **Duplicate gene inflation (2026-03-04):** Ebola listed both "NP" and
   "nucleoprotein" as separate genes (same gene). RSV listed "N protein" +
   "nucleoprotein", "F protein" + "fusion protein", "G protein" + "attachment
   protein". MERS listed "spike" + "S protein". Result: Gap counts inflated
   by 7 (53→46 when fixed). **Root cause:** No deduplication check; copy-paste
   gene lists from different sources without reconciliation.

---

## Rules

### Rule 1: Every PubMed Query MUST Be Tested

Before any PubMed scan code is merged, the query builder MUST have tests covering:

- [ ] `_name_clause()` correctly handles single names, multi-word names, and `OR`-separated synonyms
- [ ] Gene synonym expansion produces valid PubMed boolean (not quoted phrases)
- [ ] No duplicate genes per pathogen (test: `len(genes) == len(set(genes))`)
- [ ] No synonym pairs listed as separate genes (test: known synonym pairs must not coexist)
- [ ] No non-pathogen entries in pathogen lists (e.g., Human GRCh38)

**Test file:** `tests/test_pubmed_queries.py` — run with `python3 -m pytest tests/test_pubmed_queries.py -v`

### Rule 2: Never Claim "Zero Research Exists"

The correct phrasing is always:

> "Zero results returned by our PubMed query [exact query string] as of [date]"

NEVER say:
- "No research exists"
- "Completely unexplored"
- "No one has studied this"
- "This has never been done"

**Why:** PubMed is one database. It misses:
- Chinese databases (CNKI, Wanfang, VIP)
- Russian databases (eLIBRARY.ru, CyberLeninka)
- Arabic/Persian databases
- Preprints not indexed in PubMed (bioRxiv, medRxiv partial coverage)
- Patents (Google Patents, Espacenet)
- Conference proceedings not indexed
- Dissertations and theses
- Papers with different terminology than our query

### Rule 3: Every New Gene MUST Have Synonyms Checked

When adding a gene to `PATHOGEN_GENES`:

1. Search UniProt for the gene name → note all synonyms
2. Search PubMed manually with the primary name AND each synonym
3. Add all synonyms that return different results to `GENE_SYNONYMS`
4. Verify the gene is not a synonym of an already-listed gene

**Synonym trap examples:**
- NP = nucleoprotein (same protein)
- S protein = spike = spike protein (same protein)
- F protein = fusion protein (same protein)
- G protein = attachment protein (RSV context)
- inhA = enoyl-ACP reductase = Rv1484 (same enzyme)

### Rule 4: Every Scan Result MUST Be Spot-Checked

After running `pubmed_scan.py`, manually verify at least:
- 3 claims that return zero (are they genuinely zero, or query bugs?)
- 3 claims that return low counts (1-5; could be synonym issues)
- All newly-added pathogens or genes (first scan of a new entity)

Use `scripts/spot_check_claims.py` or manual PubMed queries.

### Rule 5: Numbers Flow One Direction

The single source of truth for gap counts is `web/data/pubmed-scan-results.json`.
When the scan runs, all documents MUST be updated from this file:

```
pubmed_scan.py → pubmed-scan-results.json → all papers/briefs
```

NEVER manually edit a number in a paper. Always re-run the scan and propagate.

After updating, grep for stale numbers:
```bash
grep -rn "OLD_NUMBER" publications/ docs/ --include="*.md"
```

### Rule 6: Pathogens ≠ Reference Genomes

Human GRCh38 is in the CRISPR target database for off-target screening.
It is NOT a pathogen. It MUST be excluded from gap analysis.
The scan script MUST filter out non-pathogen entries before generating claims.

### Rule 7: Document Non-English Database Limitations

Every paper MUST contain a limitation statement about
PubMed's English-language bias. Template:

> **Limitation: English-language database bias.** Our literature analysis queries
> PubMed exclusively, which has limited coverage of research published in Chinese
> (CNKI, Wanfang), Russian (eLIBRARY.ru), Arabic, and other non-English databases.
> Gap claims reflect PubMed coverage as of [date], not the totality of global
> CRISPR research. Researchers in non-English-speaking countries may have published
> relevant work not captured by our scan.

---

## Verification Checklist

Use this checklist before publishing any claim:

```
□ Query string logged with exact PubMed syntax
□ Date of query recorded
□ OR-separated names produce boolean, not phrases
□ Gene synonyms expanded (checked against UniProt)
□ No duplicate genes in gene list
□ Non-pathogen entries excluded
□ At least 3 zero-result claims manually spot-checked
□ Phrasing says "our query returned zero" not "no research exists"
□ Non-English database limitation documented
□ All documents updated from single source (JSON)
□ grep confirms no stale numbers remain
□ Tests pass: pytest tests/test_pubmed_queries.py
```

---

## Recovery Procedure

If a false claim is discovered:

1. **Stop.** Do not publish, present, or share anything until fixed.
2. Run `scripts/spot_check_claims.py` on the specific claim
3. If confirmed false: fix the query bug in `pubmed_scan.py`
4. Add a regression test to `tests/test_pubmed_queries.py`
5. Re-run full scan: `python3 scripts/pubmed_scan.py`
6. Update ALL documents (grep for old numbers)
7. Rebuild ALL PDFs: `bash scripts/build-preprints.sh`
8. Redeploy: `bash scripts/deploy-azure.sh --site-only`
9. Record lesson in `docs/lessons/`
10. Commit with `fix:` prefix explaining the false claim
