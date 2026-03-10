# Incomplete-Analysis Sections Undermine Papers

Date: 2026-03-07
Tags: [paper, review, patent, completeness, S7]

## Context

The S7 patent landscape section included Table S7c — an exact-match patent search for
10 out of 52 guide sequences. This was the only incomplete analysis in the paper.

## What Happened

The incompleteness (10/52) required a caveat that was so heavy it undermined the entire
section. A reviewer would rightly ask "why only 10?" and the answer (no API, manual
process) is not publishable. The section took up space without landing a punch.

## The Lesson

1. **Never publish partial computational results when the full set is feasible.** If you
   can't compute all 52, don't publish 10 with a caveat — it reads as sloppy.
2. **If the full computation is infeasible, restructure the argument.** Instead of the
   incomplete exact-match search, we reframed around data we *do* have completely:
   sequence novelty (≥8 mm from all 83 published guides) and gene-level patent density.
3. **Caveats should be one sentence, not a paragraph.** If a caveat needs a paragraph,
   the analysis it's hedging probably shouldn't be in the paper.
4. **API availability check first.** Before building a data pipeline, verify that the
   source has a programmatic API. Google Patents doesn't — this should have been caught
   before the 10-guide manual search was started.

## Evidence

- Original S7c: 10/52 guides tested, 0 matches, multi-sentence caveat
- Fix: Cut S7c, kept S7a/S7b landscape tables, strengthened argument from complete data
- Commit: cc602d9
