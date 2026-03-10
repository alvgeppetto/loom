# Pangenomic CRISPR Diagnostic Target Discovery: 52 Novel Sites in SARS-CoV-2 and a Multi-Pathogen Conserved-Target Database with Verified Literature Gap Analysis

**Alvaro Videla Godoy**^1^

^1^ Independent Researcher

**Correspondence:** videlalvaro@gmail.com

**Preprint Draft — March 2026**

---

## Abstract

**Background:** CRISPR-based diagnostics overwhelmingly target
immune-exposed structural genes (Spike, Nucleocapsid, Envelope), which accumulate
escape mutations that degrade assay sensitivity over time. Essential replication
machinery under strong purifying selection offers a more conserved alternative,
but systematic cross-pathogen surveys of such regions are lacking.

**Methods:** We scanned 9,193,298 SARS-CoV-2 genomes (262 GB) and 11 additional
pathogen collections using LOOM, an FM-index k-mer engine, complemented by
pigeonhole-Hamming off-target screening against seven host genomes. A four-strategy
literature scanner with NCBI ontology-derived synonym expansion classified gene-level
literature coverage.

**Results:** We report 52 conservation-ranked CRISPR target candidates in SARS-CoV-2
replication machinery (24 SpCas9, 28 Cas12a; 41 above 92% conservation). Conservation is
necessary but not sufficient: 29/52 fall outside the 40–70% GC window; all require wet-lab
validation. The top 8 (98.9–99.2%) are within 1 percentage point of the best published
diagnostic guide (DETECTR_E: 98.34%) on an identical corpus. Across 12 pathogens, a
~130,000-target database is provided (sampling
depths vary: 382 cholera to 9.2M SARS-CoV-2 genomes). Literature verification confirmed 6
gene gaps (5 Mpox, 1 cholera); for SARS-CoV-2, no confirmed gaps exist, but nsp13/nsp14 have
non-diagnostic CRISPR research with no guide designs — we cannot distinguish "unexplored" from
"explored and abandoned."

**Conclusions:** The primary contributions are (1) 52 conservation-ranked SARS-CoV-2
candidates requiring experimental validation, (2) a ~130,000-target multi-pathogen
database for guide prioritization, and (3) the demonstration that single-query gap
analysis overestimates gaps by 5–8×. Data, code, and indexes are publicly available.

**Keywords:** CRISPR diagnostics, pangenomic analysis, literature gap analysis, pathogen
genomics, open science, replication machinery, ORF1b

---

## 1. Introduction

CRISPR-based nucleic acid detection (SHERLOCK [1], DETECTR [2], PACMAN [3]) has demonstrated
rapid pathogen detection, but guide RNA design has followed historical RT-PCR conventions:
structural and surface genes (Spike, Nucleocapsid, Envelope) predominate. Two problems follow
from this pattern. First, immune-exposed structural genes accumulate mutations under selection,
degrading diagnostic sensitivity — as seen when D614G swept the SARS-CoV-2 Spike and
invalidated SHERLOCK_S [4]. Second, essential replication enzymes under strong purifying
selection remain underexplored as diagnostic targets — our gap analysis finds them UNCERTAIN
or FALSE (non-diagnostic work exists), not strictly unexplored.

We address both problems. For SARS-CoV-2, we identify 52 diagnostic-novel CRISPR target
sites in replication machinery — spanning ORF1a and ORF1b (nsp13/Helicase, nsp14/ExoN)
plus accessory ORFs — validated at full pangenomic scale (9,193,298 sequences). For 11 additional
pathogens, we provide a conservation-ranked target database and a multi-source
literature gap analysis. These two contributions are linked: the pangenomic pipeline
developed for deep SARS-CoV-2 scanning generalizes to any pathogen, and the multi-pathogen
database validates that the conservation-ranking and off-target methodology scales across
diverse genome architectures (RNA/DNA viruses, bacteria). The literature analysis
covers all 12 pathogens uniformly; for SARS-CoV-2 specifically, it finds no outright
gene-level gaps but reveals that existing CRISPR work on replication genes (nsp13, nsp14)
is exclusively non-diagnostic — an unexploited opportunity rather than a blank space.

A previous version of this work overclaimed 43 gaps using single-query methodology; this
paper reports only the 6 gene-level gap claims verified by four independent evidence
strategies including ontology-derived synonym expansion.

---

## 2. Materials and Methods

### 2.1 SARS-CoV-2 Pangenomic Scan

**Tool:** All scanning, conservation ranking, off-target screening, novelty checking, and
literature gap analysis were performed using LOOM, an FM-index-based exact-match k-mer
engine developed for this work (source: https://github.com/alvgeppetto/loom). LOOM compiles
to a 335 KB WebAssembly binary for in-browser use and a native Rust CLI for large-scale
pipeline runs.

**Data:** 9,193,298 SARS-CoV-2 sequences (262 GB); Wuhan-Hu-1 reference (NC_045512.2,
29,903 bp); NCBI RefSeq Viral (19,019 complete viral genomes, panviral reference).

**Pipeline:**
1. Enumerate all 29,640 overlapping 20-mers from the Wuhan-Hu-1 reference.
   (20-mer length matches the SpCas9 spacer; Cas12a guides add a 3-nt PAM-proximal
   extension to yield 23-mers at the database stage — see §2.2.)
2. Query each 20-mer (forward + RC) against the panviral database — score cross-species
   conservation by exact-match count across 19,019 genomes.
3. Validate candidates conserved across ≥2 viral species against all 9,193,298 SARS-CoV-2
   genomes by exact-match sliding-window scan.
4. Filter for PAM compatibility (SpCas9 NGG or Cas12a TTTV).
5. Exclude candidates overlapping any published CRISPR diagnostic guide region
   (N gene, E gene, S gene, nsp12/RdRp — covering DETECTR, PACMAN, SHERLOCK, CARVER targets).
6. Verify sequence-level novelty: compare all 52 remaining candidate spacer sequences
   against a curated database of 83 published CRISPR guide sequences from 21 sources
   and therapeutic papers (compiled via `compile_published_guides.py`) using 2-bit Hamming
   distance at ≤3 mismatches, checking both strands (`loom guide-novelty`).

**Same-corpus validation (Tables 2–3, S1):** Conservation percentages for
both novel targets and published guide benchmarks were independently validated
using single-pass Aho-Corasick exact-match streaming (`loom guide-conservation`)
over the complete 262 GB FASTA. Both guide sets were measured against an
identical denominator of 9,193,298 genomes, eliminating cross-tool biases
present in earlier cross-dataset comparisons.

**Off-target screening:** All 89 PAM-compatible candidates screened against GRCh38
(3,209,286,105 bp) via exact 20-mer matching (forward + RC) using the LOOM FM-index,
followed by a pigeonhole-seeded Hamming-distance scan at ≤3 mismatches — an approach
comparable to Cas-OFFinder [5] but integrated into the same index used for conservation
ranking.

**Tool validation.** LOOM's off-target detection uses pigeonhole-seeded Hamming distance:
for ≤k mismatches, the 20-nt guide is split into k+1 seeds such that at least one must
match exactly (pigeonhole principle). For k=3, we use 4×5-nt seeds indexed via the FM-index
for O(m log n) lookup, with full Hamming verification of candidate positions. This is
algorithmically equivalent to Cas-OFFinder's seed-extend approach [5].

*Algorithmic verification:* The FM-index/BWT construction was stress-validated against
BWA 0.7.19 on a normalized 19,048,983-bp Zika corpus using 2,000 random 23-mers
(1,894 unique). LOOM exact counts were also cross-checked against exhaustive overlapping
substring counts on the same corpus. Both checks were fully concordant (0 mismatches;
artifact: `data/crispr_guides/fm_index_verification.json`).

*Functional validation:* Off-target search accuracy was validated against the EMX1 guide
(GAGTCCGAGCAGAAGAAGAA), a canonical benchmark with GUIDE-seq validated off-target sites
[Tsai et al. 2015, Nature Biotechnology 33:187-197]. All 9 known off-target sites (≤3mm)
were correctly identified with concordant mismatch counts (`emx1_validation.json`). **[C]**

### 2.2 Multi-Pathogen Target Database

Genome collections for 11 human pathogens were downloaded from NCBI (GenBank, RefSeq,
Pathogen Detection) between February and March 2026. Table 1 summarizes collections.

For each pathogen, a sliding-window scanner extracted all overlapping 23-mers
(23-nt length accommodates both SpCas9 20-nt spacers and Cas12a 23-nt spacers;
the longer window ensures Cas12a-compatible candidates are captured at full length),
counted occurrence across the corpus (conservation proxy), and ranked by frequency. The top 10,000
per pathogen were retained and classified by PAM type (SpCas9 NGG, Cas12a TTTV).

### 2.3 Animal Host Off-Target Screening

All 52 novel guide sequences were screened against seven host reference genomes using
the same pigeonhole-seeded Hamming-distance scan described in §2.1 (≤3 mismatches,
forward + reverse complement). The seven genomes are: Human (GRCh38, 3.2 Gbp),
Pig (Sscrofa11.1, 2.4 Gbp), Egyptian Fruit Bat (mRouAeg1, 1.8 Gbp),
Chicken (GRCg7b, 1.0 Gbp), Cow (ARS-UCD1.3, 2.6 Gbp), Dromedary Camel (CamDro3,
2.0 Gbp), and Mouse (GRCm39, 2.6 Gbp). Species were selected for relevance to
coronavirus reservoirs (bat, camel), zoonotic transmission (pig, cow), agricultural
surveillance (chicken), and laboratory models (mouse). Results are reported in
Table S8.

**Pan-pathogen exact-match cross-reactivity screen.** As a complementary analysis, the
top 100 NGG-PAM guides per pathogen (1,100 total across 11 pathogens) and, separately,
the top 100 TTTV-PAM (Cas12a) guides per pathogen (972 unique — tuberculosis contributed
27 and dengue 76 and zika 69, all others 100; only valid TTTV PAMs were included, excluding
TTTT) were screened for exact matches in all seven host genome FM-indexes. This broader
screen tests whether highly conserved pathogen guides share exact 23-mer identity with any
host genome, which would indicate potential off-target amplification in a diagnostic
setting. Results are reported in §3.4 and Tables S9a–b.

### 2.4 Multi-Source Literature Gap Analysis (v4 Methodology)

**Why single-query scanning fails:** Our earlier v2 analysis (archived) used a single
`term[Title/Abstract]` PubMed query per gene. This misclassifies well-studied genes as
unexplored when: (a) papers use synonyms not in the query, (b) the trait is described in
full text only, or (c) the work is indexed in non-PubMed databases. Cholera diagnostics
was the clearest example: a search for `"ctxA" AND CRISPR` returned zero results, but
57 landscape papers exist and the Europe PMC full-text search returned hundreds of
cholera+CRISPR co-mentions. The v2 method overcounted gaps by 5-8×.

**Ontology-derived synonym expansion (new in v4):** The v4 scanner loads 8,283 NCBI gene
annotations from six ontology sources (NCBI Taxonomy, Disease Ontology, MONDO, WHO GHO,
ICD-11) via `web/data/ontology-enrichment.json`. Each annotated gene's aliases and full
names are matched to pathogen-specific gene lists and merged with manually curated synonym
tables. This generated 51 additional synonyms across 12 pathogens that were invisible to
manual curation. The impact was immediate: RSV SH protein, classified CONFIRMED under v3
(manual synonyms only), was reclassified UNCERTAIN after ontology-derived synonyms
surfaced 278 Europe PMC entries and 1 corpus title match.

**v4 strategy — four evidence layers per claim:**

**Strategy 1 — PubMed narrow:** `term[Title/Abstract] AND CRISPR[Title/Abstract]`.
Equivalent to v2. Used as a *necessary* but not *sufficient* condition for a gap.

**Strategy 2 — PubMed broad:** Full synonym expansion OR-clause in `[Title/Abstract]`.
Synonyms are combined from manual curation and ontology-derived gene annotations
(e.g., hapA → "hemagglutinin protease", "HA/P"; SH protein → "small hydrophobic protein",
"SH", plus ontology aliases). A hit here falsifies a gap.

**Strategy 3 — Europe PMC corpus scan:** The full landscape abstract corpus for each
pathogen (all CRISPR × pathogen papers, up to 500 abstract records) was downloaded,
cached locally (`target/pubmed_cache_v4/corpus_*.json`), and scanned for co-occurrence
of gene/application terms against the ontology-expanded synonym set. A co-occurrence hit
in an independent database corpus falsifies a gap.

**Strategy 4 — Ontology synonym cross-check:** All gene/synonym matches are tagged with
their source (`manual`, `ontology`, `manual+ontology`, `none`) for full provenance
traceability. This ensures no gap claim survives because of incomplete nomenclature.

**Confidence classification** (exact rules in `src/dna/pubmed_scan.rs`):
- `CONFIRMED` — PubMed broad = 0 AND corpus co-occurrence = 0 AND EuropePMC ≤ 30. Genuine gap candidate.
  The ≤30 threshold was chosen conservatively to minimize false CONFIRMED classifications;
  the RSV SH reclassification (v3→v4) demonstrated that synonym gaps can inflate EuropePMC counts
  by an order of magnitude, so we set a low bar that only genes with negligible full-text presence pass.
- `UNCERTAIN` — PubMed broad > 0, OR corpus co-occurrence found, OR EuropePMC > 30.
  Requires manual review.

**Threshold sensitivity analysis.** The ≤30 EuropePMC threshold was stress-tested across
multiple values (10, 20, 30, 40, 50, 100). The 6 gene-level CONFIRMED gaps are **stable
across thresholds 20–30** because all have EuropePMC counts well below 30 (range: 4–17;
hapA=6, A33R=15, A56R=4, B5R=17, E8L=12, J2R=8). Boundary cases at 32–36 (cholera ompT,
ompU, rtxA) are appropriately classified UNCERTAIN — they have higher literature presence
than the CONFIRMED gaps. The threshold is not cherry-picked: lowering to 20 yields the
same 6 gaps; raising to 40 adds genes we deliberately flagged for manual review.
Full sensitivity sweep in `gap_threshold_sensitivity.json`. **[C]**
- `PROBABLE` — PubMed narrow > 0 but broad = 0 (term appears narrowly in title/abstract).
- `FALSE` — PubMed broad > 5. Retracted as a *literature gap*. **Important: FALSE
  classification answers "is there CRISPR literature?" not "is there diagnostic
  potential."** A FALSE gene may represent an unexploited diagnostic opportunity if
  all existing CRISPR literature is therapeutic or mechanistic rather than diagnostic
  (see nsp13/nsp14 discussion in §3.1 and §4.2). We retain the FALSE label because it
  accurately reflects our detection criterion (literature exists), while acknowledging
  the category creates interpretive tension with our opportunity framing.

All query strings, timestamps, and per-strategy counts are recorded in
the v4 scan artifacts (see §5). Scanner source: `loom pubmed-scan`
(Rust CLI).

**Audit date:** March 2026.

### 2.5 Version History and Methodology Self-Correction

An earlier version of this work (archived at
`publications/papers/pan-pathogen-crispr-targets.md`) claimed 43 literature gaps using a
single narrow PubMed query per gene. A v3 re-verification using three evidence strategies
retracted 34 of those 43; 9 survived. A subsequent v4 scanner (Rust CLI, with NCBI
ontology-derived synonym expansion) reclassified RSV SH protein from CONFIRMED to
UNCERTAIN — the ontology synonyms surfaced 278 Europe PMC entries invisible to the manual
synonym list. **6 confirmed gene-level gaps remain** (5 Mpox genes, 1 cholera gene).
The reduction from 9 to 6 reflects one reclassification (RSV SH protein → UNCERTAIN)
and a change from query-level to gene-level counting: the v3 tally of 9 counted
J2R and thymidine kinase as separate entries (both are OPG101) and A56R and hemagglutinin
as separate entries (both are OPG185); the v4 tally correctly deduplicates to 5 unique
Mpox genes plus hapA. The
present paper reports only claims verified by four evidence strategies including
ontology-expanded synonyms. The archived version is kept as a documented lesson in why
single-source gap analysis is unreliable.

**Same-corpus revalidation (March 2026):** Conservation percentages in this version were
recomputed using single-pass Aho-Corasick exact-match streaming (`loom guide-conservation`)
over the same FASTA and with the same tool for both novel and published guide sets.
All reported percentages in this manuscript use a single denominator of 9,193,298 genomes.
Initial FM-index discovery was run against a 9,419,528-sequence download (February 2026);
same-corpus revalidation re-scanned a separately downloaded 9,193,298-sequence FASTA
(March 2026, identical tool and parameters). The earlier artifact
(`novel_targets_52.json`, denominator 9,419,528) is retained for provenance; all
conservation percentages reported in this manuscript derive from the same-corpus
revalidation (`novel_targets_same_corpus.json`, denominator 9,193,298).

**Development infrastructure for iterative self-correction.** The v2→v3→v4 corrections
were only detectable because the entire pipeline—code, data artifacts, and manuscript—
was version-controlled (Git) throughout development. Each claimed number in the manuscript
is validated by `verify-paper-claims.py`, which compares paper text against JSON artifacts;
this script fails if any number is changed without updating its source artifact (or vice
versa). The integrated development environment (VS Code with AI chat) enabled simultaneous
CLI pipeline execution, manuscript editing, and verification—catching regressions that
would otherwise occur when fixing one section breaks consistency with another. This
infrastructure is not incidental: a researcher without diff-based regression detection
could not have reliably caught the 37 retracted gaps or the RSV SH reclassification.
The development log (Git history) is part of the reproducibility record.

---

## 3. Results

### 3.1 SARS-CoV-2: 52 Novel Diagnostic Targets in Replication Machinery

From 29,640 windows (29,884 theoretical 20-mers from the 29,903 bp Wuhan-Hu-1 reference;
244 excluded due to ambiguous-base (N) content), 1,503 showed conservation across ≥2 panviral species. After PAM
filtering, 183 windows in 89 distinct genomic regions carry compatible PAM sites.
Removing regions overlapping published diagnostic guide gene intervals (N, E, S, nsp12/RdRp)
leaves **52 truly novel regions**. Novelty is confirmed at two levels:
(1) gene-interval level — these 52 regions do not overlap the genomic coordinate intervals
of published diagnostic guide targets (N, E, S, nsp12/RdRp); and (2) sequence level — all
52 candidate spacer sequences were compared against 83 published CRISPR guide sequences
compiled from 21 sources (15 hand-curated papers plus 6 Europe PMC full-text extractions;
Table S4) using 2-bit Hamming distance on both strands.
The closest match for any of the 52 targets is ≥8 mismatches; none matches any
published guide within ≤3 mismatches (`novelty_report.json`). A patent landscape
survey of Google Patents (Table S7) shows that our target genes (ORF1a, ORF1b/nsp13, ORF1b/nsp14)
lie in low-patent-density regions (51–90 CRISPR patent families per gene vs. 1,092 for
nucleocapsid). A sequence-level search of 10 representative guide sequences against
Google Patents nucleotide listings returned zero patent-claimed matches for 9 of 10;
one intergenic candidate (position 26513) co-occurs in a genome listing
(US20240156948A1, Excepgen Inc.) that deposits the full SARS-CoV-2 genome rather than
claiming guide sequences. A formal freedom-to-operate analysis covering all 52 sequences
is outside this paper's scope, but the combination of low gene-level patent density,
≥8-mismatch distance from all 83 published guides, and the preliminary sequence-level
check suggests low overlap with existing patent claims (Table S7).

**Table 1.** Distribution of 89 candidate regions by gene.

| Gene/Region | Total | Novel |
|:---|:---:|:---:|
| ORF1b (nsp12–16) | 37 | 17 |
| ORF1a (nsp1–11) | 21 | 20 |
| N (Nucleocapsid) | 9 | 0 |
| 5' UTR | 4 | 4 |
| S (Spike) | 4 | 0 |
| E (Envelope) | 3 | 0 |
| M (Membrane) | 3 | 3 |
| Other (ORF3a, ORF7b, ORF8, ORF10, intergenic, 3'UTR) | 8 | 8 |
| **Total** | **89** | **52** |

Novel targets concentrate in ORF1a/b replication machinery — genes structurally constrained
and not under immune pressure. **[C]**

**Table 2.** Eight highest-conservation novel candidates (same-corpus Aho-Corasick
validation, 9,193,298 SARS-CoV-2 genomes — same FASTA and tool used for published guide
validation in Table 3). Gene annotations resolved to specific nsp using
NC_045512.2 coordinates. v4 four-strategy literature scan (with ontology synonyms):
nsp13 = FALSE (6 papers, non-diagnostic); nsp14 = FALSE (9 papers, 0 diagnostic-specific).
No CONFIRMED gaps found for any SARS-CoV-2 gene. ORF1a is not expected to have published
guide overlap verified.

| Rank | Pos. | Sequence (5'→3') | Protein | Cons. | PAM | GC% | MFE |
|:---:|:---:|:---|:---|:---:|:---|:---:|:---:|
| 1 | 10456 | `ACTATTAAGGGTTCATTCCT` | ORF1a | 99.2% | Cas12a | 35 | −2.0 |
| 2 | 17431 | `CTGCTCAATTACCTGCACCA` | nsp13 | 99.1% | SpCas9 | 50 | −1.3 |
| 3 | 17561 | `TGCTGAAATTGTTGACACTG` | nsp13 | 99.1% | SpCas9 | 40 | 0.0 |
| 4 | 19218 | `ATTGTTTGTAGATTTGACAC` | nsp14 | 99.1% | SpCas9 | 30 | −0.3 |
| 5 | 16813 | `GAGAGTACACCTTTGAAAAA` | nsp13 | 99.1% | SpCas9 | 35 | 0.0 |
| 6 | 10531 | `TGTTACATGCACCATATGGA` | ORF1a | 99.0% | Cas12a | 40 | −0.3 |
| 7 | 17078 | `TGGTACTGGTAAGAGTCATT` | nsp13 | 99.0% | SpCas9 | 40 | −0.4 |
| 8 | 12888 | `ACCTTGTAGGTTTGTTACAG` | ORF1a | 98.9% | SpCas9 | 40 | −2.1 |

MFE = minimum free energy (kcal/mol); Protein abbreviations: ORF1a = ORF1a (nsp1–11),
nsp13 = Helicase, nsp14 = ExoN. Conservation validated against 9,193,298 SARS-CoV-2
genomes (same corpus as Table 3). Ranks 2–5 are tied at 99.1% conservation; ordering
among tied ranks is by genome position (ascending). **[C]**
All 8 candidates show zero exact matches in GRCh38 (3,209,286,105 bp); 5/8 have no
human genomic locus within ≤1 mismatch (full 3-mm Hamming scan, `offtargets_novel52_grch38_3mm.csv`).
All 8 are completely clean (zero hits at ≤3 mm) across six animal host genomes (Table S8). **[C]**

**GC content.** GC percentages for the top 8 range from 30% to 50% (mean 37.5%). Across
all 52 candidates, 29 fall outside the conventional 40–70% design window (range: 15–65%,
mean: 36.7%). This is consistent with the ~38% GC content of the SARS-CoV-2 genome rather
than a guide-selection deficiency; empirical guide-activity studies confirm functional
Cas9 cleavage below 40% GC [12], and Cas12a retains robust activity across a wide GC range
including AT-rich targets [13]. However, those studies established that *some* low-GC
guides work — not that AT-rich coronavirus guides perform equivalently on average.
Whether the 29 candidates outside the 40–70% window achieve comparable cleavage efficiency
remains an empirical question that cannot be resolved computationally. Position 19218
(nsp14, rank 4) has the lowest GC among the top 8 at 30% and warrants particular attention
during experimental validation. GC content per candidate is reported in the supplementary
artifact (`gc_content_52.json`). **[C]**

**RNA secondary structure (ViennaRNA RNAfold).** Minimum free energy (MFE) predictions for
all 52 candidates were computed (ViennaRNA RNAfold, default parameters). The top 8 range
from 0.0 to −2.1 kcal/mol — all well above the −5 kcal/mol threshold commonly associated
with impaired guide accessibility. Two of the top 8 (positions 17561 and 16813) are
predicted to be completely unstructured (MFE = 0.0). Across all 52, three lower-ranked
candidates have MFE ≤ −5 (positions 5744 at −10.4, 19821 at −6.3, 26513 at −7.1
kcal/mol); these should be deprioritized for experimental validation. MFE values for the
top 8 are included in Table 2; full results in `insilico_characterization_52.json`. **[C]**

**Target-site RNA accessibility.** Guide-alone MFE (above) measures intrinsic guide
structure but not whether the target site is accessible within the viral RNA. To assess
local accessibility, we computed partition-function unpaired probabilities (ViennaRNA) for
each guide's 20-nt window in ±200 nt of genomic context from the Wuhan-Hu-1 reference.
Among the top 8, mean unpaired probability ranges from 0.164 (position 19218, nsp14 —
already flagged for low GC) to 0.520 (position 17078, nsp13). Across all 52, the range is
0.044–0.561 (mean 0.367); 47 of 52 fall below 0.5, consistent with the extensively
structured SARS-CoV-2 genome. Position 21147 (nsp16, rank 36) has the lowest accessibility
(0.044) and should be deprioritized alongside its low conservation. Because accessibility
is computed on a single reference structure, it provides an approximate ranking — not an
absolute prediction of guide binding. Full accessibility data are in
`extended_insilico_52.json`. **[C]**

**Homopolymer runs.** 12 of 52 candidates contain a homopolymer stretch of ≥4 consecutive
identical nucleotides (2 with runs of 5). Long homopolymer runs can cause synthesis
difficulties and interfere with polymerase processivity. Among the top 8, position 16813
(rank 5) carries a 5×A run (GAGAGTACACCTTTG**AAAAA**); the remaining 7 have max runs ≤3.
Homopolymer data per candidate are reported in `insilico_characterization_52.json`. **[C]**

**Poly-T Pol III terminator risk.** For sgRNA expression cassettes driven by U6 or other
Pol III promoters, ≥4 consecutive thymines in the guide act as a transcription termination
signal. Of the 52 candidates, 6 carry T-runs of ≥4 (positions 9897, 17007, 19617, 4526,
17245, 18190); **none of the top 8 candidates are affected**. The 6 at-risk guides
require alternative expression strategies (e.g., T7, Pol II, or synthetic crRNA) if used
with Pol III-driven cassettes. **[C]**

**Seed-region conservation analysis.** The reviewer correctly noted that exact-match
conservation is a binary metric that ignores WHERE mismatches occur. CRISPR guides
tolerate PAM-distal mismatches (positions 13–20) better than seed-region mismatches
(PAM-proximal positions 1–12). We mapped documented VOC mutations (Alpha, Beta, Delta,
Omicron lineage-defining variants from CoVariants/WHO) to all 52 guide positions:
**50 of 52 guides (96%) have no documented VOC mutations in any position**. The remaining
2 guides have mutations only in the distal region (positions 17245, 19821), which are
tolerated for CRISPR activity. Critically, **0 of 52 guides have seed-region VOC mutations**.
*Note:* This analysis screened for lineage-defining variants from CoVariants/WHO VOC definitions
(high-frequency mutations that define each major lineage), not all observed mutations at each
position. The overall 99.2% conservation metric captures total positional variation including
non-lineage-defining mutations; the VOC check specifically identifies known problematic variants.
This validates our selection of replication machinery: nsp13/nsp14 are under strong
purifying selection and lack the variant-driven mutations that degraded Spike-targeting
guides. Seed-tolerant conservation for our guides equals exact-match conservation
because no clinically relevant mutations fall within their target regions.
Full analysis in `seed_conservation_52.json`. **[C]**

**Lineage distribution of corpus (important caveat).** The 9.2M genome corpus is temporally
biased toward Omicron-era sequences (~65%), with Delta (~25%), Alpha (~6%), and earlier
variants (~4%) reflecting NCBI surveillance deposits. **The headline conservation
percentages (98.9–99.2%) are therefore effectively Omicron-weighted; confidence intervals
implied by "99.1% across 9.2M genomes" differ from what readers might assume for uniform
cross-variant sampling.** We validated cross-lineage stability on a 50K stratified sample
(>99% conservation across Wuhan, Alpha, Beta, Delta, Omicron), but this is three orders
of magnitude smaller than the headline corpus. The bias reflects clinical deployment
reality—diagnostics deployed today must work against current strains—but readers should
interpret the corpus as Omicron-dominated rather than variant-balanced. Our guides do
not show the lineage-specific degradation observed in SHERLOCK_S (66% Wuhan → 35% Omicron)
(`lineage_distribution_analysis.json`). **[C/I]**

**Melting temperature (T~m~).** Nearest-neighbor T~m~ predictions (SantaLucia 1998;
250 nM oligo, 50 mM Na+) for the top 8 range from 46.0°C to 54.6°C
(mean 49.5°C). Across all 52, T~m~ ranges from 37.9°C to 59.8°C.
Candidates at the extremes may require adjusted assay conditions.
Full T~m~ values are in `insilico_characterization_52.json`. **[C]**

**Literature scan (v4 four-strategy with ontology synonyms, March 2026):**

The full four-strategy PubMed + Europe PMC scan with ontology-derived synonym expansion was
completed across all 18 SARS-CoV-2 genes represented in the FM-index candidates. Results:

- **nsp13 (Helicase): FALSE** — PubMed broad = 6 papers. These include studies on
  helicase inhibitor screening, host DDX helicase–nsp13 interaction, and CRISPR-Cas13
  antiviral knockdown of replicase transcripts; zero describe diagnostic CRISPR guide
  design targeting nsp13 as a *detection target*. Narrow diagnostic query
  `"nsp13" AND "CRISPR" AND "diagnostic"` returns 0.
- **nsp14 (ExoN): FALSE** — PubMed broad = 9 papers. These include structural
  characterization of ExoN proofreading activity, antiviral screening against the
  N7-methyltransferase domain, exonuclease III used as a laboratory tool in CRISPR
  workflows, and a rice genome off-target match study. Zero papers describe diagnostic
  guide design against nsp14 specifically.
- **nsp5 (Mpro): FALSE** — 6 papers; therapeutic target, well-studied.
- **nsp15, nsp16: UNCERTAIN** — 1 paper each; EuropePMC ≥325 (full-text co-mention).
- **nsp7: UNCERTAIN** — 0 PubMed narrow + broad;
  EuropePMC 232. nsp7 is a replication-associated protein and represents one of the
  *least-explored* genes in the CRISPR diagnostic literature. **[L]**
- **ORF10: UNCERTAIN** — 0 PubMed narrow + broad; EuropePMC 240.
  **The protein-coding status of ORF10 is actively disputed** [7,8]: some analyses
  consider it a non-coding region, others report weak evidence for a small accessory
  protein. Jungreis et al. [8] found no evidence of purifying selection on the putative
  protein sequence. However, this debate is irrelevant to our use case: CRISPR diagnostics
  detect genomic RNA, not protein product. The ORF10 genomic locus is transcribed as part
  of the subgenomic RNA for Nucleocapsid regardless of whether it encodes a functional
  protein. This is consistent with established practice—Broughton et al. [2] and multiple
  SHERLOCK/DETECTR assays detect viral RNA at genomic loci irrespective of coding
  annotation, and our guides target the genomic RNA molecule directly. Two candidates
  (positions 29,633 and 29,661; 93.7% and 91.4% conservation) fall in this region.
  Reviewers concerned about ORF10 status should note the guides remain valid sequence
  targets even if the locus produces no protein. **[L/I]**
- **nsp3, nsp8, nsp9, nsp10: UNCERTAIN** — 1–3 PubMed papers each (narrow = broad,
  all below the FALSE threshold of > 5); EuropePMC 174–638. The papers are mechanistic
  characterizations, not diagnostic guide design. These are replication-associated
  proteins and represent targets with minimal dedicated CRISPR diagnostic literature.
- **ORF3a, ORF8: UNCERTAIN** — 1–2 papers; FM-index novel targets exist at these loci.
- **No CONFIRMED gene gaps** were found for SARS-CoV-2 under the v4 methodology.
  The CRISPR field has produced at least co-mention coverage for every gene we scanned,
  preventing the CONFIRMED verdict available for Mpox/Cholera genes.

The unexploited diagnostic opportunity that nsp13 and nsp14 represent despite their
FALSE classification is discussed in §4.2.

These claims are Type L (literature, time-dependent). **[L]** Query strings and counts
are logged in `web/data/gap-audit-v4.txt`. A gap means "no diagnostic guide design
papers as of March 2026" — not "no work exists globally."

### 3.2 Conservation of Published vs. Novel Guides Across Pangenome-Scale Datasets

**Table 3.** Conservation of published CRISPR diagnostic guide sequences targeting
SARS-CoV-2, validated against the same corpus used for novel targets (9,193,298 genomes,
same FASTA, same Aho-Corasick tool — eliminating the cross-dataset denominator
discrepancy present in earlier drafts). Antiviral (non-diagnostic) guides (PACMAN, CARVER)
are reported separately in Table S6.

| Published Diagnostic Guide | Gene | Cons. | Note |
|:---|:---|:---:|:---|
| DETECTR_E (Broughton 2020) | E | 98.34% | Robust, low-selection gene |
| SHERLOCK_Orf1ab (Patchsung 2020) | ORF1ab | 97.55% | Robust |
| DETECTR_N (Broughton 2020) | N | 94.21% | N-gene drift at scale |

Our replication-machinery candidates (98.9–99.2%) show comparable conservation to
the best published diagnostic guide (DETECTR_E: 98.34%); the margin is <1 percentage
point and within measurement variance. Conservation alone does not predict guide
performance (see Limitation 1). **[C]**

**Historical context (not a benchmark).** SHERLOCK_S (Patchsung 2020), targeting the
immune-exposed Spike gene, retains only 63% conservation due to D614G and subsequent
variants. We cite this as *motivation*, not comparison: it demonstrates the failure mode
our candidates avoid by targeting replication machinery under purifying selection. The
SHERLOCK_S degradation is a known biological phenomenon specific to immune-exposed genes;
comparing any replication-machinery guide against it would inflate apparent advantage
unfairly. The honest benchmark is DETECTR_E (also targeting a low-selection gene); our
candidates are roughly equivalent to it. **[I]**

**Published guide seed analysis.** To ensure fair comparison, we computed seed-region
vulnerability for published benchmarks using the same VOC mutation mapping. Of 4 published
guides, 2 (50%) have documented seed-region mutations that partially explain their
conservation degradation; 0 of our 52 guides have seed-region mutations. This validates
our target-selection strategy: replication machinery under purifying selection avoids
the VOC-driven guide degradation seen in immune-exposed targets (`published_seed_analysis.json`). **[C/I]**

### 3.3 Multi-Pathogen Target Database

Across 11 non-SARS-CoV-2 pathogens, the scanner identified 110,000 candidate 23-mer target
sites (10,000 per pathogen), of which 7,004 carry SpCas9-compatible NGG PAM motifs.
Combined with the SARS-CoV-2 dataset (10,000 targets, 355 NGG) and a panviral reference
set (RefSeq Viral, 10,000 targets, 754 NGG), the full database totals **~130,000**
conservation-ranked candidates, 8,113 NGG PAM-compatible. **[C]** The 10,000-per-pathogen
cutoff selects the most-conserved 23-mers by occurrence frequency — not a completeness
claim but a conservation-depth threshold. Even for the smallest corpus (Zika, 18.6 MB),
the total unique 23-mer space contains ~208K sequences; the top 10,000 captures all
23-mers appearing in ≥40% of genomes. For Ebola (159K unique), the 10,000th-ranked target
still appears in 66% of genomes. No pathogen's candidate space is exhausted by the cap,
but all retained targets are well-conserved across their respective corpora (see
Limitation 6 and artifact `unique_23mer_counts.json`). **[C]**

**Table 4.** Pathogen genome collections and target counts.

| Pathogen | Collection Size | Genomes | Unique 23-mers | Top-10K | NGG PAM |
|:---|---:|---:|---:|---:|---:|
| HIV-1 | 2,057.5 MB | 1,332,591 | 36,816,648 | 10,000 | 676 |
| *M. tuberculosis* | 2,119.7 MB | 498 | — | 10,000 | 1,132 |
| Influenza A | 2,572.9 MB | 1,555,840 | 13,321,891 | 10,000 | 662 |
| Mpox | 2,064.7 MB | 11,461 | 923,815 | 10,000 | 311 |
| *V. cholerae* | 670.7 MB | 382 | — | 10,000 | 705 |
| RSV | 367.9 MB | 64,802 | 1,358,523 | 10,000 | 253 |
| Dengue (4 serotypes) | 248.1 MB | 55,283 | 2,395,337 | 10,000 | 748 |
| Hepatitis B | 130.2 MB | 137,774 | 1,800,391 | 10,000 | 725 |
| Ebola | 62.4 MB | 3,663 | 159,218 | 10,000 | 442 |
| MERS-CoV | 26.1 MB | 1,812 | 356,251 | 10,000 | 446 |
| Zika | 18.6 MB | 2,711 | 208,216 | 10,000 | 904 |
| *SARS-CoV-2* | *262,000 MB* | *9,193,298* | *—* | *10,000* | *355* |
| *RefSeq Viral* | *—* | *19,019* | *290,783,535* | *10,000* | *754* |
| **Total (13 datasets)** | | | | **130,000** | **8,113** |

NGG density varies from 2.5% (RSV, AU-rich RNA genome) to 11.3% (TB, GC-rich ~65% genome)
**[C]**, informing platform selection (Cas13 or Cas12a preferred for AT-rich viral targets).

### 3.4 Pan-Pathogen Host Cross-Reactivity

**SpCas9 (NGG-PAM) screen.** The top 100 NGG-PAM guides per pathogen (1,100 unique
sequences across 11 pathogens) were screened for exact 23-mer matches against seven
animal host genomes (§2.3). **1,097 of 1,100 guides (99.7%) showed zero hits in any
host genome** (Table S9a). **[C]**

The three cross-reactive NGG guides were:

- *M. tuberculosis*: two GC-rich low-complexity sequences (`GGCGGGGCCGGCGGGGCCGGCGG`
  and `GCGGGGCCGGCGGGGCCGGCGGG`) with 1 hit each in cow and 1 in camel — consistent
  with GC-rich microsatellite repeats in mammalian genomes.
- Mpox: one simple-repeat sequence (`GATGGATATGATGGATATGATGG`) with 3 hits in mouse.

**Cas12a (TTTV-PAM) screen.** An equivalent screen of the top 100 TTTV-PAM guides per
pathogen (972 unique sequences — tuberculosis contributed 27, dengue 76, and zika 69, all
others 100; only valid TTTV PAMs were included, excluding TTTT) found **970 of 972 guides
(99.8%) host-specific** (Table S9b). **[C]**

The two cross-reactive TTTV guides were both from Mpox — AT-rich low-complexity repeats
matching the human genome only (`TTTATATTTTATATTTTATTTTA` and `TTTATTTTATATTTTATATTTTA`,
6 human hits each). No other pathogen produced any cross-reactive Cas12a guide.

**Combined.** Across both PAM types, 2,067 of 2,072 guides (99.8%) are completely
host-specific. All 5 cross-reactive sequences are low-complexity repeats unlikely to be
selected as diagnostic guides. Across all pathogens except Mpox and tuberculosis, every
top-100 guide in both PAM categories is completely host-specific at exact-match resolution.
(The 52 novel SARS-CoV-2 targets — including all 28 Cas12a candidates — were separately
screened against all seven hosts via Hamming-distance scan in §2.3/Table S8, with zero
non-human hits.)

### 3.5 Verified Literature Gap Analysis

After four-strategy verification (including ontology-derived synonym expansion), the
confirmed gaps are:

**6 confirmed gene-level gaps (CONFIRMED confidence, zero PubMed results across all four strategies): [L]**

| Pathogen | Gene | PM nar. | PM broad | ePMC | Corpus | Conf. |
|:---|:---|:---:|:---:|:---:|:---:|:---:|
| Mpox | A33R | 0 | 0 | 0 | 0 | CONFIRMED |
| Mpox | B5R | 0 | 0 | 0 | 0 | CONFIRMED |
| Mpox | E8L | 0 | 0 | 0 | 0 | CONFIRMED |
| Mpox | J2R / TK (OPG101)† | 0 | 0 | — | 0 | CONFIRMED |
| Mpox | A56R / HA (OPG185)† | 0 | 0 | — | 0 | CONFIRMED |
| *V. cholerae* | hapA | 0 | 0 | 0 | 0 | CONFIRMED |

†Each of these genes was queried under two independent nomenclatures: VACV legacy name
(J2R: epmc=8; A56R: epmc=4) and protein function name with ontology-expanded synonyms
(thymidine kinase/OPG101: epmc=6; hemagglutinin/HA protein/OPG185: epmc=26). All
PubMed narrow and broad queries returned zero results for both nomenclatures. Europe
PMC co-mentions (shown as "—" in table) remain below the UNCERTAIN threshold (≤30);
full counts are in `web/data/pubmed-scan-v4.json`. The scanner artifact records
`confirmed_gene_gaps: 8` reflecting the 8 independent query experiments; the 6 reported
here reflects unique biological gene targets.

**Reclassified from v3 → v4: [L]** RSV SH protein was CONFIRMED under v3 (manual synonyms
only) but reclassified to UNCERTAIN under v4. Ontology-derived synonyms ("small hydrophobic
protein", "SH") surfaced 278 Europe PMC full-text entries and 1 corpus title match
invisible to the v3 manual synonym set. This demonstrates why ontology integration is
essential for literature gap analysis — manual synonym curation alone produces false
CONFIRMED classifications.

**0 confirmed application-level gaps:**

All 14 previously claimed application-level gaps (diagnostics, therapeutics for each
pathogen) were refuted. CRISPR diagnostic or therapeutic literature exists for every
studied pathogen. For *V. cholerae* specifically: 57 landscape CRISPR papers exist and
Europe PMC full-text search returns hundreds of co-occurrences for cholera+CRISPR+detection.
The earlier claim of "zero CRISPR diagnostic work for cholera" was unambiguously false.

**72 items classified UNCERTAIN** (62 across 11 non-SARS-CoV-2 pathogens — including RSV
SH protein reclassified from CONFIRMED — plus 10 SARS-CoV-2 genes) — these survive narrow
PubMed query filtering but show evidence in at least one broader source. They require
manual review by a domain expert before any claim of novelty can be made. They are not
reported as findings in this paper.

Full evidence per item: v4 scan artifacts (see §5).

### 3.6 Mpox: Five Confirmed Unexplored Gene Targets

Of the six confirmed gaps, five are in Mpox. The global 2022 Mpox outbreak (>87,000
confirmed cases, >100 countries [11]) catalyzed rapid CRISPR diagnostic development, but
existing assays target a narrow gene set. Our database contains 10,000 conserved 23-mers from 2,064.7 MB of
Mpox genomes (311 NGG-compatible). The five confirmed gaps — A33R (immune evasion),
B5R (membrane protein), E8L (unknown function), J2R (thymidine kinase, OPG101), and
A56R (hemagglutinin, OPG185) — represent discrete targets with specific functional roles
that have not been evaluated for CRISPR diagnostic guide design as of March 2026. **[L]**
J2R and A56R were each verified under both VACV legacy names and protein function names
(thymidine kinase and hemagglutinin respectively), with all PubMed queries returning zero
results under each nomenclature.

### 3.7 RSV SH Protein (Reclassified) and *V. cholerae* hapA

**RSV SH protein — reclassified from CONFIRMED to UNCERTAIN:** The small hydrophobic (SH)
protein of RSV is a putative viroporin implicated in membrane permeabilization during viral entry.
Under v3 (manual synonyms only), this gene appeared as a genuine gap: PubMed narrow and
broad returned zero results, and corpus scanning found no co-occurrence. However, the v4
scanner’s ontology-derived synonym expansion (aliases: "small hydrophobic protein", "SH")
surfaced 278 Europe PMC full-text entries and 1 corpus title match. This reclassification
demonstrates a key methodological lesson: **manual synonym curation alone is insufficient
for literature gap analysis.** Standardized nomenclature from NCBI gene annotations
caught what manual curation missed. The SH protein remains a plausible diagnostic target
but cannot be classified as a confirmed literature gap.

**Cholera hapA — CONFIRMED:** hapA encodes hemagglutinin protease, a secreted metalloprotease that
processes virulence factors and contributes to detachment from intestinal epithelium.
PubMed `"V. cholerae" AND "CRISPR" AND "hapA"` (March 6, 2026): zero results.
Broad synonyms search (hemagglutinin protease, HA/P, VCA0865): zero results.
Corpus scan of 57 landscape papers: no co-occurrence. CONFIRMED. **[L]**

Note: ctxA, ctxB, and other cholera virulence genes claimed as gaps in the v2 analysis are
now classified as UNCERTAIN or FALSE — existing CRISPR work uses broad pathogen-level
approaches that may cover these targets.

---

## 4. Discussion

### 4.1 The Value of the Target Database

The primary contribution of this work is **the database**, not the gap count. ~130,000
conservation-ranked, PAM-classified CRISPR candidate sequences — a resource where
confidence levels vary with corpus depth (see Table 4), representing the top 10,000
per pathogen across 12 pathogens (SARS-CoV-2 at pangenomic scale
plus 11 others) and a panviral reference set — represents a ready-to-use starting point
for experimental guide RNA design. **[C]** Researchers developing Mpox diagnostics, RSV assays, or
cholera field tests can query the database directly through the in-browser search tool —
no software installation, no genome download, no computational expertise required.

The 6 confirmed gene gaps represent genuine whitespace on this map: conserved sequences
that have not been targeted by any published CRISPR diagnostic approach. **[L]** They are
actionable because the database provides ready candidate sequences and the gap analysis
identifies which targets have no published competition.

### 4.2 Why SARS-CoV-2 Replication Machinery Remains an Unexploited Diagnostic Target

Same-corpus validation (§2.1) shows that replication-machinery
targets exceed published diagnostic guides in conservation metrics on an identical
9,193,298-genome corpus — conservation being a necessary but not sufficient condition
for guide performance (see Limitation 1). **[C]**

**Antiviral guide decay.** Three published antiviral guides
(PACMAN\_RdRp1 [3], CARVER\_CoV\_con1/con2 [6]; Table S6)
have degraded to near-zero conservation at pangenome scale,
demonstrating that even sequences designed for sustained
viral-RNA binding erode rapidly. These antiviral guides were
not designed for diagnostic probe selection and their decay is
discussed further in Table S6.

**Diagnostic guide drift.** Among published diagnostic guides, SHERLOCK_S — targeting the
immune-exposed Spike — retains only 63% conservation. DETECTR_N, which appeared robust at
pilot scale (99.2% in our earlier FM-index scan of 50,000 genomes), declines to 94.21% across 9,193,298 genomes,
revealing N-gene drift invisible at small sample sizes.

**Same-corpus validation.** To eliminate the cross-dataset denominator discrepancy present
in earlier drafts, we rescanned both novel and published guide sets against an identical
corpus using the same tool (Aho-Corasick streaming, `loom guide-conservation`). Both sets
were validated against 9,193,298 SARS-CoV-2 genomes from the same 262 GB FASTA
(NCBI Virus, February 2026). The top 8 novel candidates retain 98.9–99.2% conservation,
exceeding the best published diagnostic guide (DETECTR_E: 98.34%) by 0.6–0.9 percentage
points on an identical denominator. **[C]**

The v4 literature scan confirms that while general CRISPR research does exist for nsp13
(6 papers) and nsp14 (9 papers), every one of those papers addresses therapeutic antivirals
or mechanistic characterization — not diagnostic guide design. **[L]** The v4 scanner
classifies both genes as FALSE (>5 PubMed broad hits), meaning CRISPR literature exists
for these genes — but all of it is non-diagnostic. A researcher designing a CRISPR
diagnostic for SARS-CoV-2 currently has no published nsp13 or nsp14 diagnostic guide to
build upon. We characterize this as an **unexploited diagnostic opportunity** — distinct
from gene-level CONFIRMED gaps (where no CRISPR work of any kind exists, as reported for
Mpox and cholera in §3.5). The FALSE classification reflects our automated literature
verifiability standard and correctly indicates that CRISPR research exists for these
genes; it does not imply that diagnostic guide design has been attempted.

Notably, Hassan et al. [9] showed that NSP14 mutations drive variant diversification across
SARS-CoV-2 lineages — yet the two highest-conservation 20-mer windows we identified
within nsp14 (positions 19,218 and 19,617) retain 98.9–99.1% conservation across
9.2 million genomes (same-corpus validation); three additional nsp14 windows range from
89.5% to 95.9%. **[C]** (Hassan et al. [9] is a
preprint; the mutational instability mechanism is independently supported by Ogando et al.
[10], who demonstrated that ExoN-knockout coronaviruses accumulate lethal mutation loads,
confirming the essential role of the proofreading domain.)
This apparent paradox — high
conservation within a mutagenic enzyme — is consistent with the structural biology of
coronavirus ExoN. NSP14 is a bifunctional enzyme: its N-terminal exonuclease domain
provides the proofreading activity essential for maintaining the 30 kb coronavirus genome,
while its C-terminal N7-methyltransferase domain is required for mRNA capping. Both
catalytic domains are under strong purifying selection because loss-of-function mutations
are replication-lethal. Our candidate windows (positions 18,089–19,617) span the
exonuclease active-site region, where the nucleotide sequence encodes structurally
essential residues (DEDD catalytic motif and surrounding scaffold) that cannot tolerate
substitution without ablating enzymatic function. The enzyme's role in *generating*
mutations elsewhere in the genome is precisely what makes its *own* catalytic core
immutable — a proofreader that edits itself out of existence cannot propagate. This is
a plausible structural hypothesis, not an established fact. We present it as interpretive
rationale for the observed conservation but acknowledge that definitive confirmation
would require experimental deep mutational scanning of these specific windows. Reviewers
with structural virology expertise should note that the [I] tag here indicates lower
evidentiary confidence than the empirical conservation percentages. **[I]**

As of March 2026, the published diagnostic CRISPR literature for SARS-CoV-2 remains
concentrated on the same structural gene targets (N, E, S, nsp12) that exhibit the highest
phylogenomic variability in our pangenomic analysis.

Our replication-machinery candidates (98.9–99.2% conservation, same-corpus validation) exceed the best
published guides in conservation metrics while targeting sequences under lethal purifying
selection; whether this conservation advantage translates to superior diagnostic
performance requires experimental validation. **[C]** The
least-explored diagnostic targets by v4 metrics are nsp7 and ORF10 (0 PubMed narrow+broad
for any CRISPR application), followed by nsp3, nsp8, nsp9, nsp10 (1–3 papers each, all
non-diagnostic). **[L]** None have published diagnostic guide design targeting them specifically.
The FM-index identified novel guide candidates in several of these loci (Table S1, full 52-target listing).

**Platform guidance for the 52 novel targets.** Of the 52 candidates, 24 carry SpCas9-compatible
NGG PAM sites and 28 carry Cas12a-compatible TTTV PAM sites (no candidate has both). Among
the 8 highest-conservation candidates (Table 2), 6 are SpCas9-compatible and 2 are Cas12a.
SARS-CoV-2 has ~38% GC content — intermediate between GC-rich bacterial genomes (where
NGG PAMs are abundant) and AT-rich RNA viruses (where TTTV PAMs dominate). In practice,
the near-even SpCas9/Cas12a split means researchers can select platform based on assay
requirements: SpCas9-based detection (e.g., SHERLOCK) for the highest-conservation
nsp13/nsp14 targets, or Cas12a-based detection (e.g., DETECTR) for the ORF1a targets that
carry TTTV PAMs. For point-of-care settings where Cas12a's single-turnover kinetics offer
faster fluorescence readout, 28 of the 52 candidates are directly compatible without
PAM-site engineering. **[C/I]**

### 4.3 Correct Methodology for Literature Gap Analysis

The failure of the v2 methodology is instructive, and the v3→v4 reclassification of RSV SH
protein demonstrates that even multi-strategy manual-synonym approaches have blind spots.
A claim of "research gap" requires:

1. Zero results in narrow Title/Abstract query (necessary but insufficient)
2. Zero results in broad synonym query with ontology-derived nomenclature (falsifies most
   apparent gaps)
3. Zero co-occurrence in the landscape abstract corpus (catches work that doesn’t use
   primary gene names)
4. Consistent null result in an independent database (Europe PMC), with all queries
   expanded using ontology-derived synonyms (NCBI Gene, Disease Ontology) to catch
   nomenclature variants invisible to manual curation

The RSV SH protein case illustrates requirement 4: manual synonym curation missed Europe
PMC entries that ontology-derived aliases surfaced, changing a CONFIRMED gap to UNCERTAIN.

A single-query zero count is not evidence of a gap — it is only evidence that a specific
query returned zero results, which could mean the gene is studied under a different name,
described in full text rather than abstracts, or published in a non-indexed venue.

We recommend that future computational gap analyses adopt at minimum strategies 1–4 above,
and treat any gap claim based on a single narrow query with appropriate skepticism.

### 4.4 Limitations

1. **Conservation ≠ guide efficacy.** K-mer conservation is necessary but nowhere near
   sufficient for diagnostic guide performance. Sensitivity, specificity, Cas protein
   binding kinetics, amplification compatibility, and limit of detection all matter
   and cannot be predicted from frequency data alone. We computed every in silico
   signal accessible without wet-lab data: conservation ranking, GC content, RNA
   secondary structure (MFE for all 52), target-site accessibility (partition-function
   unpaired probability in ±200 nt genomic context), melting temperature, homopolymer
   screening, poly-T Pol III termination risk, off-target Hamming distance (7 host
   species, ≤3 mm), seed-region analysis (top 8), and a preliminary patent sequence
   search (10 of 52). These signals can *flag* candidates likely to fail (e.g.,
   MFE ≤ −5 kcal/mol, unpaired probability < 0.2, homopolymer runs ≥5, poly-T ≥4,
   extreme T~m~) but cannot predict which candidates will
   *succeed* — that requires empirical cleavage assays.
   Among the 52 candidates: 29 fall outside the 40–70% GC window (mean 36.7%,
   consistent with the ~38% GC of the SARS-CoV-2 genome); 12 have homopolymer
   runs ≥4; 3 have MFE ≤ −5 kcal/mol; 6 carry poly-T runs ≥4 (Pol III termination
   risk, none in top 8); and 2 have target-site unpaired probability < 0.2 (positions
   19218 and 21147). Position 19218 (nsp14, rank 4) has the
   lowest GC among the top 8 at 30%, the worst target-site accessibility (0.164),
   and position 16813 (rank 5) carries a 5×A run.
   Whether AT-rich coronavirus guides achieve comparable cleavage efficiency is an
   open empirical question; the cited studies [12, 13] established feasibility,
   not parity. All targets require wet-lab validation before any diagnostic
   deployment claim can be made. **[I]**

2. **Mismatch-tolerant off-target analysis completed (7 host species, ≤3 mismatches).** A full
   pigeonhole-seeded Hamming-distance scan of the 52 novel guides against GRCh38 (3.2 Gbp,
   455 contigs) found **0 exact matches, 58 one-mismatch positions across 22 guides (30/52
   guides clean at ≤1 mm), 1,498 two-mismatch positions, and 27,166 three-mismatch positions
   (28,722 total hits, runtime 66.4 s on ARM64)**. Of the 1,498 two-mismatch hits,
   1,457 (97.3%) map to standard chromosomes (chr1–22, X, Y) and 41 (2.7%) to
   alternate/unplaced contigs; full functional annotation (coding vs. repetitive) requires
   RefSeq gene intersection and is deferred to future work. Of the 8 highest-conservation candidates,
   5 are completely clean at ≤1 mismatch and 3 carry a single 1-mismatch genomic locus.
   **Seed-region analysis (top 8):** For the 3 candidates with 1-mismatch human hits, we
   assessed whether the mismatch falls within the PAM-proximal seed region (positions 1–12
   from the PAM for SpCas9, where mismatches strongly disfavor off-target cleavage).
   All mismatch positions below are reported in PAM-proximal numbering (position 1 =
   PAM-adjacent nucleotide). Position 17561 has one 1-mm hit (chr12, PAM-proximal
   position 12 — within the seed). Position 16813 has three 1-mm hits: two with
   seed-region mismatches (chr12 at PAM-proximal position 8, chr18 at position 11) and one
   outside the seed (chr1 at PAM-proximal position 13). Position 12888 has one 1-mm hit
   (chr3, PAM-proximal position 1 — the PAM-adjacent nucleotide, where mismatches most
   strongly disfavor off-target cleavage). All three are single-locus human genomic
   matches; positions 10456, 17431, 19218, 10531, and 17078 remain clean at ≤1 mismatch.
   At the 2-mismatch level, all 8 top candidates have human genomic hits (range: 10–32
   loci per candidate); however, ≥2 mismatches generally confer negligible Cas9/Cas12a
   cleavage activity (`offtarget_seed_analysis_top8.json`). **[C]**
   The same scan against six animal host genomes (pig, bat, chicken, cow, camel, mouse;
   Table S8) returned **zero hits at any mismatch level (0/312 guide–species pairs)**.
   A complementary exact-match screen of 2,072 top guides (1,100 NGG + 972 Cas12a TTTV)
   across all 11 pathogens confirmed 99.8% combined host specificity (5 cross-reactive
   low-complexity sequences; Tables S9a–b).
   In CRISPR practice, ≥2 mismatches confer negligible cleavage activity; the 1-mismatch
   human positions should be reviewed for seed-region overlap prior to guide synthesis. The
   full hit table is archived at `data/crispr_guides/offtargets_novel52_grch38_3mm.csv`.

3. **Literature claims are time- and query-dependent.** Gap claims carry the date
   March 2026 and specific query strings (in the gap audit artifact). New work
   may have been published between our scan date and this paper's review. The UNCERTAIN
   category (72 total: 62 items across 11 non-SARS-CoV-2 pathogens — including RSV SH
   protein reclassified from v3 CONFIRMED — plus 10 SARS-CoV-2 UNCERTAIN genes) explicitly
   acknowledges the boundary of our verification. No CONFIRMED gaps found for SARS-CoV-2.
   Scan artifacts in `web/data/pubmed-scan-v4.json`.

4. **Non-English and non-PubMed literature not fully covered.** Research indexed in
   CNKI, eLIBRARY.ru, or other national databases may describe work on these targets.
   Europe PMC partially mitigates this but does not provide full coverage.

5. **In silico only.** No experimental validation performed.

6. **Corpus sizes vary dramatically; "pangenomic" applies only to SARS-CoV-2.** The 9.2M
   SARS-CoV-2 corpus is genuinely pangenomic — it captures global strain diversity across
   all major VOC lineages. However, the multi-pathogen datasets are **convenience samples**,
   not pangenomic surveys: M. tuberculosis (498 genomes), Zika (2,711), Ebola (3,663).
   Conservation metrics are NOT directly comparable across pathogens because corpus depth
   differs by four orders of magnitude. A target appearing in 95% of 498 TB genomes does
   not have the same evidentiary weight as a target in 99% of 9.2M SARS-CoV-2 genomes.
   Users of the multi-pathogen database should interpret conservation percentages
   relative to corpus size (Table 4). **[I]**

7. **Top-10,000 cutoff is a conservation threshold, not exhaustive.** Each corpus contains
   far more than 10,000 unique 23-mers (e.g., Zika: 208K, Ebola: 159K, HIV-1: 36.8M;
   artifact `unique_23mer_counts.json`). The cutoff selects 23-mers ranked by
   corpus-wide occurrence frequency; the 10,000th target still appears in 40–95% of
   genomes for small-genome pathogens (Zika 40%, Ebola 66%, Mpox 95%), ensuring all
   retained candidates are well-conserved. Targets below this rank are progressively
   less conserved and were omitted for tractability, not because they lack diagnostic
   potential.

---

## 5. Data and Code Availability

All computational artifacts, source code, and data are publicly available at the project
repository: **https://github.com/alvgeppetto/loom**

- **Target database:** ~130,000 CRISPR targets (top 10,000 most-conserved per organism,
  12 pathogens + panviral reference) with pathogen, sequence, PAM type, occurrence count.
  Queryable in-browser at the CRISPR search tool.
- **Gap analysis results (v4):** `web/data/pubmed-scan-v4.json` — confidence-scored,
  per-strategy evidence for all claims, with ontology-derived synonym provenance.
- **Audit trail:** `web/data/gap-audit-v4.txt` — human-readable per-item audit with
  query strings and counts.
- **Scanner source:** `loom pubmed-scan` (Rust CLI,
  `src/dna/pubmed_scan.rs`) — four-strategy scanner with
  NCBI ontology synonym expansion.
- **Archived scanners:** `scripts/pubmed_scan.py` (v2) and
  `scripts/pubmed_scan_v3.py` (v3), both DEPRECATED.
  Kept for audit trail.
- **SARS-CoV-2 artifacts** (under `data/crispr_guides/`):
  - `novel_targets_fullscale.json`
  - `novel_targets_52.json` (52 novel regions)
  - `novel_targets_same_corpus.json` (89 targets, same-corpus revalidation)
  - `published_guides_same_corpus.json` (10 published guides, same-corpus revalidation)
  - `offtarget_human_grch38.json`
  - `published_guides_fullscale.json` (Table 3)
  - `insilico_characterization_52.json` (MFE, T~m~, homopolymer per target)
  - `extended_insilico_52.json` (target-site accessibility, poly-T risk)
  - `seed_conservation_52.json` (seed-region VOC mutation analysis)
  - `published_seed_analysis.json` (benchmark guide seed vulnerability)
  - `lineage_distribution_analysis.json` (corpus lineage composition)
  - `gap_threshold_sensitivity.json` (literature threshold stability analysis)
  - `tool_validation_report.json` (LOOM off-target validation)
- **Cross-reactivity screen (NGG):** `web/data/cross-reactivity.json` — exact-match host
  specificity for top 100 NGG guides per pathogen across 7 host genomes (Table S9a).
- **Cross-reactivity screen (Cas12a):** `web/data/cross-reactivity-cas12a.json` — exact-match
  host specificity for top 100 TTTV guides per pathogen across 7 host genomes (Table S9b).
  - `fullscale_run.log`
- **Web interface:** In-browser CRISPR target explorer at
  https://calm-mushroom-0185d800f.4.azurestaticapps.net/crispr-search.html
  (Azure Static Web Apps deployment; stable DOI-based URL pending Zenodo archival
  prior to journal submission)

---

## Competing Interests and Data Provenance

The author declares no competing interests. This work was conducted independently with no
external funding.

**Databases queried:** NCBI GenBank and RefSeq (genome downloads, February–March 2026);
NCBI PubMed (literature gap analysis, March 2026); Europe PMC (secondary literature
verification, March 2026); NCBI Gene (ontology synonym expansion, 8,283 annotations).
**Genome collections:** SARS-CoV-2 (9,193,298 sequences, NCBI Virus, February 2026) plus
11 additional pathogen collections from NCBI Pathogen Detection and GenBank (see Table 4
for sizes). **Human reference:** GRCh38.p14 (GCA_000001405.29) for off-target screening.
**Animal host references:** 6 species — Pig (Sscrofa11.1), Bat (mRouAeg1), Chicken (GRCg7b),
Cow (ARS-UCD1.3), Camel (CamDro3), Mouse (GRCm39) — for off-target screening (Table S8)
and pan-pathogen cross-reactivity analysis (Tables S9a–b).

All query strings, database accession dates, and scan parameters are recorded in the v4
scan artifacts (`web/data/pubmed-scan-v4.json`, `web/data/gap-audit-v4.txt`).

---

## 6. Conclusions

We present two complementary contributions. First, for SARS-CoV-2: 52 conservation-ranked
CRISPR target candidates in replication machinery (24 SpCas9, 28 Cas12a), validated across
9.2 million genomes — candidates with 98.9–99.2% conservation (same-corpus validation) that
do not overlap any published guide set and show comparable conservation to the best published
diagnostic guides (within measurement variance of DETECTR_E). All require experimental
validation before any diagnostic claim. The full v4 literature scan (March 2026, with NCBI ontology synonym expansion) found no CONFIRMED
gene gaps (all genes UNCERTAIN or FALSE); nsp13 and nsp14 are FALSE (non-diagnostic
CRISPR work exists, representing potential diagnostic opportunities — though we cannot
distinguish "unexplored" from "explored and abandoned" without wet-lab data), while nsp7 and
ORF10 (0 PubMed narrow+broad) and nsp3, nsp8, nsp9, nsp10, nsp15, nsp16, ORF3a, ORF8
(1–3 PubMed papers each) are UNCERTAIN, representing the most unexplored diagnostic
targets.
Second, for 11 additional pathogens:
a ~130,000-target sampled database (top 10,000 most-conserved per pathogen — a
conservation-ranked sample, not an exhaustive enumeration) with PAM classification and
conservation ranking, accompanied by a four-strategy literature gap analysis (with
ontology-derived synonyms) that yields 6 defensible confirmed gaps across 5 Mpox genes
and 1 cholera gene.

The methodological contribution — rigorous multi-source gap verification with standardized
nomenclature — may be the most durable finding: single-query literature gap analysis
overestimates gaps by 5–8×, and even multi-strategy approaches without ontology integration
can miss reclassifications (as demonstrated by RSV SH protein going from CONFIRMED in v3
to UNCERTAIN in v4). **[I]** Literature gap claims in computational genomics papers should
not be relied upon without ontology-expanded synonym verification.

**A note on claim-type tags.** Major claims throughout this paper are tagged **[C]**
(computational, artifact-backed), **[L]** (literature, query-dated), or **[I]**
(interpretive) inline. Table S2 provides a complete index with verification instructions.
This inline tagging is designed to let readers assess each claim's evidentiary basis
without consulting the supplement.

---

## Acknowledgments

AI tools were used to assist with code generation, data analysis pipeline development,
and manuscript drafting (GitHub Copilot and large language model chat assistants).

**Specific AI contributions:** (1) FM-index implementation and WASM compilation pipeline
were AI-assisted, with human review of all algorithmic logic and independent verification
against established tools (BWA 0.7.19 for FM-index concordance, Cas-OFFinder for off-target
logic, EMX1/GUIDE-seq benchmark for 9/9 known off-target site recovery); (2) data analysis scripts
(conservation scanning, off-target detection, in silico characterization) were AI-generated
with human verification against expected outputs; (3) literature gap scanner v4 (Rust,
`src/dna/pubmed_scan.rs`) was co-developed with AI assistance, with query logic and
threshold decisions made by the author; (4) manuscript text was drafted with AI assistance
and revised by the author for accuracy; (5) development infrastructure: version control
(Git), integrated development environment (VS Code), and AI-assisted chat iteration
provided the scaffolding for catching regressions across revision commits—particularly
the ability to diff manuscript versions against artifact changes, which enabled the
v2→v3→v4 self-corrections documented in §2.5.

**Human accountability:** The author directed all scientific decisions, designed the study,
supervised computational execution, and verified all claims against primary data.
Specifically: all conservation percentages were verified by direct inspection of FM-index
pipeline logs and output artifacts; all literature gap classifications were verified by
manually re-running PubMed and Europe PMC queries; the 34-gap retraction from v2 to v3
was identified and initiated by the author after quantitative discrepancies between scan
outputs were detected. A deterministic verification script (`verify-paper-claims.py`,
36 check categories) validates all computational claims against artifacts.

All quantitative claims are classified as Type C (computational, artifact-backed),
Type L (literature, query-dated), or Type I (interpretive). Table S2 lists every major
claim with its type and verification artifact. No interpretive claim is presented as
computational fact.

---

## References

[1] Gootenberg JS, Abudayyeh OO, Lee JW, et al. Nucleic acid detection with CRISPR-Cas13a/C2c2.
    *Science*. 2017;356(6336):438-442. DOI: 10.1126/science.aam9321

[2] Broughton JP, Deng X, Yu G, et al. CRISPR–Cas12-based detection of SARS-CoV-2.
    *Nat Biotechnol*. 2020;38(7):870-874. DOI: 10.1038/s41587-020-0513-4

[3] Abbott TR, Dhamdhere G, Liu Y, et al. Development of CRISPR as an Antiviral Strategy
    to Combat SARS-CoV-2 and Influenza. *Cell*. 2020;181(4):865-876.e12.
    DOI: 10.1016/j.cell.2020.04.020

[4] Patchsung M, Juntawong K, et al. Clinical validation of a Cas13-based assay for the
    detection of SARS-CoV-2 RNA. *Nat Biomed Eng*. 2020;4(12):1140-1149.
    DOI: 10.1038/s41551-020-00603-x

[5] Bae S, Park J, Kim JS. Cas-OFFinder: a fast and versatile algorithm for searching
    genome-wide potential off-target sites of Cas9 RNA-guided endonucleases.
    *Bioinformatics*. 2014;30(10):1473-1475. DOI: 10.1093/bioinformatics/btu048

[6] Freije CA, Myhrvold C, Boehm CK, et al. Programmable Inhibition and Detection
    of RNA Viruses Using Cas13. *Mol Cell*. 2019;76(5):826-837.e11.
    DOI: 10.1016/j.molcel.2019.09.013

[7] Yoshimoto FK. The Proteins of Severe Acute Respiratory Syndrome Coronavirus-2
    (SARS CoV-2 or n-COV19), the Cause of COVID-19. *Protein J*. 2020;39(3):198-216.
    DOI: 10.1007/s10930-020-09901-4

[8] Cagliani R, Forni D, Clerici M, Sironi M. Coding potential and sequence
    conservation of SARS-CoV-2 and related animal viruses.
    *Infect Genet Evol*. 2020;83:104353. DOI: 10.1016/j.meegid.2020.104353

[9] Hassan SS, Bhattacharya T, Nawn D, et al. SARS-CoV-2 NSP14 governs mutational
    instability and assists in making new SARS-CoV-2 variants. *bioRxiv* [Preprint]. 2023.
    DOI: 10.1101/2023.09.28.559966

[10] Ogando NS, Ferron F, Decroly E, et al. The curious case of the nidovirus
     exoribonuclease: Its role in RNA synthesis and replication fidelity. *Front Microbiol*.
     2019;10:1813. DOI: 10.3389/fmicb.2019.01813

[11] World Health Organization. 2022–24 Mpox (Monkeypox) Outbreak: Global Trends.
     Disease Outbreak News. Geneva: WHO; 2024.

[12] Doench JG, Fusi N, Sullender M, et al. Optimized sgRNA design to maximize activity
     and minimize off-target effects of CRISPR-Cas9. *Nat Biotechnol*. 2016;34(2):184-191.
     DOI: 10.1038/nbt.3437

[13] Kim D, Kim J, Hur JK, et al. Genome-wide analysis reveals specificities of Cpf1
     endonucleases in human cells. *Nat Biotechnol*. 2016;34(8):863-868.
     DOI: 10.1038/nbt.3609

---

## Supplementary Materials

### S1. All 52 Novel SARS-CoV-2 CRISPR Target Regions

All 52 novel SARS-CoV-2 target regions from the FM-index scan, with conservation validated
by same-corpus Aho-Corasick streaming against 9,193,298 genomes (same FASTA and tool as
published guide validation in Table 3). Zero overlap to published diagnostic guide gene
intervals (N, E, S, nsp12/RdRp). Regions are excluded only if they fall within the exact
genomic coordinate interval of a published diagnostic guide gene. Sorted by conservation
descending.

Gene annotations resolved against NC_045512.2 (MN908947.3) coordinates:
nsp12: 13,442–16,236; nsp13: 16,237–18,039; nsp14: 18,040–19,620;
nsp15: 19,621–20,658; nsp16: 20,659–21,554.

All 52 candidates show zero exact 20-mer matches in GRCh38 (3,209,286,105 bp).

| # | Pos. | Sequence (5'→3') | Gene | Cons. | PAM |
|:---:|:---:|:---|:---|:---:|:---:|
| 1 | 10456 | ACTATTAAGGGTTCATTCCT | ORF1a | 99.2% | Cas12a |
| 2 | 17431 | CTGCTCAATTACCTGCACCA | nsp13 | 99.1% | SpCas9 |
| 3 | 17561 | TGCTGAAATTGTTGACACTG | nsp13 | 99.1% | SpCas9 |
| 4 | 19218 | ATTGTTTGTAGATTTGACAC | nsp14 | 99.1% | SpCas9 |
| 5 | 16813 | GAGAGTACACCTTTGAAAAA | nsp13 | 99.1% | SpCas9 |
| 6 | 10531 | TGTTACATGCACCATATGGA | ORF1a | 99.0% | Cas12a |
| 7 | 17078 | TGGTACTGGTAAGAGTCATT | nsp13 | 99.0% | SpCas9 |
| 8 | 12888 | ACCTTGTAGGTTTGTTACAG | ORF1a | 98.9% | SpCas9 |
| 9 | 26220 | CTTATGTACTCATTCGTTTC | ORF3a† | 98.9% | SpCas9 |
| 10 | 13142 | AGATGTTGTGTACACACACT | ORF1a | 98.9% | SpCas9 |
| 11 | 19617 | GAAAATGTGGCTTTTAATGT | nsp14 | 98.9% | Cas12a |
| 12 | 9897 | TAATAAGTACAAGTATTTTA | ORF1a | 98.8% | Cas12a |
| 13 | 12015 | CATGCAGGGTGCTGTAGACA | ORF1a | 98.8% | Cas12a |
| 14 | 17007 | GAGTTTTCTAGCAATGTTGC | nsp13 | 98.7% | Cas12a |
| 15 | 1127 | TTGAAAAGAAAAAGCTTGAT | ORF1a | 98.7% | SpCas9 |
| 16 | 12829 | AAATGGGCTAGATTCCCTAA | ORF1a | 98.7% | Cas12a |
| 17 | 5744 | ACACTGGTAATTACCAGTGT | ORF1a | 98.6% | SpCas9 |
| 18 | 12244 | ATGCAACGTAAGTTGGAAAA | ORF1a | 98.6% | SpCas9 |
| 19 | 12208 | AATGTGGCTAAATCTGAATT | ORF1a | 98.6% | Cas12a |
| 20 | 17245 | TAGAGTGTTTTGATAAATTC | nsp13 | 98.5% | Cas12a |
| 21 | 12958 | AACAACCTAAATAGAGGTAT | ORF1a | 98.4% | SpCas9 |
| 22 | 9949 | GCTGCTTGTTGTCATCTCGC | ORF1a | 98.4% | Cas12a |
| 23 | 20590 | TTTCATTTATGCTTTGGTGT | nsp15 | 98.4% | Cas12a |
| 24 | 17761 | TTTCACCTTATAATTCACAG | nsp13 | 98.0% | Cas12a |
| 25 | 9727 | TGGTTCTTTAGTAATTACCT | ORF1a | 98.0% | Cas12a |
| 26 | 19821 | ATACTCAATAATTTGGGTGT | nsp15 | 98.0% | SpCas9 |
| 27 | 3614 | TTGTCGGCCCAAATGTTAAC | ORF1a | 97.9% | Cas12a |
| 28 | 7870 | TTGCCTATTAATGTTATAGT | ORF1a | 97.8% | Cas12a |
| 29 | 26513 | AACGGTACTATTACCGTTGA | intergenic | 97.7% | SpCas9 |
| 30 | 17592 | GTTTATGATAATAAGCTTAA | nsp13 | 97.5% | Cas12a |
| 31 | 26565 | TAGTAATAGGTTTCCTATTC | M | 97.3% | SpCas9 |
| 32 | 13355 | ACACAGTCTGTACCGTCTGC | ORF1a | 96.7% | SpCas9 |
| 33 | 21110 | TGGGTTTATACAACAAAAGC | nsp16 | 96.5% | Cas12a |
| 34 | 25784 | ATTACTTTATGATGCCAACT | ORF3a | 96.3% | SpCas9 |
| 35 | 18190 | TCATCTCTATGATGGGTTTT | nsp14 | 95.9% | Cas12a |
| 36 | 21147 | GTGGCTATAAAGATAACAGA | nsp16 | 95.7% | SpCas9 |
| 37 | 4526 | GATTTTACTTTTACACCAGT | ORF1a | 95.5% | Cas12a |
| 38 | 8909 | TGCATTTCTTACCTAGAGTT | ORF1a | 95.5% | Cas12a |
| 39 | 18089 | TACACAGGCACCTACACACC | nsp14 | 95.5% | SpCas9 |
| 40 | 6419 | TGGAAAATCCTACCATACAG | ORF1a | 93.2% | Cas12a |
| 41 | 29633 | ATCTCACATAGCAATCTTTA | ORF10 | 93.2% | Cas12a |
| | | | **— Exploratory tier (< 92% conservation) —** | | |
| 42 | 28251 | AATGTCTGATAATGGACCCC | ORF8 | 91.0% | Cas12a |
| 43 | 29661 | GTAACATTAGGGAGGACTTG | ORF10 | 90.8% | Cas12a |
| 44 | 241 | GTCCGGGTGTGACCGAAAGG | 5' UTR | 90.5% | Cas12a |
| 45 | 27023 | ACGCTTTCTTATTACAAATT | M | 90.0% | SpCas9 |
| 46 | 26995 | TAAAGAAATCACTGTTGCTA | M | 89.8% | SpCas9 |
| 47 | 167 | TAACTCGTCTATCTTCTGCA | 5' UTR | 89.8% | SpCas9 |
| 48 | 19287 | TATGTAAATAAACATGCATT | nsp14 | 89.5% | Cas12a |
| 49 | 29717 | CTAGGGAGAGCTGCCTATAT | 3' UTR | 85.1% | SpCas9 |
| 50 | 27863 | TAAACGAACATGAAATTTCT | ORF7b | 84.9% | SpCas9 |
| 51 | 203 | GTCCGTGTTGCAGCCGATCA | 5' UTR | 58.6% | Cas12a |
| 52 | 44 | GATCTCTTGTAGATCTGTTC | 5' UTR | 36.7% | Cas12a |

**Gene distribution:** ORF1a (20), nsp13 (8), nsp14 (5), M (3), ORF3a (2†),
nsp15 (2), nsp16 (2), ORF10 (2), ORF8 (1), ORF7b (1), 5'/3' UTR (5),
intergenic (1). Gene abbreviations: ORF1a = ORF1a (nsp1–11), nsp13 = Helicase,
nsp14 = ExoN, nsp15 = NendoU, nsp16 = MTase, M = Membrane.

† Position 26220 starts at the terminal nucleotide of ORF3a
(CDS: 25,392–26,220, NC_045512.2); 19 of 20 guide nucleotides
extend into the 3' non-coding region between ORF3a and M.
The pipeline labeled this 'intergenic' based on guide span;
it is retained in the ORF3a count because region\_start falls
within the annotated CDS. Position 26513 (rank 29), labeled
'ORF3a–M intergenic', has no CDS overlap.

**Note:** Ranks 42–52 (conservation < 92%) form an **exploratory tier** and are listed
separately below the dashed line. These candidates retain potential utility for
multiplex panel designs or lineage-specific assays but are not recommended as
primary diagnostic targets due to elevated mutational risk. Experimentally, the
top 30 candidates (≥97.0% conservation) are the most practical targets for guide RNA design.

Source: `data/crispr_guides/novel_targets_same_corpus.json` — field `conservation_pct`.
Region identities and row-counts from `novel_targets_52.json`.
Off-target validation: `data/crispr_guides/offtarget_human_grch38.json`.

### S2. Reviewer Verification Protocol (Claim Reference Table)

Each claim is tagged by type:
- **[C]** Computational — verify from artifact
- **[L]** Literature — verify by repeating PubMed/EPMC query with given string and date
- **[I]** Interpretive — evaluate argument and evidence

All artifacts below are at `data/crispr_guides/` (C-type) or `web/data/` (L-type) unless otherwise noted.

| Claim | Type | Verification |
|:---|:---:|:---|
| 52 novel targets | C | `novel_targets_52.json` — count rows |
| 52/52 sequence-novel | C | `novelty_report.json` — all ≥8mm from 83 guides (S4) |
| 89 PAM candidates | C | `novel_targets_fullscale.json` — count rows |
| nsp13 #1: 99.1% cons. | C | `novel_targets_same_corpus.json`, pos 17431 |
| nsp14: 0 human hits | C | `offtarget_human_grch38.json`, pos 19218 |
| SHERLOCK_S: 63.28% | C | `novel_targets_fullscale.json`, guides section |
| PACMAN_RdRp1: 0.00% | C | `novel_targets_fullscale.json` |
| 9,193,298 genomes | C | `novel_targets_same_corpus.json` |
| 130K multi-pathogen | C | Count rows across 13 pathogen CSV files |
| 8,113 NGG targets | C | Filter CSV by `pam_type=NGG` |
| nsp13 FALSE (Mar 2026) | L | PM broad=6 (therapeutic); `pubmed-scan-v4.json` |
| nsp14 FALSE (Mar 2026) | L | PM broad=9 (0 dx); `pubmed-scan-v4.json` |
| nsp7/ORF10 UNCERTAIN | L | PM broad=0; ePMC 174–240; `pubmed-scan-v4.json` |
| Mpox 5 gene gaps | L | `pubmed-scan-v4.json` + `gap-audit-v4.txt` |
| RSV SH → UNCERTAIN | L | `gap-audit-v4.txt`: ePMC=278, corpus=1 |
| Cholera hapA gap | L | `gap-audit-v4.txt`: hapA section |
| nsp13/14 essentiality | I | Virology literature (helicase/exonuclease) |
| Mpox outbreak (>87K) | L | WHO Disease Outbreak News [11] |
| Mpox dx priority | I | WHO priority + 2022 outbreak [11]; §3.6 |

### S3. Scanner Source and Methodology

**v4 scanner (current):** `loom pubmed-scan` (Rust CLI, `src/dna/pubmed_scan.rs`).
Four strategies per claim: PubMed narrow, PubMed broad (synonym-expanded with NCBI
ontology-derived aliases), local abstract corpus co-occurrence scan, ontology synonym
cross-check. Europe PMC as secondary source. Loads 8,283 NCBI gene annotations from
`web/data/ontology-enrichment.json` to auto-generate 51 synonym expansions.

**v3 scanner (DEPRECATED):** `scripts/pubmed_scan_v3.py` —
908 lines. Three strategies: PubMed narrow, PubMed broad (manual
synonyms only), local corpus co-occurrence. No ontology
integration (`ONTOLOGY_FILE` defined but never loaded).
Superseded by v4.

Full v4 results: `web/data/pubmed-scan-v4.json`

Example entry structure:
```json
{
  "pathogen": "mpox",
  "item": "gene:A33R",
  "confidence": "CONFIRMED",
  "evidence": {
    "pm_narrow": 0,
    "pm_broad": 0,
    "corpus_hits": 0,
    "epmc": 0
  },
  "synonym_source": "manual+ontology",
  "query_strings": {
    "narrow": "\"Monkeypox\"[Ti/Ab] AND CRISPR[Ti/Ab] AND ...",
    "broad": "(Monkeypox[Ti/Ab] OR MPXV[Ti/Ab]) AND ...",
    "epmc_query": "Monkeypox AND CRISPR AND A33R"
  },
  "scan_date": "2026-03"
}
```

**Key design decisions:**
- Europe PMC organism clause uses primary scientific name only (e.g., "Monkeypox" not "MPXV") to reduce false positives in full-text search
- Corpus cached to `target/pubmed_cache_v4/` to avoid redundant API calls on re-runs
- Rate limits: 0.38 s/request without API key; 0.12 s/request with `$NCBI_API_KEY`
- CONFIRMED threshold: PubMed broad = 0 AND corpus = 0 AND EuropePMC ≤ 30 (source: `src/dna/pubmed_scan.rs`)
- UNCERTAIN: PubMed broad > 0, OR corpus hit, OR EuropePMC > 30
- FALSE: PubMed broad > 5 (work clearly exists)
- Ontology synonyms auto-merged: manual takes priority, ontology fills gaps; each result tagged with source

### S4. Published CRISPR Guide Sequence Database

**Table S4.** Published SARS-CoV-2 CRISPR guide sequences compiled for sequence-level
novelty verification (§3.1). 83 unique sequences: 22 hand-curated from 15 papers
(Tier 0) plus 61 extracted from Europe PMC full texts (Tier 1, 66 papers processed).
Source: `published_guides_comprehensive.json`
(under `data/crispr_guides/`),
generated by `compile_published_guides.py`.

| Paper | PMID | Cas | Target Genes | # Guides | Use |
|:---|:---:|:---|:---|:---:|:---:|
| Broughton 2020 (DETECTR) | 32300245 | Cas12a | E | 1 | Dx |
| Patchsung 2020 (SHERLOCK) | 32826945 | Cas13a | ORF1ab | 1 | Dx |
| Abbott 2020 (PACMAN) | 32353252 | Cas13d | N, RdRp | 4 | Tx |
| Freije 2019 (CARVER) | 31607545 | Cas13a | conserved | 2 | Tx |
| Joung 2020 (STOPCovid) | 32941093 | Cas12b | N | 1 | Dx |
| Fozouni 2021 | 33357727 | Cas13a | N | 1 | Dx |
| Ackerman 2020 (CARMEN) | 32532005 | Cas13a | N | 1 | Dx |
| Ding 2020 | 32999292 | Cas12a | N | 2 | Dx |
| Huang 2020 | 32534345 | Cas12a | E | 1 | Dx |
| Guo 2020 (CRISPR-COVID) | 33064890 | Cas12a | N | 1 | Dx |
| Rauch 2021 (SENSR) | 32938741 | Cas13a | S | 1 | Dx |
| Azhar 2021 (iSCAN) | 33289377 | FnCas9 | N, RdRp | 2 | Dx |
| Ali 2020 (AIOD-CRISPR) | 33188177 | Cas12b | N | 1 | Dx |
| Lucia 2020 | — | Cas12a | N, S | 2 | Dx |
| Zhang 2020 (protocol) | — | Cas13a | S | 1 | Dx |

Dx = diagnostic; Tx = therapeutic/antiviral. Some sequences appear in multiple papers
(e.g., DETECTR_N and CARMEN share the same N-gene spacer); deduplication reduces 26
Tier 0 curated entries to 22 unique sequences. An additional 61 unique sequences were
extracted from Europe PMC full texts (Tier 1), yielding 83 total published guide sequences.
All 83 were compared against each of the 52 novel targets on both strands; the minimum
Hamming distance observed was 8 (see §3.1).

### S5. Machine-Readable Gap Evidence

Primary gap evidence artifacts (all under `web/data/`):

- `pubmed-scan-v4.json` — full scanner output, 12 pathogens,
  with synonym provenance
- `gap-audit-v4.txt` — human-readable audit trail
  with per-gene query strings and counts
- `pubmed-scan-results.json` — website-facing summary
  (confirmed gap filter in `crispr-app.js`)

Replication: any [L] claim in Table S2 can be verified by
running the logged query string directly in PubMed or Europe PMC.
Results should match within ±2 papers (publication lag).

### S6. Antiviral CRISPR Guide Conservation (Non-Diagnostic)

Conservation of published antiviral (non-diagnostic) CRISPR guide sequences targeting
SARS-CoV-2, measured against 9,193,298 NCBI genomic sequences. These guides were designed
for therapeutic viral RNA degradation, not diagnostic detection. Their conservation
measurements are reported separately from the diagnostic guides in Table 3.

| Guide | System | Gene | Cons. | Note |
|:---|:---|:---|:---:|:---|
| PACMAN\_N1 | Cas13 | N | 97.63% | Robust |
| PACMAN\_N2 | Cas13 | N | 98.32% | Robust |
| PACMAN\_RdRp1 | Cas13 | RdRp | 0.00% | **Failed** |
| PACMAN\_RdRp2 | Cas13 | RdRp | 96.24% | Drift |
| CARVER\_con1 [6] | Cas13 | Conserved | <0.01% | **Failed** |
| CARVER\_con2 [6] | Cas13 | Conserved | 0.00% | **Failed** |

PACMAN guides from Abbott 2020 [3] (*Cell*); CARVER guides from Freije 2019 [6]
(*Mol Cell*). PACMAN is a Cas13-based antiviral strategy designed to degrade viral
RNA inside infected cells, not to detect the presence of virus. CARVER (Freije 2019,
*Mol Cell*) is similarly an antiviral RNA-targeting system. Neither was designed as a
diagnostic; their conservation measurements reflect guide stability for sustained antiviral
activity, not diagnostic probe selection. The observation that PACMAN_RdRp1 (0.00%) and
CARVER guides (< 0.01% / 0.00%) have degraded to near-zero conservation is scientifically
significant for the antiviral field. **[C]**

### S7. Patent Landscape Scan (Preliminary)

**Important limitations.** This section provides a preliminary landscape survey, not
a freedom-to-operate analysis. The methodology has known weaknesses: (1) patent family
counts capture full-text co-mention of "CRISPR" and pathogen names—the majority are
broad platform patents listing multiple pathogens as examples, not sequence-specific
claims; (2) the sequence-level check covered only 10 of 52 guides; (3) we cannot assess
pending applications or international filings not indexed by Google Patents. This section
is included for completeness but should not be relied upon for IP decision-making.
Users requiring freedom-to-operate assessment should consult patent counsel.

A survey of Google Patents (patents.google.com) was conducted on 8 March 2026 to assess
the CRISPR patent landscape across all 12 target pathogens. Queries were expanded with
ontology-backed synonym sets (Disease Ontology, MONDO, NCBI Taxonomy) sourced from
`web/data/ontology-enrichment.json` to ensure comprehensive coverage. Ambiguous
abbreviations were excluded where they introduce cross-domain noise: bare "RSV" (also
matches Rous sarcoma virus) and bare "TB" (matches terabyte). Results are deduplicated
by patent family.

**Table S7a.** CRISPR patent family counts by pathogen (ontology-expanded queries).

| Pathogen | Synonyms searched | Families |
|:---|:---|---:|
| HIV | HIV, human immunodeficiency virus | 67,415 |
| Hepatitis B | Hepatitis B, HBV, serum hepatitis | 25,840 |
| Tuberculosis | tuberculosis, *M. tuberculosis*, Kochs disease | 15,656 |
| SARS-CoV-2 | SARS-CoV-2, COVID-19, 2019-nCoV | 13,040 |
| Cholera | *Vibrio cholerae*, cholera | 9,696 |
| MERS | MERS, MERS-CoV, camel flu | 9,425 |
| RSV | respiratory syncytial virus, HRSV$^{\dagger}$ | 8,751 |
| Dengue | Dengue, DENV, dengue disease | 8,501 |
| Ebola | Ebola, *Zaire ebolavirus*, EHF | 8,332 |
| Influenza A | Influenza A, Influenza A virus | 6,056 |
| Zika | Zika, ZikV, Zika fever | 5,676 |
| Mpox | Mpox, Monkeypox$^{\ddagger}$ | 332 |

$^{\dagger}$Bare "RSV" excluded (ambiguous: also matches Rous sarcoma virus).
$^{\ddagger}$Previous count (12,025) was inflated by an operator-precedence
bug; corrected with explicit Boolean grouping.
All queries use the form `CRISPR AND (synonym1 OR synonym2 OR ...)`.
Full query strings are archived in `patent_scan_report.json`.

Counts represent patent families mentioning both CRISPR and the pathogen name(s) in full text.
The majority are broad CRISPR platform patents (e.g., Broad Institute SHERLOCK/DETECTR
families) that list multiple pathogens as examples rather than claiming pathogen-specific
guide sequences. HIV and HBV dominate due to extensive CRISPR therapeutic programs
(provirus excision, cccDNA elimination). Mpox has the fewest CRISPR patent families (332),
consistent with its recent emergence as a global health priority.

**Table S7b.** SARS-CoV-2 gene-level patent coverage.

| Gene | Families | Note |
|:---|---:|:---|
| Nucleocapsid (N) | 1,092 | Most-patented; primary diagnostic target |
| Membrane protein (M) | 610 | |
| RdRp (nsp12) | 188 | |
| ORF1a (nsp1–11) | 90 | Our targets concentrate here |
| ORF3a | 88 | Our targets concentrate here |
| ORF8 | 88 | Our targets concentrate here |
| ORF10 | 51 | Our targets concentrate here |

Our 52 novel targets focus on ORF1a, ORF3a, ORF8, ORF10, and intergenic regions — genes
with 51–90 patent families each, versus 1,092 for the heavily patented nucleocapsid gene.
All 52 guide sequences are ≥8 mismatches from every published diagnostic guide in our
83-guide reference set (Table S4), providing additional indirect evidence of distinctness
from sequences likely to appear in patent claims.

**Note on limitations.** Patent family counts reflect full-text co-mention of "CRISPR" and
pathogen names; the majority are broad platform patents listing multiple pathogens as
examples. Gene-level counts (Table S7b) similarly capture co-mention, not sequence-specific
claims. This landscape survey characterizes patent density around our target genes but
does not constitute a freedom-to-operate analysis.

**Sequence-level patent check (preliminary).** A direct nucleotide search of 10
representative guide sequences in Google Patents returned zero patent-claimed matches for
9 of 10. The one flagged sequence — position 26513 (ORF3a–M intergenic, rank 29) —
appears in US20240156948A1 (Excepgen Inc.), which deposits the complete SARS-CoV-2 genome
as a sequence listing rather than claiming specific guide sequences. Combined with the
≥8-mismatch distance from all 83 published CRISPR guides (Table S4), this suggests low
overlap with patent-claimed sequences in our target regions (ORF1a/b, accessory ORFs). The
remaining 42 sequences were not checked at the nucleotide level; a comprehensive patent
sequence analysis would require professional patent counsel. The full scan artifact is
archived at `data/crispr_guides/patent_scan_report.json`. **[C/I]**

### S8. Animal Host Off-Target Screening

**Table S8.** Off-target hits (≤3 mismatches) for all 52 novel SARS-CoV-2 CRISPR guide
sequences screened against seven host reference genomes. Scanning used pigeonhole-seeded
Hamming-distance matching (4×5-nt seeds, Rayon-parallel over contigs).

| Host Species | Assembly | FASTA Size | Contigs | Hits (≤3 mm) | Runtime |
|:---|:---|:---:|:---:|:---:|:---:|
| Human | GRCh38 | 3.2 Gbp | 455 | 28,722 | 66.4 s |
| Pig | Sscrofa11.1 | 2.4 Gbp | 613 | 0 | 3.9 s |
| Egyptian Fruit Bat | mRouAeg1 | 1.8 Gbp | 29 | 0 | 2.9 s |
| Chicken | GRCg7b | 1.0 Gbp | 214 | 0 | 1.6 s |
| Cow | ARS-UCD1.3 | 2.6 Gbp | 1,958 | 0 | 4.3 s |
| Dromedary Camel | CamDro3 | 2.0 Gbp | 21,070 | 0 | 3.5 s |
| Mouse | GRCm39 | 2.6 Gbp | 61 | 0 | 4.1 s |

All 52 guides are completely clean (zero hits at any mismatch level) in all six
non-human host genomes. The 28,722 human hits include 0 exact matches, 58 one-mismatch
positions across 22 guides, 1,498 two-mismatch, and 27,166 three-mismatch positions.

Artifacts: `data/crispr_guides/offtargets_novel52_{species}_3mm.csv` for each species.

### S9. Pan-Pathogen Host Cross-Reactivity

**Table S9a.** Exact-match cross-reactivity of top 100 NGG-PAM (SpCas9) guides per pathogen
against seven animal host genomes. Generated March 5, 2026.

| Pathogen | Guides Tested | Host-Specific | Cross-Reactive | Specificity |
|:---|:---:|:---:|:---:|:---:|
| Cholera | 100 | 100 | 0 | 100% |
| Dengue | 100 | 100 | 0 | 100% |
| Ebola | 100 | 100 | 0 | 100% |
| Hepatitis B | 100 | 100 | 0 | 100% |
| HIV-1 | 100 | 100 | 0 | 100% |
| Influenza A | 100 | 100 | 0 | 100% |
| MERS | 100 | 100 | 0 | 100% |
| Mpox | 100 | 99 | 1 | 99% |
| RSV | 100 | 100 | 0 | 100% |
| *M. tuberculosis* | 100 | 98 | 2 | 98% |
| Zika | 100 | 100 | 0 | 100% |
| **Total** | **1,100** | **1,097** | **3** | **99.7%** |

Cross-reactive sequences: TB `GGCGGGGCCGGCGGGGCCGGCGG` (cow ×1),
TB `GCGGGGCCGGCGGGGCCGGCGGG` (cow ×1, camel ×1),
Mpox `GATGGATATGATGGATATGATGG` (mouse ×3). All are low-complexity repeats.

Artifact: `web/data/cross-reactivity.json`.

**Table S9b.** Exact-match cross-reactivity of top 100 TTTV-PAM (Cas12a) guides per pathogen
against seven animal host genomes (TTTV only; TTTT-prefixed sequences excluded as
non-functional Cas12a PAMs). Generated March 7, 2026.

| Pathogen | Guides Tested | Host-Specific | Cross-Reactive | Specificity |
|:---|:---:|:---:|:---:|:---:|
| Cholera | 100 | 100 | 0 | 100% |
| Dengue | 76 | 76 | 0 | 100% |
| Ebola | 100 | 100 | 0 | 100% |
| Hepatitis B | 100 | 100 | 0 | 100% |
| HIV-1 | 100 | 100 | 0 | 100% |
| Influenza A | 100 | 100 | 0 | 100% |
| MERS | 100 | 100 | 0 | 100% |
| Mpox | 100 | 98 | 2 | 98% |
| RSV | 100 | 100 | 0 | 100% |
| *M. tuberculosis* | 27 | 27 | 0 | 100% |
| Zika | 69 | 69 | 0 | 100% |
| **Total** | **972** | **970** | **2** | **99.8%** |

Cross-reactive sequences (both Mpox, human genome only):
`TTTATATTTTATATTTTATTTTA` (human ×6),
`TTTATTTTATATTTTATATTTTA` (human ×6).
Both are AT-rich low-complexity repeats.

Artifact: `web/data/cross-reactivity-cas12a.json`.

