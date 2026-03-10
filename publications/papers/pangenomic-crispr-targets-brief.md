# Pangenomic CRISPR Target Discovery: 52 Novel Sites in SARS-CoV-2 Replication Machinery and a Multi-Pathogen Conserved-Target Database

**Alvaro Videla Godoy**^1^

^1^ Independent Researcher

**Correspondence:** videlalvaro@gmail.com

**Preprint Brief — March 2026**

---

## Abstract

CRISPR diagnostics overwhelmingly target immune-exposed structural genes (Spike, Nucleocapsid), which accumulate escape mutations that degrade assay sensitivity. We scanned 9,193,298 SARS-CoV-2 genomes using an FM-index k-mer engine and identified 52 novel CRISPR target sites in replication machinery — regions under strong purifying selection and absent from all published guide sets. The top 8 candidates retain 98.9–99.2% conservation, comparable to the best published diagnostic guide (DETECTR_E: 98.34%) on an identical corpus. Across 12 pathogens, we provide a ~130,000-target conservation-ranked database and confirm 6 literature gaps (5 Mpox genes, 1 cholera gene) via four-strategy verification with ontology-derived synonym expansion. Conservation is necessary but not sufficient for guide performance: all 52 candidates require wet-lab validation. Data, code, and an in-browser search tool are publicly available.

---

## Introduction

CRISPR-based detection platforms — SHERLOCK [1], DETECTR [2], PACMAN [3] — have proven rapid pathogen detection, but guide RNA design follows historical RT-PCR conventions: structural and surface genes predominate. This creates two problems. First, immune-exposed genes accumulate mutations under selection pressure on the pathogen, degrading diagnostic sensitivity — as demonstrated when SHERLOCK_S targeting the SARS-CoV-2 Spike degraded to 63% conservation across 9.2 million genomes [4]. Second, essential replication enzymes under strong purifying selection remain largely unexamined as diagnostic targets.

We address both problems. For SARS-CoV-2, we report 52 diagnostic-novel CRISPR target sites in replication machinery (ORF1a, nsp13/Helicase, nsp14/ExoN, and accessory ORFs), validated at full pangenomic scale. For 11 additional pathogens, we provide a conservation-ranked target database and a multi-source literature gap analysis. The pangenomic pipeline developed for SARS-CoV-2 generalizes to any pathogen, and the multi-pathogen database validates that conservation-ranking methodology scales across diverse genome architectures.

---

## Results

### SARS-CoV-2: 52 Novel Targets in Replication Machinery

From 29,640 candidate windows, 183 carry compatible PAM sites across 89 genomic regions. After excluding regions overlapping published diagnostic targets (N, E, S, nsp12/RdRp), 52 truly novel regions remain — all ≥8 mismatches from any of 83 published guide sequences. Targets concentrate in ORF1a/b replication machinery under purifying selection.

**Table 1.** Eight highest-conservation novel candidates (same-corpus validation, 9,193,298 genomes).

| Rank | Pos. | Sequence (5'→3') | Protein | Cons. | PAM | GC% | MFE |
|:---:|:---:|:---|:---|:---:|:---|:---:|:---:|
| 1 | 10456 | ACTATTAAGGGTTCATTCCT | ORF1a | 99.2% | Cas12a | 35 | −2.0 |
| 2 | 17431 | CTGCTCAATTACCTGCACCA | nsp13 | 99.1% | SpCas9 | 50 | −1.3 |
| 3 | 17561 | TGCTGAAATTGTTGACACTG | nsp13 | 99.1% | SpCas9 | 40 | 0.0 |
| 4 | 19218 | ATTGTTTGTAGATTTGACAC | nsp14 | 99.1% | SpCas9 | 30 | −0.3 |
| 5 | 16813 | GAGAGTACACCTTTGAAAAA | nsp13 | 99.1% | SpCas9 | 35 | 0.0 |
| 6 | 10531 | TGTTACATGCACCATATGGA | ORF1a | 99.0% | Cas12a | 40 | −0.3 |
| 7 | 17078 | TGGTACTGGTAAGAGTCATT | nsp13 | 99.0% | SpCas9 | 40 | −0.4 |
| 8 | 12888 | ACCTTGTAGGTTTGTTACAG | ORF1a | 98.9% | SpCas9 | 40 | −2.1 |

MFE = minimum free energy (kcal/mol, ViennaRNA RNAfold). Cons. = conservation across 9,193,298 SARS-CoV-2 genomes, same corpus and tool as published guide benchmarks (Table 2).

The top 8 are within 1 percentage point of the best published diagnostic guide. All 8 show zero exact matches in GRCh38 (3.2 Gbp) and zero hits at any mismatch level across six non-human host genomes. Of the 52 candidates: 29 fall outside the 40–70% GC window (mean 36.7%, consistent with SARS-CoV-2's ~38% GC content); 3 have MFE ≤ −5 kcal/mol; 0 of 52 have seed-region VOC mutations.

**Published guide comparison.** Same-corpus validation against the three published diagnostic guides:

**Table 2.** Published CRISPR diagnostic guide conservation (9,193,298 genomes, same corpus as Table 1).

| Guide | Gene | Cons. |
|:---|:---|:---:|
| DETECTR_E (Broughton 2020) | E | 98.34% |
| SHERLOCK_Orf1ab (Patchsung 2020) | ORF1ab | 97.55% |
| DETECTR_N (Broughton 2020) | N | 94.21% |

SHERLOCK_S (Spike) retains only 63% conservation — demonstrating the failure mode our replication-machinery targets avoid.

### Multi-Pathogen Target Database

Across 12 pathogens plus a panviral reference, the database contains ~130,000 conservation-ranked CRISPR candidate sequences (top 10,000 per pathogen), of which 8,113 carry SpCas9 NGG PAMs.

**Table 3.** Pathogen genome collections and target counts.

| Pathogen | Genomes | Top-10K | NGG PAM |
|:---|---:|---:|---:|
| HIV-1 | 1,332,591 | 10,000 | 676 |
| *M. tuberculosis* | 498 | 10,000 | 1,132 |
| Influenza A | 1,555,840 | 10,000 | 662 |
| Mpox | 11,461 | 10,000 | 311 |
| *V. cholerae* | 382 | 10,000 | 705 |
| RSV | 64,802 | 10,000 | 253 |
| Dengue | 55,283 | 10,000 | 748 |
| Hepatitis B | 137,774 | 10,000 | 725 |
| Ebola | 3,663 | 10,000 | 442 |
| MERS-CoV | 1,812 | 10,000 | 446 |
| Zika | 2,711 | 10,000 | 904 |
| SARS-CoV-2 | 9,193,298 | 10,000 | 355 |
| RefSeq Viral | 19,019 | 10,000 | 754 |
| **Total** | | **130,000** | **8,113** |

Corpus sizes vary by four orders of magnitude; conservation metrics are not directly comparable across pathogens. Pan-pathogen cross-reactivity screening (2,072 top guides across both PAM types against 7 host genomes) confirmed 99.8% host specificity at exact-match resolution.

### Verified Literature Gaps

After four-strategy verification with ontology-derived synonyms, 6 confirmed gene-level gaps remain:

**Table 4.** Confirmed literature gaps (zero CRISPR publications across all evidence strategies, March 2026).

| Pathogen | Gene | Function |
|:---|:---|:---|
| Mpox | A33R | Immune evasion |
| Mpox | B5R | Membrane protein |
| Mpox | E8L | Unknown function |
| Mpox | J2R (OPG101) | Thymidine kinase |
| Mpox | A56R (OPG185) | Hemagglutinin |
| *V. cholerae* | hapA | Hemagglutinin protease |

These 6 genes have zero CRISPR publications across all four evidence strategies as of March 2026. For SARS-CoV-2, no confirmed gene-level gaps exist — but nsp13 (6 papers) and nsp14 (9 papers) have exclusively non-diagnostic CRISPR research, representing an unexploited diagnostic opportunity.

---

## Discussion

The 52 replication-machinery targets are the headline result: CRISPR guide candidates in regions under lethal purifying selection, absent from every published guide set, and conserved at 98.9–99.2% across 9.2 million genomes — matching the best published diagnostic guide on an identical corpus. Replication enzymes avoid the variant-driven degradation that crippled Spike-targeting guides (SHERLOCK_S: 63%). Whether this conservation advantage translates to superior diagnostic performance remains an empirical question.

The accompanying ~130,000-target database across 12 pathogens provides a ready-to-use starting point for experimental guide design, queryable in-browser without software installation. The 6 confirmed gene gaps (5 Mpox, 1 cholera) mark actionable whitespace: conserved sequences with no published CRISPR diagnostic work.

A methodological note: literature gap claims require multi-source verification with ontology-derived synonym expansion. Single-query PubMed scans overestimate gaps by 5–8×; ontology integration catches nomenclature mismatches that manual curation misses.

---

## Limitations

- **Conservation ≠ efficacy.** All 52 candidates require wet-lab validation; k-mer conservation cannot predict cleavage efficiency, binding kinetics, or limit of detection.
- **In silico only.** No experimental validation was performed. 29 of 52 candidates fall outside the 40–70% GC window.
- **Corpus bias.** The 9.2M SARS-CoV-2 corpus is Omicron-weighted (~65%); multi-pathogen datasets are convenience samples (382 to 1.5M genomes), not pangenomic surveys.
- **Literature claims are time-dependent.** Gap claims carry the date March 2026 and specific query strings. New publications may have appeared since our scan.
- **Non-English and non-PubMed literature not fully covered.** Research indexed in national databases may describe work on these targets.

---

## Data Availability and Methods

All artifacts, source code, target database, and an in-browser CRISPR target search tool are publicly available at: **https://github.com/alvgeppetto/loom**

**Annotated walkthrough:** A plain-language guide to this brief is available at: <https://calm-mushroom-0185d800f.4.azurestaticapps.net/annotated-brief.html>

All scanning was performed using LOOM, an FM-index k-mer engine. SARS-CoV-2 candidates were identified by scoring all 29,640 overlapping 20-mers from Wuhan-Hu-1 (NC_045512.2) for conservation across 9,193,298 genomes (262 GB), filtered for PAM compatibility, screened for off-targets against 7 host genomes (≤3 mismatches), and verified for sequence-level novelty against 83 published guides (minimum Hamming distance ≥8). Literature gaps were classified by a four-strategy scanner with NCBI ontology-derived synonym expansion; a gene is confirmed as a gap only when all four strategies return zero results. Full methodology is described in the companion manuscript.

---

## References

[1] Gootenberg JS, Abudayyeh OO, Lee JW, et al. Nucleic acid detection with CRISPR-Cas13a/C2c2. *Science*. 2017;356(6336):438-442. DOI: 10.1126/science.aam9321

[2] Broughton JP, Deng X, Yu G, et al. CRISPR–Cas12-based detection of SARS-CoV-2. *Nat Biotechnol*. 2020;38(7):870-874. DOI: 10.1038/s41587-020-0513-4

[3] Abbott TR, Dhamdhere G, Liu Y, et al. Development of CRISPR as an Antiviral Strategy to Combat SARS-CoV-2 and Influenza. *Cell*. 2020;181(4):865-876.e12. DOI: 10.1016/j.cell.2020.04.020

[4] Patchsung M, Juntawong K, et al. Clinical validation of a Cas13-based assay for the detection of SARS-CoV-2 RNA. *Nat Biomed Eng*. 2020;4(12):1140-1149. DOI: 10.1038/s41551-020-00603-x
