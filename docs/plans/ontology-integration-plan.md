# Ontology Integration Plan — LOOM CRISPR Database

**Date:** March 4, 2026
**Status:** Draft / Research

---

## Current State

LOOM's CRISPR target database has **minimal metadata per pathogen**: a display name, virus family, genome count, and flat target records (sequence, position, PAM, conservation score). There are no standardized identifiers linking pathogens to diseases, symptoms, host organisms, gene functions, or clinical classification systems.

The existing `ontology.rs` module in brenda is a **co-occurrence graph for search term expansion** — not a biomedical ontology. It's unrelated to this effort.

---

## Why Ontologies Matter for This Project

Linking LOOM's 140,000 CRISPR targets to established biomedical ontologies would:

1. **Make the database interoperable** — researchers can query by disease (not just pathogen name), link to clinical workflows, and connect to other databases
2. **Enable richer research gap analysis** — cross-reference gaps against WHO priority pathogens, outbreak status, diagnostic availability
3. **Strengthen the papers** — reviewers expect standardized identifiers (NCBI TaxID at minimum); ontology-backed metadata upgrades credibility from "search tool" to "research resource"
4. **Open new queries** — "show me all CRISPR targets for WHO priority pathogens with no existing point-of-care diagnostic" requires linking disease ontology → pathogen → our targets

---

## Relevant Ontologies (Ranked by Relevance)

### Tier 1 — Direct & High-Value

| Ontology | What It Covers | Relation to LOOM | Effort |
|----------|---------------|-----------------|--------|
| **NCBI Taxonomy** | Organism classification (TaxID) | Every pathogen already has a TaxID. Adding it links our data to GenBank, RefSeq, UniProt, and every major bioinformatics DB. | **Trivial** — just add 14 TaxIDs to the JSON |
| **Disease Ontology (DO/DOID)** | Human diseases | Maps each pathogen → disease(s) it causes. Cholera → DOID:1498, COVID-19 → DOID:0080600, Ebola → DOID:4325. Connects our targets to clinical context. | **Low** — 14 manual mappings |
| **Sequence Ontology (SO)** | Genomic features | Classifies what our targets *are*: guide RNA (SO:0001998), PAM site, protospacer. Makes target records self-describing. | **Low** — add 3-4 SO terms to schema |
| **MONDO** | Unified disease ontology | Superset of DOID + OMIM + Orphanet. If we pick one disease ID system, MONDO is the most interoperable. | **Low** — 14 mappings, overlaps with DO |

### Tier 2 — Valuable for Enrichment

| Ontology | What It Covers | Relation to LOOM | Effort |
|----------|---------------|-----------------|--------|
| **ICD-11** | WHO disease classification | The global standard clinicians use. Adding ICD-11 codes means our database speaks the same language as hospitals and surveillance systems. | **Low** — 14 codes |
| **Infectious Disease Ontology (IDO)** | Infectious disease processes | Formal model of infection, transmission, host response. Could classify *why* each pathogen matters (transmission route, mortality, outbreak potential). | **Medium** — requires selecting relevant IDO classes |
| **VIDO** (Virus Infectious Disease Ontology) | Virus-specific extension of IDO | Covers viral replication, tropism, immune evasion — directly relevant for 12/14 of our pathogens (all except cholera and TB). | **Medium** |
| **Gene Ontology (GO)** | Gene/protein function | If we annotate which *genes* our targets fall in (not just positions), GO terms describe their function (e.g., "RNA-dependent RNA polymerase activity"). | **Medium-High** — requires gene-level mapping |
| **Human Phenotype Ontology (HPO)** | Symptoms & clinical features | Links diseases to their symptoms. "What CRISPR targets exist for diseases causing hemorrhagic fever?" becomes answerable. | **Medium** — needs disease → phenotype mapping |

### Tier 3 — Future / Specialized

| Ontology | What It Covers | Relation to LOOM | Effort |
|----------|---------------|-----------------|--------|
| **GENEPIO** | Genomic epidemiology | Sample metadata standards for outbreak genomics. Relevant if LOOM targets are used in real surveillance. | Medium |
| **Vaccine Ontology (VO)** | Vaccines | Links pathogens to existing/in-development vaccines. "Targets for pathogens with no vaccine" is a powerful query. | Medium |
| **PHI-base** | Pathogen-host interactions | Which host genes interact with which pathogen genes. Relevant for off-target analysis context. | High |
| **ChEBI** | Chemical entities | Relevant only if we start annotating CRISPR reagents (Cas proteins, tracrRNAs). | Low priority |
| **SNOMED CT** | Clinical terminology | Broader than ICD but license-restricted in some jurisdictions. ICD-11 is the open equivalent. | Low priority |

---

## Concrete Mapping: LOOM Pathogens → Ontology IDs

| LOOM Key | Species | NCBI TaxID | DOID | MONDO | ICD-11 |
|----------|---------|-----------|------|-------|--------|
| cholera | *V. cholerae* | 666 | DOID:1498 | MONDO:0015904 | 1A00 |
| dengue | Dengue virus | 12637 | DOID:11205 | MONDO:0005502 | 1D20 |
| ebola | Ebola virus | 186538 | DOID:4325 | MONDO:0005737 | 1D60 |
| hepatitis-b | Hepatitis B virus | 10407 | DOID:2043 | MONDO:0005344 | 1E50.1 |
| hiv-1 | HIV-1 | 11676 | DOID:526 | MONDO:0005109 | 1C60 |
| influenza-a | Influenza A virus | 11320 | DOID:8469 | MONDO:0005812 | 1E30 |
| mers | MERS-CoV | 1335626 | DOID:0080642 | MONDO:0100313 | 1D65 |
| mpox | Monkeypox virus | 10244 | DOID:12594 | MONDO:0002594 | 1E70 |
| rsv | RSV | 12814 | DOID:0050140 | MONDO:0005603 | 1E31 |
| sars-cov-2 | SARS-CoV-2 | 2697049 | DOID:0080600 | MONDO:0100096 | RA01 |
| tuberculosis | *M. tuberculosis* | 1773 | DOID:399 | MONDO:0018076 | 1B10 |
| zika | Zika virus | 64320 | DOID:0060478 | MONDO:0018087 | 1D4D |
| refseq-viral | Mixed | — | — | — | — |
| human-grch38 | *H. sapiens* | 9606 | — | — | — |

> **Note:** TaxIDs, DOIDs, and ICD-11 codes above need verification against current releases before publishing. The MONDO IDs in particular should be confirmed via the MONDO browser (https://mondo.monarchinitiative.org/).

---

## What This Unlocks

### Immediate (Tier 1 integration)

1. **Cross-database linking** — a researcher finds a target in LOOM, clicks the TaxID, lands in NCBI with full genome context
2. **Machine-readable disease labels** — DOID/MONDO IDs let other tools (e.g., Open Targets, DisGeNET) reference our data
3. **Proper SO annotation** — our JSON becomes self-documenting: each target is explicitly tagged as a CRISPR guide RNA candidate

### With Tier 2

4. **Symptom-based queries** — "show targets for hemorrhagic fever pathogens" (HPO → DOID → our pathogens)  
5. **Gene function context** — "this target sits in the RNA polymerase gene (GO:0003723)" gives researchers biological meaning
6. **Surveillance alignment** — ICD-11 codes connect our database to WHO/national disease surveillance systems
7. **WHO priority matching** — cross-reference our gap analysis against WHO priority pathogen lists programmatically

### With Tier 3

8. **Vaccine gap analysis** — "pathogens in our DB with no licensed vaccine" (from VO)
9. **Epidemiological metadata** — GENEPIO-compliant sample metadata for outbreak integration

---

## Implementation Plan

### Phase 1: Structural IDs (1-2 hours)

Add ontology identifiers to the CRISPR targets JSON schema:

```json
{
  "cholera": {
    "name": "V. cholerae",
    "family": "Vibrionaceae",
    "genomes": 382,
    "ncbi_taxid": 666,
    "disease_ontology_id": "DOID:1498",
    "mondo_id": "MONDO:0015904",
    "icd11": "1A00",
    "targets": [...]
  }
}
```

- Verify all IDs against current ontology releases
- Update `crispr-targets.json` and `animal-targets.json`
- Update the web UI to display/link these IDs

### Phase 2: Sequence Ontology for Targets (1 hour)

Add SO terms to the target schema:

```json
{
  "s": "GCCTCAAGAGGGACTGTCAACGC",
  "p": 3850358,
  "m": "NGG",
  "o": 12915,
  "so": "SO:0001998"  
}
```

This is probably overkill for the JSON (every target is the same type) — better as a schema-level declaration:
```json
{ "target_type": "SO:0001998", "target_type_label": "sgRNA" }
```

### Phase 3: Gene-Level Annotation (days, not hours)

Map target positions to gene annotations:
- Download GFF3 annotation files for each pathogen reference genome
- For each target, determine which gene it falls in
- Add gene name + GO terms to each target record
- This is the highest-value enrichment but requires real bioinformatics work

### Phase 4: Disease-Phenotype Links (optional)

- Pull HPO annotations for each disease
- Enable symptom-based browsing in the web UI
- Cross-reference with vaccine availability (VO) for gap analysis

---

## Relation to Existing Papers

| Paper | How Ontologies Help |
|-------|-------------------|
| **Paper 1: Pan-Pathogen Targets** | Reviewers expect NCBI TaxIDs at minimum. DOID/MONDO IDs position the database as interoperable with the OBO Foundry ecosystem. Gene-level GO annotation would be a strong differentiator. |
| **Paper 2: LOOM Tool Paper** | Mention SO terms in the data model section. Shows the tool speaks standard vocabularies. |
| **Paper 3: Verifiable AI Science** | Ontology IDs are machine-verifiable claims. A target linked to TaxID:666 + DOID:1498 is independently checkable — this directly supports the verifiability framework. |

---

## Recommendation

**Do Phase 1 now** (add TaxID, DOID, MONDO, ICD-11 to the JSON). It's 1-2 hours of work, dramatically improves credibility with reviewers, and makes the database interoperable. The mapping table above is a starting point — verify each ID before publishing.

Phase 3 (gene-level GO annotation) is the richest improvement but is real work. Defer unless paper reviewers request it.

---

## Key Ontology Resources

- NCBI Taxonomy: https://www.ncbi.nlm.nih.gov/taxonomy
- Disease Ontology: https://disease-ontology.org/
- MONDO: https://mondo.monarchinitiative.org/
- Sequence Ontology: http://www.sequenceontology.org/
- Gene Ontology: http://geneontology.org/
- ICD-11: https://icd.who.int/
- OBO Foundry (hub for all bio-ontologies): http://obofoundry.org/
- IDO: https://infectiousdiseaseontology.org/
- VIDO: part of IDO family
