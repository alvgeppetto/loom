# Runbook: Same-Corpus CRISPR Revalidation

Re-run the CRISPR guide conservation scans against the identical FASTA corpus
to confirm denominator and hit-rate consistency.

## Prerequisites

- LOOM binary built: `cargo build --release`
- FASTA corpus mounted at the expected path (see below)
- ~40 min wall-clock for both scans on Mac Studio M3 Ultra

## Corpus

```
/Volumes/T9/loom-openscience/raw/sars-cov-2/ncbi_dataset/data/genomic.fna
```

Canonical genome count: **9,193,298** (confirmed 2026-03-07).

## Scans

### Novel guides (52)

```bash
tmux new-session -d -s rerun-novel \
  "./target/release/loom guide-conservation \
     --guides data/crispr_guides/novel_52_guides.tsv \
     --fasta /Volumes/T9/loom-openscience/raw/sars-cov-2/ncbi_dataset/data/genomic.fna \
     --output data/crispr_guides/novel_targets_same_corpus.json \
     2>&1 | tee logs/same-corpus-rerun-novel.log"
```

Expected: ~23 min, ~6,500 genomes/s, denominator = 9,193,298.

### Published guides (10)

```bash
tmux new-session -d -s rerun-published \
  "./target/release/loom guide-conservation \
     --guides data/crispr_guides/published_guides.tsv \
     --fasta /Volumes/T9/loom-openscience/raw/sars-cov-2/ncbi_dataset/data/genomic.fna \
     --output data/crispr_guides/published_guides_same_corpus.json \
     2>&1 | tee logs/same-corpus-rerun-published.log"
```

Expected: ~14 min, ~11,000 genomes/s, denominator = 9,193,298.

## Verification

After both scans finish:

```bash
# Confirm denominator in logs
grep -o 'Done: [0-9]* genomes' logs/same-corpus-rerun-*.log

# Confirm denominator in JSON output
jq '.[0].total_genomes' data/crispr_guides/novel_targets_same_corpus.json
jq '.[0].total_genomes' data/crispr_guides/published_guides_same_corpus.json
```

Both must report **9,193,298**.

## Build guard

`scripts/build-preprints.sh` contains a `check_regressions()` gate that blocks
PDF builds if the merged manuscript contains the stale denominator (`9,419,528`).

## Cross-paper audit

```bash
bash scripts/audit-paper-consistency.sh
```

The audit script checks all papers for canonical vs legacy denominator tokens.
Only `pangenomic-crispr-targets-merged.md` should be in `publications/papers/`.
All others are archived in `publications/archive/papers/`.

## Reference

| Artifact | Path |
|---|---|
| Novel scan output | `data/crispr_guides/novel_targets_same_corpus.json` |
| Published scan output | `data/crispr_guides/published_guides_same_corpus.json` |
| Scan logs (2026-03-07) | `logs/same-corpus-rerun-20260307/` |
| Merged manuscript | `publications/papers/pangenomic-crispr-targets-merged.md` |
| Build script | `scripts/build-preprints.sh` |
| Audit script | `scripts/audit-paper-consistency.sh` |
