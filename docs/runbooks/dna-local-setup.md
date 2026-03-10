# DNA Local Setup Runbook

Reproducible commands for running the full DNA ingestion and search pipeline
on a Mac Studio M3 Ultra (96 GB RAM, macOS).  No network access required after
data download.

## Prerequisites

```bash
# Build the loom binary (release for benchmarks, debug for dev)
cargo build --release

# Verify the build
./target/release/loom --help
```

## Step 1 — Prepare FASTA/FASTQ data

Place your FASTA or FASTQ files in a local directory, e.g. `data/dna/raw/`.
Supported extensions: `.fasta`, `.fa`, `.fna`, `.ffn`, `.fastq`, `.fq`

For a quick sanity test, create a synthetic file:

```bash
mkdir -p data/dna/raw
cat > data/dna/raw/sample.fasta << 'EOF'
>seq1 Homo sapiens chromosome 1 fragment
ACGTACGTACGTNNNACGTRYSWKMBDHV
>seq2 Pathogen reference
GCTAGCTAGCTAGCTAACGTTTTTAAAA
EOF
```

## Step 2 — Generate DNA manifest (checksums + metadata)

```bash
./target/release/loom dna-manifest \
  --input  data/dna/raw \
  --output data/dna/manifest.json \
  --split  train \
  --source-url "https://your-dataset-source.org" \
  --license "CC-BY-4.0"

# Inspect
cat data/dna/manifest.json
```

## Step 3 — Build normalized JSONL corpus

```bash
./target/release/loom dna-corpus \
  --input   data/dna/raw \
  --output  data/dna/corpus.jsonl \
  --min-len 32

# Inspect first record
head -1 data/dna/corpus.jsonl | python3 -m json.tool
```

Expected JSONL record format:

```json
{
  "id": "seq1",
  "sequence": "ACGTACGTACGTNNNNNNNNNN",
  "source_path": "data/dna/raw/sample.fasta",
  "length": 22
}
```

## Step 4 — Build FM-index from FASTA files

```bash
./target/release/loom index \
  --path   data/dna/raw \
  --output data/dna/dna.idx

# Index size sanity check
ls -lh data/dna/dna.idx
```

## Step 5 — Search for a DNA pattern

```bash
# Count occurrences
./target/release/loom search ACGTACGT \
  --index data/dna/dna.idx \
  --count

# Show locations with context
./target/release/loom search ACGTACGT \
  --index   data/dna/dna.idx \
  --context 2
```

## Step 6 — Generate BWT training samples (Python bridge)

```bash
# Requires: pip install -e python/ (or PYTHONPATH=python)
python3 - << 'EOF'
from pathlib import Path
from loom.training.dna_dataset import build_samples_from_dna_jsonl, deterministic_split

samples = build_samples_from_dna_jsonl(
    jsonl_path=Path("data/dna/corpus.jsonl"),
    window_size=128,
    max_samples=5_000,
    seed=42,
)
train, val, test = deterministic_split(samples)
print(f"train={len(train)}  val={len(val)}  test={len(test)}")
EOF
```

## Quick end-to-end smoke test

```bash
cargo test --test integration test_dna_corpus_and_index_search_end_to_end
```

## Troubleshooting

| Symptom | Fix |
|---|---|
| `No FASTA/FASTQ files found` | Check file extensions match `.fasta/.fa/.fna/.ffn/.fastq/.fq` |
| `sequence length < min_len` | Lower `--min-len` or check file contains real sequences |
| Search returns 0 results | Confirm FASTA and index were built from the same directory |
| Python `ModuleNotFoundError` | `export PYTHONPATH=python` before running Python commands |
