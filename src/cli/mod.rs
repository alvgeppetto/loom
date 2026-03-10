//! CLI module - command-line interface

use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

/// Build-mode selector accepted by `loom index --build-mode`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum BuildModeArg {
    /// Load full corpus into RAM; fastest but requires ~6× corpus bytes of free RAM.
    InMemory,
    /// Use bounded-memory streaming path; safe for corpora larger than available RAM.
    Streaming,
    /// Choose automatically using the preflight RAM estimator (default).
    Auto,
}

#[derive(Parser)]
#[command(name = "loom")]
#[command(about = "BWT-augmented code search", long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Build an index from a directory
    Index {
        /// Path to the source directory
        #[arg(short, long)]
        path: PathBuf,

        /// File extensions to include (e.g., "erl,hrl")
        #[arg(short, long, default_value = "")]
        extensions: String,

        /// Output path for the index file
        #[arg(short, long, default_value = "loom.idx")]
        output: PathBuf,

        /// Also build and embed a vocabulary (VCAB) section
        #[arg(long)]
        vocab: bool,

        /// Encrypt the output index with AES-GCM-256 (writes <output>.enc)
        #[arg(long)]
        encrypt: bool,

        /// Environment variable holding the encryption password (used with --encrypt)
        #[arg(long, default_value = "LOOM_KEY")]
        key_env: String,

        /// Environment variable holding the 32-byte blind-index key (builds permuted corpus)
        #[arg(long)]
        blind_key_env: Option<String>,

        /// Build strategy: in-memory (fast, RAM-bound), streaming (low-memory), or auto (default)
        #[arg(long, default_value = "auto", value_name = "MODE")]
        build_mode: BuildModeArg,

        /// Run memory estimator preflight check and print JSON, then exit without building
        #[arg(long)]
        estimate_memory: bool,

        /// Split corpus into N shards (writes shard_NNNN.idx + shard_manifest.json to output dir)
        #[arg(long)]
        shards: Option<usize>,

        /// Target bytes per shard when using --shards (default: 64 MiB)
        #[arg(long, default_value = "67108864")]
        shard_size: usize,

    },

    /// Search the index for a pattern
    Search {
        /// The pattern to search for
        pattern: String,

        /// Path to the index file (.idx or .enc)
        #[arg(short, long, default_value = "loom.idx")]
        index: PathBuf,

        /// Number of context lines to show
        #[arg(short, long, default_value = "3")]
        context: usize,

        /// Only show count of matches
        #[arg(long)]
        count: bool,

        /// Environment variable holding the decryption password (for .enc files)
        #[arg(long, default_value = "LOOM_KEY")]
        key_env: String,

        /// Environment variable holding the 32-byte blind-index key (for .blind files)
        #[arg(long)]
        blind_key_env: Option<String>,
    },

    /// Search all shards in a manifest (federated fan-out + merge)
    ShardSearch {
        /// The pattern to search for
        pattern: String,

        /// Path to shard_manifest.json
        #[arg(short, long, default_value = "loom.shards/shard_manifest.json")]
        manifest: PathBuf,

        /// Context characters around each match
        #[arg(short, long, default_value = "80")]
        context: usize,

        /// Only show total count of matches across all shards
        #[arg(long)]
        count: bool,

        /// Maximum parallel worker threads (0 = one per shard)
        #[arg(long, default_value = "0")]
        workers: usize,
    },

    /// Ask a natural language question (requires LLM)
    Ask {
        /// The question to ask
        question: String,

        /// Path to the index file
        #[arg(short, long, default_value = "loom.idx")]
        index: PathBuf,
    },

    /// Look up vocabulary terms by prefix (or exact match with --exact)
    VocabLookup {
        /// Prefix to search for (or exact term with --exact)
        prefix: String,

        /// Path to the index file
        #[arg(short, long, default_value = "loom.idx")]
        index: PathBuf,

        /// Check exact existence instead of prefix search
        #[arg(long)]
        exact: bool,
    },

    /// Look up related terms via co-occurrence ontology
    Related {
        /// The term to look up
        term: String,

        /// Path to the index file
        #[arg(short, long, default_value = "loom.idx")]
        index: PathBuf,

        /// Maximum number of related terms to show
        #[arg(short, long, default_value = "10")]
        limit: usize,
    },

    /// Build a normalized DNA corpus JSONL from FASTA/FASTQ inputs
    DnaCorpus {
        /// Input FASTA/FASTQ file or directory
        #[arg(short, long)]
        input: PathBuf,

        /// Output JSONL file path
        #[arg(short, long, default_value = "dna_corpus.jsonl")]
        output: PathBuf,

        /// Minimum normalized sequence length to keep
        #[arg(long, default_value = "32")]
        min_len: usize,
    },

    /// Build a DNA dataset manifest with checksums and metadata
    DnaManifest {
        /// Input FASTA/FASTQ file or directory
        #[arg(short, long)]
        input: PathBuf,

        /// Output JSON manifest path
        #[arg(short, long, default_value = "dna_manifest.json")]
        output: PathBuf,

        /// Split tag to assign to all discovered files (train/val/test)
        #[arg(long, default_value = "train")]
        split: String,

        /// Source dataset URL to annotate each entry
        #[arg(long)]
        source_url: Option<String>,

        /// Dataset license label to annotate each entry
        #[arg(long)]
        license: Option<String>,
    },

    /// Benchmark DNA corpus parsing latency (cold, warm median, p95, p99).
    ///
    /// Parses the given FASTA/FASTQ file or directory repeatedly and reports
    /// latency statistics.  The first round is the cold measurement; all
    /// subsequent rounds feed the warm distribution.
    DnaLatency {
        /// Input FASTA/FASTQ file or directory to parse
        #[arg(short, long)]
        input: PathBuf,

        /// Minimum normalized sequence length to keep (mirrors dna-corpus)
        #[arg(long, default_value = "32")]
        min_len: usize,

        /// Total number of parse rounds (first = cold, rest = warm)
        #[arg(long, default_value = "20")]
        rounds: usize,

        /// Write JSON stats to this file instead of stdout
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Compare two DNA benchmark reports (e.g. plain-index vs transformer).
    ///
    /// Loads two JSON reports produced by `dna-benchmark`, computes metric
    /// deltas and KPI gate checks, and outputs a comparison report in JSON
    /// and optional Markdown.
    DnaCompare {
        /// JSON report for the baseline run (e.g. plain-index)
        #[arg(long)]
        baseline: PathBuf,

        /// Label for the baseline run
        #[arg(long, default_value = "plain-index")]
        baseline_label: String,

        /// JSON report for the candidate run (e.g. transformer)
        #[arg(long)]
        candidate: PathBuf,

        /// Label for the candidate run
        #[arg(long, default_value = "transformer")]
        candidate_label: String,

        /// Write JSON comparison report to this file (stdout if omitted)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Also write a Markdown report to this file
        #[arg(long)]
        markdown: Option<PathBuf>,
    },

    /// Run the plain-index DNA retrieval benchmark.
    ///
    /// Loads a JSONL corpus produced by `dna-corpus`, generates query
    /// substrings, runs exact substring search, and outputs a JSON report
    /// with hit@k, precision@k, recall@k, NDCG@k, and latency statistics.
    DnaBenchmark {
        /// Input JSONL corpus file (produced by dna-corpus)
        #[arg(short, long)]
        corpus: PathBuf,

        /// Number of top results for @k metrics
        #[arg(long, default_value = "10")]
        k: usize,

        /// PRNG seed for deterministic query generation
        #[arg(long, default_value = "42")]
        seed: u64,

        /// Number of queries to generate
        #[arg(long, default_value = "100")]
        n_queries: usize,

        /// Length in bases of each query substring
        #[arg(long, default_value = "30")]
        query_len: usize,

        /// Minimum sequence length to include in index
        #[arg(long, default_value = "30")]
        min_seq_len: usize,

        /// Write JSON report to this file instead of stdout
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Scan FASTA files for CRISPR targets (fast 2-bit k-mer sliding window).
    CrisprScan {
        /// Input FASTA file or directory
        #[arg(short, long)]
        input: PathBuf,

        /// Output CSV file path
        #[arg(short, long, default_value = "crispr_targets.csv")]
        output: PathBuf,

        /// K-mer size (default 23 = 20bp guide + 3bp PAM context)
        #[arg(short, long, default_value = "23")]
        kmer_size: usize,

        /// Maximum targets to output (sorted by occurrences desc)
        #[arg(long, default_value = "10000")]
        max_targets: usize,

        /// Pathogen name for CSV output
        #[arg(short, long)]
        name: String,
    },

    /// Scan FASTA for conservation of specific CRISPR guide sequences.
    ///
    /// Streams a large FASTA file, checks each genome for the presence of
    /// guide sequences (+ reverse complements) using Aho-Corasick, and
    /// reports per-guide conservation rates.
    GuideConservation {
        /// Input FASTA file (can be very large, streamed record-by-record)
        #[arg(short, long)]
        fasta: PathBuf,

        /// TSV file with guide_id and sequence columns
        #[arg(short, long)]
        guides: PathBuf,

        /// Output JSON file path
        #[arg(short, long, default_value = "guide_conservation.json")]
        output: PathBuf,
    },

    /// Batch search: load one index, search many patterns from a file.
    ///
    /// Reads one pattern per line from the input file, counts occurrences
    /// of each in the loaded index, and writes CSV to stdout: pattern,count
    BatchSearch {
        /// Path to the FM-index file
        #[arg(short, long)]
        index: PathBuf,

        /// File with one search pattern per line
        #[arg(short, long)]
        patterns: PathBuf,
    },

    /// Print the LOOM operational SLO profile (capacity planning + runbook) as JSON.
    OpsProfile {
        /// Write JSON profile to this file instead of stdout
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Scan a reference genome for CRISPR off-targets (pigeonhole seeding + 2-bit Hamming).
    ///
    /// For each guide in --guides, finds every position in --reference where the
    /// guide (or its reverse complement) matches with at most --max-mismatches edits.
    /// Uses the pigeonhole principle to seed-filter candidates before Hamming verify:
    /// splitting 20 nt into 4 × 5-nt seeds guarantees ≥ 1 exact seed for ≤ 3 mismatches.
    /// Parallel over chromosome sub-chunks via Rayon.
    OffTargetScan {
        /// CSV or plain-text file containing 20-nt guide spacer sequences.
        /// CSV must have a column named sequence_23mer, sequence_20mer, guide_sequence, or sequence.
        #[arg(short, long)]
        guides: PathBuf,

        /// Reference genome FASTA (plain or .gz, e.g. data/human/GRCh38.fa).
        #[arg(short, long)]
        reference: PathBuf,

        /// Output CSV file path.
        #[arg(short, long, default_value = "offtargets.csv")]
        output: PathBuf,

        /// Maximum mismatches to report (1–3, default 3).
        #[arg(long, default_value = "3")]
        max_mismatches: u8,

        /// Only use the first N guides from the input (0 = all).
        #[arg(long, default_value = "0")]
        max_guides: usize,

        /// Parallel sub-chunk size in guide-start positions (default 4 000 000).
        #[arg(long, default_value = "4000000")]
        chunk_size: usize,
    },

    /// Sequence-level novelty verification: compare novel CRISPR targets
    /// against a comprehensive published-guide database using 2-bit Hamming distance.
    ///
    /// For each of the novel targets, checks both strands of every published guide.
    /// A target is "sequence-novel" if no published guide matches within --max-mismatches.
    GuideNovelty {
        /// JSON file with novel target regions (novel_targets_52.json format).
        #[arg(short, long)]
        novel_targets: PathBuf,

        /// JSON file with published guide database (from compile_published_guides.py).
        #[arg(short, long)]
        published_guides: PathBuf,

        /// Output JSON report path.
        #[arg(short, long, default_value = "novelty_report.json")]
        output: PathBuf,

        /// Maximum Hamming mismatches to consider a "match" (default 3).
        #[arg(long, default_value = "3")]
        max_mismatches: u8,
    },

    /// Add a pre-serialized FM-index (FMIX) section to an .idx file so
    /// WASM can load it without rebuilding the suffix array.
    PrepareWasm {
        /// Path to the .idx file to upgrade
        #[arg(short, long)]
        index: PathBuf,

        /// Output path (defaults to overwriting input)
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Split an index into WASM-ready shards for Web Worker loading.
    PrepareWasmShards {
        /// Path to the .idx file to split
        #[arg(short, long)]
        index: PathBuf,

        /// Directory to write shard files into
        #[arg(short, long, default_value = "wasm_shards")]
        output_dir: PathBuf,

        /// Target bytes per shard (default: 200 MiB)
        #[arg(long, default_value = "209715200")]
        target_shard_bytes: usize,
    },

    /// PubMed literature gap scanner v4 — ontology-integrated, multi-strategy.
    ///
    /// Queries PubMed + Europe PMC across all configured pathogens, auto-generates
    /// gene synonyms from NCBI ontology annotations, and produces confidence-scored
    /// gap classifications (CONFIRMED / PROBABLE / UNCERTAIN / FALSE).
    PubmedScan {
        /// Scan only this pathogen key (e.g. sars-cov-2). Merges with existing output.
        #[arg(long)]
        only: Option<String>,

        /// When using --only, overwrite output instead of merging with existing results.
        #[arg(long)]
        no_merge: bool,

        /// Clear corpus cache before scanning (forces fresh API calls).
        #[arg(long)]
        no_cache: bool,

        /// Output JSON path (default: web/data/pubmed-scan-v4.json)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Path to ontology-enrichment.json (default: web/data/ontology-enrichment.json)
        #[arg(long)]
        ontology: Option<PathBuf>,
    },
}
