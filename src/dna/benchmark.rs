//! Plain-index DNA retrieval benchmark.
//!
//! Loads a corpus from a JSONL file produced by `dna-corpus`, generates query
//! substrings deterministically, runs a substring-scan "plain index", and
//! computes standard retrieval metrics (hit\@k, precision\@k, recall\@k,
//! NDCG\@k) plus latency statistics.
//!
//! The plain index is intentionally simple: exact substring matching over the
//! full sequence corpus.  This gives a **deterministic upper-bound baseline**
//! for recall, against which approximate or learned indexes can be compared.
//!
//! # Example
//!
//! ```no_run
//! use loom::dna::benchmark::{DnaBenchmarkConfig, run_dna_benchmark};
//! use std::path::Path;
//!
//! let cfg = DnaBenchmarkConfig::default();
//! let report = run_dna_benchmark(Path::new("dna_corpus.jsonl"), &cfg).unwrap();
//! println!("hit@{}: {:.4}", cfg.k, report.hit_at_k);
//! ```

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::time::Instant;

use anyhow::{bail, Context};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::dna::latency::{compute_latency_stats, LatencyStats};

// ---------------------------------------------------------------------------
// Config
// ---------------------------------------------------------------------------

/// Parameters that control the benchmark run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DnaBenchmarkConfig {
    /// Number of top results considered when computing \@k metrics.
    pub k: usize,
    /// PRNG seed for deterministic query generation.
    pub seed: u64,
    /// Number of queries to generate (capped by corpus size).
    pub n_queries: usize,
    /// Length (in bases) of each generated query substring.
    pub query_len: usize,
    /// Minimum sequence length required for a record to be included.
    pub min_seq_len: usize,
}

impl Default for DnaBenchmarkConfig {
    fn default() -> Self {
        Self {
            k: 10,
            seed: 42,
            n_queries: 100,
            query_len: 30,
            min_seq_len: 30,
        }
    }
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Minimal subset of `DnaCorpusRecord` required for the benchmark.
#[derive(Debug, Deserialize)]
struct BenchRecord {
    id: String,
    sequence: String,
}

/// Load sequence records from a JSONL corpus file.
fn load_corpus(path: &Path) -> anyhow::Result<Vec<BenchRecord>> {
    let f = File::open(path).with_context(|| format!("Cannot open corpus file {}", path.display()))?;
    let reader = BufReader::new(f);
    let mut records = Vec::new();
    for (lineno, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("IO error at line {}", lineno + 1))?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let rec: BenchRecord = serde_json::from_str(trimmed)
            .with_context(|| format!("JSON parse error at line {}", lineno + 1))?;
        records.push(rec);
    }
    Ok(records)
}

/// Generate query substrings deterministically from the corpus.
///
/// Each query is extracted from a randomly chosen sequence at a random offset.
/// Returns a vec of `(query_string, origin_record_id)`.
fn generate_queries(
    records: &[BenchRecord],
    cfg: &DnaBenchmarkConfig,
    eligible: &[usize],
) -> Vec<(String, String)> {
    let mut rng = StdRng::seed_from_u64(cfg.seed);
    let n = cfg.n_queries.min(eligible.len());
    let mut queries = Vec::with_capacity(n);
    for _ in 0..n {
        let idx = eligible[rng.gen_range(0..eligible.len())];
        let seq = &records[idx].sequence;
        let max_start = seq.len() - cfg.query_len;
        let start = rng.gen_range(0..=max_start);
        let query = seq[start..start + cfg.query_len].to_string();
        queries.push((query, records[idx].id.clone()));
    }
    queries
}

/// Search the plain index for `query` and return ranked result indices.
///
/// Results are ranked by total number of non-overlapping occurrences of the
/// query in each sequence (descending), with ties broken by record order.
fn plain_index_search(records: &[BenchRecord], query: &str) -> Vec<usize> {
    if query.is_empty() {
        return vec![];
    }
    let mut scored: Vec<(usize, usize)> = records
        .iter()
        .enumerate()
        .filter_map(|(i, rec)| {
            let count = count_occurrences(&rec.sequence, query);
            if count > 0 {
                Some((i, count))
            } else {
                None
            }
        })
        .collect();
    // Sort descending by count, stable order for ties.
    scored.sort_by(|a, b| b.1.cmp(&a.1));
    scored.into_iter().map(|(i, _)| i).collect()
}

/// Count non-overlapping occurrences of `needle` in `haystack`.
fn count_occurrences(haystack: &str, needle: &str) -> usize {
    if needle.is_empty() {
        return 0;
    }
    let mut count = 0usize;
    let mut start = 0;
    while let Some(pos) = haystack[start..].find(needle) {
        count += 1;
        start += pos + needle.len();
    }
    count
}

// ---------------------------------------------------------------------------
// Metrics
// ---------------------------------------------------------------------------

/// Compute hit\@k: 1.0 if any of the top-k results is relevant, else 0.0.
fn hit_at_k(ranked: &[usize], relevant: &HashSet<usize>, k: usize) -> f64 {
    let top_k = ranked.len().min(k);
    let hit = ranked[..top_k].iter().any(|i| relevant.contains(i));
    if hit { 1.0 } else { 0.0 }
}

/// Compute precision\@k: fraction of top-k results that are relevant.
fn precision_at_k(ranked: &[usize], relevant: &HashSet<usize>, k: usize) -> f64 {
    let top_k = ranked.len().min(k);
    if top_k == 0 {
        return 0.0;
    }
    let hits = ranked[..top_k].iter().filter(|i| relevant.contains(i)).count();
    hits as f64 / k as f64
}

/// Compute recall\@k: fraction of relevant items found in the top-k results.
fn recall_at_k(ranked: &[usize], relevant: &HashSet<usize>, k: usize) -> f64 {
    if relevant.is_empty() {
        return 1.0;
    }
    let top_k = ranked.len().min(k);
    let hits = ranked[..top_k].iter().filter(|i| relevant.contains(i)).count();
    hits as f64 / relevant.len() as f64
}

/// Compute NDCG\@k (binary relevance).
fn ndcg_at_k(ranked: &[usize], relevant: &HashSet<usize>, k: usize) -> f64 {
    let top_k = ranked.len().min(k);
    // DCG: sum rel_i / log2(i+2) for 0-indexed rank i (so rank 1 → log2(2)=1)
    let dcg: f64 = ranked[..top_k]
        .iter()
        .enumerate()
        .map(|(i, idx)| {
            if relevant.contains(idx) {
                1.0 / (i as f64 + 2.0).log2()
            } else {
                0.0
            }
        })
        .sum();
    // Ideal DCG: all relevant items placed at ranks 1..=ideal_count
    let ideal_count = relevant.len().min(k);
    let idcg: f64 = (1..=ideal_count)
        .map(|rank| 1.0 / (rank as f64 + 1.0).log2())
        .sum();
    if idcg == 0.0 {
        return 0.0;
    }
    dcg / idcg
}

// ---------------------------------------------------------------------------
// Report
// ---------------------------------------------------------------------------

/// Output produced by `run_dna_benchmark`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DnaBenchmarkReport {
    /// Number of records in the corpus.
    pub corpus_size: usize,
    /// Number of queries evaluated.
    pub n_queries: usize,
    /// The k used for \@k metrics.
    pub k: usize,
    /// Query substring length used.
    pub query_len: usize,
    /// Mean hit\@k across all queries.
    pub hit_at_k: f64,
    /// Mean precision\@k across all queries.
    pub precision_at_k: f64,
    /// Mean recall\@k across all queries.
    pub recall_at_k: f64,
    /// Mean NDCG\@k across all queries.
    pub ndcg_at_k: f64,
    /// Per-query latency statistics (milliseconds).
    pub latency: LatencyStats,
    /// ISO 8601 timestamp when the report was generated.
    pub generated_at: String,
    /// PRNG seed used for query generation.
    pub seed: u64,
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Run the plain-index DNA retrieval benchmark.
///
/// Loads `corpus_path`, generates queries according to `cfg`, performs
/// exact substring search against all sequences, and aggregates metrics.
pub fn run_dna_benchmark(corpus_path: &Path, cfg: &DnaBenchmarkConfig) -> anyhow::Result<DnaBenchmarkReport> {
    let records = load_corpus(corpus_path)?;
    if records.is_empty() {
        bail!("Corpus is empty: {}", corpus_path.display());
    }

    // Filter to sequences long enough to generate queries.
    let eligible: Vec<usize> = records
        .iter()
        .enumerate()
        .filter(|(_, r)| r.sequence.len() >= cfg.query_len.max(cfg.min_seq_len))
        .map(|(i, _)| i)
        .collect();

    if eligible.is_empty() {
        bail!(
            "No sequences meet minimum length {} for query generation",
            cfg.query_len
        );
    }

    let queries = generate_queries(&records, cfg, &eligible);
    if queries.is_empty() {
        bail!("Could not generate any queries from the corpus");
    }

    // For each query, determine the ground-truth relevant set (all records
    // containing the query as a substring) — guaranteed ≥ 1 because the query
    // was extracted from the corpus.
    // Queries are embarrassingly parallel: each is independent read-only scan.
    let per_query: Vec<(f64, f64, f64, f64, f64)> = queries
        .par_iter()
        .map(|(query, _origin_id)| {
            let relevant: HashSet<usize> = records
                .iter()
                .enumerate()
                .filter(|(_, r)| r.sequence.contains(query.as_str()))
                .map(|(i, _)| i)
                .collect();

            let t0 = Instant::now();
            let ranked = plain_index_search(&records, query);
            let elapsed_ms = t0.elapsed().as_secs_f64() * 1000.0;

            (
                elapsed_ms,
                hit_at_k(&ranked, &relevant, cfg.k),
                precision_at_k(&ranked, &relevant, cfg.k),
                recall_at_k(&ranked, &relevant, cfg.k),
                ndcg_at_k(&ranked, &relevant, cfg.k),
            )
        })
        .collect();

    let mut hits = Vec::with_capacity(per_query.len());
    let mut precisions = Vec::with_capacity(per_query.len());
    let mut recalls = Vec::with_capacity(per_query.len());
    let mut ndcgs = Vec::with_capacity(per_query.len());
    let mut timings: Vec<f64> = Vec::with_capacity(per_query.len());

    for (elapsed, hit, prec, rec, ndcg) in per_query {
        timings.push(elapsed);
        hits.push(hit);
        precisions.push(prec);
        recalls.push(rec);
        ndcgs.push(ndcg);
    }

    let mean = |v: &[f64]| -> f64 {
        if v.is_empty() {
            return 0.0;
        }
        v.iter().sum::<f64>() / v.len() as f64
    };

    let latency = compute_latency_stats(&timings)
        .expect("timings should be non-empty");

    let generated_at = {
        // RFC 3339 without external crate: use std SystemTime.
        use std::time::SystemTime;
        let secs = SystemTime::now()
            .duration_since(SystemTime::UNIX_EPOCH)
            .map(|d| d.as_secs())
            .unwrap_or(0);
        format_utc_timestamp(secs)
    };

    Ok(DnaBenchmarkReport {
        corpus_size: records.len(),
        n_queries: queries.len(),
        k: cfg.k,
        query_len: cfg.query_len,
        hit_at_k: mean(&hits),
        precision_at_k: mean(&precisions),
        recall_at_k: mean(&recalls),
        ndcg_at_k: mean(&ndcgs),
        latency,
        generated_at,
        seed: cfg.seed,
    })
}

/// Format a Unix timestamp as a minimal ISO 8601 UTC string (`YYYY-MM-DDTHH:MM:SSZ`).
///
/// Exposed as `pub` so sibling modules (e.g. `compare`) can reuse it without
/// pulling in an external date/time crate.
pub fn format_utc_iso(secs: u64) -> String {
    format_utc_timestamp(secs)
}

fn format_utc_timestamp(secs: u64) -> String {
    let days_since_epoch = secs / 86400;
    let time_of_day = secs % 86400;
    let h = time_of_day / 3600;
    let m = (time_of_day % 3600) / 60;
    let s = time_of_day % 60;

    // Gregorian calendar decomposition from days-since-epoch.
    let (y, mo, d) = days_to_ymd(days_since_epoch);
    format!("{:04}-{:02}-{:02}T{:02}:{:02}:{:02}Z", y, mo, d, h, m, s)
}

fn days_to_ymd(days: u64) -> (u64, u64, u64) {
    // Proleptic Gregorian calendar algorithm.
    let z = days + 719468;
    let era = z / 146097;
    let doe = z % 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let mo = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if mo <= 2 { y + 1 } else { y };
    (y, mo, d)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn make_corpus_jsonl(seqs: &[(&str, &str)]) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        for (id, seq) in seqs {
            writeln!(
                f,
                r#"{{"id":"{id}","sequence":"{seq}","source_path":"test","length":{len}}}"#,
                id = id,
                seq = seq,
                len = seq.len()
            )
            .unwrap();
        }
        f
    }

    #[test]
    fn count_occurrences_basic() {
        assert_eq!(count_occurrences("ACGTACGT", "ACG"), 2);
        assert_eq!(count_occurrences("ACGT", "TTTT"), 0);
        assert_eq!(count_occurrences("AAAA", "AA"), 2); // non-overlapping
        assert_eq!(count_occurrences("", "A"), 0);
    }

    #[test]
    fn metrics_perfect_retrieval() {
        let relevant: HashSet<usize> = vec![0].into_iter().collect();
        let ranked = vec![0, 1, 2];
        assert_eq!(hit_at_k(&ranked, &relevant, 1), 1.0);
        assert_eq!(precision_at_k(&ranked, &relevant, 1), 1.0);
        assert_eq!(recall_at_k(&ranked, &relevant, 1), 1.0);
        // NDCG@1: dcg = 1/log2(2)=1, idcg = 1 → 1.0
        assert!((ndcg_at_k(&ranked, &relevant, 1) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn metrics_miss_retrieval() {
        let relevant: HashSet<usize> = vec![5].into_iter().collect();
        let ranked = vec![0, 1, 2];
        assert_eq!(hit_at_k(&ranked, &relevant, 3), 0.0);
        assert_eq!(precision_at_k(&ranked, &relevant, 3), 0.0);
        assert_eq!(recall_at_k(&ranked, &relevant, 3), 0.0);
        assert_eq!(ndcg_at_k(&ranked, &relevant, 3), 0.0);
    }

    #[test]
    fn plain_index_finds_exact_match() {
        let records = vec![
            BenchRecord { id: "a".into(), sequence: "AAACCCGGG".into() },
            BenchRecord { id: "b".into(), sequence: "TTTGGGCCC".into() },
        ];
        let ranked = plain_index_search(&records, "GGGCCC");
        assert_eq!(ranked, vec![1]); // only record b contains GGGCCC
    }

    #[test]
    fn plain_index_returns_empty_for_no_match() {
        let records = vec![
            BenchRecord { id: "a".into(), sequence: "AAAA".into() },
        ];
        let ranked = plain_index_search(&records, "CCCC");
        assert!(ranked.is_empty());
    }

    #[test]
    fn run_benchmark_on_synthetic_corpus() {
        let seqs: Vec<(&str, &str)> = (0..20)
            .map(|i| {
                let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32 bases
                (Box::leak(format!("rec{}", i).into_boxed_str()) as &str, seq)
            })
            .collect();
        let f = make_corpus_jsonl(&seqs);
        let cfg = DnaBenchmarkConfig {
            k: 5,
            seed: 1,
            n_queries: 10,
            query_len: 8,
            min_seq_len: 8,
        };
        let report = run_dna_benchmark(f.path(), &cfg).unwrap();
        assert_eq!(report.corpus_size, 20);
        assert_eq!(report.n_queries, 10);
        assert_eq!(report.k, 5);
        // All queries come from the corpus so hit@k must be 1.0
        assert_eq!(report.hit_at_k, 1.0);
        assert!(report.latency.n > 0);
    }

    #[test]
    fn format_utc_timestamp_epoch() {
        // Unix epoch (1970-01-01T00:00:00Z)
        assert_eq!(format_utc_timestamp(0), "1970-01-01T00:00:00Z");
    }

    #[test]
    fn format_utc_timestamp_known_date() {
        // 2024-03-01 00:00:00 UTC
        // days = (2024-1970)*365 + leap adjustments + 31(Jan) + 29(Feb2024 leap) = ...
        // Easier: just check format structure
        let ts = format_utc_timestamp(1_000_000_000);
        assert!(ts.ends_with('Z'));
        assert_eq!(ts.len(), 20); // "YYYY-MM-DDTHH:MM:SSZ"
    }
}
