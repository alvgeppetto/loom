//! Federated shard search — parallel query fan-out + deterministic merge.
//!
//! Loads a [`ShardManifest`] and fans a query out to all shards in parallel
//! using a bounded worker pool built on `std::thread` + channels.  Results
//! are merged with a deterministic stable sort:
//!   primary key:   `shard_id` ascending   (preserves corpus / manifest order)
//!   secondary key: `position` ascending   (stable within a shard)

use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::thread;

use serde::{Deserialize, Serialize};

use crate::indexer::shards::ShardManifest;
use crate::indexer::{IndexBuilder, IndexError};

// ── Public types ─────────────────────────────────────────────────────────────

/// One match result from a federated search across shards.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FederatedMatch {
    /// Zero-based shard index from which this result was retrieved.
    pub shard_id: usize,
    /// Byte offset within the shard corpus.
    pub position: u64,
    /// Source file path, if resolvable from the shard (line, column).
    pub file: Option<String>,
    /// Line number within the file (1-based), if resolvable.
    pub line: Option<usize>,
    /// Column number within the line (1-based), if resolvable.
    pub col: Option<usize>,
    /// Context snippet around the match.
    pub context: String,
}

/// Summary stats produced alongside a federated search.
#[derive(Debug, Serialize, Deserialize)]
pub struct FederatedStats {
    /// Total number of shards searched.
    pub shards_searched: usize,
    /// Total match count across all shards.
    pub total_matches: usize,
    /// Wall-clock time for the fan-out + merge (milliseconds).
    pub elapsed_ms: f64,
}

// ── FederatedSearcher ─────────────────────────────────────────────────────────

/// Searches across multiple index shards with parallel fan-out.
pub struct FederatedSearcher {
    /// Absolute paths to each shard index file, in manifest order.
    shard_paths: Vec<PathBuf>,
    /// Maximum number of worker threads to spawn concurrently.
    max_workers: usize,
}

impl FederatedSearcher {
    /// Construct from a [`ShardManifest`] JSON file.
    ///
    /// `max_workers` caps concurrent worker threads.  Pass `0` to spawn one
    /// thread per shard (bounded only by shard count).
    pub fn from_manifest(manifest_path: &Path, max_workers: usize) -> Result<Self, IndexError> {
        let manifest = ShardManifest::load(manifest_path)?;
        let base_dir = manifest_path
            .parent()
            .unwrap_or(Path::new("."))
            .to_path_buf();

        let shard_paths: Vec<PathBuf> = manifest
            .shards
            .iter()
            .map(|e| base_dir.join(&e.path))
            .collect();

        let workers = if max_workers == 0 || max_workers > shard_paths.len() {
            shard_paths.len().max(1)
        } else {
            max_workers
        };

        Ok(Self {
            shard_paths,
            max_workers: workers,
        })
    }

    /// Search all shards for `pattern` and return a deterministically-ordered
    /// merged result list.
    ///
    /// Results are sorted by `(shard_id asc, position asc)` so output is
    /// reproducible regardless of thread scheduling.  Returns both the matches
    /// and a [`FederatedStats`] summary.
    pub fn search(
        &self,
        pattern: &str,
        context_chars: usize,
    ) -> (Vec<FederatedMatch>, FederatedStats) {
        let t0 = std::time::Instant::now();
        let pattern = pattern.to_string();
        let n_shards = self.shard_paths.len();

        let results: Vec<Vec<FederatedMatch>> = self.fan_out(move |idx, path| {
            let builder = match IndexBuilder::load(path) {
                Ok(b) => b,
                Err(_) => return vec![],
            };
            let positions = builder.search(&pattern);
            positions
                .into_iter()
                .map(|pos| {
                    let (file, line, col) = match builder.line_col_at_offset(pos) {
                        Some((f, l, c)) => (Some(f), Some(l), Some(c)),
                        None => (None, None, None),
                    };
                    let context =
                        builder.context_at_offset(pos, context_chars, context_chars);
                    FederatedMatch {
                        shard_id: idx,
                        position: pos,
                        file,
                        line,
                        col,
                        context,
                    }
                })
                .collect()
        });

        let mut all: Vec<FederatedMatch> = results.into_iter().flatten().collect();
        // Deterministic stable sort: shard_id ascending, then position ascending.
        all.sort_by_key(|m| (m.shard_id, m.position));

        let total = all.len();
        let stats = FederatedStats {
            shards_searched: n_shards,
            total_matches: total,
            elapsed_ms: t0.elapsed().as_secs_f64() * 1000.0,
        };
        (all, stats)
    }

    /// Count total occurrences of `pattern` across all shards.
    pub fn count(&self, pattern: &str) -> (u64, FederatedStats) {
        let t0 = std::time::Instant::now();
        let n_shards = self.shard_paths.len();
        let pattern = pattern.to_string();

        let counts: Vec<u64> = self.fan_out(move |_idx, path| {
            IndexBuilder::load(path)
                .map(|b| b.count(&pattern))
                .unwrap_or(0)
        });

        let total: u64 = counts.into_iter().sum();
        let stats = FederatedStats {
            shards_searched: n_shards,
            total_matches: total as usize,
            elapsed_ms: t0.elapsed().as_secs_f64() * 1000.0,
        };
        (total, stats)
    }

    // ── Internal fan-out ──────────────────────────────────────────────────────

    /// Run `f(shard_id, path)` on each shard concurrently with a bounded
    /// worker pool and collect results preserving manifest (shard) order.
    fn fan_out<T, F>(&self, f: F) -> Vec<T>
    where
        T: Send + 'static,
        F: Fn(usize, &Path) -> T + Send + Sync + 'static,
    {
        let n = self.shard_paths.len();
        if n == 0 {
            return vec![];
        }

        let f = Arc::new(f);
        // Pre-allocate result slots indexed by shard_id.
        let results: Arc<Mutex<Vec<Option<T>>>> =
            Arc::new(Mutex::new((0..n).map(|_| None).collect()));

        // Work queue: (shard_id, path).
        let work: Arc<Mutex<std::collections::VecDeque<(usize, PathBuf)>>> = Arc::new(
            Mutex::new(
                self.shard_paths
                    .iter()
                    .enumerate()
                    .map(|(i, p)| (i, p.clone()))
                    .collect(),
            ),
        );

        let num_threads = self.max_workers.min(n);
        let mut handles = Vec::with_capacity(num_threads);

        for _ in 0..num_threads {
            let work = Arc::clone(&work);
            let results = Arc::clone(&results);
            let f = Arc::clone(&f);
            let handle = thread::spawn(move || loop {
                let item = work.lock().unwrap().pop_front();
                let (idx, path) = match item {
                    Some(t) => t,
                    None => break,
                };
                let val = f(idx, &path);
                results.lock().unwrap()[idx] = Some(val);
            });
            handles.push(handle);
        }

        for h in handles {
            let _ = h.join();
        }

        Arc::try_unwrap(results)
            .ok()
            .expect("BUG: Arc still held after all threads joined")
            .into_inner()
            .expect("Mutex poisoned")
            .into_iter()
            .flatten()
            .collect()
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::indexer::{Corpus, SourceFile, FILE_DELIMITER};
    use crate::indexer::shards::build_shards;
    use std::path::PathBuf;
    use tempfile::TempDir;

    fn tiny_corpus(pairs: &[(&str, &str)]) -> Corpus {
        let mut concatenated: Vec<u8> = Vec::new();
        let mut files = Vec::new();
        for (name, content) in pairs {
            let start = concatenated.len();
            concatenated.push(FILE_DELIMITER);
            concatenated.extend_from_slice(name.as_bytes());
            concatenated.push(FILE_DELIMITER);
            concatenated.extend_from_slice(content.as_bytes());
            let end = concatenated.len();
            files.push(SourceFile {
                path: PathBuf::from(name),
                content: content.to_string(),
                start_offset: start,
                end_offset: end,
            });
        }
        Corpus { files, concatenated }
    }

    fn build_test_shards(dir: &TempDir) -> std::path::PathBuf {
        let corpus = tiny_corpus(&[
            ("alpha.rs", "fn alpha() { let x = 42; }"),
            ("beta.rs", "fn beta() { let y = 99; }"),
            ("gamma.rs", "fn gamma() { let x = 0; }"),
        ]);
        // Force 3 shards (one file each).
        build_shards(&corpus, dir.path(), 1).unwrap();
        dir.path().join("shard_manifest.json")
    }

    #[test]
    fn federated_count_sums_all_shards() {
        let dir = TempDir::new().unwrap();
        let manifest_path = build_test_shards(&dir);

        let searcher = FederatedSearcher::from_manifest(&manifest_path, 0).unwrap();
        let (count, stats) = searcher.count("fn");
        // Each of the 3 shards contains "fn" once → total 3.
        assert_eq!(count, 3, "expected 3 'fn' occurrences across all shards");
        assert_eq!(stats.shards_searched, 3);
    }

    #[test]
    fn federated_search_returns_sorted_results() {
        let dir = TempDir::new().unwrap();
        let manifest_path = build_test_shards(&dir);

        let searcher = FederatedSearcher::from_manifest(&manifest_path, 2).unwrap();
        let (matches, stats) = searcher.search("fn", 20);

        assert_eq!(matches.len(), 3);
        assert_eq!(stats.total_matches, 3);

        // Verify sorted by (shard_id, position).
        for w in matches.windows(2) {
            let a = &w[0];
            let b = &w[1];
            assert!(
                (a.shard_id, a.position) <= (b.shard_id, b.position),
                "results must be sorted by (shard_id, position)"
            );
        }
    }

    #[test]
    fn federated_search_result_has_context() {
        let dir = TempDir::new().unwrap();
        let manifest_path = build_test_shards(&dir);

        let searcher = FederatedSearcher::from_manifest(&manifest_path, 0).unwrap();
        let (matches, _) = searcher.search("alpha", 40);

        assert!(!matches.is_empty(), "expected at least one match for 'alpha'");
        assert!(
            !matches[0].context.is_empty(),
            "context should not be empty"
        );
    }

    #[test]
    fn federated_search_cross_shard_term() {
        let dir = TempDir::new().unwrap();
        let manifest_path = build_test_shards(&dir);

        // "let x" appears in alpha.rs and gamma.rs (shards 0 and 2).
        let searcher = FederatedSearcher::from_manifest(&manifest_path, 0).unwrap();
        let (matches, _) = searcher.search("let x", 20);

        assert_eq!(matches.len(), 2, "'let x' should appear in 2 shards");
        assert_eq!(matches[0].shard_id, 0);
        assert_eq!(matches[1].shard_id, 2);
    }

    #[test]
    fn federated_stats_elapsed_is_non_negative() {
        let dir = TempDir::new().unwrap();
        let manifest_path = build_test_shards(&dir);

        let searcher = FederatedSearcher::from_manifest(&manifest_path, 1).unwrap();
        let (_, stats) = searcher.search("fn", 10);
        assert!(stats.elapsed_ms >= 0.0);
    }

    #[test]
    fn from_manifest_respects_worker_cap() {
        let dir = TempDir::new().unwrap();
        let manifest_path = build_test_shards(&dir);

        // max_workers = 1 → should still find all results.
        let searcher = FederatedSearcher::from_manifest(&manifest_path, 1).unwrap();
        let (count, _) = searcher.count("fn");
        assert_eq!(count, 3);
    }
}
