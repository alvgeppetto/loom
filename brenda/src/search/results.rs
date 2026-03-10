//! Search results representation

use std::collections::HashMap;
use std::path::PathBuf;

/// A single search result with location and context
#[derive(Debug, Clone)]
pub struct SearchResult {
    pub file: PathBuf,
    pub line: usize,
    pub column: usize,
    pub byte_offset: usize,
    pub context_before: String,
    pub matched_text: String,
    pub context_after: String,
}

/// Collection of search results
#[derive(Debug)]
pub struct SearchResults {
    pub results: Vec<SearchResult>,
    pub total_count: usize,
    pub truncated: bool,
}

impl SearchResults {
    pub fn empty() -> Self {
        Self {
            results: Vec::new(),
            total_count: 0,
            truncated: false,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.results.is_empty()
    }
}

/// Returns a scope-based weight: source files rank higher than test/generated files.
fn scope_weight(file: &str) -> f64 {
    if file.contains("generated") || file.contains("gen_") {
        0.5
    } else if file.contains("/test") || file.ends_with("_test.rs") || file.contains("tests/") {
        0.8
    } else {
        1.0
    }
}

/// A file-level scored result combining counts across all matching query terms.
#[derive(Debug, Clone)]
pub struct ScoredFile {
    pub file:        String,
    pub line:        usize,
    pub score:       f64,
    pub term_count:  usize,   // distinct query terms matching this file
    pub hit_count:   usize,   // total occurrence count across all terms
}

impl ScoredFile {
    pub fn aggregate_score(&self, alpha: f64) -> f64 {
        self.score + alpha * self.term_count as f64
    }
}

/// Rank a list of raw FM-index hits into scored, deduplicated file results.
///
/// `hits` is a slice of `(file, line, count)` tuples where `count` is the
/// number of occurrences of a single query term in that file.
pub fn rank_results(
    hits: &[(String, usize, usize)],  // (file, line, count)
    alpha: f64,
) -> Vec<ScoredFile> {
    let mut file_map: HashMap<String, ScoredFile> = HashMap::new();
    for (file, line, count) in hits {
        let freq_score = (1.0 + *count as f64).ln();
        let scope_w    = scope_weight(file);
        let entry = file_map.entry(file.clone()).or_insert(
            ScoredFile { file: file.clone(), line: *line,
                         score: 0.0, term_count: 0, hit_count: 0 }
        );
        entry.score     += freq_score * scope_w;
        entry.term_count += 1;
        entry.hit_count  += count;
    }
    let mut ranked: Vec<ScoredFile> = file_map.into_values().collect();
    ranked.sort_by(|a, b|
        b.aggregate_score(alpha)
         .partial_cmp(&a.aggregate_score(alpha))
         .unwrap()
    );
    ranked
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rank_results_orders_by_aggregate_score() {
        let hits = vec![
            ("src/channel.rs".to_string(), 10, 5),
            ("src/queue.rs".to_string(),   3,  1),
            ("src/channel.rs".to_string(), 10, 3),  // second term hit on same file
        ];
        let ranked = rank_results(&hits, 0.5);
        assert!(!ranked.is_empty());
        // channel.rs has two term hits and higher total count, should be first
        assert_eq!(ranked[0].file, "src/channel.rs");
        assert_eq!(ranked[0].term_count, 2);
        assert_eq!(ranked[0].hit_count, 8);
    }

    #[test]
    fn test_rank_results_scope_weight_depresses_test_files() {
        let hits = vec![
            ("src/core.rs".to_string(),    1, 10),
            ("tests/core_test.rs".to_string(), 1, 10),
        ];
        let ranked = rank_results(&hits, 0.5);
        assert_eq!(ranked.len(), 2);
        // source file should rank first due to scope_weight = 1.0 vs 0.8
        assert_eq!(ranked[0].file, "src/core.rs");
    }
}
