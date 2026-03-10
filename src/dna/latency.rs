//! DNA query latency tracking and statistics.
//!
//! Provides infrastructure for measuring search/parse latency with cold-start
//! and warm (steady-state) breakdowns, plus standard percentile stats.
//!
//! # Usage
//!
//! ```no_run
//! use loom::dna::latency::{time_call, compute_latency_stats};
//!
//! let mut timings: Vec<f64> = Vec::new();
//! for _ in 0..20 {
//!     let (_result, elapsed_ms) = time_call(|| expensive_operation());
//!     timings.push(elapsed_ms);
//! }
//! let stats = compute_latency_stats(&timings).unwrap();
//! println!("p99: {:.2} ms", stats.p99_ms);
//! # fn expensive_operation() {}
//! ```

use std::time::Instant;

use serde::{Deserialize, Serialize};

/// Latency statistics derived from a series of millisecond timings.
///
/// The first timing is considered the **cold** measurement (index/data not yet
/// in CPU cache or OS page cache).  All subsequent timings feed the **warm**
/// distribution used for median, p95, and p99.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LatencyStats {
    /// Total number of timed samples.
    pub n: usize,
    /// Cold (first call) latency in milliseconds.
    pub cold_ms: f64,
    /// Median warm latency in milliseconds.
    pub warm_median_ms: f64,
    /// 95th-percentile warm latency in milliseconds.
    pub p95_ms: f64,
    /// 99th-percentile warm latency in milliseconds.
    pub p99_ms: f64,
    /// Mean across **all** samples (cold + warm) in milliseconds.
    pub mean_ms: f64,
    /// Minimum across all samples.
    pub min_ms: f64,
    /// Maximum across all samples.
    pub max_ms: f64,
}

/// Compute latency statistics from a slice of millisecond timings.
///
/// Returns `None` when the slice is empty.  The first element is treated as
/// the cold-start sample; the remaining elements form the warm distribution.
/// If only one sample exists it is used for all stats.
pub fn compute_latency_stats(timings_ms: &[f64]) -> Option<LatencyStats> {
    if timings_ms.is_empty() {
        return None;
    }

    let n = timings_ms.len();
    let cold_ms = timings_ms[0];
    let mean_ms = timings_ms.iter().sum::<f64>() / n as f64;
    let min_ms = timings_ms
        .iter()
        .cloned()
        .fold(f64::INFINITY, f64::min);
    let max_ms = timings_ms
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);

    // Warm distribution excludes the first (cold) sample when n > 1.
    let warm: Vec<f64> = if n > 1 {
        timings_ms[1..].to_vec()
    } else {
        timings_ms.to_vec()
    };

    let mut warm_sorted = warm.clone();
    warm_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let warm_median_ms = percentile_sorted(&warm_sorted, 50.0);
    let p95_ms = percentile_sorted(&warm_sorted, 95.0);
    let p99_ms = percentile_sorted(&warm_sorted, 99.0);

    Some(LatencyStats {
        n,
        cold_ms,
        warm_median_ms,
        p95_ms,
        p99_ms,
        mean_ms,
        min_ms,
        max_ms,
    })
}

/// Compute a percentile from an already-sorted slice using nearest-rank.
///
/// `p` must be in [0, 100].  Returns 0.0 for an empty slice.
pub fn percentile_sorted(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let idx = ((p / 100.0) * (sorted.len() - 1) as f64).round() as usize;
    sorted[idx.min(sorted.len() - 1)]
}

/// Time a single function call.
///
/// Returns the call's result and the elapsed wall-clock time in milliseconds.
pub fn time_call<T>(f: impl FnOnce() -> T) -> (T, f64) {
    let start = Instant::now();
    let result = f();
    let elapsed_ms = start.elapsed().as_secs_f64() * 1_000.0;
    (result, elapsed_ms)
}

#[cfg(test)]
mod tests {
    use super::{compute_latency_stats, percentile_sorted, time_call};

    #[test]
    fn returns_none_for_empty_input() {
        assert!(compute_latency_stats(&[]).is_none());
    }

    #[test]
    fn single_sample_populates_all_fields() {
        let stats = compute_latency_stats(&[42.0]).unwrap();
        assert_eq!(stats.n, 1);
        assert!((stats.cold_ms - 42.0).abs() < 1e-9);
        assert!((stats.warm_median_ms - 42.0).abs() < 1e-9);
        assert!((stats.p95_ms - 42.0).abs() < 1e-9);
        assert!((stats.p99_ms - 42.0).abs() < 1e-9);
        assert!((stats.mean_ms - 42.0).abs() < 1e-9);
        assert!((stats.min_ms - 42.0).abs() < 1e-9);
        assert!((stats.max_ms - 42.0).abs() < 1e-9);
    }

    #[test]
    fn cold_is_first_sample_warm_excludes_it() {
        // Timings: cold=100, warm=[1,2,3,4,5]
        let timings = vec![100.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        let stats = compute_latency_stats(&timings).unwrap();
        assert_eq!(stats.n, 6);
        assert!((stats.cold_ms - 100.0).abs() < 1e-9);
        // Warm median of [1,2,3,4,5] = 3
        assert!((stats.warm_median_ms - 3.0).abs() < 1e-9);
        // Mean across all = (100+1+2+3+4+5)/6
        let expected_mean = 115.0 / 6.0;
        assert!((stats.mean_ms - expected_mean).abs() < 1e-6);
        assert!((stats.min_ms - 1.0).abs() < 1e-9);
        assert!((stats.max_ms - 100.0).abs() < 1e-9);
    }

    #[test]
    fn p95_and_p99_on_100_sample_range() {
        let timings: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let stats = compute_latency_stats(&timings).unwrap();
        // Warm = [1..99], sorted. p95 = index round(0.95*98)=93 → value 94
        assert!(stats.p95_ms >= 90.0);
        assert!(stats.p99_ms >= 95.0);
        assert!(stats.p99_ms <= 99.0);
    }

    #[test]
    fn percentile_sorted_edge_cases() {
        assert!((percentile_sorted(&[5.0], 0.0) - 5.0).abs() < 1e-9);
        assert!((percentile_sorted(&[5.0], 100.0) - 5.0).abs() < 1e-9);
        assert!((percentile_sorted(&[], 50.0) - 0.0).abs() < 1e-9);

        let s = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert!((percentile_sorted(&s, 0.0) - 1.0).abs() < 1e-9);
        assert!((percentile_sorted(&s, 50.0) - 3.0).abs() < 1e-9);
        assert!((percentile_sorted(&s, 100.0) - 5.0).abs() < 1e-9);
    }

    #[test]
    fn time_call_returns_result_and_positive_elapsed() {
        let (value, elapsed_ms) = time_call(|| 2 + 2);
        assert_eq!(value, 4);
        assert!(elapsed_ms >= 0.0);
    }
}
