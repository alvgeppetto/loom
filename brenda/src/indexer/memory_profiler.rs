//! Preflight RAM estimator for FM-index construction.
//!
//! Estimates peak memory required to build an FM-index in-memory from a corpus
//! of `corpus_bytes` bytes, and recommends an appropriate build mode.

use serde::Serialize;

/// Multiplier for projected peak RAM relative to corpus size.
///
/// Empirically, building a sampled FM-index (level-2 suffix-order sampler)
/// requires approximately:
///   - 1× corpus text
///   - 1× BWT construction scratch
///   - 0.5× sampled suffix array (level 2 ≈ half density)
///   - 0.5× rank/occ table
///   - 2× SA construction peak (induced-sort scratch)
/// ≈ 5× corpus bytes, rounded up to 6 for safety.
const PEAK_RAM_MULTIPLIER: u64 = 6;

/// Fraction of available RAM (0–100) that is considered safe for in-memory
/// index construction.  Staying at 70 % leaves headroom for OS, page cache,
/// and other processes.
const SAFE_HEADROOM_PCT: u64 = 70;

/// Extra safety margin (fraction of estimated peak) added to the reported
/// `safety_margin_bytes` field so operators can plan conservatively.
const MARGIN_PCT: u64 = 25;

/// Recommended build strategy.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum BuildMode {
    /// Whole corpus fits in memory with headroom — fast path.
    InMemory,
    /// Corpus is too large for safe in-memory construction; use streaming/sharded path.
    Streaming,
}

impl std::fmt::Display for BuildMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BuildMode::InMemory => write!(f, "in-memory"),
            BuildMode::Streaming => write!(f, "streaming"),
        }
    }
}

/// Machine-readable output of the preflight RAM estimator.
#[derive(Debug, Serialize)]
pub struct MemoryEstimate {
    /// Raw corpus size in bytes (input).
    pub corpus_bytes: u64,
    /// Projected peak RAM for FM-index construction (bytes).
    pub estimated_peak_ram_bytes: u64,
    /// Extra headroom added on top of the peak estimate (bytes).
    pub safety_margin_bytes: u64,
    /// Total recommended available RAM (peak + margin).
    pub recommended_available_bytes: u64,
    /// Available RAM on this machine (bytes), or 0 if it could not be queried.
    pub available_ram_bytes: u64,
    /// Recommended build mode based on available RAM.
    pub recommended_mode: BuildMode,
    /// Human-readable rationale string.
    pub rationale: String,
}

/// Estimate peak RAM for building an FM-index from `corpus_bytes` bytes and
/// return a [`MemoryEstimate`] with a recommended [`BuildMode`].
pub fn estimate_memory(corpus_bytes: u64) -> MemoryEstimate {
    let estimated_peak = corpus_bytes.saturating_mul(PEAK_RAM_MULTIPLIER);
    let safety_margin = estimated_peak.saturating_mul(MARGIN_PCT) / 100;
    let recommended_available = estimated_peak.saturating_add(safety_margin);

    let available_ram = query_available_ram_bytes();

    let safe_limit = if available_ram > 0 {
        available_ram.saturating_mul(SAFE_HEADROOM_PCT) / 100
    } else {
        0
    };

    let recommended_mode = if available_ram == 0 || estimated_peak > safe_limit {
        BuildMode::Streaming
    } else {
        BuildMode::InMemory
    };

    let rationale = if available_ram == 0 {
        format!(
            "Available RAM could not be determined. Corpus {}  MiB → estimated peak {} MiB. \
             Defaulting to streaming mode for safety.",
            corpus_bytes / 1_048_576,
            estimated_peak / 1_048_576,
        )
    } else {
        format!(
            "Corpus {} MiB → estimated peak {} MiB ({}× multiplier). \
             Available RAM {} MiB, safe limit {} MiB ({} % headroom). \
             Recommended: {}.",
            corpus_bytes / 1_048_576,
            estimated_peak / 1_048_576,
            PEAK_RAM_MULTIPLIER,
            available_ram / 1_048_576,
            safe_limit / 1_048_576,
            SAFE_HEADROOM_PCT,
            recommended_mode,
        )
    };

    MemoryEstimate {
        corpus_bytes,
        estimated_peak_ram_bytes: estimated_peak,
        safety_margin_bytes: safety_margin,
        recommended_available_bytes: recommended_available,
        available_ram_bytes: available_ram,
        recommended_mode,
        rationale,
    }
}

/// Query total available physical RAM from the OS.
///
/// Returns 0 if the query fails (callers treat 0 as "unknown").
fn query_available_ram_bytes() -> u64 {
    #[cfg(target_os = "linux")]
    {
        query_linux_memavailable().unwrap_or(0)
    }
    #[cfg(target_os = "macos")]
    {
        query_macos_hw_memsize().unwrap_or(0)
    }
    #[cfg(not(any(target_os = "linux", target_os = "macos")))]
    {
        0
    }
}

#[cfg(target_os = "linux")]
fn query_linux_memavailable() -> Option<u64> {
    let content = std::fs::read_to_string("/proc/meminfo").ok()?;
    for line in content.lines() {
        // Prefer MemAvailable (free + reclaimable) over MemTotal.
        if line.starts_with("MemAvailable:") {
            let kb: u64 = line
                .split_whitespace()
                .nth(1)?
                .parse()
                .ok()?;
            return Some(kb * 1024);
        }
    }
    // Fallback to MemTotal.
    for line in content.lines() {
        if line.starts_with("MemTotal:") {
            let kb: u64 = line.split_whitespace().nth(1)?.parse().ok()?;
            return Some(kb * 1024);
        }
    }
    None
}

#[cfg(target_os = "macos")]
fn query_macos_hw_memsize() -> Option<u64> {
    let output = std::process::Command::new("sysctl")
        .args(["-n", "hw.memsize"])
        .output()
        .ok()?;
    let s = std::str::from_utf8(&output.stdout).ok()?.trim();
    s.parse().ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_corpus_returns_in_memory_or_streaming_not_panics() {
        let est = estimate_memory(0);
        assert_eq!(est.corpus_bytes, 0);
        assert_eq!(est.estimated_peak_ram_bytes, 0);
    }

    #[test]
    fn small_corpus_recommends_in_memory_when_ram_known() {
        // 1 MiB corpus → peak 6 MiB.  Any machine with > 9 MiB free should
        // recommend in-memory.  We only assert this when available_ram is known.
        let est = estimate_memory(1_048_576);
        if est.available_ram_bytes > 0 {
            assert_eq!(est.recommended_mode, BuildMode::InMemory);
        }
    }

    #[test]
    fn safety_margin_is_25_pct_of_peak() {
        let est = estimate_memory(100_000);
        let expected_peak = 100_000 * PEAK_RAM_MULTIPLIER;
        assert_eq!(est.estimated_peak_ram_bytes, expected_peak);
        assert_eq!(est.safety_margin_bytes, expected_peak * 25 / 100);
    }

    #[test]
    fn huge_corpus_recommends_streaming() {
        // Simulate a corpus that no real machine has RAM for: 128 TiB.
        let corpus_bytes = 128u64 * 1024 * 1024 * 1024 * 1024;
        let est = estimate_memory(corpus_bytes);
        assert_eq!(est.recommended_mode, BuildMode::Streaming);
    }

    #[test]
    fn estimate_is_serializable_to_json() {
        let est = estimate_memory(512 * 1024);
        let json = serde_json::to_string(&est).expect("serialize");
        assert!(json.contains("corpus_bytes"));
        assert!(json.contains("recommended_mode"));
        assert!(json.contains("rationale"));
    }
}
