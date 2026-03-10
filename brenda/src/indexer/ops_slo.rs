//! Operational SLO profile — capacity planning table + runbook helpers.
//!
//! Provides a static, machine-readable sizing table and a human-readable
//! runbook stub for operators deploying LOOM in production environments.

use serde::{Deserialize, Serialize};

// ── Sizing tier table ─────────────────────────────────────────────────────────

/// One row in the capacity-planning sizing table.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SizingTier {
    /// Human-readable tier label (e.g. "XS – demo").
    pub tier: &'static str,
    /// Corpus upper bound (bytes).
    pub max_corpus_bytes: u64,
    /// Corpus upper bound (human-readable).
    pub max_corpus_label: &'static str,
    /// Recommended available RAM (bytes).
    pub recommended_ram_bytes: u64,
    /// Recommended available RAM (human-readable).
    pub recommended_ram_label: &'static str,
    /// Recommended build mode for this tier.
    pub build_mode: &'static str,
    /// Recommended number of index shards (1 = monolithic).
    pub recommended_shards: usize,
    /// Target bytes per shard (0 = N/A for in-memory single-index).
    pub shard_size_bytes: u64,
    /// Approximate cold-search p99 latency target (milliseconds).
    pub search_p99_ms: u64,
    /// Notes for operators.
    pub notes: &'static str,
}

/// Static sizing table covering XS through XXL corpus tiers.
///
/// Based on the `PEAK_RAM_MULTIPLIER = 6` constant in `memory_profiler.rs`
/// plus a 25 % safety margin.
pub static SIZING_TABLE: &[SizingTier] = &[
    SizingTier {
        tier: "XS – demo",
        max_corpus_bytes: 10 * 1024 * 1024,        // 10 MiB
        max_corpus_label: "10 MiB",
        recommended_ram_bytes: 1 * 1024 * 1024 * 1024, // 1 GiB
        recommended_ram_label: "1 GiB",
        build_mode: "in-memory",
        recommended_shards: 1,
        shard_size_bytes: 0,
        search_p99_ms: 5,
        notes: "Laptop / CI; single monolithic index; no sharding required.",
    },
    SizingTier {
        tier: "S – small project",
        max_corpus_bytes: 100 * 1024 * 1024,       // 100 MiB
        max_corpus_label: "100 MiB",
        recommended_ram_bytes: 8 * 1024 * 1024 * 1024, // 8 GiB
        recommended_ram_label: "8 GiB",
        build_mode: "in-memory",
        recommended_shards: 1,
        shard_size_bytes: 0,
        search_p99_ms: 20,
        notes: "Developer workstation; in-memory path viable; keep swap off.",
    },
    SizingTier {
        tier: "M – medium codebase",
        max_corpus_bytes: 512 * 1024 * 1024,       // 512 MiB
        max_corpus_label: "512 MiB",
        recommended_ram_bytes: 32 * 1024 * 1024 * 1024, // 32 GiB
        recommended_ram_label: "32 GiB",
        build_mode: "in-memory or 2 shards",
        recommended_shards: 2,
        shard_size_bytes: 256 * 1024 * 1024,        // 256 MiB
        search_p99_ms: 50,
        notes: "Mono-repo or large library; 2 shards of 256 MiB each; \
                federated search adds ~5 ms fan-out overhead.",
    },
    SizingTier {
        tier: "L – large mono-repo",
        max_corpus_bytes: 4 * 1024 * 1024 * 1024,  // 4 GiB
        max_corpus_label: "4 GiB",
        recommended_ram_bytes: 128 * 1024 * 1024 * 1024, // 128 GiB
        recommended_ram_label: "128 GiB",
        build_mode: "streaming / 8 shards",
        recommended_shards: 8,
        shard_size_bytes: 512 * 1024 * 1024,        // 512 MiB
        search_p99_ms: 200,
        notes: "Cloud VM (r6i.4xlarge or equiv.); 8 shards; parallel fan-out \
                with ≤8 workers; SSD-backed index storage recommended.",
    },
    SizingTier {
        tier: "XL – organisation-wide",
        max_corpus_bytes: 32 * 1024 * 1024 * 1024, // 32 GiB
        max_corpus_label: "32 GiB",
        recommended_ram_bytes: 512 * 1024 * 1024 * 1024, // 512 GiB (u-series)
        recommended_ram_label: "512 GiB",
        build_mode: "streaming / 64 shards",
        recommended_shards: 64,
        shard_size_bytes: 512 * 1024 * 1024,        // 512 MiB
        search_p99_ms: 500,
        notes: "Bare-metal or memory-optimised cloud instance; distribute shards \
                across NVMe array; LOOM process pinned to NUMA node.",
    },
    SizingTier {
        tier: "XXL – DNA / genomic corpus",
        max_corpus_bytes: 256 * 1024 * 1024 * 1024, // 256 GiB
        max_corpus_label: "256 GiB",
        recommended_ram_bytes: 2 * 1024 * 1024 * 1024 * 1024, // 2 TiB
        recommended_ram_label: "2 TiB",
        build_mode: "streaming / 512 shards",
        recommended_shards: 512,
        shard_size_bytes: 512 * 1024 * 1024,        // 512 MiB
        search_p99_ms: 2000,
        notes: "Genomic or telemetry scale; streaming ingest via dna-corpus; \
                shards distributed across storage nodes; requires cluster coordination.",
    },
];

// ── Runbook ───────────────────────────────────────────────────────────────────

/// Short operational runbook sections for common failure scenarios.
#[derive(Debug, Serialize, Deserialize)]
pub struct RunbookEntry {
    /// Failure scenario label.
    pub scenario: &'static str,
    /// Symptoms to look for.
    pub symptoms: &'static str,
    /// Recommended recovery steps (newline-separated).
    pub recovery_steps: &'static str,
    /// Prevention guidance.
    pub prevention: &'static str,
}

/// Static runbook for common LOOM operational failures.
pub static RUNBOOK: &[RunbookEntry] = &[
    RunbookEntry {
        scenario: "OOM during index build",
        symptoms: "Process killed (exit 137 / SIGKILL); `dmesg` shows OOM killer invocation.",
        recovery_steps: "\
1. Confirm corpus size with `du -sb <corpus_dir>`.\n\
2. Run `loom index --estimate-memory` to check projected peak RAM.\n\
3. If build-mode was `in-memory`, switch to `--build-mode streaming`.\n\
4. If corpus > 100 MiB, rebuild with `--shards` and a target `--shard-size 268435456` (256 MiB).\n\
5. Re-run the build; each shard is independently resumable — delete partial shard files and rerun.",
        prevention: "Always run `loom index --estimate-memory` before indexing corpora > 1 GiB. \
                     Set `--build-mode auto` (default) so the estimator gates the build path.",
    },
    RunbookEntry {
        scenario: "Corrupt or incomplete shard file",
        symptoms: "FederatedSearcher returns fewer results than expected; \
                   `loom shard-search --count` disagrees with single-index baseline.",
        recovery_steps: "\
1. Open `shard_manifest.json` and note the SHA-256 for each shard entry.\n\
2. Re-compute checksums: `sha256sum <shard_dir>/shard_*.idx`.\n\
3. Identify shards where checksum mismatches manifest.\n\
4. Delete the affected shard files.\n\
5. Re-run `loom index --shards` — the builder overwrites only the missing/deleted shards.",
        prevention: "Store `shard_manifest.json` in version control or object storage so the \
                     authoritative checksum set is always available. Validate checksums as part \
                     of deployment CI.",
    },
    RunbookEntry {
        scenario: "Federated search latency spike",
        symptoms: "p99 latency > 2× baseline; fan-out threads serialise on a single slow shard.",
        recovery_steps: "\
1. Identify the slow shard: time `loom search --index shard_NNNN.idx <pattern>`.\n\
2. Check if the slow shard is on a degraded or cold storage device.\n\
3. Reduce `--workers` to avoid saturating I/O bandwidth.\n\
4. If one shard is disproportionately large, rebuild with a smaller `--shard-size`.",
        prevention: "Keep shard sizes uniform (use a fixed `--shard-size`). \
                     Pre-warm shard files with `fadvise WILLNEED` or equivalent before query bursts.",
    },
    RunbookEntry {
        scenario: "Disk full during shard build",
        symptoms: "Build exits with `IO error: No space left on device`; partial shards on disk.",
        recovery_steps: "\
1. Free disk space or mount a larger volume.\n\
2. Delete incomplete shard files in the output directory (`shard_*.idx` with size 0).\n\
3. Re-run `loom index --shards`; completed shards will be overwritten only if their file is \
   missing or corrupt (validate via manifest checksum).",
        prevention: "Ensure `2 × corpus_size` free space on the target volume before building. \
                     Rule of thumb: each shard file ≈ 1.2× its shard corpus bytes.",
    },
];

// ── Machine-readable profile output ──────────────────────────────────────────

/// Full ops SLO profile: sizing table + runbook.
#[derive(Debug, Serialize, Deserialize)]
pub struct OpsProfile {
    pub version: &'static str,
    pub description: &'static str,
    pub sizing_table: Vec<SizingTierOwned>,
    pub runbook: Vec<RunbookEntryOwned>,
}

/// Owned (heap-allocated) version of [`SizingTier`] for JSON serialization.
#[derive(Debug, Serialize, Deserialize)]
pub struct SizingTierOwned {
    pub tier: String,
    pub max_corpus_label: String,
    pub recommended_ram_label: String,
    pub build_mode: String,
    pub recommended_shards: usize,
    pub shard_size_bytes: u64,
    pub search_p99_ms: u64,
    pub notes: String,
}

/// Owned (heap-allocated) version of [`RunbookEntry`] for JSON serialization.
#[derive(Debug, Serialize, Deserialize)]
pub struct RunbookEntryOwned {
    pub scenario: String,
    pub symptoms: String,
    pub recovery_steps: String,
    pub prevention: String,
}

/// Build and return the full [`OpsProfile`] struct.
pub fn ops_profile() -> OpsProfile {
    OpsProfile {
        version: "v1",
        description: "LOOM operational SLO profile — capacity planning + runbook",
        sizing_table: SIZING_TABLE
            .iter()
            .map(|t| SizingTierOwned {
                tier: t.tier.to_string(),
                max_corpus_label: t.max_corpus_label.to_string(),
                recommended_ram_label: t.recommended_ram_label.to_string(),
                build_mode: t.build_mode.to_string(),
                recommended_shards: t.recommended_shards,
                shard_size_bytes: t.shard_size_bytes,
                search_p99_ms: t.search_p99_ms,
                notes: t.notes.to_string(),
            })
            .collect(),
        runbook: RUNBOOK
            .iter()
            .map(|r| RunbookEntryOwned {
                scenario: r.scenario.to_string(),
                symptoms: r.symptoms.to_string(),
                recovery_steps: r.recovery_steps.to_string(),
                prevention: r.prevention.to_string(),
            })
            .collect(),
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sizing_table_is_non_empty() {
        assert!(!SIZING_TABLE.is_empty());
    }

    #[test]
    fn sizing_table_tiers_are_ascending() {
        let bytes: Vec<u64> = SIZING_TABLE.iter().map(|t| t.max_corpus_bytes).collect();
        for w in bytes.windows(2) {
            assert!(w[0] < w[1], "sizing table must be sorted by ascending corpus size");
        }
    }

    #[test]
    fn runbook_has_at_least_three_entries() {
        assert!(RUNBOOK.len() >= 3);
    }

    #[test]
    fn ops_profile_serializes_to_json() {
        let profile = ops_profile();
        let json = serde_json::to_string(&profile).expect("serialize");
        assert!(json.contains("sizing_table"));
        assert!(json.contains("runbook"));
        assert!(json.contains("version"));
    }

    #[test]
    fn sharded_tiers_have_shard_size() {
        for tier in SIZING_TABLE {
            if tier.recommended_shards > 1 {
                assert!(
                    tier.shard_size_bytes > 0,
                    "tier '{}' recommends {} shards but shard_size_bytes is 0",
                    tier.tier,
                    tier.recommended_shards
                );
            }
        }
    }

    #[test]
    fn search_p99_latencies_are_monotonic() {
        let latencies: Vec<u64> = SIZING_TABLE.iter().map(|t| t.search_p99_ms).collect();
        for w in latencies.windows(2) {
            assert!(w[0] <= w[1], "p99 latency must be non-decreasing with tier size");
        }
    }
}
