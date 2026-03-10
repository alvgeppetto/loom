//! Sharded index builder — large-file-safe index construction.
//!
//! Splits a corpus into fixed-size chunks, builds one FM-index shard per
//! chunk, persists each shard to `<output_dir>/shard_NNNN.idx`, and writes a
//! JSON manifest (`<output_dir>/shard_manifest.json`) with deterministic
//! ordering and SHA-256 checksums.

use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use super::corpus::{Corpus, SourceFile, FILE_DELIMITER};
use super::{IndexBuilder, IndexError};

// ── Public types ─────────────────────────────────────────────────────────────

/// One entry in the shard manifest.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ShardEntry {
    /// Zero-based shard index (stable ordering).
    pub id: usize,
    /// Shard file name, relative to the manifest directory.
    pub path: String,
    /// Byte offset into the *original* corpus where this shard begins.
    pub corpus_start_offset: u64,
    /// Byte offset into the *original* corpus where this shard ends (exclusive).
    pub corpus_end_offset: u64,
    /// Number of source files in this shard.
    pub file_count: usize,
    /// SHA-256 hex digest of the shard `.idx` file.
    pub sha256: String,
}

/// Manifest describing all shards produced for a corpus.
#[derive(Debug, Serialize, Deserialize)]
pub struct ShardManifest {
    /// Ordered shard entries (ascending shard id).
    pub shards: Vec<ShardEntry>,
    /// Total corpus size in bytes (sum of all shard corpora).
    pub total_corpus_bytes: u64,
    /// Target maximum bytes per shard (actual shards may be smaller).
    pub target_shard_bytes: u64,
    /// Total number of source files across all shards.
    pub total_files: usize,
    /// ISO-8601 UTC creation timestamp.
    pub created_at: String,
}

impl ShardManifest {
    /// Load a manifest from a JSON file.
    pub fn load(path: &Path) -> Result<Self, IndexError> {
        let data = fs::read(path)?;
        serde_json::from_slice(&data).map_err(|e| IndexError::Serialization(e.to_string()))
    }

    /// Persist this manifest as a JSON file.
    pub fn save(&self, path: &Path) -> Result<(), IndexError> {
        let json =
            serde_json::to_vec_pretty(self).map_err(|e| IndexError::Serialization(e.to_string()))?;
        fs::write(path, json)?;
        Ok(())
    }

    /// Return the absolute path to shard `id`'s index file, given that the
    /// manifest lives in `base_dir`.
    pub fn shard_path(&self, id: usize, base_dir: &Path) -> Option<PathBuf> {
        self.shards.get(id).map(|e| base_dir.join(&e.path))
    }
}

// ── Core builder ─────────────────────────────────────────────────────────────

/// Build sharded indices from `corpus`, placing shard files and manifest in
/// `output_dir`.
///
/// Files from the corpus are grouped so that the concatenated content per
/// shard does not exceed `target_shard_bytes`.  Each shard is saved to
/// `shard_NNNN.idx`; the manifest is saved to `shard_manifest.json` in the
/// same directory.
pub fn build_shards(
    corpus: &Corpus,
    output_dir: &Path,
    target_shard_bytes: usize,
) -> Result<ShardManifest, IndexError> {
    fs::create_dir_all(output_dir)?;

    let groups = partition_files(&corpus.files, target_shard_bytes);
    let mut shard_entries: Vec<ShardEntry> = Vec::with_capacity(groups.len());

    for (id, files) in groups.iter().enumerate() {
        let (sub_corpus, corpus_start, corpus_end) = build_sub_corpus(files);
        let builder = IndexBuilder::new(&sub_corpus)?;

        let shard_filename = format!("shard_{:04}.idx", id);
        let shard_path = output_dir.join(&shard_filename);
        builder.save(&shard_path)?;

        let sha256 = sha256_file(&shard_path)?;

        shard_entries.push(ShardEntry {
            id,
            path: shard_filename,
            corpus_start_offset: corpus_start as u64,
            corpus_end_offset: corpus_end as u64,
            file_count: files.len(),
            sha256,
        });
    }

    let manifest = ShardManifest {
        total_corpus_bytes: corpus.size() as u64,
        target_shard_bytes: target_shard_bytes as u64,
        total_files: corpus.files.len(),
        created_at: utc_now_rfc3339(),
        shards: shard_entries,
    };

    let manifest_path = output_dir.join("shard_manifest.json");
    manifest.save(&manifest_path)?;

    Ok(manifest)
}

// ── Internal helpers ──────────────────────────────────────────────────────────

/// Partition corpus files into groups with total content ≤ target_shard_bytes.
///
/// Preserves the original file order.  A single file larger than the target
/// is placed in its own shard rather than being split.
fn partition_files<'a>(
    files: &'a [SourceFile],
    target_shard_bytes: usize,
) -> Vec<Vec<&'a SourceFile>> {
    let mut groups: Vec<Vec<&'a SourceFile>> = Vec::new();
    let mut current: Vec<&'a SourceFile> = Vec::new();
    let mut current_bytes: usize = 0;

    for file in files {
        let file_size = file.end_offset.saturating_sub(file.start_offset);
        if !current.is_empty() && current_bytes + file_size > target_shard_bytes {
            groups.push(std::mem::take(&mut current));
            current_bytes = 0;
        }
        current.push(file);
        current_bytes += file_size;
    }

    if !current.is_empty() {
        groups.push(current);
    }

    groups
}

/// Rebuild a [`Corpus`] from a subset of files with contiguous offsets from 0.
///
/// Returns `(sub_corpus, original_start_offset, original_end_offset)`.
fn build_sub_corpus(files: &[&SourceFile]) -> (Corpus, usize, usize) {
    let mut concatenated: Vec<u8> = Vec::new();
    let mut sub_files: Vec<SourceFile> = Vec::with_capacity(files.len());

    let original_start = files.first().map(|f| f.start_offset).unwrap_or(0);
    let original_end = files.last().map(|f| f.end_offset).unwrap_or(0);

    for file in files {
        let start_offset = concatenated.len();
        concatenated.push(FILE_DELIMITER);
        concatenated.extend_from_slice(file.path.to_string_lossy().as_bytes());
        concatenated.push(FILE_DELIMITER);
        concatenated.extend_from_slice(file.content.as_bytes());
        let end_offset = concatenated.len();

        sub_files.push(SourceFile {
            path: file.path.clone(),
            content: file.content.clone(),
            start_offset,
            end_offset,
        });
    }

    let sub = Corpus {
        files: sub_files,
        concatenated,
    };
    (sub, original_start, original_end)
}

/// Compute the SHA-256 hex digest of a file on disk.
fn sha256_file(path: &Path) -> Result<String, IndexError> {
    let mut file = fs::File::open(path)?;
    let mut hasher = Sha256::new();
    let mut buf = [0u8; 65536];
    loop {
        let n = file.read(&mut buf)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

/// Return a naive RFC-3339 UTC timestamp string derived from `SystemTime`.
fn utc_now_rfc3339() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);

    let s = secs;
    let (year, month, day) = days_to_ymd(s / 86400);
    let hour = (s % 86400) / 3600;
    let min = (s % 3600) / 60;
    let sec = s % 60;
    format!("{:04}-{:02}-{:02}T{:02}:{:02}:{:02}Z", year, month, day, hour, min, sec)
}

fn days_to_ymd(mut days: u64) -> (u64, u64, u64) {
    let mut year = 1970u64;
    loop {
        let dy = if is_leap(year) { 366 } else { 365 };
        if days < dy {
            break;
        }
        days -= dy;
        year += 1;
    }
    let month_days: [u64; 12] = if is_leap(year) {
        [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    } else {
        [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    };
    let mut month = 1u64;
    for &md in &month_days {
        if days < md {
            break;
        }
        days -= md;
        month += 1;
    }
    (year, month, days + 1)
}

fn is_leap(year: u64) -> bool {
    (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    use tempfile::TempDir;

    fn make_file(name: &str, content: &str, start: usize) -> SourceFile {
        let path = PathBuf::from(name);
        let marker = format!("\x01{}\x01{}", name, content);
        let end = start + marker.len();
        SourceFile {
            path,
            content: content.to_string(),
            start_offset: start,
            end_offset: end,
        }
    }

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

    #[test]
    fn build_shards_single_shard_small_corpus() {
        let corpus = tiny_corpus(&[("a.rs", "fn a() {}"), ("b.rs", "fn b() {}")]);
        let dir = TempDir::new().unwrap();
        let manifest = build_shards(&corpus, dir.path(), 1_000_000).unwrap();

        assert_eq!(manifest.shards.len(), 1);
        assert_eq!(manifest.total_files, 2);
        assert_eq!(manifest.shards[0].file_count, 2);

        let shard_file = dir.path().join("shard_0000.idx");
        assert!(shard_file.exists(), "shard_0000.idx should exist");

        let manifest_file = dir.path().join("shard_manifest.json");
        assert!(manifest_file.exists(), "shard_manifest.json should exist");
    }

    #[test]
    fn build_shards_multiple_shards_small_limit() {
        let corpus = tiny_corpus(&[
            ("a.rs", "fn a() {}"),
            ("b.rs", "fn b() {}"),
            ("c.rs", "fn c() {}"),
        ]);
        let dir = TempDir::new().unwrap();
        // Set a tiny shard limit so each file gets its own shard.
        let manifest = build_shards(&corpus, dir.path(), 1).unwrap();

        assert_eq!(manifest.shards.len(), 3);
        assert_eq!(manifest.total_files, 3);
        for entry in &manifest.shards {
            assert_eq!(entry.file_count, 1);
            assert!(!entry.sha256.is_empty());
        }
    }

    #[test]
    fn manifest_roundtrip() {
        let corpus = tiny_corpus(&[("x.rs", "let x = 1;")]);
        let dir = TempDir::new().unwrap();
        let original = build_shards(&corpus, dir.path(), 512 * 1024).unwrap();

        let manifest_path = dir.path().join("shard_manifest.json");
        let loaded = ShardManifest::load(&manifest_path).unwrap();

        assert_eq!(loaded.shards.len(), original.shards.len());
        assert_eq!(loaded.total_files, original.total_files);
        assert_eq!(loaded.shards[0].sha256, original.shards[0].sha256);
    }

    #[test]
    fn shard_is_searchable() {
        let corpus = tiny_corpus(&[("main.rs", "fn hello_world() {}")]);
        let dir = TempDir::new().unwrap();
        build_shards(&corpus, dir.path(), 1_000_000).unwrap();

        let shard_path = dir.path().join("shard_0000.idx");
        let builder = IndexBuilder::load(&shard_path).unwrap();
        let hits = builder.search("hello_world");
        assert!(!hits.is_empty(), "expected at least one hit for 'hello_world'");
    }

    #[test]
    fn shard_ordering_is_deterministic() {
        let corpus = tiny_corpus(&[
            ("one.rs", "one"),
            ("two.rs", "two"),
            ("three.rs", "three"),
        ]);
        let dir = TempDir::new().unwrap();
        let manifest = build_shards(&corpus, dir.path(), 1).unwrap();

        let ids: Vec<usize> = manifest.shards.iter().map(|e| e.id).collect();
        let sorted: Vec<usize> = {
            let mut v = ids.clone();
            v.sort_unstable();
            v
        };
        assert_eq!(ids, sorted, "shard ids should be in ascending order");
    }

    #[test]
    fn sha256_is_hex_64_chars() {
        let corpus = tiny_corpus(&[("f.rs", "content")]);
        let dir = TempDir::new().unwrap();
        let manifest = build_shards(&corpus, dir.path(), 1_000_000).unwrap();
        let hex = &manifest.shards[0].sha256;
        assert_eq!(hex.len(), 64, "SHA-256 hex digest should be 64 characters");
        assert!(hex.chars().all(|c| c.is_ascii_hexdigit()));
    }

    #[test]
    fn partition_files_respects_target() {
        let files: Vec<SourceFile> = vec![
            make_file("a.rs", &"a".repeat(100), 0),
            make_file("b.rs", &"b".repeat(100), 110),
            make_file("c.rs", &"c".repeat(100), 220),
        ];
        // Target of 150 bytes → groups of 1 file each since 110 + 110 > 150.
        let groups = partition_files(&files, 150);
        assert_eq!(groups.len(), 3);
    }
}
