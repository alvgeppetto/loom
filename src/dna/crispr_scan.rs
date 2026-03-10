//! Fast CRISPR target scanner.
//!
//! Single-pass sliding-window extraction of k-mers with 2-bit encoding,
//! PAM classification, and parallel chunk processing via rayon.
//! Matches the algorithm described in the CRISPR paper:
//!   publications/papers/crispr-novel-targets-preprint.tex

use rayon::prelude::*;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::time::Instant;

/// PAM site classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PamType {
    /// SpCas9: last 2 bases of 23-mer are GG
    Ngg,
    /// Cas12a: first 3 bases of 23-mer are TTT
    Tttn,
    /// No recognized PAM
    None,
}

impl PamType {
    fn label(self) -> &'static str {
        match self {
            PamType::Ngg => "NGG",
            PamType::Tttn => "TTTN",
            PamType::None => "none",
        }
    }
}

/// A discovered CRISPR target.
pub struct CrisprTarget {
    pub kmer: u64,
    pub first_pos: u64,
    pub occurrences: u32,
    pub pam: PamType,
}

// ── 2-bit DNA encoding ──────────────────────────────────

/// Encode a single base to 2 bits. Returns None for non-ACGT.
#[inline(always)]
fn encode_base(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => Option::None,
    }
}

/// Decode a 2-bit encoded k-mer back to a DNA string.
fn decode_kmer(kmer: u64, k: usize) -> String {
    let mut s = String::with_capacity(k);
    for i in (0..k).rev() {
        let bits = (kmer >> (i * 2)) & 3;
        s.push(match bits {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => unreachable!(),
        });
    }
    s
}

/// Classify PAM site from a 2-bit encoded k-mer.
#[inline(always)]
fn classify_pam(kmer: u64, k: usize) -> PamType {
    // SpCas9 NGG: last 2 bases are GG (G = 2)
    let last = kmer & 3;
    let second_last = (kmer >> 2) & 3;
    if second_last == 2 && last == 2 {
        return PamType::Ngg;
    }

    // Cas12a TTTN: first 3 bases are TTT (T = 3)
    let shift_first = (k as u32 - 1) * 2;
    let first = (kmer >> shift_first) & 3;
    let second = (kmer >> (shift_first - 2)) & 3;
    let third = (kmer >> (shift_first - 4)) & 3;
    if first == 3 && second == 3 && third == 3 {
        return PamType::Tttn;
    }

    PamType::None
}

// ── FASTA reading ───────────────────────────────────────

/// Read FASTA file(s) and return concatenated uppercase DNA bytes (ACGT only, N preserved).
/// Only use for files that fit in memory.
fn read_fasta_sequences(input: &Path) -> anyhow::Result<Vec<u8>> {
    let paths = collect_fasta_paths(input)?;
    let mut seq = Vec::new();
    for path in &paths {
        let file = std::fs::File::open(path)?;
        let reader = BufReader::with_capacity(1 << 20, file); // 1 MB buffer
        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('>') || trimmed.starts_with(';') {
                continue;
            }
            seq.extend(trimmed.bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    Ok(seq)
}

/// Collect sorted FASTA file paths from a file or directory.
fn collect_fasta_paths(input: &Path) -> anyhow::Result<Vec<std::path::PathBuf>> {
    if input.is_dir() {
        let mut ps: Vec<_> = std::fs::read_dir(input)?
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .filter(|p| {
                matches!(
                    p.extension().and_then(|e| e.to_str()),
                    Some("fa" | "fasta" | "fna")
                )
            })
            .collect();
        ps.sort();
        Ok(ps)
    } else {
        Ok(vec![input.to_path_buf()])
    }
}

/// Total size of all FASTA files in bytes.
fn total_fasta_size(input: &Path) -> anyhow::Result<u64> {
    let paths = collect_fasta_paths(input)?;
    let mut total = 0u64;
    for p in &paths {
        total += std::fs::metadata(p)?.len();
    }
    Ok(total)
}

// ── Streaming scanner for large files ───────────────────

/// Stream-scan FASTA files in fixed-size chunks to limit memory usage.
/// Returns the merged global k-mer map.
fn stream_scan_fasta(
    input: &Path,
    k: usize,
    chunk_capacity: usize,
) -> anyhow::Result<HashMap<u64, (u64, u32)>> {
    let paths = collect_fasta_paths(input)?;
    let num_threads = rayon::current_num_threads();

    let mut global_map: HashMap<u64, (u64, u32)> = HashMap::new();
    let mut buf = Vec::with_capacity(chunk_capacity + 1024);
    let mut global_offset: u64 = 0;
    let mut chunk_num: usize = 0;
    let mut total_seq_bytes: u64 = 0;

    let flush = |buf: &[u8],
                 global_map: &mut HashMap<u64, (u64, u32)>,
                 base_offset: u64,
                 num_threads: usize,
                 k: usize,
                 chunk_num: usize| {
        let seq_len = buf.len();
        if seq_len < k {
            return;
        }
        let t = Instant::now();
        let partial_maps: Vec<_> = (0..num_threads)
            .into_par_iter()
            .map(|i| {
                let start = i * (seq_len / num_threads);
                let end = if i == num_threads - 1 {
                    seq_len
                } else {
                    ((i + 1) * (seq_len / num_threads) + k - 1).min(seq_len)
                };
                scan_chunk(&buf[start..end], k, base_offset + start as u64)
            })
            .collect();

        for partial in partial_maps {
            for (kmer, (pos, count)) in partial {
                let entry = global_map.entry(kmer).or_insert((pos, 0u32));
                entry.1 += count;
                if pos < entry.0 {
                    entry.0 = pos;
                }
            }
        }
        eprintln!(
            "[crispr-scan]   chunk {} ({:.1} MB, {} unique so far) scanned in {:.2}s",
            chunk_num,
            seq_len as f64 / 1_048_576.0,
            global_map.len(),
            t.elapsed().as_secs_f64(),
        );

        // Prune singletons to bound memory. Safe for top-N queries:
        // k-mers appearing only once can't be in the top 10,000.
        // Cross-chunk count-2 k-mers may be lost, but won't rank highly.
        if global_map.len() > 20_000_000 {
            let before = global_map.len();
            global_map.retain(|_, (_, count)| *count > 1);
            eprintln!(
                "[crispr-scan]   pruned singletons: {} → {} entries",
                before,
                global_map.len(),
            );
        }
    };

    for path in &paths {
        let file = std::fs::File::open(path)?;
        let reader = BufReader::with_capacity(1 << 20, file);
        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('>') || trimmed.starts_with(';') {
                continue;
            }
            buf.extend(trimmed.bytes().map(|b| b.to_ascii_uppercase()));

            if buf.len() >= chunk_capacity {
                flush(&buf, &mut global_map, global_offset, num_threads, k, chunk_num);
                total_seq_bytes += buf.len() as u64;
                global_offset += buf.len() as u64;
                chunk_num += 1;
                // Keep last k-1 bytes as overlap for next chunk
                let overlap = buf[buf.len() - (k - 1)..].to_vec();
                buf.clear();
                buf.extend_from_slice(&overlap);
            }
        }
    }
    // Flush remaining
    if buf.len() >= k {
        flush(&buf, &mut global_map, global_offset, num_threads, k, chunk_num);
        total_seq_bytes += buf.len() as u64;
    }

    eprintln!(
        "[crispr-scan] Streamed {:.1} GB total sequence, {} unique {}-mers",
        total_seq_bytes as f64 / 1_073_741_824.0,
        global_map.len(),
        k,
    );

    Ok(global_map)
}

/// Count k-mers in a byte slice. Returns HashMap<kmer_u64, (first_position, count)>.
fn scan_chunk(seq: &[u8], k: usize, base_offset: u64) -> HashMap<u64, (u64, u32)> {
    let mask: u64 = (1u64 << (k as u64 * 2)) - 1;
    let mut map: HashMap<u64, (u64, u32)> = HashMap::new();
    let mut kmer: u64 = 0;
    let mut valid: usize = 0;

    for (i, &base) in seq.iter().enumerate() {
        if let Some(enc) = encode_base(base) {
            kmer = ((kmer << 2) | enc) & mask;
            valid += 1;
            if valid >= k {
                let pos = base_offset + i as u64 + 1 - k as u64;
                let entry = map.entry(kmer).or_insert((pos, 0));
                entry.1 += 1;
                if pos < entry.0 {
                    entry.0 = pos;
                }
            }
        } else {
            // N or ambiguous base: reset window
            valid = 0;
            kmer = 0;
        }
    }

    map
}

/// Merge partial k-mer maps, keeping minimum position and summing counts.
fn merge_maps(
    maps: Vec<HashMap<u64, (u64, u32)>>,
) -> HashMap<u64, (u64, u32)> {
    let total_keys: usize = maps.iter().map(|m| m.len()).sum();
    let mut merged = HashMap::with_capacity(total_keys);
    for partial in maps {
        for (kmer, (pos, count)) in partial {
            let entry = merged.entry(kmer).or_insert((pos, 0u32));
            entry.1 += count;
            if pos < entry.0 {
                entry.0 = pos;
            }
        }
    }
    merged
}

// ── Public API ──────────────────────────────────────────

/// Run a full CRISPR target scan.
///
/// Reads FASTA from `input`, scans for k-mers of size `k`,
/// classifies PAM sites, and writes results to CSV at `output`.
pub fn run_scan(
    input: &Path,
    output: &Path,
    k: usize,
    max_targets: usize,
    pathogen_name: &str,
) -> anyhow::Result<()> {
    anyhow::ensure!(k >= 20 && k <= 32, "k-mer size must be 20..32 (got {k})");

    let t0 = Instant::now();
    let file_size = total_fasta_size(input)?;
    let streaming_threshold = 500 * 1_048_576u64; // 500 MB

    let (merged, scan_time) = if file_size > streaming_threshold {
        // Streaming mode: process in 256 MB chunks to stay within RAM
        let chunk_capacity = 256 * 1_048_576usize; // 256 MB per chunk
        eprintln!(
            "[crispr-scan] Streaming mode ({:.1} GB input, {} MB chunks)",
            file_size as f64 / 1_073_741_824.0,
            chunk_capacity / 1_048_576,
        );
        let t1 = Instant::now();
        let map = stream_scan_fasta(input, k, chunk_capacity)?;
        let elapsed = t1.elapsed().as_secs_f64();
        (map, elapsed)
    } else {
        // In-memory mode: load everything, parallel scan
        eprintln!("[crispr-scan] Reading FASTA from {}...", input.display());
        let seq = read_fasta_sequences(input)?;
        let seq_len = seq.len();
        eprintln!(
            "[crispr-scan] Loaded {:.1} MB of sequence in {:.2}s",
            seq_len as f64 / 1_048_576.0,
            t0.elapsed().as_secs_f64()
        );

        let t1 = Instant::now();
        let num_threads = rayon::current_num_threads();
        let partial_maps: Vec<_> = (0..num_threads)
            .into_par_iter()
            .map(|i| {
                let start = i * (seq_len / num_threads);
                let end = if i == num_threads - 1 {
                    seq_len
                } else {
                    ((i + 1) * (seq_len / num_threads) + k - 1).min(seq_len)
                };
                scan_chunk(&seq[start..end], k, start as u64)
            })
            .collect();

        let map = merge_maps(partial_maps);
        let elapsed = t1.elapsed().as_secs_f64();
        (map, elapsed)
    };

    let unique_kmers = merged.len();
    eprintln!(
        "[crispr-scan] Scanned {} unique {}-mers in {:.2}s",
        unique_kmers, k, scan_time,
    );

    // 3. Classify PAM and build target list
    let mut targets: Vec<CrisprTarget> = merged
        .into_iter()
        .map(|(kmer, (first_pos, occurrences))| CrisprTarget {
            kmer,
            first_pos,
            occurrences,
            pam: classify_pam(kmer, k),
        })
        .collect();

    // Sort by occurrences descending
    targets.sort_unstable_by(|a, b| b.occurrences.cmp(&a.occurrences));

    // Cap at max_targets
    if targets.len() > max_targets {
        eprintln!(
            "[crispr-scan] Capping from {} to {} targets",
            targets.len(),
            max_targets
        );
        targets.truncate(max_targets);
    }

    // 4. Write CSV
    let t2 = Instant::now();
    if let Some(parent) = output.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut f = std::io::BufWriter::new(std::fs::File::create(output)?);
    writeln!(
        f,
        "pathogen,gene,position,sequence_23mer,pam_type,source_motif,occurrences"
    )?;
    let gene_label = format!("{pathogen_name}.fna");
    for t in &targets {
        let seq_str = decode_kmer(t.kmer, k);
        writeln!(
            f,
            "{},{},{},{},{},sliding_window,{}",
            pathogen_name,
            gene_label,
            t.first_pos,
            seq_str,
            t.pam.label(),
            t.occurrences,
        )?;
    }
    f.flush()?;

    let total = t0.elapsed();
    eprintln!(
        "[crispr-scan] Done: {} targets → {} in {:.2}s (write: {:.2}s)",
        targets.len(),
        output.display(),
        total.as_secs_f64(),
        t2.elapsed().as_secs_f64(),
    );
    let pam_ngg = targets.iter().filter(|t| t.pam == PamType::Ngg).count();
    let pam_tttn = targets.iter().filter(|t| t.pam == PamType::Tttn).count();
    let pam_none = targets.iter().filter(|t| t.pam == PamType::None).count();
    eprintln!(
        "[crispr-scan] PAM breakdown: NGG={}, TTTN={}, none={}",
        pam_ngg, pam_tttn, pam_none,
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encode_decode_roundtrip() {
        let seq = "ACGTACGTACGTACGTACGTACG"; // 23-mer
        let k = seq.len();
        let mask: u64 = (1u64 << (k as u64 * 2)) - 1;
        let mut kmer: u64 = 0;
        for b in seq.bytes() {
            kmer = ((kmer << 2) | encode_base(b).unwrap()) & mask;
        }
        assert_eq!(decode_kmer(kmer, k), seq);
    }

    #[test]
    fn pam_ngg_detected() {
        // 23-mer ending in GG
        let seq = b"ACGTACGTACGTACGTACGTACGG";
        let k = 23;
        let mask: u64 = (1u64 << (k as u64 * 2)) - 1;
        let mut kmer: u64 = 0;
        // only use last 23 chars
        for &b in &seq[seq.len() - k..] {
            kmer = ((kmer << 2) | encode_base(b).unwrap()) & mask;
        }
        assert_eq!(classify_pam(kmer, k), PamType::Ngg);
    }

    #[test]
    fn pam_tttn_detected() {
        // 23-mer starting with TTT
        let seq = b"TTTACGTACGTACGTACGTACGT";
        let k = 23;
        let mask: u64 = (1u64 << (k as u64 * 2)) - 1;
        let mut kmer: u64 = 0;
        for &b in seq.iter() {
            kmer = ((kmer << 2) | encode_base(b).unwrap()) & mask;
        }
        assert_eq!(classify_pam(kmer, k), PamType::Tttn);
    }

    #[test]
    fn scan_small_sequence() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGT"; // 28 bp
        let map = scan_chunk(seq, 23, 0);
        // Should have 6 unique 23-mers (28 - 23 + 1 = 6)
        assert_eq!(map.values().map(|(_, c)| *c as usize).sum::<usize>(), 6);
    }

    #[test]
    fn ambiguous_bases_reset_window() {
        let seq = b"ACGTACGTACGTNACGTACGTACGT"; // N at position 12
        let map = scan_chunk(seq, 23, 0);
        // The N splits the sequence: 13 bases before (too short for 23-mer),
        // 12 bases after (too short for 23-mer)
        assert!(map.is_empty());
    }
}
