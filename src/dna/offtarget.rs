//! Off-target analysis: pigeonhole seeding + 2-bit Hamming verification.
//!
//! # Algorithm
//!
//! For a 20-nt guide with at most `max_mismatches` edits against a target,
//! the **pigeonhole principle** guarantees that if we split the guide into
//! `max_mismatches + 1` non-overlapping seeds, at least one seed must be an
//! exact match. With max_mismatches = 3 and 4 seeds of 5 bases each:
//!
//!   0         5        10        15       20
//!   ┣━seed₀━━┫━seed₁━━┫━seed₂━━┫━seed₃━━┫
//!
//! We scan the genome once, extracting every 5-mer. For each 5-mer at
//! position `p`, we look it up in a seed hash keyed by (seed_value,
//! seed_offset). A hit at (seed_val, soff) means a guide could start at
//! `p − soff`; we then verify the full 20-mer with a 2-bit Hamming check —
//! 4 ARM64 instructions (XOR → OR-shift → AND → CNT).
//!
//! Forward and reverse-complement of every guide are indexed separately so
//! both strands are covered in a single genome pass.
//!
//! Parallelism: Rayon partition over fixed-size sub-chunks within each
//! chromosome sequence (no cross-chunk boundary issue since guide starts
//! are uniquely assigned to chunks).

use rayon::prelude::*;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::time::Instant;

/// Length of each CRISPR spacer we consider (no PAM).
pub const GUIDE_LEN: usize = 20;
/// Seed length for pigeonhole filtering.
const SEED_LEN: usize = 5;
/// Number of non-overlapping seeds covering the full guide.
const NUM_SEEDS: usize = GUIDE_LEN / SEED_LEN; // 4

// ── 2-bit DNA encoding ─────────────────────────────────────────────────────

/// Encode base to 2 bits (A=0, C=1, G=2, T=3). Returns None for non-ACGT.
#[inline(always)]
fn enc(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// 2-bit encode exactly `k` bases from `seq[start..]`.
/// Returns None if any base is non-ACGT or the slice is too short.
#[inline]
pub(crate) fn encode_slice(seq: &[u8], start: usize, k: usize) -> Option<u64> {
    if start + k > seq.len() {
        return None;
    }
    let mut val: u64 = 0;
    for &b in &seq[start..start + k] {
        val = (val << 2) | enc(b)?;
    }
    Some(val)
}

/// Decode a 2-bit encoded k-mer back to a DNA string.
#[cfg(test)]
fn decode_kmer(kmer: u64, k: usize) -> String {
    let mut s = String::with_capacity(k);
    for i in (0..k).rev() {
        s.push(match (kmer >> (i * 2)) & 3 {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => unreachable!(),
        });
    }
    s
}

/// Count mismatching base positions between two 2-bit encoded k-mers.
///
/// Each 2-bit group encodes one base. A mismatch contributes ≥1 bit to XOR.
/// This exploits the 2-bit structure: `(xor | (xor >> 1)) & MASK` has exactly
/// one 1-bit per mismatching position → `count_ones()` = Hamming distance.
/// Compiles to ~4 instructions on ARM64: EOR, ORR-shift, AND, CNT.
#[inline(always)]
pub(crate) fn hamming_2bit(a: u64, b: u64, k: u32) -> u32 {
    let xor = a ^ b;
    let diff = (xor | (xor >> 1)) & 0x5555_5555_5555_5555u64;
    let used_bits = k * 2;
    let mask: u64 = if used_bits >= 64 {
        u64::MAX
    } else {
        (1u64 << used_bits) - 1
    };
    (diff & mask).count_ones()
}

/// Return the list of 0-based mismatch positions (guide-relative, position 0 = first base).
fn mismatch_positions(guide: u64, refseq: u64, k: usize) -> Vec<u8> {
    let xor = guide ^ refseq;
    let mut positions = Vec::new();
    for i in 0..k {
        let shift = (k - 1 - i) * 2;
        if (xor >> shift) & 3 != 0 {
            positions.push(i as u8);
        }
    }
    positions
}

/// Reverse complement of a 2-bit encoded k-mer packed at the LSBs.
///
/// Complement flips all bits (A=00↔T=11, C=01↔G=10). Reverse swaps the order
/// of the k 2-bit groups. The reversal is done by iterating pairs — the
/// compiler unrolls short fixed-length loops well.
pub(crate) fn reverse_complement_2bit(x: u64, k: usize) -> u64 {
    // Complement all bits within the used region
    let bits = k * 2;
    let mask: u64 = if bits >= 64 { u64::MAX } else { (1u64 << bits) - 1 };
    let comp = (!x) & mask;
    // Reverse the order of the k 2-bit groups.
    // Read from LSB (pair 0) upward and build new value MSB-first so that
    // the last base ends up at the highest bit position.
    let mut out: u64 = 0;
    for i in 0..k {
        let pair = (comp >> (i * 2)) & 3;
        out = (out << 2) | pair;
    }
    out
}

// ── Guide indexing ──────────────────────────────────────────────────────────

/// A single indexed guide (forward or RC strand).
#[derive(Clone)]
struct Guide {
    /// 2-bit encoded spacer sequence (GUIDE_LEN = 20 bases, packed at LSBs).
    encoded: u64,
    /// Original forward-strand spacer sequence (always).
    seq: String,
    /// User-visible guide number (0-based index into the input list).
    idx: u32,
    /// '+' for forward, '−' for reverse complement.
    strand: char,
}

/// Packed hash key: upper 32 bits = 5-mer value, lower 8 bits = seed offset.
type SeedKey = u64;

#[inline(always)]
fn make_key(seed_val: u64, seed_offset: u8) -> SeedKey {
    (seed_val << 8) | seed_offset as u64
}

/// Build the guide Vec and seed hash from a list of raw 20-nt sequences.
///
/// Returns `(guides, hash)` where hash maps seed key → indices into `guides`.
fn build_guide_index(raw_guides: &[String]) -> (Vec<Guide>, HashMap<SeedKey, Vec<u32>>) {
    let mut guides: Vec<Guide> = Vec::with_capacity(raw_guides.len() * 2);

    for (idx, seq) in raw_guides.iter().enumerate() {
        let bytes = seq.as_bytes();
        if let Some(fwd_enc) = encode_slice(bytes, 0, GUIDE_LEN) {
            let rc_enc = reverse_complement_2bit(fwd_enc, GUIDE_LEN);
            guides.push(Guide {
                encoded: fwd_enc,
                seq: seq.clone(),
                idx: idx as u32,
                strand: '+',
            });
            guides.push(Guide {
                encoded: rc_enc,
                seq: seq.clone(),
                idx: idx as u32,
                strand: '-',
            });
        }
    }

    let seed_mask: u64 = (1u64 << (SEED_LEN * 2)) - 1;
    let mut hash: HashMap<SeedKey, Vec<u32>> = HashMap::new();

    for (gi, guide) in guides.iter().enumerate() {
        for s in 0..NUM_SEEDS {
            let offset = s * SEED_LEN;
            // Extract seed at guide-relative offset `offset`, packed at LSBs.
            // In LSB-packed representation: base i is at bit offset (k-1-i)*2.
            let shift = (GUIDE_LEN - offset - SEED_LEN) * 2;
            let seed_val = (guide.encoded >> shift) & seed_mask;
            let key = make_key(seed_val, offset as u8);
            hash.entry(key).or_default().push(gi as u32);
        }
    }

    (guides, hash)
}

// ── One-pass parallel off-target scan ──────────────────────────────────────

/// One detected off-target hit.
pub struct OtHit {
    pub guide_idx: u32,
    pub guide_seq: String,
    pub chrom: String,
    pub pos: u64,
    pub strand: char,
    pub mismatches: u32,
    pub mm_positions: Vec<u8>,
}

/// Scan one sub-chunk of raw genome sequence for off-target hits.
///
/// `seq` must contain all bases needed for guide starts in `[0, seq.len() - GUIDE_LEN]`.
/// `base_offset` is the absolute genomic coordinate of `seq[0]`.
fn scan_subchunk(
    seq: &[u8],
    base_offset: u64,
    chrom: &str,
    guides: &[Guide],
    seed_hash: &HashMap<SeedKey, Vec<u32>>,
    max_mm: u8,
) -> Vec<OtHit> {
    let n = seq.len();
    if n < GUIDE_LEN {
        return Vec::new();
    }

    let max_guide_starts = n - GUIDE_LEN + 1;

    // Deduplication: a guide can be triggered by multiple seed offsets at the
    // same position. Only report each (guide_vec_idx, abs_start) pair once.
    // Off-target counts are typically tiny so this set stays small.
    let mut seen: HashMap<(u32, u64), ()> = HashMap::new();
    let mut hits: Vec<OtHit> = Vec::new();

    // For each genome position `p`, extract the 5-mer and check all 4 seed offsets.
    // (p, soff) → potential guide start at p − soff.
    for p in 0..=(n.saturating_sub(SEED_LEN)) {
        if let Some(seed_val) = encode_slice(seq, p, SEED_LEN) {
            for &soff in &[0u8, 5, 10, 15] {
                let key = make_key(seed_val, soff);
                if let Some(candidates) = seed_hash.get(&key) {
                    // Guide starts at genome position p − soff
                    let guide_start = match p.checked_sub(soff as usize) {
                        Some(gs) if gs + GUIDE_LEN <= n => gs,
                        _ => continue,
                    };
                    // Only process guide starts within this sub-chunk's responsibility.
                    if guide_start >= max_guide_starts {
                        continue;
                    }
                    if let Some(ref_enc) = encode_slice(seq, guide_start, GUIDE_LEN) {
                        for &gi in candidates {
                            let mm = hamming_2bit(guides[gi as usize].encoded, ref_enc, GUIDE_LEN as u32);
                            if mm <= max_mm as u32 {
                                let abs_start = base_offset + guide_start as u64;
                                let key2 = (gi, abs_start);
                                if seen.insert(key2, ()).is_none() {
                                    let mm_pos = mismatch_positions(
                                        guides[gi as usize].encoded,
                                        ref_enc,
                                        GUIDE_LEN,
                                    );
                                    hits.push(OtHit {
                                        guide_idx: guides[gi as usize].idx,
                                        guide_seq: guides[gi as usize].seq.clone(),
                                        chrom: chrom.to_string(),
                                        pos: abs_start,
                                        strand: guides[gi as usize].strand,
                                        mismatches: mm,
                                        mm_positions: mm_pos,
                                    });
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    hits
}

/// Scan a full chromosome sequence in parallel sub-chunks.
///
/// The chromosome sequence is divided into blocks of `chunk_size` guide starts.
/// Each block is processed independently by Rayon, then hits are merged.
fn scan_chromosome(
    seq: &[u8],
    chrom: &str,
    guides: &[Guide],
    seed_hash: &HashMap<SeedKey, Vec<u32>>,
    max_mm: u8,
    chunk_size: usize,
) -> Vec<OtHit> {
    let n = seq.len();
    if n < GUIDE_LEN {
        return Vec::new();
    }

    let total_starts = n - GUIDE_LEN + 1;
    let num_chunks = (total_starts + chunk_size - 1) / chunk_size;

    // Build per-chunk (start, end) ranges of guide start positions.
    let ranges: Vec<(usize, usize)> = (0..num_chunks)
        .map(|i| {
            let gs_start = i * chunk_size;
            let gs_end = ((i + 1) * chunk_size).min(total_starts);
            // The sub-chunk slice spans [gs_start, gs_end + GUIDE_LEN - 1]
            let slice_end = (gs_end + GUIDE_LEN - 1).min(n);
            (gs_start, slice_end)
        })
        .collect();

    ranges
        .into_par_iter()
        .flat_map(|(slice_start, slice_end)| {
            let subseq = &seq[slice_start..slice_end];
            scan_subchunk(
                subseq,
                slice_start as u64,
                chrom,
                guides,
                seed_hash,
                max_mm,
            )
        })
        .collect()
}

// ── FASTA streaming ─────────────────────────────────────────────────────────

/// Iterate over chromosomes in a FASTA, yielding (name, sequence).
///
/// Sequences are uppercased; non-DNA bytes (N, other IUPAC) are kept as-is
/// since `encode_slice` will return None for them, naturally skipping any
/// window that contains them.
fn iter_fasta_chroms(
    path: &Path,
) -> anyhow::Result<Vec<(String, Vec<u8>)>> {
    let file = std::fs::File::open(path)?;
    let reader: Box<dyn BufRead> = if path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e == "gz")
        .unwrap_or(false)
    {
        Box::new(BufReader::new(flate2::read::GzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(4 * 1_048_576, file))
    };

    let mut chroms: Vec<(String, Vec<u8>)> = Vec::new();
    let mut current_name = String::new();
    let mut current_seq: Vec<u8> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim_end();
        if trimmed.starts_with('>') {
            if !current_name.is_empty() {
                chroms.push((current_name.clone(), current_seq.clone()));
                current_seq.clear();
            }
            // Take only the first word of the header as the chromosome name.
            current_name = trimmed[1..]
                .split_ascii_whitespace()
                .next()
                .unwrap_or("unknown")
                .to_string();
        } else if !trimmed.is_empty() {
            current_seq.extend(trimmed.bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    if !current_name.is_empty() {
        chroms.push((current_name, current_seq));
    }

    Ok(chroms)
}

// ── Guide loading ───────────────────────────────────────────────────────────

/// Load guide spacer sequences (exactly 20 nt) from a CSV or plain-text file.
///
/// CSV: looks for a column named `sequence_23mer`, `sequence_20mer`,
/// `guide_sequence`, or `sequence`. Takes the first 20 bases.
/// Plain text: one sequence per non-empty, non-comment line.
pub fn load_guides(path: &Path) -> anyhow::Result<Vec<String>> {
    let raw = std::fs::read_to_string(path)?;
    // Detect format: if it looks like a CSV (has commas and a header matching
    // known column names), parse as CSV. Otherwise treat as plain text.
    if raw.lines().next().map(|l| l.contains(',')).unwrap_or(false) {
        load_guides_csv(path)
    } else {
        load_guides_text(&raw)
    }
}

fn load_guides_csv(path: &Path) -> anyhow::Result<Vec<String>> {
    let mut rdr = csv::ReaderBuilder::new()
        .flexible(true)
        .from_path(path)?;
    let headers = rdr.headers()?.clone();
    let col = headers
        .iter()
        .position(|h| {
            matches!(
                h,
                "sequence_23mer"
                    | "sequence_20mer"
                    | "guide_sequence"
                    | "sequence"
                    | "spacer"
            )
        })
        .ok_or_else(|| {
            anyhow::anyhow!(
                "CSV must have a column named 'sequence_23mer', 'sequence_20mer', \
                 'guide_sequence', 'spacer', or 'sequence'"
            )
        })?;

    let mut guides = Vec::new();
    for result in rdr.records() {
        let rec = result?;
        if let Some(seq) = rec.get(col) {
            let seq = seq.trim();
            if seq.len() >= GUIDE_LEN {
                guides.push(seq[..GUIDE_LEN].to_ascii_uppercase());
            }
        }
    }
    Ok(guides)
}

fn load_guides_text(raw: &str) -> anyhow::Result<Vec<String>> {
    let guides: Vec<String> = raw
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .filter(|l| l.len() >= GUIDE_LEN)
        .map(|l| l[..GUIDE_LEN].to_ascii_uppercase())
        .collect();
    Ok(guides)
}

// ── Public API ──────────────────────────────────────────────────────────────

/// Configuration for the off-target scan.
pub struct OffTargetConfig {
    /// Maximum allowed mismatches (1–3, default 3).
    pub max_mismatches: u8,
    /// Parallel sub-chunk size in guide-start positions (default 4 M).
    pub chunk_size: usize,
    /// Only keep the first N guides from the input (0 = all).
    pub max_guides: usize,
}

impl Default for OffTargetConfig {
    fn default() -> Self {
        Self {
            max_mismatches: 3,
            chunk_size: 4_000_000,
            max_guides: 0,
        }
    }
}

/// Run the full off-target scan.
///
/// * `guides_path` — CSV or plain-text file containing 20-nt spacer sequences.
/// * `reference_path` — Reference genome FASTA (plain or `.gz`).
/// * `output_path` — Output CSV path.
/// * `cfg` — Scan configuration.
pub fn run_offtarget_scan(
    guides_path: &Path,
    reference_path: &Path,
    output_path: &Path,
    cfg: &OffTargetConfig,
) -> anyhow::Result<()> {
    let t0 = Instant::now();

    // 1. Load guide sequences.
    let mut raw_guides = load_guides(guides_path)?;
    if cfg.max_guides > 0 && raw_guides.len() > cfg.max_guides {
        eprintln!(
            "[offtarget] Capping input from {} to {} guides",
            raw_guides.len(),
            cfg.max_guides
        );
        raw_guides.truncate(cfg.max_guides);
    }
    // Deduplicate while preserving order.
    let mut seen_seqs = std::collections::HashSet::new();
    raw_guides.retain(|s| seen_seqs.insert(s.clone()));
    eprintln!(
        "[offtarget] Loaded {} unique guides, max_mismatches={}",
        raw_guides.len(),
        cfg.max_mismatches
    );

    // 2. Build guide index (forward + RC).
    let (guides, seed_hash) = build_guide_index(&raw_guides);
    eprintln!(
        "[offtarget] Guide index: {} entries ({} fwd + RC), {} seed buckets",
        guides.len(),
        raw_guides.len(),
        seed_hash.len(),
    );

    // 3. Create output CSV.
    if let Some(parent) = output_path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let out_file = std::fs::File::create(output_path)?;
    let mut wtr = csv::Writer::from_writer(std::io::BufWriter::new(out_file));
    wtr.write_record(&[
        "guide_idx",
        "guide_seq",
        "chrom",
        "pos",
        "strand",
        "mismatches",
        "mm_positions",
    ])?;

    // 4. Stream reference genome chromosome by chromosome.
    eprintln!(
        "[offtarget] Scanning reference genome: {}",
        reference_path.display()
    );
    let chroms = iter_fasta_chroms(reference_path)?;
    eprintln!("[offtarget] {} chromosomes/contigs found", chroms.len());

    let mut total_hits: u64 = 0;

    for (chrom_name, chrom_seq) in &chroms {
        let t_chrom = Instant::now();
        let hits = scan_chromosome(
            chrom_seq,
            chrom_name,
            &guides,
            &seed_hash,
            cfg.max_mismatches,
            cfg.chunk_size,
        );
        let n_hits = hits.len();
        total_hits += n_hits as u64;

        for hit in &hits {
            wtr.write_record(&[
                hit.guide_idx.to_string(),
                hit.guide_seq.clone(),
                hit.chrom.clone(),
                hit.pos.to_string(),
                hit.strand.to_string(),
                hit.mismatches.to_string(),
                hit.mm_positions
                    .iter()
                    .map(|p| p.to_string())
                    .collect::<Vec<_>>()
                    .join(";"),
            ])?;
        }

        eprintln!(
            "[offtarget]   {} — {:.1} Mb, {} hits in {:.2}s",
            chrom_name,
            chrom_seq.len() as f64 / 1_048_576.0,
            n_hits,
            t_chrom.elapsed().as_secs_f64(),
        );
    }

    wtr.flush()?;

    eprintln!(
        "[offtarget] Done. {} total hits across {} chroms in {:.1}s → {}",
        total_hits,
        chroms.len(),
        t0.elapsed().as_secs_f64(),
        output_path.display(),
    );

    Ok(())
}

// ── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encode_decode_round_trip() {
        let seq = b"ACGTACGTACGTACGTACGT"; // 20 bases
        let enc = encode_slice(seq, 0, 20).unwrap();
        assert_eq!(decode_kmer(enc, 20), "ACGTACGTACGTACGTACGT");
    }

    #[test]
    fn hamming_identical_is_zero() {
        let seq = b"AAAACCCCGGGGTTTTAAAA";
        let a = encode_slice(seq, 0, 20).unwrap();
        assert_eq!(hamming_2bit(a, a, 20), 0);
    }

    #[test]
    fn hamming_one_mismatch() {
        let a = encode_slice(b"AAAACCCCGGGGTTTTAAAA", 0, 20).unwrap();
        // Change first base A→C
        let b = encode_slice(b"CAAACCCCGGGGTTTTAAAA", 0, 20).unwrap();
        assert_eq!(hamming_2bit(a, b, 20), 1);
    }

    #[test]
    fn hamming_all_mismatch() {
        let a = encode_slice(b"AAAAAAAAAAAAAAAAAAAA", 0, 20).unwrap(); // all A
        let b = encode_slice(b"TTTTTTTTTTTTTTTTTTTT", 0, 20).unwrap(); // all T (complement)
        assert_eq!(hamming_2bit(a, b, 20), 20);
    }

    #[test]
    fn rc_complement_and_reverse() {
        // RC of AAAA (k=4) should be TTTT
        let a = encode_slice(b"AAAA", 0, 4).unwrap();
        let rc = reverse_complement_2bit(a, 4);
        assert_eq!(decode_kmer(rc, 4), "TTTT");

        // RC of ACGT (k=4) should be ACGT (palindrome)
        let acgt = encode_slice(b"ACGT", 0, 4).unwrap();
        let rc2 = reverse_complement_2bit(acgt, 4);
        assert_eq!(decode_kmer(rc2, 4), "ACGT");
    }

    #[test]
    fn exact_match_found() {
        // Build a tiny guide and scan a genome that contains it exactly.
        let guide = "AAAACCCCGGGGTTTTAAAA"; // 20 nt
        let genome = format!("NNNNN{}NNNNN", guide);
        let raw_guides = vec![guide.to_string()];
        let (guides, hash) = build_guide_index(&raw_guides);
        let hits = scan_subchunk(genome.as_bytes(), 0, "test", &guides, &hash, 0);
        // Exactly one hit (the forward strand), at offset 5.
        let fwd_hits: Vec<_> = hits.iter().filter(|h| h.strand == '+').collect();
        assert_eq!(fwd_hits.len(), 1);
        assert_eq!(fwd_hits[0].pos, 5);
        assert_eq!(fwd_hits[0].mismatches, 0);
    }

    #[test]
    fn one_mismatch_found() {
        let guide = "AAAACCCCGGGGTTTTAAAA";
        // Introduce one mismatch at position 0 (A→T)
        let target = "TAAACCCCGGGGTTTTAAAA";
        let genome = format!("NNNNN{}NNNNN", target);
        let raw_guides = vec![guide.to_string()];
        let (guides, hash) = build_guide_index(&raw_guides);
        let hits = scan_subchunk(genome.as_bytes(), 0, "test", &guides, &hash, 1);
        let fwd_hits: Vec<_> = hits.iter().filter(|h| h.strand == '+').collect();
        assert_eq!(fwd_hits.len(), 1);
        assert_eq!(fwd_hits[0].mismatches, 1);
        assert_eq!(fwd_hits[0].mm_positions, vec![0u8]);
    }

    #[test]
    fn no_hit_beyond_threshold() {
        let guide = "AAAACCCCGGGGTTTTAAAA";
        // 4 mismatches — should not be found with max_mm = 3
        let target = "TTTTCCCCGGGGTTTTAAAA";
        let genome = format!("NNNNN{}NNNNN", target);
        let raw_guides = vec![guide.to_string()];
        let (guides, hash) = build_guide_index(&raw_guides);
        let hits = scan_subchunk(genome.as_bytes(), 0, "test", &guides, &hash, 3);
        let fwd_hits: Vec<_> = hits.iter().filter(|h| h.strand == '+').collect();
        assert_eq!(fwd_hits.len(), 0);
    }

    #[test]
    fn rc_strand_detected() {
        let guide = "AAAACCCCGGGGTTTTAAAA";
        // RC of guide embedded in genome (RC of each base, reversed)
        let rc = reverse_complement_2bit(
            encode_slice(guide.as_bytes(), 0, 20).unwrap(),
            20,
        );
        let rc_str = decode_kmer(rc, 20);
        let genome = format!("NNNNN{}NNNNN", rc_str);
        let raw_guides = vec![guide.to_string()];
        let (guides, hash) = build_guide_index(&raw_guides);
        let hits = scan_subchunk(genome.as_bytes(), 0, "test", &guides, &hash, 0);
        let rc_hits: Vec<_> = hits.iter().filter(|h| h.strand == '-').collect();
        assert_eq!(rc_hits.len(), 1);
        assert_eq!(rc_hits[0].mismatches, 0);
    }

    #[test]
    fn no_duplicate_hits() {
        // A guide with a palindrome-like sequence hit by multiple seeds should
        // appear exactly once in the output.
        let guide = "AAAAACGGGCAAAAACGGGC"; // seeds at 0,5,10,15 differ
        let genome = format!("NNNNN{}NNNNN", guide);
        let raw_guides = vec![guide.to_string()];
        let (guides, hash) = build_guide_index(&raw_guides);
        let hits = scan_subchunk(genome.as_bytes(), 0, "test", &guides, &hash, 0);
        let fwd_hits: Vec<_> = hits.iter().filter(|h| h.strand == '+').collect();
        assert_eq!(fwd_hits.len(), 1, "expected exactly 1 forward hit, got {}", fwd_hits.len());
    }
}
