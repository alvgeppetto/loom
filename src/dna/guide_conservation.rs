//! Guide conservation scanner.
//!
//! Streams a large FASTA file record-by-record, checks each genome for
//! the presence of specific CRISPR guide sequences (and their reverse
//! complements) using an Aho-Corasick automaton, and reports per-guide
//! conservation rates across the full pangenome.

use aho_corasick::AhoCorasick;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::time::Instant;

/// One guide loaded from the TSV.
#[derive(Debug, Clone)]
struct Guide {
    id: String,
    sequence: String,
    revcomp: String,
}

/// Result for a single guide.
#[derive(Debug, serde::Serialize)]
pub struct GuideResult {
    pub guide_id: String,
    pub sequence: String,
    pub hits: u64,
    pub total_genomes: u64,
    pub conservation_pct: f64,
}

/// Reverse complement a DNA string.
fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'C' | 'c' => 'G',
            'G' | 'g' => 'C',
            _ => 'N',
        })
        .collect()
}

/// Load guides from a TSV file (header: guide_id, sequence, ...).
fn load_guides(path: &Path) -> anyhow::Result<Vec<Guide>> {
    let file = std::fs::File::open(path)?;
    let reader = BufReader::new(file);
    let mut guides = Vec::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if i == 0 {
            continue; // skip header
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 2 {
            continue;
        }
        let id = fields[0].to_string();
        let seq = fields[1].to_uppercase();
        let rc = reverse_complement(&seq);
        guides.push(Guide {
            id,
            sequence: seq,
            revcomp: rc,
        });
    }
    Ok(guides)
}

/// Run the guide conservation scan.
///
/// Reads guides from `guides_tsv`, streams FASTA from `fasta_path`,
/// and writes JSON results to `output`.
pub fn run(fasta_path: &Path, guides_tsv: &Path, output: &Path) -> anyhow::Result<()> {
    let t0 = Instant::now();

    // Load guides
    let guides = load_guides(guides_tsv)?;
    anyhow::ensure!(!guides.is_empty(), "No guides loaded from TSV");
    eprintln!(
        "[guide-conservation] Loaded {} guides from {}",
        guides.len(),
        guides_tsv.display()
    );

    // Build pattern list: for each guide, fwd + revcomp.
    // Pattern index i*2 = fwd of guide i, i*2+1 = revcomp of guide i.
    let mut patterns: Vec<String> = Vec::with_capacity(guides.len() * 2);
    for g in &guides {
        patterns.push(g.sequence.clone());
        patterns.push(g.revcomp.clone());
    }
    let ac = AhoCorasick::new(&patterns)?;
    eprintln!(
        "[guide-conservation] Built Aho-Corasick automaton for {} patterns ({} guides × 2)",
        patterns.len(),
        guides.len()
    );

    // Per-guide hit counters
    let num_guides = guides.len();
    let mut hits = vec![0u64; num_guides];
    let mut total_genomes: u64 = 0;

    // Stream FASTA record by record
    let file = std::fs::File::open(fasta_path)?;
    let reader = BufReader::with_capacity(4 << 20, file); // 4 MB buffer
    let mut current_seq = Vec::<u8>::with_capacity(32_000); // typical ~30 KB genome
    let mut in_record = false;

    eprintln!(
        "[guide-conservation] Scanning {}...",
        fasta_path.display()
    );

    let report_interval = 1_000_000u64;

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim_ascii();

        if trimmed.starts_with(">") {
            // Process previous record
            if in_record && !current_seq.is_empty() {
                scan_genome(&ac, &current_seq, num_guides, &mut hits);
                total_genomes += 1;
                if total_genomes % report_interval == 0 {
                    let elapsed = t0.elapsed().as_secs_f64();
                    let rate = total_genomes as f64 / elapsed;
                    eprintln!(
                        "[guide-conservation] {:>10} genomes scanned  ({:.0} genomes/s, {:.1} min elapsed)",
                        total_genomes, rate, elapsed / 60.0
                    );
                }
                current_seq.clear();
            }
            in_record = true;
        } else if in_record {
            // Accumulate sequence lines (uppercase)
            for b in trimmed.bytes() {
                current_seq.push(b.to_ascii_uppercase());
            }
        }
    }
    // Process last record
    if in_record && !current_seq.is_empty() {
        scan_genome(&ac, &current_seq, num_guides, &mut hits);
        total_genomes += 1;
    }

    let elapsed = t0.elapsed().as_secs_f64();
    eprintln!(
        "[guide-conservation] Done: {} genomes in {:.1}s ({:.0} genomes/s)",
        total_genomes,
        elapsed,
        total_genomes as f64 / elapsed
    );

    // Build results
    let results: Vec<GuideResult> = guides
        .iter()
        .enumerate()
        .map(|(i, g)| {
            let pct = if total_genomes > 0 {
                hits[i] as f64 / total_genomes as f64 * 100.0
            } else {
                0.0
            };
            GuideResult {
                guide_id: g.id.clone(),
                sequence: g.sequence.clone(),
                hits: hits[i],
                total_genomes,
                conservation_pct: (pct * 100.0).round() / 100.0, // 2 decimal places
            }
        })
        .collect();

    // Print summary table
    eprintln!("\n{:<25} {:>12} {:>12} {:>10}", "Guide", "Hits", "Total", "Rate");
    eprintln!("{}", "-".repeat(62));
    for r in &results {
        eprintln!(
            "{:<25} {:>12} {:>12} {:>9.2}%",
            r.guide_id, r.hits, r.total_genomes, r.conservation_pct
        );
    }

    // Write JSON
    if let Some(parent) = output.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut f = std::fs::File::create(output)?;
    let json = serde_json::to_string_pretty(&results)?;
    f.write_all(json.as_bytes())?;
    f.write_all(b"\n")?;
    eprintln!(
        "\n[guide-conservation] Results written to {}",
        output.display()
    );

    Ok(())
}

/// Check one genome against the Aho-Corasick automaton.
/// Sets hit flags for each guide that matches (fwd or revcomp).
#[inline]
fn scan_genome(ac: &AhoCorasick, seq: &[u8], num_guides: usize, hits: &mut [u64]) {
    // Track which guides hit this genome (avoid double-counting fwd + revcomp)
    let mut seen = vec![false; num_guides];

    for mat in ac.find_overlapping_iter(seq) {
        let guide_idx = mat.pattern().as_usize() / 2; // fwd and revcomp map to same guide
        if !seen[guide_idx] {
            seen[guide_idx] = true;
            hits[guide_idx] += 1;
        }
        // Early exit if all guides found in this genome
        if seen.iter().all(|&s| s) {
            break;
        }
    }
}
