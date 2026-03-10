//! Sequence-level novelty verification for CRISPR guide candidates.
//!
//! Compares novel guide sequences against a comprehensive database of published
//! CRISPR guide sequences using 2-bit Hamming distance (reusing the same ARM64-
//! optimised core from `offtarget.rs`).
//!
//! For each novel target, we check both strands of every published guide and
//! report the closest match. A target is "sequence-novel" if no published guide
//! matches within `max_mismatches` (default 3).

use crate::dna::offtarget::{encode_slice, hamming_2bit, reverse_complement_2bit};
use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::path::Path;
use std::time::Instant;

// ── Input schemas ───────────────────────────────────────────────────────────

/// One novel target from `novel_targets_52.json`.
#[derive(Debug, Deserialize)]
struct NovelTargetsFile {
    regions: Vec<NovelRegion>,
}

#[derive(Debug, Deserialize)]
struct NovelRegion {
    best_seq: String,
    best_pos: u64,
    gene: String,
    #[serde(default)]
    novelty_score: f64,
}

/// The comprehensive published-guide database from the Python curation script.
#[derive(Debug, Deserialize)]
struct PublishedDb {
    total_guides: usize,
    guides: Vec<PublishedGuide>,
}

#[derive(Debug, Deserialize)]
struct PublishedGuide {
    id: String,
    sequence: String,
    #[serde(default)]
    paper: String,
    #[serde(default)]
    pmid: String,
    #[serde(default)]
    gene: String,
    #[serde(default)]
    wuhan_hu1_pos: Option<i64>,
    length: usize,
}

// ── Output schema ───────────────────────────────────────────────────────────

#[derive(Debug, Serialize)]
struct NoveltyReport {
    novel_count: usize,
    total_targets: usize,
    published_guides_checked: usize,
    max_mismatches: u8,
    targets: Vec<TargetVerdict>,
}

#[derive(Debug, Serialize)]
struct TargetVerdict {
    sequence: String,
    position: u64,
    gene: String,
    is_novel: bool,
    closest_distance: u32,
    closest_published_id: String,
    closest_published_seq: String,
    closest_published_paper: String,
}

// ── Core comparison ─────────────────────────────────────────────────────────

/// Compare two sequences of potentially different lengths.
/// If lengths differ, compare only the overlapping prefix (shorter length)
/// and add the length difference as extra mismatches.
fn compare_guides(novel_enc: u64, novel_len: usize, pub_enc: u64, pub_rc: u64, pub_len: usize) -> u32 {
    let min_len = novel_len.min(pub_len);
    let extra = (novel_len.max(pub_len) - min_len) as u32;

    // Truncate both to the shorter length for comparison
    let novel_trunc = if novel_len > min_len {
        novel_enc >> ((novel_len - min_len) * 2)
    } else {
        novel_enc
    };

    let pub_fwd_trunc = if pub_len > min_len {
        pub_enc >> ((pub_len - min_len) * 2)
    } else {
        pub_enc
    };

    let pub_rc_trunc = if pub_len > min_len {
        pub_rc >> ((pub_len - min_len) * 2)
    } else {
        pub_rc
    };

    let d_fwd = hamming_2bit(novel_trunc, pub_fwd_trunc, min_len as u32) + extra;
    let d_rc = hamming_2bit(novel_trunc, pub_rc_trunc, min_len as u32) + extra;

    d_fwd.min(d_rc)
}

// ── Public entry point ──────────────────────────────────────────────────────

pub fn run(
    novel_path: &Path,
    published_path: &Path,
    output_path: &Path,
    max_mismatches: u8,
) -> Result<()> {
    let start = Instant::now();

    // Load novel targets
    let novel_text = std::fs::read_to_string(novel_path)
        .with_context(|| format!("reading {}", novel_path.display()))?;
    let novel_file: NovelTargetsFile = serde_json::from_str(&novel_text)
        .with_context(|| "parsing novel targets JSON")?;
    println!("Loaded {} novel targets", novel_file.regions.len());

    // Load published guide database
    let pub_text = std::fs::read_to_string(published_path)
        .with_context(|| format!("reading {}", published_path.display()))?;
    let pub_db: PublishedDb = serde_json::from_str(&pub_text)
        .with_context(|| "parsing published guides JSON")?;
    println!(
        "Loaded {} published guides from database",
        pub_db.guides.len()
    );

    // 2-bit encode all published guides (forward + RC)
    struct EncodedPub {
        fwd: u64,
        rc: u64,
        len: usize,
        idx: usize,
    }

    let encoded_pubs: Vec<EncodedPub> = pub_db
        .guides
        .iter()
        .enumerate()
        .filter_map(|(i, g)| {
            let seq = g.sequence.as_bytes();
            let len = seq.len();
            let fwd = encode_slice(seq, 0, len)?;
            let rc = reverse_complement_2bit(fwd, len);
            Some(EncodedPub { fwd, rc, len, idx: i })
        })
        .collect();

    println!(
        "Encoded {} published guides (2-bit), {} skipped (ambiguous bases)",
        encoded_pubs.len(),
        pub_db.guides.len() - encoded_pubs.len()
    );

    // Compare each novel target against all published guides
    let mut verdicts = Vec::with_capacity(novel_file.regions.len());
    let mut novel_count = 0usize;
    let mut total_comparisons = 0u64;

    for region in &novel_file.regions {
        let seq = region.best_seq.to_uppercase().replace('U', "T");
        let novel_len = seq.len();
        let novel_enc = match encode_slice(seq.as_bytes(), 0, novel_len) {
            Some(e) => e,
            None => {
                eprintln!(
                    "  WARN: cannot encode novel target at pos {} (non-ACGT base), skipping",
                    region.best_pos
                );
                continue;
            }
        };

        let mut best_dist = u32::MAX;
        let mut best_idx = 0usize;

        for ep in &encoded_pubs {
            let d = compare_guides(novel_enc, novel_len, ep.fwd, ep.rc, ep.len);
            total_comparisons += 1;
            if d < best_dist {
                best_dist = d;
                best_idx = ep.idx;
                if d == 0 {
                    break; // Exact match, can't do better
                }
            }
        }

        let is_novel = best_dist > max_mismatches as u32;
        if is_novel {
            novel_count += 1;
        }

        let closest = &pub_db.guides[best_idx];
        let status = if is_novel { "NOVEL" } else { "MATCH" };
        println!(
            "  {} pos={:<6} gene={:<10} closest_dist={} closest=\"{}\" ({})",
            status, region.best_pos, region.gene, best_dist, closest.id, closest.paper
        );

        verdicts.push(TargetVerdict {
            sequence: seq,
            position: region.best_pos,
            gene: region.gene.clone(),
            is_novel,
            closest_distance: best_dist,
            closest_published_id: closest.id.clone(),
            closest_published_seq: closest.sequence.clone(),
            closest_published_paper: closest.paper.clone(),
        });
    }

    // Summary
    println!("\n── Summary ──────────────────────────────────────────");
    println!("Targets checked:       {}", verdicts.len());
    println!("Published guides:      {}", encoded_pubs.len());
    println!("Total comparisons:     {}", total_comparisons);
    println!("Max mismatches:        {}", max_mismatches);
    println!("Sequence-novel:        {}/{}", novel_count, verdicts.len());
    println!("Elapsed:               {:.3}s", start.elapsed().as_secs_f64());

    // Write report
    let report = NoveltyReport {
        novel_count,
        total_targets: verdicts.len(),
        published_guides_checked: encoded_pubs.len(),
        max_mismatches,
        targets: verdicts,
    };

    if let Some(parent) = output_path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let json = serde_json::to_string_pretty(&report)?;
    std::fs::write(output_path, json.as_bytes())
        .with_context(|| format!("writing {}", output_path.display()))?;
    println!("\nReport written to {}", output_path.display());

    Ok(())
}
