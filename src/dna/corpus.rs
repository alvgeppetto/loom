use anyhow::{Context, bail};
use serde::Serialize;
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use crate::dna::normalize_dna_sequence;

/// Maximum sequence length accepted per record (guard against runaway allocations).
const MAX_SEQUENCE_BYTES: usize = 256 * 1024 * 1024; // 256 MiB

#[derive(Debug, Clone, Serialize)]
pub struct DnaCorpusRecord {
    pub id: String,
    pub sequence: String,
    pub source_path: String,
    pub length: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DnaFileFormat {
    Fasta,
    Fastq,
}

/// Compression wrapper detected from filename.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Compression {
    None,
    Gzip,
    Xz,
}

/// Build a normalized DNA corpus JSONL from `input_path` using bounded-memory
/// streaming I/O.  Supports plain, `.gz`, and `.xz` compressed FASTA/FASTQ files.
/// Returns the number of records written.
pub fn build_dna_corpus_jsonl(input_path: &Path, output_path: &Path, min_len: usize) -> anyhow::Result<usize> {
    let files = collect_sequence_files(input_path)?;
    if files.is_empty() {
        bail!("No FASTA/FASTQ files found under {}", input_path.display());
    }

    if let Some(parent) = output_path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)
                .with_context(|| format!("Failed to create output directory {}", parent.display()))?;
        }
    }

    let mut out = fs::File::create(output_path)
        .with_context(|| format!("Failed to create output file {}", output_path.display()))?;

    let mut written = 0usize;
    for file in &files {
        let (format, compression) = detect_format_and_compression(file)?;
        let source = file.to_string_lossy().to_string();

        let result = stream_sequence_records(file, format, compression, min_len, |record| {
            writeln!(out, "{}", serde_json::to_string(&DnaCorpusRecord {
                id: record.id.clone(),
                length: record.sequence.len(),
                sequence: record.sequence.clone(),
                source_path: source.clone(),
            })?)
            .with_context(|| format!("Failed writing to {}", output_path.display()))?;
            written += 1;
            Ok(())
        });

        if let Err(e) = result {
            // Emit a warning for malformed files but continue with remaining files.
            eprintln!("[loom] WARNING: skipping {}: {}", file.display(), e);
        }
    }

    Ok(written)
}

/// Internal streaming state for a single record being accumulated.
struct RecordAccumulator {
    id: String,
    sequence: String,
}

/// Stream records from a single FASTA/FASTQ file, calling `callback` for each
/// valid record that passes `min_len`.  Supports plain, gz, and xz inputs.
fn stream_sequence_records<F>(
    path: &Path,
    format: DnaFileFormat,
    compression: Compression,
    min_len: usize,
    mut callback: F,
) -> anyhow::Result<()>
where
    F: FnMut(&RecordAccumulator) -> anyhow::Result<()>,
{
    let file = fs::File::open(path)
        .with_context(|| format!("Cannot open {}", path.display()))?;

    // Build a type-erased boxed BufReader over the appropriate decompressor.
    let reader: Box<dyn BufRead> = match compression {
        Compression::None => Box::new(BufReader::new(file)),
        Compression::Gzip => {
            use flate2::read::GzDecoder;
            Box::new(BufReader::new(GzDecoder::new(file)))
        }
        Compression::Xz => {
            #[cfg(not(target_arch = "wasm32"))]
            {
                use xz2::read::XzDecoder;
                Box::new(BufReader::new(XzDecoder::new(file)))
            }
            #[cfg(target_arch = "wasm32")]
            {
                return Err(anyhow::anyhow!("xz decompression is not supported in WASM builds"));
            }
        }
    };

    match format {
        DnaFileFormat::Fasta => stream_fasta(reader, min_len, &mut callback),
        DnaFileFormat::Fastq => stream_fastq(reader, min_len, &mut callback),
    }
}

/// Stream-parse FASTA using a line-by-line bounded reader.
fn stream_fasta<R: BufRead, F>(
    reader: R,
    min_len: usize,
    callback: &mut F,
) -> anyhow::Result<()>
where
    F: FnMut(&RecordAccumulator) -> anyhow::Result<()>,
{
    let mut current: Option<RecordAccumulator> = None;

    for (line_no, line_result) in reader.lines().enumerate() {
        let line = line_result
            .with_context(|| format!("I/O error reading FASTA at line {}", line_no + 1))?;
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with(';') {
            continue;
        }

        if let Some(rest) = trimmed.strip_prefix('>') {
            // Flush the previous record.
            if let Some(prev) = current.take() {
                if prev.sequence.len() >= min_len {
                    callback(&prev)?;
                }
            }
            let id = rest
                .split_whitespace()
                .next()
                .unwrap_or("unknown")
                .to_string();
            current = Some(RecordAccumulator { id, sequence: String::new() });
        } else if let Some(ref mut acc) = current {
            let chunk = normalize_dna_sequence(trimmed);
            if acc.sequence.len().saturating_add(chunk.len()) > MAX_SEQUENCE_BYTES {
                bail!(
                    "Sequence '{}' exceeds maximum allowed size of {} bytes",
                    acc.id,
                    MAX_SEQUENCE_BYTES
                );
            }
            acc.sequence.push_str(&chunk);
        }
        // Lines before the first '>' header are silently ignored.
    }

    // Flush the last record.
    if let Some(prev) = current {
        if prev.sequence.len() >= min_len {
            callback(&prev)?;
        }
    }

    Ok(())
}

/// Stream-parse FASTQ using a 4-line-per-record bounded reader.
fn stream_fastq<R: BufRead, F>(
    reader: R,
    min_len: usize,
    callback: &mut F,
) -> anyhow::Result<()>
where
    F: FnMut(&RecordAccumulator) -> anyhow::Result<()>,
{
    let mut lines = reader.lines().enumerate();
    let mut record_idx = 0usize;

    loop {
        // Line 1: header (@...)
        let header_line = match lines.next() {
            None => break,
            Some((n, Err(e))) => bail!("I/O error at FASTQ line {}: {}", n + 1, e),
            Some((_, Ok(l))) => l,
        };
        let header = header_line.trim();
        if header.is_empty() {
            continue;
        }
        if !header.starts_with('@') {
            bail!(
                "Malformed FASTQ: expected '@' at record {} header, got {:?}",
                record_idx + 1,
                &header[..header.len().min(40)]
            );
        }
        let id = header[1..]
            .split_whitespace()
            .next()
            .unwrap_or("unknown")
            .to_string();

        // Line 2: sequence
        let seq_line = match lines.next() {
            None => bail!("Truncated FASTQ: missing sequence line for record '{}'", id),
            Some((n, Err(e))) => bail!("I/O error at FASTQ line {}: {}", n + 1, e),
            Some((_, Ok(l))) => l,
        };
        let raw_seq = seq_line.trim();

        // Line 3: '+' separator
        let plus_line = match lines.next() {
            None => bail!("Truncated FASTQ: missing '+' line for record '{}'", id),
            Some((n, Err(e))) => bail!("I/O error at FASTQ line {}: {}", n + 1, e),
            Some((_, Ok(l))) => l,
        };
        let plus = plus_line.trim();
        if !plus.starts_with('+') {
            bail!(
                "Malformed FASTQ: expected '+' separator at record '{}', got {:?}",
                id,
                &plus[..plus.len().min(40)]
            );
        }

        // Line 4: quality scores
        let qual_line = match lines.next() {
            None => bail!("Truncated FASTQ: missing quality line for record '{}'", id),
            Some((n, Err(e))) => bail!("I/O error at FASTQ line {}: {}", n + 1, e),
            Some((_, Ok(l))) => l,
        };
        let qual = qual_line.trim();

        // Guard: quality length must match sequence length.
        if qual.len() != raw_seq.len() {
            bail!(
                "Malformed FASTQ: sequence length {} != quality length {} for record '{}'",
                raw_seq.len(),
                qual.len(),
                id
            );
        }

        let sequence = normalize_dna_sequence(raw_seq);
        if sequence.len() > MAX_SEQUENCE_BYTES {
            bail!(
                "Sequence '{}' exceeds maximum allowed size of {} bytes",
                id,
                MAX_SEQUENCE_BYTES
            );
        }

        if !sequence.is_empty() && sequence.len() >= min_len {
            callback(&RecordAccumulator { id, sequence })?;
        }

        record_idx += 1;
    }

    Ok(())
}

fn detect_format_and_compression(path: &Path) -> anyhow::Result<(DnaFileFormat, Compression)> {
    let name = path
        .file_name()
        .and_then(|v| v.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();

    // Strip compression suffix first.
    let (base, compression) = if name.ends_with(".gz") {
        (&name[..name.len() - 3], Compression::Gzip)
    } else if name.ends_with(".xz") {
        (&name[..name.len() - 3], Compression::Xz)
    } else {
        (name.as_str(), Compression::None)
    };

    let format = if base.ends_with(".fasta") || base.ends_with(".fa") || base.ends_with(".fna") || base.ends_with(".ffn") {
        DnaFileFormat::Fasta
    } else if base.ends_with(".fastq") || base.ends_with(".fq") {
        DnaFileFormat::Fastq
    } else {
        bail!("Unsupported sequence file extension: {}", path.display())
    };

    Ok((format, compression))
}

fn collect_sequence_files(input_path: &Path) -> anyhow::Result<Vec<PathBuf>> {
    if input_path.is_file() {
        return Ok(vec![input_path.to_path_buf()]);
    }

    if !input_path.is_dir() {
        bail!("Input path does not exist or is not a file/directory: {}", input_path.display());
    }

    let mut files = Vec::new();
    for entry in walkdir::WalkDir::new(input_path).follow_links(true) {
        let entry = entry?;
        if !entry.file_type().is_file() {
            continue;
        }
        let path = entry.path();
        let lowercase = path.to_string_lossy().to_ascii_lowercase();
        if is_sequence_filename(&lowercase) {
            files.push(path.to_path_buf());
        }
    }

    files.sort();
    Ok(files)
}

fn is_sequence_filename(lower: &str) -> bool {
    // Plain formats
    lower.ends_with(".fasta") || lower.ends_with(".fa") || lower.ends_with(".fna")
        || lower.ends_with(".ffn") || lower.ends_with(".fastq") || lower.ends_with(".fq")
    // Gzip-compressed
        || lower.ends_with(".fasta.gz") || lower.ends_with(".fa.gz")
        || lower.ends_with(".fna.gz") || lower.ends_with(".ffn.gz")
        || lower.ends_with(".fastq.gz") || lower.ends_with(".fq.gz")
    // XZ-compressed
        || lower.ends_with(".fasta.xz") || lower.ends_with(".fa.xz")
        || lower.ends_with(".fna.xz") || lower.ends_with(".ffn.xz")
        || lower.ends_with(".fastq.xz") || lower.ends_with(".fq.xz")
}

// ── Legacy helpers kept for unit-test compatibility ──────────────────────────

/// Parse all FASTA records from an in-memory string (used in tests).
pub fn parse_fasta_records(input: &str) -> Vec<(String, String)> {
    let mut records: Vec<(String, String)> = Vec::new();
    let mut current_id: Option<String> = None;
    let mut current_seq = String::new();

    for line in input.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with(';') {
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix('>') {
            if let Some(id) = current_id.take() {
                if !current_seq.is_empty() {
                    records.push((id, current_seq.clone()));
                }
                current_seq.clear();
            }
            let id = rest.split_whitespace().next().unwrap_or("unknown").to_string();
            current_id = Some(id);
            continue;
        }
        current_seq.push_str(&normalize_dna_sequence(trimmed));
    }
    if let Some(id) = current_id {
        if !current_seq.is_empty() {
            records.push((id, current_seq));
        }
    }
    records
}

/// Parse all FASTQ records from an in-memory string (used in tests).
pub fn parse_fastq_records(input: &str) -> Vec<(String, String)> {
    let mut records = Vec::new();
    let lines: Vec<&str> = input.lines().collect();
    let mut index = 0usize;

    while index + 3 < lines.len() {
        let header = lines[index].trim();
        let sequence_line = lines[index + 1].trim();
        let plus = lines[index + 2].trim();
        let _quality = lines[index + 3].trim();
        index += 4;

        if !header.starts_with('@') || !plus.starts_with('+') {
            continue;
        }
        let id = header[1..].split_whitespace().next().unwrap_or("unknown").to_string();
        let sequence = normalize_dna_sequence(sequence_line);
        if sequence.is_empty() {
            continue;
        }
        records.push((id, sequence));
    }
    records
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    // ── Legacy parser unit tests (in-memory helpers) ─────────────────────────

    #[test]
    fn parses_fasta_records_with_normalization() {
        let data = ">seq1\nacgtryswk\n>seq2 something\nTTUU";
        let records = parse_fasta_records(data);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].0, "seq1");
        assert_eq!(records[0].1, "ACGTNNNNN");
        assert_eq!(records[1].0, "seq2");
        assert_eq!(records[1].1, "TTTT");
    }

    #[test]
    fn parses_fastq_records_with_normalization() {
        let data = "@r1\nacgtnu\n+\n!!!!!!\n@r2 extra\nRYSW\n+\n####";
        let records = parse_fastq_records(data);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].0, "r1");
        assert_eq!(records[0].1, "ACGTNT");
        assert_eq!(records[1].0, "r2");
        assert_eq!(records[1].1, "NNNN");
    }

    // ── Streaming FASTA tests ─────────────────────────────────────────────────

    #[test]
    fn stream_fasta_basic() {
        let data = ">s1\nACGT\n>s2 desc\nNNNN\n";
        let reader = std::io::BufReader::new(data.as_bytes());
        let mut out = Vec::new();
        stream_fasta(reader, 1, &mut |r| { out.push((r.id.clone(), r.sequence.clone())); Ok(()) }).unwrap();
        assert_eq!(out.len(), 2);
        assert_eq!(out[0], ("s1".to_string(), "ACGT".to_string()));
        assert_eq!(out[1], ("s2".to_string(), "NNNN".to_string()));
    }

    #[test]
    fn stream_fasta_filters_by_min_len() {
        let data = ">short\nAC\n>long\nACGTACGT\n";
        let reader = std::io::BufReader::new(data.as_bytes());
        let mut out = Vec::new();
        stream_fasta(reader, 4, &mut |r| { out.push(r.id.clone()); Ok(()) }).unwrap();
        assert_eq!(out, vec!["long".to_string()]);
    }

    #[test]
    fn stream_fasta_multiline_sequence() {
        let data = ">multi\nACGT\nTGCA\nACGT\n";
        let reader = std::io::BufReader::new(data.as_bytes());
        let mut out = Vec::new();
        stream_fasta(reader, 1, &mut |r| { out.push(r.sequence.clone()); Ok(()) }).unwrap();
        assert_eq!(out[0], "ACGTTGCAACGT");
    }

    #[test]
    fn stream_fasta_ignores_semicolon_comments() {
        let data = "; comment\n>s1\nACGT\n; another\n";
        let reader = std::io::BufReader::new(data.as_bytes());
        let mut out = Vec::new();
        stream_fasta(reader, 1, &mut |r| { out.push(r.id.clone()); Ok(()) }).unwrap();
        assert_eq!(out, vec!["s1".to_string()]);
    }

    // ── Streaming FASTQ tests ─────────────────────────────────────────────────

    #[test]
    fn stream_fastq_basic() {
        let data = "@r1\nACGT\n+\n!!!!\n@r2\nNNNN\n+\n####\n";
        let reader = std::io::BufReader::new(data.as_bytes());
        let mut out = Vec::new();
        stream_fastq(reader, 1, &mut |r| { out.push((r.id.clone(), r.sequence.clone())); Ok(()) }).unwrap();
        assert_eq!(out.len(), 2);
        assert_eq!(out[0], ("r1".to_string(), "ACGT".to_string()));
    }

    #[test]
    fn stream_fastq_rejects_quality_length_mismatch() {
        // quality string length (3) != sequence length (4)
        let data = "@r1\nACGT\n+\n!!!\n";
        let reader = std::io::BufReader::new(data.as_bytes());
        let result = stream_fastq(reader, 1, &mut |_| Ok(()));
        assert!(result.is_err(), "Expected error for quality/seq length mismatch");
    }

    #[test]
    fn stream_fastq_rejects_missing_plus_separator() {
        let data = "@r1\nACGT\nNOTAPLUS\n!!!!\n";
        let reader = std::io::BufReader::new(data.as_bytes());
        let result = stream_fastq(reader, 1, &mut |_| Ok(()));
        assert!(result.is_err(), "Expected error for missing '+' separator");
    }

    #[test]
    fn stream_fastq_rejects_bad_header() {
        let data = "NOTAT\nACGT\n+\n!!!!\n";
        let reader = std::io::BufReader::new(data.as_bytes());
        let result = stream_fastq(reader, 1, &mut |_| Ok(()));
        assert!(result.is_err(), "Expected error for missing '@' header");
    }

    // ── Format/compression detection ─────────────────────────────────────────

    #[test]
    fn detect_plain_fasta() {
        let (fmt, cmp) = detect_format_and_compression(Path::new("genome.fasta")).unwrap();
        assert_eq!(fmt, DnaFileFormat::Fasta);
        assert_eq!(cmp, Compression::None);
    }

    #[test]
    fn detect_gzip_fastq() {
        let (fmt, cmp) = detect_format_and_compression(Path::new("reads.fastq.gz")).unwrap();
        assert_eq!(fmt, DnaFileFormat::Fastq);
        assert_eq!(cmp, Compression::Gzip);
    }

    #[test]
    fn detect_xz_fasta() {
        let (fmt, cmp) = detect_format_and_compression(Path::new("ref.fa.xz")).unwrap();
        assert_eq!(fmt, DnaFileFormat::Fasta);
        assert_eq!(cmp, Compression::Xz);
    }

    #[test]
    fn detect_unknown_extension_errors() {
        let result = detect_format_and_compression(Path::new("data.bam"));
        assert!(result.is_err());
    }

    // ── Gzip round-trip integration test ─────────────────────────────────────

    #[test]
    fn stream_gz_fasta_round_trip() {
        use flate2::write::GzEncoder;
        use flate2::Compression as GzCompression;
        use std::io::Write as _;

        let fasta_content = b">gz1\nACGTACGT\n>gz2\nTTTTAAAA\n";
        let mut gz_tmp = NamedTempFile::with_suffix(".fasta.gz").unwrap();
        {
            let mut encoder = GzEncoder::new(&mut gz_tmp, GzCompression::default());
            encoder.write_all(fasta_content).unwrap();
            encoder.finish().unwrap();
        }

        let out_tmp = NamedTempFile::new().unwrap();
        let n = build_dna_corpus_jsonl(gz_tmp.path(), out_tmp.path(), 1).unwrap();
        assert_eq!(n, 2, "expected 2 gz FASTA records");
    }
}

