use anyhow::{Context, bail};
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Serialize)]
pub struct DnaManifestFile {
    pub path: String,
    pub format: String,
    pub size_bytes: u64,
    pub sha256: String,
    pub split: String,
    pub source_url: Option<String>,
    pub license: Option<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct DnaManifest {
    pub version: String,
    pub input_root: String,
    pub total_files: usize,
    pub total_bytes: u64,
    pub files: Vec<DnaManifestFile>,
}

pub fn build_dna_manifest(
    input_path: &Path,
    split: &str,
    source_url: Option<&str>,
    license: Option<&str>,
) -> anyhow::Result<DnaManifest> {
    let files = collect_sequence_files(input_path)?;
    if files.is_empty() {
        bail!("No FASTA/FASTQ files found under {}", input_path.display());
    }

    let mut entries = Vec::with_capacity(files.len());
    let mut total_bytes = 0u64;

    for file in files {
        let format = detect_format(&file)?;
        let metadata = fs::metadata(&file)
            .with_context(|| format!("Failed to stat {}", file.display()))?;
        let size_bytes = metadata.len();
        total_bytes += size_bytes;

        let digest = sha256_file(&file)?;

        entries.push(DnaManifestFile {
            path: file.to_string_lossy().to_string(),
            format,
            size_bytes,
            sha256: digest,
            split: split.to_string(),
            source_url: source_url.map(ToString::to_string),
            license: license.map(ToString::to_string),
        });
    }

    Ok(DnaManifest {
        version: "v1".to_string(),
        input_root: input_path.to_string_lossy().to_string(),
        total_files: entries.len(),
        total_bytes,
        files: entries,
    })
}

fn sha256_file(path: &Path) -> anyhow::Result<String> {
    let mut file = fs::File::open(path)
        .with_context(|| format!("Failed to open {}", path.display()))?;
    let mut hasher = Sha256::new();
    let mut buffer = [0u8; 64 * 1024];

    loop {
        let read = file
            .read(&mut buffer)
            .with_context(|| format!("Failed to read {}", path.display()))?;
        if read == 0 {
            break;
        }
        hasher.update(&buffer[..read]);
    }

    Ok(format!("{:x}", hasher.finalize()))
}

pub fn write_dna_manifest_json(path: &Path, manifest: &DnaManifest) -> anyhow::Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("Failed to create output directory {}", parent.display()))?;
    }
    let serialized = serde_json::to_string_pretty(manifest)?;
    fs::write(path, serialized).with_context(|| format!("Failed to write {}", path.display()))
}

fn detect_format(path: &Path) -> anyhow::Result<String> {
    let lowercase = path
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();

    if lowercase.ends_with(".fasta") || lowercase.ends_with(".fa") || lowercase.ends_with(".fna") || lowercase.ends_with(".ffn") {
        return Ok("fasta".to_string());
    }
    if lowercase.ends_with(".fastq") || lowercase.ends_with(".fq") {
        return Ok("fastq".to_string());
    }

    bail!("Unsupported sequence file extension: {}", path.display())
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
        if lowercase.ends_with(".fasta")
            || lowercase.ends_with(".fa")
            || lowercase.ends_with(".fna")
            || lowercase.ends_with(".ffn")
            || lowercase.ends_with(".fastq")
            || lowercase.ends_with(".fq")
        {
            files.push(path.to_path_buf());
        }
    }

    files.sort();
    Ok(files)
}

#[cfg(test)]
mod tests {
    use super::{build_dna_manifest, write_dna_manifest_json};
    use tempfile::tempdir;

    #[test]
    fn builds_manifest_with_checksums_and_metadata() {
        let dir = tempdir().expect("tempdir");
        let fasta = dir.path().join("a.fasta");
        let fastq = dir.path().join("b.fastq");

        std::fs::write(&fasta, ">a\nACGT\n").expect("write fasta");
        std::fs::write(&fastq, "@r1\nACGT\n+\n!!!!\n").expect("write fastq");

        let manifest = build_dna_manifest(
            dir.path(),
            "train",
            Some("https://example.org/dataset"),
            Some("CC-BY-4.0"),
        )
        .expect("build manifest");

        assert_eq!(manifest.total_files, 2);
        assert!(manifest.total_bytes > 0);
        assert_eq!(manifest.files[0].split, "train");
        assert!(manifest.files[0].sha256.len() == 64);
        assert_eq!(manifest.files[0].source_url.as_deref(), Some("https://example.org/dataset"));
        assert_eq!(manifest.files[0].license.as_deref(), Some("CC-BY-4.0"));
    }

    #[test]
    fn writes_manifest_json() {
        let dir = tempdir().expect("tempdir");
        let fasta = dir.path().join("only.fa");
        std::fs::write(&fasta, ">s\nACGT\n").expect("write fasta");

        let manifest = build_dna_manifest(dir.path(), "val", None, None).expect("build manifest");
        let out = dir.path().join("manifest.json");
        write_dna_manifest_json(&out, &manifest).expect("write manifest");

        let raw = std::fs::read_to_string(out).expect("read manifest");
        assert!(raw.contains("\"version\": \"v1\""));
        assert!(raw.contains("\"split\": \"val\""));
    }
}
