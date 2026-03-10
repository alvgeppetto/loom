use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

use fm_index::converter::RangeConverter;
use fm_index::suffix_array::SuffixOrderSampler;
use fm_index::{BackwardSearchIndex, FMIndex};
use serde::{Deserialize, Serialize};

use super::vocab::{deserialize_vocab, extract_vocab, serialize_vocab, VocabIndex};
use super::ontology::OntologySection;
use super::{Corpus, IndexError};

/// 4-byte section tags for optional data after the core corpus.
const SECTION_TAG_ONTOLOGY: &[u8; 4] = b"ONTO";
/// 4-byte tag identifying the vocabulary section.
const TAG_VOCAB: [u8; 4] = *b"VCAB";
/// Sentinel tag (all-zeros) indicating no more sections follow.
const TAG_END: [u8; 4] = [0u8; 4];
/// Section tag marking a blind (alphabet-permuted) index.
const SECTION_TAG_BLIND: &[u8; 4] = b"BLIN";
/// Section tag for embedded MicroGPT weights (tokenizer + model).
const TAG_MGPT: [u8; 4] = *b"MGPT";
/// Section tag for pre-serialized FM-index. Allows WASM to
/// deserialize the index without rebuilding from the raw corpus.
pub const TAG_FMIX: [u8; 4] = *b"FMIX";

/// 4-byte magic prefix inside FMIX payloads written with postcard.
/// Used to distinguish postcard payloads from stale bincode ones.
const FMIX_MAGIC: &[u8; 4] = b"PCv1";

/// Serializable index metadata
#[derive(Serialize, Deserialize)]
pub struct IndexMetadata {
    files: Vec<FileMetadata>,
    corpus_size: usize,
}

#[derive(Serialize, Deserialize, Clone)]
pub struct FileMetadata {
    path: String,
    start_offset: usize,
    end_offset: usize,
}

/// Builds and manages the FM-index
pub struct IndexBuilder {
    index: FMIndex<
        u8,
        RangeConverter<u8>,
        fm_index::suffix_array::SuffixOrderSampledArray,
    >,
    metadata: IndexMetadata,
    corpus_text: Vec<u8>,
    /// Optional tagged sections appended after the corpus.
    /// Key = 4-byte ASCII tag (e.g. "ONTO"), Value = raw bytes.
    sections: HashMap<[u8; 4], Vec<u8>>,
    /// Typed vocabulary index, derived from the VCAB section when present.
    pub vocab: Option<VocabIndex>,
}

impl IndexBuilder {
    /// Create a new index builder from a corpus (no vocab section).
    pub fn new(corpus: &Corpus) -> Result<Self, IndexError> {
        Self::new_with_vocab(corpus, false)
    }

    /// Create a new index builder from a corpus, optionally building a vocab section.
    pub fn new_with_vocab(corpus: &Corpus, build_vocab: bool) -> Result<Self, IndexError> {
        let text = corpus.concatenated.clone();

        let converter = RangeConverter::new(1u8, 255u8);
        let sampler = SuffixOrderSampler::new().level(2);
        let index = FMIndex::new(text.clone(), converter, sampler);

        let metadata = IndexMetadata {
            files: corpus
                .files
                .iter()
                .map(|f| FileMetadata {
                    path: f.path.to_string_lossy().to_string(),
                    start_offset: f.start_offset,
                    end_offset: f.end_offset,
                })
                .collect(),
            corpus_size: corpus.size(),
        };

        // Build vocabulary from all file content if requested.
        let vocab = if build_vocab {
            let all_text: String = corpus.files.iter().map(|f| f.content.as_str()).collect::<Vec<_>>().join("\n");
            let terms = extract_vocab(&all_text);
            Some(VocabIndex::build(terms))
        } else {
            None
        };

        Ok(Self {
            index,
            metadata,
            corpus_text: text,
            sections: HashMap::new(),
            vocab,
        })
    }

    /// Create a new blind index builder from a permuted corpus.
    ///
    /// Identical to [`Self::new`] but stamps a `BLIN` section tag so searchers
    /// know to apply the trapdoor permutation before querying.
    pub fn new_blind(corpus: &Corpus) -> Result<Self, IndexError> {
        let mut builder = Self::new(corpus)?;
        builder.set_section(*SECTION_TAG_BLIND, b"1".to_vec());
        Ok(builder)
    }

    /// Returns `true` if this index was built with alphabet permutation (blind index).
    pub fn is_blind(&self) -> bool {
        self.sections.contains_key(SECTION_TAG_BLIND)
    }

    /// Count occurrences of a raw byte pattern (used for blind/permuted queries).
    pub fn count_bytes(&self, pattern: &[u8]) -> u64 {
        self.index.search_backward(pattern).count()
    }

    /// Search for a raw byte pattern and return positions (used for blind/permuted queries).
    pub fn search_bytes(&self, pattern: &[u8]) -> Vec<u64> {
        self.index.search_backward(pattern).locate()
    }

    /// Count occurrences of a pattern
    pub fn count(&self, pattern: &str) -> u64 {
        self.index.search_backward(pattern.as_bytes()).count()
    }

    /// Search for a pattern and return positions
    pub fn search(&self, pattern: &str) -> Vec<u64> {
        self.index.search_backward(pattern.as_bytes()).locate()
    }

    /// Get the file containing a given offset
    pub fn file_at_offset(&self, offset: u64) -> Option<&FileMetadata> {
        self.metadata
            .files
            .iter()
            .find(|f| offset >= f.start_offset as u64 && offset < f.end_offset as u64)
    }

    /// Get line and column for an offset within a file
    pub fn line_col_at_offset(&self, offset: u64) -> Option<(String, usize, usize)> {
        let file = self.file_at_offset(offset)?;
        let offset = offset as usize;

        let content_start = file.start_offset + file.path.len() + 2;
        if offset < content_start {
            return None;
        }

        let relative_offset = offset - content_start;
        let content = &self.corpus_text[content_start..offset];
        let line = content.iter().filter(|&&b| b == b'\n').count() + 1;

        let last_newline = content.iter().rposition(|&b| b == b'\n');
        let column = match last_newline {
            Some(pos) => relative_offset - pos,
            None => relative_offset + 1,
        };

        Some((file.path.clone(), line, column))
    }

    /// Get context around an offset
    pub fn context_at_offset(&self, offset: u64, chars_before: usize, chars_after: usize) -> String {
        let offset = offset as usize;
        let start = offset.saturating_sub(chars_before);
        let end = (offset + chars_after).min(self.corpus_text.len());
        String::from_utf8_lossy(&self.corpus_text[start..end]).to_string()
    }

    /// Get index metadata
    pub fn metadata(&self) -> &IndexMetadata {
        &self.metadata
    }

    /// Get number of indexed files
    pub fn file_count(&self) -> usize {
        self.metadata.files.len()
    }

    /// Get corpus size
    pub fn corpus_size(&self) -> usize {
        self.metadata.corpus_size
    }

    // ── Tagged-section I/O helpers ───────────────────────────────────────────

    fn write_section<W: Write>(writer: &mut W, tag: &[u8; 4], data: &[u8]) -> Result<(), IndexError> {
        writer.write_all(tag)?;
        writer.write_all(&(data.len() as u64).to_le_bytes())?;
        writer.write_all(data)?;
        Ok(())
    }

    /// Save index to disk (with optional VCAB section if `self.vocab` is set).
    pub fn save(&self, path: &Path) -> Result<(), IndexError> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // Write metadata as JSON
        let metadata_json = serde_json::to_vec(&self.metadata)
            .map_err(|e| IndexError::Serialization(e.to_string()))?;
        let metadata_len = metadata_json.len() as u64;
        writer.write_all(&metadata_len.to_le_bytes())?;
        writer.write_all(&metadata_json)?;

        // Write corpus text
        let corpus_len = self.corpus_text.len() as u64;
        writer.write_all(&corpus_len.to_le_bytes())?;
        writer.write_all(&self.corpus_text)?;

        // Write optional tagged sections: [4-byte tag][8-byte len][data...]
        // First, emit generic sections (e.g. ONTO from the ontology builder).
        for (tag, data) in &self.sections {
            writer.write_all(tag)?;
            let data_len = data.len() as u64;
            writer.write_all(&data_len.to_le_bytes())?;
            writer.write_all(data)?;
        }

        // Write optional VCAB section (typed vocab index).
        if let Some(ref v) = self.vocab {
            let mut payload: Vec<u8> = Vec::new();
            serialize_vocab(v, &mut payload)
                .map_err(|e| IndexError::Serialization(e.to_string()))?;
            Self::write_section(&mut writer, &TAG_VOCAB, &payload)?;
        }

        // Write end sentinel so readers know no more sections follow.
        writer.write_all(&TAG_END)?;

        writer.flush()?;
        Ok(())
    }

    /// Load index from disk (reads VCAB section if present; backward-compatible).
    pub fn load(path: &Path) -> Result<Self, IndexError> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);

        // Read metadata
        let mut len_buf = [0u8; 8];
        reader.read_exact(&mut len_buf)?;
        let metadata_len = u64::from_le_bytes(len_buf) as usize;

        let mut metadata_buf = vec![0u8; metadata_len];
        reader.read_exact(&mut metadata_buf)?;
        let metadata: IndexMetadata = serde_json::from_slice(&metadata_buf)
            .map_err(|e| IndexError::Serialization(e.to_string()))?;

        // Read corpus text
        reader.read_exact(&mut len_buf)?;
        let corpus_len = u64::from_le_bytes(len_buf) as usize;

        let mut corpus_text = vec![0u8; corpus_len];
        reader.read_exact(&mut corpus_text)?;

        // Read optional tagged sections (backward compatible: old files simply end here).
        // A new-format file ends with the TAG_END sentinel [0,0,0,0].
        let mut sections = HashMap::new();
        let mut vocab: Option<VocabIndex> = None;
        loop {
            let mut tag_buf = [0u8; 4];
            match reader.read_exact(&mut tag_buf) {
                Ok(()) => {}
                Err(ref e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(IndexError::Io(e)),
            }
            // End-of-sections sentinel (new format).
            if tag_buf == TAG_END {
                break;
            }
            match reader.read_exact(&mut len_buf) {
                Ok(()) => {}
                Err(ref e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(IndexError::Io(e)),
            }
            let section_len = u64::from_le_bytes(len_buf) as usize;
            let mut section_data = vec![0u8; section_len];
            reader.read_exact(&mut section_data)?;

            // Decode typed sections we know about; store all for generic access.
            if tag_buf == TAG_VOCAB {
                vocab = Some(deserialize_vocab(&section_data));
            }
            sections.insert(tag_buf, section_data);
        }

        // Use pre-serialized FM-index (FMIX section) when available;
        // otherwise rebuild from the corpus.
        let index = if let Some(fmix_data) = sections.get(&TAG_FMIX) {
            if fmix_data.len() > 4 && &fmix_data[..4] == FMIX_MAGIC {
                postcard::from_bytes(&fmix_data[4..])
                    .map_err(|e| IndexError::Serialization(format!("FM-index deserialize: {e}")))?
            } else {
                // Stale bincode payload — rebuild from corpus.
                let converter = RangeConverter::new(1u8, 255u8);
                let sampler = SuffixOrderSampler::new().level(2);
                FMIndex::new(corpus_text.clone(), converter, sampler)
            }
        } else {
            let converter = RangeConverter::new(1u8, 255u8);
            let sampler = SuffixOrderSampler::new().level(2);
            FMIndex::new(corpus_text.clone(), converter, sampler)
        };

        Ok(Self {
            index,
            metadata,
            corpus_text,
            sections,
            vocab,
        })
    }

    /// Serialize the FM-index into a `FMIX` section so WASM can
    /// deserialize it directly without rebuilding from the corpus.
    pub fn prepare_wasm(&mut self) -> Result<(), IndexError> {
        let payload = postcard::to_stdvec(&self.index)
            .map_err(|e| IndexError::Serialization(format!("FM-index serialize: {e}")))?;
        let mut fmix_data = Vec::with_capacity(4 + payload.len());
        fmix_data.extend_from_slice(FMIX_MAGIC);
        fmix_data.extend_from_slice(&payload);
        self.sections.insert(TAG_FMIX, fmix_data);
        Ok(())
    }

    /// Check whether this index has a pre-built FMIX section.
    pub fn has_fmix(&self) -> bool {
        self.sections.contains_key(&TAG_FMIX)
    }

    /// Get an optional section by its 4-byte tag.
    pub fn get_section(&self, tag: &[u8; 4]) -> Option<&[u8]> {
        self.sections.get(tag).map(|v| v.as_slice())
    }

    /// Set an optional section by its 4-byte tag.
    pub fn set_section(&mut self, tag: [u8; 4], data: Vec<u8>) {
        self.sections.insert(tag, data);
    }

    /// Get the raw ontology section bytes, if present.
    pub fn ontology_bytes(&self) -> Option<&[u8]> {
        self.get_section(SECTION_TAG_ONTOLOGY)
    }

    /// Set the ontology section from raw bytes.
    pub fn set_ontology_bytes(&mut self, data: Vec<u8>) {
        self.set_section(*SECTION_TAG_ONTOLOGY, data);
    }

    /// Deserialize and return the typed ontology section, if present.
    pub fn ontology(&self) -> Option<OntologySection> {
        self.ontology_bytes()
            .and_then(|b| OntologySection::from_bytes(b).ok())
    }

    /// Store a typed ontology section (serialized to bytes).
    pub fn set_ontology(&mut self, onto: OntologySection) {
        self.set_ontology_bytes(onto.to_bytes());
    }

    /// Extract (path, content) pairs from the corpus for ontology building.
    pub fn extract_files_for_ontology(&self) -> Vec<(String, String)> {
        self.metadata
            .files
            .iter()
            .filter_map(|f| {
                let start = f.start_offset;
                let end = f.end_offset.min(self.corpus_text.len());
                if start >= end {
                    return None;
                }
                let content = String::from_utf8_lossy(&self.corpus_text[start..end]).into_owned();
                Some((f.path.clone(), content))
            })
            .collect()
    }

    /// Get the raw MGPT section bytes, if present.
    pub fn microgpt_bytes(&self) -> Option<&[u8]> {
        self.get_section(&TAG_MGPT)
    }

    /// Returns `true` if this index has an embedded MicroGPT model.
    pub fn has_microgpt(&self) -> bool {
        self.sections.contains_key(&TAG_MGPT)
    }

    /// Set the MGPT section from raw bytes.
    pub fn set_microgpt_bytes(&mut self, data: Vec<u8>) {
        self.set_section(TAG_MGPT, data);
    }
}

/// Manifest describing a set of WASM-ready shards.
#[derive(Debug, Serialize, Deserialize)]
pub struct WasmShardManifest {
    pub shards: Vec<WasmShardEntry>,
    pub total_corpus_bytes: usize,
    pub total_files: usize,
    pub overlap_bytes: usize,
}

/// One shard in a WASM shard manifest.
#[derive(Debug, Serialize, Deserialize)]
pub struct WasmShardEntry {
    pub id: usize,
    pub path: String,
    pub corpus_start: usize,
    pub corpus_end: usize,
    pub file_count: usize,
    pub size_bytes: u64,
}

/// Domain selector for embedded MicroGPT training.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MicrogptDomain {
    /// DNA/genomic corpus (character-level tokenizer, small vocab)
    Dna,
    /// Source code corpus (term-level tokenizer, ~4096 vocab)
    Code,
    /// Technical manual corpus (term-level tokenizer, ~2048 vocab)
    Manual,
    /// Auto-detect from corpus content
    Auto,
}

impl IndexBuilder {
    /// Split this index into WASM-ready shards, each ≤ `target_shard_bytes`.
    ///
    /// Each shard gets its own FM-index + FMIX section so a Web Worker can
    /// load it via `LoomWasmIndex.fromSerialized()` independently.
    pub fn prepare_wasm_shards(
        &self,
        output_dir: &Path,
        target_shard_bytes: usize,
    ) -> Result<WasmShardManifest, IndexError> {
        use std::fs;
        fs::create_dir_all(output_dir)
            .map_err(|e| IndexError::Io(e))?;

        let corpus = &self.corpus_text;
        let total = corpus.len();
        let overlap: usize = 1024; // bytes of overlap between adjacent shards

        // Compute shard boundaries
        let mut boundaries: Vec<(usize, usize)> = Vec::new();
        let mut start = 0usize;
        while start < total {
            let end = (start + target_shard_bytes).min(total);
            boundaries.push((start, end));
            start = if end == total { total } else { end.saturating_sub(overlap) };
        }

        let mut entries = Vec::new();
        for (id, &(shard_start, shard_end)) in boundaries.iter().enumerate() {
            let shard_corpus = &corpus[shard_start..shard_end];

            // Count files whose start_offset falls within this shard
            let file_count = self.metadata.files.iter()
                .filter(|f| f.start_offset >= shard_start && f.start_offset < shard_end)
                .count();

            // Build per-shard FM-index
            let converter = RangeConverter::new(1u8, 255u8);
            let sampler = SuffixOrderSampler::new().level(2);
            let shard_index = FMIndex::new(shard_corpus.to_vec(), converter, sampler);

            // Serialize FM-index with postcard + PCv1 magic
            let payload = postcard::to_stdvec(&shard_index)
                .map_err(|e| IndexError::Serialization(format!("shard FM-index serialize: {e}")))?;
            let mut fmix_data = Vec::with_capacity(4 + payload.len());
            fmix_data.extend_from_slice(FMIX_MAGIC);
            fmix_data.extend_from_slice(&payload);

            // Build shard metadata
            let shard_files: Vec<FileMetadata> = self.metadata.files.iter()
                .filter(|f| f.start_offset >= shard_start && f.start_offset < shard_end)
                .map(|f| FileMetadata {
                    path: f.path.clone(),
                    start_offset: f.start_offset - shard_start,
                    end_offset: f.end_offset.min(shard_end) - shard_start,
                })
                .collect();
            let shard_meta = IndexMetadata {
                files: shard_files,
                corpus_size: shard_corpus.len(),
            };

            // Build a minimal IndexBuilder and save it
            let mut shard_builder = IndexBuilder {
                index: shard_index,
                metadata: shard_meta,
                corpus_text: shard_corpus.to_vec(),
                sections: HashMap::new(),
                vocab: None,
            };
            shard_builder.sections.insert(TAG_FMIX, fmix_data);

            let shard_name = format!("shard_{:04}.idx", id);
            let shard_path = output_dir.join(&shard_name);
            shard_builder.save(&shard_path)?;

            let size_bytes = fs::metadata(&shard_path)
                .map(|m| m.len())
                .unwrap_or(0);

            entries.push(WasmShardEntry {
                id,
                path: shard_name,
                corpus_start: shard_start,
                corpus_end: shard_end,
                file_count,
                size_bytes,
            });
        }

        let manifest = WasmShardManifest {
            shards: entries,
            total_corpus_bytes: total,
            total_files: self.metadata.files.len(),
            overlap_bytes: overlap,
        };

        // Write manifest JSON
        let manifest_path = output_dir.join("shard_manifest.json");
        let manifest_json = serde_json::to_string_pretty(&manifest)
            .map_err(|e| IndexError::Serialization(format!("manifest JSON: {e}")))?;
        fs::write(&manifest_path, manifest_json)
            .map_err(|e| IndexError::Io(e))?;

        Ok(manifest)
    }

    /// Extract vocabulary terms from the corpus for tokenizer construction.
    pub fn extract_vocab_terms(&self) -> Vec<String> {
        if let Some(ref v) = self.vocab {
            v.terms().to_vec()
        } else {
            let all_text: String = self.extract_files_for_ontology()
                .iter()
                .map(|(_, c)| c.as_str())
                .collect::<Vec<_>>()
                .join("\n");
            extract_vocab(&all_text)
        }
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::corpus::SourceFile;
    use crate::indexer::Corpus;
    use std::path::PathBuf;
    use tempfile::NamedTempFile;

    /// Build a minimal IndexBuilder from a small corpus string.
    fn tiny_builder(text: &str) -> IndexBuilder {
        let path = PathBuf::from("test.rs");
        let marker = format!("\x01test.rs\x01{}", text);
        let concatenated = marker.as_bytes().to_vec();
        let end = concatenated.len();
        let corpus = Corpus {
            files: vec![SourceFile {
                path,
                content: text.to_string(),
                start_offset: 0,
                end_offset: end,
            }],
            concatenated,
        };
        IndexBuilder::new(&corpus).expect("builder")
    }

    #[test]
    fn save_and_load_no_sections() {
        let builder = tiny_builder("hello world");
        let tmp = NamedTempFile::new().unwrap();
        builder.save(tmp.path()).unwrap();

        let loaded = IndexBuilder::load(tmp.path()).unwrap();
        assert_eq!(loaded.corpus_text, builder.corpus_text);
        assert!(loaded.sections.is_empty());
    }

    #[test]
    fn save_and_load_single_section() {
        let mut builder = tiny_builder("hello world");
        builder.set_section(*b"TEST", b"section-data".to_vec());

        let tmp = NamedTempFile::new().unwrap();
        builder.save(tmp.path()).unwrap();

        let loaded = IndexBuilder::load(tmp.path()).unwrap();
        assert_eq!(loaded.get_section(b"TEST"), Some(b"section-data".as_slice()));
    }

    #[test]
    fn save_and_load_multiple_sections() {
        let mut builder = tiny_builder("hello world");
        builder.set_section(*b"AAA\x00", b"aaa-data".to_vec());
        builder.set_section(*b"BBB\x00", b"bbb-data".to_vec());

        let tmp = NamedTempFile::new().unwrap();
        builder.save(tmp.path()).unwrap();

        let loaded = IndexBuilder::load(tmp.path()).unwrap();
        assert_eq!(loaded.get_section(b"AAA\x00"), Some(b"aaa-data".as_slice()));
        assert_eq!(loaded.get_section(b"BBB\x00"), Some(b"bbb-data".as_slice()));
        assert_eq!(loaded.sections.len(), 2);
    }

    #[test]
    fn ontology_section_round_trip() {
        let mut builder = tiny_builder("fn deliver() {}");
        let onto_json = br#"{"deliver":["ack","confirm"]}"#;
        builder.set_ontology_bytes(onto_json.to_vec());

        let tmp = NamedTempFile::new().unwrap();
        builder.save(tmp.path()).unwrap();

        let loaded = IndexBuilder::load(tmp.path()).unwrap();
        assert_eq!(loaded.ontology_bytes(), Some(onto_json.as_slice()));
    }

    #[test]
    fn backward_compat_old_index_no_sections() {
        // Simulate reading a file written by old code that has no section data.
        // We save normally (no sections set) and verify load succeeds with empty sections.
        let builder = tiny_builder("old format");
        let tmp = NamedTempFile::new().unwrap();
        builder.save(tmp.path()).unwrap();

        let loaded = IndexBuilder::load(tmp.path()).unwrap();
        assert!(loaded.sections.is_empty());
        assert_eq!(loaded.get_section(b"ONTO"), None);
    }

    #[test]
    fn test_save_load_without_vocab() {
        let builder = tiny_builder("fn hello_world() { let myVar = 42; }");
        let tmp = NamedTempFile::new().unwrap();
        builder.save(tmp.path()).unwrap();

        let loaded = IndexBuilder::load(tmp.path()).unwrap();
        assert_eq!(loaded.file_count(), 1);
        assert!(loaded.vocab.is_none());
    }

    #[test]
    fn test_save_load_with_vocab_section() {
        let path = PathBuf::from("test.rs");
        let content = "fn hello_world() { let myVar = 42; }".to_string();
        let marker = format!("\x01test.rs\x01{}", content);
        let concatenated = marker.as_bytes().to_vec();
        let end_offset = concatenated.len();
        let corpus = Corpus {
            files: vec![SourceFile {
                path,
                content,
                start_offset: 0,
                end_offset,
            }],
            concatenated,
        };
        let builder = IndexBuilder::new_with_vocab(&corpus, true).unwrap();
        assert!(builder.vocab.is_some());

        let tmp = NamedTempFile::new().unwrap();
        builder.save(tmp.path()).unwrap();

        let loaded = IndexBuilder::load(tmp.path()).unwrap();
        let v = loaded.vocab.expect("VCAB section should be present after round-trip");
        assert!(v.exists("hello"), "expected 'hello' in vocab");
        assert!(v.exists("world"), "expected 'world' in vocab");
        assert!(v.exists("my"), "expected 'my' in vocab");
        assert!(v.exists("var"), "expected 'var' in vocab");
    }

    #[test]
    fn mgpt_backward_compat_no_mgpt() {
        // Index without MGPT section should load fine.
        let builder = tiny_builder("fn main() {}");
        assert!(!builder.has_microgpt());

        let tmp = NamedTempFile::new().unwrap();
        builder.save(tmp.path()).unwrap();

        let loaded = IndexBuilder::load(tmp.path()).unwrap();
        assert!(!loaded.has_microgpt());
        assert!(loaded.microgpt_bytes().is_none());
    }
}

