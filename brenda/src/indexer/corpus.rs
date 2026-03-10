//! Corpus - collection of source files for indexing

use std::path::{Path, PathBuf};
use walkdir::WalkDir;

use super::IndexError;

/// Delimiter byte for file boundaries (SOH - Start of Header)
/// We use 0x01 instead of 0x00 because fm-index reserves 0x00 for sentinel
pub const FILE_DELIMITER: u8 = 0x01;

/// Represents a file in the corpus with its content and metadata
#[derive(Debug, Clone)]
pub struct SourceFile {
    pub path: PathBuf,
    pub content: String,
    pub start_offset: usize,
    pub end_offset: usize,
}

/// A corpus of source files ready for BWT transformation
#[derive(Debug)]
pub struct Corpus {
    pub files: Vec<SourceFile>,
    pub concatenated: Vec<u8>,
}

impl Corpus {
    /// Build a corpus from a directory, filtering by file extensions
    pub fn from_directory(root: &Path, extensions: &[&str]) -> Result<Self, IndexError> {
        let mut files = Vec::new();
        let mut concatenated = Vec::new();

        for entry in WalkDir::new(root)
            .follow_links(true)
            .into_iter()
            .filter_map(|e| e.ok())
        {
            let path = entry.path();

            // Skip non-files
            if !path.is_file() {
                continue;
            }

            // Check extension
            let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");
            if !extensions.is_empty() && !extensions.contains(&ext) {
                continue;
            }

            // Read file content
            let content = match std::fs::read_to_string(path) {
                Ok(c) => c,
                Err(_) => continue, // Skip binary files or read errors
            };

            // Skip files containing our delimiter byte (very rare in source code)
            if content.as_bytes().contains(&FILE_DELIMITER) {
                continue;
            }

            let relative_path = path.strip_prefix(root).unwrap_or(path);
            let relative_path_str = relative_path.to_string_lossy();
            let start_offset = concatenated.len();

            // Add file marker: \x01<filepath>\x01<content>
            concatenated.push(FILE_DELIMITER);
            concatenated.extend_from_slice(relative_path_str.as_bytes());
            concatenated.push(FILE_DELIMITER);
            concatenated.extend_from_slice(content.as_bytes());

            let end_offset = concatenated.len();

            files.push(SourceFile {
                path: relative_path.to_path_buf(),
                content,
                start_offset,
                end_offset,
            });
        }

        if files.is_empty() {
            return Err(IndexError::NoFilesFound);
        }

        Ok(Self { files, concatenated })
    }

    /// Build a corpus with pre-permuted concatenated bytes and the same file metadata.
    ///
    /// Used by blind-index construction: file offsets are unchanged; only the
    /// concatenated bytes have been alphabet-permuted.
    pub fn from_permuted(original: &Corpus, permuted_bytes: Vec<u8>) -> Self {
        Self {
            files: original.files.clone(),
            concatenated: permuted_bytes,
        }
    }

    /// Get total size of the corpus in bytes
    pub fn size(&self) -> usize {
        self.concatenated.len()
    }

    /// Get number of files in the corpus
    pub fn file_count(&self) -> usize {
        self.files.len()
    }

    /// Find which file contains a given byte offset
    pub fn file_at_offset(&self, offset: usize) -> Option<&SourceFile> {
        self.files
            .iter()
            .find(|f| offset >= f.start_offset && offset < f.end_offset)
    }
}
