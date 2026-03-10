//! Indexer module - BWT/FM-index construction from source files

mod builder;
mod corpus;
pub mod memory_profiler;
pub mod ops_slo;
pub mod shards;
pub mod vocab;
pub mod ontology;
pub mod sections;
pub mod encrypt;
pub mod alphabet;

pub use builder::IndexBuilder;
pub use builder::MicrogptDomain;
pub use builder::{WasmShardManifest, WasmShardEntry};
pub use corpus::{Corpus, SourceFile, FILE_DELIMITER};
pub use shards::{build_shards, ShardManifest, ShardEntry};
pub use ops_slo::ops_profile;
pub use vocab::{VocabIndex, extract_vocab};
pub use memory_profiler::{BuildMode, MemoryEstimate, estimate_memory};

use thiserror::Error;

#[derive(Error, Debug)]
pub enum IndexError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("No files found matching pattern")]
    NoFilesFound,

    #[error("Serialization error: {0}")]
    Serialization(String),
}

/// Main indexer that orchestrates corpus building and FM-index construction
pub struct Indexer {
    corpus: Corpus,
}

impl Indexer {
    /// Create a new indexer from a directory path
    pub fn from_directory(path: &std::path::Path, extensions: &[&str]) -> Result<Self, IndexError> {
        let corpus = Corpus::from_directory(path, extensions)?;
        Ok(Self { corpus })
    }

    /// Get the corpus
    pub fn corpus(&self) -> &Corpus {
        &self.corpus
    }
}
