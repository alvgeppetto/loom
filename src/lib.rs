//! LOOM - BWT-Powered Exact-Retrieval Engine
//!
//! Open-science CRISPR target discovery and pathogen genomics search.
//! Core indexing and search functionality is provided by the `brenda` crate.

pub mod cli;
pub mod dna;
pub mod wasm;

// Re-export brenda's public API so existing code works unchanged.
pub use brenda::indexer;
pub use brenda::search;
pub use brenda::Indexer;
pub use brenda::{VocabIndex, extract_vocab};
pub use brenda::Searcher;
