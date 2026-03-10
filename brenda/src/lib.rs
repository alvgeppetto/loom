//! Brenda — BWT/FM-index construction and search engine
//!
//! A standalone library for building and querying FM-index based text indices.
//! Supports corpus construction, vocabulary extraction, ontology building,
//! sharded indexing, encrypted indices, and federated search.

pub mod indexer;
pub mod search;

pub use indexer::Indexer;
pub use indexer::{VocabIndex, extract_vocab};
pub use indexer::{IndexBuilder, IndexError, Corpus, SourceFile, FILE_DELIMITER};
pub use indexer::{build_shards, ShardManifest, ShardEntry};
pub use indexer::{BuildMode, MemoryEstimate, estimate_memory};
pub use indexer::MicrogptDomain;
pub use search::Searcher;
pub use search::{FederatedSearcher, FederatedMatch, FederatedStats};
pub use search::{Query, SearchResult, SearchResults, ScoredFile, rank_results};
