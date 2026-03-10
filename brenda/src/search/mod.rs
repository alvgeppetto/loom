//! Search module - query execution against FM-index

mod query;
mod results;
pub mod federated;

pub use query::Query;
pub use results::{rank_results, ScoredFile, SearchResult, SearchResults};
pub use federated::{FederatedSearcher, FederatedMatch, FederatedStats};

use thiserror::Error;

#[derive(Error, Debug)]
pub enum SearchError {
    #[error("Index not loaded")]
    IndexNotLoaded,

    #[error("Invalid query: {0}")]
    InvalidQuery(String),
}

/// Searcher for executing queries against an FM-index
pub struct Searcher {
    // Will hold reference to FM-index and corpus metadata
}

impl Searcher {
    /// Execute a search query
    pub fn search(&self, _query: &str) -> Result<SearchResults, SearchError> {
        // TODO: Implement FM-index search
        Ok(SearchResults::empty())
    }

    /// Count occurrences of a pattern
    pub fn count(&self, _pattern: &str) -> Result<usize, SearchError> {
        // TODO: Implement FM-index count
        Ok(0)
    }
}
