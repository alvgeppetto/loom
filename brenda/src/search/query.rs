//! Query representation

/// A search query with optional context
#[derive(Debug, Clone)]
pub struct Query {
    pub pattern: String,
    pub context_lines: usize,
    pub max_results: Option<usize>,
}

impl Query {
    pub fn new(pattern: impl Into<String>) -> Self {
        Self {
            pattern: pattern.into(),
            context_lines: 0,
            max_results: None,
        }
    }

    pub fn with_context(mut self, lines: usize) -> Self {
        self.context_lines = lines;
        self
    }

    pub fn with_max_results(mut self, max: usize) -> Self {
        self.max_results = Some(max);
        self
    }
}
