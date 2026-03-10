//! Vocabulary index — extracts identifiers from corpus, builds secondary FM-index.
//!
//! V1.1: token extraction from raw corpus text.
//! V1.2: VocabIndex wrapping an FM-index for prefix/exact lookups.

use fm_index::converter::RangeConverter;
use fm_index::suffix_array::SuffixOrderSampler;
use fm_index::{FMIndex};

// ── V1.1: token extractor ────────────────────────────────────────────────────

/// Split a camelCase or PascalCase identifier into lowercase subwords.
///
/// `fooBar`  → `["foo", "bar"]`
/// `XMLParser` → `["x", "m", "l", "parser"]`  (naive – each upper starts new word)
fn split_camel(s: &str) -> Vec<String> {
    let mut words: Vec<String> = Vec::new();
    let mut current = String::new();

    for ch in s.chars() {
        if ch.is_uppercase() && !current.is_empty() {
            words.push(current.to_lowercase());
            current = String::new();
        }
        current.push(ch);
    }
    if !current.is_empty() {
        words.push(current.to_lowercase());
    }
    words
}

/// Split a snake_case identifier into lowercase subwords.
///
/// `foo_bar` → `["foo", "bar"]`
fn split_snake(s: &str) -> Vec<String> {
    s.split('_')
        .filter(|p| !p.is_empty())
        .map(|p| p.to_lowercase())
        .collect()
}

/// Extract, split, deduplicate and sort all code identifiers from `text`.
///
/// Returns a sorted `Vec<String>` of unique lowercase subwords.
pub fn extract_vocab(text: &str) -> Vec<String> {
    // Simple hand-rolled regex equivalent: [A-Za-z_][A-Za-z0-9_]*
    let mut tokens: Vec<String> = Vec::new();
    let bytes = text.as_bytes();
    let mut i = 0;
    let n = bytes.len();

    while i < n {
        let b = bytes[i];
        if b.is_ascii_alphabetic() || b == b'_' {
            let start = i;
            i += 1;
            while i < n && (bytes[i].is_ascii_alphanumeric() || bytes[i] == b'_') {
                i += 1;
            }
            // SAFETY: only ASCII bytes selected above
            let token = &text[start..i];

            // Emit subwords from both splits
            let has_upper = token.chars().any(|c| c.is_uppercase());
            let has_under = token.contains('_');

            if has_under {
                tokens.extend(split_snake(token));
            } else if has_upper {
                tokens.extend(split_camel(token));
            } else {
                tokens.push(token.to_lowercase());
            }
        } else {
            i += 1;
        }
    }

    // Deduplicate and sort
    tokens.sort_unstable();
    tokens.dedup();
    // Filter single-char noise and pure underscores
    tokens.retain(|t| t.len() > 1);
    tokens
}

// ── V1.2: VocabIndex ────────────────────────────────────────────────────────

type VocabFMIndex = FMIndex<u8, RangeConverter<u8>, fm_index::suffix_array::SuffixOrderSampledArray>;

/// Secondary FM-index built over the sorted vocabulary.
pub struct VocabIndex {
    /// Sorted vocabulary terms (for binary-search prefix/exists).
    sorted_terms: Vec<String>,
    /// FM-index over newline-separated sorted terms (for future pattern queries).
    _fm_index: VocabFMIndex,
    /// The raw text fed to the FM-index (`term\n...`).
    _text: Vec<u8>,
}

impl VocabIndex {
    /// Build a `VocabIndex` from a list of already-sorted, deduplicated terms.
    pub fn build(terms: Vec<String>) -> Self {
        // Concatenate as `term\n` — each term ends with newline so boundaries are clear.
        let mut text: Vec<u8> = Vec::new();
        for t in &terms {
            text.extend_from_slice(t.as_bytes());
            text.push(b'\n');
        }
        // Ensure non-empty (FM-index panics on empty input)
        if text.is_empty() {
            text.push(b'\n');
        }

        let converter = RangeConverter::new(1u8, 255u8);
        let sampler = SuffixOrderSampler::new().level(2);
        let fm_index = FMIndex::new(text.clone(), converter, sampler);

        Self {
            sorted_terms: terms,
            _fm_index: fm_index,
            _text: text,
        }
    }

    /// Return all terms whose prefix matches `prefix` (at most `limit` results).
    pub fn prefix_search(&self, prefix: &str, limit: usize) -> Vec<String> {
        if prefix.is_empty() {
            return self.sorted_terms.iter().take(limit).cloned().collect();
        }
        // Binary-search for the start of the prefix range.
        let start = self
            .sorted_terms
            .partition_point(|t| t.as_str() < prefix);
        self.sorted_terms[start..]
            .iter()
            .take_while(|t| t.starts_with(prefix))
            .take(limit)
            .cloned()
            .collect()
    }

    /// Check whether an exact term exists in the vocabulary.
    pub fn exists(&self, term: &str) -> bool {
        self.sorted_terms.binary_search_by(|t| t.as_str().cmp(term)).is_ok()
    }

    /// Number of unique terms.
    pub fn len(&self) -> usize {
        self.sorted_terms.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sorted_terms.is_empty()
    }

    /// Borrow the sorted term list.
    pub fn terms(&self) -> &[String] {
        &self.sorted_terms
    }
}

// ── Serialization helpers (used by V1.3) ────────────────────────────────────

use std::io::Write;

/// Serialize a `VocabIndex` into `writer` (newline-separated terms, no framing).
pub fn serialize_vocab<W: Write>(idx: &VocabIndex, writer: &mut W) -> std::io::Result<()> {
    // Store sorted terms as newline-separated UTF-8.
    let payload: Vec<u8> = idx
        .sorted_terms
        .iter()
        .flat_map(|t| {
            let mut v = t.as_bytes().to_vec();
            v.push(b'\n');
            v
        })
        .collect();
    writer.write_all(&payload)
}

/// Deserialize a `VocabIndex` from a byte slice (newline-separated terms).
pub fn deserialize_vocab(data: &[u8]) -> VocabIndex {
    let text = String::from_utf8_lossy(data);
    let terms: Vec<String> = text
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| l.to_string())
        .collect();
    VocabIndex::build(terms)
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── V1.5 token extraction tests ──────────────────────────────────────────

    #[test]
    fn test_snake_case_split() {
        let mut v = split_snake("foo_bar_baz");
        v.sort_unstable();
        assert_eq!(v, vec!["bar", "baz", "foo"]);
    }

    #[test]
    fn test_camel_case_split() {
        let words = split_camel("fooBarBaz");
        assert_eq!(words, vec!["foo", "bar", "baz"]);
    }

    #[test]
    fn test_pascal_case_split() {
        let words = split_camel("MyStructName");
        assert_eq!(words, vec!["my", "struct", "name"]);
    }

    #[test]
    fn test_extract_vocab_basic() {
        let text = "fn foo_bar() { let myValue = 42; }";
        let vocab = extract_vocab(text);
        assert!(vocab.contains(&"foo".to_string()));
        assert!(vocab.contains(&"bar".to_string()));
        assert!(vocab.contains(&"my".to_string()));
        assert!(vocab.contains(&"value".to_string()));
        assert!(vocab.contains(&"fn".to_string()));
        assert!(vocab.contains(&"let".to_string()));
    }

    #[test]
    fn test_extract_vocab_dedup_and_sort() {
        let text = "foo foo bar foo";
        let vocab = extract_vocab(text);
        assert_eq!(vocab.iter().filter(|t| *t == "foo").count(), 1);
        let is_sorted: bool = vocab.windows(2).all(|w| w[0] <= w[1]);
        assert!(is_sorted, "vocab should be sorted");
    }

    #[test]
    fn test_extract_vocab_no_single_chars() {
        let text = "a b c fn let";
        let vocab = extract_vocab(text);
        // 'a', 'b', 'c' should be filtered (len <= 1)
        for t in &vocab {
            assert!(t.len() > 1, "single-char token leaked: {t}");
        }
    }

    // ── V1.5 VocabIndex tests ────────────────────────────────────────────────

    fn sample_index() -> VocabIndex {
        let terms: Vec<String> = vec![
            "bar", "baz", "foo", "foobar", "foobaz", "qux",
        ]
        .into_iter()
        .map(String::from)
        .collect();
        VocabIndex::build(terms)
    }

    #[test]
    fn test_exists_true() {
        let idx = sample_index();
        assert!(idx.exists("foo"));
        assert!(idx.exists("bar"));
    }

    #[test]
    fn test_exists_false() {
        let idx = sample_index();
        assert!(!idx.exists("nothere"));
    }

    #[test]
    fn test_prefix_search_basic() {
        let idx = sample_index();
        let results = idx.prefix_search("foo", 20);
        assert_eq!(results, vec!["foo", "foobar", "foobaz"]);
    }

    #[test]
    fn test_prefix_search_limit() {
        let idx = sample_index();
        let results = idx.prefix_search("foo", 1);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0], "foo");
    }

    #[test]
    fn test_prefix_search_no_match() {
        let idx = sample_index();
        let results = idx.prefix_search("zzz", 20);
        assert!(results.is_empty());
    }

    #[test]
    fn test_prefix_search_empty_prefix() {
        let idx = sample_index();
        let results = idx.prefix_search("", 100);
        assert_eq!(results.len(), 6);
    }

    // ── V1.5 round-trip serialize/deserialize ────────────────────────────────

    #[test]
    fn test_vocab_roundtrip() {
        let idx = sample_index();
        let mut buf: Vec<u8> = Vec::new();
        serialize_vocab(&idx, &mut buf).unwrap();
        let idx2 = deserialize_vocab(&buf);
        assert_eq!(idx.sorted_terms, idx2.sorted_terms);
    }
}
