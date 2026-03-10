//! Co-occurrence ontology — section-level and file-level term relationships
//!
//! Builds a weighted graph of term co-occurrences from indexed corpus files.
//! Two scopes: section-level (high weight α=3) and file-level (lower weight β=1).

use std::collections::HashMap;
use std::io::{Read, Write};

use serde::{Deserialize, Serialize};

use super::sections::{detect_sections, extract_tokens};

/// A single co-occurrence edge between two terms
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CooccurrenceEdge {
    pub term_a: String,
    pub term_b: String,
    pub weight: f32,
    pub section_count: u32,
    pub file_count: u32,
}

/// The full ontology: a list of co-occurrence edges
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Ontology {
    pub edges: Vec<CooccurrenceEdge>,
    /// Number of files processed
    pub file_count: usize,
    /// Number of sections detected
    pub section_count: usize,
}

/// Configuration for ontology building
pub struct OntologyConfig {
    /// Weight for section-level co-occurrence
    pub alpha: f32,
    /// Weight for file-level co-occurrence
    pub beta: f32,
    /// Maximum edges to keep per term
    pub top_k: usize,
    /// Minimum weight to keep an edge
    pub min_weight: f32,
}

impl Default for OntologyConfig {
    fn default() -> Self {
        Self {
            alpha: 3.0,
            beta: 1.0,
            top_k: 10,
            min_weight: 1.0,
        }
    }
}

/// Build an ontology from corpus files.
///
/// Each file is represented as (path, content) pairs.
pub fn build_ontology(
    files: &[(&str, &str)],
    config: &OntologyConfig,
) -> Ontology {
    // Accumulate pair counts at both scopes
    let mut section_pairs: HashMap<(String, String), u32> = HashMap::new();
    let mut file_pairs: HashMap<(String, String), u32> = HashMap::new();
    let mut total_sections = 0;

    for &(path, content) in files {
        // Detect sections within the file
        let sections = detect_sections(content, path);
        total_sections += sections.len();

        // Section-level co-occurrence
        for section in &sections {
            let tokens = &section.tokens;
            for i in 0..tokens.len() {
                for j in (i + 1)..tokens.len() {
                    let key = ordered_pair(&tokens[i], &tokens[j]);
                    *section_pairs.entry(key).or_insert(0) += 1;
                }
            }
        }

        // File-level co-occurrence (all tokens in the file)
        let file_tokens = extract_tokens(content);
        for i in 0..file_tokens.len() {
            for j in (i + 1)..file_tokens.len() {
                let key = ordered_pair(&file_tokens[i], &file_tokens[j]);
                *file_pairs.entry(key).or_insert(0) += 1;
            }
        }
    }

    // Combine into weighted edges
    let mut all_pairs: HashMap<(String, String), (u32, u32)> = HashMap::new();

    for (pair, count) in &section_pairs {
        all_pairs.entry(pair.clone()).or_insert((0, 0)).0 = *count;
    }
    for (pair, count) in &file_pairs {
        all_pairs.entry(pair.clone()).or_insert((0, 0)).1 = *count;
    }

    let mut edges: Vec<CooccurrenceEdge> = all_pairs
        .into_iter()
        .map(|((a, b), (sec, file))| {
            let weight = config.alpha * sec as f32 + config.beta * file as f32;
            CooccurrenceEdge {
                term_a: a,
                term_b: b,
                weight,
                section_count: sec,
                file_count: file,
            }
        })
        .filter(|e| e.weight >= config.min_weight)
        .collect();

    // Top-k pruning: for each term, keep only its top_k strongest edges
    let pruned = top_k_prune(&mut edges, config.top_k);

    Ontology {
        edges: pruned,
        file_count: files.len(),
        section_count: total_sections,
    }
}

/// Ensure pair ordering is deterministic (alphabetical)
fn ordered_pair(a: &str, b: &str) -> (String, String) {
    if a <= b {
        (a.to_string(), b.to_string())
    } else {
        (b.to_string(), a.to_string())
    }
}

/// Keep only top_k edges per term
fn top_k_prune(edges: &mut [CooccurrenceEdge], top_k: usize) -> Vec<CooccurrenceEdge> {
    // Sort by weight descending
    edges.sort_by(|a, b| b.weight.partial_cmp(&a.weight).unwrap_or(std::cmp::Ordering::Equal));

    // Count per term
    let mut term_counts: HashMap<&str, usize> = HashMap::new();
    let mut kept = Vec::new();

    for edge in edges.iter() {
        let count_a = term_counts.get(edge.term_a.as_str()).copied().unwrap_or(0);
        let count_b = term_counts.get(edge.term_b.as_str()).copied().unwrap_or(0);

        if count_a < top_k && count_b < top_k {
            *term_counts.entry(&edge.term_a).or_insert(0) += 1;
            *term_counts.entry(&edge.term_b).or_insert(0) += 1;
            kept.push(edge.clone());
        }
    }

    kept
}

impl Ontology {
    /// Look up the top related terms for a given term
    pub fn related(&self, term: &str, limit: usize) -> Vec<(&CooccurrenceEdge, &str)> {
        let lower = term.to_lowercase();
        let mut matches: Vec<(&CooccurrenceEdge, &str)> = self
            .edges
            .iter()
            .filter_map(|e| {
                if e.term_a == lower {
                    Some((e, e.term_b.as_str()))
                } else if e.term_b == lower {
                    Some((e, e.term_a.as_str()))
                } else {
                    None
                }
            })
            .collect();

        matches.sort_by(|a, b| b.0.weight.partial_cmp(&a.0.weight).unwrap_or(std::cmp::Ordering::Equal));
        matches.truncate(limit);
        matches
    }

    /// Serialize to bytes (JSON)
    pub fn to_bytes(&self) -> Result<Vec<u8>, String> {
        serde_json::to_vec(self).map_err(|e| e.to_string())
    }

    /// Deserialize from bytes (JSON)
    pub fn from_bytes(data: &[u8]) -> Result<Self, String> {
        serde_json::from_slice(data).map_err(|e| e.to_string())
    }

    /// Write to a writer
    pub fn write_to<W: Write>(&self, writer: &mut W) -> Result<(), String> {
        let data = self.to_bytes()?;
        let len = data.len() as u64;
        writer.write_all(&len.to_le_bytes()).map_err(|e| e.to_string())?;
        writer.write_all(&data).map_err(|e| e.to_string())?;
        Ok(())
    }

    /// Read from a reader
    pub fn read_from<R: Read>(reader: &mut R) -> Result<Self, String> {
        let mut len_buf = [0u8; 8];
        reader.read_exact(&mut len_buf).map_err(|e| e.to_string())?;
        let len = u64::from_le_bytes(len_buf) as usize;

        let mut data = vec![0u8; len];
        reader.read_exact(&mut data).map_err(|e| e.to_string())?;
        Self::from_bytes(&data)
    }

    /// Summary stats
    pub fn stats(&self) -> String {
        let unique_terms: std::collections::HashSet<&str> = self
            .edges
            .iter()
            .flat_map(|e| [e.term_a.as_str(), e.term_b.as_str()])
            .collect();

        format!(
            "{} edges, {} unique terms, {} files, {} sections",
            self.edges.len(),
            unique_terms.len(),
            self.file_count,
            self.section_count,
        )
    }
}

/// Tag for the ONTO section in .idx files
pub const ONTO_TAG: &[u8; 4] = b"ONTO";

/// A single co-occurrence triple for binary serialization
#[derive(Debug, Clone)]
pub struct OntologyTriple {
    pub term_a: String,
    pub term_b: String,
    pub weight: f32,
    /// 0 = file-level only, 1 = has section-level co-occurrences
    pub scope: u8,
}

/// Binary representation of the ontology for embedding in .idx
#[derive(Debug, Clone)]
pub struct OntologySection {
    pub triples: Vec<OntologyTriple>,
}

impl From<&Ontology> for OntologySection {
    fn from(ont: &Ontology) -> Self {
        let mut triples: Vec<OntologyTriple> = ont
            .edges
            .iter()
            .map(|e| OntologyTriple {
                term_a: e.term_a.clone(),
                term_b: e.term_b.clone(),
                weight: e.weight,
                scope: if e.section_count > 0 { 1 } else { 0 },
            })
            .collect();
        // Sort deterministically by (term_a, term_b)
        triples.sort_by(|a, b| a.term_a.cmp(&b.term_a).then(a.term_b.cmp(&b.term_b)));
        Self { triples }
    }
}

impl OntologySection {
    /// Serialize to binary: num_triples(u64) + per-triple records
    /// Each triple: term_a_len(u16) + term_a_bytes + term_b_len(u16) + term_b_bytes + weight(f32 LE) + scope(u8)
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::new();
        buf.extend_from_slice(&(self.triples.len() as u64).to_le_bytes());
        for t in &self.triples {
            let ta = t.term_a.as_bytes();
            buf.extend_from_slice(&(ta.len() as u16).to_le_bytes());
            buf.extend_from_slice(ta);
            let tb = t.term_b.as_bytes();
            buf.extend_from_slice(&(tb.len() as u16).to_le_bytes());
            buf.extend_from_slice(tb);
            buf.extend_from_slice(&t.weight.to_le_bytes());
            buf.push(t.scope);
        }
        buf
    }

    /// Deserialize from binary bytes
    pub fn from_bytes(data: &[u8]) -> Result<Self, String> {
        let mut pos = 0;

        macro_rules! need {
            ($n:expr) => {
                if pos + $n > data.len() {
                    return Err(format!("truncated at offset {}", pos));
                }
            };
        }

        need!(8);
        let n = u64::from_le_bytes(data[pos..pos + 8].try_into().unwrap()) as usize;
        pos += 8;

        // Defensive bound: each triple needs at least 2 + 0 + 2 + 0 + 4 + 1 = 9 bytes.
        // Reject impossible counts early to avoid oversized allocations on malformed input.
        let min_bytes_per_triple = 9usize;
        if n > (data.len().saturating_sub(pos)) / min_bytes_per_triple {
            return Err(format!(
                "invalid ontology triple count {} for payload size {}",
                n,
                data.len()
            ));
        }

        let mut triples = Vec::with_capacity(n);
        for _ in 0..n {
            need!(2);
            let ta_len = u16::from_le_bytes(data[pos..pos + 2].try_into().unwrap()) as usize;
            pos += 2;
            need!(ta_len);
            let term_a = String::from_utf8(data[pos..pos + ta_len].to_vec())
                .map_err(|e| e.to_string())?;
            pos += ta_len;

            need!(2);
            let tb_len = u16::from_le_bytes(data[pos..pos + 2].try_into().unwrap()) as usize;
            pos += 2;
            need!(tb_len);
            let term_b = String::from_utf8(data[pos..pos + tb_len].to_vec())
                .map_err(|e| e.to_string())?;
            pos += tb_len;

            need!(5); // 4 (f32) + 1 (u8)
            let weight = f32::from_le_bytes(data[pos..pos + 4].try_into().unwrap());
            pos += 4;
            let scope = data[pos];
            pos += 1;

            triples.push(OntologyTriple { term_a, term_b, weight, scope });
        }

        Ok(Self { triples })
    }

    /// Look up top related terms for a given term
    pub fn related(&self, term: &str, limit: usize) -> Vec<(&OntologyTriple, &str)> {
        let lower = term.to_lowercase();
        let mut matches: Vec<(&OntologyTriple, &str)> = self
            .triples
            .iter()
            .filter_map(|t| {
                if t.term_a == lower {
                    Some((t, t.term_b.as_str()))
                } else if t.term_b == lower {
                    Some((t, t.term_a.as_str()))
                } else {
                    None
                }
            })
            .collect();
        matches.sort_by(|a, b| {
            b.0.weight
                .partial_cmp(&a.0.weight)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        matches.truncate(limit);
        matches
    }

    /// Number of triples
    pub fn len(&self) -> usize {
        self.triples.len()
    }

    /// Whether empty
    pub fn is_empty(&self) -> bool {
        self.triples.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_simple_ontology() {
        let files = vec![
            ("test.rs", "fn handle_message(msg: &Message) {\n    deliver(msg);\n    acknowledge(msg);\n}\n\nfn process_queue(q: &Queue) {\n    let msg = q.pop();\n    deliver(msg);\n}\n"),
        ];
        let config = OntologyConfig::default();
        let ont = build_ontology(&files, &config);

        assert!(!ont.edges.is_empty());
        assert_eq!(ont.file_count, 1);

        // "deliver" and "msg" should co-occur
        let deliver_related = ont.related("deliver", 10);
        assert!(!deliver_related.is_empty());
    }

    #[test]
    fn test_section_vs_file_weighting() {
        // Two terms in same section should have higher weight than
        // two terms in different sections of the same file
        let files = vec![
            ("test.erl", "-module(test).\n\nstart(Config) ->\n    channel_open(Config),\n    deliver(Config).\n\n%%--- Other ---\n\nstop(Config) ->\n    cleanup(Config).\n"),
        ];
        let config = OntologyConfig {
            alpha: 3.0,
            beta: 1.0,
            top_k: 20,
            min_weight: 0.0, // keep all for testing
        };
        let ont = build_ontology(&files, &config);

        // channel_open and deliver are in same section → section_count > 0
        let edge = ont.edges.iter().find(|e| {
            (e.term_a.contains("channel") && e.term_b.contains("deliver"))
                || (e.term_a.contains("deliver") && e.term_b.contains("channel"))
        });
        if let Some(e) = edge {
            assert!(e.section_count > 0, "section co-occurrence expected");
        }
    }

    #[test]
    fn test_related_lookup() {
        let ont = Ontology {
            edges: vec![
                CooccurrenceEdge {
                    term_a: "ack".to_string(),
                    term_b: "deliver".to_string(),
                    weight: 15.0,
                    section_count: 3,
                    file_count: 5,
                },
                CooccurrenceEdge {
                    term_a: "ack".to_string(),
                    term_b: "message".to_string(),
                    weight: 10.0,
                    section_count: 2,
                    file_count: 4,
                },
            ],
            file_count: 10,
            section_count: 50,
        };

        let related = ont.related("ack", 10);
        assert_eq!(related.len(), 2);
        assert_eq!(related[0].1, "deliver"); // highest weight first
        assert_eq!(related[1].1, "message");
    }

    #[test]
    fn test_top_k_pruning() {
        let content = "fn a() { b; c; d; e; f; g; h; i; j; k; l; m; n; }\n\
                        fn p() { q; r; s; t; u; v; w; x; y; z; }\n";
        let files = vec![("big.rs", content)];
        let config = OntologyConfig {
            top_k: 5,
            ..OntologyConfig::default()
        };
        let ont = build_ontology(&files, &config);

        // Each term should have at most top_k edges
        let mut term_edge_counts: HashMap<&str, usize> = HashMap::new();
        for e in &ont.edges {
            *term_edge_counts.entry(&e.term_a).or_insert(0) += 1;
            *term_edge_counts.entry(&e.term_b).or_insert(0) += 1;
        }
        for (_, count) in &term_edge_counts {
            assert!(*count <= 5, "term has {} edges, expected <= 5", count);
        }
    }

    #[test]
    fn test_serialize_roundtrip() {
        let ont = Ontology {
            edges: vec![CooccurrenceEdge {
                term_a: "foo".to_string(),
                term_b: "bar".to_string(),
                weight: 7.0,
                section_count: 1,
                file_count: 2,
            }],
            file_count: 1,
            section_count: 3,
        };

        let bytes = ont.to_bytes().unwrap();
        let restored = Ontology::from_bytes(&bytes).unwrap();

        assert_eq!(restored.edges.len(), 1);
        assert_eq!(restored.edges[0].term_a, "foo");
        assert_eq!(restored.edges[0].weight, 7.0);
        assert_eq!(restored.file_count, 1);
        assert_eq!(restored.section_count, 3);
    }

    #[test]
    fn test_ontology_section_roundtrip() {
        let ont = Ontology {
            edges: vec![
                CooccurrenceEdge {
                    term_a: "ack".to_string(),
                    term_b: "deliver".to_string(),
                    weight: 9.0,
                    section_count: 2,
                    file_count: 3,
                },
                CooccurrenceEdge {
                    term_a: "channel".to_string(),
                    term_b: "queue".to_string(),
                    weight: 4.0,
                    section_count: 0,
                    file_count: 4,
                },
            ],
            file_count: 5,
            section_count: 20,
        };
        let section = OntologySection::from(&ont);
        assert_eq!(section.triples.len(), 2);
        // Sorted: ack < channel
        assert_eq!(section.triples[0].term_a, "ack");
        assert_eq!(section.triples[0].scope, 1); // section_count > 0
        assert_eq!(section.triples[1].scope, 0); // file-level only

        let bytes = section.to_bytes();
        let restored = OntologySection::from_bytes(&bytes).unwrap();
        assert_eq!(restored.triples.len(), 2);
        assert_eq!(restored.triples[0].term_a, "ack");
        assert_eq!(restored.triples[0].term_b, "deliver");
        assert!((restored.triples[0].weight - 9.0).abs() < 1e-6);
        assert_eq!(restored.triples[0].scope, 1);
        assert_eq!(restored.triples[1].term_a, "channel");
        assert_eq!(restored.triples[1].scope, 0);

        let related = restored.related("ack", 10);
        assert_eq!(related.len(), 1);
        assert_eq!(related[0].1, "deliver");
    }
}
