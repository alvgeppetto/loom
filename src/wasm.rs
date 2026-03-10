use crate::indexer::ontology::{Ontology, OntologySection, ONTO_TAG};
use fm_index::converter::RangeConverter;
use fm_index::suffix_array::SuffixOrderSampler;
use fm_index::{BackwardSearchIndex, FMIndex};
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

const TAG_END: [u8; 4] = [0u8; 4];
const DEFAULT_ONTOLOGY_HINTS: [&str; 6] = [
    "rabbit",
    "erlang",
    "gen_server",
    "handle_call",
    "handle_cast",
    "handle_info",
];
const STOPWORDS: [&str; 24] = [
    "how", "does", "where", "what", "which", "when", "why", "is", "are", "the", "a",
    "an", "to", "and", "or", "in", "on", "for", "with", "from", "of", "it", "this",
    "that",
];

#[derive(Deserialize, Clone)]
struct FileMetadata {
    path: String,
    start_offset: usize,
    end_offset: usize,
}

#[derive(Deserialize)]
struct IndexMetadata {
    files: Vec<FileMetadata>,
    corpus_size: usize,
}

#[derive(Serialize)]
struct SearchHit {
    pattern: String,
    path: String,
    line: usize,
    column: usize,
    offset: u64,
    context: String,
}

#[derive(Serialize, Debug, Clone)]
struct OntologyRelatedHit {
    term: String,
    weight: f32,
    scope: u8,
}

#[derive(Serialize)]
struct OntologyEnrichmentRun {
    #[serde(rename = "type")]
    run_type: String,
    #[serde(rename = "rowCount")]
    row_count: usize,
    #[serde(rename = "translated_sparql")]
    translated_sparql: Option<String>,
    #[serde(rename = "translation_mode")]
    translation_mode: String,
    #[serde(rename = "timing_ms")]
    timing_ms: Option<u64>,
    error: Option<String>,
    seed: Option<String>,
}

#[derive(Serialize)]
struct OntologyEnrichmentResult {
    runs: Vec<OntologyEnrichmentRun>,
    terms: Vec<String>,
    summary: String,
}

#[wasm_bindgen]
pub struct LoomWasmIndex {
    index: FMIndex<
        u8,
        RangeConverter<u8>,
        fm_index::suffix_array::SuffixOrderSampledArray,
    >,
    metadata: IndexMetadata,
    corpus_text: Vec<u8>,
    ontology: Option<OntologySection>,
    /// Raw MGPT section bytes — parsed to keep index-format compat, unused post-pivot.
    #[allow(dead_code)]
    mgpt_bytes: Option<Vec<u8>>,
}

#[wasm_bindgen]
impl LoomWasmIndex {
    #[wasm_bindgen(js_name = fromSerialized)]
    pub fn from_serialized(bytes: Vec<u8>) -> Result<LoomWasmIndex, JsValue> {
        if bytes.len() < 16 {
            return Err(JsValue::from_str("Invalid index: too small"));
        }

        let mut cursor = 0usize;

        let metadata_len = read_u64(&bytes, &mut cursor)? as usize;
        let metadata_end = cursor.saturating_add(metadata_len);
        if metadata_end > bytes.len() {
            return Err(JsValue::from_str("Invalid index: metadata truncated"));
        }

        let metadata: IndexMetadata = serde_json::from_slice(&bytes[cursor..metadata_end])
            .map_err(|e| JsValue::from_str(&format!("Invalid metadata JSON: {e}")))?;
        cursor = metadata_end;

        let corpus_len = read_u64(&bytes, &mut cursor)? as usize;
        let corpus_end = cursor.saturating_add(corpus_len);
        if corpus_end > bytes.len() {
            return Err(JsValue::from_str("Invalid index: corpus truncated"));
        }

        // Remember corpus location but do NOT copy yet — keep peak memory low.
        let corpus_start = cursor;
        cursor = corpus_end;

        let sections = parse_sections(&bytes, &mut cursor)?;

        // PCv1 magic prefix for postcard-encoded FMIX payloads.
        const FMIX_MAGIC: &[u8; 4] = b"PCv1";

        // Deserialize pre-built FM-index when FMIX section is present;
        // otherwise fall back to rebuilding from corpus (small indices only).
        let (index, corpus_text) = if let Some((start, end)) = sections.fmix_range {
            let fmix_data = &bytes[start..end];
            if fmix_data.len() > 4 && &fmix_data[..4] == FMIX_MAGIC {
                let fm: FMIndex<u8, RangeConverter<u8>, fm_index::suffix_array::SuffixOrderSampledArray> =
                    postcard::from_bytes(&fmix_data[4..])
                        .map_err(|e| JsValue::from_str(&format!("FMIX deserialize failed: {e}")))?;
                // Copy corpus AFTER deserializing FM-index, then drop input bytes.
                let ct = bytes[corpus_start..corpus_end].to_vec();
                drop(bytes);
                (fm, ct)
            } else {
                // Stale bincode payload — rebuild from corpus.
                let ct = bytes[corpus_start..corpus_end].to_vec();
                drop(bytes);
                let converter = RangeConverter::new(1u8, 255u8);
                let sampler = SuffixOrderSampler::new().level(2);
                let fm = FMIndex::new(ct.clone(), converter, sampler);
                (fm, ct)
            }
        } else {
            let ct = bytes[corpus_start..corpus_end].to_vec();
            drop(bytes);
            let converter = RangeConverter::new(1u8, 255u8);
            let sampler = SuffixOrderSampler::new().level(2);
            let fm = FMIndex::new(ct.clone(), converter, sampler);
            (fm, ct)
        };

        Ok(LoomWasmIndex {
            index,
            metadata,
            corpus_text,
            ontology: sections.ontology,
            mgpt_bytes: sections.mgpt_bytes,
        })
    }

    #[wasm_bindgen(js_name = fileCount)]
    pub fn file_count(&self) -> usize {
        self.metadata.files.len()
    }

    #[wasm_bindgen(js_name = corpusSize)]
    pub fn corpus_size(&self) -> usize {
        self.metadata.corpus_size
    }

    pub fn count(&self, pattern: &str) -> u64 {
        self.index.search_backward(pattern.as_bytes()).count()
    }

    pub fn search(
        &self,
        pattern: &str,
        limit: usize,
        context_chars: usize,
    ) -> Result<JsValue, JsValue> {
        let mut hits = Vec::new();
        let positions = self.index.search_backward(pattern.as_bytes()).locate();

        for offset in positions.into_iter().take(limit) {
            if let Some((path, line, column)) = self.line_col_at_offset(offset) {
                let context = self.context_at_offset(offset, context_chars, context_chars);
                hits.push(SearchHit {
                    pattern: pattern.to_string(),
                    path,
                    line,
                    column,
                    offset,
                    context,
                });
            }
        }

        serde_wasm_bindgen::to_value(&hits)
            .map_err(|e| JsValue::from_str(&format!("Failed to serialize hits: {e}")))
    }

    #[wasm_bindgen(js_name = relatedTerms)]
    pub fn related_terms(&self, term: &str, limit: usize) -> Result<JsValue, JsValue> {
        let hits = self.related_terms_internal(term, limit);
        serde_wasm_bindgen::to_value(&hits)
            .map_err(|e| JsValue::from_str(&format!("Failed to serialize related terms: {e}")))
    }

    #[wasm_bindgen(js_name = enrichTermsWithOntology)]
    pub fn enrich_terms_with_ontology(
        &self,
        question: &str,
        max_terms: usize,
        related_per_term: usize,
    ) -> Result<JsValue, JsValue> {
        let result = self.enrich_terms_with_ontology_internal(question, max_terms, related_per_term);
        serde_wasm_bindgen::to_value(&result)
            .map_err(|e| JsValue::from_str(&format!("Failed to serialize ontology enrichment: {e}")))
    }

    /// Enrich search terms using a custom set of domain hint strings.
    ///
    /// Identical to `enrichTermsWithOntology` but replaces the built-in
    /// `DEFAULT_ONTOLOGY_HINTS` with the caller-supplied `hints_json` array.
    /// This is the WASM counterpart of the Python `--prompt` CLI flag: pass
    /// domain-specific seed terms (e.g. DNA motifs for a genome task) so the
    /// ontology expansion starts from the right vocabulary.
    ///
    /// # Parameters
    /// - `hints_json`: JSON array of strings, e.g. `["TATAAA","ATG","AATAAA"]`.
    ///   Pass `"[]"` or `null` to use no domain hints (question tokens only).
    #[wasm_bindgen(js_name = enrichTermsWithCustomHints)]
    pub fn enrich_terms_with_custom_hints(
        &self,
        question: &str,
        hints_json: &str,
        max_terms: usize,
        related_per_term: usize,
    ) -> Result<JsValue, JsValue> {
        let hints: Vec<String> = serde_json::from_str(hints_json)
            .map_err(|e| JsValue::from_str(&format!("Invalid hints JSON: {e}")))?;
        let result = self.enrich_terms_with_custom_hints_internal(
            question,
            &hints,
            max_terms,
            related_per_term,
        );
        serde_wasm_bindgen::to_value(&result)
            .map_err(|e| JsValue::from_str(&format!("Failed to serialize enrichment: {e}")))
    }

    /// Search for multiple patterns at once and return all hits combined.
    ///
    /// This is the multi-term search endpoint for the JS/WASM frontend.
    /// Pass a JSON array of search terms (from any source — LLM expansion,
    /// custom prompt, or direct user input) and receive a flat array of hits.
    ///
    /// # Parameters
    /// - `terms_json`: JSON array of pattern strings, e.g. `["TATAAA","ATG"]`.
    /// - `limit_per_term`: Maximum hits to return per term.
    /// - `context_chars`: Characters of context on each side of a hit.
    #[wasm_bindgen(js_name = searchMulti)]
    pub fn search_multi(
        &self,
        terms_json: &str,
        limit_per_term: usize,
        context_chars: usize,
    ) -> Result<JsValue, JsValue> {
        let terms: Vec<String> = serde_json::from_str(terms_json)
            .map_err(|e| JsValue::from_str(&format!("Invalid terms JSON: {e}")))?;
        let mut all_hits: Vec<SearchHit> = Vec::new();
        for term in &terms {
            let positions = self.index.search_backward(term.as_bytes()).locate();
            for offset in positions.into_iter().take(limit_per_term) {
                if let Some((path, line, column)) = self.line_col_at_offset(offset) {
                    let context = self.context_at_offset(offset, context_chars, context_chars);
                    all_hits.push(SearchHit {
                        pattern: term.clone(),
                        path,
                        line,
                        column,
                        offset,
                        context,
                    });
                }
            }
        }
        serde_wasm_bindgen::to_value(&all_hits)
            .map_err(|e| JsValue::from_str(&format!("Failed to serialize hits: {e}")))
    }
}

impl LoomWasmIndex {
    fn related_terms_internal(&self, term: &str, limit: usize) -> Vec<OntologyRelatedHit> {
        let Some(ontology) = self.ontology.as_ref() else {
            return Vec::new();
        };

        ontology
            .related(&term.to_lowercase(), limit)
            .into_iter()
            .map(|(triple, related)| OntologyRelatedHit {
                term: related.to_string(),
                weight: triple.weight,
                scope: triple.scope,
            })
            .collect()
    }

    fn enrich_terms_with_ontology_internal(
        &self,
        question: &str,
        max_terms: usize,
        related_per_term: usize,
    ) -> OntologyEnrichmentResult {
        let seeds = build_seed_terms(question);
        if seeds.is_empty() {
            return OntologyEnrichmentResult {
                runs: vec![OntologyEnrichmentRun {
                    run_type: "embedded_ontology".to_string(),
                    row_count: 0,
                    translated_sparql: None,
                    translation_mode: "embedded-index".to_string(),
                    timing_ms: None,
                    error: None,
                    seed: None,
                }],
                terms: Vec::new(),
                summary: String::new(),
            };
        }

        let mut runs = Vec::new();
        let mut ranked_terms: Vec<(String, f32)> = Vec::new();

        for seed in &seeds {
            ranked_terms.push((seed.clone(), 10_000.0));
            let related = self.related_terms_internal(seed, related_per_term);
            runs.push(OntologyEnrichmentRun {
                run_type: "embedded_related".to_string(),
                row_count: related.len(),
                translated_sparql: None,
                translation_mode: "embedded-index".to_string(),
                timing_ms: None,
                error: None,
                seed: Some(seed.clone()),
            });
            for related_hit in related {
                ranked_terms.push((related_hit.term, related_hit.weight));
            }
        }

        if self.ontology.is_none() {
            runs.push(OntologyEnrichmentRun {
                run_type: "embedded_ontology_unavailable".to_string(),
                row_count: 0,
                translated_sparql: None,
                translation_mode: "embedded-index".to_string(),
                timing_ms: None,
                error: Some("ONTO section not present in index".to_string()),
                seed: None,
            });
        }

        ranked_terms.sort_by(|a, b| {
            b.1.partial_cmp(&a.1)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.0.cmp(&b.0))
        });

        let mut seen = std::collections::HashSet::new();
        let mut terms = Vec::new();
        for (term, _) in ranked_terms {
            let key = term.to_lowercase();
            if key.len() < 3 || !seen.insert(key) {
                continue;
            }
            terms.push(term);
            if terms.len() >= max_terms {
                break;
            }
        }

        let summary = if terms.is_empty() {
            String::new()
        } else {
            let preview = terms.iter().take(12).cloned().collect::<Vec<_>>().join(", ");
            format!("Embedded ontology terms: {preview}")
        };

        OntologyEnrichmentResult { runs, terms, summary }
    }

    fn enrich_terms_with_custom_hints_internal(
        &self,
        question: &str,
        hints: &[String],
        max_terms: usize,
        related_per_term: usize,
    ) -> OntologyEnrichmentResult {
        // Build seed terms from the question, but swap in caller-supplied hints
        // instead of DEFAULT_ONTOLOGY_HINTS.
        let stopwords: std::collections::HashSet<&'static str> = STOPWORDS.into_iter().collect();
        let mut seeds: Vec<String> = question
            .split(|ch: char| !(ch.is_ascii_alphanumeric() || ch == '_'))
            .map(|token| token.trim().to_lowercase())
            .filter(|token| token.len() > 2 && !stopwords.contains(token.as_str()))
            .collect();
        for hint in hints {
            seeds.push(hint.to_lowercase());
        }
        let mut seen_s = std::collections::HashSet::new();
        seeds.retain(|t| seen_s.insert(t.clone()));

        if seeds.is_empty() {
            return OntologyEnrichmentResult {
                runs: vec![OntologyEnrichmentRun {
                    run_type: "custom_hints".to_string(),
                    row_count: 0,
                    translated_sparql: None,
                    translation_mode: "custom-hints".to_string(),
                    timing_ms: None,
                    error: None,
                    seed: None,
                }],
                terms: Vec::new(),
                summary: String::new(),
            };
        }

        let mut runs = Vec::new();
        let mut ranked_terms: Vec<(String, f32)> = Vec::new();

        for seed in &seeds {
            ranked_terms.push((seed.clone(), 10_000.0));
            let related = self.related_terms_internal(seed, related_per_term);
            runs.push(OntologyEnrichmentRun {
                run_type: "custom_hints_related".to_string(),
                row_count: related.len(),
                translated_sparql: None,
                translation_mode: "custom-hints".to_string(),
                timing_ms: None,
                error: None,
                seed: Some(seed.clone()),
            });
            for hit in related {
                ranked_terms.push((hit.term, hit.weight));
            }
        }

        ranked_terms.sort_by(|a, b| {
            b.1.partial_cmp(&a.1)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.0.cmp(&b.0))
        });

        let mut seen = std::collections::HashSet::new();
        let mut terms = Vec::new();
        for (term, _) in ranked_terms {
            let key = term.to_lowercase();
            if key.len() < 2 || !seen.insert(key) {
                continue;
            }
            terms.push(term);
            if terms.len() >= max_terms {
                break;
            }
        }

        let summary = if terms.is_empty() {
            String::new()
        } else {
            let preview = terms.iter().take(12).cloned().collect::<Vec<_>>().join(", ");
            format!("Custom-hints terms: {preview}")
        };

        OntologyEnrichmentResult { runs, terms, summary }
    }

    fn file_at_offset(&self, offset: u64) -> Option<&FileMetadata> {
        self.metadata
            .files
            .iter()
            .find(|file| offset >= file.start_offset as u64 && offset < file.end_offset as u64)
    }

    fn line_col_at_offset(&self, offset: u64) -> Option<(String, usize, usize)> {
        let file = self.file_at_offset(offset)?;
        let offset = offset as usize;
        let content_start = file.start_offset + file.path.len() + 2;

        if offset < content_start || offset >= self.corpus_text.len() {
            return None;
        }

        let relative_offset = offset - content_start;
        let content = &self.corpus_text[content_start..offset];
        let line = content.iter().filter(|&&byte| byte == b'\n').count() + 1;

        let last_newline = content.iter().rposition(|&byte| byte == b'\n');
        let column = match last_newline {
            Some(position) => relative_offset - position,
            None => relative_offset + 1,
        };

        Some((file.path.clone(), line, column))
    }

    fn context_at_offset(&self, offset: u64, chars_before: usize, chars_after: usize) -> String {
        let offset = offset as usize;
        let start = offset.saturating_sub(chars_before);
        let end = (offset + chars_after).min(self.corpus_text.len());
        String::from_utf8_lossy(&self.corpus_text[start..end]).to_string()
    }
}

const TAG_MGPT: [u8; 4] = *b"MGPT";
const TAG_FMIX: [u8; 4] = *b"FMIX";

/// Parsed tagged sections from the .idx byte stream.
struct ParsedSections {
    ontology: Option<OntologySection>,
    mgpt_bytes: Option<Vec<u8>>,
    /// Byte range (start, end) of the pre-serialized FM-index in the input buffer.
    fmix_range: Option<(usize, usize)>,
}

fn parse_sections(bytes: &[u8], cursor: &mut usize) -> Result<ParsedSections, JsValue> {
    let mut ontology: Option<OntologySection> = None;
    let mut mgpt_bytes: Option<Vec<u8>> = None;
    let mut fmix_range: Option<(usize, usize)> = None;
    while *cursor < bytes.len() {
        let end_tag = cursor.saturating_add(4);
        if end_tag > bytes.len() {
            return Err(JsValue::from_str("Invalid index: truncated section tag"));
        }

        let mut tag = [0u8; 4];
        tag.copy_from_slice(&bytes[*cursor..end_tag]);
        *cursor = end_tag;

        if tag == TAG_END {
            break;
        }

        let len = read_u64(bytes, cursor)? as usize;
        let end = cursor.saturating_add(len);
        if end > bytes.len() {
            return Err(JsValue::from_str("Invalid index: section truncated"));
        }

        if tag == *ONTO_TAG {
            let onto_bytes = &bytes[*cursor..end];

            // Prefer compact binary ONTO format; fallback to legacy JSON ontology payloads.
            // Never fail index boot solely due to ONTO parse issues.
            if let Ok(parsed) = OntologySection::from_bytes(onto_bytes) {
                ontology = Some(parsed);
            } else if let Ok(legacy) = Ontology::from_bytes(onto_bytes) {
                ontology = Some(OntologySection::from(&legacy));
            }
        } else if tag == TAG_MGPT {
            mgpt_bytes = Some(bytes[*cursor..end].to_vec());
        } else if tag == TAG_FMIX {
            // Store range only — avoid copying potentially-GB blob.
            fmix_range = Some((*cursor, end));
        }

        *cursor = end;
    }

    Ok(ParsedSections { ontology, mgpt_bytes, fmix_range })
}

fn build_seed_terms(question: &str) -> Vec<String> {
    let stopwords: std::collections::HashSet<&'static str> = STOPWORDS.into_iter().collect();
    let mut terms: Vec<String> = question
        .split(|ch: char| !(ch.is_ascii_alphanumeric() || ch == '_'))
        .map(|token| token.trim().to_lowercase())
        .filter(|token| token.len() > 2 && !stopwords.contains(token.as_str()))
        .collect();

    for hint in DEFAULT_ONTOLOGY_HINTS {
        terms.push(hint.to_string());
    }

    let mut seen = std::collections::HashSet::new();
    let mut deduped = Vec::new();
    for term in terms {
        if seen.insert(term.clone()) {
            deduped.push(term);
        }
    }
    deduped
}

fn read_u64(bytes: &[u8], cursor: &mut usize) -> Result<u64, JsValue> {
    let end = cursor.saturating_add(8);
    if end > bytes.len() {
        return Err(JsValue::from_str("Invalid index: unexpected EOF"));
    }

    let mut buffer = [0u8; 8];
    buffer.copy_from_slice(&bytes[*cursor..end]);
    *cursor = end;
    Ok(u64::from_le_bytes(buffer))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::indexer::ontology::OntologyTriple;
    use serde_json::json;

    fn make_serialized_index(include_ontology: bool) -> Vec<u8> {
        let corpus_text = b"\x01test.erl\x01-module(test).\nhandle_call(A,B,C)->ok.\nqueue_durable_test() -> durable.\n".to_vec();

        let metadata = json!({
            "files": [
                {
                    "path": "test.erl",
                    "start_offset": 0,
                    "end_offset": corpus_text.len()
                }
            ],
            "corpus_size": corpus_text.len()
        });
        let metadata_bytes = serde_json::to_vec(&metadata).unwrap();

        let mut bytes = Vec::new();
        bytes.extend_from_slice(&(metadata_bytes.len() as u64).to_le_bytes());
        bytes.extend_from_slice(&metadata_bytes);
        bytes.extend_from_slice(&(corpus_text.len() as u64).to_le_bytes());
        bytes.extend_from_slice(&corpus_text);

        if include_ontology {
            let section = OntologySection {
                triples: vec![
                    OntologyTriple {
                        term_a: "queue".to_string(),
                        term_b: "durable".to_string(),
                        weight: 1.0,
                        scope: 1,
                    },
                    OntologyTriple {
                        term_a: "queue".to_string(),
                        term_b: "rabbit_amqqueue".to_string(),
                        weight: 0.8,
                        scope: 1,
                    },
                    OntologyTriple {
                        term_a: "handle_call".to_string(),
                        term_b: "gen_server".to_string(),
                        weight: 0.7,
                        scope: 1,
                    },
                ],
            };

            let onto_bytes = section.to_bytes();
            bytes.extend_from_slice(ONTO_TAG);
            bytes.extend_from_slice(&(onto_bytes.len() as u64).to_le_bytes());
            bytes.extend_from_slice(&onto_bytes);
        }

        bytes.extend_from_slice(&TAG_END);
        bytes
    }

    fn make_serialized_index_with_malformed_onto() -> Vec<u8> {
        let corpus_text = b"\x01test.erl\x01-module(test).\nhandle_call(A,B,C)->ok.\n".to_vec();

        let metadata = json!({
            "files": [
                {
                    "path": "test.erl",
                    "start_offset": 0,
                    "end_offset": corpus_text.len()
                }
            ],
            "corpus_size": corpus_text.len()
        });
        let metadata_bytes = serde_json::to_vec(&metadata).unwrap();

        let mut bytes = Vec::new();
        bytes.extend_from_slice(&(metadata_bytes.len() as u64).to_le_bytes());
        bytes.extend_from_slice(&metadata_bytes);
        bytes.extend_from_slice(&(corpus_text.len() as u64).to_le_bytes());
        bytes.extend_from_slice(&corpus_text);

        // Corrupt ONTO payload that used to risk allocator trap in binary parser.
        let malformed = b"not-a-valid-onto!";
        bytes.extend_from_slice(ONTO_TAG);
        bytes.extend_from_slice(&(malformed.len() as u64).to_le_bytes());
        bytes.extend_from_slice(malformed);

        bytes.extend_from_slice(&TAG_END);
        bytes
    }

    #[test]
    fn parses_ontology_and_returns_related_terms() {
        let bytes = make_serialized_index(true);
        let index = LoomWasmIndex::from_serialized(bytes).expect("index should parse");
        let related = index.related_terms_internal("queue", 5);
        let terms: Vec<String> = related.into_iter().map(|r| r.term).collect();
        assert!(terms.contains(&"durable".to_string()));
        assert!(terms.contains(&"rabbit_amqqueue".to_string()));
    }

    #[test]
    fn enrichment_uses_embedded_ontology_terms() {
        let bytes = make_serialized_index(true);
        let index = LoomWasmIndex::from_serialized(bytes).expect("index should parse");
        let enriched = index.enrich_terms_with_ontology_internal(
            "How does queue durability work?",
            12,
            4,
        );

        assert!(enriched.terms.iter().any(|t| t == "queue"));
        assert!(enriched.terms.iter().any(|t| t == "durable"));
        assert!(enriched
            .runs
            .iter()
            .any(|r| r.run_type == "embedded_related" && r.row_count > 0));
    }

    #[test]
    fn enrichment_falls_back_when_ontology_missing() {
        let bytes = make_serialized_index(false);
        let index = LoomWasmIndex::from_serialized(bytes).expect("index should parse");
        let enriched = index.enrich_terms_with_ontology_internal("queue durability", 8, 4);

        assert!(enriched.terms.iter().any(|t| t == "queue"));
        assert!(enriched
            .runs
            .iter()
            .any(|r| r.run_type == "embedded_ontology_unavailable"));
    }

    #[test]
    fn malformed_ontology_payload_does_not_break_boot() {
        let bytes = make_serialized_index_with_malformed_onto();
        let index = LoomWasmIndex::from_serialized(bytes).expect("index should still parse");
        assert_eq!(index.file_count(), 1);
        assert_eq!(index.related_terms_internal("queue", 5).len(), 0);
    }
}
