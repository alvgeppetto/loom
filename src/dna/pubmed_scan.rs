//! PubMed Literature Scanner v4 — Ontology-integrated, multi-strategy, Rust-native.
//!
//! Replaces Python v2/v3 scanners with a Rust CLI that:
//!   1. Loads gene synonyms from NCBI ontology annotations (ontology-enrichment.json)
//!   2. Runs multi-strategy PubMed + Europe PMC queries per (pathogen × gene/application)
//!   3. Performs abstract corpus scanning for refutation
//!   4. Assigns confidence scores: CONFIRMED / PROBABLE / UNCERTAIN / FALSE
//!
//! Ontology integration: for each pathogen, the scanner loads all NCBI gene annotations
//! from web/data/ontology-enrichment.json and auto-generates synonym expansions from
//! official gene symbols and names. Manual synonyms in PATHOGEN_CONFIG are merged as
//! overrides (manual takes priority, ontology fills gaps).

use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

// ── Configuration ────────────────────────────────────────────────────────────

/// Rate limits (seconds between requests)
const PUBMED_DELAY_WITH_KEY: f64 = 0.12;
const PUBMED_DELAY_NO_KEY: f64 = 0.38;
const EPMC_DELAY: f64 = 0.38;
const MAX_LANDSCAPE_FETCH: usize = 500;
const MAX_DETAIL_FETCH: usize = 5;
const MAX_RETRIES: u32 = 4;

const EUTILS_BASE: &str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils";
const EUROPEPMC_BASE: &str = "https://www.ebi.ac.uk/europepmc/webservices/rest";

const EXCLUDE_KEYS: &[&str] = &["human-grch38", "human", "refseq-viral"];

// ── Ontology types ───────────────────────────────────────────────────────────

#[derive(Debug, Deserialize)]
struct OntologyFile {
    pathogens: HashMap<String, OntologyPathogen>,
}

#[derive(Debug, Deserialize)]
struct OntologyPathogen {
    species: Option<String>,
    genes: Option<OntologyGenes>,
    disease_ontology: Option<DiseaseOntology>,
    taxonomy: Option<Taxonomy>,
}

#[derive(Debug, Deserialize)]
struct OntologyGenes {
    reference_accession: Option<String>,
    gene_count: Option<usize>,
    genes: Vec<OntologyGene>,
}

#[derive(Debug, Deserialize)]
struct OntologyGene {
    symbol: Option<String>,
    name: Option<String>,
    gene_id: Option<u64>,
    #[serde(rename = "type")]
    gene_type: Option<String>,
}

#[derive(Debug, Deserialize)]
struct DiseaseOntology {
    name: Option<String>,
    synonyms: Option<Vec<String>>,
}

#[derive(Debug, Deserialize)]
struct Taxonomy {
    scientific_name: Option<String>,
}

// ── Pathogen config ──────────────────────────────────────────────────────────

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PathogenConfig {
    search_name: String,
    alt_names: Vec<String>,
    mesh_term: Option<String>,
    genes: Vec<String>,
    /// Manual synonyms (override ontology)
    gene_synonyms: HashMap<String, Vec<String>>,
    applications: BTreeMap<String, Vec<String>>,
}

// ── Result types ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "UPPERCASE")]
pub enum Confidence {
    Confirmed,
    Probable,
    Uncertain,
    False,
}

impl std::fmt::Display for Confidence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Confidence::Confirmed => write!(f, "CONFIRMED"),
            Confidence::Probable => write!(f, "PROBABLE"),
            Confidence::Uncertain => write!(f, "UNCERTAIN"),
            Confidence::False => write!(f, "FALSE"),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct QueryCounts {
    pubmed_narrow: i64,
    pubmed_broad: i64,
    europepmc: i64,
    corpus_title_matches: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct TopPaper {
    pmid: String,
    title: String,
    year: String,
    journal: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ClaimResult {
    label: String,
    confidence: Confidence,
    counts: QueryCounts,
    queries: BTreeMap<String, String>,
    corpus_hits: Vec<String>,
    top_papers: Vec<TopPaper>,
    /// Synonyms actually used in queries (for transparency)
    synonyms_used: Vec<String>,
    /// Source of synonyms: "manual", "ontology", or "manual+ontology"
    synonym_source: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PathogenResult {
    pathogen_key: String,
    search_name: String,
    landscape_total: i64,
    corpus_fetched: usize,
    genes: BTreeMap<String, ClaimResult>,
    applications: BTreeMap<String, ClaimResult>,
    confirmed_gene_gaps: Vec<String>,
    confirmed_app_gaps: Vec<String>,
    uncertain: Vec<String>,
    false_claims: Vec<String>,
}

#[derive(Debug, Serialize, Deserialize)]
struct ScanOutput {
    generated: String,
    version: u32,
    methodology: String,
    ontology_source: String,
    summary: ScanSummary,
    confirmed_gene_gaps: Vec<GapEntry>,
    confirmed_app_gaps: Vec<AppGapEntry>,
    uncertain: Vec<ItemEntry>,
    false_claims: Vec<ItemEntry>,
    pathogens: BTreeMap<String, PathogenResult>,
}

#[derive(Debug, Serialize, Deserialize)]
struct ScanSummary {
    pathogens_scanned: usize,
    confirmed_gene_gaps: usize,
    confirmed_app_gaps: usize,
    confirmed_gaps: usize,
    uncertain_count: usize,
    false_count: usize,
    ontology_genes_loaded: usize,
    ontology_synonyms_added: usize,
}

#[derive(Debug, Serialize, Deserialize)]
struct GapEntry {
    pathogen: String,
    gene: String,
}

#[derive(Debug, Serialize, Deserialize)]
struct AppGapEntry {
    pathogen: String,
    application: String,
}

#[derive(Debug, Serialize, Deserialize)]
struct ItemEntry {
    pathogen: String,
    item: String,
}

// ── HTTP helpers ─────────────────────────────────────────────────────────────

#[cfg(not(target_arch = "wasm32"))]
fn http_get_json(url: &str, retries: u32) -> Result<serde_json::Value, String> {
    for attempt in 0..retries {
        match ureq::get(url).timeout(Duration::from_secs(25)).call() {
            Ok(resp) => {
                let body: serde_json::Value = resp
                    .into_json()
                    .map_err(|e| format!("JSON parse error: {}", e))?;
                return Ok(body);
            }
            Err(e) => {
                if attempt + 1 < retries {
                    let wait = Duration::from_secs_f64(2.0 * (attempt as f64 + 1.0));
                    eprintln!(
                        "  [retry {}/{}] {}, sleeping {:.1}s...",
                        attempt + 1,
                        retries - 1,
                        e,
                        wait.as_secs_f64()
                    );
                    std::thread::sleep(wait);
                } else {
                    return Err(format!(
                        "HTTP fetch failed after {} attempts: {}",
                        retries,
                        &url[..url.len().min(80).min({
                        let mut e = url.len().min(80);
                        while e > 0 && !url.is_char_boundary(e) { e -= 1; }
                        e
                    })]
                    ));
                }
            }
        }
    }
    Err("unreachable".into())
}

fn rate_limit(delay_secs: f64) {
    std::thread::sleep(Duration::from_secs_f64(delay_secs));
}

// ── PubMed + Europe PMC API ──────────────────────────────────────────────────

#[cfg(not(target_arch = "wasm32"))]
fn pubmed_search(query: &str, max_ids: usize, api_key: &str, delay: f64) -> (i64, Vec<String>) {
    let mut params = vec![
        ("db", "pubmed"),
        ("retmode", "json"),
    ];
    let max_str = max_ids.to_string();
    params.push(("retmax", &max_str));
    let key_owned: String;
    if !api_key.is_empty() {
        key_owned = api_key.to_string();
        params.push(("api_key", &key_owned));
    }
    // URL-encode the query
    let encoded_query = urlencoded(query);
    let param_str: String = params
        .iter()
        .map(|(k, v)| format!("{}={}", k, urlencoded(v)))
        .collect::<Vec<_>>()
        .join("&");
    let url = format!("{}/esearch.fcgi?term={}&{}", EUTILS_BASE, encoded_query, param_str);

    match http_get_json(&url, MAX_RETRIES) {
        Ok(data) => {
            let r = &data["esearchresult"];
            let count = r["count"]
                .as_str()
                .and_then(|s| s.parse::<i64>().ok())
                .unwrap_or(0);
            let ids: Vec<String> = r["idlist"]
                .as_array()
                .map(|arr| {
                    arr.iter()
                        .filter_map(|v| v.as_str().map(String::from))
                        .collect()
                })
                .unwrap_or_default();
            rate_limit(delay);
            (count, ids)
        }
        Err(e) => {
            eprintln!("  [PubMed error] {}", e);
            rate_limit(delay);
            (0, vec![])
        }
    }
}

#[cfg(not(target_arch = "wasm32"))]
fn pubmed_summaries(pmids: &[String], api_key: &str, delay: f64) -> Vec<TopPaper> {
    if pmids.is_empty() {
        return vec![];
    }
    let ids_joined = pmids[..pmids.len().min(MAX_LANDSCAPE_FETCH)].join(",");
    let mut url = format!(
        "{}/esummary.fcgi?db=pubmed&id={}&retmode=json",
        EUTILS_BASE, ids_joined
    );
    if !api_key.is_empty() {
        url.push_str(&format!("&api_key={}", api_key));
    }

    match http_get_json(&url, MAX_RETRIES) {
        Ok(data) => {
            let result = &data["result"];
            let mut papers = Vec::new();
            for pid in pmids {
                if let Some(r) = result.get(pid) {
                    if r.get("error").is_some() {
                        continue;
                    }
                    papers.push(TopPaper {
                        pmid: pid.clone(),
                        title: r["title"]
                            .as_str()
                            .unwrap_or("")
                            .chars()
                            .take(120)
                            .collect(),
                        year: r["pubdate"]
                            .as_str()
                            .unwrap_or("")
                            .split_whitespace()
                            .next()
                            .unwrap_or("")
                            .to_string(),
                        journal: r["source"].as_str().unwrap_or("").to_string(),
                    });
                }
            }
            rate_limit(delay);
            papers
        }
        Err(e) => {
            eprintln!("  [PubMed summary error] {}", e);
            rate_limit(delay);
            vec![]
        }
    }
}

#[cfg(not(target_arch = "wasm32"))]
fn europepmc_search(query: &str, delay: f64) -> i64 {
    let encoded = urlencoded(query);
    let url = format!(
        "{}/search?query={}&format=json&pageSize=1&resultType=lite",
        EUROPEPMC_BASE, encoded
    );
    match http_get_json(&url, MAX_RETRIES) {
        Ok(data) => {
            let count = data["hitCount"].as_i64().unwrap_or(0);
            rate_limit(delay);
            count
        }
        Err(e) => {
            eprintln!("  [EuropePMC error] {}", e);
            rate_limit(delay);
            -1
        }
    }
}

fn urlencoded(s: &str) -> String {
    let mut result = String::with_capacity(s.len() * 2);
    for b in s.bytes() {
        match b {
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'-' | b'_' | b'.' | b'~' => {
                result.push(b as char);
            }
            b' ' => result.push('+'),
            _ => {
                result.push('%');
                result.push(char::from(b"0123456789ABCDEF"[(b >> 4) as usize]));
                result.push(char::from(b"0123456789ABCDEF"[(b & 0x0F) as usize]));
            }
        }
    }
    result
}

// ── Ontology loading ─────────────────────────────────────────────────────────

/// Load ontology-enrichment.json and return auto-generated synonyms per pathogen.
fn load_ontology_synonyms(
    ontology_path: &Path,
    configs: &BTreeMap<String, PathogenConfig>,
) -> (HashMap<String, HashMap<String, Vec<String>>>, usize, usize) {
    let mut total_genes = 0usize;
    let mut total_synonyms_added = 0usize;
    let mut result: HashMap<String, HashMap<String, Vec<String>>> = HashMap::new();

    let data = match std::fs::read_to_string(ontology_path) {
        Ok(s) => s,
        Err(e) => {
            eprintln!(
                "WARNING: Could not read ontology file {}: {}",
                ontology_path.display(),
                e
            );
            return (result, 0, 0);
        }
    };

    let ontology: OntologyFile = match serde_json::from_str(&data) {
        Ok(o) => o,
        Err(e) => {
            eprintln!("WARNING: Could not parse ontology JSON: {}", e);
            return (result, 0, 0);
        }
    };

    for (pathogen_key, cfg) in configs {
        let ont_pathogen = match ontology.pathogens.get(pathogen_key) {
            Some(p) => p,
            None => continue,
        };

        let gene_list = match &ont_pathogen.genes {
            Some(g) => &g.genes,
            None => continue,
        };

        total_genes += gene_list.len();

        // Build lookups: lowercase symbol → (symbol, name)
        let mut by_symbol: HashMap<String, (&str, &str)> = HashMap::new();
        let mut by_name: HashMap<String, (&str, &str)> = HashMap::new();
        for g in gene_list {
            let sym = g.symbol.as_deref().unwrap_or("").trim();
            let name = g.name.as_deref().unwrap_or("").trim();
            if !sym.is_empty() {
                by_symbol.insert(sym.to_lowercase(), (sym, name));
            }
            if !name.is_empty() {
                by_name.insert(name.to_lowercase(), (sym, name));
            }
        }

        let mut pathogen_synonyms: HashMap<String, Vec<String>> = HashMap::new();

        for gene in &cfg.genes {
            let gene_lower = gene.to_lowercase();
            let matched: Option<(&str, &str)> = by_symbol
                .get(&gene_lower)
                .copied()
                // Strategy 2: strip " protein" suffix
                .or_else(|| {
                    gene_lower
                        .strip_suffix(" protein")
                        .and_then(|stem| by_symbol.get(stem).copied())
                })
                // Strategy 3: match by name
                .or_else(|| by_name.get(&gene_lower).copied());

            if let Some((sym, name)) = matched {
                let mut additions: HashSet<String> = HashSet::new();
                if !name.is_empty() && name.to_lowercase() != gene_lower {
                    additions.insert(name.to_string());
                }
                if !sym.is_empty() && sym.to_lowercase() != gene_lower {
                    additions.insert(sym.to_string());
                }
                if !additions.is_empty() {
                    total_synonyms_added += additions.len();
                    let mut sorted: Vec<String> = additions.into_iter().collect();
                    sorted.sort();
                    pathogen_synonyms.insert(gene.clone(), sorted);
                }
            }
        }

        if !pathogen_synonyms.is_empty() {
            result.insert(pathogen_key.clone(), pathogen_synonyms);
        }
    }

    (result, total_genes, total_synonyms_added)
}

/// Merge ontology synonyms into manual synonyms. Manual entries take priority;
/// ontology adds additional aliases not already present.
fn merge_synonyms(
    configs: &mut BTreeMap<String, PathogenConfig>,
    ontology_synonyms: &HashMap<String, HashMap<String, Vec<String>>>,
) {
    for (pathogen_key, ont_syns) in ontology_synonyms {
        if let Some(cfg) = configs.get_mut(pathogen_key) {
            for (gene, aliases) in ont_syns {
                let existing: HashSet<String> = cfg
                    .gene_synonyms
                    .get(gene)
                    .map(|v| v.iter().map(|s| s.to_lowercase()).collect())
                    .unwrap_or_default();

                let new_aliases: Vec<String> = aliases
                    .iter()
                    .filter(|a| !existing.contains(&a.to_lowercase()))
                    .cloned()
                    .collect();

                if !new_aliases.is_empty() {
                    cfg.gene_synonyms
                        .entry(gene.clone())
                        .or_default()
                        .extend(new_aliases);
                }
            }
        }
    }
}

// ── Query builders ───────────────────────────────────────────────────────────

fn organism_clause_pm(cfg: &PathogenConfig) -> String {
    let mut names = vec![cfg.search_name.clone()];
    names.extend(cfg.alt_names.iter().cloned());
    let parts: Vec<String> = names
        .iter()
        .map(|n| format!("\"{}\"[Title/Abstract]", n))
        .collect();
    format!("({})", parts.join(" OR "))
}

fn organism_clause_epmc(cfg: &PathogenConfig) -> String {
    let name = cfg.search_name.split(" OR ").next().unwrap_or(&cfg.search_name).trim();
    format!("\"{}\"", name)
}

fn crispr_clause_pm() -> &'static str {
    "(CRISPR[Title/Abstract] OR \"guide RNA\"[Title/Abstract] \
     OR Cas9[Title/Abstract] OR Cas12[Title/Abstract] \
     OR Cas13[Title/Abstract] OR \"gene editing\"[Title/Abstract])"
}

fn crispr_clause_epmc() -> &'static str {
    "(CRISPR OR \"guide RNA\" OR Cas9 OR Cas12 OR Cas13 OR \"gene editing\")"
}

fn gene_clause_pm(gene: &str, synonyms: &HashMap<String, Vec<String>>) -> String {
    let mut all = vec![gene.to_string()];
    if let Some(syns) = synonyms.get(gene) {
        all.extend(syns.iter().cloned());
    }
    let parts: Vec<String> = all
        .iter()
        .map(|g| format!("\"{}\"[Title/Abstract]", g))
        .collect();
    format!("({})", parts.join(" OR "))
}

fn gene_clause_epmc(gene: &str, synonyms: &HashMap<String, Vec<String>>) -> String {
    let mut all = vec![gene.to_string()];
    if let Some(syns) = synonyms.get(gene) {
        all.extend(syns.iter().cloned());
    }
    let parts: Vec<String> = all
        .iter()
        .map(|g| format!("\"{}\"", g))
        .collect();
    format!("({})", parts.join(" OR "))
}

fn app_clause_pm_narrow(term: &str) -> String {
    format!("\"{}\"[Title/Abstract]", term)
}

fn app_clause_pm_broad(terms: &[String]) -> String {
    let parts: Vec<String> = terms
        .iter()
        .map(|t| format!("\"{}\"[Title/Abstract]", t))
        .collect();
    format!("({})", parts.join(" OR "))
}

fn app_clause_epmc_broad(terms: &[String]) -> String {
    let parts: Vec<String> = terms
        .iter()
        .map(|t| format!("\"{}\"", t))
        .collect();
    format!("({})", parts.join(" OR "))
}

// ── Corpus scan ──────────────────────────────────────────────────────────────

fn corpus_matches(corpus: &[TopPaper], terms: &[String]) -> Vec<String> {
    let norm_terms: Vec<String> = terms.iter().map(|t| t.to_lowercase()).collect();
    let mut hits = Vec::new();
    for paper in corpus {
        let hay = format!("{} {}", paper.title, paper.journal).to_lowercase();
        if norm_terms.iter().any(|t| hay.contains(t.as_str())) {
            hits.push(paper.title.chars().take(100).collect());
        }
    }
    hits
}

// ── Landscape + corpus fetch (with disk cache) ──────────────────────────────

#[cfg(not(target_arch = "wasm32"))]
fn fetch_landscape_corpus(
    key: &str,
    cfg: &PathogenConfig,
    cache_dir: &Path,
    api_key: &str,
    pm_delay: f64,
) -> (i64, Vec<TopPaper>) {
    let cache_file = cache_dir.join(format!("corpus_{}.json", key));

    // Check cache
    if cache_file.exists() {
        if let Ok(data) = std::fs::read_to_string(&cache_file) {
            if let Ok(cached) = serde_json::from_str::<serde_json::Value>(&data) {
                let total = cached["total"].as_i64().unwrap_or(0);
                let papers: Vec<TopPaper> = cached["papers"]
                    .as_array()
                    .map(|arr| {
                        arr.iter()
                            .filter_map(|v| serde_json::from_value(v.clone()).ok())
                            .collect()
                    })
                    .unwrap_or_default();
                eprintln!("  (corpus cache hit: {} papers)", papers.len());
                return (total, papers);
            }
        }
    }

    let org = organism_clause_pm(cfg);
    let crispr = crispr_clause_pm();
    let query = format!("{} AND {}", org, crispr);
    let (total, pmids) = pubmed_search(&query, MAX_LANDSCAPE_FETCH, api_key, pm_delay);

    let mut papers: Vec<TopPaper> = Vec::new();
    // Batch fetch summaries: max 200/request
    let chunks: Vec<&[String]> = pmids.chunks(200).collect();
    for chunk in chunks {
        let chunk_vec: Vec<String> = chunk.to_vec();
        papers.extend(pubmed_summaries(&chunk_vec, api_key, pm_delay));
        rate_limit(pm_delay * 2.0);
    }

    // Write cache
    let cache_data = serde_json::json!({
        "total": total,
        "papers": papers,
    });
    if let Ok(json_str) = serde_json::to_string(&cache_data) {
        let _ = std::fs::write(&cache_file, json_str);
    }

    (total, papers)
}

// ── Per-claim multi-strategy scoring ─────────────────────────────────────────

#[cfg(not(target_arch = "wasm32"))]
fn query_and_score(
    label: &str,
    pm_narrow_query: &str,
    pm_broad_query: &str,
    epmc_query: &str,
    corpus: &[TopPaper],
    corpus_terms: &[String],
    synonyms_used: Vec<String>,
    synonym_source: &str,
    api_key: &str,
    pm_delay: f64,
    epmc_delay: f64,
) -> ClaimResult {
    // Strategy 1: narrow
    let (n_narrow, narrow_ids) = pubmed_search(pm_narrow_query, MAX_DETAIL_FETCH, api_key, pm_delay);

    // Strategy 2: broad
    let (n_broad, broad_ids) = if pm_narrow_query == pm_broad_query {
        (n_narrow, narrow_ids.clone())
    } else {
        pubmed_search(pm_broad_query, MAX_DETAIL_FETCH, api_key, pm_delay)
    };

    // Strategy 3: Europe PMC
    let n_epmc = europepmc_search(epmc_query, epmc_delay);

    // Strategy 4: corpus scan
    let corpus_hits = corpus_matches(corpus, corpus_terms);

    // Fetch top paper summaries
    let mut all_ids: Vec<String> = Vec::new();
    let mut seen = HashSet::new();
    for id in narrow_ids.iter().take(3).chain(broad_ids.iter().take(3)) {
        if seen.insert(id.clone()) {
            all_ids.push(id.clone());
        }
    }
    let top_papers = if all_ids.is_empty() {
        vec![]
    } else {
        pubmed_summaries(&all_ids, api_key, pm_delay)
    };

    // Confidence scoring (same rules as v3):
    //   FALSE     — PubMed broad > 5
    //   UNCERTAIN — PubMed broad > 0, OR corpus title match, OR EuropePMC > 30
    //   PROBABLE  — PubMed narrow > 0 but PubMed broad = 0
    //   CONFIRMED — all counts == 0 and corpus == 0
    let confidence = if n_broad > 5 {
        Confidence::False
    } else if n_broad > 0 || !corpus_hits.is_empty() || n_epmc > 30 {
        Confidence::Uncertain
    } else if n_narrow > 0 {
        Confidence::Probable
    } else {
        Confidence::Confirmed
    };

    let mut queries = BTreeMap::new();
    queries.insert("pubmed_narrow".into(), pm_narrow_query.to_string());
    queries.insert("pubmed_broad".into(), pm_broad_query.to_string());
    queries.insert("europepmc".into(), epmc_query.to_string());

    ClaimResult {
        label: label.to_string(),
        confidence,
        counts: QueryCounts {
            pubmed_narrow: n_narrow,
            pubmed_broad: n_broad,
            europepmc: n_epmc,
            corpus_title_matches: corpus_hits.len(),
        },
        queries,
        corpus_hits: corpus_hits.into_iter().take(5).collect(),
        top_papers: top_papers.into_iter().take(5).collect(),
        synonyms_used,
        synonym_source: synonym_source.to_string(),
    }
}

// ── Pathogen scanning ────────────────────────────────────────────────────────

#[cfg(not(target_arch = "wasm32"))]
fn scan_pathogen(
    key: &str,
    cfg: &PathogenConfig,
    cache_dir: &Path,
    api_key: &str,
    pm_delay: f64,
    epmc_delay: f64,
    ontology_synonyms: &HashMap<String, HashMap<String, Vec<String>>>,
) -> PathogenResult {
    eprintln!("\n{}", "=".repeat(60));
    eprintln!("  {} ({})", cfg.search_name, key);
    eprintln!("{}", "=".repeat(60));

    // Fetch landscape + corpus
    eprintln!("  Fetching landscape + corpus...");
    let (total_papers, corpus) = fetch_landscape_corpus(key, cfg, cache_dir, api_key, pm_delay);
    eprintln!(
        "  Landscape: {} total CRISPR papers, {} abstracts fetched",
        total_papers,
        corpus.len()
    );

    let org_pm = organism_clause_pm(cfg);
    let org_epmc = organism_clause_epmc(cfg);
    let crispr_pm = crispr_clause_pm();
    let crispr_epmc = crispr_clause_epmc();

    let ont_syns = ontology_synonyms.get(key);

    let mut gene_results = BTreeMap::new();
    let mut app_results = BTreeMap::new();

    // ── Gene-level scan ──────────────────────────────────────────────────
    eprintln!("\n  Gene scan ({} genes):", cfg.genes.len());
    for gene in &cfg.genes {
        let gene_cl_pm = gene_clause_pm(gene, &cfg.gene_synonyms);
        let gene_cl_epmc = gene_clause_epmc(gene, &cfg.gene_synonyms);
        let mut all_gene_terms: Vec<String> = vec![gene.clone()];
        if let Some(syns) = cfg.gene_synonyms.get(gene) {
            all_gene_terms.extend(syns.iter().cloned());
        }

        let pm_narrow = format!("{} AND {} AND {}", org_pm, crispr_pm, gene_cl_pm);
        let pm_broad = pm_narrow.clone(); // gene query already uses full synonym expansion
        let epmc_q = format!("{} AND {} AND {}", org_epmc, crispr_epmc, gene_cl_epmc);

        // Determine synonym source
        let has_manual = cfg.gene_synonyms.contains_key(gene);
        let has_ontology = ont_syns.map_or(false, |s| s.contains_key(gene));
        let syn_source = match (has_manual, has_ontology) {
            (true, true) => "manual+ontology",
            (true, false) => "manual",
            (false, true) => "ontology",
            (false, false) => "none",
        };

        let result = query_and_score(
            &format!("{}:gene:{}", key, gene),
            &pm_narrow,
            &pm_broad,
            &epmc_q,
            &corpus,
            &all_gene_terms,
            cfg.gene_synonyms.get(gene).cloned().unwrap_or_default(),
            syn_source,
            api_key,
            pm_delay,
            epmc_delay,
        );

        let marker = match result.confidence {
            Confidence::Confirmed => "○",
            Confidence::Probable => "?",
            Confidence::Uncertain => "⚠",
            Confidence::False => "●",
        };
        eprintln!(
            "    {} {:25}  pm={:>3} pm_broad={:>3} epmc={:>4} corpus={:>2}  → {} [{}]",
            marker,
            gene,
            result.counts.pubmed_narrow,
            result.counts.pubmed_broad,
            result.counts.europepmc,
            result.counts.corpus_title_matches,
            result.confidence,
            syn_source
        );

        gene_results.insert(gene.clone(), result);
    }

    // ── Application-level scan ───────────────────────────────────────────
    eprintln!(
        "\n  Application scan ({} categories):",
        cfg.applications.len()
    );
    for (app_cat, app_terms) in &cfg.applications {
        let primary = &app_terms[0];

        let pm_narrow = format!(
            "{} AND {} AND {}",
            org_pm,
            crispr_pm,
            app_clause_pm_narrow(primary)
        );
        let pm_broad = format!(
            "{} AND {} AND {}",
            org_pm,
            crispr_pm,
            app_clause_pm_broad(app_terms)
        );
        let epmc_q = format!(
            "{} AND {} AND {}",
            org_epmc,
            crispr_epmc,
            app_clause_epmc_broad(app_terms)
        );

        let result = query_and_score(
            &format!("{}:app:{}", key, app_cat),
            &pm_narrow,
            &pm_broad,
            &epmc_q,
            &corpus,
            app_terms,
            app_terms.clone(),
            "manual",
            api_key,
            pm_delay,
            epmc_delay,
        );

        let marker = match result.confidence {
            Confidence::Confirmed => "○",
            Confidence::Probable => "?",
            Confidence::Uncertain => "⚠",
            Confidence::False => "●",
        };
        eprintln!(
            "    {} {:22}  pm={:>3} pm_broad={:>3} epmc={:>4} corpus={:>2}  → {}",
            marker,
            app_cat,
            result.counts.pubmed_narrow,
            result.counts.pubmed_broad,
            result.counts.europepmc,
            result.counts.corpus_title_matches,
            result.confidence
        );

        app_results.insert(app_cat.clone(), result);
    }

    // Summarise
    let confirmed_gene_gaps: Vec<String> = gene_results
        .iter()
        .filter(|(_, r)| r.confidence == Confidence::Confirmed)
        .map(|(g, _)| g.clone())
        .collect();
    let confirmed_app_gaps: Vec<String> = app_results
        .iter()
        .filter(|(_, r)| r.confidence == Confidence::Confirmed)
        .map(|(a, _)| a.clone())
        .collect();
    let uncertain: Vec<String> = gene_results
        .iter()
        .filter(|(_, r)| matches!(r.confidence, Confidence::Probable | Confidence::Uncertain))
        .map(|(g, _)| format!("gene:{}", g))
        .chain(
            app_results
                .iter()
                .filter(|(_, r)| {
                    matches!(r.confidence, Confidence::Probable | Confidence::Uncertain)
                })
                .map(|(a, _)| format!("app:{}", a)),
        )
        .collect();
    let false_claims: Vec<String> = gene_results
        .iter()
        .filter(|(_, r)| r.confidence == Confidence::False)
        .map(|(g, _)| format!("gene:{}", g))
        .chain(
            app_results
                .iter()
                .filter(|(_, r)| r.confidence == Confidence::False)
                .map(|(a, _)| format!("app:{}", a)),
        )
        .collect();

    eprintln!("\n  Summary for {}:", cfg.search_name);
    eprintln!("    CONFIRMED gaps:  genes={:?}, apps={:?}", confirmed_gene_gaps, confirmed_app_gaps);
    eprintln!("    UNCERTAIN:       {:?}", uncertain);
    eprintln!("    FALSE claims:    {:?}", false_claims);

    PathogenResult {
        pathogen_key: key.to_string(),
        search_name: cfg.search_name.clone(),
        landscape_total: total_papers,
        corpus_fetched: corpus.len(),
        genes: gene_results,
        applications: app_results,
        confirmed_gene_gaps,
        confirmed_app_gaps,
        uncertain,
        false_claims,
    }
}

// ── Audit file ───────────────────────────────────────────────────────────────

fn write_audit(output: &ScanOutput, path: &Path) -> std::io::Result<()> {
    use std::io::Write;
    let mut f = std::fs::File::create(path)?;

    // Safe UTF-8 truncation: never split in the middle of a char
    fn trunc(s: &str, max: usize) -> &str {
        if s.len() <= max { return s; }
        let mut end = max;
        while end > 0 && !s.is_char_boundary(end) { end -= 1; }
        &s[..end]
    }

    writeln!(f, "{}", "=".repeat(70))?;
    writeln!(f, "LOOM PUBMED SCAN v4 — AUDIT TRAIL (Rust, ontology-integrated)")?;
    writeln!(f, "Generated: {}", output.generated)?;
    writeln!(f, "Ontology: {}", output.ontology_source)?;
    writeln!(
        f,
        "Total confirmed gaps: {}",
        output.summary.confirmed_gaps
    )?;
    writeln!(
        f,
        "Uncertain/needs review: {}",
        output.summary.uncertain_count
    )?;
    writeln!(f, "Refuted (FALSE) claims: {}", output.summary.false_count)?;
    writeln!(
        f,
        "Ontology genes loaded: {}, synonyms added: {}",
        output.summary.ontology_genes_loaded, output.summary.ontology_synonyms_added
    )?;
    writeln!(f, "{}", "=".repeat(70))?;
    writeln!(f)?;
    writeln!(f, "CONFIDENCE KEY:")?;
    writeln!(
        f,
        "  CONFIRMED — 0 results across all query strategies + sources + corpus"
    )?;
    writeln!(
        f,
        "  PROBABLE  — 0 in narrow PubMed query only; not verified elsewhere"
    )?;
    writeln!(
        f,
        "  UNCERTAIN — 0 in narrow, but >0 in broad/EuropePMC or corpus hit"
    )?;
    writeln!(
        f,
        "  FALSE     — >5 results in broad query or EuropePMC; NOT a gap"
    )?;
    writeln!(f)?;

    for (key, pr) in &output.pathogens {
        writeln!(f, "{}", "─".repeat(60))?;
        writeln!(f, "PATHOGEN: {} ({})", pr.search_name, key)?;
        writeln!(f, "  Landscape: {} total CRISPR papers", pr.landscape_total)?;
        writeln!(f, "  Corpus fetched: {} abstracts scanned", pr.corpus_fetched)?;
        writeln!(f)?;

        writeln!(f, "  GENE GAPS:")?;
        for (gene, r) in &pr.genes {
            writeln!(
                f,
                "    [{:9}] {:25} pm_narrow={:>3}  pm_broad={:>3}  epmc={:>4}  corpus={:>2}  syns={}",
                r.confidence.to_string(),
                gene,
                r.counts.pubmed_narrow,
                r.counts.pubmed_broad,
                r.counts.europepmc,
                r.counts.corpus_title_matches,
                r.synonym_source
            )?;
            if !r.synonyms_used.is_empty() {
                writeln!(f, "               synonyms: {:?}", &r.synonyms_used[..r.synonyms_used.len().min(8)])?;
            }
            for h in r.corpus_hits.iter().take(3) {
                writeln!(f, "               corpus: {}", trunc(h, 80))?;
            }
            if matches!(r.confidence, Confidence::Uncertain | Confidence::False) {
                for p in r.top_papers.iter().take(2) {
                    writeln!(
                        f,
                        "               paper {}: {}",
                        p.year,
                        trunc(&p.title, 80)
                    )?;
                }
            }
        }

        writeln!(f)?;
        writeln!(f, "  APPLICATION GAPS:")?;
        for (app, r) in &pr.applications {
            writeln!(
                f,
                "    [{:9}] {:22} pm_narrow={:>3}  pm_broad={:>3}  epmc={:>4}  corpus={:>2}",
                r.confidence.to_string(),
                app,
                r.counts.pubmed_narrow,
                r.counts.pubmed_broad,
                r.counts.europepmc,
                r.counts.corpus_title_matches
            )?;
            for h in r.corpus_hits.iter().take(3) {
                writeln!(f, "               corpus: {}", &h[..h.len().min(80)])?;
            }
        }

        writeln!(f)?;
        if !pr.false_claims.is_empty() {
            writeln!(
                f,
                "  *** FALSE CLAIMS (must be removed/corrected): {:?}",
                pr.false_claims
            )?;
        }
        if !pr.uncertain.is_empty() {
            writeln!(f, "  *** NEEDS MANUAL REVIEW: {:?}", pr.uncertain)?;
        }
        writeln!(f)?;
    }

    Ok(())
}

// ── Pathogen config (same data as v3, but without gene_synonyms — those come from ontology) ──

fn build_pathogen_configs() -> BTreeMap<String, PathogenConfig> {
    let mut m = BTreeMap::new();

    // Helper macro for conciseness
    macro_rules! cfg {
        ($key:expr, $search:expr, $alt:expr, $mesh:expr, $genes:expr,
         $syns:expr, $apps:expr) => {
            m.insert(
                $key.to_string(),
                PathogenConfig {
                    search_name: $search.to_string(),
                    alt_names: $alt.iter().map(|s: &&str| s.to_string()).collect(),
                    mesh_term: if $mesh.is_empty() {
                        None
                    } else {
                        Some($mesh.to_string())
                    },
                    genes: $genes.iter().map(|s: &&str| s.to_string()).collect(),
                    gene_synonyms: $syns,
                    applications: $apps,
                },
            );
        };
    }

    fn syns(pairs: &[(&str, &[&str])]) -> HashMap<String, Vec<String>> {
        pairs
            .iter()
            .map(|(k, v)| (k.to_string(), v.iter().map(|s| s.to_string()).collect()))
            .collect()
    }

    fn apps(pairs: &[(&str, &[&str])]) -> BTreeMap<String, Vec<String>> {
        pairs
            .iter()
            .map(|(k, v)| (k.to_string(), v.iter().map(|s| s.to_string()).collect()))
            .collect()
    }

    cfg!(
        "cholera",
        "Vibrio cholerae",
        &["cholera", "V. cholerae"],
        "Vibrio cholerae[MeSH]",
        &["ctxA", "ctxB", "tcpA", "toxR", "ompU", "ompT", "hapA", "hlyA", "rtxA"],
        syns(&[
            ("ctxA", &["cholera toxin A", "CT-A", "cholera toxin subunit A", "ctx operon", "ctxAB"]),
            ("ctxB", &["cholera toxin B", "CT-B", "cholera toxin subunit B", "ctxAB"]),
            ("tcpA", &["toxin-coregulated pilus", "TCP pilus", "tcp operon"]),
            ("toxR", &["ToxR", "virulence regulator"]),
            ("ompU", &["outer membrane protein U"]),
            ("rtxA", &["RTX toxin"]),
            ("hlyA", &["haemolysin A", "hemolysin"]),
        ]),
        apps(&[
            ("diagnostics", &[
                "diagnostic", "diagnostics", "detection assay", "detection", "identify",
                "identification", "rapid test", "rapid detection", "test strip", "biosensor",
                "lateral flow", "POC", "point-of-care", "SHERLOCK", "DETECTR", "Cas12", "Cas13",
                "isothermal", "LAMP", "field-deployable", "serology", "serogroup", "serovar",
                "serotyping", "environmental surveillance", "water testing", "wastewater",
            ]),
            ("therapeutics", &[
                "antimicrobial", "antibacterial", "therapy", "treatment", "kill", "inhibit",
                "antibiotic",
            ]),
        ])
    );

    cfg!(
        "dengue",
        "Dengue",
        &["dengue virus", "DENV"],
        "Dengue virus[MeSH]",
        &["NS1", "NS2A", "NS2B", "NS3", "NS4A", "NS4B", "NS5", "capsid", "prM", "envelope"],
        syns(&[
            ("NS5", &["RNA-dependent RNA polymerase", "RdRp", "methyltransferase"]),
            ("NS3", &["helicase", "protease", "NS3 helicase"]),
            ("NS1", &["non-structural protein 1"]),
            ("prM", &["pre-membrane", "premembrane"]),
            ("envelope", &["E protein", "E gene"]),
        ]),
        apps(&[
            ("diagnostics", &[
                "diagnostic", "detection", "assay", "rapid test", "SHERLOCK", "DETECTR",
                "point-of-care", "Cas12", "Cas13", "isothermal", "LAMP", "biosensor",
                "lateral flow", "serotype", "serotyping",
            ]),
            ("therapeutics", &["antiviral", "therapy", "inhibit", "treatment"]),
        ])
    );

    cfg!(
        "ebola",
        "Ebola",
        &["Ebola virus", "EBOV", "ebolavirus", "filovirus"],
        "Ebolavirus[MeSH]",
        &["NP", "VP35", "VP40", "GP", "VP30", "VP24", "L protein"],
        syns(&[
            ("NP", &["nucleoprotein"]),
            ("GP", &["glycoprotein", "surface glycoprotein"]),
            ("L protein", &["RNA-dependent RNA polymerase", "RdRp", "L gene"]),
            ("VP40", &["matrix protein"]),
        ]),
        apps(&[
            ("diagnostics", &[
                "diagnostic", "detection", "assay", "rapid test", "SHERLOCK", "DETECTR",
                "point-of-care", "Cas12", "Cas13", "isothermal", "field-deployable",
                "biosensor", "outbreak",
            ]),
            ("therapeutics", &["antiviral", "therapy", "treatment", "inhibit", "knockdown"]),
        ])
    );

    cfg!(
        "mpox",
        "Mpox OR Monkeypox",
        &["MPXV", "mpox virus", "orthopoxvirus", "monkeypox virus"],
        "Monkeypox virus[MeSH]",
        &["A33R", "B5R", "B6R", "E8L", "J2R", "A56R", "thymidine kinase", "hemagglutinin"],
        syns(&[
            ("A33R", &["EEV-specific protein", "extracellular enveloped virus"]),
            ("B5R", &["EEV membrane antigen"]),
            ("B6R", &["B6R gene"]),
            ("thymidine kinase", &["TK gene", "J2R thymidine kinase"]),
            ("hemagglutinin", &["HA protein", "A56R"]),
        ]),
        apps(&[
            ("diagnostics", &[
                "diagnostic", "detection", "assay", "rapid test", "SHERLOCK", "DETECTR",
                "Cas12", "Cas13", "point-of-care", "isothermal", "PCR-free", "surveillance",
                "biosensor",
            ]),
            ("therapeutics", &["antiviral", "therapy", "treatment", "inhibit"]),
        ])
    );

    cfg!(
        "tuberculosis",
        "Mycobacterium tuberculosis",
        &["M. tuberculosis", "MTB", "TB"],
        "Mycobacterium tuberculosis[MeSH]",
        &["rpoB", "katG", "inhA", "gyrA", "gyrB", "embB", "pncA", "ethA", "rrs", "eis",
          "IS6110", "Ag85B", "ESAT-6", "CFP-10"],
        syns(&[
            ("inhA", &["enoyl-ACP reductase", "isoniazid resistance", "Rv1484"]),
            ("gyrB", &["gyrase subunit B", "fluoroquinolone resistance"]),
            ("pncA", &["pyrazinamidase", "pyrazinamide resistance", "Rv2043c"]),
            ("rrs", &["16S rRNA", "16S ribosomal RNA", "aminoglycoside resistance"]),
            ("eis", &["enhanced intracellular survival", "kanamycin resistance", "Rv2416c"]),
            ("CFP-10", &["Rv3874", "ESAT-6/CFP-10", "culture filtrate protein 10", "EsxB"]),
            ("ethA", &["ethionamide monooxygenase", "ethionamide resistance", "Rv3854c"]),
            ("embB", &["arabinosyltransferase", "ethambutol resistance"]),
            ("rpoB", &["RNA polymerase beta", "rifampicin resistance", "rifampin resistance"]),
            ("katG", &["catalase-peroxidase", "isoniazid resistance KatG"]),
            ("IS6110", &["insertion sequence IS6110", "repeat element TB"]),
            ("ESAT-6", &["early secreted antigenic target", "EsxA", "Rv3875"]),
            ("Ag85B", &["antigen 85B", "Rv1886c", "fbpB"]),
        ]),
        apps(&[
            ("diagnostics", &[
                "diagnostic", "detection", "drug resistance", "MDR-TB", "XDR-TB",
                "rifampicin resistance", "isoniazid resistance", "rapid test", "SHERLOCK",
                "DETECTR", "Cas12", "Cas13", "point-of-care", "isothermal", "LAMP", "sputum",
            ]),
            ("therapeutics", &[
                "antimicrobial", "antibacterial", "treatment", "inhibit", "gene silencing",
                "CRISPRi",
            ]),
        ])
    );

    cfg!(
        "zika",
        "Zika",
        &["Zika virus", "ZIKV"],
        "Zika Virus[MeSH]",
        &["NS1", "NS2A", "NS2B", "NS3", "NS4A", "NS4B", "NS5", "capsid", "prM", "envelope"],
        syns(&[
            ("NS5", &["RNA-dependent RNA polymerase", "RdRp", "methyltransferase"]),
            ("NS3", &["helicase", "serine protease"]),
            ("NS2B", &["NS2B cofactor", "protease cofactor"]),
            ("prM", &["pre-membrane protein"]),
            ("envelope", &["E protein", "E gene"]),
        ]),
        apps(&[
            ("diagnostics", &[
                "diagnostic", "detection", "assay", "rapid test", "SHERLOCK", "DETECTR",
                "Cas12", "Cas13", "point-of-care", "isothermal", "LAMP",
                "microcephaly screening", "congenital",
            ]),
            ("therapeutics", &["antiviral", "therapy", "treatment"]),
        ])
    );

    cfg!(
        "rsv",
        "RSV",
        &["respiratory syncytial virus", "hRSV", "human RSV"],
        "Respiratory Syncytial Viruses[MeSH]",
        &["F protein", "G protein", "N protein", "L protein", "M protein", "M2-1", "SH protein"],
        syns(&[
            ("F protein", &["fusion protein", "F gene", "RSV F"]),
            ("G protein", &["attachment glycoprotein", "G gene"]),
            ("N protein", &["nucleoprotein", "N gene"]),
            ("L protein", &["RNA-dependent RNA polymerase", "L gene", "RdRp"]),
            ("SH protein", &["small hydrophobic protein"]),
        ]),
        apps(&[
            ("diagnostics", &[
                "diagnostic", "detection", "assay", "rapid test", "SHERLOCK", "DETECTR",
                "Cas12", "Cas13", "point-of-care", "isothermal", "LAMP",
                "pediatric", "infant", "bronchiolitis",
            ]),
            ("therapeutics", &["antiviral", "therapy", "inhibit", "palivizumab"]),
        ])
    );

    cfg!(
        "influenza-a",
        "Influenza A",
        &["influenza virus", "flu", "H1N1", "H3N2", "H5N1", "avian influenza"],
        "Influenza A virus[MeSH]",
        &["HA", "NA", "PB1", "PB2", "PA", "NP", "M1", "M2", "NS1", "NEP"],
        syns(&[
            ("HA", &["hemagglutinin", "H1", "H3", "H5"]),
            ("NA", &["neuraminidase", "N1", "N2"]),
            ("PB1", &["polymerase basic 1", "PB1-F2"]),
            ("PB2", &["polymerase basic 2"]),
            ("PA", &["polymerase acidic"]),
            ("NP", &["nucleoprotein"]),
            ("NS1", &["non-structural protein 1", "interferon antagonist"]),
            ("M2", &["M2 ion channel", "amantadine target"]),
        ]),
        apps(&[
            ("diagnostics", &[
                "diagnostic", "detection", "assay", "rapid test", "subtyping", "SHERLOCK",
                "DETECTR", "Cas12", "Cas13", "point-of-care", "isothermal", "surveillance",
                "pandemic",
            ]),
            ("therapeutics", &["antiviral", "therapy", "inhibit", "oseltamivir"]),
        ])
    );

    cfg!(
        "mers",
        "MERS",
        &["MERS-CoV", "Middle East respiratory syndrome",
          "Middle East respiratory syndrome coronavirus"],
        "Middle East Respiratory Syndrome Coronavirus[MeSH]",
        &["spike", "RdRp", "ORF1a", "ORF1b", "nucleocapsid", "envelope", "membrane"],
        syns(&[
            ("spike", &["S protein", "spike glycoprotein", "receptor-binding domain", "RBD"]),
            ("RdRp", &["nsp12", "RNA-dependent RNA polymerase"]),
            ("ORF1a", &["nsp1", "nsp3", "papain-like protease", "PLpro"]),
            ("ORF1b", &["nsp13", "nsp14", "helicase"]),
            ("nucleocapsid", &["N protein", "N gene"]),
        ]),
        apps(&[
            ("diagnostics", &[
                "diagnostic", "detection", "assay", "rapid test", "SHERLOCK", "DETECTR",
                "Cas12", "Cas13", "point-of-care", "isothermal", "camel",
                "zoonotic surveillance",
            ]),
            ("therapeutics", &["antiviral", "therapy", "treatment", "inhibit"]),
        ])
    );

    cfg!(
        "hepatitis-b",
        "Hepatitis B",
        &["HBV", "hepatitis B virus"],
        "Hepatitis B virus[MeSH]",
        &["HBsAg", "HBcAg", "HBx", "polymerase", "precore", "cccDNA", "pgRNA"],
        syns(&[
            ("HBsAg", &["surface antigen", "HBs antigen", "HBsAg gene", "HBV surface"]),
            ("cccDNA", &["covalently closed circular DNA", "episomal DNA", "nuclear HBV"]),
            ("HBx", &["HBx protein", "X gene", "hepatitis B x antigen"]),
            ("polymerase", &["HBV polymerase", "reverse transcriptase", "RT gene"]),
            ("pgRNA", &["pregenomic RNA", "pg RNA", "3.5 kb RNA"]),
        ]),
        apps(&[
            ("diagnostics", &[
                "diagnostic", "detection", "assay", "surface antigen", "HBsAg detection",
                "viral load", "quantification",
            ]),
            ("therapeutics", &[
                "cure", "cccDNA disruption", "functional cure", "antiviral", "knockdown",
                "silencing", "gene editing", "excision",
            ]),
        ])
    );

    cfg!(
        "hiv-1",
        "HIV",
        &["HIV-1", "human immunodeficiency virus"],
        "HIV-1[MeSH]",
        &["gag", "pol", "env", "tat", "rev", "nef", "vif", "LTR", "integrase", "gp120"],
        syns(&[
            ("gag", &["capsid", "matrix", "p24", "gag gene"]),
            ("pol", &["reverse transcriptase", "integrase", "protease"]),
            ("env", &["gp120", "gp41", "envelope glycoprotein"]),
            ("tat", &["transactivator", "Tat protein"]),
            ("LTR", &["long terminal repeat", "HIV LTR"]),
            ("integrase", &["IN gene", "HIV integrase"]),
        ]),
        apps(&[
            ("diagnostics", &[
                "diagnostic", "detection", "viral load", "quantification", "proviral DNA",
            ]),
            ("therapeutics", &[
                "cure", "latency reversal", "reservoir elimination", "excision",
                "gene editing", "knockdown", "CRISPRi", "CRISPR excision",
            ]),
        ])
    );

    cfg!(
        "sars-cov-2",
        "SARS-CoV-2",
        &["COVID-19", "SARS CoV 2", "2019-nCoV",
          "coronavirus disease 2019", "severe acute respiratory syndrome coronavirus 2"],
        "SARS-CoV-2[MeSH]",
        &["spike", "RdRp", "nucleocapsid", "envelope", "membrane",
          "nsp13", "nsp14", "nsp15", "nsp16",
          "nsp3", "nsp5", "nsp7", "nsp8", "nsp9", "nsp10",
          "ORF3a", "ORF8", "ORF10"],
        syns(&[
            ("spike", &["S protein", "S gene", "spike glycoprotein",
                        "receptor-binding domain", "RBD", "furin cleavage site", "D614G"]),
            ("RdRp", &["nsp12", "RNA-dependent RNA polymerase", "replicase",
                       "polymerase", "ORF1b polymerase"]),
            ("nucleocapsid", &["N protein", "N gene", "nucleocapsid protein"]),
            ("envelope", &["E protein", "E gene", "envelope protein"]),
            ("membrane", &["M protein", "M gene", "membrane protein"]),
            ("nsp13", &["helicase", "nsp-13", "NSP13", "SARS-CoV-2 helicase",
                        "coronavirus helicase", "nsp13 helicase"]),
            ("nsp14", &["exonuclease", "ExoN", "nsp-14", "NSP14",
                        "proofreading exonuclease", "N7-methyltransferase",
                        "nsp14 exonuclease"]),
            ("nsp15", &["endoribonuclease", "NendoU", "nsp-15", "NSP15",
                        "uridylate-specific endoribonuclease"]),
            ("nsp16", &["2'-O-methyltransferase", "nsp-16", "NSP16", "cap methyltransferase"]),
            ("nsp3", &["papain-like protease", "PLpro", "PL2pro", "nsp-3", "NSP3",
                       "ADP-ribose phosphatase", "macro domain"]),
            ("nsp5", &["main protease", "3CLpro", "Mpro", "3C-like protease", "nsp-5", "NSP5"]),
            ("nsp7", &["nsp-7", "NSP7", "primase cofactor"]),
            ("nsp8", &["nsp-8", "NSP8", "primase"]),
            ("nsp9", &["nsp-9", "NSP9", "ssRNA-binding protein"]),
            ("nsp10", &["nsp-10", "NSP10", "cap-0 methyltransferase activator"]),
            ("ORF3a", &["ORF 3a", "3a protein", "viroporin 3a"]),
            ("ORF8", &["ORF 8", "8 protein", "ORF8 protein"]),
            ("ORF10", &["ORF 10", "10 protein"]),
        ]),
        apps(&[
            ("diagnostics", &[
                "diagnostic", "diagnostics", "detection", "assay", "rapid test",
                "SHERLOCK", "DETECTR", "Cas12", "Cas13", "point-of-care", "POC",
                "isothermal", "LAMP", "biosensor", "lateral flow", "RT-LAMP",
                "CRISPR-based detection", "guide RNA", "field-deployable",
                "wastewater surveillance", "clinical testing",
            ]),
            ("therapeutics", &[
                "antiviral", "therapy", "treatment", "inhibit", "knockdown",
                "silencing", "gene editing", "CRISPRi", "CRISPR knockout",
            ]),
        ])
    );

    m
}

// ── Public entry point ───────────────────────────────────────────────────────

#[cfg(not(target_arch = "wasm32"))]
pub fn run(
    only: Option<&str>,
    no_merge: bool,
    no_cache: bool,
    output_path: Option<&Path>,
    ontology_path: &Path,
    cache_dir: &Path,
    audit_path: &Path,
) -> anyhow::Result<()> {
    let api_key = std::env::var("NCBI_API_KEY").unwrap_or_default();
    let pm_delay = if api_key.is_empty() {
        PUBMED_DELAY_NO_KEY
    } else {
        PUBMED_DELAY_WITH_KEY
    };

    let mut configs = build_pathogen_configs();

    // Load ontology and merge synonyms
    eprintln!("Loading ontology from {}...", ontology_path.display());
    let (ontology_synonyms, total_genes, total_syns_added) =
        load_ontology_synonyms(ontology_path, &configs);
    eprintln!(
        "  Ontology: {} genes loaded, {} synonyms added across {} pathogens",
        total_genes,
        total_syns_added,
        ontology_synonyms.len()
    );
    merge_synonyms(&mut configs, &ontology_synonyms);

    // Clear cache if requested
    if no_cache && cache_dir.exists() {
        eprintln!("Clearing cache directory...");
        let _ = std::fs::remove_dir_all(cache_dir);
    }
    std::fs::create_dir_all(cache_dir)?;

    // Determine which pathogens to scan
    let available_keys: Vec<String> = configs.keys().cloned().collect();
    let keys_to_scan: Vec<String> = if let Some(only_key) = only {
        if !configs.contains_key(only_key) {
            anyhow::bail!("No config for '{}' in PATHOGEN_CONFIG", only_key);
        }
        vec![only_key.to_string()]
    } else {
        available_keys.clone()
    };

    eprintln!("Scanning {} pathogens: {:?}", keys_to_scan.len(), keys_to_scan);
    eprintln!(
        "API key: {}",
        if api_key.is_empty() {
            "not set (3/s)"
        } else {
            "set (10/s)"
        }
    );

    // Estimate API calls
    let total_calls: usize = keys_to_scan
        .iter()
        .map(|k| {
            let c = &configs[k];
            1 + c.genes.len() * 3 + c.applications.len() * 3
        })
        .sum();
    let est_min = total_calls as f64 * pm_delay.max(EPMC_DELAY) / 60.0;
    eprintln!("Estimated: ~{} API calls, ~{:.0} min\n", total_calls, est_min);

    let start = Instant::now();
    let mut all_results: BTreeMap<String, PathogenResult> = BTreeMap::new();
    for key in &keys_to_scan {
        let cfg = &configs[key];
        let result = scan_pathogen(key, cfg, cache_dir, &api_key, pm_delay, EPMC_DELAY, &ontology_synonyms);
        all_results.insert(key.clone(), result);
    }

    // If --only mode, merge new result into existing output
    let default_output = PathBuf::from("web/data/pubmed-scan-v4.json");
    let out_file = output_path.unwrap_or(&default_output);
    if only.is_some() && !no_merge && out_file.exists() {
        eprintln!("\nMerging result into existing {} ...", out_file.display());
        if let Ok(data) = std::fs::read_to_string(out_file) {
            if let Ok(existing) = serde_json::from_str::<ScanOutput>(&data) {
                for (k, v) in existing.pathogens {
                    all_results.entry(k).or_insert(v);
                }
            }
        }
    }

    // Compile summary
    let all_confirmed_gene: Vec<GapEntry> = all_results
        .iter()
        .flat_map(|(k, r)| {
            r.confirmed_gene_gaps
                .iter()
                .map(move |g| GapEntry {
                    pathogen: k.clone(),
                    gene: g.clone(),
                })
        })
        .collect();
    let all_confirmed_app: Vec<AppGapEntry> = all_results
        .iter()
        .flat_map(|(k, r)| {
            r.confirmed_app_gaps
                .iter()
                .map(move |a| AppGapEntry {
                    pathogen: k.clone(),
                    application: a.clone(),
                })
        })
        .collect();
    let all_uncertain: Vec<ItemEntry> = all_results
        .iter()
        .flat_map(|(k, r)| {
            r.uncertain.iter().map(move |it| ItemEntry {
                pathogen: k.clone(),
                item: it.clone(),
            })
        })
        .collect();
    let all_false: Vec<ItemEntry> = all_results
        .iter()
        .flat_map(|(k, r)| {
            r.false_claims.iter().map(move |it| ItemEntry {
                pathogen: k.clone(),
                item: it.clone(),
            })
        })
        .collect();

    let output = ScanOutput {
        generated: chrono::Utc::now().to_rfc3339(),
        version: 4,
        methodology: format!(
            "Multi-strategy gap analysis v4 (Rust, ontology-integrated). \
             Each claimed gap is tested with (1) narrow PubMed query [Title/Abstract], \
             (2) broad PubMed query with full synonym expansion [Title/Abstract], \
             (3) Europe PMC independent search, \
             (4) local corpus scan of landscape abstracts. \
             Gene synonyms are auto-generated from NCBI ontology annotations \
             ({} genes across {} pathogens) and merged with manual synonyms. \
             Only claims with CONFIRMED confidence (0 in all strategies) are \
             reported as gaps. UNCERTAIN and FALSE claims are flagged for review.",
            total_genes, ontology_synonyms.len()
        ),
        ontology_source: format!(
            "web/data/ontology-enrichment.json ({} genes, {} auto-synonyms added)",
            total_genes, total_syns_added
        ),
        summary: ScanSummary {
            pathogens_scanned: all_results.len(),
            confirmed_gene_gaps: all_confirmed_gene.len(),
            confirmed_app_gaps: all_confirmed_app.len(),
            confirmed_gaps: all_confirmed_gene.len() + all_confirmed_app.len(),
            uncertain_count: all_uncertain.len(),
            false_count: all_false.len(),
            ontology_genes_loaded: total_genes,
            ontology_synonyms_added: total_syns_added,
        },
        confirmed_gene_gaps: all_confirmed_gene,
        confirmed_app_gaps: all_confirmed_app,
        uncertain: all_uncertain,
        false_claims: all_false,
        pathogens: all_results,
    };

    // Write output JSON
    let json = serde_json::to_string_pretty(&output)?;
    if let Some(parent) = out_file.parent() {
        std::fs::create_dir_all(parent)?;
    }
    std::fs::write(out_file, &json)?;
    eprintln!("\nResults → {}", out_file.display());

    // Write audit trail
    write_audit(&output, audit_path)?;
    eprintln!("Audit trail → {}", audit_path.display());

    let elapsed = start.elapsed();
    eprintln!("\n{}", "=".repeat(70));
    eprintln!("  PUBMED SCAN v4 — SUMMARY (Rust, ontology-integrated)");
    eprintln!("{}", "=".repeat(70));
    eprintln!("  Confirmed gene gaps:  {}", output.summary.confirmed_gene_gaps);
    eprintln!("  Confirmed app gaps:   {}", output.summary.confirmed_app_gaps);
    eprintln!("  Needs manual review:  {}", output.summary.uncertain_count);
    eprintln!("  False (retracted):    {}", output.summary.false_count);
    eprintln!("  Ontology synonyms:    {}", total_syns_added);
    eprintln!("  Total scan time:      {:.1}s", elapsed.as_secs_f64());
    eprintln!("{}", "-".repeat(70));
    if !output.false_claims.is_empty() {
        eprintln!("  FALSE CLAIMS (must be removed from paper):");
        for f in &output.false_claims {
            eprintln!("    {:15} | {}", f.pathogen, f.item);
        }
    }
    eprintln!("{}", "=".repeat(70));

    Ok(())
}
