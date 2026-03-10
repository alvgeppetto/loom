use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use loom::cli::{Cli, Commands, BuildModeArg};
use loom::dna::corpus::build_dna_corpus_jsonl;
use loom::dna::latency::{compute_latency_stats, time_call};
use loom::dna::manifest::{build_dna_manifest, write_dna_manifest_json};
use loom::dna::benchmark::{DnaBenchmarkConfig, run_dna_benchmark};
use loom::dna::compare::{compare_reports_labeled, load_report};
use loom::indexer::{IndexBuilder, alphabet::AlphabetPermuter, encrypt};
use loom::indexer::{estimate_memory, BuildMode};
use loom::search::FederatedSearcher;
use std::time::Instant;

/// Read a secret string from the environment variable named `env_var`.
fn secret_from_env(env_var: &str) -> anyhow::Result<String> {
    std::env::var(env_var)
        .map_err(|_| anyhow::anyhow!("Environment variable {} is not set", env_var))
}

/// Read a 32-byte key from the environment variable named `env_var`.
/// The value is interpreted as raw bytes (truncated/padded to 32).
fn key_from_env(env_var: &str) -> anyhow::Result<[u8; 32]> {
    let val = secret_from_env(env_var)?;
    // Derive key: treat the env var value as a password and use a fixed zero salt.
    // For production use, store a proper salt; here we use Argon2id with a zero salt
    // as a deterministic KDF for the blind key (the password itself is the secret).
    let salt = [0u8; 32];
    let key = encrypt::derive_key(&val, &salt)?;
    Ok(key)
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Index {
            path,
            extensions,
            output,
            vocab,
            encrypt: do_encrypt,
            key_env,
            blind_key_env,
            build_mode,
            estimate_memory: do_estimate,
            shards: shards_count,
            shard_size,
        } => {
            let exts: Vec<&str> = if extensions.is_empty() {
                vec![]
            } else {
                extensions.split(',').collect()
            };

            println!("Indexing {} with extensions {:?}...", path.display(), exts);

            let start = Instant::now();
            let indexer = loom::Indexer::from_directory(&path, &exts)?;
            let corpus = indexer.corpus();

            println!(
                "Found {} files, {} bytes total",
                corpus.file_count(),
                corpus.size()
            );

            // --- Preflight RAM estimator ---
            let mem_est = estimate_memory(corpus.size() as u64);
            let effective_mode: BuildMode = match build_mode {
                BuildModeArg::InMemory => BuildMode::InMemory,
                BuildModeArg::Streaming => BuildMode::Streaming,
                BuildModeArg::Auto => {
                    eprintln!("[loom] auto build-mode: {}", mem_est.rationale);
                    mem_est.recommended_mode.clone()
                }
            };

            if do_estimate {
                println!("{}", serde_json::to_string_pretty(&mem_est)?);
                return Ok(());
            }

            // --- Sharded index build path ---
            if let Some(_n_shards) = shards_count {
                let shard_dir = output.with_extension("shards");
                println!(
                    "Building sharded index (target {} bytes/shard) → {}",
                    shard_size,
                    shard_dir.display()
                );
                let manifest = loom::indexer::build_shards(corpus, &shard_dir, shard_size)?;
                let elapsed = start.elapsed();
                println!(
                    "Done! {} shards, {} files, {:.2}s",
                    manifest.shards.len(),
                    manifest.total_files,
                    elapsed.as_secs_f64()
                );
                println!(
                    "Manifest: {}",
                    shard_dir.join("shard_manifest.json").display()
                );
                return Ok(());
            }

            println!("[loom] build-mode: {}", effective_mode);
            if effective_mode == BuildMode::Streaming {
                eprintln!(
                    "[loom] WARNING: streaming mode selected but the current indexer path \
                     loads the corpus into memory. For corpora larger than available RAM, \
                     use the streaming-corpus-ingest path (dna-corpus command) and rebuild."
                );
            }
            // ---------------------------------

            // If blind-index key is requested, permute the corpus before building the FM-index.
            let builder = if let Some(ref bke) = blind_key_env {
                let bkey = key_from_env(bke)?;
                let permuter = AlphabetPermuter::new(&bkey);
                let mut permuted = corpus.concatenated.clone();
                permuter.permute_bytes(&mut permuted);

                let pb = ProgressBar::new_spinner();
                pb.set_style(ProgressStyle::default_spinner().template("{spinner:.green} {msg}")?);
                pb.set_message("Building blind FM-index (permuted corpus)...");
                pb.enable_steady_tick(std::time::Duration::from_millis(100));

                // Build a corpus with permuted text by temporarily mutating corpus data.
                let blind_corpus = loom::indexer::Corpus::from_permuted(corpus, permuted);
                let b = IndexBuilder::new_blind(&blind_corpus)?;
                pb.finish_with_message("Blind FM-index built!");
                b
            } else {
                let pb = ProgressBar::new_spinner();
                pb.set_style(ProgressStyle::default_spinner().template("{spinner:.green} {msg}")?);
                pb.set_message("Building FM-index...");
                pb.enable_steady_tick(std::time::Duration::from_millis(100));

                let b = IndexBuilder::new_with_vocab(corpus, vocab)?;
                pb.finish_with_message(if vocab { "FM-index + vocab built!" } else { "FM-index built!" });
                b
            };

            // Determine output extension.
            let idx_path = if blind_key_env.is_some() {
                output.with_extension("blind")
            } else {
                output.clone()
            };

            println!("Saving index to {}...", idx_path.display());
            builder.save(&idx_path)?;

            // Optionally encrypt.
            if do_encrypt {
                let password = secret_from_env(&key_env)?;
                encrypt::encrypt_idx_with_password(&idx_path, &password)?;
                println!("Encrypted index written to {}.enc", idx_path.display());
            }

            let elapsed = start.elapsed();
            println!(
                "Done! Indexed {} files in {:.2}s",
                builder.file_count(),
                elapsed.as_secs_f64()
            );
            if vocab {
                if let Some(v) = &builder.vocab {
                    println!("Vocab section: {} unique terms", v.len());
                }
            }

            Ok(())
        }

        Commands::Search {
            pattern,
            index,
            context,
            count,
            key_env,
            blind_key_env,
        } => {
            let start = Instant::now();

            // Auto-detect encrypted format by file extension.
            let is_enc = index.extension().and_then(|e| e.to_str()) == Some("enc");
            let is_blind_ext = index.extension().and_then(|e| e.to_str()) == Some("blind");

            println!("Loading index from {}...", index.display());
            let builder = if is_enc {
                let password = secret_from_env(&key_env)?;
                encrypt::load_encrypted_with_password(&index, &password)?
            } else {
                IndexBuilder::load(&index)?
            };
            println!("Index loaded ({} files, {} bytes)", builder.file_count(), builder.corpus_size());
            let is_blind = is_blind_ext || builder.is_blind();

            // If blind index, permute the search pattern before querying.
            let effective_pattern: Vec<u8> = if is_blind {
                let bkey = key_from_env(
                    blind_key_env
                        .as_deref()
                        .unwrap_or("LOOM_BLIND_KEY"),
                )?;
                loom::indexer::alphabet::Trapdoor::permute_term(&pattern, &bkey)
            } else {
                pattern.as_bytes().to_vec()
            };

            if count {
                let n = builder.count_bytes(&effective_pattern);
                println!("'{}' occurs {} times", pattern, n);
            } else {
                let positions = builder.search_bytes(&effective_pattern);
                println!("Found {} occurrences of '{}':\n", positions.len(), pattern);

                for (i, pos) in positions.iter().take(20).enumerate() {
                    if let Some((file, line, col)) = builder.line_col_at_offset(*pos) {
                        println!("{}. {}:{}:{}", i + 1, file, line, col);

                        if context > 0 {
                            let ctx = builder.context_at_offset(*pos, context * 40, context * 40);
                            let ctx_lines: Vec<&str> = ctx.lines().collect();
                            for line in ctx_lines.iter().take(context * 2 + 1) {
                                println!("   {}", line);
                            }
                            println!();
                        }
                    }
                }

                if positions.len() > 20 {
                    println!("... and {} more results", positions.len() - 20);
                }
            }

            let elapsed = start.elapsed();
            println!("\nSearch completed in {:.3}s", elapsed.as_secs_f64());
            Ok(())
        }

        Commands::Ask { question, index: _ } => {
            println!("Question: {}", question);
            println!("(LLM integration not yet implemented)");
            Ok(())
        }

        Commands::BatchSearch { index, patterns } => {
            let start = Instant::now();
            eprintln!("Loading index from {}...", index.display());
            let builder = IndexBuilder::load(&index)?;
            eprintln!("Index loaded ({} bytes) in {:.2}s",
                      builder.corpus_size(), start.elapsed().as_secs_f64());

            let content = std::fs::read_to_string(&patterns)?;
            let lines: Vec<&str> = content.lines()
                .map(|l| l.trim())
                .filter(|l| !l.is_empty() && !l.starts_with('#'))
                .collect();
            eprintln!("Searching {} patterns...", lines.len());
            println!("pattern,count");
            for pattern in &lines {
                let n = builder.count_bytes(pattern.as_bytes());
                println!("{},{}", pattern, n);
            }
            eprintln!("Batch search completed in {:.3}s", start.elapsed().as_secs_f64());
            Ok(())
        }

        Commands::ShardSearch {
            pattern,
            manifest,
            context,
            count: count_only,
            workers,
        } => {
            let searcher = FederatedSearcher::from_manifest(&manifest, workers)?;

            if count_only {
                let (total, stats) = searcher.count(&pattern);
                println!(
                    "'{}' occurs {} times across {} shards ({:.2}ms)",
                    pattern, total, stats.shards_searched, stats.elapsed_ms
                );
            } else {
                let (matches, stats) = searcher.search(&pattern, context);
                println!(
                    "Found {} occurrences of '{}' across {} shards ({:.2}ms)\n",
                    matches.len(), pattern, stats.shards_searched, stats.elapsed_ms
                );
                for (i, m) in matches.iter().take(20).enumerate() {
                    let loc = match (&m.file, m.line, m.col) {
                        (Some(f), Some(l), Some(c)) => format!("{}:{}:{}", f, l, c),
                        (Some(f), _, _) => f.clone(),
                        _ => format!("shard {} offset {}", m.shard_id, m.position),
                    };
                    println!("{}. [shard {}] {}", i + 1, m.shard_id, loc);
                    if !m.context.is_empty() {
                        for line in m.context.lines().take(3) {
                            println!("   {}", line);
                        }
                        println!();
                    }
                }
                if matches.len() > 20 {
                    println!("... and {} more results", matches.len() - 20);
                }
            }
            Ok(())
        }

        Commands::VocabLookup { prefix, index, exact } => {
            let builder = IndexBuilder::load(&index)?;
            let vocab = builder.vocab.as_ref().ok_or_else(|| {
                anyhow::anyhow!("Index has no vocab section. Rebuild with `loom index --vocab`.")
            })?;

            if exact {
                if vocab.exists(&prefix) {
                    println!("✓ '{}' exists in vocabulary ({} terms)", prefix, vocab.len());
                } else {
                    println!("✗ '{}' not found in vocabulary", prefix);
                }
            } else {
                let matches = vocab.prefix_search(&prefix, 20);
                if matches.is_empty() {
                    println!("No terms found with prefix '{}'", prefix);
                } else {
                    println!("Terms matching '{}' ({} shown, {} total in vocab):", prefix, matches.len(), vocab.len());
                    for term in &matches {
                        println!("  {}", term);
                    }
                }
            }
            Ok(())
        }

        Commands::Related { term, index, limit } => {
            let start = Instant::now();
            println!("Loading index from {}...", index.display());
            let builder = IndexBuilder::load(&index)?;
            println!("Index loaded ({} files, {} bytes)", builder.file_count(), builder.corpus_size());

            let elapsed = start.elapsed();

            if let Some(onto) = builder.ontology() {
                println!("Using embedded ontology ({} triples)", onto.len());
                println!();
                let related = onto.related(&term, limit);
                if related.is_empty() {
                    println!("No related terms found for '{}'", term);
                } else {
                    println!("Top {} terms related to '{}':", related.len(), term);
                    println!("{:<30} {:>8} {:>6}", "Term", "Weight", "Scope");
                    println!("{}", "-".repeat(46));
                    for (triple, other) in &related {
                        let scope_label = if triple.scope == 1 { "sect" } else { "file" };
                        println!("{:<30} {:>8.1} {:>6}", other, triple.weight, scope_label);
                    }
                }
            } else {
                println!("No embedded ontology; building on the fly...");
                let files = builder.extract_files_for_ontology();
                if files.is_empty() {
                    println!("No files found in index");
                    return Ok(());
                }
                let config = loom::indexer::ontology::OntologyConfig::default();
                let ontology = loom::indexer::ontology::build_ontology(
                    &files.iter().map(|(p, c)| (p.as_str(), c.as_str())).collect::<Vec<_>>(),
                    &config,
                );
                println!("Ontology: {}", ontology.stats());
                println!();
                let related = ontology.related(&term, limit);
                if related.is_empty() {
                    println!("No related terms found for '{}'", term);
                } else {
                    println!("Top {} terms related to '{}':", related.len(), term);
                    println!("{:<30} {:>8} {:>8} {:>8}", "Term", "Weight", "Section", "File");
                    println!("{}", "-".repeat(58));
                    for (edge, other) in &related {
                        println!("{:<30} {:>8.1} {:>8} {:>8}", other, edge.weight, edge.section_count, edge.file_count);
                    }
                }
            }

            println!("\nCompleted in {:.3}s", elapsed.as_secs_f64());
            Ok(())
        }

        Commands::DnaCorpus {
            input,
            output,
            min_len,
        } => {
            let start = Instant::now();
            println!(
                "Building normalized DNA corpus from {} -> {} (min_len={})",
                input.display(),
                output.display(),
                min_len
            );
            let written = build_dna_corpus_jsonl(&input, &output, min_len)?;
            println!("Wrote {} DNA sequence records", written);
            println!("Completed in {:.3}s", start.elapsed().as_secs_f64());
            Ok(())
        }

        Commands::DnaManifest {
            input,
            output,
            split,
            source_url,
            license,
        } => {
            let start = Instant::now();
            println!(
                "Building DNA manifest from {} -> {} (split={})",
                input.display(),
                output.display(),
                split
            );

            let manifest = build_dna_manifest(
                &input,
                &split,
                source_url.as_deref(),
                license.as_deref(),
            )?;
            write_dna_manifest_json(&output, &manifest)?;

            println!(
                "Manifest written: {} files, {} bytes",
                manifest.total_files, manifest.total_bytes
            );
            println!("Completed in {:.3}s", start.elapsed().as_secs_f64());
            Ok(())
        }

        Commands::DnaLatency {
            input,
            min_len,
            rounds,
            output,
        } => {
            use std::fs;
            use std::io::Write as _;

            let rounds = rounds.max(1);
            println!(
                "Running DNA parse latency benchmark: {} rounds on {}",
                rounds,
                input.display()
            );

            let mut timings: Vec<f64> = Vec::with_capacity(rounds);
            for round in 0..rounds {
                let (result, elapsed_ms) =
                    time_call(|| build_dna_corpus_jsonl(&input, &std::path::PathBuf::from("/dev/null"), min_len));
                match result {
                    Ok(n) => {
                        if round == 0 {
                            println!("  cold round: {:.3} ms ({} records)", elapsed_ms, n);
                        }
                    }
                    Err(e) => {
                        eprintln!("Parse error on round {}: {}", round, e);
                    }
                }
                timings.push(elapsed_ms);
            }

            match compute_latency_stats(&timings) {
                None => eprintln!("No timing data collected."),
                Some(stats) => {
                    let json = serde_json::to_string_pretty(&stats)?;
                    if let Some(out_path) = output {
                        if let Some(parent) = out_path.parent() {
                            fs::create_dir_all(parent)?;
                        }
                        let mut f = fs::File::create(&out_path)?;
                        f.write_all(json.as_bytes())?;
                        println!("Latency stats written to {}", out_path.display());
                    } else {
                        println!("{}", json);
                    }
                }
            }

            Ok(())
        }

        Commands::DnaBenchmark {
            corpus,
            k,
            seed,
            n_queries,
            query_len,
            min_seq_len,
            output,
        } => {
            use std::fs;
            use std::io::Write as _;

            let cfg = DnaBenchmarkConfig {
                k,
                seed,
                n_queries,
                query_len,
                min_seq_len,
            };

            println!(
                "Running DNA plain-index benchmark on {} (n_queries={}, k={}, query_len={}, seed={})...",
                corpus.display(),
                n_queries,
                k,
                query_len,
                seed,
            );

            let report = run_dna_benchmark(&corpus, &cfg)?;

            println!(
                "  corpus_size={} n_queries={} hit@{}={:.4} precision@{}={:.4} recall@{}={:.4} ndcg@{}={:.4}",
                report.corpus_size,
                report.n_queries,
                k, report.hit_at_k,
                k, report.precision_at_k,
                k, report.recall_at_k,
                k, report.ndcg_at_k,
            );
            println!(
                "  latency: cold={:.3}ms warm_median={:.3}ms p95={:.3}ms p99={:.3}ms",
                report.latency.cold_ms,
                report.latency.warm_median_ms,
                report.latency.p95_ms,
                report.latency.p99_ms,
            );

            let json = serde_json::to_string_pretty(&report)?;
            if let Some(out_path) = output {
                if let Some(parent) = out_path.parent() {
                    fs::create_dir_all(parent)?;
                }
                let mut f = fs::File::create(&out_path)?;
                f.write_all(json.as_bytes())?;
                println!("Benchmark report written to {}", out_path.display());
            } else {
                println!("{}", json);
            }

            Ok(())
        }

        Commands::DnaCompare {
            baseline,
            baseline_label,
            candidate,
            candidate_label,
            output,
            markdown,
        } => {
            use std::fs;
            use std::io::Write as _;

            println!(
                "Comparing '{}' ({}) vs '{}' ({})",
                baseline_label,
                baseline.display(),
                candidate_label,
                candidate.display(),
            );

            let baseline_report = load_report(&baseline)?;
            let candidate_report = load_report(&candidate)?;

            let cmp = compare_reports_labeled(
                &baseline_report,
                &baseline_label,
                &candidate_report,
                &candidate_label,
            );

            println!("  Verdict: {}", cmp.verdict);
            println!(
                "  hit@k: baseline={:.4} candidate={:.4} delta={:+.4}",
                cmp.hit_at_k.baseline, cmp.hit_at_k.candidate, cmp.hit_at_k.delta
            );
            println!(
                "  recall@k: baseline={:.4} candidate={:.4} delta={:+.4}",
                cmp.recall_at_k.baseline, cmp.recall_at_k.candidate, cmp.recall_at_k.delta
            );
            println!(
                "  p95 latency: baseline={:.3}ms candidate={:.3}ms delta={:+.3}ms",
                cmp.p95_latency_ms.baseline, cmp.p95_latency_ms.candidate, cmp.p95_latency_ms.delta
            );

            let json = serde_json::to_string_pretty(&cmp)?;
            if let Some(out_path) = output {
                if let Some(parent) = out_path.parent() {
                    fs::create_dir_all(parent)?;
                }
                let mut f = fs::File::create(&out_path)?;
                f.write_all(json.as_bytes())?;
                println!("Comparison JSON written to {}", out_path.display());
            } else {
                println!("{}", json);
            }

            if let Some(md_path) = markdown {
                if let Some(parent) = md_path.parent() {
                    fs::create_dir_all(parent)?;
                }
                let mut f = fs::File::create(&md_path)?;
                f.write_all(cmp.render_markdown().as_bytes())?;
                println!("Comparison Markdown written to {}", md_path.display());
            }

            Ok(())
        }

        Commands::CrisprScan {
            input,
            output,
            kmer_size,
            max_targets,
            name,
        } => {
            loom::dna::crispr_scan::run_scan(&input, &output, kmer_size, max_targets, &name)?;
            Ok(())
        }

        Commands::OffTargetScan {
            guides,
            reference,
            output,
            max_mismatches,
            max_guides,
            chunk_size,
        } => {
            anyhow::ensure!(
                max_mismatches <= 3,
                "--max-mismatches must be 1–3 (got {})",
                max_mismatches
            );
            let cfg = loom::dna::offtarget::OffTargetConfig {
                max_mismatches,
                chunk_size,
                max_guides,
            };
            loom::dna::offtarget::run_offtarget_scan(&guides, &reference, &output, &cfg)?;
            Ok(())
        }

        Commands::GuideConservation {
            fasta,
            guides,
            output,
        } => {
            loom::dna::guide_conservation::run(&fasta, &guides, &output)?;
            Ok(())
        }

        Commands::GuideNovelty {
            novel_targets,
            published_guides,
            output,
            max_mismatches,
        } => {
            loom::dna::guide_novelty::run(&novel_targets, &published_guides, &output, max_mismatches)?;
            Ok(())
        }

        Commands::OpsProfile { output } => {
            use std::fs;
            use std::io::Write as _;

            let profile = loom::indexer::ops_profile();
            let json = serde_json::to_string_pretty(&profile)?;

            if let Some(out_path) = output {
                if let Some(parent) = out_path.parent() {
                    fs::create_dir_all(parent)?;
                }
                let mut f = fs::File::create(&out_path)?;
                f.write_all(json.as_bytes())?;
                println!("Ops SLO profile written to {}", out_path.display());
            } else {
                println!("{}", json);
            }
            Ok(())
        }

        Commands::PrepareWasm { index, output } => {
            let start = Instant::now();
            println!("Loading index from {}...", index.display());
            let mut builder = IndexBuilder::load(&index)?;

            println!("Serializing FM-index into FMIX section...");
            builder.prepare_wasm()?;
            let out_path = output.as_ref().unwrap_or(&index);
            builder.save(out_path)?;
            println!("Saved WASM-ready index to {}", out_path.display());

            println!("Completed in {:.3}s", start.elapsed().as_secs_f64());
            Ok(())
        }

        Commands::PrepareWasmShards { index, output_dir, target_shard_bytes } => {
            let start = Instant::now();
            println!("Loading index from {}...", index.display());
            let builder = IndexBuilder::load(&index)?;

            println!("Splitting into shards (target {} bytes each)...", target_shard_bytes);
            let manifest = builder.prepare_wasm_shards(&output_dir, target_shard_bytes)?;

            println!("Created {} shards in {}:", manifest.shards.len(), output_dir.display());
            for s in &manifest.shards {
                println!("  {} — {} files, {} bytes, corpus [{}, {})",
                    s.path, s.file_count, s.size_bytes, s.corpus_start, s.corpus_end);
            }
            println!("Completed in {:.3}s", start.elapsed().as_secs_f64());
            Ok(())
        }

        Commands::PubmedScan {
            only,
            no_merge,
            no_cache,
            output,
            ontology,
        } => {
            let repo_root = std::env::current_dir()?;
            let ontology_path = ontology.unwrap_or_else(|| repo_root.join("web/data/ontology-enrichment.json"));
            let cache_dir = repo_root.join("target/pubmed_cache_v4");
            let audit_path = repo_root.join("web/data/gap-audit-v4.txt");

            loom::dna::pubmed_scan::run(
                only.as_deref(),
                no_merge,
                no_cache,
                output.as_deref(),
                &ontology_path,
                &cache_dir,
                &audit_path,
            )?;
            Ok(())
        }
    }
}
