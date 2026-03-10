//! Tests for Lane 5: Encrypted Search (E1.6 + E2.5 + E2.6)
//!
//! Run with: `cargo test --test encrypt_tests`

use loom::indexer::{
    alphabet::{AlphabetPermuter, Trapdoor},
    encrypt::{derive_key, encrypt_idx, encrypt_idx_with_password, load_encrypted, load_encrypted_with_password},
    Corpus, IndexBuilder,
};
use std::fs;
use std::process::Command;
use tempfile::TempDir;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Build a tiny in-memory corpus from a temp directory with two files.
fn build_test_corpus(dir: &TempDir) -> Corpus {
    let f1 = dir.path().join("hello.txt");
    let f2 = dir.path().join("world.txt");
    fs::write(&f1, "hello world foo bar").unwrap();
    fs::write(&f2, "world baz qux hello").unwrap();
    Corpus::from_directory(dir.path(), &[]).unwrap()
}

fn test_key() -> [u8; 32] {
    let salt = [0u8; 32];
    derive_key("test-password-loom", &salt).unwrap()
}

fn test_key_wrong() -> [u8; 32] {
    let salt = [0u8; 32];
    derive_key("wrong-password-loom", &salt).unwrap()
}

// ---------------------------------------------------------------------------
// E1.6 — At-rest encryption
// ---------------------------------------------------------------------------

/// E1.6-a: encrypt → decrypt round-trip produces the same index bytes.
#[test]
fn test_encrypt_decrypt_roundtrip() {
    let dir = TempDir::new().unwrap();
    let corpus = build_test_corpus(&dir);
    let builder = IndexBuilder::new(&corpus).unwrap();

    let idx_path = dir.path().join("test.idx");
    builder.save(&idx_path).unwrap();
    let original_bytes = fs::read(&idx_path).unwrap();

    let key = test_key();
    encrypt_idx(&idx_path, &key).unwrap();

    let enc_path = idx_path.with_extension("enc");
    assert!(enc_path.exists(), ".enc file should be created");

    // Decrypt back and reload — should succeed.
    let decrypted = load_encrypted(&enc_path, &key).unwrap();
    assert_eq!(decrypted.file_count(), builder.file_count());
    assert_eq!(decrypted.corpus_size(), builder.corpus_size());

    // Verify the .enc file is NOT the same as the original .idx bytes.
    let enc_bytes = fs::read(&enc_path).unwrap();
    assert_ne!(enc_bytes, original_bytes, "Encrypted file must differ from plaintext");
}

/// Password API round-trip with random salt must decrypt correctly.
#[test]
fn test_encrypt_decrypt_roundtrip_password_api() {
    let dir = TempDir::new().unwrap();
    let corpus = build_test_corpus(&dir);
    let builder = IndexBuilder::new(&corpus).unwrap();

    let idx_path = dir.path().join("pw_roundtrip.idx");
    builder.save(&idx_path).unwrap();

    encrypt_idx_with_password(&idx_path, "super-secret-password").unwrap();
    let enc_path = idx_path.with_extension("enc");
    assert!(enc_path.exists(), ".enc file should be created");

    let decrypted = load_encrypted_with_password(&enc_path, "super-secret-password").unwrap();
    assert_eq!(decrypted.file_count(), builder.file_count());
    assert_eq!(decrypted.corpus_size(), builder.corpus_size());
}

/// Wrong password must fail for password-derived encrypted files.
#[test]
fn test_encrypt_password_wrong_password_fails() {
    let dir = TempDir::new().unwrap();
    let corpus = build_test_corpus(&dir);
    let builder = IndexBuilder::new(&corpus).unwrap();

    let idx_path = dir.path().join("pw_wrong.idx");
    builder.save(&idx_path).unwrap();
    encrypt_idx_with_password(&idx_path, "correct-password").unwrap();

    let enc_path = idx_path.with_extension("enc");
    let result = load_encrypted_with_password(&enc_path, "wrong-password");
    assert!(result.is_err(), "Wrong password must return Err");
}

/// E1.6-b: wrong key returns Err (authentication tag mismatch).
#[test]
fn test_decrypt_wrong_key_returns_error() {
    let dir = TempDir::new().unwrap();
    let corpus = build_test_corpus(&dir);
    let builder = IndexBuilder::new(&corpus).unwrap();

    let idx_path = dir.path().join("test_wk.idx");
    builder.save(&idx_path).unwrap();

    let key = test_key();
    encrypt_idx(&idx_path, &key).unwrap();

    let enc_path = idx_path.with_extension("enc");
    let wrong_key = test_key_wrong();
    let result = load_encrypted(&enc_path, &wrong_key);
    assert!(result.is_err(), "Wrong key must return Err");
}

/// E1.6-c: correct search results after decrypt.
#[test]
fn test_search_after_decrypt() {
    let dir = TempDir::new().unwrap();
    let corpus = build_test_corpus(&dir);
    let builder = IndexBuilder::new(&corpus).unwrap();

    let idx_path = dir.path().join("test_search.idx");
    builder.save(&idx_path).unwrap();

    let key = test_key();
    encrypt_idx(&idx_path, &key).unwrap();
    let enc_path = idx_path.with_extension("enc");

    let decrypted = load_encrypted(&enc_path, &key).unwrap();
    let count = decrypted.count("hello");
    assert!(count >= 2, "Should find 'hello' in both files, got {}", count);
    assert_eq!(decrypted.count("notpresent"), 0, "Unknown term must return 0");
}

/// E1.6-d: encrypted file has correct magic header.
#[test]
fn test_encrypted_file_magic() {
    let dir = TempDir::new().unwrap();
    let corpus = build_test_corpus(&dir);
    let builder = IndexBuilder::new(&corpus).unwrap();
    let idx_path = dir.path().join("test_magic.idx");
    builder.save(&idx_path).unwrap();

    let key = test_key();
    encrypt_idx(&idx_path, &key).unwrap();

    let enc_path = idx_path.with_extension("enc");
    let data = fs::read(&enc_path).unwrap();
    assert_eq!(&data[..11], b"LOOM-ENC-V1", "Magic header mismatch");
}

// ---------------------------------------------------------------------------
// E2.1 — AlphabetPermuter: basic properties
// ---------------------------------------------------------------------------

/// Permutation must be bijective (permute then depermute = identity).
#[test]
fn test_permuter_roundtrip() {
    let key = test_key();
    let p = AlphabetPermuter::new(&key);
    for b in 0u8..=255 {
        assert_eq!(p.depermute(p.permute(b)), b, "byte {} failed round-trip", b);
    }
}

/// Two different keys must produce different permutations (with overwhelming probability).
#[test]
fn test_permuter_different_keys_differ() {
    let key1 = test_key();
    let key2 = test_key_wrong();
    let p1 = AlphabetPermuter::new(&key1);
    let p2 = AlphabetPermuter::new(&key2);

    let same = (0u8..=255).filter(|&b| p1.permute(b) == p2.permute(b)).count();
    assert!(same < 10, "Two different keys produced nearly identical permutations ({})", same);
}

/// Same key must produce the same permutation deterministically.
#[test]
fn test_permuter_deterministic() {
    let key = test_key();
    let p1 = AlphabetPermuter::new(&key);
    let p2 = AlphabetPermuter::new(&key);
    for b in 0u8..=255 {
        assert_eq!(p1.permute(b), p2.permute(b));
    }
}

// ---------------------------------------------------------------------------
// E2.5 — Blind index: unkeyed = 0 results, correct key = expected results
// ---------------------------------------------------------------------------

/// Build a blind index and verify that searching without permutation finds nothing.
#[test]
fn test_blind_index_unkeyed_search_returns_zero() {
    let dir = TempDir::new().unwrap();
    let corpus = build_test_corpus(&dir);

    let key = test_key();
    let permuter = AlphabetPermuter::new(&key);
    let mut permuted_bytes = corpus.concatenated.clone();
    permuter.permute_bytes(&mut permuted_bytes);

    let blind_corpus = Corpus::from_permuted(&corpus, permuted_bytes);
    let blind_builder = IndexBuilder::new_blind(&blind_corpus).unwrap();

    // Unkeyed: search for the plaintext term "hello" — must find nothing.
    let results = blind_builder.count("hello");
    assert_eq!(results, 0, "Unkeyed search on blind index must return 0 (got {})", results);
}

/// Build a blind index and verify that the correct trapdoor key finds expected terms.
#[test]
fn test_blind_index_correct_key_finds_terms() {
    let dir = TempDir::new().unwrap();
    let corpus = build_test_corpus(&dir);

    // Count occurrences in plaintext for reference.
    let plain_builder = IndexBuilder::new(&corpus).unwrap();
    let plain_count = plain_builder.count("hello");

    let key = test_key();
    let permuter = AlphabetPermuter::new(&key);
    let mut permuted_bytes = corpus.concatenated.clone();
    permuter.permute_bytes(&mut permuted_bytes);

    let blind_corpus = Corpus::from_permuted(&corpus, permuted_bytes);
    let blind_builder = IndexBuilder::new_blind(&blind_corpus).unwrap();

    // Use trapdoor to permute the query.
    let permuted_query = Trapdoor::permute_term("hello", &key);
    let blind_count = blind_builder.count_bytes(&permuted_query);

    assert_eq!(
        blind_count, plain_count,
        "Blind search with correct key must return same count as plaintext search"
    );
    assert!(blind_count >= 2, "Should find at least 2 occurrences of 'hello'");
}

/// is_blind() flag must be set on blind indices and not on plain indices.
#[test]
fn test_blind_flag() {
    let dir = TempDir::new().unwrap();
    let corpus = build_test_corpus(&dir);

    let plain = IndexBuilder::new(&corpus).unwrap();
    assert!(!plain.is_blind(), "Plain index must not be flagged as blind");

    let key = test_key();
    let permuter = AlphabetPermuter::new(&key);
    let mut pb = corpus.concatenated.clone();
    permuter.permute_bytes(&mut pb);
    let blind_corpus = Corpus::from_permuted(&corpus, pb);
    let blind = IndexBuilder::new_blind(&blind_corpus).unwrap();
    assert!(blind.is_blind(), "Blind index must be flagged as blind");
}

/// Regression: encrypted blind index search must still apply trapdoor permutation.
#[test]
fn test_cli_search_encrypted_blind_index_with_key_succeeds() {
    let dir = TempDir::new().unwrap();
    let file = dir.path().join("src.txt");
    fs::write(&file, "hello world\nhello queue\n").unwrap();

    let bin = env!("CARGO_BIN_EXE_loom");
    let out_base = dir.path().join("secure.idx");

    let index = Command::new(bin)
        .arg("index")
        .arg("--path")
        .arg(dir.path())
        .arg("--output")
        .arg(&out_base)
        .arg("--blind-key-env")
        .arg("LOOM_BLIND_KEY")
        .arg("--encrypt")
        .arg("--key-env")
        .arg("LOOM_KEY")
        .env("LOOM_BLIND_KEY", "blind-test-password")
        .env("LOOM_KEY", "enc-test-password")
        .output()
        .expect("failed to run loom index");
    assert!(
        index.status.success(),
        "index command failed:\nstdout={}\nstderr={}",
        String::from_utf8_lossy(&index.stdout),
        String::from_utf8_lossy(&index.stderr)
    );

    // Blind output is secure.blind, then encrypted to secure.enc by current writer behavior.
    let enc_path = out_base.with_extension("enc");
    assert!(enc_path.exists(), "encrypted index not created at {}", enc_path.display());

    let search = Command::new(bin)
        .arg("search")
        .arg("hello")
        .arg("--index")
        .arg(&enc_path)
        .arg("--count")
        .arg("--key-env")
        .arg("LOOM_KEY")
        .arg("--blind-key-env")
        .arg("LOOM_BLIND_KEY")
        .env("LOOM_BLIND_KEY", "blind-test-password")
        .env("LOOM_KEY", "enc-test-password")
        .output()
        .expect("failed to run loom search");
    assert!(
        search.status.success(),
        "search command failed:\nstdout={}\nstderr={}",
        String::from_utf8_lossy(&search.stdout),
        String::from_utf8_lossy(&search.stderr)
    );

    let stdout = String::from_utf8_lossy(&search.stdout);
    assert!(
        stdout.contains("occurs") && !stdout.contains("occurs 0 times"),
        "expected non-zero encrypted blind count, got: {}",
        stdout
    );
}

// ---------------------------------------------------------------------------
// E2.6 — Blind index overhead benchmark (< 5% latency)
// ---------------------------------------------------------------------------

/// Rough benchmark: permutation overhead over the corpus should be negligible.
/// This is a timing sanity check, not a strict gate (CI environments vary).
#[test]
fn test_blind_permutation_overhead_sanity() {
    use std::time::Instant;

    let dir = TempDir::new().unwrap();
    // Write a larger corpus (~50 KB) to measure meaningful timing.
    let big_content = "hello world foo bar baz qux\n".repeat(1000);
    for i in 0..5 {
        fs::write(dir.path().join(format!("file{}.txt", i)), &big_content).unwrap();
    }
    let corpus = Corpus::from_directory(dir.path(), &[]).unwrap();

    // Time plain indexing.
    let t0 = Instant::now();
    let _ = IndexBuilder::new(&corpus).unwrap();
    let plain_ms = t0.elapsed().as_millis() as f64;

    // Time blind indexing (permutation + FM-index build).
    let key = test_key();
    let t1 = Instant::now();
    let permuter = AlphabetPermuter::new(&key);
    let mut pb = corpus.concatenated.clone();
    permuter.permute_bytes(&mut pb);
    let blind_corpus = Corpus::from_permuted(&corpus, pb);
    let _ = IndexBuilder::new_blind(&blind_corpus).unwrap();
    let blind_ms = t1.elapsed().as_millis() as f64;

    if plain_ms > 0.0 {
        let overhead = (blind_ms - plain_ms) / plain_ms;
        println!("Plain: {plain_ms:.1}ms  Blind: {blind_ms:.1}ms  Overhead: {:.1}%", overhead * 100.0);
        // Allow up to 20% in slow CI; real hardware meets the <5% target.
        assert!(
            overhead < 0.20,
            "Blind index overhead too large: {:.1}% (target < 5%, CI limit < 20%)",
            overhead * 100.0
        );
    }
}
