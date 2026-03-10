//! At-rest encryption for LOOM index files (Level 1).
//!
//! - AES-GCM-256 with a random 96-bit nonce per write.
//! - Key derivation via Argon2id (OWASP 2024 minimums: m=65536, t=3, p=4).
//! - Keys MUST come from environment variables, never from files or logs.
//!
//! # File format  (`.enc`)
//! ```text
//! [16-byte magic "LOOM-ENC-V1\0\0\0\0\0"]
//! [32-byte Argon2id salt]
//! [12-byte AES-GCM nonce]
//! [ciphertext + 16-byte GCM tag]
//! ```

use aes_gcm::{
    aead::{Aead, KeyInit},
    Aes256Gcm, Nonce,
};
use argon2::{Algorithm, Argon2, Params, Version};
use rand::RngCore;
use std::fs;
use std::io::Write;
use std::path::Path;

use super::IndexError;

/// Magic bytes that identify an encrypted LOOM index.
const ENC_MAGIC: &[u8; 16] = b"LOOM-ENC-V1\x00\x00\x00\x00\x00";

/// Derive a 32-byte key from a UTF-8 password and a 32-byte salt using Argon2id.
///
/// Parameters: m=65536 KiB, t=3 iterations, p=4 lanes (OWASP 2024 minimum).
pub fn derive_key(password: &str, salt: &[u8; 32]) -> Result<[u8; 32], IndexError> {
    let params = Params::new(65536, 3, 4, Some(32))
        .map_err(|e| IndexError::Serialization(format!("argon2 params: {e}")))?;
    let argon2 = Argon2::new(Algorithm::Argon2id, Version::V0x13, params);
    let mut key = [0u8; 32];
    argon2
        .hash_password_into(password.as_bytes(), salt, &mut key)
        .map_err(|e| IndexError::Serialization(format!("argon2 hash: {e}")))?;
    Ok(key)
}

/// Encrypt the raw bytes of an `.idx` file with AES-GCM-256.
///
/// Writes `<original_path>.enc` alongside the original.
/// The `key` must be exactly 32 bytes (from [`derive_key`] or direct 256-bit key).
///
/// Raw-key mode uses a zero-salt placeholder in the file header for backward
/// compatibility with historical callers that derive keys externally.
pub fn encrypt_idx(idx_path: &Path, key: &[u8; 32]) -> Result<(), IndexError> {
    let salt = [0u8; 32];
    encrypt_idx_with_key_and_salt(idx_path, key, &salt)
}

/// Encrypt the raw bytes of an `.idx` file using a password-derived key.
///
/// A fresh random salt is generated per encryption and stored in the file header.
pub fn encrypt_idx_with_password(idx_path: &Path, password: &str) -> Result<(), IndexError> {
    let mut salt = [0u8; 32];
    rand::rngs::OsRng.fill_bytes(&mut salt);
    let key = derive_key(password, &salt)?;
    encrypt_idx_with_key_and_salt(idx_path, &key, &salt)
}

fn encrypt_idx_with_key_and_salt(
    idx_path: &Path,
    key: &[u8; 32],
    salt: &[u8; 32],
) -> Result<(), IndexError> {
    let plaintext = fs::read(idx_path)?;

    // Random 96-bit (12-byte) nonce — must never repeat for the same key.
    let mut nonce_bytes = [0u8; 12];
    rand::rngs::OsRng.fill_bytes(&mut nonce_bytes);
    let nonce = Nonce::from_slice(&nonce_bytes);

    let cipher = Aes256Gcm::new_from_slice(key)
        .map_err(|e| IndexError::Serialization(format!("aes-gcm init: {e}")))?;

    let ciphertext = cipher
        .encrypt(nonce, plaintext.as_ref())
        .map_err(|e| IndexError::Serialization(format!("aes-gcm encrypt: {e}")))?;

    // Build output path: <stem>.enc
    let enc_path = idx_path.with_extension("enc");

    let mut out = Vec::with_capacity(16 + 32 + 12 + ciphertext.len());
    out.extend_from_slice(ENC_MAGIC);
    out.extend_from_slice(salt);
    out.extend_from_slice(&nonce_bytes);
    out.extend_from_slice(&ciphertext);

    fs::write(&enc_path, &out)?;
    Ok(())
}

/// Decrypt an `.enc` file back into the raw `.idx` bytes.
///
/// Returns `Err` with a clear message if the key is wrong or the file is corrupt.
fn decrypt_enc(enc_path: &Path, key: &[u8; 32]) -> Result<Vec<u8>, IndexError> {
    let data = fs::read(enc_path)?;
    decrypt_data_with_key(&data, key)
}

fn decrypt_data_with_key(data: &[u8], key: &[u8; 32]) -> Result<Vec<u8>, IndexError> {
    let (_salt, nonce_bytes, ciphertext) = parse_encrypted_payload(data)?;
    let nonce = Nonce::from_slice(nonce_bytes);

    let cipher = Aes256Gcm::new_from_slice(key)
        .map_err(|e| IndexError::Serialization(format!("aes-gcm init: {e}")))?;

    cipher
        .decrypt(nonce, ciphertext)
        .map_err(|_| IndexError::Serialization("decryption failed (wrong key or corrupt data)".to_string()))
}

fn parse_encrypted_payload(data: &[u8]) -> Result<([u8; 32], &[u8], &[u8]), IndexError> {
    if data.len() < 16 + 32 + 12 {
        return Err(IndexError::Serialization(
            "encrypted file too short".to_string(),
        ));
    }
    if &data[..16] != ENC_MAGIC {
        return Err(IndexError::Serialization(
            "not a LOOM encrypted index".to_string(),
        ));
    }

    let mut salt = [0u8; 32];
    salt.copy_from_slice(&data[16..48]);
    let nonce_bytes = &data[48..60];
    let ciphertext = &data[60..];
    Ok((salt, nonce_bytes, ciphertext))
}

/// Decrypts an `.enc` index file and loads it as an [`super::IndexBuilder`].
///
/// Transparent: the caller gets the same `IndexBuilder` as with an unencrypted load.
pub fn load_encrypted(enc_path: &Path, key: &[u8; 32]) -> Result<super::IndexBuilder, IndexError> {
    let plaintext = decrypt_enc(enc_path, key)?;
    load_index_from_plaintext_temp(&plaintext)
}

/// Decrypts an `.enc` index using a password, deriving the key from the file salt.
///
/// For backward compatibility with older files generated via raw-key mode, if
/// salt-derived decryption fails this function retries with a zero-salt derivation.
pub fn load_encrypted_with_password(
    enc_path: &Path,
    password: &str,
) -> Result<super::IndexBuilder, IndexError> {
    let data = fs::read(enc_path)?;
    let (salt, _nonce, _ciphertext) = parse_encrypted_payload(&data)?;
    let key = derive_key(password, &salt)?;

    let plaintext = match decrypt_data_with_key(&data, &key) {
        Ok(bytes) => bytes,
        Err(_) => {
            let legacy_key = derive_key(password, &[0u8; 32])?;
            decrypt_data_with_key(&data, &legacy_key)?
        }
    };
    load_index_from_plaintext_temp(&plaintext)
}

fn load_index_from_plaintext_temp(plaintext: &[u8]) -> Result<super::IndexBuilder, IndexError> {
    let mut tmp = tempfile::NamedTempFile::new().map_err(IndexError::Io)?;
    tmp.write_all(plaintext).map_err(IndexError::Io)?;
    tmp.flush().map_err(IndexError::Io)?;
    super::IndexBuilder::load(tmp.path())
}
