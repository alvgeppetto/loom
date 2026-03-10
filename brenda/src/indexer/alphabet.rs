//! Blind-index alphabet permutation (Level 2 searchable encryption).
//!
//! A 32-byte key → deterministic Fisher-Yates shuffle of the full byte alphabet [0, 255].
//! The CSPRNG seed is derived via HKDF-SHA256 so no raw key material is used as entropy.
//!
//! # Usage
//! ```rust,ignore
//! let p = AlphabetPermuter::new(&key);
//! let encrypted_byte = p.permute(b'r');
//! let original_byte  = p.depermute(encrypted_byte);
//! ```

use hkdf::Hkdf;
use sha2::Sha256;

/// Deterministic, reversible byte-alphabet permutation derived from a 32-byte key.
pub struct AlphabetPermuter {
    /// `forward[b]` = the permuted byte for original byte `b`.
    forward: [u8; 256],
    /// `inverse[b]` = the original byte for permuted byte `b`.
    inverse: [u8; 256],
}

impl AlphabetPermuter {
    /// Create a new permuter from a 32-byte key.
    ///
    /// Uses HKDF-SHA256 to expand the key into a deterministic byte stream,
    /// then performs Fisher-Yates shuffle over `[0, 255]`.
    pub fn new(key: &[u8; 32]) -> Self {
        // Derive 256 bytes of pseudorandom material via HKDF-SHA256.
        let hk = Hkdf::<Sha256>::new(None, key);
        let mut prk = [0u8; 256];
        hk.expand(b"LOOM-BLIND-PERMUTE-V1", &mut prk)
            .expect("HKDF expand: output length ≤ 255*32 bytes — always valid");

        // Initialise identity permutation.
        let mut forward: [u8; 256] = core::array::from_fn(|i| i as u8);

        // Fisher-Yates shuffle using the HKDF output bytes as selection indices.
        // For index i we need a value in [0, i] — use prk[i] mod (i+1) as the swap index.
        for i in (1..=255usize).rev() {
            let j = (prk[i] as usize) % (i + 1);
            forward.swap(i, j);
        }

        // Build inverse map.
        let mut inverse = [0u8; 256];
        for (orig, &perm) in forward.iter().enumerate() {
            inverse[perm as usize] = orig as u8;
        }

        Self { forward, inverse }
    }

    /// Permute a single byte.
    #[inline]
    pub fn permute(&self, b: u8) -> u8 {
        self.forward[b as usize]
    }

    /// Reverse a single permuted byte back to the original byte.
    #[inline]
    pub fn depermute(&self, b: u8) -> u8 {
        self.inverse[b as usize]
    }

    /// Permute a slice of bytes in-place.
    pub fn permute_bytes(&self, data: &mut [u8]) {
        for b in data.iter_mut() {
            *b = self.forward[*b as usize];
        }
    }

    /// Reverse-permute a slice of bytes in-place.
    pub fn depermute_bytes(&self, data: &mut [u8]) {
        for b in data.iter_mut() {
            *b = self.inverse[*b as usize];
        }
    }
}

/// Trapdoor: permute a UTF-8 query term so it can be searched against a blind index.
///
/// The permuted term is returned as a `Vec<u8>` (may not be valid UTF-8 after permutation).
pub struct Trapdoor;

impl Trapdoor {
    /// Permute every byte of `term` with the given 32-byte key.
    ///
    /// The result can be used as the pattern for FM-index backward search against a
    /// corpus that was built with the same key via [`AlphabetPermuter`].
    pub fn permute_term(term: &str, key: &[u8; 32]) -> Vec<u8> {
        let p = AlphabetPermuter::new(key);
        term.bytes().map(|b| p.permute(b)).collect()
    }
}
