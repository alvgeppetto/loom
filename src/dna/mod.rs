pub mod benchmark;
pub mod compare;
pub mod corpus;
pub mod crispr_scan;
pub mod guide_conservation;
pub mod guide_novelty;
pub mod latency;
pub mod manifest;
pub mod offtarget;
pub mod pubmed_scan;

pub fn normalize_dna_sequence(input: &str) -> String {
    let mut output = String::with_capacity(input.len());
    for line in input.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('>') || trimmed.starts_with(';') {
            continue;
        }
        for base in trimmed.chars().filter_map(normalize_base) {
            output.push(base);
        }
    }
    output
}

fn normalize_base(base: char) -> Option<char> {
    match base.to_ascii_uppercase() {
        'A' | 'C' | 'G' | 'T' | 'N' => Some(base.to_ascii_uppercase()),
        'U' => Some('T'),
        'R' | 'Y' | 'S' | 'W' | 'K' | 'M' | 'B' | 'D' | 'H' | 'V' => Some('N'),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::normalize_dna_sequence;

    #[test]
    fn keeps_canonical_bases() {
        assert_eq!(normalize_dna_sequence("ACGTN"), "ACGTN");
    }

    #[test]
    fn uppercases_and_maps_u_to_t() {
        assert_eq!(normalize_dna_sequence("acguu"), "ACGTT");
    }

    #[test]
    fn maps_iupac_ambiguity_codes_to_n() {
        assert_eq!(normalize_dna_sequence("RYSWKMBDHV"), "NNNNNNNNNN");
    }

    #[test]
    fn removes_non_sequence_characters() {
        assert_eq!(normalize_dna_sequence(">chr1\nAC GT-12"), "ACGT");
    }
}
