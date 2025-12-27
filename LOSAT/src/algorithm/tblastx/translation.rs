//! Translation and frame generation for TBLASTX
//!
//! This module handles translation of DNA sequences to amino acid sequences
//! and generation of all six reading frames.

use bio::alphabets::dna;
use crate::utils::genetic_code::GeneticCode;

/// A query frame containing translated sequence and metadata
#[derive(Debug, Clone)]
pub struct QueryFrame {
    /// Frame number: 1..3 for forward, -1..-3 for reverse
    pub frame: i8,
    /// Translated amino acid sequence (encoded as 0-24 for A-Z)
    pub aa_seq: Vec<u8>,
    /// Original DNA sequence length
    pub orig_len: usize,
}

/// Convert a codon to amino acid index (0-24 for A-Z, 25 for stop/other)
fn codon_to_aa_idx(codon: &[u8], table: &GeneticCode) -> u8 {
    let aa = table.get(codon);
    // '*' などA-Z以外は25にマップされる
    if aa >= b'A' && aa <= b'Z' {
        aa - b'A'
    } else {
        super::constants::STOP_CODON
    }
}

/// Translate a DNA sequence to an amino acid sequence
fn translate_sequence(seq: &[u8], table: &GeneticCode) -> Vec<u8> {
    let len_aa = seq.len() / 3;
    let mut aa_seq = Vec::with_capacity(len_aa);
    for chunk in seq.chunks(3) {
        if chunk.len() == 3 {
            aa_seq.push(codon_to_aa_idx(chunk, table));
        }
    }
    aa_seq
}

/// Generate all six reading frames (3 forward + 3 reverse) for a DNA sequence
///
/// # Arguments
/// * `seq` - DNA sequence to translate
/// * `table` - Genetic code table to use for translation
///
/// # Returns
/// Vector of QueryFrame structures, one for each valid reading frame
pub fn generate_frames(seq: &[u8], table: &GeneticCode) -> Vec<QueryFrame> {
    let mut frames = Vec::with_capacity(6);
    let rev_seq = dna::revcomp(seq);
    let seq_len = seq.len();

    // Forward (Frames 1,2,3)
    for i in 0..3 {
        if i + 3 <= seq_len {
            frames.push(QueryFrame {
                frame: (i as i8) + 1,
                aa_seq: translate_sequence(&seq[i..], table),
                orig_len: seq_len,
            });
        }
    }
    // Reverse (Frames -1,-2,-3)
    for i in 0..3 {
        if i + 3 <= rev_seq.len() {
            frames.push(QueryFrame {
                frame: -((i as i8) + 1),
                aa_seq: translate_sequence(&rev_seq[i..], table),
                orig_len: seq_len,
            });
        }
    }
    frames
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_frames() {
        let code = GeneticCode::from_id(1);
        let seq = b"ATGCGATCGATCGATCG";
        let frames = generate_frames(seq, &code);
        
        // Should generate up to 6 frames
        assert!(frames.len() <= 6);
        
        // Check that frames have valid frame numbers
        for frame in &frames {
            assert!(frame.frame >= -3 && frame.frame <= 3 && frame.frame != 0);
            assert_eq!(frame.orig_len, seq.len());
        }
    }
}

