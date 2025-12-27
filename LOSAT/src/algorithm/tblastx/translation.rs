//! Translation and frame generation for TBLASTX
//!
//! This module handles translation of DNA sequences to amino acid sequences
//! and generation of all six reading frames.
//!
//! Amino acid encoding follows NCBI BLAST's order: ARNDCQEGHILKMFPSTWYVBJZX*
//! This is different from ASCII order (A-Z) to match NCBI BLAST's internal representation.

use bio::alphabets::dna;
use crate::utils::genetic_code::GeneticCode;

/// NCBI BLAST amino acid order: ARNDCQEGHILKMFPSTWYVBJZX*
/// This is the canonical order used by NCBI BLAST for amino acid encoding.
/// Reference: ncbi-blast/c++/src/util/tables/sm_blosum62.c
pub const NCBI_AA_ORDER: &[u8] = b"ARNDCQEGHILKMFPSTWYVBJZX*";

/// Number of standard amino acids in NCBI BLAST order (25: ARNDCQEGHILKMFPSTWYVBJZX*)
pub const NCBI_AA_COUNT: u8 = 25;

/// Stop codon index in NCBI BLAST order (24 = '*' position)
pub const NCBI_STOP_CODON_IDX: u8 = 24;

/// Lookup table: ASCII amino acid character -> NCBI BLAST index
/// Maps A-Z and * to indices 0-24 in NCBI order (ARNDCQEGHILKMFPSTWYVBJZX*)
/// Invalid characters map to 255
static ASCII_TO_NCBI_AA: [u8; 256] = {
    let mut table = [255u8; 256];
    // NCBI order: ARNDCQEGHILKMFPSTWYVBJZX*
    table[b'A' as usize] = 0;   // A
    table[b'R' as usize] = 1;    // R
    table[b'N' as usize] = 2;    // N
    table[b'D' as usize] = 3;    // D
    table[b'C' as usize] = 4;    // C
    table[b'Q' as usize] = 5;    // Q
    table[b'E' as usize] = 6;    // E
    table[b'G' as usize] = 7;    // G
    table[b'H' as usize] = 8;    // H
    table[b'I' as usize] = 9;    // I
    table[b'L' as usize] = 10;   // L
    table[b'K' as usize] = 11;   // K
    table[b'M' as usize] = 12;   // M
    table[b'F' as usize] = 13;   // F
    table[b'P' as usize] = 14;   // P
    table[b'S' as usize] = 15;   // S
    table[b'T' as usize] = 16;   // T
    table[b'W' as usize] = 17;   // W
    table[b'Y' as usize] = 18;   // Y
    table[b'V' as usize] = 19;   // V
    table[b'B' as usize] = 20;   // B (Asn or Asp)
    table[b'J' as usize] = 21;   // J (Leu or Ile)
    table[b'Z' as usize] = 22;   // Z (Gln or Glu)
    table[b'X' as usize] = 23;   // X (any amino acid)
    table[b'*' as usize] = 24;   // * (stop codon)
    table
};

/// A query frame containing translated sequence and metadata
#[derive(Debug, Clone)]
pub struct QueryFrame {
    /// Frame number: 1..3 for forward, -1..-3 for reverse
    pub frame: i8,
    /// Translated amino acid sequence (encoded as NCBI BLAST indices 0-24)
    /// 0-23: ARNDCQEGHILKMFPSTWYVBJZX
    /// 24: stop codon (*)
    pub aa_seq: Vec<u8>,
    /// Original DNA sequence length
    pub orig_len: usize,
}

/// Convert a codon to NCBI BLAST amino acid index (0-24)
/// 
/// # Arguments
/// * `codon` - 3-nucleotide codon
/// * `table` - Genetic code table
/// 
/// # Returns
/// NCBI BLAST amino acid index (0-24) where:
/// - 0-23: ARNDCQEGHILKMFPSTWYVBJZX
/// - 24: stop codon (*)
/// 
/// Reference: NCBI BLAST's amino acid encoding in ARNDCQEGHILKMFPSTWYVBJZX* order
fn codon_to_aa_idx(codon: &[u8], table: &GeneticCode) -> u8 {
    let aa_char = table.get(codon);
    let ncbi_idx = ASCII_TO_NCBI_AA[aa_char as usize];
    
    // If the amino acid is not in NCBI order, treat as stop codon
    if ncbi_idx == 255 {
        NCBI_STOP_CODON_IDX
    } else {
        ncbi_idx
    }
}

/// Translate a DNA sequence to an amino acid sequence using NCBI BLAST encoding
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

    #[test]
    fn test_ascii_to_ncbi_mapping() {
        // Test that ASCII characters map correctly to NCBI indices
        assert_eq!(ASCII_TO_NCBI_AA[b'A' as usize], 0);
        assert_eq!(ASCII_TO_NCBI_AA[b'R' as usize], 1);
        assert_eq!(ASCII_TO_NCBI_AA[b'*' as usize], 24);
        assert_eq!(ASCII_TO_NCBI_AA[b'X' as usize], 23);
    }
}
