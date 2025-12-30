//! Translation and frame generation for TBLASTX
//!
//! This module handles translation of DNA sequences to amino acid sequences
//! and generation of all six reading frames.
//!
//! NCBI BLAST style: Each translated sequence is wrapped with sentinel bytes
//! at the beginning and end to ensure proper extension termination.
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:438-453

use bio::alphabets::dna;
use crate::utils::genetic_code::GeneticCode;
use crate::utils::matrix::aa_char_to_ncbi_index;
use super::constants::SENTINEL_BYTE;

/// A query frame containing translated sequence and metadata
#[derive(Debug, Clone)]
pub struct QueryFrame {
    /// Frame number: 1..3 for forward, -1..-3 for reverse
    pub frame: i8,
    /// Translated amino acid sequence with NCBI-style sentinel bytes.
    /// 
    /// Layout: [SENTINEL_BYTE, aa_0, aa_1, ..., aa_n-1, SENTINEL_BYTE]
    /// - aa_seq[0] = SENTINEL_BYTE (255)
    /// - aa_seq[1..len-1] = actual amino acids (encoded in NCBI matrix order 0-24)
    /// - aa_seq[len-1] = SENTINEL_BYTE (255)
    /// 
    /// When accessing amino acid at logical position i, use aa_seq[i + 1].
    /// When reporting coordinates, subtract 1 from the raw position.
    /// 
    /// Reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:438-453
    pub aa_seq: Vec<u8>,
    /// Number of actual amino acids (excluding sentinels)
    pub aa_len: usize,
    /// Original DNA sequence length
    pub orig_len: usize,
    /// SEG-masked AA intervals (logical positions, frame-specific)
    /// Seeds should not be generated from these positions
    pub seg_masks: Vec<(usize, usize)>,
}

/// Convert a codon to amino acid index in NCBI matrix order (0-24)
/// NCBI order: ARNDCQEGHILKMFPSTWYVBJZX*
fn codon_to_aa_idx(codon: &[u8], table: &GeneticCode) -> u8 {
    let aa = table.get(codon);
    aa_char_to_ncbi_index(aa)
}

/// Translate a DNA sequence to an amino acid sequence with NCBI-style sentinels.
/// 
/// NCBI BLAST places NULLB sentinel bytes at the beginning and end of each
/// translated sequence to ensure proper extension termination.
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:438-453
///   prot_seq[0] = NULLB;  // First character is sentinel
///   prot_seq[index_prot] = NULLB;  // Last character is sentinel
///
/// Returns (aa_seq_with_sentinels, actual_aa_count)
fn translate_sequence(seq: &[u8], table: &GeneticCode) -> (Vec<u8>, usize) {
    let len_aa = seq.len() / 3;
    // +2 for sentinel bytes at beginning and end
    let mut aa_seq = Vec::with_capacity(len_aa + 2);
    
    // NCBI BLAST style: sentinel at the beginning
    aa_seq.push(SENTINEL_BYTE);
    
    for chunk in seq.chunks(3) {
        if chunk.len() == 3 {
            aa_seq.push(codon_to_aa_idx(chunk, table));
        }
    }
    
    // NCBI BLAST style: sentinel at the end
    aa_seq.push(SENTINEL_BYTE);
    
    (aa_seq, len_aa)
}

/// Generate all six reading frames (3 forward + 3 reverse) for a DNA sequence
///
/// # Arguments
/// * `seq` - DNA sequence to translate
/// * `table` - Genetic code table to use for translation
///
/// # Returns
/// Vector of QueryFrame structures, one for each valid reading frame.
/// Each frame's aa_seq contains sentinel bytes at positions 0 and len-1.
pub fn generate_frames(seq: &[u8], table: &GeneticCode) -> Vec<QueryFrame> {
    let mut frames = Vec::with_capacity(6);
    let rev_seq = dna::revcomp(seq);
    let seq_len = seq.len();

    // Forward (Frames 1,2,3)
    for i in 0..3 {
        if i + 3 <= seq_len {
            let (aa_seq, aa_len) = translate_sequence(&seq[i..], table);
            frames.push(QueryFrame {
                frame: (i as i8) + 1,
                aa_seq,
                aa_len,
                orig_len: seq_len,
                seg_masks: Vec::new(),
            });
        }
    }
    // Reverse (Frames -1,-2,-3)
    for i in 0..3 {
        if i + 3 <= rev_seq.len() {
            let (aa_seq, aa_len) = translate_sequence(&rev_seq[i..], table);
            frames.push(QueryFrame {
                frame: -((i as i8) + 1),
                aa_seq,
                aa_len,
                orig_len: seq_len,
                seg_masks: Vec::new(),
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
        
        // Check that frames have valid frame numbers and sentinel structure
        for frame in &frames {
            assert!(frame.frame >= -3 && frame.frame <= 3 && frame.frame != 0);
            assert_eq!(frame.orig_len, seq.len());
            // Check sentinels
            assert_eq!(frame.aa_seq[0], SENTINEL_BYTE);
            assert_eq!(frame.aa_seq[frame.aa_seq.len() - 1], SENTINEL_BYTE);
            // Check aa_len matches actual amino acids
            assert_eq!(frame.aa_len, frame.aa_seq.len() - 2);
        }
    }
    
    #[test]
    fn test_sentinel_positions() {
        let code = GeneticCode::from_id(1);
        let seq = b"ATGATGATGATG"; // 4 codons = 4 amino acids
        let frames = generate_frames(seq, &code);
        
        let frame = &frames[0]; // Frame 1
        // aa_seq should be: [SENTINEL, aa0, aa1, aa2, aa3, SENTINEL]
        assert_eq!(frame.aa_seq.len(), 6); // 4 AAs + 2 sentinels
        assert_eq!(frame.aa_len, 4);
        assert_eq!(frame.aa_seq[0], SENTINEL_BYTE);
        assert_eq!(frame.aa_seq[5], SENTINEL_BYTE);
        // Middle positions should be valid amino acids (not sentinels)
        for i in 1..5 {
            assert_ne!(frame.aa_seq[i], SENTINEL_BYTE);
        }
    }
}
