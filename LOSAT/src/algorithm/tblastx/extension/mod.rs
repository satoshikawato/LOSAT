//! Extension algorithms for TBLASTX
//!
//! This module implements ungapped and gapped extension algorithms for protein sequences,
//! following NCBI BLAST's approach with X-drop termination and affine gap penalties.
//!
//! Sequences are encoded in NCBISTDAA (0-27), and scoring uses BLOSUM62 via
//! the ncbistdaa_to_blosum62 conversion table.
//!
//! # Module Structure
//!
//! - `ungapped` - Ungapped extension with SIMD optimization
//! - `two_hit` - Two-hit extension algorithm (NCBI style)
//! - `gapped` - Gapped extension with affine gap penalties

mod ungapped;
mod two_hit;
mod gapped;

// Re-export public functions
pub use ungapped::extend_hit_ungapped;
pub use two_hit::extend_hit_two_hit;
pub use gapped::extend_gapped_protein;

use crate::stats::KarlinParams;
use crate::algorithm::common::evalue::calculate_evalue_alignment_length;
use crate::utils::matrix::blosum62_score;
use super::constants::{SENTINEL_BYTE, SENTINEL_PENALTY};

/// Get the substitution matrix score for two amino acids in NCBISTDAA encoding.
///
/// If either character is a sentinel byte (SENTINEL_BYTE = 0, NCBI NULLB), returns
/// SENTINEL_PENALTY (-4) to trigger X-drop termination at sequence boundaries.
///
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:120
///   const Uint1 kProtSentinel = NULLB;
/// Reference: ncbi-blast/c++/include/algo/blast/core/ncbi_std.h:181
///   #define NULLB '\0'
#[inline(always)]
pub fn get_score(a: u8, b: u8) -> i32 {
    // Check for sentinel bytes (NCBI BLAST style sequence boundary markers)
    if a == SENTINEL_BYTE || b == SENTINEL_BYTE {
        return SENTINEL_PENALTY;
    }
    blosum62_score(a, b)
}

/// Convert amino acid coordinates to nucleotide coordinates.
///
/// NCBI reference (link_hsps.c:1247-1262):
/// For positive frames:
///   start = CODON_LENGTH*(segment->offset) + segment->frame - 1;
///   end = CODON_LENGTH*segment->end + segment->frame - 2;
/// For negative frames:
///   start = seq_length - CODON_LENGTH*segment->offset + segment->frame;
///   end = seq_length - CODON_LENGTH*segment->end + segment->frame + 1;
pub fn convert_coords(aa_start: usize, aa_end: usize, frame: i8, dna_len: usize) -> (usize, usize) {
    let f_abs = frame.abs() as usize;
    let shift = f_abs - 1;

    if frame > 0 {
        // NCBI: start = 3*offset + frame - 1 = 3*offset + shift
        // NCBI: end = 3*end + frame - 2 = 3*end + shift - 1
        let start_bp = aa_start * 3 + shift;
        let end_bp = aa_end * 3 + shift - 1;
        (start_bp, end_bp)
    } else {
        // NCBI: start = len - 3*offset + frame = len - 3*offset - f_abs
        //             = len - (3*offset + f_abs) = len - (3*offset + shift + 1)
        // NCBI: end = len - 3*end + frame + 1 = len - 3*end - f_abs + 1
        //           = len - (3*end + f_abs - 1) = len - (3*end + shift)
        let start_bp = dna_len - (aa_start * 3 + shift + 1);
        let end_bp_calc = dna_len - (aa_end * 3 + shift);
        (start_bp, end_bp_calc)
    }
}

/// Calculate bit score and E-value for a protein alignment
///
/// This function uses alignment length as the effective search space,
/// which is appropriate for protein alignments.
pub fn calculate_statistics(score: i32, aln_len: usize, params: &KarlinParams) -> (f64, f64) {
    calculate_evalue_alignment_length(score, aln_len, params)
}
