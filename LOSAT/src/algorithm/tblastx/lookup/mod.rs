//! Lookup table construction for TBLASTX - EXACT NCBI BLAST port
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c
//!            ncbi-blast/c++/include/algo/blast/core/blast_aalookup.h
//!
//! # Module Structure
//!
//! - `backbone` - Backbone lookup table (BackboneCell, BlastAaLookupTable, build_ncbi_lookup)
//! - `neighbor_map` - Neighbor word map for neighbor_map mode

mod backbone;
mod neighbor_map;

// Re-export public types and functions
pub use backbone::{
    BackboneCell, BlastAaLookupTable, build_ncbi_lookup, build_direct_lookup,
    AA_HITS_PER_CELL,
};
pub use neighbor_map::{NeighborMap, NeighborLookup};

use crate::utils::matrix::BLASTAA_SIZE;
use crate::stats::KarlinParams;

/// NCBI BLAST parameters for protein lookup indexing.
///
/// We use the NCBI-style bit-shift + mask indexing scheme from:
/// - `blast_lookup.h`: `ComputeTableIndex` / `ComputeTableIndexIncremental`
/// - `blast_aalookup.c`: `BlastAaLookupTableNew`
///
/// NCBI uses BLASTAA_SIZE = 28 for all amino acid lookup operations.
pub const LOOKUP_WORD_LENGTH: usize = 3;
pub const LOOKUP_ALPHABET_SIZE: usize = BLASTAA_SIZE; // 28 (NCBI standard)

/// Get the charsize (bits per residue) for NCBI-style indexing
/// ilog2(28) + 1 = 5
pub const fn get_charsize() -> usize {
    5
}

/// Get the mask for k-mer indexing
pub const fn get_mask() -> usize {
    (1usize << (LOOKUP_WORD_LENGTH * get_charsize())) - 1
}

// Presence Vector - NCBI style
pub(crate) const PV_ARRAY_BTS: usize = 6;
pub(crate) const PV_ARRAY_MASK: usize = 63;

#[inline(always)]
pub(crate) fn ilog2(mut x: usize) -> usize {
    let mut r = 0usize;
    while x > 1 {
        x >>= 1;
        r += 1;
    }
    r
}

#[inline(always)]
pub(crate) fn compute_backbone_size(word_length: usize, alphabet_size: usize, charsize: usize) -> usize {
    // NCBI: for (i=0;i<word_length;i++) backbone_size |= (alphabet_size-1) << (i*charsize);
    //       backbone_size++;
    let mut backbone_size: usize = 0;
    for i in 0..word_length {
        backbone_size |= (alphabet_size - 1) << (i * charsize);
    }
    backbone_size + 1
}

#[inline(always)]
pub(crate) fn compute_mask(word_length: usize, charsize: usize) -> usize {
    (1usize << (word_length * charsize)) - 1
}

#[inline(always)]
pub(crate) fn encode_kmer_3(aa0: usize, aa1: usize, aa2: usize, charsize: usize) -> usize {
    (aa0 << (2 * charsize)) | (aa1 << charsize) | aa2
}

/// Encode a 3-mer starting at `pos` in `seq` using the NCBI-style bit-shift index.
/// Uses BLASTAA_SIZE = 28 alphabet.
pub fn encode_aa_kmer(seq: &[u8], pos: usize) -> Option<usize> {
    if pos + (LOOKUP_WORD_LENGTH - 1) >= seq.len() {
        return None;
    }
    let a0 = seq[pos] as usize;
    let a1 = seq[pos + 1] as usize;
    let a2 = seq[pos + 2] as usize;
    // NCBI accepts all residues 0-27
    if a0 >= LOOKUP_ALPHABET_SIZE || a1 >= LOOKUP_ALPHABET_SIZE || a2 >= LOOKUP_ALPHABET_SIZE {
        return None;
    }
    let charsize = ilog2(LOOKUP_ALPHABET_SIZE) + 1; // 5
    Some(encode_kmer_3(a0, a1, a2, charsize))
}

/// Test if a k-mer index is set in the presence vector
#[inline(always)]
pub fn pv_test(pv: &[u64], index: usize) -> bool {
    (pv[index >> PV_ARRAY_BTS] >> (index & PV_ARRAY_MASK)) & 1 != 0
}

#[inline(always)]
pub(crate) fn pv_set(pv: &mut [u64], index: usize) {
    pv[index >> PV_ARRAY_BTS] |= 1u64 << (index & PV_ARRAY_MASK);
}

/// Encode a k-mer from three amino acid indices (0-27 each).
/// Uses NCBI-style bit-shift indexing with BLASTAA_SIZE=28.
pub fn encode_kmer(aa0: usize, aa1: usize, aa2: usize) -> usize {
    let charsize = ilog2(LOOKUP_ALPHABET_SIZE) + 1; // 5
    encode_kmer_3(aa0, aa1, aa2, charsize)
}

/// Decode a k-mer index back to three amino acid indices.
pub fn decode_kmer(index: usize) -> (usize, usize, usize) {
    let charsize = ilog2(LOOKUP_ALPHABET_SIZE) + 1; // 5
    let mask = (1 << charsize) - 1;
    let aa2 = index & mask;
    let aa1 = (index >> charsize) & mask;
    let aa0 = (index >> (2 * charsize)) & mask;
    (aa0, aa1, aa2)
}

/// Query context - stores per-frame information for a query
#[derive(Clone)]
pub struct QueryContext {
    pub q_idx: u32,
    pub f_idx: u8,
    pub frame: i8,
    pub aa_seq: Vec<u8>,
    pub aa_seq_nomask: Option<Vec<u8>>,
    pub aa_len: usize,
    pub orig_len: usize,
    pub frame_base: i32,
    /// NCBI: kbp[context] - Karlin parameters for this context
    /// Reference: link_hsps.c line 750-752, 866-867
    pub karlin_params: KarlinParams,
}

/// Compute unmasked intervals from SEG mask ranges
pub(crate) fn compute_unmasked_intervals(seg_masks: &[(usize, usize)], aa_len: usize) -> Vec<(usize, usize)> {
    if seg_masks.is_empty() {
        return vec![(0, aa_len)];
    }
    let mut result = Vec::new();
    let mut pos = 0;
    let mut sorted = seg_masks.to_vec();
    sorted.sort_by_key(|m| m.0);
    for &(s, e) in &sorted {
        if pos < s {
            result.push((pos, s));
        }
        pos = pos.max(e);
    }
    if pos < aa_len {
        result.push((pos, aa_len));
    }
    result
}
