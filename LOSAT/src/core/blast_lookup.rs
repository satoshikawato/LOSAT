//! Common Lookup Table Infrastructure
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c
//!
//! This module provides common lookup table infrastructure shared between
//! amino acid (blast_aalookup) and nucleotide (blast_nalookup) lookup tables.
//!
//! During migration, this module re-exports commonly used types from both
//! AA and NA lookup modules.

// Common constants
pub const PV_ARRAY_BTS: usize = 6;
pub const PV_ARRAY_MASK: usize = 63;
pub const PV_BUCKET_BITS: usize = 64;

/// Test if a bit is set in the presence vector
#[inline(always)]
pub fn pv_test(pv: &[u64], index: usize) -> bool {
    (pv[index >> PV_ARRAY_BTS] & (1u64 << (index & PV_ARRAY_MASK))) != 0
}

/// Set a bit in the presence vector
#[inline(always)]
pub fn pv_set(pv: &mut [u64], index: usize) {
    pv[index >> PV_ARRAY_BTS] |= 1u64 << (index & PV_ARRAY_MASK);
}

/// Calculate the size of a presence vector array for a given backbone size
#[inline]
pub fn pv_array_size(backbone_size: usize) -> usize {
    (backbone_size + PV_BUCKET_BITS - 1) / PV_BUCKET_BITS
}
