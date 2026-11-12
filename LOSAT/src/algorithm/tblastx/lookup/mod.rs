//! Lookup table construction for TBLASTX - EXACT NCBI BLAST port
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c
//!            ncbi-blast/c++/include/algo/blast/core/blast_aalookup.h
//!
//! # Module Structure
//!
//! - `backbone` - Backbone lookup table (BackboneCell, BlastAaLookupTable, build_ncbi_lookup)
//! - `compressed` - compressed amino-acid lookup primitives for `blastp-fast`

mod backbone;
// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_aalookup.h:186-267
// ```c
// /* ------------ compressed alphabet protein blast defines ---------------*/
// typedef struct BlastCompressedAaLookupTable { ... } BlastCompressedAaLookupTable;
// ```
pub mod compressed;

// Re-export public types and functions
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(crate) use backbone::x_prepare_blosum62_lookup_query_in_pool;
// TBLASTN's query set setup on the batch's pool (`tblastn/stage_d_pipeline.rs`); no NCBI
// counterpart for the re-export itself.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(crate) use backbone::x_lookup_contexts_for_profile_in_pool;
// EXPERIMENT (LOSAT_X_BPLKPFLAT): the lookup query with the contexts' Karlin blocks left to the
// caller (`blastp/x_query_blocks.rs`); no NCBI counterpart for the re-export itself.
#[cfg(test)]
pub use backbone::build_ncbi_lookup;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(crate) use backbone::XBlosum62LookupQuery;
pub use backbone::{
    BackboneCell, BackboneCell16, BlastAaLookupTable, XCell, XHit, AA_HITS_PER_CELL,
};
// LOSAT_X_TBNLUTLEAN: the summary line `main` prints at exit (re-exported by `tblastn`). No NCBI
// counterpart: it prints counters and changes no value NCBI computes.
pub use backbone::x_tbn_lut_print_summary;
pub(crate) use backbone::{
    build_ncbi_lookup_for_profile, ideal_karlin_params_of, matrix_score_bounds,
    prepare_blosum62_lookup_query_for_word_size, prepare_matrix_lookup_query_for_word_size,
    x_compressed_contexts_for_profile, x_lookup_contexts_for_profile, x_lookup_query_for_profile,
    x_lookup_table_for_profile, x_lookup_table_for_profile_in,
};

use crate::stats::KarlinParams;
use crate::utils::matrix::BLASTAA_SIZE;

/// NCBI BLAST parameters for protein lookup indexing.
///
/// We use the NCBI-style bit-shift + mask indexing scheme from:
/// - `blast_lookup.h`: `ComputeTableIndex` / `ComputeTableIndexIncremental`
/// - `blast_aalookup.c`: `BlastAaLookupTableNew`
///
/// NCBI uses BLASTAA_SIZE = 28 for all amino acid lookup operations.
// NCBI reference (NCBISTDAA alphabet includes stop codon):
// ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:115-118
// ```c
// const char NCBISTDAA_TO_AMINOACID[BLASTAA_SIZE] = {
// '-','A','B','C','D','E','F','G','H','I','K','L','M',
// 'N','P','Q','R','S','T','V','W','X','Y','Z','U','*',
// 'O', 'J'};
// ```
pub const LOOKUP_WORD_LENGTH: usize = 3;
pub const LOOKUP_ALPHABET_SIZE: usize = BLASTAA_SIZE; // 28 (NCBI standard)

// Presence Vector - NCBI style (PV_ARRAY_TYPE = Uint4, 32-bit).
// Reference: ncbi-blast/c++/include/algo/blast/core/blast_lookup.h:41-44
pub(crate) const PV_ARRAY_BTS: usize = 5;
pub(crate) const PV_ARRAY_MASK: usize = 31;

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
pub(crate) fn compute_backbone_size(
    word_length: usize,
    alphabet_size: usize,
    charsize: usize,
) -> usize {
    // NCBI: for (i=0;i<word_length;i++) backbone_size |= (alphabet_size-1) << (i*charsize);
    //       backbone_size++;
    let mut backbone_size: usize = 0;
    for i in 0..word_length {
        backbone_size |= (alphabet_size - 1) << (i * charsize);
    }
    backbone_size + 1
}

#[cfg(test)]
#[inline(always)]
pub(crate) fn encode_kmer_3(aa0: usize, aa1: usize, aa2: usize, charsize: usize) -> usize {
    (aa0 << (2 * charsize)) | (aa1 << charsize) | aa2
}

#[inline(always)]
pub(crate) fn compute_mask(word_length: usize, charsize: usize) -> usize {
    (1usize << (word_length * charsize)) - 1
}

/// Test if a k-mer index is set in the presence vector
#[cfg(test)]
#[inline(always)]
pub fn pv_test(pv: &[u32], index: usize) -> bool {
    // NCBI PV_TEST macro:
    // lookup[(index) >> PV_ARRAY_BTS] & (1U << (index & PV_ARRAY_MASK))
    // Reference: ncbi-blast/c++/include/algo/blast/core/blast_lookup.h:55-57
    (pv[index >> PV_ARRAY_BTS] & (1u32 << (index & PV_ARRAY_MASK))) != 0
}

#[inline(always)]
pub(crate) fn pv_set(pv: &mut [u32], index: usize) {
    // NCBI PV_SET macro:
    // lookup[(index) >> PV_ARRAY_BTS] |= (1U << (index & PV_ARRAY_MASK))
    // Reference: ncbi-blast/c++/include/algo/blast/core/blast_lookup.h:49-50
    pv[index >> PV_ARRAY_BTS] |= 1u32 << (index & PV_ARRAY_MASK);
}

/// Query context - stores per-frame information for a query
#[derive(Clone)]
pub struct QueryContext {
    pub q_idx: u32,
    pub frame: i8,
    pub aa_seq: Vec<u8>,
    pub aa_seq_nomask: Option<Vec<u8>>,
    pub aa_len: usize,
    pub orig_len: usize,
    pub frame_base: i32,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_query_info.h:60-69
    // ```c
    // typedef struct BlastContextInfo {
    //     Int4 query_offset;
    //     Int4 query_length;
    //     ...
    //     Boolean is_valid; /**< Determine if this context is valid or not. */
    // } BlastContextInfo;
    // ```
    pub is_valid: bool,
    /// NCBI: `kbp[context]` - Karlin parameters for this context
    /// Reference: link_hsps.c line 750-752, 866-867
    pub karlin_params: KarlinParams,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:446-453
// ```c++
//     /* create an empty backbone */
//
//     exact_backbone = (Int4 **) calloc(lookup->backbone_size, sizeof(Int4 *));
//
//     /* find all the exact matches, grouping together all offsets of identical
//        query words. The query bias is not used here, since the next stage
//        will need real offsets into the query sequence */
//
// ```
pub(crate) use backbone::{build_ncbi_lookup_from_prepared, build_ncbi_lookup_from_prepared_of};
// BLASTP's table built in place (lookup/x_lut_direct.rs); no NCBI counterpart for the re-export
// itself.
pub(crate) use backbone::{x_bp_build_ncbi_lookup, x_bp_build_ncbi_lookup_as};
// EXPERIMENT (LOSAT_X_BPLKPAHEAD): exact blastp's table built ahead of its block
// (lookup/x_lut_direct.rs `XAaSpares`); no NCBI counterpart for the re-export itself.
pub(crate) use backbone::{XAaSpares, XBpLookup};
// The comparison of two tables (lookup/x_lut_direct.rs), for blastp's unit tests with
// `parallel` threads (blast_engine/x_lookup_ahead.rs).
#[cfg(all(
    test,
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(crate) use backbone::{x_assert_same_table_as, x_same_table};
// The BLASTX table's vectors given back for the next build (LOSAT_X_BXPOOLN1:
// by searches of two threads or more; x_lut_direct.rs), and the table of either width a BLASTX
// build returns.
pub(crate) use backbone::{
    x_bx_pooled_for, x_give_back_pooled, x_give_back_pooled16, x_hits_per_cell_hundredths, XBxTable,
};

// EXPERIMENT (LOSAT_X_TBXLUT): TBLASTX's full table built in place
// (lookup/x_lut_direct.rs); no NCBI counterpart for the re-export itself.
pub(crate) use backbone::{x_tbx_build_ncbi_lookup, x_tbxlut_reach_stub};

// EXPERIMENT (PAIRLUT): the direct builder of the TBLASTX unit
// path's group tables (lookup/x_pair_lut.rs); no NCBI counterpart for the re-export itself.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:454-456
// ```c
//     BlastLookupIndexQueryExactMatches(exact_backbone, lookup->word_length,
//                                       lookup->charsize, lookup->word_length,
//                                       query, location);
// ```
pub(crate) use backbone::x_pair_lut;
// The reference builder's table of a BLASTX lookup query, for blastx::seed's tests (item SC9 of line
// NP, a test of threaded builds); no NCBI counterpart for the re-export itself.
#[cfg(all(
    test,
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(crate) use backbone::x_reference_table_from_prepared;

#[cfg(test)]
mod lookup_tests;
