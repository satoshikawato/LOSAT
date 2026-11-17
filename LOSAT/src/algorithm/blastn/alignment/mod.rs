//! Alignment algorithms for BLASTN
//!
//! This module contains:
//! - Greedy alignment algorithms for high-identity sequences
//! - Gapped extension algorithms with affine gap penalties
//! - Alignment statistics calculation
//! - Utility functions (GCD, coordinate conversion, etc.)

pub mod gapped;
pub mod greedy;
pub mod statistics;
pub mod utilities;
// No NCBI counterpart: the fences that a thread's alignments met (LOSAT_X_BNTBRANGE traceback,
// preliminary greedy windows).
pub(crate) mod x_fence;
// No NCBI counterpart: the preliminary greedy alignment on a window of the packed subject, one
// base per byte, or on a one-byte copy of the whole subject (see its doc for the NCBI references).
pub mod x_greedy_prelim;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:385-389
// ```c
//           query = query_blk->sequence +
//               query_info->contexts[hsp->context].query_offset;
//           query_nomask = query_blk->sequence_nomask +
//               query_info->contexts[hsp->context].query_offset;
//           query_length = query_info->contexts[hsp->context].query_length;
// ```
/// NCBI's query pointer of one context with its length: the letters of the concatenated query
/// `concat` (every context, with the sentinel between two of them; NCBI's `sequence` or
/// `sequence_nomask`) from the context's first letter `start`, and the context's `len`. NCBI's
/// `Int4` offsets relative to that pointer (`at`) read the concatenated query outside the context
/// too, before it for a negative offset and past its end for one at `len` or more, as NCBI's
/// pointer arithmetic does for an HSP from a hit that passes an end of its context
/// (remaining-gaps X1b-o).
#[derive(Clone, Copy, Debug)]
pub struct ContextQuery<'a> {
    pub concat: &'a [u8],
    pub start: usize,
    pub len: usize,
}

impl ContextQuery<'_> {
    /// The place in `concat` of NCBI's offset `off` from the context's start (`query + off`).
    #[inline(always)]
    pub fn at(&self, off: i32) -> usize {
        (self.start as isize + off as isize) as usize
    }
}

// Re-export public APIs
pub use greedy::{
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2762-2936 (BLAST_GreedyGappedAlignment)
    greedy_gapped_alignment_with_traceback,
    GreedyAlignScratch,
};
// NCBI reference: blast_gapalign.c:3248-3305 (BlastGetOffsetsForGappedAlignment)
// NCBI reference: blast_gapalign.c:3323-3389 (BlastGetStartForGappedAlignmentNucl)
// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
// The traceback extension reads its matrix with `get_unchecked` for letters below 16 (BLASTNA,
// which `run_local` encodes first), so it is the crate's only.
pub use gapped::{
    blast_get_offsets_for_gapped_alignment,
    blast_get_start_for_gapped_alignment_nucl,
    build_blastna_matrix,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    extend_gapped_heuristic_with_scratch,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    GapAlignScratch,
};
pub(crate) use gapped::{extend_gapped_heuristic_with_traceback_with_scratch, stats_from_edit_ops};
