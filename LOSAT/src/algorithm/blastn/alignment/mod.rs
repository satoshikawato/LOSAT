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

// Re-export public APIs
pub use greedy::{
    greedy_align_one_direction,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2762-2936 (BLAST_GreedyGappedAlignment)
    greedy_gapped_alignment_with_traceback,
    GreedyAlignMem,
    GreedyAlignScratch,
    GreedyOffset,
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
    extend_gapped_heuristic,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    extend_gapped_heuristic_with_scratch,
    extend_gapped_one_direction,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    GapAlignScratch,
};
pub(crate) use gapped::{extend_gapped_heuristic_with_traceback_with_scratch, stats_from_edit_ops};
pub use statistics::align_region;
pub use utilities::gdb3;
