//! Alignment algorithms for BLASTN
//!
//! This module contains:
//! - Greedy alignment algorithms for high-identity sequences
//! - Gapped extension algorithms with affine gap penalties
//! - Alignment statistics calculation
//! - Utility functions (GCD, coordinate conversion, etc.)

pub mod utilities;
pub mod greedy;
pub mod gapped;
pub mod statistics;

// Re-export public APIs
pub use greedy::{
    GreedyAlignMem,
    GreedyOffset,
    greedy_align_one_direction,
    greedy_align_one_direction_ex,
};
// NCBI reference: blast_gapalign.c:3248-3305 (BlastGetOffsetsForGappedAlignment)
// NCBI reference: blast_gapalign.c:3323-3389 (BlastGetStartForGappedAlignmentNucl)
pub use gapped::{
    blast_get_offsets_for_gapped_alignment,
    blast_get_start_for_gapped_alignment_nucl,
    extend_gapped_heuristic,
    extend_gapped_one_direction,
    extend_final_traceback,
    extend_gapped_heuristic_with_traceback,
    extend_gapped_one_direction_with_traceback,
    AlnStats,
};
pub use utilities::gdb3;
pub use statistics::align_region;
