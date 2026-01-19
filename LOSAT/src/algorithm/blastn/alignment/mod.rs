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
    GreedyAlignScratch,
    GreedyOffset,
    greedy_align_one_direction,
    greedy_align_one_direction_ex,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2762-2936 (BLAST_GreedyGappedAlignment)
    greedy_gapped_alignment_score_only,
    greedy_gapped_alignment_with_traceback,
};
// NCBI reference: blast_gapalign.c:3248-3305 (BlastGetOffsetsForGappedAlignment)
// NCBI reference: blast_gapalign.c:3323-3389 (BlastGetStartForGappedAlignmentNucl)
pub use gapped::{
    build_blastna_matrix,
    blast_get_offsets_for_gapped_alignment,
    blast_get_start_for_gapped_alignment_nucl,
    extend_gapped_heuristic,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    extend_gapped_heuristic_with_scratch,
    extend_gapped_one_direction,
    extend_gapped_heuristic_with_traceback,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    extend_gapped_heuristic_with_traceback_with_scratch,
    extend_gapped_one_direction_with_traceback,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    GapAlignScratch,
    AlnStats,
};
pub use utilities::gdb3;
pub use statistics::align_region;
