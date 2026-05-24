//! NCBI-style sum-statistics HSP linking with lh_helper optimization
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/link_hsps.c
//!
//! Key optimizations from NCBI:
//! 1. `changed` flag: skip recomputation when previous link is still valid
//! 2. `linked_to` counter: track how many HSPs link to this one
//! 3. `path_changed` flag: skip recomputation when no chains were affected
//! 4. `use_current_max`: reuse max chains when they weren't affected by removal
//! 5. `next_larger`: skip HSPs with too-small sum
//! 6. Dual-index approach: index=0 (small gaps), index=1 (large gaps)
//!
//! # Module Structure
//!
//! - `params` - Parameter structures (LinkingParams, LinkHspCutoffs) and helpers
//! - `cutoffs` - NCBI cutoff calculation (calculate_link_hsp_cutoffs_ncbi)
//! - `linking` - Main linking algorithm (apply_sum_stats_even_gap_linking)

mod cutoffs;
mod linking;
mod params;

// Re-export parameter types and functions
pub use params::{
    compute_avg_query_length_ncbi, find_smallest_lambda, find_smallest_lambda_params,
    LinkHspCutoffs, LinkingParams, BLAST_GAP_DECAY_RATE,
};

// Re-export cutoff calculation
pub use cutoffs::calculate_link_hsp_cutoffs_ncbi;

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/link_hsps.c:553-558
// ```c
// for (frame_index=0; frame_index<num_query_frames; frame_index++)
// {
//    hp_start->next = hp_frame_start[frame_index];
//    number_of_hsps = hp_frame_number[frame_index];
// }
// ```
// Re-export main linking functions
pub use linking::{
    apply_sum_stats_even_gap_linking, apply_sum_stats_even_gap_linking_with_parallel,
};
