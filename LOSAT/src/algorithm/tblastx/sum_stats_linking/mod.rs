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

mod params;
mod cutoffs;
mod linking;

// Re-export parameter types and functions
pub use params::{
    LinkHspCutoffs,
    LinkingParams,
    BLAST_GAP_DECAY_RATE,
    find_smallest_lambda,
    find_smallest_lambda_params,
    compute_avg_query_length_ncbi,
};

// Re-export cutoff calculation
pub use cutoffs::calculate_link_hsp_cutoffs_ncbi;

// Re-export main linking function
pub use linking::apply_sum_stats_even_gap_linking;
