//! NCBI-style sum-statistics HSP linking with lh_helper optimization
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c
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
//! - `linking` - Main linking algorithm (apply_sum_stats_even_gap_linking_with_parallel)
//! - `linking_incr` - The per-group kernel `linking` runs by default (incremental recomputation)

mod cutoffs;
mod linking;
mod linking_incr;
mod params;

// Re-export parameter types and functions
pub use params::{compute_avg_query_length_ncbi, find_smallest_lambda_params, LinkingParams};

// Re-export cutoff calculation
pub use cutoffs::calculate_link_hsp_cutoffs_ncbi;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
// Re-export main linking functions
pub use linking::apply_sum_stats_even_gap_linking_with_parallel;
pub(crate) use linking::link_one_hit_with_int4_ends;

// The in-place HSP list path
// (tblastx/x_hsp_inplace.rs, a child module of `linking`). `ungapped_hit_sort.rs` calls its sort
// dispatch, `main` its summary.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:485-486
pub(crate) use linking::x_hsp_inplace;
pub use linking::x_hsp_inplace_print_summary;

// The summary of the index sorts of the
// in-place path (tblastx/x_link_sort.rs, a child module of `linking`), for `main`.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:993-994
// No NCBI counterpart for the summary: it prints counters.
pub use linking::x_link_sort_print_summary;

// The summary of the groups linked as
// sub-slices of the sorted list (tblastx/x_link_slices.rs, a child module of `linking`), for
// `main`.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
// No NCBI counterpart for the summary: it prints counters.
pub use linking::x_link_slices_print_summary;
pub(crate) use linking::{x_link_helped, XHelped};
