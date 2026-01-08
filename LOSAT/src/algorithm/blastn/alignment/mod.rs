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
pub use gapped::{
    extend_gapped_heuristic,
    extend_gapped_one_direction,
    extend_final_traceback,
    AlnStats,
};
pub use utilities::gdb3;
pub use statistics::align_region;


