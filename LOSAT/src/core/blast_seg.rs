//! SEG filter implementation
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1625-1638
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2282-2325
//
// Kappa subject SEG and query SEG must use the same `blast_seg.c` port.
// Re-export the canonical implementation from `utils::seg` so there is a
// single NCBI-faithful SEG code path in LOSAT.
pub use crate::utils::seg::{SegMasker, SegParams};
