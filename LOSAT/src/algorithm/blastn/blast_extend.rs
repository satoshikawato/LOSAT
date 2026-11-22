//! Diagonal table management for BLASTN
//!
//! NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:57-60
//!
//! This module contains the DiagStruct for tracking hits on diagonals
//! during the two-hit word finding algorithm.

/// Diagonal structure for tracking hits (NCBI DiagStruct equivalent)
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:57-60
#[derive(Clone, Copy, Default)]
pub struct DiagStruct {
    /// Last hit position on this diagonal (with diag_offset added)
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:96-105
    /// LOSAT: Initialize to 0 (equivalent behavior for first hit)
    pub last_hit: i32,
    /// Flag indicating if hit was saved/extended (1 = extended, 0 = new hit)
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:666
    pub flag: u8,
}
