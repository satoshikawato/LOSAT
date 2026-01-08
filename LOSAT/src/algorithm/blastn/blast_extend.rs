//! Diagonal table management for BLASTN
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.h:57-60
//!
//! This module contains the DiagStruct for tracking hits on diagonals
//! during the two-hit word finding algorithm.

/// Diagonal structure for tracking hits (NCBI DiagStruct equivalent)
///
/// Reference: blast_extend.h:57-60
/// ```c
/// typedef struct DiagStruct {
///     Int4 last_hit;  /**< last hit position on this diagonal */
///     Uint4 flag;     /**< flag indicating if hit was saved/extended */
/// } DiagStruct;
/// ```
#[derive(Clone, Copy, Default)]
pub struct DiagStruct {
    /// Last hit position on this diagonal (with diag_offset added)
    /// NCBI reference: blast_extend.c:103: diag_struct_array[i].last_hit = -diag->window;
    /// LOSAT: Initialize to 0 (equivalent behavior for first hit)
    pub last_hit: usize,
    /// Flag indicating if hit was saved/extended (1 = extended, 0 = new hit)
    /// NCBI reference: na_ungapped.c:666: hit_saved = hit_level_array[real_diag].flag;
    pub flag: u8,
}
