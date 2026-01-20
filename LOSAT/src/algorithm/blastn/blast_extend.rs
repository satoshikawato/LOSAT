//! Diagonal table management for BLASTN
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.h:57-60
//!
//! This module contains the DiagStruct for tracking hits on diagonals
//! during the two-hit word finding algorithm.

/// Diagonal structure for tracking hits (NCBI DiagStruct equivalent)
///
/// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_extend.h:57-60
/// ```c
/// typedef struct DiagStruct {
///    signed int last_hit   : 31; /**< Offset of the last hit */
///    unsigned int flag      : 1 ; /**< Reset the next extension? */
/// } DiagStruct;
/// ```
#[derive(Clone, Copy, Default)]
pub struct DiagStruct {
    /// Last hit position on this diagonal (with diag_offset added)
    /// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:96-105
    /// ```c
    /// for (i = 0; i < n; i++) {
    ///     diag_struct_array[i].flag = 0;
    ///     diag_struct_array[i].last_hit = -diag->window;
    ///     if (diag->hit_len_array) diag->hit_len_array[i] = 0;
    /// }
    /// ```
    /// LOSAT: Initialize to 0 (equivalent behavior for first hit)
    pub last_hit: i32,
    /// Flag indicating if hit was saved/extended (1 = extended, 0 = new hit)
    /// NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:666
    /// ```c
    /// hit_saved = hit_level_array[real_diag].flag;
    /// ```
    pub flag: u8,
}
