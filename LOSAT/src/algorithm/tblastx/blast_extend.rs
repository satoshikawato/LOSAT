//! Diagonal table management for TBLASTX
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:52-180
//!
//! This module contains the DiagStruct for tracking hits on diagonals
//! during the two-hit word finding algorithm.

/// Diagonal structure for tracking hits (NCBI DiagStruct equivalent)
///
/// Reference: blast_extend.c:52-61
/// ```c
/// typedef struct DiagStruct {
///     Int4 last_hit;  /**< Offset of the last hit on this diagonal */
///     Uint1 flag;     /**< Whether this diagonal has been extended */
/// } DiagStruct;
/// ```
#[derive(Clone, Copy)]
pub struct DiagStruct {
    /// Last hit position on this diagonal
    /// NCBI reference: blast_extend.c:53
    pub last_hit: i32,
    /// Flag indicating if hit was saved/extended (1 = extended, 0 = new hit)
    /// NCBI reference: blast_extend.c:54
    pub flag: u8,
}

impl Default for DiagStruct {
    /// NCBI: calloc zeros memory for initial allocation
    /// Reference: blast_extend.c:109-180 BlastExtendWordNew
    fn default() -> Self {
        Self { last_hit: 0, flag: 0 }
    }
}

impl DiagStruct {
    /// Clear diagonal state for a new scan
    ///
    /// NCBI s_BlastDiagClear (blast_extend.c:101-103):
    /// ```c
    /// diag_struct_array[i].flag = 0;
    /// diag_struct_array[i].last_hit = -diag->window;
    /// ```
    #[inline]
    pub fn clear(window: i32) -> Self {
        Self { last_hit: -window, flag: 0 }
    }
}
