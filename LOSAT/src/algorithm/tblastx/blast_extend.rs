//! Diagonal table management for TBLASTX
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:52-180
//!
//! This module contains the DiagStruct for tracking hits on diagonals
//! during the two-hit word finding algorithm.

/// Diagonal structure for tracking hits (NCBI DiagStruct equivalent)
///
/// NCBI uses a bitfield: signed last_hit:31 + unsigned flag:1.
/// Reference: ncbi-blast/c++/include/algo/blast/core/blast_extend.h:57-61
#[derive(Clone, Copy)]
pub struct DiagStruct {
    /// Packed bitfield storage mirroring NCBI layout.
    /// Reference: ncbi-blast/c++/include/algo/blast/core/blast_extend.h:57-61
    raw: u32,
}

impl Default for DiagStruct {
    /// NCBI: calloc zeros memory for initial allocation
    /// Reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:109-180
    fn default() -> Self {
        Self { raw: 0 }
    }
}

impl DiagStruct {
    // Bitfield masks (last_hit: bits 0-30, flag: bit 31).
    // Reference: ncbi-blast/c++/include/algo/blast/core/blast_extend.h:57-61
    const LAST_HIT_MASK: u32 = 0x7fff_ffff;
    const FLAG_MASK: u32 = 0x8000_0000;

    /// Return the flag (0/1).
    #[inline]
    pub fn flag(&self) -> u32 {
        (self.raw >> 31) & 1
    }

    /// Set the flag (0/1).
    #[inline]
    pub fn set_flag(&mut self, flag: u32) {
        let bit = (flag & 1) << 31;
        self.raw = (self.raw & Self::LAST_HIT_MASK) | bit;
    }

    /// Return last_hit as signed 31-bit value (sign-extended).
    #[inline]
    pub fn last_hit(&self) -> i32 {
        let raw = self.raw & Self::LAST_HIT_MASK;
        // Sign-extend 31-bit to 32-bit signed.
        // Reference: ncbi-blast/c++/include/algo/blast/core/blast_extend.h:57-61
        ((raw << 1) as i32) >> 1
    }

    /// Set last_hit, truncating to 31-bit signed.
    #[inline]
    pub fn set_last_hit(&mut self, value: i32) {
        let masked = (value as u32) & Self::LAST_HIT_MASK;
        self.raw = (self.raw & Self::FLAG_MASK) | masked;
    }

    /// Clear diagonal state for a new scan
    ///
    /// NCBI s_BlastDiagClear
    /// Reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:101-103
    /// ```c
    /// diag_struct_array[i].flag = 0;
    /// diag_struct_array[i].last_hit = -diag->window;
    /// ```
    #[inline]
    pub fn clear(window: i32) -> Self {
        let mut diag = Self::default();
        diag.set_flag(0);
        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:101-103
        diag.set_last_hit(-window);
        diag
    }
}
