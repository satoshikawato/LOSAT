//! The diagonal table's entry of the protein word finders (NCBI's `DiagStruct`; TBLASTX, BLASTP,
//! BLASTX and the programs built on them) and the two-hit diagonal test of one hit on it
//! ([`x_test_hit`]), shared by BLASTP's scans and the one word of a hard chunk shorter than the
//! word (`common::hard_chunk_scan::OneWord`, BLASTP and BLASTX).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-180

/// Diagonal structure for tracking hits (NCBI DiagStruct equivalent)
///
/// NCBI uses a bitfield: signed last_hit:31 + unsigned flag:1.
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:57-61
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct DiagStruct {
    /// Packed bitfield storage mirroring NCBI layout.
    /// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:57-61
    raw: u32,
}

impl Default for DiagStruct {
    /// NCBI: calloc zeros memory for initial allocation
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:109-180
    fn default() -> Self {
        Self { raw: 0 }
    }
}

impl DiagStruct {
    // Bitfield masks (last_hit: bits 0-30, flag: bit 31).
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:57-61
    const LAST_HIT_MASK: u32 = 0x7fff_ffff;
    const FLAG_MASK: u32 = 0x8000_0000;

    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:57-60
    // No computation: it returns the packed (last_hit, flag) word of one cell so that the hit loops
    // and the shadow modes of other switches can read and compare diagonal tables cell by cell. It
    // does not change any value NCBI computes.
    /// The packed cell, for comparisons.
    #[inline]
    pub fn raw_bits(&self) -> u32 {
        self.raw
    }

    // No NCBI counterpart: it stores a packed (last_hit, flag) word that the caller computed with
    // the bit-field's rules (last_hit truncated to 31 bits, flag in bit 31), so that
    // LOSAT_X_BPDIAG's batched two-hit test (`blastp/x_seed_fused.rs`) writes an entry in one
    // store. It does not change any value NCBI computes.
    /// Store the packed cell.
    #[inline]
    pub fn set_raw_bits(&mut self, raw: u32) {
        self.raw = raw;
    }

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
        // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:57-61
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
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:101-103
    #[inline]
    pub fn clear(window: i32) -> Self {
        let mut diag = Self::default();
        diag.set_flag(0);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:101-103
        diag.set_last_hit(-window);
        diag
    }
}

/// The diagonal test of one word hit as NCBI applies it before an extension, on the hit's own
/// entry. Returns `Some((last_hit, diff))` when the hit needs the rest of the test.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552
//
// The sums and differences wrap (two's complement i32): that is what `+` and `-` do in release
// builds (no overflow checks), and what NCBI's arithmetic on its `subject_offset` (a `Uint4`) gives.
// LOSAT_X_BPDIAG's batches compute the same wrapping values per lane, and their tests compare
// with this function on offsets whose sum passes 2^31 - 1.
#[inline(always)]
pub(crate) fn x_test_hit(
    entry: &mut DiagStruct,
    subject_offset: i32,
    diag_offset: i32,
    window: i32,
    word_size: i32,
) -> Option<(i32, i32)> {
    if entry.flag() != 0 {
        let subject_plus_offset = subject_offset.wrapping_add(diag_offset);
        if subject_plus_offset >= entry.last_hit() {
            entry.set_last_hit(subject_plus_offset);
            entry.set_flag(0);
        }
        return None;
    }
    let last_hit = entry.last_hit().wrapping_sub(diag_offset);
    let diff = subject_offset.wrapping_sub(last_hit);
    if diff >= window {
        entry.set_last_hit(subject_offset.wrapping_add(diag_offset));
        return None;
    }
    if diff < word_size {
        return None;
    }
    Some((last_hit, diff))
}
