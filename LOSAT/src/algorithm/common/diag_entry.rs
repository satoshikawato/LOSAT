//! The diagonal table's entry of the protein word finders (NCBI's `DiagStruct`; TBLASTX, BLASTP,
//! BLASTX and the programs built on them) and the two-hit diagonal test of one hit on it
//! ([`x_test_hit`]), shared by BLASTP's scans and the one word of a hard chunk shorter than the
//! word (`common::hard_chunk_scan::OneWord`, BLASTP and BLASTX); the entry of a hit's diagonal
//! in a search's table ([`DiagTable`]: one ring, [`DiagRing`], or NCBI's ring of each query batch,
//! `common::query_batch_tables::QueryBatchTables`).
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
        diag.set_last_hit(window.wrapping_neg());
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-518, 776-779
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:46-65
// Each search of a query batch makes its own table, a ring of the power of two just at least the
// batch's concatenated length plus the window, and finds a hit's entry by the hit's diagonal in
// that batch's offsets: the query offset less the subject offset (two-hit), or the reverse
// (one-hit), masked by the ring's length.
/// The diagonal table of a word finder as a scan reads it: the entry of a hit's diagonal and the
/// first offset of the hit's query batch in the concatenated query, from which NCBI's search of
/// that batch counts its query offsets (the context test reads them,
/// `aa_ungapped.c:566-573`). One ring for the queries of one NCBI batch ([`DiagRing`], batch start
/// 0), or NCBI's ring of each batch where a search holds several
/// (`common::query_batch_tables::QueryBatchTables`).
pub(crate) trait DiagTable {
    /// The entry of the two-hit hit at query offset `q` and subject offset `s`, and the first
    /// offset of its query batch.
    fn two_hit_entry(&mut self, q: i32, s: i32) -> (&mut DiagStruct, i32);
    /// The entry of the one-hit hit at query offset `q` and subject offset `s`.
    fn one_hit_entry(&mut self, q: i32, s: i32) -> &mut DiagStruct;
}

/// One ring of diagonal entries (`mask`: its length less one) whose queries are one NCBI query
/// batch from offset 0: the table of a search of one batch (BLASTX's unit, BLASTP's call where the
/// scan shares it).
pub(crate) struct DiagRing<'a> {
    pub(crate) entries: &'a mut [DiagStruct],
    pub(crate) mask: i32,
}

impl DiagTable for DiagRing<'_> {
    #[inline(always)]
    fn two_hit_entry(&mut self, q: i32, s: i32) -> (&mut DiagStruct, i32) {
        (
            &mut self.entries[(q.wrapping_sub(s) & self.mask) as usize],
            0,
        )
    }

    #[inline(always)]
    fn one_hit_entry(&mut self, q: i32, s: i32) -> &mut DiagStruct {
        &mut self.entries[(s.wrapping_sub(q) & self.mask) as usize]
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-510, 566-573
// The scan's query offset is unsigned and counts from the start of the hit's query batch: the
// context test compares the offset less the distance with the context's first offset in unsigned
// arithmetic, so a distance past the batch's start passes it.
/// The distance `diff` of a two-hit hit at query offset `q` as the context test of a search of all
/// the queries reads it (`q - diff` in unsigned arithmetic below the context's first offset: a new
/// hit), so that the test gives the result of NCBI's search of the hit's query batch, which
/// starts at `batch_start`: where `q - diff` falls at or after 0 and before the batch's start,
/// NCBI's offset less the distance lies below 0 and passes the test, so the distance becomes one
/// past the query offset (the offset less it lies below 0 too). Every other distance is the same in
/// both (a context lies at or after the start of its batch).
#[inline(always)]
pub(crate) fn batch_context_distance(q: i32, diff: i32, batch_start: i32) -> i32 {
    if (q.wrapping_sub(diff) as u32) < batch_start as u32 {
        q.wrapping_add(1)
    } else {
        diff
    }
}
