//! EXPERIMENT (LOSAT_X_MB24; in the strict set from the pull request that added it, 2026-10-26: on
//! unless it is `0`): the head table (`hashtable`) and the chains (`next_pos`) of a contiguous
//! megablast table of 4^12 cells (`build_mb_lookup`) with elements of 3 bytes instead of 4.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1081-1091
//! ```c
//!          if (mb_lt->hashtable[ecode] == 0) {
//!             PV_SET(pv_array, ecode, pv_array_bts);
//!          }
//!          else {
//!             helper_array[ecode/kCompressionFactor]++;
//!          }
//!          mb_lt->next_pos[index] = mb_lt->hashtable[ecode];
//!          mb_lt->hashtable[ecode] = index;
//! ```
//! NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1406-1418
//! ```c
//! Int4 q_off = lookup->hashtable[index];
//!
//! while (q_off) {
//!     offset_pairs[i].qs_offsets.q_off   = q_off - 1;
//!     offset_pairs[i++].qs_offsets.s_off = s_off;
//!     q_off = lookup->next_pos[q_off];
//! }
//! ```
//!
//! What it changes. NCBI's `hashtable` and `next_pos` are `Int4` arrays: one head per cell and one
//! link per query offset, each a 1-based query offset or 0 (the end of a chain). A table whose
//! `next_pos` has fewer than 2^24 elements holds only values below 2^24, which fit in 3 bytes. Here
//! both arrays are kept in the bytes of a `Vec<u32>`: element `i` is the low 24 bits of the
//! little-endian `u32` at byte `3 i`, and the vector has `words_for(n)` = `(3 n + 4) / 4` words for
//! `n` elements, at least `3 n + 1` bytes, so that the 4-byte read of the last element stays
//! inside it. A table of 4^12 cells takes 48 MiB of heads instead of 64, and the chains 3 bytes a
//! query position instead of 4 (megablast of Sakai, 10.9 M words: 31.5 MiB instead of 42). The head
//! table is the batch loop's kept table (`x_mb_lut_keep.rs`) read as bytes, as the 4-byte table is:
//! the build takes `words_for(4^12)` words of it.
//!
//! Why the table is NCBI's. Every element holds the value NCBI's `Int4` holds: the fills apply the
//! same update to the same words in the same order (each fill's own argument) and only the width of
//! the stores changes. The serial fills (`x_mb_fill.rs`, `x_fill_mb_queued`) read a head as 4 bytes
//! and write it back with its fourth byte (the first byte of the next cell) unchanged, and write a
//! link as one 4-byte store whose fourth byte, the first byte of the next element, is 0: the words
//! come in increasing order of their offsets (the queries' offsets ascend without overlap, `plan`,
//! and each query's words are read from its start), so the next element is not written yet, and an
//! element no word writes is never read (a chain only visits offsets of words). On a pool
//! (`x_mb_lookup_par.rs`, `x_mb_lookup_once.rs`) two workers may own neighbouring cells (at the
//! edge of a block) or neighbouring offsets, so there every head and link is read and written as
//! exactly its 3 bytes (`load_exact`, `store_exact`), and no worker touches another's bytes. The
//! readers (`MbLookupTable::for_each_hit`, `contains_hit`, `count_hits`) walk the same chains.
//! No floating-point operation is involved.
//!
//! Where (`plan`): a contiguous table of 4^12 cells (the megablast task's table for queries of about
//! 300,000 table entries or more), newest first (not the blastn task's re-threaded order), without
//! database word counts or `BLEMIR_DEBUG`, with fewer than 2^24 elements of `next_pos` and queries
//! whose offsets ascend without overlap (as the callers place them). Elsewhere the 4-byte table is
//! built as before; the fills of the smaller tables are not changed, nor are their readers' paths:
//! a 3-byte table leaves `hashtable` empty, so the readers find it behind the bounds check of the
//! head that the CSR tables of LOSAT_X_BNCSR already take.
//!
//! The variable is read once. `LOSAT_X_MB24` (a switch of the strict set, `utils/xswitch.rs`: on
//! unless it is `0`) runs this module.

use super::x_mb_lookup_par::Words;
use std::sync::OnceLock;

/// The cells of the tables this module takes (4^12).
pub(super) const CELLS: usize = 1 << 24;

/// The bits of an element.
pub(super) const MASK: u32 = (1 << 24) - 1;

/// No NCBI counterpart: reads LOSAT_X_MB24 once (a switch of the strict set: on unless it is `0`;
/// off builds the 4-byte table). It chooses a LOSAT layout and changes no value NCBI computes.
pub(super) fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_MB24"))
}

/// The `u32` words that hold `n` elements.
pub(super) const fn words_for(n: usize) -> usize {
    (3 * n + 4) / 4
}

/// The elements that `words` words hold (the largest `n` with `words_for(n) <= words`).
pub(super) const fn capacity(words: usize) -> usize {
    (4 * words).saturating_sub(1) / 3
}

/// Dispatch test of `build_mb_lookup`: the table of `cells` cells is built with 3-byte elements
/// (module documentation, "Where"); `next_pos_len`: the elements of its `next_pos`.
pub(super) fn plan(
    words: &Words<'_>,
    cells: usize,
    db_word_counts: bool,
    debug: bool,
    next_pos_len: usize,
) -> bool {
    if !on() {
        return false;
    }
    let taken = in_scope(words, cells, db_word_counts, debug, next_pos_len);
    crate::utils::xreach::x_reach!("LOSAT_X_MB24", taken);
    taken
}

fn in_scope(
    words: &Words<'_>,
    cells: usize,
    db_word_counts: bool,
    debug: bool,
    next_pos_len: usize,
) -> bool {
    if cells != CELLS
        || words.lut_word_length != 12
        || words.kmer_mask != (CELLS - 1) as u64
        || words.ascending_cells
        || db_word_counts
        || debug
        || next_pos_len > CELLS
    {
        return false;
    }
    ascending_offsets(words)
}

/// The queries the fills read (those of `lut_word_length` letters or more) have offsets that ascend
/// without overlap, and no mask ends before it starts (`build_unmasked_ranges` then gives ranges
/// that ascend without overlap; an inverted mask would give overlapping ranges), so their words come
/// in increasing order of offset.
fn ascending_offsets(words: &Words<'_>) -> bool {
    let mut end = 0usize;
    for (q, (seq, &offset)) in words
        .queries_blastna
        .iter()
        .zip(words.query_offsets)
        .enumerate()
    {
        if seq.len() < words.lut_word_length {
            continue;
        }
        let start = offset.max(0) as usize;
        let masks = words.query_masks.get(q).map_or(&[][..], Vec::as_slice);
        if start < end || masks.iter().any(|m| m.start > m.end) {
            return false;
        }
        end = start + seq.len();
    }
    true
}

/// A table's 3-byte heads and chains (`MbLookupTable::x_narrow`).
pub(super) struct Narrow {
    /// The heads, one per cell (the batch loop's kept table or a fresh one; it may be longer).
    heads: Vec<u32>,
    /// The links, one per element of NCBI's `next_pos`.
    chain: Vec<u32>,
}

impl Narrow {
    pub(super) fn new(heads: Vec<u32>, chain: Vec<u32>) -> Self {
        Self { heads, chain }
    }

    /// The head table, for the batch loop's keep (`MbLookupTable`'s `Drop`).
    pub(super) fn take_heads(&mut self) -> Vec<u32> {
        std::mem::take(&mut self.heads)
    }

    /// EXPERIMENT (LOSAT_X_BNSHRINKFREE): the heads (unless the keep took them) and the chain,
    /// shrunk before they are freed (`MbLookupTable`'s `Drop`, `x_shrink_free`).
    pub(super) fn x_shrink_free(&mut self) {
        crate::algorithm::blastn::x_shrink_free::free(std::mem::take(&mut self.heads));
        crate::algorithm::blastn::x_shrink_free::free(std::mem::take(&mut self.chain));
    }

    /// `for_each_hit` (blast_nascan.c:1406-1418 above): the chain of `cell`, newest first.
    #[inline(always)]
    pub(super) fn each(&self, cell: usize, mut callback: impl FnMut(u32)) {
        let mut q_off = get(&self.heads, cell);
        while q_off != 0 {
            callback(q_off);
            q_off = get(&self.chain, q_off as usize);
        }
    }

    /// `contains_hit`: whether the chain of `cell` holds `q_off_1` (the whole chain is walked,
    /// as the 4-byte table's test does).
    #[inline(always)]
    pub(super) fn contains(&self, cell: usize, q_off_1: u32) -> bool {
        let mut found = false;
        self.each(cell, |q_off| {
            if q_off == q_off_1 {
                found = true;
            }
        });
        found
    }

    /// `count_hits`: the length of the chain of `cell`.
    #[inline(always)]
    pub(super) fn count(&self, cell: usize) -> usize {
        let mut count = 0usize;
        self.each(cell, |_| count += 1);
        count
    }

    /// Tests only: the words of the heads.
    #[cfg(test)]
    pub(super) fn heads(&self) -> &[u32] {
        &self.heads
    }

    /// Tests only: the heads of `cells` cells and the `next_pos` of `next_pos_len` elements as
    /// 4-byte arrays.
    #[cfg(test)]
    pub(super) fn wide(&self, cells: usize, next_pos_len: usize) -> (Vec<u32>, Vec<u32>) {
        (widen(&self.heads, cells), widen(&self.chain, next_pos_len))
    }
}

/// Tests only (`x_mb_lookup_once`'s `OnceStats::first_pages`): the pages of `words` (at its address)
/// of the words that hold the first or the last byte of a cell other than 0 among the first `cells`:
/// the pages the read-once fill's marks give for the cells it wrote.
#[cfg(test)]
pub(super) fn written_pages(words: &[u32], cells: usize) -> std::collections::BTreeSet<usize> {
    const PAGE_BYTES: usize = 4096;
    let base = words.as_ptr() as usize;
    let mut pages = std::collections::BTreeSet::new();
    for (w, _) in words.iter().enumerate().filter(|(_, &v)| v != 0) {
        for cell in 4 * w / 3..=((4 * w + 3) / 3).min(cells.saturating_sub(1)) {
            if get(words, cell) != 0 {
                for byte in [3 * cell, 3 * cell + 2] {
                    pages.insert((base + 4 * (byte >> 2)) / PAGE_BYTES);
                }
            }
        }
    }
    pages
}

/// Tests only: the first `n` elements of `words` as 4-byte values.
#[cfg(test)]
pub(super) fn widen(words: &[u32], n: usize) -> Vec<u32> {
    (0..n).map(|i| get(words, i)).collect()
}

#[inline(always)]
fn bytes(words: &[u32]) -> &[u8] {
    // SAFETY: the bytes of the words (every byte is a valid `u8`), borrowed as long as `words`.
    unsafe { std::slice::from_raw_parts(words.as_ptr().cast::<u8>(), std::mem::size_of_val(words)) }
}

#[inline(always)]
fn bytes_mut(words: &mut [u32]) -> &mut [u8] {
    // SAFETY: as `bytes`, borrowed mutably as long as `words`.
    unsafe {
        std::slice::from_raw_parts_mut(
            words.as_mut_ptr().cast::<u8>(),
            std::mem::size_of_val(words),
        )
    }
}

/// Element `i` of `words`.
#[inline(always)]
pub(super) fn get(words: &[u32], i: usize) -> u32 {
    let b = &bytes(words)[3 * i..3 * i + 4];
    u32::from_le_bytes([b[0], b[1], b[2], b[3]]) & MASK
}

/// Element `i` of `words` set to `value` (below 2^24), the byte after it kept.
#[inline(always)]
pub(super) fn set(words: &mut [u32], i: usize, value: u32) {
    debug_assert!(value <= MASK);
    let b = &mut bytes_mut(words)[3 * i..3 * i + 3];
    b.copy_from_slice(&value.to_le_bytes()[..3]);
}

/// Element `i` of `words` set to `value` (below 2^24) by one 4-byte store that also writes 0 to the
/// first byte of element `i + 1` (the serial fills' links: that element is not written yet).
#[inline(always)]
pub(super) fn set_spill(words: &mut [u32], i: usize, value: u32) {
    debug_assert!(value <= MASK);
    bytes_mut(words)[3 * i..3 * i + 4].copy_from_slice(&value.to_le_bytes());
}

/// The address of element `i` of the words at `base` (for a prefetch; nothing is read).
#[inline(always)]
pub(super) fn at(base: *const u32, i: usize) -> *const i8 {
    base.cast::<i8>().wrapping_add(3 * i)
}

/// The 4 bytes at element `i` of the words at `base` (the element in the low 24 bits).
/// SAFETY: the 4 bytes lie in the words, and no other thread writes them.
#[inline(always)]
pub(super) unsafe fn load_word(base: *const u32, i: usize) -> u32 {
    // SAFETY: the caller's.
    u32::from_le(unsafe { base.cast::<u8>().add(3 * i).cast::<u32>().read_unaligned() })
}

/// Writes `word` as the 4 bytes at element `i` of the words at `base`.
/// SAFETY: the 4 bytes lie in the words, and no other thread reads or writes them.
#[inline(always)]
pub(super) unsafe fn store_word(base: *mut u32, i: usize, word: u32) {
    // SAFETY: the caller's.
    unsafe {
        base.cast::<u8>()
            .add(3 * i)
            .cast::<u32>()
            .write_unaligned(word.to_le())
    }
}

/// Element `i` of the words at `base`, reading its 3 bytes only.
/// SAFETY: the 3 bytes lie in the words, and no other thread writes them.
#[inline(always)]
pub(super) unsafe fn load_exact(base: *const u32, i: usize) -> u32 {
    // SAFETY: the caller's.
    unsafe {
        let p = base.cast::<u8>().add(3 * i);
        u32::from(u16::from_le(p.cast::<u16>().read_unaligned())) | u32::from(*p.add(2)) << 16
    }
}

/// Element `i` of the words at `base` set to `value` (below 2^24), writing its 3 bytes only.
/// SAFETY: the 3 bytes lie in the words, and no other thread reads or writes them.
#[inline(always)]
pub(super) unsafe fn store_exact(base: *mut u32, i: usize, value: u32) {
    debug_assert!(value <= MASK);
    // SAFETY: the caller's.
    unsafe {
        let p = base.cast::<u8>().add(3 * i);
        p.cast::<u16>().write_unaligned((value as u16).to_le());
        *p.add(2) = (value >> 16) as u8;
    }
}

#[cfg(test)]
mod tests {
    use super::super::MaskedInterval;
    use super::*;

    fn words<'a>(
        refs: &'a [&'a [u8]],
        offsets: &'a [i32],
        masks: &'a [Vec<MaskedInterval>],
        lut: usize,
    ) -> Words<'a> {
        Words {
            queries_blastna: refs,
            query_offsets: offsets,
            query_masks: masks,
            word_length: 16,
            lut_word_length: lut,
            kmer_mask: (1u64 << (2 * lut)) - 1,
            ascending_cells: false,
        }
    }

    // The scope: 4^12 cells and `lut_word_length` 12, newest first, no database word counts or
    // BLEMIR_DEBUG, `next_pos` of at most 2^24 elements, query offsets that ascend without overlap
    // (queries shorter than a word do not count) and no inverted mask.
    #[test]
    fn the_scope_takes_ascending_queries_of_4_12_tables_only() {
        let queries: Vec<Vec<u8>> = vec![vec![0; 100], vec![1; 5], vec![2; 50]];
        let refs: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
        let plain = vec![Vec::new(); 3];
        let ok = words(&refs, &[0, 101, 107], &plain, 12);
        assert!(in_scope(&ok, CELLS, false, false, 158));
        assert!(in_scope(&ok, CELLS, false, false, CELLS));
        assert!(!in_scope(&ok, CELLS, false, false, CELLS + 1));
        assert!(!in_scope(&ok, CELLS, true, false, 158));
        assert!(!in_scope(&ok, CELLS, false, true, 158));
        let lut11 = words(&refs, &[0, 101, 107], &plain, 11);
        assert!(!in_scope(&lut11, 1 << 22, false, false, 158));
        let ascending = Words {
            ascending_cells: true,
            ..words(&refs, &[0, 101, 107], &plain, 12)
        };
        assert!(!in_scope(&ascending, CELLS, false, false, 158));
        // a query shorter than a word does not count; overlapping or descending ones do
        let short_first = words(&refs, &[0, 0, 101], &plain, 12);
        assert!(in_scope(&short_first, CELLS, false, false, 152));
        let overlapping = words(&refs, &[0, 101, 99], &plain, 12);
        assert!(!in_scope(&overlapping, CELLS, false, false, 150));
        let descending = words(&refs, &[60, 0, 0], &plain, 12);
        assert!(!in_scope(&descending, CELLS, false, false, 161));
        let inverted = vec![
            vec![MaskedInterval { start: 70, end: 20 }],
            Vec::new(),
            Vec::new(),
        ];
        let inverted = words(&refs, &[0, 101, 107], &inverted, 12);
        assert!(!in_scope(&inverted, CELLS, false, false, 158));
    }

    #[test]
    fn the_words_hold_the_elements_and_the_last_read() {
        for n in 0..64 {
            let w = words_for(n);
            assert!(4 * w >= 3 * n + 1, "{n}");
            assert!(capacity(w) >= n, "{n}");
            assert!(words_for(capacity(w)) <= w, "{n}");
        }
        assert_eq!(words_for(CELLS), 12_582_913);
    }

    #[test]
    fn every_store_writes_its_element_only() {
        let n = 9;
        let values = [
            0,
            1,
            MASK,
            0x00AB_CDEF,
            0x0001_0000,
            0x0000_FFFF,
            7,
            MASK - 1,
            0x0080_0000,
        ];
        // set (and the exact store) keep the neighbours; set_spill writes in increasing order
        let mut a = vec![0u32; words_for(n)];
        let mut b = vec![0u32; words_for(n)];
        let mut c = vec![0u32; words_for(n)];
        for (i, &v) in values.iter().enumerate().rev() {
            set(&mut a, i, v);
            // SAFETY: element `i` lies in `b`, which no other thread uses.
            unsafe { store_exact(b.as_mut_ptr(), i, v) };
        }
        for (i, &v) in values.iter().enumerate() {
            set_spill(&mut c, i, v);
        }
        for (i, &v) in values.iter().enumerate() {
            assert_eq!(get(&a, i), v);
            assert_eq!(get(&b, i), v);
            assert_eq!(get(&c, i), v);
            // SAFETY: as above.
            unsafe {
                assert_eq!(load_exact(b.as_ptr(), i), v);
                assert_eq!(load_word(b.as_ptr(), i) & MASK, v);
            }
        }
        assert_eq!(a, b);
        assert_eq!(a, c);
        // a merged store keeps the next element's byte
        // SAFETY: as above.
        unsafe {
            let w = load_word(a.as_ptr(), 3);
            store_word(a.as_mut_ptr(), 3, (w & !MASK) | 5);
        }
        assert_eq!(get(&a, 3), 5);
        assert_eq!(get(&a, 4), values[4]);
        assert_eq!(widen(&a, n)[4..], values[4..]);
    }
}
