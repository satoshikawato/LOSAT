//! The two-hit seed stage of blastp-fast (the compressed lookup table,
//! `CompressedAaLookupTable`) with the subject scan and the diagonal test in one loop.
//!
//! The reference path copies every word hit of the subject into the offset-pair buffer
//! (`BlastCompressedAaLookupTable::scan_subject`, NCBI's `s_BlastCompressedAaScanSubject`) and then
//! reads the buffer back for the two-hit test (NCBI's `s_BlastAaWordFinder_TwoHit`). Here the scan
//! lists the subject offsets whose word passes `PV_TEST` (with the word's index), then reads those
//! lookup cells in that order (the backbone cell a few cells ahead is prefetched) and hands each
//! cell's hits, in place of the buffer, to the test of the amino-acid lookup table's fused loop
//! (`x_seed_fused::x_test_hit`, then the caller's closure for a hit that needs the context check and
//! the extension).
//!
//! Why the outcome is the reference's: the argument of `x_seed_fused.rs` holds once the hits reach
//! the test in the reference's order, because the test of a hit reads and writes only the entry of
//! its diagonal and the closure of a candidate writes only that entry and the HSP list. The scan
//! reads only the subject and the lookup table, which the test does not write, so when the scan's
//! steps run relative to the tests (all of a round first in the reference, the list of offsets
//! first here) changes nothing. The hits come in the reference's order:
//! - the subject offsets are scanned in increasing order over one range, `0..=max(length - word, 0)`,
//!   as `blastp_for_each_offset_pair` sets `scan_range` (the subject has one range, `seq_ranges[0]`,
//!   and no masked ranges in this scan).
//! - the index of each offset is formed by the reference's statements: the same priming
//!   (`compute_compressed_index_u8` over `word - 1` letters, letters past the buffer read as
//!   NULLB), the same rolling update (`compressed_reciprocal_divide`), the same flush of a letter
//!   outside the compressed alphabet (`scaled_compress_table` entry -1, the next `word - 1` letters
//!   only refill the index), and the same `PV_TEST` and `num_used` checks (`num_used` read when the
//!   listed offsets' cells are read). A word is tested only when its `word` letters are all inside
//!   the compressed alphabet.
//! - a cell's hits are taken in the order NCBI's scan copies them: `query_offset`, then the
//!   backbone's other offsets (up to 5 hits), or for more than 5 hits the two offsets kept in the
//!   backbone, then the overflow list from its head: the cell filled last first (it holds
//!   `(n - 3) & 3` offsets, or 4), then the older cells of 4, newest first (`x_cell_hits`, the walk
//!   of `copy_scan_hits`).
//! - the offset-pair buffer's capacity (`array_size` = `OFFSET_ARRAY_SIZE` + the longest cell) only
//!   cuts the reference's scan into rounds. A round ends before a cell whose hits do not fit, and
//!   the next round starts at that offset with an empty buffer, which takes any cell, so a cell is
//!   never split and the scan always moves on. The next round primes the index again at that offset:
//!   the word there was tested, so its first `word - 1` letters are inside the alphabet and the
//!   priming stops at once, with the index `sum c[s + k] * A^(k + 1)` (k < word - 1, `A` the
//!   compressed alphabet's size). The rolling scan has `index / A` there as `preshift`, and the
//!   priming's `preshift` is `(index * recip) >> 32` of the primed index; both equal the exact
//!   quotient because NCBI's reciprocal (`ceil(2^32 / A)`) gives `floor(i / A)` for every
//!   `i < A^word`, the bound of every index the scan forms (checked for all three tables in
//!   `the_reciprocal_gives_the_exact_quotient_below_every_tables_size`). From there the rounds run
//!   the statements of the one scan here on the same state (offset, index, `preshift`,
//!   `next_char`), so the hits and their order are those of a scan without rounds.
//!
//! The closure and the prefetch distance are those of `x_seed_fused.rs`.

use std::cell::Cell;

use super::x_seed_fused::{x_prefetch, x_test_hit, X_PREFETCH_DISTANCE};
use crate::algorithm::tblastx::blast_extend::DiagStruct;
#[cfg(target_arch = "x86_64")]
use crate::algorithm::tblastx::lookup::compressed::CompressedLookupBackboneCell;
use crate::algorithm::tblastx::lookup::compressed::{
    compressed_pv_test, compressed_reciprocal_divide, compute_compressed_index_u8,
    BlastCompressedAaLookupTable, CompressedBackbonePayload, COMPRESSED_HITS_CELL_MASK,
    COMPRESSED_HITS_PER_BACKBONE_CELL, COMPRESSED_HITS_PER_OVERFLOW_CELL,
};

/// No NCBI counterpart: the buffers of the fused loop, reused across subjects (NCBI copies the hits
/// into its offset-pair array).
#[derive(Default)]
struct XScratch {
    /// The scan's cells, in scan order: (subject offset, backbone index).
    cells: Vec<(i32, u32)>,
    /// The hits of one cell (`x_cell_hits`).
    hits: Vec<i32>,
}

thread_local! {
    static X_SCRATCH: Cell<Option<XScratch>> = const { Cell::new(None) };
}

// No NCBI counterpart: how many cells ahead the fused loop reads the backbone cells. Tuning only.
const X_CELL_AHEAD: usize = 8;

// No NCBI counterpart: a cache hint for one backbone cell; it reads nothing and changes no value.
#[inline(always)]
fn x_prefetch_cell(lookup: &BlastCompressedAaLookupTable, index: usize) {
    #[cfg(target_arch = "x86_64")]
    if let Some(cell) = lookup.backbone.get(index) {
        // SAFETY: a prefetch of a valid address has no effect on program state.
        unsafe {
            std::arch::x86_64::_mm_prefetch(
                (cell as *const CompressedLookupBackboneCell).cast::<i8>(),
                std::arch::x86_64::_MM_HINT_T0,
            );
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    let _ = (lookup, index);
}

/// Writes the hits of backbone cell `index` (with `numhits` > 0 hits) into `dest` in NCBI's scan
/// order and returns `numhits`; `dest` holds at least `numhits + 8` entries.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:348-396
// A list of NCBI's shape holds `numhits - 3` offsets: the head `first_cell_entries` (or 4), every
// other cell 4. The copies here are of whole cells (one move of four offsets); the offsets copied
// past the cell's hits are overwritten by the next cell or never read.
#[inline(always)]
fn x_cell_hits(
    lookup: &BlastCompressedAaLookupTable,
    index: usize,
    numhits: usize,
    dest: &mut [i32],
) -> usize {
    let cell = &lookup.backbone[index];
    dest[0] = cell.query_offset;
    if numhits <= COMPRESSED_HITS_PER_BACKBONE_CELL + 1 {
        let CompressedBackbonePayload::Inline(query_offsets) = &cell.payload else {
            unreachable!("NCBI BLAST small compressed cell uses inline storage");
        };
        dest[1..1 + COMPRESSED_HITS_PER_BACKBONE_CELL].copy_from_slice(query_offsets);
        return numhits;
    }
    let CompressedBackbonePayload::Overflow {
        query_offsets,
        head,
    } = &cell.payload
    else {
        unreachable!("NCBI BLAST large compressed cell uses overflow storage");
    };
    dest[1] = query_offsets[0];
    dest[2] = query_offsets[1];
    let overflow = &lookup.overflow_cells;
    let first_cell_entries = (numhits - 3) & COMPRESSED_HITS_CELL_MASK;
    let mut curr_cell = head.expect("NCBI BLAST compressed overflow head must exist");
    let mut written = 3;
    if first_cell_entries != 0 {
        dest[written..written + COMPRESSED_HITS_PER_OVERFLOW_CELL]
            .copy_from_slice(&overflow[curr_cell].query_offsets);
        written += first_cell_entries;
        let Some(next) = overflow[curr_cell].next else {
            debug_assert_eq!(written, numhits);
            return written;
        };
        curr_cell = next;
    }
    while written < numhits {
        dest[written..written + COMPRESSED_HITS_PER_OVERFLOW_CELL]
            .copy_from_slice(&overflow[curr_cell].query_offsets);
        written += COMPRESSED_HITS_PER_OVERFLOW_CELL;
        let Some(next) = overflow[curr_cell].next else {
            break;
        };
        curr_cell = next;
    }
    debug_assert_eq!(written, numhits);
    written.min(numhits)
}

/// Scans `subject` with the compressed `lookup` as the reference's rounds of
/// `scan_subject` do together, and calls `visit(subject_offset, hits)` for every cell found, with the
/// cell's hits in NCBI's order.
///
/// `subject` starts at the first residue and keeps the sentinel after the last one, as the
/// reference scan reads it; `subject_length` is the number of residues.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:495-503
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:274-341
// The rounds that `array_size` makes are not repeated (module header). The scan first lists the
// offsets whose index passes `PV_TEST` (with the index), then reads their cells in that order, a few
// cells ahead of the test (the cells are read, and their hits tested, in the reference's order).
#[inline(always)]
fn x_scan_cells<V>(
    lookup: &BlastCompressedAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    scratch: &mut XScratch,
    mut visit: V,
) where
    V: FnMut(i32, &[i32]),
{
    let word_length = usize::try_from(lookup.word_length)
        .expect("NCBI BLAST compressed word length must fit usize");
    let compressed_alphabet_size = usize::try_from(lookup.compressed_alphabet_size)
        .expect("NCBI BLAST compressed alphabet size must fit usize");
    let recip = i64::from(lookup.reciprocal_alphabet_size);
    let pv_array_bts = usize::try_from(lookup.pv_array_bts)
        .expect("NCBI BLAST compressed PV shift must fit usize");
    let scaled_compress_table = &lookup.scaled_compress_table_u8;
    let pv = &lookup.pv;
    let backbone = &lookup.backbone;
    let longest_chain = usize::try_from(lookup.longest_chain.max(0))
        .expect("NCBI BLAST compressed longest chain must fit usize");
    let XScratch { cells, hits } = scratch;
    if hits.len() < longest_chain + 8 {
        hits.resize(longest_chain + 8, 0);
    }
    cells.clear();

    let mut s_last =
        i32::try_from(subject_length).expect("subject length must fit Int4") - lookup.word_length;
    if s_last < 0 {
        s_last = 0;
    }
    let s_last = s_last as usize;
    // The reference reads a letter past the subject's buffer (a subject shorter than a word) as
    // NULLB, whose entry is -1.
    let byte_at = |i: usize| subject.get(i).copied().unwrap_or(0);

    let mut s = 0usize;
    let mut index = 0i32;
    let mut padded_word = [0u8; 8];
    while s <= s_last {
        let end = s + word_length - 1;
        let word = match subject.get(s..end) {
            Some(word) => word,
            None => {
                let available = subject.get(s..).unwrap_or(&[]);
                let word = &mut padded_word[..word_length - 1];
                word.fill(0);
                word[..available.len()].copy_from_slice(available);
                &padded_word[..word_length - 1]
            }
        };
        let (computed, skip) = compute_compressed_index_u8(
            word_length - 1,
            word,
            compressed_alphabet_size,
            scaled_compress_table,
        );
        index = computed;
        if skip == 0 {
            break;
        }
        s += 1;
    }

    let mut next_char = if s <= s_last {
        byte_at(s + word_length - 1)
    } else {
        0
    };
    let mut preshift = compressed_reciprocal_divide(index, recip);

    while s <= s_last {
        let mut compressed_char = scaled_compress_table[usize::from(next_char)];
        next_char = byte_at(s + word_length);

        if compressed_char < 0 {
            preshift = 0;
            s += 1;
            let mut skip = lookup.word_length - 1;
            while skip != 0 && s <= s_last {
                compressed_char = scaled_compress_table[usize::from(next_char)];
                next_char = byte_at(s + word_length);
                if compressed_char < 0 {
                    skip = lookup.word_length - 1;
                    preshift = 0;
                    s += 1;
                    continue;
                }
                index = preshift + compressed_char;
                preshift = compressed_reciprocal_divide(index, recip);
                skip -= 1;
                s += 1;
            }
            // NCBI's `s--` is balanced by the outer `for` increment; `s` already holds the
            // offset the scan resumes at.
            continue;
        }

        index = preshift + compressed_char;
        preshift = compressed_reciprocal_divide(index, recip);

        // index = preshift + compressed_char >= 0: both terms are (from) non-negative letters.
        if compressed_pv_test(pv, index as usize, pv_array_bts) {
            cells.push((s as i32, index as u32));
        }
        s += 1;
    }

    for (n, &(s, index)) in cells.iter().enumerate() {
        if let Some(&(_, later)) = cells.get(n + X_CELL_AHEAD) {
            x_prefetch_cell(lookup, later as usize);
        }
        let index = index as usize;
        let numhits = backbone[index].num_used;
        if numhits != 0 {
            let written = x_cell_hits(lookup, index, numhits as usize, hits);
            visit(s, &hits[..written]);
        }
    }
}

/// Scans `subject` with the compressed `lookup` (two-hit mode) and applies the diagonal test to every
/// word hit in the reference's order; `candidate(query_offset, subject_offset, last_hit, diff,
/// entry)` receives each hit that needs the rest of the test, with its diagonal entry. The
/// arguments are those of `x_seed_fused::x_scan_two_hit`.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:508-518
#[allow(clippy::too_many_arguments)]
pub(crate) fn x_scan_two_hit_compressed<F>(
    lookup: &BlastCompressedAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    diag_array: &mut [DiagStruct],
    diag_mask: i32,
    diag_offset: i32,
    window: i32,
    word_size: i32,
    mut candidate: F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    assert!(diag_mask >= 0 && diag_array.len() > diag_mask as usize);
    let prefetch_distance = X_PREFETCH_DISTANCE;
    let mut scratch = X_SCRATCH.with(Cell::take).unwrap_or_default();
    x_scan_cells(
        lookup,
        subject,
        subject_length,
        &mut scratch,
        |s, cell_hits| {
            for (i, &query_offset) in cell_hits.iter().enumerate() {
                // Read ahead the entry of a later hit of this cell (a cache hint; see `x_scan_two_hit`).
                if let Some(&ahead) = cell_hits.get(i + prefetch_distance) {
                    x_prefetch(diag_array, ((ahead - s) & diag_mask) as usize);
                }
                let diag_coord = ((query_offset - s) & diag_mask) as usize;
                // SAFETY: diag_coord <= diag_mask < diag_array.len() (asserted above).
                let entry = unsafe { diag_array.get_unchecked_mut(diag_coord) };
                if let Some((last_hit, diff)) = x_test_hit(entry, s, diag_offset, window, word_size)
                {
                    candidate(query_offset, s, last_hit, diff, entry);
                }
            }
        },
    );
    X_SCRATCH.with(|cell| cell.set(Some(scratch)));
}

#[cfg(test)]
#[path = "x_seed_fused_compressed_tests.rs"]
mod tests;
