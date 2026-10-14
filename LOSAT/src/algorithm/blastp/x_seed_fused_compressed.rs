//! EXPERIMENT (LOSAT_X_BPSEED): the two-hit seed stage of blastp-fast (the compressed lookup table,
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
//! The closure, the switches (`LOSAT_X_BPSEED`, `LOSAT_X_BPSEEDSHADOW`, `LOSAT_X_BPSEED_PREFETCH`) and
//! the shadow are those of `x_seed_fused.rs`. `LOSAT_X_BPSEED` is a switch of the
//! strict set (`utils/xswitch.rs`): on unless it is `0`, or `LOSAT_X_DEFAULT=0` and it is not set.

use std::cell::Cell;

use super::x_seed_fused::{x_prefetch, x_prefetch_distance, x_test_hit};
use crate::algorithm::tblastx::blast_extend::DiagStruct;
use crate::algorithm::tblastx::lookup::compressed::{
    compressed_pv_test, compressed_reciprocal_divide, compute_compressed_index_u8,
    BlastCompressedAaLookupTable, CompressedBackbonePayload, CompressedLookupBackboneCell,
    COMPRESSED_HITS_CELL_MASK, COMPRESSED_HITS_PER_BACKBONE_CELL,
    COMPRESSED_HITS_PER_OVERFLOW_CELL,
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
// ```c
//              dest[0].qs_offsets.q_off = backbone_cell->query_offset;
//              dest[0].qs_offsets.s_off = s_off;
//
//              if (numhits <= COMPRESSED_HITS_PER_BACKBONE_CELL+1) {
//             	 dest ++;
//                 /* hits all live in the backbone */
//                 for (i = 0; i < numhits-1; i++) {
//                    dest[i].qs_offsets.q_off = backbone_cell->payload.query_offsets[i];
//                    dest[i].qs_offsets.s_off = s_off;
//                 }
//              }
//              else {
//                 /* hits are in the backbone cell and in the overflow list */
//                 CompressedOverflowCell* curr_cell =
//                                     backbone_cell->payload.overflow_list.head;
//                 ...
//                 Int4 first_cell_entries = (numhits - 3) & COMPRESSED_HITS_CELL_MASK;
//
//                 /* copy hits from backbone */
//                 dest[1].qs_offsets.q_off = backbone_cell->payload.overflow_list.query_offsets[0];
//                 dest[1].qs_offsets.s_off = s_off;
//                 dest[2].qs_offsets.q_off = backbone_cell->payload.overflow_list.query_offsets[1];
//                 dest[2].qs_offsets.s_off = s_off;
//
//                 /* handle the overflow list */
//
//                 /* first cell can be partially filled */
//                 dest += 3;
//                 for (i = 0; i < first_cell_entries; i++) {
//                    dest[i].qs_offsets.q_off =  curr_cell->query_offsets[i];
//                    dest[i].qs_offsets.s_off = s_off;
//                 }
//
//                 /* handle the rest of the list */
//                 if (first_cell_entries) {
//                    curr_cell = curr_cell->next;
//                 }
//                    while (curr_cell != NULL) {
//                       dest += i;
//                       for (i = 0; i < COMPRESSED_HITS_PER_OVERFLOW_CELL; i++) {
//                          dest[i].qs_offsets.q_off =  curr_cell->query_offsets[i];
//                          dest[i].qs_offsets.s_off = s_off;
//                       }
//                       curr_cell = curr_cell->next;    /* prefetch */
//                    }
//              }
// ```
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
// ```c
//     scan_range[0] = 0;
//     scan_range[1] = subject->seq_ranges[0].left;
//     scan_range[2] = subject->seq_ranges[0].right - wordsize;
//
//     if (scan_range[2] < scan_range[1])
//         scan_range[2] = scan_range[1];
//
//     while (scan_range[1] <= scan_range[2]) {
//         /* scan the subject sequence for hits */
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:274-341
// ```c
//     /* prime the index */
//     for(s = s_first; s <= s_last; s++){
//         index = s_ComputeCompressedIndex(word_length - 1, s,
//                                          compressed_alphabet_size,
//                                          &skip, lookup);
//         if(!skip)
//           break;
//     }
//
//     next_char = ((s <= s_last)? s[word_length-1] : 0);
//     preshift = (Int4)((((Int8)index) * recip) >> 32);
//
//     /* main scanning loop */
//     for (; s <= s_last; s++) {
//        /* compute the index value */
//
//        compressed_char = scaled_compress_table[next_char];
//        next_char = s[word_length];
//
//        if(compressed_char < 0){ /* flush (rare) "bad" character(s) */
//          preshift = 0;
//          s++;
//          for(skip = word_length-1; skip && (s <= s_last) ; s++){
//            compressed_char = scaled_compress_table[next_char];
//            next_char = s[word_length];
//
//            if(compressed_char < 0){ /* not again! */
//              skip = word_length-1;
//              preshift = 0;
//              continue;
//            }
//
//            index = preshift + compressed_char;
//            preshift = (Int4)((((Int8)( index )) * recip) >> 32);
//            skip--;
//          }
//
//          s--; /*undo the following increment*/
//          continue;
//        }
//        ...
//        index = preshift + compressed_char;
//        preshift = (Int4)((((Int8)( index )) * recip) >> 32);
//
//        /* if there are hits */
//        if (PV_TEST(pv, index, pv_array_bts)) {
//           Int4 s_off = (Int4)(s - subject->sequence);
//
//           CompressedLookupBackboneCell* backbone_cell =
//                                         lookup->backbone + index;
//
//           numhits = backbone_cell->num_used;
//
//           /* and there is enough space in the destination array */
//           if (numhits) {
//           if (numhits <= (array_size - totalhits)) {
// ```
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-518
// ```c
//         /* for each hit, */
//         for (i = 0; i < hits; ++i) {
//             Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//             Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
//
//             /* calculate the diagonal associated with this query-subject pair
//              */
//
//             diag_coord = (query_offset - subject_offset) & diag_mask;
// ```
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
    let prefetch_distance = x_prefetch_distance();
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
mod tests {
    use super::*;
    use crate::algorithm::tblastx::blast_aascan::BlastOffsetPair as OffsetPair;
    use crate::algorithm::tblastx::lookup::compressed::{
        build_blosum62_compressed_lookup, compressed_lookup_parameters,
        compute_encoded_compressed_index,
    };

    // A small deterministic generator (xorshift64).
    struct XorShift(u64);

    impl XorShift {
        fn below(&mut self, n: usize) -> usize {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            (x % n as u64) as usize
        }
    }

    // NCBISTDAA letters inside the 15-letter compressed alphabet and outside it (NULLB, X, U, *, O).
    const LETTERS: [u8; 23] = [
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 27,
    ];
    const OUTSIDE: [u8; 5] = [0, 21, 24, 25, 26];

    // The reference: `blastp_for_each_offset_pair`'s rounds of `scan_subject` with an array of
    // `array_size` pairs; returns the hits in order and the number of rounds.
    fn reference_hits(
        lookup: &BlastCompressedAaLookupTable,
        subject: &[u8],
        subject_length: usize,
        array_size: usize,
    ) -> (Vec<(i32, i32)>, usize) {
        let seq_ranges = [(0, subject_length as i32)];
        let mut scan_range = [0, seq_ranges[0].0, seq_ranges[0].1 - lookup.word_length];
        if scan_range[2] < scan_range[1] {
            scan_range[2] = scan_range[1];
        }
        let mut offset_pairs = vec![OffsetPair::default(); array_size];
        let mut hits = Vec::new();
        let mut rounds = 0;
        while scan_range[1] <= scan_range[2] {
            let count = lookup.scan_subject(
                subject,
                &seq_ranges,
                &mut offset_pairs,
                array_size as i32,
                &mut scan_range,
            );
            hits.extend(
                offset_pairs[..count as usize]
                    .iter()
                    .map(|pair| (pair.s_off as i32, pair.q_off as i32)),
            );
            rounds += 1;
        }
        (hits, rounds)
    }

    fn fused_hits(
        lookup: &BlastCompressedAaLookupTable,
        subject: &[u8],
        subject_length: usize,
    ) -> Vec<(i32, i32)> {
        let mut scratch = XScratch::default();
        let mut hits = Vec::new();
        x_scan_cells(lookup, subject, subject_length, &mut scratch, |s, cell| {
            hits.extend(cell.iter().map(|&q| (s, q)));
        });
        hits
    }

    // The reference's two-hit test (`blast_engine.rs`, the closure of `blastp_for_each_offset_pair`)
    // over `hits`, with `candidate` in place of the context check and the extension.
    #[allow(clippy::too_many_arguments)]
    fn reference_two_hit(
        hits: &[(i32, i32)],
        diag_array: &mut [DiagStruct],
        diag_mask: i32,
        diag_offset: i32,
        window: i32,
        word_size: i32,
        candidate: &mut dyn FnMut(i32, i32, i32, i32, &mut DiagStruct),
    ) {
        for &(subject_offset, query_offset) in hits {
            let diag_entry =
                &mut diag_array[((query_offset - subject_offset) & diag_mask) as usize];
            if diag_entry.flag() != 0 {
                let subject_plus_offset = subject_offset + diag_offset;
                if subject_plus_offset < diag_entry.last_hit() {
                    continue;
                }
                diag_entry.set_last_hit(subject_plus_offset);
                diag_entry.set_flag(0);
                continue;
            }
            let last_hit = diag_entry.last_hit() - diag_offset;
            let diff = subject_offset - last_hit;
            if diff >= window {
                diag_entry.set_last_hit(subject_offset + diag_offset);
                continue;
            }
            if diff < word_size {
                continue;
            }
            candidate(query_offset, subject_offset, last_hit, diff, diag_entry);
        }
    }

    // A stand-in for the extension: it logs the call and writes the entry as the reference closure
    // does (a flag and a later `last_hit` after a right extension, else this hit).
    fn stand_in(
        log: &mut Vec<(i32, i32, i32, i32)>,
        diag_offset: i32,
        query_offset: i32,
        subject_offset: i32,
        last_hit: i32,
        diff: i32,
        entry: &mut DiagStruct,
    ) {
        log.push((query_offset, subject_offset, last_hit, diff));
        match (query_offset ^ subject_offset) & 3 {
            0 => {
                entry.set_flag(1);
                entry.set_last_hit(subject_offset + 7 + diag_offset);
            }
            1 => {}
            _ => entry.set_last_hit(subject_offset + diag_offset),
        }
    }

    fn random_subject(rng: &mut XorShift, length: usize, outside_one_in: usize) -> Vec<u8> {
        let mut subject: Vec<u8> = (0..length)
            .map(|_| {
                if outside_one_in != 0 && rng.below(outside_one_in) == 0 {
                    OUTSIDE[rng.below(OUTSIDE.len())]
                } else {
                    LETTERS[rng.below(LETTERS.len())]
                }
            })
            .collect();
        subject.push(0);
        subject
    }

    // Runs the subjects one after the other through both loops from the same diagonal table (the
    // offset moves on by length + window after each subject, as `finish_subject` does) and compares
    // the hits, the candidates and the table.
    fn assert_same_two_hit(
        lookup: &BlastCompressedAaLookupTable,
        subjects: &[Vec<u8>],
        diag_len: usize,
        array_sizes: &[usize],
    ) {
        let window = 40;
        let word_size = lookup.word_length;
        let diag_mask = diag_len as i32 - 1;
        let mut fused_diag = vec![DiagStruct::default(); diag_len];
        let mut reference_diag = vec![DiagStruct::default(); diag_len];
        let mut diag_offset = window;
        let mut total_candidates = 0usize;
        for (n, subject) in subjects.iter().enumerate() {
            let length = subject.len() - 1;
            let fused = fused_hits(lookup, subject, length);
            let array_size = array_sizes[n % array_sizes.len()];
            let (reference, _) = reference_hits(lookup, subject, length, array_size);
            assert_eq!(fused, reference, "subject {n}: hits");

            let mut fused_log = Vec::new();
            x_scan_two_hit_compressed(
                lookup,
                subject,
                length,
                &mut fused_diag,
                diag_mask,
                diag_offset,
                window,
                word_size,
                |q, s, last_hit, diff, entry| {
                    stand_in(&mut fused_log, diag_offset, q, s, last_hit, diff, entry)
                },
            );
            let mut reference_log = Vec::new();
            reference_two_hit(
                &reference,
                &mut reference_diag,
                diag_mask,
                diag_offset,
                window,
                word_size,
                &mut |q, s, last_hit, diff, entry| {
                    stand_in(&mut reference_log, diag_offset, q, s, last_hit, diff, entry)
                },
            );
            assert_eq!(fused_log, reference_log, "subject {n}: candidates");
            assert!(
                fused_diag
                    .iter()
                    .zip(&reference_diag)
                    .all(|(a, b)| a.raw_bits() == b.raw_bits()),
                "subject {n}: diagonal table"
            );
            total_candidates += fused_log.len();
            diag_offset += length as i32 + window;
        }
        assert!(total_candidates > 0, "the test reaches the extension");
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1288-1301
    // ```c
    //     /* set word size and threshold information. The reciprocals
    //        below are 2^32 / (compressed alphabet size) */
    //     ...
    //     if (word_size == 6 || word_size == 5) {
    //         lookup->compressed_alphabet_size = 15;
    //         lookup->reciprocal_alphabet_size = 286331154;
    //     }
    //     else {
    //         lookup->compressed_alphabet_size = 10;
    //         lookup->reciprocal_alphabet_size = 429496730;
    //     }
    // ```
    // The resumed rounds of the reference scan give the rolling index only if the reciprocal gives
    // the exact quotient for every index the scan forms (module header).
    #[test]
    fn the_reciprocal_gives_the_exact_quotient_below_every_tables_size() {
        for word_size in 5..=7usize {
            let (alphabet, recip) = compressed_lookup_parameters(word_size).expect("table");
            let bound = (alphabet as i64).pow(word_size as u32);
            for index in 0..bound {
                let index = index as i32;
                assert_eq!(
                    compressed_reciprocal_divide(index, i64::from(recip)),
                    index / alphabet as i32,
                    "word {word_size} index {index}"
                );
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:825-844
    // ```c
    //     default:
    //     {
    //         ...
    //         Int4 cell_index = (num_entries -3) & COMPRESSED_HITS_CELL_MASK ;
    //
    //         if (cell_index == 0 ) { /* can't be empty => it's full  */
    //
    //             /* fetch next vacant cell */
    //         	CompressedOverflowCell * new_cell = s_CompressedListGetNewCell(lookup);
    //
    //             /* shuffle the pointers */
    //             new_cell->next = backbone_cell->payload.overflow_list.head;
    //             backbone_cell->payload.overflow_list.head = new_cell;
    //         }
    // ```
    // Cells of 1 to 20 hits (every position of the newest group), each word once in the subject and
    // between letters outside the alphabet: the fused loop gives each cell's hits in the order of
    // NCBI's scan, with arrays that make the reference stop before many cells.
    #[test]
    fn cells_of_one_to_twenty_hits_come_in_ncbis_group_order() {
        let mut rng = XorShift(0x9e37_79b9_7f4a_7c15);
        let mut lookup = BlastCompressedAaLookupTable::new(5, 0.0).expect("compressed lookup");
        let mut subject = Vec::new();
        let mut expected = Vec::new();
        for numhits in 1..=20usize {
            let (word, index) = loop {
                let word: Vec<u8> = (0..5).map(|_| LETTERS[rng.below(LETTERS.len())]).collect();
                let encoded: Vec<u8> = word
                    .iter()
                    .map(|&letter| lookup.compress_table[letter as usize])
                    .collect();
                let index = compute_encoded_compressed_index(5, &encoded, 15);
                if lookup.backbone[index].num_used == 0 {
                    break (word, index);
                }
            };
            for k in 0..numhits {
                lookup.add_word_hit(index, (1000 * numhits + 37 * k) as i32);
            }
            let order = lookup.query_offsets_in_scan_order(index);
            let s = subject.len() as i32;
            expected.extend(order.iter().map(|&q| (s, q)));
            subject.extend_from_slice(&word);
            subject.push(OUTSIDE[rng.below(OUTSIDE.len())]);
        }
        lookup.finalize();
        assert_eq!(lookup.longest_chain, 20);
        let length = subject.len();
        subject.push(0);
        let fused = fused_hits(&lookup, &subject, length);
        assert_eq!(fused, expected);
        for array_size in [20, 21, 27, 40, 4096 + 20] {
            let (reference, rounds) = reference_hits(&lookup, &subject, length, array_size);
            assert_eq!(fused, reference, "array {array_size}");
            if array_size < 40 {
                assert!(rounds > 5, "array {array_size}: {rounds} rounds");
            }
        }
        // The 8-hit cell of NCBI's example: the first three, the newest group, the older group.
        let eight: Vec<i32> = fused
            .iter()
            .filter(|(_, q)| q / 1000 == 8)
            .map(|&(_, q)| q - 8000)
            .collect();
        assert_eq!(eight, vec![0, 37, 74, 259, 111, 148, 185, 222]);
    }

    // Queries of few letters and a low threshold give cells of many hits; subjects of random
    // letters with letters outside the alphabet (none, rare, frequent, runs), short subjects, and
    // arrays from the longest cell up (many rounds) to the reference's size. Small diagonal tables
    // make hits of one subject offset share entries, so their order shows in the table.
    #[test]
    fn the_fused_loop_matches_the_reference_scan_and_two_hit_test() {
        let mut rng = XorShift(0xbb67_ae85_84ca_a73b);
        // Word size 5 only: the tables of 6 and 7 have 15^6 and 10^7 backbone cells.
        let word_size = 5usize;
        for (threshold, letters) in [
            (11.0, &[1u8, 11, 17, 10][..]),
            (9.0, &[1u8, 11, 17][..]),
            (13.0, &LETTERS[..]),
        ] {
            let mut query = Vec::new();
            let mut locations = Vec::new();
            for length in [120usize, 75, 200, 33] {
                let left = query.len() as i32;
                query.extend((0..length).map(|_| letters[rng.below(letters.len())]));
                locations.push((left, left + length as i32 - 1));
                query.push(0);
            }
            let lookup = build_blosum62_compressed_lookup(word_size, threshold, &query, &locations)
                .expect("compressed lookup");
            assert!(lookup.longest_chain > 5, "cells with overflow lists");
            let longest = lookup.longest_chain as usize;
            let mut subjects = Vec::new();
            for outside_one_in in [0usize, 50, 6, 2] {
                for _ in 0..6 {
                    let length = 1 + rng.below(300);
                    let mut subject = random_subject(&mut rng, length, outside_one_in);
                    // Some subjects copy a stretch of the query (many hits, many extensions).
                    if length > 60 && rng.below(2) == 0 {
                        let from = rng.below(query.len() - 60);
                        let to = rng.below(length - 50);
                        subject[to..to + 50].copy_from_slice(&query[from..from + 50]);
                    }
                    subjects.push(subject);
                }
            }
            for length in 1..=word_size + 1 {
                subjects.push(random_subject(&mut rng, length, 0));
            }
            let array_sizes = [
                longest,
                longest + 1,
                longest + 5,
                2 * longest,
                4096 + longest,
            ];
            for diag_len in [64usize, 1024, 4096] {
                assert_same_two_hit(&lookup, &subjects, diag_len, &array_sizes);
            }
        }
    }
}
