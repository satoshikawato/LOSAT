//! EXPERIMENT (LOSAT_X_BPSEED with LOSAT_X_BPCLUT): the fused seed loop of blastp-fast
//! (`x_seed_fused_compressed.rs`) on the flat compressed lookup table
//! (`tblastx/lookup/x_clut_direct.rs`).
//!
//! The loop is `x_scan_two_hit_compressed` with the flat table's scan (`for_each_cell`) in place of
//! `x_scan_cells`: the same subject offsets and cells in the same order, and each cell's hits as its
//! slice, which holds them in the order of `x_cell_hits` (NCBI's scan order; the permutation is
//! proved in `x_clut_direct.rs`). The diagonal test, the prefetch of a later hit's entry and the
//! closure are those of `x_scan_two_hit_compressed`, so the argument of `x_seed_fused_compressed.rs`
//! carries over unchanged.

use std::cell::Cell;

use super::x_seed_fused::{x_bp_scan_bf, x_prefetch, x_prefetch_distance, x_test_hit};
use crate::algorithm::tblastx::blast_extend::DiagStruct;
use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::XFlatCompressedLookup;

thread_local! {
    // No NCBI counterpart: the scan's list of (subject offset, backbone index), reused across
    // subjects.
    static X_CELLS: Cell<Option<Vec<(i32, u32)>>> = const { Cell::new(None) };
}

/// `x_seed_fused_compressed::x_scan_two_hit_compressed` on the flat table (same arguments).
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
pub(crate) fn x_scan_two_hit_compressed_flat<F>(
    lookup: &XFlatCompressedLookup,
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
    // EXPERIMENT (LOSAT_X_BPSCANBF): on a table with enough cells used (LOSAT_X_BPSCANBFMIN), the
    // offsets listed without a branch on the presence bit (`for_each_cell`).
    if x_bp_scan_bf() && lookup.x_scan_bf_pays() {
        x_scan_flat_bf(
            lookup,
            subject,
            subject_length,
            diag_array,
            diag_mask,
            diag_offset,
            window,
            word_size,
            &mut candidate,
        );
        return;
    }
    x_scan_flat_in::<false, F>(
        lookup,
        subject,
        subject_length,
        diag_array,
        diag_mask,
        diag_offset,
        window,
        word_size,
        &mut candidate,
    );
}

/// EXPERIMENT (LOSAT_X_BPSCANBF): [`x_scan_two_hit_compressed_flat`] with the offsets listed
/// without a branch on the presence bit. Out of line, and with a closure of its own type, so that
/// the instance with the branch is compiled as without this one (`x_seed_fused::x_scan_two_hit_bf`).
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn x_scan_flat_bf<F>(
    lookup: &XFlatCompressedLookup,
    subject: &[u8],
    subject_length: usize,
    diag_array: &mut [DiagStruct],
    diag_mask: i32,
    diag_offset: i32,
    window: i32,
    word_size: i32,
    candidate: &mut F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let mut own = |query_offset, subject_offset, last_hit, diff, entry: &mut DiagStruct| {
        candidate(query_offset, subject_offset, last_hit, diff, entry)
    };
    x_scan_flat_in::<true, _>(
        lookup,
        subject,
        subject_length,
        diag_array,
        diag_mask,
        diag_offset,
        window,
        word_size,
        &mut own,
    );
}

/// The loop of [`x_scan_two_hit_compressed_flat`]; `BRANCH_FREE` is `for_each_cell`'s.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn x_scan_flat_in<const BRANCH_FREE: bool, F>(
    lookup: &XFlatCompressedLookup,
    subject: &[u8],
    subject_length: usize,
    diag_array: &mut [DiagStruct],
    diag_mask: i32,
    diag_offset: i32,
    window: i32,
    word_size: i32,
    candidate: &mut F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    assert!(diag_mask >= 0 && diag_array.len() > diag_mask as usize);
    let prefetch_distance = x_prefetch_distance();
    let mut cells = X_CELLS.with(Cell::take).unwrap_or_default();
    lookup.for_each_cell::<BRANCH_FREE, _>(subject, subject_length, &mut cells, |s, cell_hits| {
        for (i, &query_offset) in cell_hits.iter().enumerate() {
            // Read ahead the entry of a later hit of this cell (a cache hint; see `x_scan_two_hit`).
            if let Some(&ahead) = cell_hits.get(i + prefetch_distance) {
                x_prefetch(diag_array, ((ahead - s) & diag_mask) as usize);
            }
            let diag_coord = ((query_offset - s) & diag_mask) as usize;
            // SAFETY: diag_coord <= diag_mask < diag_array.len() (asserted above).
            let entry = unsafe { diag_array.get_unchecked_mut(diag_coord) };
            if let Some((last_hit, diff)) = x_test_hit(entry, s, diag_offset, window, word_size) {
                candidate(query_offset, s, last_hit, diff, entry);
            }
        }
    });
    X_CELLS.with(|slot| slot.set(Some(cells)));
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::blastp::x_seed_fused_compressed::x_scan_two_hit_compressed;
    use crate::algorithm::tblastx::lookup::compressed::try_build_blosum62_compressed_lookup;
    use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::{
        x_try_build_flat_compressed_lookup, XFlatBuild,
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552
    // ```c
    //             diag_coord = (query_offset - subject_offset) & diag_mask;
    //             /* If the reset bit is set, an extension just happened. */
    //             if (diag_array[diag_coord].flag) {
    //             ...
    //                 if (diff >= window) {
    //             ...
    //                 if (diff < wordsize) {
    //                     continue;
    //                 }
    // ```
    // Queries of few letters and low thresholds (cells of many hits, several lookup segments),
    // subjects with letters outside the alphabet and copies of the query, small diagonal tables: the
    // loop on the flat table makes the candidates and leaves the diagonal table of the loop of
    // `x_seed_fused_compressed.rs` on the reference's table, subject after subject.
    #[test]
    fn the_flat_loop_is_the_fused_loop_on_the_reference_table() {
        let mut rng = XorShift(0x510e_527f_ade6_82d1);
        for (threshold, letters) in [
            (1100, &[1u8, 11, 17, 10][..]),
            (900, &[1u8, 11, 17][..]),
            (1300, &LETTERS[..]),
        ] {
            let mut query = Vec::new();
            let mut locations = Vec::new();
            for length in [120usize, 75, 200, 33] {
                let left = query.len() as i32;
                query.extend((0..length).map(|_| letters[rng.below(letters.len())]));
                // Two segments (a masked run between them).
                let cut = left + length as i32 / 2;
                locations.push((left, cut - 3));
                locations.push((cut, left + length as i32 - 1));
                query.push(0);
            }
            let tables = [0..4, 6..8];
            let reference =
                try_build_blosum62_compressed_lookup(5, threshold, &query, &locations, &tables)
                    .expect("word size 5")
                    .expect("not refused");
            let flat =
                match x_try_build_flat_compressed_lookup(5, threshold, &query, &locations, &tables)
                {
                    Some(XFlatBuild::Table(table)) => table,
                    _ => panic!("flat table"),
                };
            assert!(reference.longest_chain > 5, "cells with overflow lists");
            let mut subjects = Vec::new();
            for outside_one_in in [0usize, 50, 6, 2] {
                for _ in 0..6 {
                    let length = 1 + rng.below(300);
                    let mut subject: Vec<u8> = (0..length)
                        .map(|_| {
                            if outside_one_in != 0 && rng.below(outside_one_in) == 0 {
                                OUTSIDE[rng.below(OUTSIDE.len())]
                            } else {
                                LETTERS[rng.below(LETTERS.len())]
                            }
                        })
                        .collect();
                    if length > 60 && rng.below(2) == 0 {
                        let from = rng.below(query.len() - 60);
                        let to = rng.below(length - 50);
                        subject[to..to + 50].copy_from_slice(&query[from..from + 50]);
                    }
                    subject.push(0);
                    subjects.push(subject);
                }
            }
            for diag_len in [64usize, 1024, 4096] {
                let window = 40;
                let diag_mask = diag_len as i32 - 1;
                let mut flat_diag = vec![DiagStruct::default(); diag_len];
                let mut reference_diag = vec![DiagStruct::default(); diag_len];
                let mut diag_offset = window;
                let mut total = 0usize;
                for (n, subject) in subjects.iter().enumerate() {
                    let length = subject.len() - 1;
                    let mut flat_log = Vec::new();
                    x_scan_two_hit_compressed_flat(
                        &flat,
                        subject,
                        length,
                        &mut flat_diag,
                        diag_mask,
                        diag_offset,
                        window,
                        5,
                        |q, s, last_hit, diff, entry| {
                            stand_in(&mut flat_log, diag_offset, q, s, last_hit, diff, entry)
                        },
                    );
                    let mut reference_log = Vec::new();
                    x_scan_two_hit_compressed(
                        &reference,
                        subject,
                        length,
                        &mut reference_diag,
                        diag_mask,
                        diag_offset,
                        window,
                        5,
                        |q, s, last_hit, diff, entry| {
                            stand_in(&mut reference_log, diag_offset, q, s, last_hit, diff, entry)
                        },
                    );
                    assert_eq!(flat_log, reference_log, "subject {n}: candidates");
                    assert!(
                        flat_diag
                            .iter()
                            .zip(&reference_diag)
                            .all(|(a, b)| a.raw_bits() == b.raw_bits()),
                        "subject {n}: diagonal table"
                    );
                    total += flat_log.len();
                    diag_offset += length as i32 + window;
                }
                assert!(total > 0, "the test reaches the extension");
            }
        }
    }
}
