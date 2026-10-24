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
//!
//! EXPERIMENT (LOSAT_X_BPCLUTSCANBF): on a table with at least LOSAT_X_BPCLUTSCANBFMIN used cells
//! per 10,000 (`XFlatCompressedLookup::x_scan_bf_pays`) the scan lists the offsets without a branch
//! on the presence bit (`for_each_cell::<true, _>`, the form of LOSAT_X_BPSCANBF's
//! `x_seed_fused::x_scan_cells_bf`); the list, and so the cells and their order, are the same. That
//! instance is out of line with a closure of its own type (`x_scan_flat_bf`), so the instance with
//! the branch, which the other tables take, is compiled as without it. A switch of the strict set
//! (`utils/xswitch.rs`): on unless it is `0`.
//!
//! EXPERIMENT (LOSAT_X_BPFSEED): on a subject where `XHitBf::fits` holds, the loop takes three
//! changes (`x_scan_flat_fseed`): each hit goes first through blastp's branch-free two-hit test
//! (`XHitBf::test`, LOSAT_X_BPHITBF, whose argument in `x_seed_fused.rs` holds here: the diagonal
//! table and the extension's stores are those of blastp's loop), and a flag-1 entry or a
//! candidate goes on to `x_test_hit`, unchanged; no hit reads ahead a later hit's diagonal entry
//! (a cache hint); and the scan reads each listed cell's slice bounds as it finds the cell
//! (`XFlatCompressedLookup::x_for_each_range`), so the visits read no `starts` entry. The cells,
//! their order, the hits and their order are the same, the table holds after each hit what
//! `x_test_hit` leaves there, and the extension is called for the same hits with the same values
//! and entries. The instance is out of line with a closure of its own type, so the loops without
//! it keep their code; a subject outside `XHitBf::fits` runs them. A switch of the strict set: on
//! unless it is `0`.

use std::cell::Cell;

use super::x_seed_fused::{x_prefetch, x_prefetch_distance, x_test_hit, XHitBf};
use crate::algorithm::tblastx::blast_extend::DiagStruct;
use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::XFlatCompressedLookup;

thread_local! {
    // No NCBI counterpart: the scan's list of (subject offset, backbone index), reused across
    // subjects.
    static X_CELLS: Cell<Option<Vec<(i32, u32)>>> = const { Cell::new(None) };
}

thread_local! {
    // No NCBI counterpart: LOSAT_X_BPFSEED's list of (subject offset, slice start, slice end),
    // reused across subjects.
    static X_RANGES: Cell<Option<Vec<(i32, u32, u32)>>> = const { Cell::new(None) };
}

// No NCBI counterpart: reads the LOSAT_X_BPFSEED switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPFSEED): a strict-set switch, on unless it is `0`.
pub(crate) fn x_bp_fseed() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPFSEED"))
}

// No NCBI counterpart: reads the LOSAT_X_BPCLUTSCANBF switch once; it does not change any value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_BPCLUTSCANBF): a strict-set switch, on unless it is `0`.
pub(crate) fn x_bp_clut_scan_bf() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPCLUTSCANBF"))
}

// No NCBI counterpart: chooses LOSAT's list form; both list the same offsets.
/// EXPERIMENT (LOSAT_X_BPCLUTSCANBF): whether the fused loop lists `lookup`'s offsets without a
/// branch on the presence bit: the switch is on and the table has enough cells used.
pub(crate) fn x_clut_scan_bf(lookup: &XFlatCompressedLookup) -> bool {
    x_bp_clut_scan_bf() && lookup.x_scan_bf_pays()
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
    // EXPERIMENT (LOSAT_X_BPFSEED): the branch-free two-hit test, no read-ahead of the diagonal
    // entries and the slice bounds read by the scan, where the test gives the reference's
    // outcomes (module header).
    if x_bp_fseed() && x_fseed_fits(diag_offset, subject_length, window, word_size) {
        if x_clut_scan_bf(lookup) {
            x_scan_flat_fseed::<true, F>(
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
        } else {
            x_scan_flat_fseed::<false, F>(
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
        return;
    }
    // EXPERIMENT (LOSAT_X_BPCLUTSCANBF): on a table with enough cells used, the offsets listed
    // without a branch on the presence bit (`for_each_cell`).
    if x_clut_scan_bf(lookup) {
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

/// EXPERIMENT (LOSAT_X_BPCLUTSCANBF): [`x_scan_two_hit_compressed_flat`] with the offsets listed
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

// No NCBI counterpart: chooses LOSAT's loop; both loops give the reference's outcomes.
/// EXPERIMENT (LOSAT_X_BPFSEED): whether the subject takes [`x_scan_flat_fseed`] ([`XHitBf::fits`]),
/// recorded for `LOSAT_X_REACH`.
fn x_fseed_fits(diag_offset: i32, subject_length: usize, window: i32, word_size: i32) -> bool {
    let fits = XHitBf::fits(diag_offset, subject_length, window, word_size);
    crate::utils::xreach::x_reach!("LOSAT_X_BPFSEED", fits);
    fits
}

/// EXPERIMENT (LOSAT_X_BPFSEED): [`x_scan_two_hit_compressed_flat`] on a subject where
/// [`XHitBf::fits`] holds, with the branch-free two-hit test first, no read-ahead of the diagonal
/// entries, and the cells listed with their slice bounds
/// ([`XFlatCompressedLookup::x_for_each_range`]; `BRANCH_FREE` is its list form, LOSAT_X_BPCLUTSCANBF).
/// Out of line, and with a closure of its own type, so that the loops without it are compiled as
/// without this one (see [`x_scan_flat_bf`]).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-518 (quoted at
// `x_scan_two_hit_compressed_flat`)
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn x_scan_flat_fseed<const BRANCH_FREE: bool, F>(
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
    assert!(diag_mask >= 0 && diag_array.len() > diag_mask as usize);
    // The caller checked `XHitBf::fits` for the subject.
    let bf = XHitBf::in_range(diag_offset, window, word_size);
    let mut ranges = X_RANGES.with(Cell::take).unwrap_or_default();
    lookup.x_for_each_range::<BRANCH_FREE, _>(
        subject,
        subject_length,
        &mut ranges,
        |s, cell_hits| {
            let here = bf.here(s);
            for &query_offset in cell_hits {
                let diag_coord = ((query_offset - s) & diag_mask) as usize;
                // SAFETY: diag_coord <= diag_mask < diag_array.len() (asserted above).
                let entry = unsafe { diag_array.get_unchecked_mut(diag_coord) };
                if bf.test(entry, here) {
                    continue;
                }
                if let Some((last_hit, diff)) = x_test_hit(entry, s, diag_offset, window, word_size)
                {
                    own(query_offset, s, last_hit, diff, entry);
                }
            }
        },
    );
    X_RANGES.with(|slot| slot.set(Some(ranges)));
}

/// The loop of [`x_scan_two_hit_compressed_flat`]; `BRANCH_FREE` is `for_each_cell`'s.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-518 (quoted at
// `x_scan_two_hit_compressed_flat`)
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
    // `x_seed_fused_compressed.rs` on the reference's table, subject after subject, with the list
    // of the branch and the branch-free one (LOSAT_X_BPCLUTSCANBF).
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
                let mut bf_diag = vec![DiagStruct::default(); diag_len];
                let mut reference_diag = vec![DiagStruct::default(); diag_len];
                let mut diag_offset = window;
                let mut total = 0usize;
                for (n, subject) in subjects.iter().enumerate() {
                    let length = subject.len() - 1;
                    // Both instances, whichever the switch and the table's density choose.
                    let mut flat_log = Vec::new();
                    x_scan_flat_in::<false, _>(
                        &flat,
                        subject,
                        length,
                        &mut flat_diag,
                        diag_mask,
                        diag_offset,
                        window,
                        5,
                        &mut |q, s, last_hit, diff, entry: &mut DiagStruct| {
                            stand_in(&mut flat_log, diag_offset, q, s, last_hit, diff, entry)
                        },
                    );
                    let mut bf_log = Vec::new();
                    x_scan_flat_bf(
                        &flat,
                        subject,
                        length,
                        &mut bf_diag,
                        diag_mask,
                        diag_offset,
                        window,
                        5,
                        &mut |q, s, last_hit, diff, entry: &mut DiagStruct| {
                            stand_in(&mut bf_log, diag_offset, q, s, last_hit, diff, entry)
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
                    assert_eq!(bf_log, reference_log, "subject {n}: branch-free candidates");
                    for (diag, form) in [(&flat_diag, "branch"), (&bf_diag, "branch-free")] {
                        assert!(
                            diag.iter()
                                .zip(&reference_diag)
                                .all(|(a, b)| a.raw_bits() == b.raw_bits()),
                            "subject {n}: diagonal table ({form})"
                        );
                    }
                    total += flat_log.len();
                    diag_offset += length as i32 + window;
                }
                assert!(total > 0, "the test reaches the extension");
            }
        }
    }
}
