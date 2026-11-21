use super::*;
use crate::algorithm::blastp::x_seed_fused_compressed::x_scan_two_hit_compressed;
use crate::algorithm::tblastx::lookup::compressed::try_build_blosum62_compressed_lookup;
use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::{
    x_try_build_flat_compressed_lookup, XFlatBuild, X_BLASTP,
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
// Queries of few letters and low thresholds (cells of many hits, several lookup segments),
// subjects with letters outside the alphabet and copies of the query, small diagonal tables: the
// loop on the flat table makes the candidates and leaves the diagonal table of the loop of
// `x_seed_fused_compressed.rs` on the reference's table, subject after subject, with the list
// of the branch and the branch-free one (LOSAT_X_BPCLUTSCANBF), and so do the loops after
// the list: the branch-free test (`x_scan_flat_fseed_loop`), the staged one (`x_scan_flat_lean`)
// and the chooser between them (`x_scan_flat_fseed`).
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
            match x_try_build_flat_compressed_lookup(5, threshold, &query, &locations, &tables) {
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
            let mut fseed_diag = vec![DiagStruct::default(); diag_len];
            let mut chosen_diag = vec![DiagStruct::default(); diag_len];
            let mut lean_diag = vec![DiagStruct::default(); diag_len];
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
                // The branch-free test (the subjects are in XHitBf's range).
                let bf = XHitBf::new(diag_offset, length, window, 5).expect("in range");
                let mut fseed_log = Vec::new();
                x_scan_flat_fseed_loop(
                    &flat,
                    subject,
                    length,
                    &mut fseed_diag,
                    diag_mask,
                    diag_offset,
                    window,
                    5,
                    bf,
                    &mut |q, s, last_hit, diff, entry: &mut DiagStruct| {
                        stand_in(&mut fseed_log, diag_offset, q, s, last_hit, diff, entry)
                    },
                );
                // The staged loop on the subjects of a word or more (they keep their
                // sentinel); a shorter one takes the loop of the branch-free test.
                assert_eq!(flat.x_pairs_fit(subject, length), length >= 5);
                let mut lean_log = Vec::new();
                if length >= 5 {
                    // SAFETY: x_pairs_fit holds (asserted above).
                    unsafe {
                        x_scan_flat_lean::<X_BLASTP, _>(
                            &flat,
                            subject,
                            length,
                            &mut lean_diag,
                            diag_mask,
                            diag_offset,
                            window,
                            5,
                            bf,
                            &mut |q, s, last_hit, diff, entry: &mut DiagStruct| {
                                stand_in(&mut lean_log, diag_offset, q, s, last_hit, diff, entry)
                            },
                        )
                    };
                } else {
                    x_scan_flat_fseed_loop(
                        &flat,
                        subject,
                        length,
                        &mut lean_diag,
                        diag_mask,
                        diag_offset,
                        window,
                        5,
                        bf,
                        &mut |q, s, last_hit, diff, entry: &mut DiagStruct| {
                            stand_in(&mut lean_log, diag_offset, q, s, last_hit, diff, entry)
                        },
                    );
                }
                // The chooser, on every subject.
                let mut chosen_log = Vec::new();
                x_scan_flat_fseed::<X_BLASTP, _>(
                    &flat,
                    subject,
                    length,
                    &mut chosen_diag,
                    diag_mask,
                    diag_offset,
                    window,
                    5,
                    bf,
                    &mut |q, s, last_hit, diff, entry: &mut DiagStruct| {
                        stand_in(&mut chosen_log, diag_offset, q, s, last_hit, diff, entry)
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
                assert_eq!(
                    fseed_log, reference_log,
                    "subject {n}: branch-free test candidates"
                );
                assert_eq!(lean_log, reference_log, "subject {n}: staged candidates");
                assert_eq!(chosen_log, reference_log, "subject {n}: chosen candidates");
                for (diag, form) in [
                    (&flat_diag, "branch"),
                    (&bf_diag, "branch-free"),
                    (&fseed_diag, "branch-free test"),
                    (&lean_diag, "staged"),
                    (&chosen_diag, "chosen"),
                ] {
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

// The cell whose slice ends the hits array is copied hit by hit (no run of X_RUN reads past the
// array) by the staged loop (`x_scan_flat_lean`, `x_lean_copy`), with the outcomes of the loop
// of the branch-free test.
#[test]
fn the_lean_loop_copies_the_last_cell_of_the_hits() {
    let mut rng = XorShift(0x1f83_d9ab_fb41_bd6b);
    let mut query = Vec::new();
    let mut locations = Vec::new();
    for length in [400usize, 300] {
        let left = query.len() as i32;
        query.extend((0..length).map(|_| LETTERS[rng.below(LETTERS.len())]));
        locations.push((left, left + length as i32 - 1));
        query.push(0);
    }
    let flat = match x_try_build_flat_compressed_lookup(
        5,
        1100,
        &query,
        &locations,
        &[0..locations.len()],
    ) {
        Some(XFlatBuild::Table(table)) => table,
        _ => panic!("flat table"),
    };
    // The last cell with hits, its slice at the end of `hits`, and a word that reaches it: the
    // letter of digit k (index = sum of digit k * 15^k) sits k letters into the word.
    let last = (0..flat.starts.len() - 1)
        .rev()
        .find(|&i| flat.starts[i + 1] > flat.starts[i])
        .expect("a used cell");
    assert_eq!(flat.starts[last + 1] as usize, flat.hits.len());
    // Its last run starts within X_RUN of the end, so it is copied hit by hit.
    assert_ne!(
        (flat.hits.len() - flat.starts[last] as usize) % X_RUN,
        0,
        "the last cell's runs end hit by hit"
    );
    let scale = 15i32.pow(4);
    let letter = |digit: i32| {
        (0u8..=255)
            .find(|&b| flat.scaled_compress_table_u8[usize::from(b)] == digit * scale)
            .expect("a letter for every digit")
    };
    let word: Vec<u8> = (0..5)
        .map(|k| letter(last as i32 / 15i32.pow(k) % 15))
        .collect();
    let mut subject = Vec::new();
    for _ in 0..6 {
        subject.extend((0..7).map(|_| LETTERS[rng.below(LETTERS.len())]));
        subject.extend_from_slice(&word);
    }
    let length = subject.len();
    subject.push(0);
    assert!(flat.x_pairs_fit(&subject, length));
    let mut cells = Vec::new();
    // SAFETY: x_pairs_fit holds (asserted above).
    unsafe { flat.x_list_pairs_lean::<X_BLASTP>(&subject, length, &mut cells) };
    assert!(
        cells.iter().any(|&(_, index)| index as usize == last),
        "the subject reaches it"
    );
    let (window, diag_len, diag_offset) = (40, 4096usize, 40);
    let diag_mask = diag_len as i32 - 1;
    let bf = XHitBf::new(diag_offset, length, window, 5).expect("in range");
    let mut got_diag = vec![DiagStruct::default(); diag_len];
    let mut want_diag = got_diag.clone();
    let (mut got, mut want) = (Vec::new(), Vec::new());
    // SAFETY: x_pairs_fit holds (asserted above).
    unsafe {
        x_scan_flat_lean::<X_BLASTP, _>(
            &flat,
            &subject,
            length,
            &mut got_diag,
            diag_mask,
            diag_offset,
            window,
            5,
            bf,
            &mut |q, s, last_hit, diff, entry: &mut DiagStruct| {
                stand_in(&mut got, diag_offset, q, s, last_hit, diff, entry)
            },
        )
    };
    x_scan_flat_fseed_loop(
        &flat,
        &subject,
        length,
        &mut want_diag,
        diag_mask,
        diag_offset,
        window,
        5,
        bf,
        &mut |q, s, last_hit, diff, entry: &mut DiagStruct| {
            stand_in(&mut want, diag_offset, q, s, last_hit, diff, entry)
        },
    );
    assert_eq!(got, want, "candidates");
    assert!(
        got_diag
            .iter()
            .zip(&want_diag)
            .all(|(a, b)| a.raw_bits() == b.raw_bits()),
        "diagonal table"
    );
}

// The staged loop (`x_scan_flat_lean`) has the outcomes of the loop of the branch-free test
// (`x_scan_flat_fseed_loop`) on subjects that stage several times X_STAGE_CAP hits (the copy
// stops when the stage is full and goes on from there, within a cell too) and on short ones, on
// tables of few and many hits per cell (runs whose second half is read again and runs that read
// it).
#[test]
fn the_lean_loop_is_the_branch_free_test_loop() {
    let mut rng = XorShift(0x3c6e_f372_fe94_f82b);
    for (threshold, letters) in [(1100, &LETTERS[..]), (900, &[1u8, 11, 17, 10, 5][..])] {
        let mut query = Vec::new();
        let mut locations = Vec::new();
        for length in [700usize, 500, 800] {
            let left = query.len() as i32;
            query.extend((0..length).map(|_| letters[rng.below(letters.len())]));
            locations.push((left, left + length as i32 - 1));
            query.push(0);
        }
        let flat = match x_try_build_flat_compressed_lookup(
            5,
            threshold,
            &query,
            &locations,
            &[0..locations.len()],
        ) {
            Some(XFlatBuild::Table(table)) => table,
            _ => panic!("flat table"),
        };
        // The staging of the subjects' cells, replayed: a flush within a cell, runs with 4 hits
        // or fewer left (the second half read from `k` again) and runs with more.
        let (mut within, mut few, mut more) = (false, false, false);
        let window = 40;
        let diag_len = 4096usize;
        let diag_mask = diag_len as i32 - 1;
        let mut total = 0usize;
        for length in [5usize, 9, 50, 300, 1017, 3000, 5000] {
            let mut subject: Vec<u8> = (0..length)
                .map(|_| {
                    if rng.below(40) == 0 {
                        OUTSIDE[rng.below(OUTSIDE.len())]
                    } else {
                        letters[rng.below(letters.len())]
                    }
                })
                .collect();
            if length > 50 {
                let from = rng.below(query.len() - 50);
                subject[..40].copy_from_slice(&query[from..from + 40]);
            }
            subject.push(0);
            assert!(flat.x_pairs_fit(&subject, length));
            let (mut hits_staged, mut staged) = (0, 0);
            let mut cells = Vec::new();
            flat.for_each_cell::<true, _>(&subject, length, &mut cells, |_, hits| {
                hits_staged += hits.len();
                let mut k = 0;
                while k < hits.len() {
                    if staged + X_RUN > X_STAGE_CAP {
                        within |= k > 0;
                        staged = 0;
                    }
                    let left = hits.len() - k;
                    few |= left <= 4;
                    more |= left > 4;
                    staged += left.min(X_RUN);
                    k += left.min(X_RUN);
                }
            });
            if length == 5000 {
                assert!(
                    hits_staged > 3 * X_STAGE_CAP,
                    "{hits_staged} hits: several flushes"
                );
            }
            let diag_offset = window;
            let bf = XHitBf::new(diag_offset, length, window, 5).expect("in range");
            let mut want_diag = vec![DiagStruct::default(); diag_len];
            let mut got_diag = want_diag.clone();
            let (mut want, mut got) = (Vec::new(), Vec::new());
            x_scan_flat_fseed_loop(
                &flat,
                &subject,
                length,
                &mut want_diag,
                diag_mask,
                diag_offset,
                window,
                5,
                bf,
                &mut |q, s, last_hit, diff, entry: &mut DiagStruct| {
                    stand_in(&mut want, diag_offset, q, s, last_hit, diff, entry)
                },
            );
            // SAFETY: x_pairs_fit holds (asserted above).
            unsafe {
                x_scan_flat_lean::<X_BLASTP, _>(
                    &flat,
                    &subject,
                    length,
                    &mut got_diag,
                    diag_mask,
                    diag_offset,
                    window,
                    5,
                    bf,
                    &mut |q, s, last_hit, diff, entry: &mut DiagStruct| {
                        stand_in(&mut got, diag_offset, q, s, last_hit, diff, entry)
                    },
                )
            };
            assert_eq!(
                got, want,
                "threshold {threshold}, length {length}: candidates"
            );
            assert!(
                got_diag
                    .iter()
                    .zip(&want_diag)
                    .all(|(a, b)| a.raw_bits() == b.raw_bits()),
                "threshold {threshold}, length {length}: diagonal table"
            );
            total += want.len();
        }
        assert!(
            within && few && more,
            "threshold {threshold}: a flush within a cell and runs of either kind"
        );
        assert!(
            total > 0,
            "threshold {threshold}: the test reaches the extension"
        );
    }
}

// The choice of loop: on a table with the branch-free list the fused loop
// (`x_scan_two_hit_compressed_flat`) gives the candidates and the diagonal table of
// `x_scan_flat_bf`, on subjects where `XHitBf::fits` holds (the loops of `x_scan_flat_fseed`)
// and at the first diagonal offset where it does not (`x_scan_flat_bf` itself).
#[test]
fn the_fused_loop_takes_either_loop_with_the_same_outcomes() {
    let mut rng = XorShift(0x9b05_688c_2b3e_6c1f);
    let mut query = Vec::new();
    let mut locations = Vec::new();
    for length in [700usize, 500, 800] {
        let left = query.len() as i32;
        query.extend((0..length).map(|_| LETTERS[rng.below(LETTERS.len())]));
        locations.push((left, left + length as i32 - 1));
        query.push(0);
    }
    let flat = match x_try_build_flat_compressed_lookup(
        5,
        1100,
        &query,
        &locations,
        &[0..locations.len()],
    ) {
        Some(XFlatBuild::Table(table)) => table,
        _ => panic!("flat table"),
    };
    let window = 40;
    let diag_len = 4096usize;
    let diag_mask = diag_len as i32 - 1;
    let mut total = 0usize;
    // 3,000 letters stage several times X_STAGE_CAP hits (the staged loop's flushes).
    for length in [50usize, 300, 3000] {
        let mut subject: Vec<u8> = (0..length)
            .map(|_| LETTERS[rng.below(LETTERS.len())])
            .collect();
        let from = rng.below(query.len() - 50);
        subject[..40].copy_from_slice(&query[from..from + 40]);
        subject.push(0);
        if length == 3000 {
            let mut staged = 0;
            let mut cells = Vec::new();
            flat.for_each_cell::<true, _>(&subject, length, &mut cells, |_, hits| {
                staged += hits.len();
            });
            assert!(staged > 3 * X_STAGE_CAP, "{staged} hits: several flushes");
        }
        let edge = (1 << 30) - 2 * window - length as i32;
        for (diag_offset, fits) in [(window, true), (edge - 1, true), (edge, false)] {
            assert_eq!(XHitBf::fits(diag_offset, length, window, 5), fits);
            let mut got_diag = vec![DiagStruct::default(); diag_len];
            let mut want_diag = got_diag.clone();
            let mut got = Vec::new();
            x_scan_two_hit_compressed_flat::<X_BLASTP, _>(
                &flat,
                &subject,
                length,
                &mut got_diag,
                diag_mask,
                diag_offset,
                window,
                5,
                |q, s, last_hit, diff, entry| {
                    stand_in(&mut got, diag_offset, q, s, last_hit, diff, entry)
                },
            );
            let mut want = Vec::new();
            x_scan_flat_bf(
                &flat,
                &subject,
                length,
                &mut want_diag,
                diag_mask,
                diag_offset,
                window,
                5,
                &mut |q, s, last_hit, diff, entry: &mut DiagStruct| {
                    stand_in(&mut want, diag_offset, q, s, last_hit, diff, entry)
                },
            );
            assert_eq!(
                got, want,
                "length {length}, offset {diag_offset}: candidates"
            );
            assert!(
                got_diag
                    .iter()
                    .zip(&want_diag)
                    .all(|(a, b)| a.raw_bits() == b.raw_bits()),
                "length {length}, offset {diag_offset}: diagonal table"
            );
            total += got.len();
        }
    }
    assert!(total > 0, "the test reaches the extension");
}
