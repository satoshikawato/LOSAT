use super::*;
use crate::config::ScoringMatrix;
use crate::stats::lookup_protein_params_ungapped;
use crate::utils::matrix::{blosum62_score, ncbistdaa, BLASTAA_SIZE};

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1032-1114: the complement of the masks of one context, without the
// empty intervals at the ends.
#[test]
fn lookup_segments_are_the_complement_of_the_masks() {
    assert_eq!(build_lookup_segments_for_context(0, 10, &[]), vec![(0, 9)]);
    assert_eq!(
        build_lookup_segments_for_context(0, 10, &[(2, 4)]),
        vec![(0, 1), (4, 9)]
    );
    assert_eq!(
        build_lookup_segments_for_context(0, 10, &[(0, 3)]),
        vec![(3, 9)]
    );
    assert_eq!(
        build_lookup_segments_for_context(0, 10, &[(6, 10)]),
        vec![(0, 5)]
    );
    assert_eq!(
        build_lookup_segments_for_context(0, 12, &[(0, 2), (5, 7), (9, 12)]),
        vec![(2, 4), (7, 8)]
    );
    assert_eq!(
        build_lookup_segments_for_context(4, 10, &[(2, 4)]),
        vec![(4, 5), (8, 13)]
    );
    assert!(build_lookup_segments_for_context(0, 5, &[(0, 5)]).is_empty());
    assert!(build_lookup_segments_for_context(0, 0, &[]).is_empty());
}

// No NCBI counterpart (LOSAT_X_BPDBKEEP): the unthresholded arena cuts its slabs and long chains
// from one region while it has room (one region per 32 MiB on 64-bit glibc, else one per
// allocation): the blocks it gives do not overlap and keep what was written to them.
#[test]
fn unthresholded_slabs_are_cut_from_one_region() {
    let mut arena = XChainArena::new();
    arena.x_unthresholded = true;
    let sizes = [
        XChainArena::SLAB_INTS,
        1 << 15,
        XChainArena::SLAB_INTS,
        3,
        1 << 17,
    ];
    let mut blocks = Vec::new();
    for (i, &ints) in sizes.iter().enumerate() {
        let start = arena.allocate(ints);
        // SAFETY: the arena gave `ints` ints at `start`.
        unsafe { std::slice::from_raw_parts_mut(start, ints) }.fill(i as i32 + 1);
        blocks.push((start as usize, ints));
    }
    for (i, &(start, ints)) in blocks.iter().enumerate() {
        // SAFETY: as above; the arena is alive.
        let ints_read = unsafe { std::slice::from_raw_parts(start as *const i32, ints) };
        assert!(ints_read.iter().all(|&v| v == i as i32 + 1), "block {i}");
        for &(other, other_ints) in &blocks[i + 1..] {
            let end = start + ints * 4;
            assert!(
                end <= other || other + other_ints * 4 <= start,
                "blocks {i} overlap"
            );
        }
    }
    let regions = if crate::utils::x_mapped::X_MAPPED_ALLOCATOR {
        1
    } else {
        sizes.len()
    };
    assert_eq!(arena.x_unthresholded_slabs.len(), regions);
    arena.x_reset();
    assert!(arena.x_unthresholded_slabs.is_empty());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:87-129
fn make_poly_a_query_frame() -> QueryFrame {
    QueryFrame {
        frame: 1,
        aa_seq: vec![
            0,
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
            0,
        ],
        aa_seq_nomask: None,
        aa_len: 6,
        orig_len: 6,
        seg_masks: Vec::new(),
    }
}

fn make_query_frame(seq: Vec<u8>, seg_masks: Vec<(usize, usize)>) -> QueryFrame {
    let mut aa_seq = Vec::with_capacity(seq.len() + 2);
    aa_seq.push(0);
    aa_seq.extend_from_slice(&seq);
    aa_seq.push(0);
    QueryFrame {
        frame: 1,
        aa_seq,
        aa_seq_nomask: None,
        aa_len: seq.len(),
        orig_len: seq.len(),
        seg_masks,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2778-2803 (cited at
// `x_context_karlin`)
// The contexts of a query set with their Karlin blocks computed on a pool's threads
// (`x_lookup_contexts_for_profile_in_pool`) are those of the one-thread computation: every
// field, a context whose block fails (letters X only) included, for each matrix's bounds.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn contexts_in_pool_equal_the_one_thread_contexts() {
    let mut seed = 0x5eed_u64;
    let mut residue = || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        1 + ((seed >> 33) % 22) as u8
    };
    let mut queries: Vec<Vec<QueryFrame>> = (0..37)
        .map(|i| {
            let seq: Vec<u8> = (0..(5 + i * 13)).map(|_| residue()).collect();
            vec![make_query_frame(seq, vec![(1, 3)])]
        })
        .collect();
    queries.push(vec![make_poly_a_query_frame()]);
    // Only negative scores: no Karlin block (NCBI's `is_valid` false).
    queries.push(vec![make_query_frame(vec![ncbistdaa::X; 20], vec![])]);
    for matrix in [
        ScoringMatrix::Blosum62,
        ScoringMatrix::Blosum45,
        ScoringMatrix::Pam30,
    ] {
        let one = x_lookup_contexts_for_profile(&queries, false, matrix, 3).unwrap();
        let pooled = crate::utils::threading::with_search_pool(4, "tblastn", |pool| {
            Ok(x_lookup_contexts_for_profile_in_pool(
                &queries, false, matrix, 3, pool,
            ))
        })
        .unwrap()
        .unwrap();
        assert_eq!(one.len(), pooled.len());
        assert!(one.iter().any(|context| !context.is_valid));
        for (a, b) in one.iter().zip(&pooled) {
            assert_eq!(
                (a.q_idx, a.frame, a.aa_len, a.orig_len, a.frame_base),
                (b.q_idx, b.frame, b.aa_len, b.orig_len, b.frame_base)
            );
            assert_eq!((&a.aa_seq, &a.aa_seq_nomask), (&b.aa_seq, &b.aa_seq_nomask));
            assert_eq!(a.is_valid, b.is_valid);
            assert_eq!(
                format!("{:?}", a.karlin_params),
                format!("{:?}", b.karlin_params)
            );
        }
    }
}

fn build_lookup_for_test(queries: Vec<Vec<QueryFrame>>, threshold: i32) -> BlastAaLookupTable {
    let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let (lookup, _) = build_ncbi_lookup(&queries, threshold, &ungapped, true);
    lookup
}

fn prepare_lookup_for_test(queries: Vec<Vec<QueryFrame>>) -> PreparedLookupQuery {
    let ideal_params = blosum62_ideal_karlin_params();
    let std_comp = compute_std_aa_composition();
    prepare_lookup_query(
        &queries,
        ideal_params,
        &std_comp,
        LOOKUP_WORD_LENGTH,
        true,
        ScoringMatrix::Blosum62,
        (-4, 11),
    )
}

fn debruijn(alphabet_size: usize, order: usize) -> Vec<u8> {
    fn db(
        t: usize,
        p: usize,
        alphabet_size: usize,
        order: usize,
        a: &mut [usize],
        output: &mut Vec<u8>,
    ) {
        if t > order {
            if order.is_multiple_of(p) {
                for j in 1..=p {
                    output.push(a[j] as u8);
                }
            }
            return;
        }

        a[t] = a[t - p];
        db(t + 1, p, alphabet_size, order, a, output);
        for residue in (a[t - p] + 1)..alphabet_size {
            a[t] = residue;
            db(t + 1, t, alphabet_size, order, a, output);
        }
    }

    let mut a = vec![0usize; alphabet_size * order + 1];
    let mut output = Vec::new();
    db(1, 1, alphabet_size, order, &mut a, &mut output);
    let wraparound = output
        .iter()
        .copied()
        .take(order.saturating_sub(1))
        .collect::<Vec<_>>();
    output.extend(wraparound);
    output
}

fn brute_force_lookup_count(subject_word: [u8; 3], threshold: i32) -> usize {
    let mut count = 0usize;
    for q0 in 0..BLASTAA_SIZE {
        for q1 in 0..BLASTAA_SIZE {
            for q2 in 0..BLASTAA_SIZE {
                let score = blosum62_score(q0 as u8, subject_word[0])
                    + blosum62_score(q1 as u8, subject_word[1])
                    + blosum62_score(q2 as u8, subject_word[2]);
                if score >= threshold
                    || (q0 as u8 == subject_word[0]
                        && q1 as u8 == subject_word[1]
                        && q2 as u8 == subject_word[2])
                {
                    count += 1;
                }
            }
        }
    }
    count
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:311-314
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1098-1101
#[test]
fn test_prepare_lookup_query_shares_single_boundary_sentinel_and_context_offsets() {
    let prepared = prepare_lookup_for_test(vec![
        vec![make_query_frame(
            vec![ncbistdaa::A, ncbistdaa::A, ncbistdaa::A],
            Vec::new(),
        )],
        vec![make_query_frame(
            vec![ncbistdaa::C, ncbistdaa::C, ncbistdaa::C],
            Vec::new(),
        )],
    ]);

    assert_eq!(
        prepared.concat_query,
        vec![
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
            0,
            ncbistdaa::C,
            ncbistdaa::C,
            ncbistdaa::C,
            0,
        ]
    );
    assert_eq!(prepared.frame_bases, vec![0, 4]);
    assert_eq!(prepared.lookup_locations, vec![(0, 2), (4, 6)]);
    assert_eq!(prepared.contexts[0].frame_base, 0);
    assert_eq!(prepared.contexts[1].frame_base, 4);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1045-1057
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1067-1114
#[test]
fn test_prepare_lookup_query_converts_seg_masks_to_absolute_inclusive_segments() {
    let prepared = prepare_lookup_for_test(vec![
        vec![make_query_frame(
            vec![
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
            ],
            vec![(2, 4)],
        )],
        vec![make_query_frame(
            vec![ncbistdaa::C, ncbistdaa::C, ncbistdaa::C],
            Vec::new(),
        )],
    ]);

    assert_eq!(prepared.lookup_locations, vec![(0, 1), (4, 5), (7, 9)]);
    assert_eq!(prepared.skipped_seg_mask, 4);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2778-2791
#[test]
fn test_prepare_lookup_query_marks_failed_ungapped_ka_context_invalid() {
    let prepared = prepare_lookup_for_test(vec![vec![make_query_frame(
        vec![ncbistdaa::X, ncbistdaa::X, ncbistdaa::X],
        Vec::new(),
    )]]);

    assert_eq!(prepared.contexts.len(), 1);
    assert!(!prepared.contexts[0].is_valid);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:246-250
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:112-127
#[test]
fn test_build_ncbi_lookup_query_length_matches_ncbi_sequence_length() {
    let lookup = build_lookup_for_test(
        vec![
            vec![make_query_frame(
                vec![ncbistdaa::A, ncbistdaa::A, ncbistdaa::A],
                Vec::new(),
            )],
            vec![make_query_frame(
                vec![ncbistdaa::C, ncbistdaa::C, ncbistdaa::C],
                Vec::new(),
            )],
        ],
        0,
    );

    assert_eq!(lookup.query_length, 7);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:102-128
#[test]
fn test_build_ncbi_lookup_multiple_queries_do_not_cross_shared_boundary_sentinel() {
    let lookup = build_lookup_for_test(
        vec![
            vec![make_query_frame(
                vec![ncbistdaa::A, ncbistdaa::A, ncbistdaa::A],
                Vec::new(),
            )],
            vec![make_query_frame(
                vec![ncbistdaa::A, ncbistdaa::A, ncbistdaa::A],
                Vec::new(),
            )],
        ],
        0,
    );
    let aaa_index = encode_kmer_3(
        ncbistdaa::A as usize,
        ncbistdaa::A as usize,
        ncbistdaa::A as usize,
        lookup.charsize as usize,
    ) & lookup.mask as usize;

    assert_eq!(lookup.get_hits(aaa_index), &[0, 4]);
    assert_eq!(lookup.get_context_idx(0), 0);
    assert_eq!(lookup.get_context_idx(4), 1);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:340-357
#[test]
fn test_build_ncbi_lookup_preserves_overflow_query_offset_append_order() {
    let lookup = build_lookup_for_test(vec![vec![make_poly_a_query_frame()]], 0);
    let aaa_index = encode_kmer_3(
        ncbistdaa::A as usize,
        ncbistdaa::A as usize,
        ncbistdaa::A as usize,
        lookup.charsize as usize,
    ) & lookup.mask as usize;

    assert_eq!(lookup.backbone[aaa_index].num_used, 4);
    assert_eq!(&lookup.overflow[..4], &[0, 1, 2, 3]);
    assert_eq!(lookup.get_hits(aaa_index), &[0, 1, 2, 3]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:504-543
#[test]
fn test_build_ncbi_lookup_preserves_neighbor_query_offset_append_order() {
    let lookup = build_lookup_for_test(vec![vec![make_poly_a_query_frame()]], 11);
    let aaa_index = encode_kmer_3(
        ncbistdaa::A as usize,
        ncbistdaa::A as usize,
        ncbistdaa::A as usize,
        lookup.charsize as usize,
    ) & lookup.mask as usize;

    assert_eq!(lookup.backbone[aaa_index].num_used, 4);
    assert_eq!(&lookup.overflow[..4], &[0, 1, 2, 3]);
    assert_eq!(lookup.get_hits(aaa_index), &[0, 1, 2, 3]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:96-128
#[test]
fn test_build_ncbi_lookup_skips_words_crossing_invalid_residue() {
    let lookup = build_lookup_for_test(
        vec![vec![make_query_frame(
            vec![
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
                32,
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
            ],
            Vec::new(),
        )]],
        0,
    );
    let aaa_index = encode_kmer_3(
        ncbistdaa::A as usize,
        ncbistdaa::A as usize,
        ncbistdaa::A as usize,
        lookup.charsize as usize,
    ) & lookup.mask as usize;

    assert_eq!(lookup.get_hits(aaa_index), &[0, 4]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:92-100
#[test]
fn test_build_ncbi_lookup_respects_interval_boundaries() {
    let lookup = build_lookup_for_test(
        vec![vec![make_query_frame(
            vec![
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
            ],
            vec![(3, 4)],
        )]],
        0,
    );
    let aaa_index = encode_kmer_3(
        ncbistdaa::A as usize,
        ncbistdaa::A as usize,
        ncbistdaa::A as usize,
        lookup.charsize as usize,
    ) & lookup.mask as usize;

    assert_eq!(lookup.get_hits(aaa_index), &[0]);
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_aalookup.h:57 (`overflow_cursor`, an
// `Int4`, cited at `BackboneCell`). A cursor below 2^31 is NCBI's `Int4` in
// `entries[0]` (and 0 in `entries[1]`); a cursor of 2^31 or more, which NCBI's `Int4` cannot
// hold (`-threshold 3e9`: 21,952 entries per word position, so 97,827 positions fill 2^31),
// reads back whole, as every scanner reads it (`cursor`, `x_cursor_at`), where
// `entries[0] as usize` gave 0xffff_ffff_8000_0000.
#[test]
fn cursor_holds_positions_past_int4() {
    #[allow(unused_mut)]
    let mut cursors = vec![0usize, 7, (1 << 31) - 1, 1 << 31, u32::MAX as usize];
    #[cfg(target_pointer_width = "64")]
    cursors.extend([1usize << 32, (1usize << 32) + 5, usize::MAX >> 2]);
    for cursor in cursors {
        let mut cell = BackboneCell::default();
        cell.set_count(AA_HITS_PER_CELL + 1);
        XCell::set_cursor(&mut cell, cursor);
        assert_eq!(cell.cursor(), cursor, "{cursor:#x}");
        assert_eq!(XCell::cursor(&cell), cursor, "{cursor:#x}");
        // SAFETY: `cell` is a cell that nothing writes.
        assert_eq!(
            unsafe { BackboneCell::x_cursor_at(&cell) },
            cursor,
            "{cursor:#x}"
        );
        if cursor < 1 << 31 {
            assert_eq!(cell.entries, [cursor as i32, 0, 0]);
        }
    }
    let mut narrow = BackboneCell16::default();
    narrow.set_count(AA_HITS_PER_CELL + 1);
    narrow.set_cursor((1 << 31) + 3);
    assert_eq!(narrow.cursor(), (1 << 31) + 3);
    // SAFETY: `narrow` is a cell that nothing writes.
    assert_eq!(
        unsafe { BackboneCell16::x_cursor_at(&narrow) },
        (1 << 31) + 3
    );
}

// No NCBI counterpart: `get_hits` is a safe function of a table whose fields are `pub`, so an
// index past the backbone, or a cell whose run does not fit `overflow`, panics instead of
// reading outside the table (it read with `get_unchecked` until the SWEBOK 2026-11-04 audit,
// KA4).
#[test]
#[cfg_attr(panic = "abort", ignore = "catches a panic: the build aborts on panic")]
fn get_hits_panics_outside_the_table() {
    let lookup = build_lookup_for_test(
        vec![vec![make_query_frame(
            vec![ncbistdaa::A, ncbistdaa::A, ncbistdaa::A],
            vec![],
        )]],
        0,
    );
    let past = lookup.backbone.len();
    let read = |table: &BlastAaLookupTable, index: usize| {
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| table.get_hits(index).len()))
    };
    assert!(read(&lookup, past).is_err());
    let mut broken = lookup;
    broken.backbone[0].num_used = AA_HITS_PER_CELL as i32 + 1;
    broken.backbone[0].set_cursor(broken.overflow.len());
    assert!(read(&broken, 0).is_err());
    broken.backbone[0].num_used = -1;
    assert!(read(&broken, 0).is_err());
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/aalookup_unit_test.cpp:222-244 DebruijnSequenceTest
#[test]
fn test_build_ncbi_lookup_debruijn_exact_words_match_ncbi_counts() {
    let seq = debruijn(BLASTAA_SIZE, LOOKUP_WORD_LENGTH);
    let lookup = build_lookup_for_test(vec![vec![make_query_frame(seq.clone(), Vec::new())]], 0);
    let charsize = lookup.charsize as usize;
    let residue_mask = (1usize << charsize) - 1;

    for idx in 0..lookup.backbone.len() {
        let w0 = (idx >> (2 * charsize)) & residue_mask;
        let w1 = (idx >> charsize) & residue_mask;
        let w2 = idx & residue_mask;
        if w0 < BLASTAA_SIZE && w1 < BLASTAA_SIZE && w2 < BLASTAA_SIZE {
            assert_eq!(lookup.backbone[idx].num_used, 1, "valid idx={idx}");
        } else {
            assert_eq!(lookup.backbone[idx].num_used, 0, "invalid idx={idx}");
        }
    }
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/aalookup_unit_test.cpp:253-277 NeighboringWordsTest
#[test]
fn test_build_ncbi_lookup_representative_neighbor_counts_match_bruteforce() {
    let seq = debruijn(BLASTAA_SIZE, LOOKUP_WORD_LENGTH);
    let lookup = build_lookup_for_test(vec![vec![make_query_frame(seq, Vec::new())]], 11);
    let representative_words = [
        [ncbistdaa::A, ncbistdaa::A, ncbistdaa::A],
        [ncbistdaa::W, ncbistdaa::W, ncbistdaa::W],
        [ncbistdaa::A, ncbistdaa::R, ncbistdaa::N],
    ];

    for subject_word in representative_words {
        let index = encode_kmer_3(
            subject_word[0] as usize,
            subject_word[1] as usize,
            subject_word[2] as usize,
            lookup.charsize as usize,
        ) & lookup.mask as usize;
        assert_eq!(
            lookup.backbone[index].num_used as usize,
            brute_force_lookup_count(subject_word, 11),
            "subject_word={subject_word:?}"
        );
    }
}
