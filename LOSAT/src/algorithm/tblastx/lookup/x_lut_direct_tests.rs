use super::*;

pub(super) fn frame(seq: Vec<u8>, seg_masks: Vec<(usize, usize)>, frame: i8) -> QueryFrame {
    let mut aa_seq = Vec::with_capacity(seq.len() + 2);
    aa_seq.push(0);
    aa_seq.extend_from_slice(&seq);
    aa_seq.push(0);
    QueryFrame {
        frame,
        aa_len: seq.len(),
        orig_len: seq.len(),
        aa_seq,
        aa_seq_nomask: None,
        seg_masks,
    }
}

pub(super) fn residues(seed: &mut u64, n: usize) -> Vec<u8> {
    (0..n)
        .map(|_| {
            *seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            // NCBISTDAA 1..=27 (every letter of the alphabet but the gap)
            1 + ((*seed >> 33) % 27) as u8
        })
        .collect()
}

fn prepared(queries: &[Vec<QueryFrame>], matrix: ScoringMatrix) -> PreparedLookupQuery {
    let bounds = if matrix == ScoringMatrix::Blosum62 {
        (-4, 11)
    } else {
        matrix_score_bounds(matrix)
    };
    prepare_lookup_query(
        queries,
        ideal_karlin_params_for_matrix(matrix, bounds),
        &compute_std_aa_composition(),
        LOOKUP_WORD_LENGTH,
        true,
        matrix,
        bounds,
    )
}

/// [`x_build_direct_pooled`] asked for `Int4` hits (no cutoff): the table of every program
/// before 16-bit hits.
fn pooled_wide(
    prepared: PreparedLookupQuery,
    threshold: i32,
    matrix: ScoringMatrix,
    parts: usize,
) -> std::result::Result<BlastAaLookupTable, PreparedLookupQuery> {
    x_build_direct_pooled(prepared, threshold, matrix, 3, parts, None).map(|table| match table {
        XBxTable::Wide(table) => table,
        XBxTable::Narrow(_) => panic!("a build without a cutoff takes Int4 hits"),
    })
}

fn same(queries: &[Vec<QueryFrame>], threshold: i32, matrix: ScoringMatrix) -> usize {
    let reference = build_lookup_from_prepared(prepared(queries, matrix), threshold, matrix, 3).0;
    let direct = match x_build_direct(prepared(queries, matrix), threshold, matrix, 3) {
        Ok(table) => table,
        Err(_) => panic!("in-place builder declined"),
    };
    x_assert_same_table(&direct, &reference);
    reference.overflow.len()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:504-509,580-590
#[test]
fn in_place_table_equals_reference_for_thresholds_masks_and_overflow() {
    let mut seed = 7u64;
    let mut low_complexity = vec![1u8; 300]; // poly-A: cells far above AA_HITS_PER_CELL
    low_complexity.extend(residues(&mut seed, 40));
    let queries = vec![
        vec![
            frame(residues(&mut seed, 700), vec![(10, 30), (400, 401)], 1),
            frame(residues(&mut seed, 2), vec![], 2),
            frame(residues(&mut seed, 3), vec![], 3),
            frame(low_complexity, vec![(250, 299)], -1),
            frame(residues(&mut seed, 650), vec![(0, 0)], -2),
            frame(Vec::new(), vec![], -3),
        ],
        vec![frame(residues(&mut seed, 1200), vec![], 1)],
    ];
    let mut overflow_seen = false;
    for threshold in [0, 1, 11, 12, 13, 19, 30, 1000] {
        overflow_seen |= same(&queries, threshold, ScoringMatrix::Blosum62) > 0;
    }
    assert!(overflow_seen, "the fixture must reach the overflow array");
    // A second query set on the same thread reads the neighbour lists memoised above.
    let other = vec![vec![frame(residues(&mut seed, 900), vec![(5, 50)], 1)]];
    same(&other, 12, ScoringMatrix::Blosum62);
    // Another matrix: the memo key changes and the lists are recomputed.
    same(&other, 12, ScoringMatrix::Blosum45);
    same(&queries, 11, ScoringMatrix::Pam30);
    same(&queries, 12, ScoringMatrix::Blosum62);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:461-467 (cited at
// `x_build_direct`). TBLASTN's lean build (`x_build_lean`) gives the reference table for every
// threshold, mask and matrix at word lengths 2 and 3 (LOSAT's TBLASTN: BLOSUM62 word size 3
// threshold 13, BLOSUM45 word size 2 threshold 16), with the overflow array reached at each
// setting, and for a query set without any word; it declines word length 4 (20-bit backbone
// indexes), whose table the joint builder makes (`x_full_wide_table`; the test after this one).
#[test]
fn lean_table_equals_reference_for_thresholds_masks_matrices_and_word_lengths() {
    let mut seed = 11u64;
    let mut low_complexity = vec![1u8; 300];
    low_complexity.extend(residues(&mut seed, 40));
    let queries = vec![
        vec![
            frame(residues(&mut seed, 700), vec![(10, 30), (400, 401)], 1),
            frame(residues(&mut seed, 2), vec![], 2),
            frame(residues(&mut seed, 3), vec![], 3),
            frame(low_complexity, vec![(250, 299)], -1),
            frame(residues(&mut seed, 650), vec![(0, 0)], -2),
            frame(Vec::new(), vec![], -3),
        ],
        vec![frame(residues(&mut seed, 1200), vec![], 1)],
    ];
    // No word at word length 3: frames shorter than the word, and one masked whole.
    let wordless = vec![vec![
        frame(residues(&mut seed, 2), vec![], 1),
        frame(residues(&mut seed, 40), vec![(0, 39)], 2),
        frame(Vec::new(), vec![], 3),
    ]];
    let prepared_w = |queries: &[Vec<QueryFrame>], matrix: ScoringMatrix, word_length: usize| {
        let bounds = if matrix == ScoringMatrix::Blosum62 {
            (-4, 11)
        } else {
            matrix_score_bounds(matrix)
        };
        prepare_lookup_query(
            queries,
            ideal_karlin_params_for_matrix(matrix, bounds),
            &compute_std_aa_composition(),
            word_length,
            false,
            matrix,
            bounds,
        )
    };
    let same =
        |queries: &[Vec<QueryFrame>], matrix: ScoringMatrix, word_length: usize, threshold: i32| {
            let reference = build_lookup_from_prepared(
                prepared_w(queries, matrix, word_length),
                threshold,
                matrix,
                word_length,
            )
            .0;
            let lean = match x_build_lean(
                prepared_w(queries, matrix, word_length),
                threshold,
                matrix,
                word_length,
            ) {
                Ok(table) => table,
                Err(_) => panic!("the lean builder declined word length {word_length}"),
            };
            x_assert_same_table(&lean, &reference);
            // The lists made on a pool's threads (`x_build_lean_in`), in ranges of a few sources.
            #[cfg(all(
                feature = "parallel",
                any(not(target_arch = "wasm32"), feature = "wasm-threads")
            ))]
            {
                let pooled = crate::utils::threading::with_search_pool(4, "tblastn", |pool| {
                    Ok(x_build_lean_in(
                        prepared_w(queries, matrix, word_length),
                        threshold,
                        matrix,
                        word_length,
                        Some(pool),
                    ))
                })
                .unwrap();
                match pooled {
                    Ok(table) => x_assert_same_table(&table, &reference),
                    Err(_) => panic!("the lean builder declined word length {word_length}"),
                }
            }
            reference
        };
    for (matrix, word_length, thresholds) in [
        (
            ScoringMatrix::Blosum62,
            3,
            &[0, 1, 11, 12, 13, 19, 30, 1000][..],
        ),
        (ScoringMatrix::Blosum62, 2, &[0, 8, 11, 13][..]),
        (ScoringMatrix::Blosum45, 2, &[16][..]),
        (ScoringMatrix::Blosum45, 3, &[13, 14][..]),
        (ScoringMatrix::Pam30, 3, &[11][..]),
    ] {
        let mut overflow_seen = false;
        for &threshold in thresholds {
            overflow_seen |= !same(&queries, matrix, word_length, threshold)
                .overflow
                .is_empty();
        }
        assert!(
            overflow_seen,
            "{matrix:?} word length {word_length} must reach the overflow array"
        );
    }
    let empty = same(&wordless, ScoringMatrix::Blosum62, 3, 13);
    assert!(empty.pv.iter().all(|&bits| bits == 0) && empty.overflow.is_empty());
    assert!(x_build_lean(
        prepared_w(&queries, ScoringMatrix::Blosum62, 4),
        13,
        ScoringMatrix::Blosum62,
        4
    )
    .is_err());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:461-467 (module comment:
// `s_AddWordHits` for every word of the exact backbone, then `BlastAaLookupFinalize`)
// Item X2e of line RG: the tables TBLASTX's full-table build (`x_tbx_build_lookup`: the direct
// builder, 4-letter words the joint builder) and TBLASTN's (`x_lookup_table_for_profile_in`: the
// lean builder, 4-letter words the joint builder, also on a pool of 4 and 8 threads) give for
// another matrix than BLOSUM62 and every word length are the chain builder's, and TBLASTX's
// contexts too.
#[test]
fn tblastx_and_tblastn_tables_of_other_matrices_and_word_lengths_equal_the_chain_builders() {
    let mut seed = 41u64;
    let queries = vec![
        vec![
            frame(residues(&mut seed, 330), vec![(10, 40)], 1),
            frame(residues(&mut seed, 329), vec![], 2),
            frame(Vec::new(), vec![], 3),
            frame(residues(&mut seed, 330), vec![(0, 0)], -1),
        ],
        vec![frame(residues(&mut seed, 250), vec![], 1)],
    ];
    for matrix in [
        ScoringMatrix::Blosum62,
        ScoringMatrix::Blosum45,
        ScoringMatrix::Pam30,
    ] {
        let params = crate::stats::lookup_protein_params_ungapped(matrix);
        for word_length in [2, 3, 4] {
            for threshold in [0, 7, 11, 13] {
                let label = format!("{matrix:?} word length {word_length} threshold {threshold}");
                let (reference, ref_contexts) = build_ncbi_lookup_for_profile(
                    &queries,
                    threshold,
                    &params,
                    true,
                    matrix,
                    word_length,
                );
                let (got, contexts) =
                    x_tbx_build_lookup(&queries, threshold, &params, true, matrix, word_length);
                x_assert_same_table_as(&format!("TBLASTX {label}:"), &got, &reference);
                x_assert_same_contexts(&contexts, &ref_contexts);
                let reference = build_ncbi_lookup_for_profile(
                    &queries,
                    threshold,
                    &params,
                    false,
                    matrix,
                    word_length,
                )
                .0;
                let got =
                    x_lookup_table_for_profile_in(&queries, threshold, matrix, word_length, None);
                x_assert_same_table_as(&format!("TBLASTN {label}:"), &got, &reference);
                #[cfg(all(
                    feature = "parallel",
                    any(not(target_arch = "wasm32"), feature = "wasm-threads")
                ))]
                for threads in [4, 8] {
                    let pooled =
                        crate::utils::threading::with_search_pool(threads, "tblastn", |pool| {
                            Ok(x_lookup_table_for_profile_in(
                                &queries,
                                threshold,
                                matrix,
                                word_length,
                                Some(pool),
                            ))
                        })
                        .unwrap();
                    x_assert_same_table_as(
                        &format!("TBLASTN {label}, {threads} threads:"),
                        &pooled,
                        &reference,
                    );
                }
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:338-343
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:321-336
// The pooled build: a table built on vectors another table and another build left (cleared
// and zero-filled as calloc gives them) and with the neighbour lists shared by all threads is
// the reference table: after a smaller and a larger table, with two tables alive at once, and
// after the shared lists change key.
#[test]
fn pooled_tables_equal_the_reference_on_kept_vectors() {
    let mut seed = 23u64;
    let mut low_complexity = vec![1u8; 400];
    low_complexity.extend(residues(&mut seed, 60));
    let small = vec![vec![frame(residues(&mut seed, 300), vec![], 1)]];
    let mut invalid = residues(&mut seed, 700);
    for at in [0, 5, 6, 300, 699] {
        invalid[at] = 0x40;
    }
    let large = vec![
        vec![
            frame(residues(&mut seed, 2_500), vec![(100, 140), (143, 200)], 1),
            frame(low_complexity, vec![(350, 399)], -1),
        ],
        vec![frame(residues(&mut seed, 1_800), vec![], 2)],
        vec![frame(invalid, vec![], 3)],
    ];
    let pooled = |queries: &[Vec<QueryFrame>], threshold: i32, matrix: ScoringMatrix| {
        let reference =
            build_lookup_from_prepared(prepared(queries, matrix), threshold, matrix, 3).0;
        let table = match pooled_wide(prepared(queries, matrix), threshold, matrix, 1) {
            Ok(table) => table,
            Err(_) => panic!("in-place builder declined"),
        };
        x_assert_same_table(&table, &reference);
        table
    };
    for (queries, threshold) in [(&small, 12), (&large, 11), (&small, 13), (&large, 12)] {
        let mut table = pooled(queries, threshold, ScoringMatrix::Blosum62);
        x_give_back_pooled(&mut table);
        assert!(table.backbone.is_empty() && table.overflow.is_empty() && table.pv.is_empty());
    }
    let mut first = pooled(&large, 12, ScoringMatrix::Blosum62);
    let mut second = pooled(&small, 12, ScoringMatrix::Blosum62);
    x_give_back_pooled(&mut first);
    x_give_back_pooled(&mut second);
    for (queries, threshold, matrix) in [
        (&small, 12, ScoringMatrix::Blosum45),
        (&large, 11, ScoringMatrix::Pam30),
        (&large, 12, ScoringMatrix::Blosum62),
    ] {
        let mut table = pooled(queries, threshold, matrix);
        x_give_back_pooled(&mut table);
    }
}

// No NCBI counterpart (buffer reuse only): a pooled vector taken for fewer entries than it has
// room for is given back down to them (`n + n / 4` for the overflow), one short of room grows,
// a new one has the room of its first table, and a build takes the kept overflow vector that
// fits its table best: the smallest of at least its size, else the largest.
#[test]
fn pooled_vectors_fit_the_table_they_are_taken_for() {
    let vector = x_zeroed(vec![7i32; 4_000], 1_000, true);
    assert_eq!(vector.len(), 1_000);
    assert!(vector.iter().all(|&hit| hit == 0));
    assert!((1_000..=1_250).contains(&vector.capacity()));
    let vector = x_zeroed(vector, 3_000, true);
    assert!((3_000..=3_750).contains(&vector.capacity()));
    let vector = x_zeroed(vector, 500, false);
    assert_eq!(vector.len(), 500);
    assert!((500..1_000).contains(&vector.capacity()));
    let vector = x_zeroed(vector, 0, true);
    assert!(vector.is_empty() && (1..500).contains(&vector.capacity()));
    let vector = x_zeroed(Vec::<u32>::new(), 10, false);
    assert_eq!(vector.len(), 10);
    assert!((10..1_000).contains(&vector.capacity()));

    let with_room = |room: usize| Vec::<i32>::with_capacity(room);
    let mut kept: Vec<XAaVectors> = [100, 300, 200]
        .into_iter()
        .map(|room| XAaVectors {
            overflow: with_room(room),
            ..XAaVectors::default()
        })
        .collect();
    let rooms = |kept: &[XAaVectors]| -> Vec<usize> {
        kept.iter().map(|v| v.overflow.capacity()).collect()
    };
    let (r100, r300, r200) = (rooms(&kept)[0], rooms(&kept)[1], rooms(&kept)[2]);
    let mut own = with_room(50);
    let r50 = own.capacity();
    x_fitting_overflow::<BackboneCell>(&mut kept, &mut own, 150);
    assert_eq!(
        (own.capacity(), rooms(&kept)),
        (r200, vec![r100, r300, r50])
    );
    x_fitting_overflow::<BackboneCell>(&mut kept, &mut own, 1_000);
    assert_eq!(
        (own.capacity(), rooms(&kept)),
        (r300, vec![r100, r200, r50])
    );
    x_fitting_overflow::<BackboneCell>(&mut kept, &mut own, r300);
    assert_eq!(
        (own.capacity(), rooms(&kept)),
        (r300, vec![r100, r200, r50])
    );
}

// No NCBI counterpart (scheduling only): LOSAT_X_BXBUILDPAR's passes by source ranges give the
// reference table, with 2 to 8 ranges (the most; more asked for are clamped), with overflow
// cells (the low-complexity frame), a location shorter than a word (between the first frame's
// masks) and residues outside the alphabet's bits (the third query: the exact backbone skips
// their words), a query of fewer source cells than ranges and an all-masked one, thresholds
// with and without neighbours and other matrices, on vectors other builds left.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn partitioned_tables_equal_the_reference() {
    let mut seed = 31u64;
    let mut low_complexity = vec![1u8; 400];
    low_complexity.extend(residues(&mut seed, 60));
    let small = vec![vec![frame(residues(&mut seed, 300), vec![], 1)]];
    let mut invalid = residues(&mut seed, 700);
    for at in [0, 5, 6, 300, 699] {
        invalid[at] = 0x40;
    }
    let large = vec![
        vec![
            frame(residues(&mut seed, 2_500), vec![(100, 140), (143, 200)], 1),
            frame(low_complexity, vec![(350, 399)], -1),
        ],
        vec![frame(residues(&mut seed, 1_800), vec![], 2)],
        vec![frame(invalid, vec![], 3)],
    ];
    // Fewer source cells than parts, and none at all.
    let tiny = vec![vec![frame(residues(&mut seed, 6), vec![], 1)]];
    let masked = vec![vec![frame(residues(&mut seed, 50), vec![(0, 49)], 1)]];
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap();
    pool.install(|| {
        for (queries, threshold, matrix) in [
            (&large, 11, ScoringMatrix::Blosum62),
            (&small, 12, ScoringMatrix::Blosum62),
            (&large, 13, ScoringMatrix::Blosum45),
            (&large, 0, ScoringMatrix::Blosum62),
            (&small, 12, ScoringMatrix::Pam30),
            (&tiny, 11, ScoringMatrix::Blosum62),
            (&masked, 11, ScoringMatrix::Blosum62),
        ] {
            let reference =
                build_lookup_from_prepared(prepared(queries, matrix), threshold, matrix, 3).0;
            let edge = std::ptr::eq(queries, &tiny) || std::ptr::eq(queries, &masked);
            assert!(edge || threshold == 0 || !reference.overflow.is_empty());
            for parts in [2, 3, 4, 7, 8, 1_000] {
                let mut table =
                    match pooled_wide(prepared(queries, matrix), threshold, matrix, parts) {
                        Ok(table) => table,
                        Err(_) => panic!("in-place builder declined"),
                    };
                x_assert_same_table(&table, &reference);
                x_give_back_pooled(&mut table);
            }
        }
    });
}

#[test]
fn in_place_table_of_an_all_masked_query_is_empty() {
    let mut seed = 11u64;
    let queries = vec![vec![frame(residues(&mut seed, 50), vec![(0, 49)], 1)]];
    assert_eq!(same(&queries, 12, ScoringMatrix::Blosum62), 0);
    let table = match x_build_direct(
        prepared(&queries, ScoringMatrix::Blosum62),
        12,
        ScoringMatrix::Blosum62,
        3,
    ) {
        Ok(table) => table,
        Err(_) => panic!("in-place builder declined"),
    };
    assert!(table.pv.iter().all(|&w| w == 0));
    assert_eq!(table.longest_chain, 0);
}

/// BLASTP's builder against `build_ncbi_lookup`: every table field and every context field.
/// Returns the overflow length.
fn same_bp(queries: &[Vec<QueryFrame>], threshold: i32, check_ideal: bool) -> usize {
    let params = blosum62_ideal_karlin_params();
    let (reference, ref_contexts) = build_ncbi_lookup(queries, threshold, &params, check_ideal);
    let (got, contexts) = x_bp_build_lookup(queries, threshold, check_ideal);
    x_assert_same_table(&got, &reference);
    x_assert_same_contexts(&contexts, &ref_contexts);
    reference.overflow.len()
}

/// Every field of every context the same, the Karlin blocks by their printed values.
fn x_assert_same_contexts(contexts: &[QueryContext], ref_contexts: &[QueryContext]) {
    assert_eq!(contexts.len(), ref_contexts.len());
    for (a, b) in contexts.iter().zip(ref_contexts) {
        assert!(
            a.q_idx == b.q_idx
                && a.frame == b.frame
                && a.aa_seq == b.aa_seq
                && a.aa_seq_nomask == b.aa_seq_nomask
                && a.aa_len == b.aa_len
                && a.orig_len == b.orig_len
                && a.frame_base == b.frame_base
                && a.is_valid == b.is_valid
                && format!("{:?}", a.karlin_params) == format!("{:?}", b.karlin_params),
            "context {} differs",
            a.q_idx
        );
    }
}

/// One BLASTP query: a single frame (as blastp builds them), optionally masked.
fn bp_query(seq: Vec<u8>, seg_masks: Vec<(usize, usize)>, nomask: bool) -> Vec<QueryFrame> {
    let mut f = frame(seq, seg_masks, 1);
    if nomask {
        f.aa_seq_nomask = Some(f.aa_seq.clone());
    }
    vec![f]
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:461-467
#[test]
fn bplut_table_and_contexts_equal_build_ncbi_lookup() {
    let mut seed = 23u64;
    let mut low_complexity = vec![1u8; 200];
    low_complexity.extend(residues(&mut seed, 30));
    let queries = vec![
        bp_query(residues(&mut seed, 900), vec![(40, 90)], true),
        bp_query(residues(&mut seed, 2), vec![], false),
        bp_query(low_complexity, vec![(150, 199)], true),
        bp_query(residues(&mut seed, 1), vec![], false),
        bp_query(Vec::new(), vec![], false),
        bp_query(residues(&mut seed, 1500), vec![], false),
        bp_query(residues(&mut seed, 60), vec![(0, 59)], true),
    ];
    let mut overflow_seen = false;
    for check_ideal in [false, true] {
        for threshold in [0, 11, 12, 13, 21] {
            overflow_seen |= same_bp(&queries, threshold, check_ideal) > 0;
        }
    }
    assert!(overflow_seen, "the fixture must reach the overflow array");
    // A threshold beyond Int4 (`INT_MIN`, every word a neighbour) on short queries.
    let short: Vec<Vec<QueryFrame>> = (0..3)
        .map(|_| bp_query(residues(&mut seed, 15), vec![], false))
        .collect();
    same_bp(&short, i32::MIN, false);
}

/// The query offsets of every cell of `table`, as the scan reads them.
fn cell_lists(table: &BlastAaLookupTable) -> impl Iterator<Item = &[i32]> {
    table.backbone.iter().map(|cell| {
        let n = cell.num_used as usize;
        if n <= AA_HITS_PER_CELL {
            &cell.entries[..n]
        } else {
            let cursor = cell.cursor();
            &table.overflow[cursor..cursor + n]
        }
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:504-519
/// LOSAT_X_BPDIAG (`blastp/x_seed_fused.rs`) tests the hits of a cell in batches, which is
/// exact because a cell lists a query offset at most once (the word itself is added only when
/// the neighbour enumeration does not return it). Both BLASTP builders, thresholds without
/// neighbours (0), with few and with every word as a neighbour, masks and the overflow array.
#[test]
fn every_cell_lists_a_query_offset_once() {
    let mut seed = 31u64;
    let mut low_complexity = vec![1u8; 200];
    low_complexity.extend(residues(&mut seed, 30));
    let queries = vec![
        bp_query(residues(&mut seed, 900), vec![(40, 90)], true),
        bp_query(residues(&mut seed, 2), vec![], false),
        bp_query(low_complexity, vec![(150, 199)], true),
        bp_query(residues(&mut seed, 1500), vec![], false),
        bp_query(residues(&mut seed, 60), vec![(0, 59)], true),
    ];
    let short: Vec<Vec<QueryFrame>> = (0..3)
        .map(|_| bp_query(residues(&mut seed, 15), vec![], false))
        .collect();
    let params = blosum62_ideal_karlin_params();
    let mut long_cells = 0;
    for direct in [false, true] {
        for (set, threshold) in [
            (&queries, 0),
            (&queries, 11),
            (&queries, 13),
            (&short, 1),
            (&short, i32::MIN),
        ] {
            let (table, _) = if direct {
                x_bp_build_lookup(set, threshold, false)
            } else {
                build_ncbi_lookup(set, threshold, &params, false)
            };
            for hits in cell_lists(&table) {
                let mut sorted = hits.to_vec();
                sorted.sort_unstable();
                assert!(
                    sorted.windows(2).all(|pair| pair[0] != pair[1]),
                    "direct {direct} threshold {threshold}: a cell lists an offset twice: {hits:?}"
                );
                long_cells += usize::from(hits.len() > 16);
            }
        }
    }
    assert!(long_cells > 0, "the fixture must give long cells");
}

/// BLASTP queries whose table reaches the overflow array (as in
/// `bplut_table_and_contexts_equal_build_ncbi_lookup`).
fn bp_queries(seed: &mut u64) -> Vec<Vec<QueryFrame>> {
    let mut low_complexity = vec![1u8; 200];
    low_complexity.extend(residues(seed, 30));
    vec![
        bp_query(residues(seed, 900), vec![(40, 90)], true),
        bp_query(residues(seed, 2), vec![], false),
        bp_query(low_complexity, vec![(150, 199)], true),
        bp_query(Vec::new(), vec![], false),
        bp_query(residues(seed, 1500), vec![], false),
    ]
}

/// `x_bp_build_lookup(..)`'s table and contexts are those of `got`: every table field,
/// and every context field (the Karlin blocks by their bits).
fn assert_same_bp(got: &(BlastAaLookupTable, Vec<QueryContext>), queries: &[Vec<QueryFrame>]) {
    let (want, want_contexts) = x_bp_build_lookup(queries, 11, false);
    x_assert_same_table_as("ahead:", &got.0, &want);
    assert!(x_same_table(&got.0, &want));
    assert_eq!(got.1.len(), want_contexts.len());
    for (a, b) in got.1.iter().zip(&want_contexts) {
        assert!(
            a.q_idx == b.q_idx
                && a.aa_seq == b.aa_seq
                && a.frame_base == b.frame_base
                && a.is_valid == b.is_valid
                && format!("{:?}", a.karlin_params) == format!("{:?}", b.karlin_params),
            "context {} differs",
            a.q_idx
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:321-336
/// LOSAT_X_BPLKPAHEAD: the table built ahead on vectors kept by `XAaSpares` (a table given back
/// with other values, then moved to the build by `x_make_room`) is `x_bp_build_lookup(..)`'s;
/// `x_make_room` allocates only the vectors short of their room, and keeps the calling
/// thread's neighbour lists for the builds, which give them back.
#[test]
fn aa_spares_build_the_bplut_table_on_kept_vectors() {
    let mut seed = 41u64;
    let queries = bp_queries(&mut seed);
    let params = blosum62_ideal_karlin_params();
    let (first, _) = x_bp_build_lookup(&queries, 11, false);
    let (backbone, pv, overflow) = (first.backbone.len(), first.pv.len(), first.overflow.len());
    assert!(overflow > 0, "the fixture must reach the overflow array");
    assert!(X_MEMO.with(|memo| memo.borrow().is_some()));
    let spares = XAaSpares::default();
    // Fresh spares: all three allocated here, and this thread's lists move to the spares.
    assert_eq!(spares.x_make_room(backbone, pv, overflow, overflow), 3);
    assert!(X_MEMO.with(|memo| memo.borrow().is_none()));
    assert!(x_locked(&spares.memo).is_some());
    let addresses = |v: &XAaVectors| {
        [
            v.backbone.as_ptr() as usize,
            v.overflow.as_ptr() as usize,
            v.pv.as_ptr() as usize,
        ]
    };
    let table_addresses = |t: &BlastAaLookupTable| {
        [
            t.backbone.as_ptr() as usize,
            t.overflow.as_ptr() as usize,
            t.pv.as_ptr() as usize,
        ]
    };
    let kept = addresses(&x_locked(&spares.ahead));
    let built = spares
        .x_build_ahead(
            &queries,
            11,
            &params,
            false,
            ScoringMatrix::Blosum62,
            LOOKUP_WORD_LENGTH,
        )
        .expect("room for the table");
    // Built on the kept vectors (the build allocates none of them) with the kept lists, which
    // it gives back (this thread still has none).
    assert_eq!(table_addresses(&built.0), kept);
    assert!(X_MEMO.with(|memo| memo.borrow().is_none()));
    assert!(x_locked(&spares.memo).is_some(), "the lists are given back");
    assert!(
        x_locked(&spares.ahead).is_empty(),
        "the build took the vectors"
    );
    assert_same_bp(&built, &queries);
    // A table given back with other values: its vectors serve the next build.
    let (mut table, _) = built;
    for cell in table.backbone.iter_mut() {
        cell.num_used = 99;
        cell.entries = [7; AA_HITS_PER_CELL];
    }
    table.overflow.fill(-5);
    table.pv.fill(u32::MAX);
    let given = table_addresses(&table);
    spares.give(table);
    assert_eq!(spares.x_make_room(backbone, pv, overflow, overflow), 0);
    assert_eq!(addresses(&x_locked(&spares.ahead)), given);
    assert!(x_locked(&spares.block).is_empty());
    // Its own vectors are kept (those given back meanwhile are dropped); only a vector short of
    // its room is replaced.
    let (other_table, _) = x_bp_build_lookup(&queries, 11, false);
    spares.give(other_table);
    assert_eq!(spares.x_make_room(backbone, pv, overflow, overflow), 0);
    assert_eq!(addresses(&x_locked(&spares.ahead)), given);
    assert!(x_locked(&spares.block).is_empty());
    assert_eq!(
        spares.x_make_room(backbone, pv, overflow * 4, overflow * 8),
        1
    );
    assert!(x_locked(&spares.ahead).overflow.capacity() >= overflow * 8);
    let other = bp_queries(&mut seed);
    let built = spares
        .x_build_ahead(
            &other,
            11,
            &params,
            false,
            ScoringMatrix::Blosum62,
            LOOKUP_WORD_LENGTH,
        )
        .expect("room for the table");
    assert_same_bp(&built, &other);
    // Home: the lists go back to a thread without any (not over its own), and with `release`
    // the kept vectors are dropped.
    spares.give(built.0);
    X_MEMO.with(|memo| *memo.borrow_mut() = None);
    spares.x_home(false);
    assert!(X_MEMO.with(|memo| memo.borrow().is_some()));
    assert!(x_locked(&spares.memo).is_none());
    assert!(!x_locked(&spares.block).is_empty());
    spares.x_home(true);
    assert!(x_locked(&spares.block).is_empty() && x_locked(&spares.ahead).is_empty());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:321-336
/// LOSAT_X_BPLKREUSE: each block's own build on the block spare (`x_build_on_block_spare`) is
/// `x_bp_build_lookup(true, ..)`'s table: the first on new vectors; the next ones on the
/// vectors the block before gave back with other values, the same allocations where they have
/// the room (a table as large or smaller), grown in place for a larger one. The spare is empty
/// while a table built on it is in use.
#[test]
fn aa_block_spare_builds_each_block_on_the_vectors_given_back() {
    let mut seed = 47u64;
    let queries = bp_queries(&mut seed);
    let small = vec![bp_query(residues(&mut seed, 300), vec![], false)];
    let large = vec![
        bp_query(residues(&mut seed, 4000), vec![(100, 180)], true),
        bp_query(residues(&mut seed, 3000), vec![], false),
    ];
    let params = blosum62_ideal_karlin_params();
    let table_addresses = |t: &BlastAaLookupTable| {
        [
            t.backbone.as_ptr() as usize,
            t.overflow.as_ptr() as usize,
            t.pv.as_ptr() as usize,
        ]
    };
    let spare = XAaSpares::default();
    let first = spare.x_build_on_block_spare(
        &queries,
        11,
        &params,
        false,
        ScoringMatrix::Blosum62,
        LOOKUP_WORD_LENGTH,
    );
    assert!(
        !first.0.overflow.is_empty(),
        "the fixture must reach the overflow array"
    );
    assert!(x_locked(&spare.block).is_empty());
    assert_same_bp(&first, &queries);
    // Given back with other values: the next table is built on the same allocations.
    let (mut table, _) = first;
    let overflow = table.overflow.len();
    for cell in table.backbone.iter_mut() {
        cell.num_used = 99;
        cell.entries = [7; AA_HITS_PER_CELL];
    }
    table.overflow.fill(-5);
    table.pv.fill(u32::MAX);
    let given = table_addresses(&table);
    spare.give(table);
    let again = spare.x_build_on_block_spare(
        &queries,
        11,
        &params,
        false,
        ScoringMatrix::Blosum62,
        LOOKUP_WORD_LENGTH,
    );
    assert_eq!(table_addresses(&again.0), given);
    assert!(x_locked(&spare.block).is_empty());
    assert_same_bp(&again, &queries);
    spare.give(again.0);
    let smaller = spare.x_build_on_block_spare(
        &small,
        11,
        &params,
        false,
        ScoringMatrix::Blosum62,
        LOOKUP_WORD_LENGTH,
    );
    assert_eq!(table_addresses(&smaller.0), given);
    assert_same_bp(&smaller, &small);
    // A larger table: the overflow vector grows to its size.
    spare.give(smaller.0);
    let grown = spare.x_build_on_block_spare(
        &large,
        11,
        &params,
        false,
        ScoringMatrix::Blosum62,
        LOOKUP_WORD_LENGTH,
    );
    assert!(grown.0.overflow.len() > overflow);
    assert_same_bp(&grown, &large);
}

/// LOSAT_X_BPLKPAHEAD: vectors without the room for the table give no table and are left as
/// they were; the neighbour lists are kept.
#[test]
fn aa_spares_without_the_room_give_no_table() {
    let mut seed = 43u64;
    let queries = bp_queries(&mut seed);
    let params = blosum62_ideal_karlin_params();
    let (table, _) = x_bp_build_lookup(&queries, 11, false);
    let (backbone, pv, overflow) = (table.backbone.len(), table.pv.len(), table.overflow.len());
    for (b, p, o) in [
        (backbone, pv, overflow - 1),
        (backbone - 1, pv, overflow),
        (backbone, pv - 1, overflow),
    ] {
        let spares = XAaSpares::default();
        spares.x_make_room(b, p, o, o);
        let capacities = |spares: &XAaSpares| {
            let ahead = x_locked(&spares.ahead);
            (
                ahead.backbone.capacity(),
                ahead.pv.capacity(),
                ahead.overflow.capacity(),
            )
        };
        let before = capacities(&spares);
        assert!(
            before.0 < backbone || before.1 < pv || before.2 < overflow,
            "{before:?}"
        );
        assert!(spares
            .x_build_ahead(
                &queries,
                11,
                &params,
                false,
                ScoringMatrix::Blosum62,
                LOOKUP_WORD_LENGTH
            )
            .is_none());
        assert_eq!(capacities(&spares), before);
        assert!(x_locked(&spares.memo).is_some());
    }
}

/// Item SC3: blastp's tables of word size 2 (blastp-short and `-word_size 2`). The in-place
/// build in both widths, the build on a block spare (given back and built again) and the
/// builds ahead hold the reference builder's table of the words of 2, with the same contexts,
/// for BLOSUM62 and PAM30 at thresholds that reach the overflow array.
#[test]
fn word_size_two_tables_are_the_reference_builders() {
    let mut seed = 59u64;
    let queries = bp_queries(&mut seed);
    let params = blosum62_ideal_karlin_params();
    for (matrix, threshold) in [
        (ScoringMatrix::Blosum62, 11),
        (ScoringMatrix::Blosum62, 6),
        (ScoringMatrix::Pam30, 16),
        (ScoringMatrix::Pam30, 9),
    ] {
        let name = format!("{matrix:?} threshold {threshold}:");
        let (want, want_contexts) = build_lookup_from_prepared(
            x_bp_prepared_query(&queries, false, matrix, 2),
            threshold,
            matrix,
            2,
        );
        assert_eq!(want.word_length, 2, "{name}");
        assert!(
            want.longest_chain as usize > AA_HITS_PER_CELL && !want.overflow.is_empty(),
            "{name} the fixture must reach the overflow array"
        );
        let same_contexts = |got: &[QueryContext]| {
            assert_eq!(got.len(), want_contexts.len(), "{name}");
            for (a, b) in got.iter().zip(&want_contexts) {
                assert!(
                    a.q_idx == b.q_idx
                        && a.aa_seq == b.aa_seq
                        && a.frame_base == b.frame_base
                        && a.is_valid == b.is_valid
                        && format!("{:?}", a.karlin_params) == format!("{:?}", b.karlin_params),
                    "{name} context {} differs",
                    a.q_idx
                );
            }
        };
        let wide_of = |built: XBpLookup| match built {
            XBpLookup::Wide(table) => table,
            XBpLookup::Narrow(_) => panic!("{name} a 32-bit table was asked for"),
        };
        let narrow_of = |built: XBpLookup| match built {
            XBpLookup::Narrow(table) => table,
            XBpLookup::Wide(_) => panic!("{name} the block fits 16-bit hits"),
        };
        let (built, contexts) = x_bp_build_lookup_as(&queries, threshold, false, false, matrix, 2);
        x_assert_same_table_as(&name, &wide_of(built), &want);
        same_contexts(&contexts);
        let (built, contexts) = x_bp_build_lookup_as(&queries, threshold, false, true, matrix, 2);
        x_assert_narrow_same_as_wide(&name, &narrow_of(built), &want);
        same_contexts(&contexts);
        // On a block spare: the vectors given back from a table of words of 3, then from the
        // table just built.
        let spare = XAaSpares::default();
        spare.give(x_bp_build_lookup(&queries, 11, false).0);
        for _ in 0..2 {
            let (table, contexts) =
                spare.x_build_on_block_spare(&queries, threshold, &params, false, matrix, 2);
            x_assert_same_table_as(&name, &table, &want);
            same_contexts(&contexts);
            spare.give(table);
        }
        let (built, contexts) =
            spare.x_build_on_block_spare16(&queries, threshold, false, matrix, 2);
        let narrow = narrow_of(built);
        x_assert_narrow_same_as_wide(&name, &narrow, &want);
        same_contexts(&contexts);
        // On the 16-bit vectors given back from a table of words of 3.
        let three = narrow_of(
            x_bp_build_lookup_as(
                &queries,
                11,
                false,
                true,
                ScoringMatrix::Blosum62,
                LOOKUP_WORD_LENGTH,
            )
            .0,
        );
        spare.give16(three);
        let (built, contexts) =
            spare.x_build_on_block_spare16(&queries, threshold, false, matrix, 2);
        x_assert_narrow_same_as_wide(&name, &narrow_of(built), &want);
        same_contexts(&contexts);
        // Ahead, on the room made for it.
        let spares = XAaSpares::default();
        spares.x_make_room(want.backbone.len(), want.pv.len(), want.overflow.len(), 0);
        let (table, contexts) = spares
            .x_build_ahead(&queries, threshold, &params, false, matrix, 2)
            .expect("the room made holds the table");
        x_assert_same_table_as(&name, &table, &want);
        same_contexts(&contexts);
        spares.x_make_room16(
            narrow.backbone.len(),
            narrow.pv.len(),
            narrow.overflow.len(),
            0,
        );
        let (built, contexts) = spares
            .x_build_ahead16(&queries, threshold, false, matrix, 2)
            .expect("the room made holds the table");
        x_assert_narrow_same_as_wide(&name, &narrow_of(built), &want);
        same_contexts(&contexts);
        spares.x_home(true);
    }
}

/// The comparison finds one backbone cell, one overflow entry, one pv bit and one header field
/// that differ, and names the first difference after its prefix.
#[test]
#[cfg_attr(panic = "abort", ignore = "catches a panic: the build aborts on panic")]
fn same_table_finds_one_cell_entry_bit_or_field() {
    let mut seed = 47u64;
    let queries = bp_queries(&mut seed);
    let build = || x_bp_build_lookup(&queries, 11, false).0;
    let want = build();
    assert!(x_same_table(&build(), &want));
    let cell = want
        .backbone
        .iter()
        .position(|cell| cell.num_used > 0 && cell.num_used as usize <= AA_HITS_PER_CELL)
        .expect("a cell with entries");
    let changes: [(&str, &dyn Fn(&mut BlastAaLookupTable)); 4] = [
        ("backbone cell", &|t| t.backbone[cell].entries[0] += 1),
        ("overflow array", &|t| t.overflow[0] += 1),
        ("pv", &|t| t.pv[0] ^= 1 << 31),
        ("table header", &|t| t.longest_chain += 1),
    ];
    for (what, change) in changes {
        let mut got = build();
        change(&mut got);
        assert!(!x_same_table(&got, &want), "{what}");
        let message = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            x_assert_same_table_as("PREFIX:", &got, &want)
        }))
        .expect_err(what);
        let message = message
            .downcast_ref::<String>()
            .cloned()
            .or_else(|| message.downcast_ref::<&str>().map(|s| s.to_string()))
            .unwrap_or_default();
        assert!(message.starts_with(&format!("PREFIX: {what}")), "{message}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:283-358 (cited at
// `x_build_direct`)
// LOSAT_X_BPLK16: the table of 16-bit hits holds, cell by cell, the hits of the table of
// 32-bit hits of the same block in the same order (the comparison of LOSAT_X_BPLK16SHADOW,
// against the reference builder), for short queries, a low-complexity query whose cells are far
// above AA_HITS_PER_CELL (the overflow array, with its cursor split over two entries), several
// thresholds, and the block that is just at the 65,535 letters of the choice (16-bit) and one
// letter above it (32-bit).
#[test]
fn table_of_u16_hits_equals_the_table_of_i32_hits() {
    let mut seed = 16u64;
    let mut low_complexity = vec![1u8; 300];
    low_complexity.extend(residues(&mut seed, 40));
    let queries = vec![
        vec![
            frame(residues(&mut seed, 700), vec![(10, 30), (400, 401)], 1),
            frame(residues(&mut seed, 2), vec![], 2),
            frame(low_complexity, vec![(250, 299)], -1),
            frame(Vec::new(), vec![], -3),
        ],
        vec![frame(residues(&mut seed, 1200), vec![], 1)],
    ];
    let mut overflow_seen = false;
    for threshold in [0, 11, 13, 1000] {
        let (built, contexts) = x_bp_build_lookup_as(
            &queries,
            threshold,
            false,
            true,
            ScoringMatrix::Blosum62,
            LOOKUP_WORD_LENGTH,
        );
        let XBpLookup::Narrow(narrow) = built else {
            panic!("a small block takes 16-bit hits (threshold {threshold})");
        };
        let (wide, wide_contexts) = x_bp_build_lookup(&queries, threshold, false);
        assert_eq!(contexts.len(), wide_contexts.len());
        x_assert_narrow_same_as_wide(&format!("threshold {threshold}:"), &narrow, &wide);
        overflow_seen |= !wide.overflow.is_empty();
        assert!(wide.longest_chain as usize > AA_HITS_PER_CELL);
    }
    assert!(overflow_seen, "the fixture must reach the overflow array");
    // Not asked for: the table of 32-bit hits.
    assert!(matches!(
        x_bp_build_lookup_as(
            &queries,
            11,
            false,
            false,
            ScoringMatrix::Blosum62,
            LOOKUP_WORD_LENGTH
        )
        .0,
        XBpLookup::Wide(_)
    ));
    // The choice at 65,535 letters: the lookup query's length is the last context's frame_base
    // plus its aa_len.
    let base = x_bp_prepared_query(
        &[vec![frame(residues(&mut seed, 50), vec![], 1)]],
        false,
        ScoringMatrix::Blosum62,
        LOOKUP_WORD_LENGTH,
    )
    .contexts[0]
        .frame_base as usize;
    for (length, narrow) in [(65_535 - base, true), (65_536 - base, false)] {
        let block = vec![vec![frame(residues(&mut seed, length), vec![], 1)]];
        assert_eq!(
            x_bplk16_fits(&x_bp_prepared_query(
                &block,
                false,
                ScoringMatrix::Blosum62,
                LOOKUP_WORD_LENGTH
            )),
            narrow,
            "{length} letters"
        );
        match x_bp_build_lookup_as(
            &block,
            11,
            false,
            true,
            ScoringMatrix::Blosum62,
            LOOKUP_WORD_LENGTH,
        )
        .0
        {
            XBpLookup::Narrow(table) => {
                assert!(narrow, "{length} letters took 16-bit hits");
                assert_eq!(table.query_length, 65_535);
                let reference = build_lookup_from_prepared(
                    x_bp_prepared_query(&block, false, ScoringMatrix::Blosum62, LOOKUP_WORD_LENGTH),
                    11,
                    ScoringMatrix::Blosum62,
                    LOOKUP_WORD_LENGTH,
                )
                .0;
                x_assert_narrow_same_as_wide("65,535 letters:", &table, &reference);
            }
            XBpLookup::Wide(table) => {
                assert!(!narrow, "{length} letters took 32-bit hits");
                assert_eq!(table.query_length, 65_536);
            }
        }
    }
}

// LOSAT_X_BPLK16: the vectors of a table of 16-bit hits given back (`give16`) serve the next
// block's build on the same spares (`x_build_on_block_spare16`, grown in place where short),
// and the vectors made ready for a build ahead (`x_make_room` with the narrow width, the
// padding counted in the sizes the table before reports) serve `x_build_ahead16`; a spare of the
// other width is not reused (the room is made for the width asked for) and the tables are the
// tables of 32-bit hits.
#[test]
fn spares_serve_the_builds_of_both_widths() {
    let mut seed = 160u64;
    let mut low_complexity = vec![1u8; 200];
    low_complexity.extend(residues(&mut seed, 30));
    let blocks: Vec<Vec<Vec<QueryFrame>>> = vec![
        vec![vec![frame(residues(&mut seed, 900), vec![], 1)]],
        vec![
            vec![frame(low_complexity, vec![], 1)],
            vec![frame(residues(&mut seed, 1500), vec![(7, 70)], 1)],
        ],
        vec![vec![frame(residues(&mut seed, 300), vec![], 1)]],
    ];
    let params = blosum62_ideal_karlin_params();
    let spares = XAaSpares::default();
    let narrow_of = |built: XBpLookup| match built {
        XBpLookup::Narrow(table) => table,
        XBpLookup::Wide(_) => panic!("a small block takes 16-bit hits"),
    };
    // The builds on the block spare, each on the vectors the block before gave back.
    for block in &blocks {
        let table = narrow_of(
            spares
                .x_build_on_block_spare16(
                    block,
                    11,
                    false,
                    ScoringMatrix::Blosum62,
                    LOOKUP_WORD_LENGTH,
                )
                .0,
        );
        let wide = x_bp_build_lookup(block, 11, false).0;
        x_assert_narrow_same_as_wide("on the block spare:", &table, &wide);
        spares.give16(table);
    }
    // The build ahead on the room made for the next table, from the table before it.
    let first = narrow_of(
        spares
            .x_build_on_block_spare16(
                &blocks[0],
                11,
                false,
                ScoringMatrix::Blosum62,
                LOOKUP_WORD_LENGTH,
            )
            .0,
    );
    spares.x_make_room16(first.backbone.len(), first.pv.len(), 1 << 16, 1 << 17);
    spares.give16(first);
    for block in &blocks {
        let (built, _) = spares
            .x_build_ahead16(
                block,
                11,
                false,
                ScoringMatrix::Blosum62,
                LOOKUP_WORD_LENGTH,
            )
            .expect("the room made holds the table");
        let table = narrow_of(built);
        let wide = x_bp_build_lookup(block, 11, false).0;
        x_assert_narrow_same_as_wide("built ahead:", &table, &wide);
        spares.x_make_room16(table.backbone.len(), table.pv.len(), 1 << 16, 1 << 17);
    }
    // A table of 32-bit hits built ahead takes room of its own width, not the 16-bit vectors.
    assert!(spares
        .x_build_ahead(
            &blocks[0],
            11,
            &params,
            false,
            ScoringMatrix::Blosum62,
            LOOKUP_WORD_LENGTH
        )
        .is_none());
    let wide = x_bp_build_lookup(&blocks[0], 11, false).0;
    spares.x_make_room(wide.backbone.len(), wide.pv.len(), 1 << 16, 1 << 17);
    let (table, _) = spares
        .x_build_ahead(
            &blocks[0],
            11,
            &params,
            false,
            ScoringMatrix::Blosum62,
            LOOKUP_WORD_LENGTH,
        )
        .expect("the room made holds the table");
    x_assert_same_table(&table, &wide);
}

// The tests below give narrow tables back to the pool, which is shared by the tests of this
// module: they run one at a time, so that a test finds the narrow sets it left there (no other
// test builds a narrow table).
static NARROW_POOL: std::sync::Mutex<()> = std::sync::Mutex::new(());

/// The parts of the partitioned build the tests try (0: not pooled; 1: pooled, serial).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
const PARTS: &[usize] = &[0, 1, 2, 3, 4, 7, 8, 1_000];
#[cfg(not(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
)))]
const PARTS: &[usize] = &[0, 1];

/// `f` on a pool of four threads (the partitioned build runs its passes on the pool it is
/// called on), or on this thread where the build has no pool.
fn in_pool(f: impl FnOnce() + Send) {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap()
        .install(f);
    #[cfg(not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )))]
    f();
}

/// The reference builder's table for `queries`.
fn reference(
    queries: &[Vec<QueryFrame>],
    threshold: i32,
    matrix: ScoringMatrix,
) -> BlastAaLookupTable {
    build_lookup_from_prepared(prepared(queries, matrix), threshold, matrix, 3).0
}

/// `table` equals the reference table of the same query, whatever its width.
fn assert_bx_table(what: &str, table: &XBxTable, reference: &BlastAaLookupTable) {
    match table {
        XBxTable::Wide(table) => x_assert_same_table(table, reference),
        XBxTable::Narrow(table) => x_assert_narrow_same_as_wide(what, table, reference),
    }
}

/// Gives `table`'s vectors back to the pool.
fn give_back(table: &mut XBxTable) {
    match table {
        XBxTable::Wide(table) => x_give_back_pooled(table),
        XBxTable::Narrow(table) => x_give_back_pooled16(table),
    }
}

/// The hits and used cells of `table`, and so its hits per used cell in hundredths.
fn hundredths(table: &BlastAaLookupTable) -> u64 {
    let (hits, used) = BackboneCell::x_hits_and_used(&table.backbone);
    x_hits_per_cell_hundredths(hits, used)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:283-358 (cited at
// `x_build_direct`)
// BLASTX's table with 16-bit hits: the serial build (not pooled, and pooled) and the partitioned
// build (2-8 parts, and more parts than cells) give, cell by cell, the table of the reference
// builder (count, hit list in order, pv, `longest_chain`, the zero padding after the cells and
// the overflow), for tables with cells in the overflow array (its cursor over two entries) and
// several thresholds; the width is chosen from the table's hits per used cell, on both sides of
// the cutoff (`Choose(cutoff)` at the table's own ratio and one hundredth above it), and from
// the lookup query's length, at 65,535 letters (16 bits) and 65,536 (32 bits).
#[test]
fn bx_tables_of_u16_hits_equal_the_reference_tables() {
    let _one_at_a_time = NARROW_POOL.lock().unwrap_or_else(|e| e.into_inner());
    let mut seed = 18u64;
    let mut low_complexity = vec![1u8; 400];
    low_complexity.extend(residues(&mut seed, 60));
    let small = vec![vec![frame(residues(&mut seed, 300), vec![], 1)]];
    let large = vec![
        vec![
            frame(residues(&mut seed, 2_500), vec![(100, 140), (143, 200)], 1),
            frame(low_complexity, vec![(350, 399)], -1),
        ],
        vec![frame(residues(&mut seed, 1_800), vec![], 2)],
    ];
    let tiny = vec![vec![frame(residues(&mut seed, 6), vec![], 1)]];
    let masked = vec![vec![frame(residues(&mut seed, 50), vec![(0, 49)], 1)]];
    let build_with = |queries: &[Vec<QueryFrame>],
                      threshold: i32,
                      parts: usize,
                      cutoff,
                      matrix: ScoringMatrix| {
        let built = if parts == 0 {
            x_bx_build_direct(prepared(queries, matrix), threshold, matrix, 3, cutoff)
        } else {
            x_build_direct_pooled(
                prepared(queries, matrix),
                threshold,
                matrix,
                3,
                parts,
                cutoff,
            )
        };
        match built {
            Ok(table) => table,
            Err(_) => panic!("in-place builder declined"),
        }
    };
    let mut overflow_seen = false;
    let mut narrow_built = 0;
    // BLOSUM62's tables, then other matrices' (item SC6b of line NP: BLASTX builds them the same
    // way), each after the pooled builds of the matrix before it (the shared neighbour lists are
    // kept per matrix).
    for matrix in [
        ScoringMatrix::Blosum62,
        ScoringMatrix::Blosum45,
        ScoringMatrix::Pam30,
    ] {
        in_pool(|| {
            for (queries, threshold) in [
                (&large, 11),
                (&small, 12),
                (&large, 13),
                (&large, 0),
                (&tiny, 11),
                (&masked, 11),
            ] {
                let reference = reference(queries, threshold, matrix);
                overflow_seen |= !reference.overflow.is_empty();
                let ratio = hundredths(&reference);
                // `parts` 0: not pooled; 1: pooled, serial; above: partitioned.
                for &parts in PARTS {
                    // A cutoff of 0 chooses 16 bits wherever the query fits; none, 32 bits.
                    for (cutoff, narrow) in [
                        (Some(0), true),
                        (Some(ratio), true),
                        (Some(ratio + 1), false),
                        (None, false),
                    ] {
                        let mut table = build_with(queries, threshold, parts, cutoff, matrix);
                        let what = format!(
                            "{matrix:?} threshold {threshold} parts {parts} cutoff {cutoff:?}:"
                        );
                        assert_eq!(
                            matches!(table, XBxTable::Narrow(_)),
                            narrow,
                            "{what} the width"
                        );
                        narrow_built += usize::from(narrow);
                        assert_bx_table(&what, &table, &reference);
                        if parts > 0 {
                            give_back(&mut table);
                        }
                    }
                }
            }
        });
    }
    assert!(overflow_seen, "the fixture must reach the overflow array");
    assert!(narrow_built > 50, "{narrow_built} tables of 16-bit hits");
    let matrix = ScoringMatrix::Blosum62;
    // The width at 65,535 letters of lookup query: its length is the last context's frame_base
    // plus its aa_len.
    let base = prepared(&[vec![frame(residues(&mut seed, 50), vec![], 1)]], matrix).contexts[0]
        .frame_base as usize;
    for (length, narrow) in [(65_535 - base, true), (65_536 - base, false)] {
        let block = vec![vec![frame(residues(&mut seed, length), vec![], 1)]];
        let reference = reference(&block, 13, matrix);
        assert_eq!(reference.query_length, if narrow { 65_535 } else { 65_536 });
        in_pool(|| {
            for parts in [0, 1, 3] {
                let mut table = build_with(&block, 13, parts, Some(0), matrix);
                let what = format!("{length} letters, parts {parts}:");
                assert_eq!(matches!(table, XBxTable::Narrow(_)), narrow, "{what}");
                assert_bx_table(&what, &table, &reference);
                if parts > 0 {
                    give_back(&mut table);
                }
            }
        });
    }
}

// No NCBI counterpart (buffer reuse only): a pooled build takes the kept set of its own width,
// the last given back, and leaves the sets of the other width where they are (a free of one of
// their mapped vectors would raise glibc's mmap threshold, `x_zeroed`).
#[test]
fn kept_sets_are_taken_by_width() {
    let set = |wide: usize, narrow: usize| XAaVectors {
        backbone: Vec::with_capacity(wide),
        backbone16: Vec::with_capacity(narrow),
        pv: Vec::with_capacity(1),
        ..XAaVectors::default()
    };
    let mut pool = vec![set(3, 0), set(0, 5), set(2, 0), set(0, 7)];
    let taken = x_take_kept_from::<BackboneCell>(&mut pool);
    assert_eq!(taken.backbone.capacity(), 2);
    let taken = x_take_kept_from::<BackboneCell16>(&mut pool);
    assert_eq!(taken.backbone16.capacity(), 7);
    let taken = x_take_kept_from::<BackboneCell16>(&mut pool);
    assert_eq!(taken.backbone16.capacity(), 5);
    // No narrow set left: an empty one, and the wide set stays.
    let taken = x_take_kept_from::<BackboneCell16>(&mut pool);
    assert!(taken.is_empty());
    assert_eq!(pool.len(), 1);
    assert_eq!(pool[0].backbone.capacity(), 3);
    let taken = x_take_kept_from::<BackboneCell>(&mut pool);
    assert_eq!(taken.backbone.capacity(), 3);
    assert!(pool.is_empty());
}

// No NCBI counterpart (buffer reuse only): the pool serves builds of both widths in turn: a
// table of 16-bit hits given back, builds of 32-bit hits, then a build of 16 bits on the kept
// vectors, each table the reference's. The pool is the process's, which other tests' searches
// may use at the same time, so the sets it holds are not counted here
// (`kept_sets_are_taken_by_width` takes sets by width from a pool of its own).
#[test]
fn pools_serve_the_builds_of_both_widths() {
    let _one_at_a_time = NARROW_POOL.lock().unwrap_or_else(|e| e.into_inner());
    let mut seed = 19u64;
    let queries = vec![vec![frame(residues(&mut seed, 1_500), vec![], 1)]];
    let matrix = ScoringMatrix::Blosum62;
    let reference = reference(&queries, 11, matrix);
    let build = |cutoff| match x_build_direct_pooled(
        prepared(&queries, matrix),
        11,
        matrix,
        3,
        1,
        cutoff,
    ) {
        Ok(table) => table,
        Err(_) => panic!("in-place builder declined"),
    };
    let XBxTable::Narrow(mut first) = build(Some(0)) else {
        panic!("a small query takes 16-bit hits");
    };
    x_assert_narrow_same_as_wide("first:", &first, &reference);
    x_give_back_pooled16(&mut first);
    // Builds of 32-bit hits, which give theirs back too.
    for _ in 0..3 {
        let mut wide = build(None);
        assert_bx_table("wide:", &wide, &reference);
        give_back(&mut wide);
    }
    // A build of 16 bits on a kept set (the one given back above, unless another test's
    // search took it).
    let XBxTable::Narrow(mut second) = build(Some(0)) else {
        panic!("a small query takes 16-bit hits");
    };
    x_assert_narrow_same_as_wide("second:", &second, &reference);
    x_give_back_pooled16(&mut second);
}
