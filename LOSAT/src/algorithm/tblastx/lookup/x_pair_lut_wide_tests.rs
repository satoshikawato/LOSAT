use super::super::super::tests::{frame, residues};
use super::*;

const FRAMES: [i8; 6] = [1, 2, 3, -1, -2, -3];

/// The table of the query offsets `bounds[g] <= q < bounds[g + 1]` filtered out of the full
/// table, order kept, in `BlastAaLookupFinalize`'s layout (the rule of `x_pair_par.rs`
/// `build_groups`), with the bounds of its entries.
fn filtered(full: &BlastAaLookupTable, bounds: &[i32], g: usize) -> XPairLutGroup {
    let lo = bounds[g];
    let hi = bounds.get(g + 1).copied().unwrap_or(i32::MAX);
    let mut table = BlastAaLookupTable {
        backbone: vec![BackboneCell::default(); full.backbone.len()],
        overflow: Vec::new(),
        pv: vec![0u32; full.pv.len()],
        frame_bases: full.frame_bases.clone(),
        num_contexts: full.num_contexts,
        query_length: full.query_length,
        word_length: full.word_length,
        alphabet_size: full.alphabet_size,
        charsize: full.charsize,
        mask: full.mask,
        longest_chain: 0,
        threshold: full.threshold,
        row_max: full.row_max.clone(),
    };
    let (mut q_min, mut q_max) = (u32::MAX, 0u32);
    for idx in 0..full.backbone.len() {
        let chain: Vec<i32> = full
            .get_hits(idx)
            .iter()
            .copied()
            .filter(|&q| lo <= q && q < hi)
            .collect();
        if chain.is_empty() {
            continue;
        }
        pv_set(&mut table.pv, idx);
        table.backbone[idx].num_used = chain.len() as i32;
        if chain.len() <= AA_HITS_PER_CELL {
            table.backbone[idx].entries[..chain.len()].copy_from_slice(&chain);
        } else {
            table.backbone[idx].set_cursor(table.overflow.len());
            table.overflow.extend_from_slice(&chain);
        }
        table.longest_chain = table.longest_chain.max(chain.len() as i32);
        for &q in &chain {
            q_min = q_min.min(q as u32);
            q_max = q_max.max(q as u32);
        }
    }
    XPairLutGroup {
        q_min: if table.longest_chain > 0 { q_min } else { 0 },
        q_max,
        table,
        rank: Vec::new(),
    }
}

/// Asserts that `wide` equals `reference` (a table of the full backbone), compacting the
/// reference where the wide table is compacted.
fn assert_same(g: usize, wide: &XPairLutGroup, mut reference: XPairLutGroup) {
    if !wide.rank.is_empty() {
        let rank = crate::algorithm::tblastx::x_scan_bf::x_compact(&mut reference.table);
        assert_eq!(wide.rank, rank, "table {g}: rank differs");
    } else {
        assert!(
            !crate::algorithm::tblastx::x_scan_bf::x_scan_bf_choice(&reference.table),
            "table {g}: the wide table is not compacted where the scan takes it"
        );
    }
    x_assert_same_group(
        g,
        &wide.table,
        (wide.q_min, wide.q_max),
        &reference.table,
        (reference.q_min, reference.q_max),
    );
}

/// The counts of `check`: the full table's entries, the tables with entries, the entries in
/// the tables' overflow arrays and the compacted tables.
#[derive(Default, Debug, Clone, Copy, PartialEq, Eq)]
struct Seen {
    entries: usize,
    with_entries: usize,
    in_overflow: usize,
    compacted: usize,
}

/// Builds the wide tables of `queries` (`matrix`, `word`-letter words, `threshold`) for the
/// group starts `starts` (None: one table per context) on `pool`, and compares each, field by
/// field, with the filter of the chain builder's full table (`build_lookup_from_prepared`) and,
/// for words of 2 and 3 letters, with `x_build_group`'s table; with one group, the table is the
/// full table itself. The stub (`x_into_stub`) carries the full table's header and
/// `longest_chain`.
fn check(
    queries: &[Vec<QueryFrame>],
    threshold: i32,
    matrix: ScoringMatrix,
    word: usize,
    starts: Option<&[usize]>,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> Seen {
    let params = crate::stats::lookup_protein_params_ungapped(matrix);
    let (full, full_contexts) = crate::algorithm::tblastx::lookup::build_ncbi_lookup_for_profile(
        queries, threshold, &params, true, matrix, word,
    );
    assert_eq!(full.word_length, word as i32);
    let q = XPairLutQuery::x_prepare_for(queries, threshold, true, matrix, word)
        .expect("a builder covers the word length");
    x_assert_same_contexts(q.contexts(), &full_contexts);
    let starts: Vec<usize> = match starts {
        Some(starts) => starts.to_vec(),
        None => (0..q.contexts().len()).collect(),
    };
    let bounds: Vec<i32> = starts.iter().map(|&c| full.frame_bases[c]).collect();
    let locations = q
        .x_group_locations(&starts)
        .expect("every location inside one group");
    let wide = q.x_build_tables_wide(&locations, pool, true);
    assert_eq!(wide.len(), starts.len());
    // Stripes of 32 cells (as many as the backbone holds) and of 2^16 cells (one for words of up
    // to 3 letters, 14 for 4): the same tables.
    let stripe_bits: &[u32] = &[5, 16];
    for &bits in stripe_bits {
        let striped =
            q.x_build_tables_wide_with(&locations, XWideOn::Pool(pool), true, Some(bits), true);
        for (g, (a, b)) in wide.iter().zip(&striped).enumerate() {
            assert_eq!(
                a.rank, b.rank,
                "table {g}: rank differs with stripes of 2^{bits} cells"
            );
            x_assert_same_group(
                g,
                &b.table,
                (b.q_min, b.q_max),
                &a.table,
                (a.q_min, a.q_max),
            );
        }
    }
    let mut seen = Seen {
        entries: full.backbone.iter().map(|c| c.num_used as usize).sum(),
        ..Seen::default()
    };
    for (g, table) in wide.iter().enumerate() {
        assert_same(g, table, filtered(&full, &bounds, g));
        if x_direct_covers(word) {
            let direct = q.x_build_group(&locations[g], true);
            assert_eq!(
                direct.rank, table.rank,
                "table {g}: rank differs from x_build_group"
            );
            x_assert_same_group(
                g,
                &table.table,
                (table.q_min, table.q_max),
                &direct.table,
                (direct.q_min, direct.q_max),
            );
        }
        seen.with_entries += usize::from(table.table.longest_chain > 0);
        seen.in_overflow += table.table.overflow.len();
        seen.compacted += usize::from(!table.rank.is_empty());
    }
    if starts == [0] {
        // One group: the full table itself.
        let entries: Vec<i32> = (0..full.backbone.len())
            .flat_map(|idx| full.get_hits(idx).iter().copied())
            .collect();
        let reference = XPairLutGroup {
            q_min: entries.iter().map(|&q| q as u32).min().unwrap_or(0),
            q_max: entries.iter().map(|&q| q as u32).max().unwrap_or(0),
            rank: Vec::new(),
            table: BlastAaLookupTable {
                backbone: full.backbone.clone(),
                overflow: full.overflow.clone(),
                pv: full.pv.clone(),
                frame_bases: full.frame_bases.clone(),
                row_max: full.row_max.clone(),
                ..full
            },
        };
        assert_same(0, &wide[0], reference);
    }
    let (stub, contexts) = q.x_into_stub(wide.iter().map(|g| &g.table));
    x_assert_stub_scalars(&stub, &full);
    x_assert_same_contexts(&contexts, &full_contexts);
    seen
}

/// Frames of `len` random letters (1 to 27: every letter but the gap) with a SEG-like mask
/// now and then.
fn frames(seed: &mut u64, n: usize, len: usize) -> Vec<QueryFrame> {
    (0..n)
        .map(|c| {
            let len = len + 13 * c;
            let masks = match c % 3 {
                0 => vec![(5, 12)],
                1 => vec![],
                _ => vec![(0, 3), (len / 2, len / 2 + 9)],
            };
            frame(residues(seed, len), masks, FRAMES[c])
        })
        .collect()
}

/// `f` on a one-thread pool and on pools of 4 and 8 threads.
fn on_pools(f: impl Fn(&crate::utils::threading::SearchPool<'_>) + Sync) {
    f(&crate::utils::threading::SearchPool::x_serial());
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    for threads in [4usize, 8] {
        crate::utils::threading::with_search_pool(threads, "tblastx", |pool| {
            assert_eq!(x_wide_threads(XWideOn::Pool(pool)), threads);
            f(pool);
            Ok(())
        })
        .expect("search pool");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:461-467,504-509,580-590
// (cited in the module comment): the chain of every cell, then BlastAaLookupFinalize's layout.
// Words of 2, 3 and 4 letters, four matrices, thresholds 0 (no neighbours), 11 and 13, one
// query of six contexts with SEG masks, a context shorter than a word, an empty one, a context
// of a repeated letter (cells far above AA_HITS_PER_CELL), and three queries of six contexts in
// groups of one, three and all contexts, on pools of 1, 4 and 8 threads.
#[test]
fn wide_tables_equal_the_filtered_full_table() {
    let mut seed = 0x3b17_5eed_0000_0004u64;
    let mut repeated = vec![1u8; 120];
    repeated.extend(residues(&mut seed, 40));
    let one_query = vec![vec![
        frame(residues(&mut seed, 260), vec![(10, 30), (200, 201)], 1),
        frame(residues(&mut seed, 3), vec![], 2),
        frame(repeated, vec![(70, 99)], 3),
        frame(residues(&mut seed, 230), vec![(0, 5)], -1),
        frame(Vec::new(), vec![], -2),
        frame(residues(&mut seed, 240), vec![], -3),
    ]];
    let three: Vec<Vec<QueryFrame>> = (0..3).map(|_| frames(&mut seed, 6, 60)).collect();
    let matrices = [
        ScoringMatrix::Blosum62,
        ScoringMatrix::Pam30,
        ScoringMatrix::Blosum45,
        ScoringMatrix::Pam250,
    ];
    let totals = std::sync::Mutex::new(Seen::default());
    on_pools(|pool| {
        for word in [2usize, 3, 4] {
            for (m, &matrix) in matrices.iter().enumerate() {
                for threshold in [0, 11, 13] {
                    // The larger thresholds at 4 letters only with the matrices of few
                    // neighbours there (BLOSUM45 and PAM250 at 11 make millions of entries).
                    if word == 4 && threshold == 11 && m >= 2 {
                        continue;
                    }
                    let seen = check(&one_query, threshold, matrix, word, None, pool);
                    // The 3-letter context has words of 2 and 3 letters, none of 4.
                    let filled = if word <= 3 { 5 } else { 4 };
                    assert_eq!(
                        seen.with_entries, filled,
                        "word {word} threshold {threshold}"
                    );
                    let mut t = totals.lock().expect("totals");
                    t.entries += seen.entries;
                    t.in_overflow += seen.in_overflow;
                    t.compacted += seen.compacted;
                    drop(t);
                    for starts in [None, Some(&[0, 3, 6, 9, 12, 15][..]), Some(&[0, 1, 7][..])] {
                        check(&three, threshold, matrix, word, starts, pool);
                    }
                    let seen = check(&three, threshold, matrix, word, Some(&[0]), pool);
                    assert_eq!(seen.with_entries, 1);
                }
            }
        }
    });
    let t = totals.into_inner().expect("totals");
    assert!(
        t.entries > 0 && t.in_overflow > 0 && t.compacted > 0,
        "{t:?}"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:504-509,580-582 (cited at
// `x_wide_targets`). The count of `XWideWordCount` is the number of cells `x_wide_targets` lists
// for the word, and the same for the word's letters in another order (the kept count of their
// set): words of 2, 3 and 4 letters, four matrices, thresholds 0, 11 and 13 (and 1 and INT_MIN
// at 2 and 3 letters), with letters outside the 20 amino acids; a batch's count is the sum of
// its words'.
#[test]
fn wide_word_count_is_the_cells_of_the_word() {
    let mut seed = 0x3b17_5eed_0000_0011u64;
    let matrices = [
        ScoringMatrix::Blosum62,
        ScoringMatrix::Pam30,
        ScoringMatrix::Blosum45,
        ScoringMatrix::Pam250,
    ];
    for word_length in [2usize, 3, 4] {
        for &matrix in &matrices {
            for threshold in [0, 11, 13, 1, i32::MIN] {
                if word_length == 4 && (threshold == 1 || threshold == i32::MIN) {
                    continue;
                }
                let counter = XWideWordCount::new(matrix, word_length, threshold)
                    .expect("a counter of every word length of the standard table");
                let mut query = residues(&mut seed, 60);
                // B, Z, X, U, * and the gap letter in the words too.
                query.extend_from_slice(&[2, 23, 21, 24, 25, 0, 1, 2]);
                let mut sum = 0u64;
                for word in query.windows(word_length) {
                    let mut cells = Vec::new();
                    x_wide_targets(word, word_length, &counter.scores, threshold, &mut cells);
                    let counted = counter.words(word);
                    assert_eq!(
                        counted,
                        cells.len() as u64,
                        "{matrix:?} word {word:?} threshold {threshold}"
                    );
                    let mut reversed = word.to_vec();
                    reversed.reverse();
                    assert_eq!(counter.words(&reversed), counted, "{word:?} reversed");
                    assert_eq!(counter.count(&reversed), counted, "{word:?} recounted");
                    sum += counted;
                }
                assert_eq!(
                    counter.batch_words(&mut [&query[..]].into_iter()),
                    sum,
                    "{matrix:?} word length {word_length} threshold {threshold}"
                );
            }
        }
    }
    assert_eq!(
        XWideWordCount::new(ScoringMatrix::Blosum62, 4, 11)
            .expect("counter")
            .backbone_cells(),
        compute_backbone_size(4, LOOKUP_ALPHABET_SIZE, X_WIDE_CHARSIZE) as u64
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:580-590 (cited in the
// module comment). Threshold 1 (most words are neighbours) and INT_MIN (every word, NCBI's
// conversion of a threshold beyond Int4): tables of the full backbone (more than half of the
// cells used), cells with long chains, at each word length.
#[test]
fn wide_tables_with_most_or_every_word_a_neighbour() {
    let mut seed = 0x3b17_5eed_0000_0005u64;
    let query = vec![vec![
        frame(residues(&mut seed, 9), vec![], 1),
        frame(residues(&mut seed, 8), vec![(2, 4)], 2),
        frame(residues(&mut seed, 3), vec![], 3),
        frame(residues(&mut seed, 7), vec![], -1),
    ]];
    on_pools(|pool| {
        for word in [2usize, 3, 4] {
            for matrix in [ScoringMatrix::Blosum62, ScoringMatrix::Pam30] {
                let seen = check(&query, 1, matrix, word, None, pool);
                assert!(seen.with_entries >= 3);
                let seen = check(&query, i32::MIN, matrix, word, None, pool);
                // Every valid word (28 letters) of every query word.
                let cells = 28usize.pow(word as u32);
                let words: usize = [9usize, 8, 3, 7]
                    .iter()
                    .enumerate()
                    .map(|(c, &len)| {
                        if c == 1 {
                            // The mask (2, 4) (letters 2 and 3) leaves [0, 1] and [4, 7].
                            [2usize, 4]
                                .iter()
                                .map(|&l| (l + 1).saturating_sub(word))
                                .sum()
                        } else {
                            (len + 1).saturating_sub(word)
                        }
                    })
                    .sum();
                assert_eq!(seen.entries, words * cells, "word {word}");
                // Every table with entries holds every valid word: of the full backbone; the
                // context of 3 letters has no 4-letter word (an empty table, compacted).
                let filled = if word <= 3 { 4 } else { 3 };
                assert_eq!(
                    (seen.with_entries, seen.compacted),
                    (filled, 4 - filled),
                    "word {word}: every table with entries of the full backbone"
                );
                check(&query, i32::MIN, matrix, word, Some(&[0]), pool);
            }
        }
    });
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1045-1047 (cited at
// `x_group_locations`). Contexts masked from end to end, groups without a location, and a
// query without any word: the tables without entries are empty in both builders.
#[test]
fn wide_tables_of_masked_and_empty_groups() {
    let mut seed = 0x3b17_5eed_0000_0006u64;
    let masked: Vec<Vec<QueryFrame>> = vec![(0..6)
        .map(|c| {
            let masks = if c % 2 == 0 {
                vec![(0, 50)]
            } else {
                vec![(3, 9)]
            };
            frame(residues(&mut seed, 50), masks, FRAMES[c])
        })
        .collect()];
    let none: Vec<Vec<QueryFrame>> = vec![(0..6)
        .map(|c| frame(residues(&mut seed, 30), vec![(0, 30)], FRAMES[c]))
        .collect()];
    on_pools(|pool| {
        for word in [2usize, 3, 4] {
            let seen = check(&masked, 13, ScoringMatrix::Blosum62, word, None, pool);
            assert_eq!(seen.with_entries, 3);
            let seen = check(&none, 13, ScoringMatrix::Blosum62, word, None, pool);
            assert_eq!(
                (seen.entries, seen.with_entries, seen.in_overflow),
                (0, 0, 0)
            );
        }
    });
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:461-467 (cited in the
// module comment).
// Item SC9 of line NP: the table of 4-letter words of one list of all the lookup locations
// (`x_build_full_wide`, from a lookup query handed over, `x_from_prepared`; the query here is
// TBLASTX's, BLASTX's own is `blastx::seed`'s test) is the reference builder's full table field by
// field, never compacted, built on this thread and in 4 parts on pools of 1, 4 and 8 threads,
// also with stripes of 32 and of 2^16 cells; with a query of SEG
// masks, a context shorter than a word, an empty one and a repeated letter, three matrices and
// thresholds 0, 11 and 13.
#[test]
fn blastx_s_full_table_of_four_letters_is_the_reference_table() {
    let mut seed = 0x3b17_5eed_0000_0009u64;
    let mut repeated = vec![1u8; 120];
    repeated.extend(residues(&mut seed, 40));
    let query = vec![vec![
        frame(residues(&mut seed, 260), vec![(10, 30), (200, 201)], 1),
        frame(residues(&mut seed, 3), vec![], 2),
        frame(repeated, vec![(70, 99)], 3),
        frame(residues(&mut seed, 230), vec![(0, 5)], -1),
        frame(Vec::new(), vec![], -2),
        frame(residues(&mut seed, 240), vec![], -3),
    ]];
    on_pools(|pool| {
        for matrix in [
            ScoringMatrix::Blosum62,
            ScoringMatrix::Pam30,
            ScoringMatrix::Blosum45,
        ] {
            for threshold in [0, 11, 13] {
                let params = crate::stats::lookup_protein_params_ungapped(matrix);
                let (full, _) = crate::algorithm::tblastx::lookup::build_ncbi_lookup_for_profile(
                    &query, threshold, &params, true, matrix, 4,
                );
                let prepared = XPairLutQuery::x_prepare_for(&query, threshold, true, matrix, 4)
                    .expect("the joint builder covers 4 letters")
                    .prepared;
                let q = XPairLutQuery::x_from_prepared(prepared, threshold, matrix, 4);
                for parts in [1, 4] {
                    let table = pool.install(|| q.x_build_full_wide(parts));
                    x_assert_same_table_as(
                        &format!("{matrix:?} threshold {threshold} parts {parts}:"),
                        &table,
                        &full,
                    );
                    for bits in [5, 16] {
                        let tables = [q.prepared.lookup_locations.clone()];
                        let mut built = pool.install(|| {
                            q.x_build_tables_wide_with(
                                &tables,
                                XWideOn::Here(parts),
                                false,
                                Some(bits),
                                false,
                            )
                        });
                        x_assert_same_table_as(
                            &format!(
                                "{matrix:?} threshold {threshold} parts {parts} stripes of \
                                     2^{bits} cells:"
                            ),
                            &built.pop().expect("one table").table,
                            &full,
                        );
                    }
                }
            }
        }
    });
}

// No NCBI counterpart: the word lengths each builder takes.
#[test]
fn the_builders_cover_words_of_two_to_four_letters() {
    assert!(!x_wide_covers(0));
    assert!((1..=4).all(x_wide_covers));
    assert!(!x_wide_covers(5));
    assert!((2..=4).all(x_tables_cover));
    assert!(x_direct_covers(3) && !x_direct_covers(4));
}
