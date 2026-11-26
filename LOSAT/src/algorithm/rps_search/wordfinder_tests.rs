use super::*;

#[test]
fn the_table_is_the_power_of_two_past_the_concatenation_and_the_window() {
    assert_eq!(diag_array_length(5035, 40), 8192);
    assert_eq!(diag_array_length(1130, 40), 2048);
    assert_eq!(diag_array_length(984, 40), 1024);
    assert_eq!(diag_array_length(985, 40), 2048);
    assert_eq!(diag_array_length(1024, 0), 1024);
    // `Lc` 0 and below: 1 entry where the sum with the window is 1 or less, else the power
    // of 2 past the sum (`Lc` -2 and the window 40: 64).
    assert_eq!(diag_array_length(-2, 40), 64);
    assert_eq!(diag_array_length(-2, 0), 1);
    assert_eq!(diag_array_length(1, 0), 1);
}

// Review B of D1a-4, SF2: NCBI's length, or the power of 2 past the rows and the longest
// query where that is shorter.
#[test]
fn the_table_is_ncbis_or_one_whose_diagonals_are_apart() {
    // `Lc` 2^30 - 40 (an alias file's LENGTH over the 563 rows of rps_mini) and q40.faa's
    // longest query of 2,351 letters: 4,096 entries in place of 2^30.
    assert_eq!(diag_table_length((1 << 30) - 40, 40, 563, 2351), 4096);
    // NCBI's is shorter (the diagonals 4,096 apart share its entries): kept.
    assert_eq!(diag_table_length(2100, 40, 2100, 4200), 4096);
    // `Lc` below the rows (LENGTH 400 over 563 rows): NCBI's 512.
    assert_eq!(diag_table_length(403, 40, 563, 2351), 512);
    // Queries shorter than the window: the rows and the query need 1,024 of NCBI's 2,048.
    assert_eq!(diag_table_length(1000, 40, 1000, 20), 1024);
    assert_eq!(diag_table_length(1000, 40, 1000, 60), 2048);
}

#[test]
fn the_entry_keeps_ncbis_31_bit_signed_field_and_flag() {
    let mut entry = 0u32;
    entry.set(-40, true);
    assert_eq!((entry.last_hit(), entry.flag()), (-40, true));
    entry.set(i32::MAX / 4 + 10, false);
    assert_eq!((entry.last_hit(), entry.flag()), (i32::MAX / 4 + 10, false));
}

// Diagonals `L` apart share an entry (`(q - s) & mask`): across buckets, a hit of a later
// bucket at a smaller `s` reads the entry that an earlier bucket's hit at a larger `s` wrote
// (exactness point 1 of the design). A hit of bucket 0 on diagonal -4180 makes the two hits of
// bucket 1 on diagonal -84 (4,096 apart) read as overlapping it; with a table of 8,192 entries
// the second of them extends.
#[test]
fn diagonals_a_table_length_apart_share_an_entry_across_buckets() {
    use crate::seqdb::rps::RPS_BACKBONE_SIZE;
    let concat = 2100usize;
    let mut rows = vec![[3i32; 28]; concat];
    rows[concat - 1] = [-32767; 28];
    let word = |l: usize| (l << 10) | (l << 5) | l;
    let mut cells = vec![[0i32; 4]; RPS_BACKBONE_SIZE + 1];
    let mut pv = vec![0u64; RPS_BACKBONE_SIZE / 64];
    // The last letter's rows of the words A A A (q 10), C C C (q 2060), D D D (q 2065).
    for (letter, entry) in [(1usize, 12), (3, 2062), (4, 2067)] {
        cells[word(letter)] = [1, entry, 0, 0];
        pv[word(letter) / 64] |= 1 << (word(letter) % 64);
    }
    let mut query = vec![20u8; 4200];
    query[4190..4193].fill(1);
    query[2144..2147].fill(3);
    query[2149..2152].fill(4);
    let params = RpsWordParams {
        rows: RpsRows {
            rows: &rows,
            length: rows.len(),
            past: None,
        },
        table: RpsTable {
            cells: &cells,
            overflow: &[],
            pv: &pv,
        },
        window: 40,
        cutoff: 1,
        dropoff: 20,
    };
    let hits_with_table = |length: i32| {
        let mut scratch = RpsWordScratch {
            diag: RpsDiagTable::new(),
            buckets: RpsBuckets::new(concat, Default::default()),
        };
        scratch.diag.reset_for(diag_array_length(length, 40), 40);
        let mut hits = Vec::new();
        word_finder(&params, &mut scratch, &query, &mut hits);
        (scratch.diag.mask, hits)
    };
    let (mask, ncbi) = hits_with_table(concat as i32);
    assert_eq!(mask, 4095);
    assert!(ncbi.is_empty(), "{ncbi:?}");
    let (wide_mask, wide) = hits_with_table(8000);
    assert_eq!(wide_mask, 8191);
    assert_eq!(wide.len(), 1);
    assert_eq!((wide[0].q_off, wide[0].s_off), (2065, 2149));
}

// Stale entries of earlier queries read as no hit: a table that searched another query
// gives the hits of a fresh one (exactness point 2 of the design).
#[test]
fn a_used_table_gives_the_hits_of_a_fresh_one() {
    use crate::seqdb::rps::RPS_BACKBONE_SIZE;
    // One profile of 8 positions that scores +5 on letter 1 (A) and -4 elsewhere, then its
    // sentinel row; every word A A A hits the rows of the profile.
    let mut rows = vec![[-4i32; 28]; 9];
    for row in rows.iter_mut().take(8) {
        row[1] = 5;
    }
    rows[8] = [-32767; 28];
    let mut cells = vec![[0i32; 4]; RPS_BACKBONE_SIZE + 1];
    let mut pv = vec![0u64; RPS_BACKBONE_SIZE / 64];
    let w = (1 << 10) | (1 << 5) | 1;
    let overflow: Vec<i32> = vec![3, 4, 5, 6, 7];
    cells[w] = [6, 2, 0, 0];
    pv[w / 64] |= 1 << (w % 64);
    let params = RpsWordParams {
        rows: RpsRows {
            rows: &rows,
            length: rows.len(),
            past: None,
        },
        table: RpsTable {
            cells: &cells,
            overflow: &overflow,
            pv: &pv,
        },
        window: 40,
        cutoff: 1,
        dropoff: 20,
    };
    let query = [1u8, 1, 1, 1, 1, 1, 20, 1, 1, 1, 1];
    let other = [1u8; 30];
    let run = |before: &[&[u8]]| {
        let mut scratch = RpsWordScratch {
            diag: RpsDiagTable::new(),
            buckets: RpsBuckets::new(rows.len(), Default::default()),
        };
        scratch
            .diag
            .reset_for(diag_array_length(rows.len() as i32, 40), 40);
        let mut ignored = Vec::new();
        for q in before {
            word_finder(&params, &mut scratch, q, &mut ignored);
        }
        let mut hits = Vec::new();
        word_finder(&params, &mut scratch, &query, &mut hits);
        hits
    };
    let fresh = run(&[]);
    assert!(!fresh.is_empty());
    assert_eq!(run(&[&other]), fresh);
    assert_eq!(run(&[&other, &query, &other]), fresh);
    // NCBI's clear at a quarter of `INT4_MAX` (`RpsDiagTable::exit`) on a kept table (D3 of PS12s):
    // a worker's offset reaches it after about 536 M letters of frames. A table whose offset is
    // just below it searches `other` (its entries then hold values past the limit), the offset
    // passes the limit, `query` still gets a fresh table's hits, and its exit clears the table
    // (every entry `-window`, the offset the window), after which `query` gets them again.
    let mut scratch = RpsWordScratch {
        diag: RpsDiagTable::new(),
        buckets: RpsBuckets::new(rows.len(), Default::default()),
    };
    scratch
        .diag
        .reset_for(diag_array_length(rows.len() as i32, 40), 40);
    scratch.diag.offset = i32::MAX / 4 - 5;
    let mut ignored = Vec::new();
    word_finder(&params, &mut scratch, &other, &mut ignored);
    assert!(scratch.diag.offset >= i32::MAX / 4);
    assert!(scratch
        .diag
        .array
        .iter()
        .any(|entry| entry.last_hit() > i32::MAX / 4));
    let mut hits = Vec::new();
    word_finder(&params, &mut scratch, &query, &mut hits);
    assert_eq!(hits, fresh);
    assert_eq!(scratch.diag.offset, 40);
    assert!(scratch
        .diag
        .array
        .iter()
        .all(|entry| (entry.last_hit(), entry.flag()) == (-40, false)));
    let mut hits = Vec::new();
    word_finder(&params, &mut scratch, &query, &mut hits);
    assert_eq!(hits, fresh);
}

// D3 of PS12s: a worker's table kept from one volume to the next of the same table length
// gives a query the hits of a fresh table (the offset rule), and a volume of another length
// gets a zeroed table of its own length.
#[test]
fn a_table_kept_from_another_volume_gives_the_hits_of_a_fresh_one() {
    use crate::seqdb::rps::RPS_BACKBONE_SIZE;
    // Two volumes whose tables have 1,024 entries (984 and 900 rows), each one profile scoring
    // +5 on letter `letter` and -4 elsewhere; the word of three such letters hits every row.
    let volume = |rows_count: usize, letter: usize| {
        let mut rows = vec![[-4i32; 28]; rows_count];
        for row in rows.iter_mut().take(rows_count - 1) {
            row[letter] = 5;
        }
        rows[rows_count - 1] = [-32767; 28];
        let mut cells = vec![[0i32; 4]; RPS_BACKBONE_SIZE + 1];
        let mut pv = vec![0u64; RPS_BACKBONE_SIZE / 64];
        let w = (letter << 10) | (letter << 5) | letter;
        let entries: Vec<i32> = (2..rows_count as i32 - 1).step_by(37).collect();
        cells[w] = [entries.len() as i32, entries[0], 0, 0];
        pv[w / 64] |= 1 << (w % 64);
        (rows, cells, pv, entries[1..].to_vec())
    };
    let (rows_a, cells_a, pv_a, overflow_a) = volume(984, 1);
    let (rows_b, cells_b, pv_b, overflow_b) = volume(900, 3);
    let a = RpsWordParams {
        rows: RpsRows {
            rows: &rows_a,
            length: rows_a.len(),
            past: None,
        },
        table: RpsTable {
            cells: &cells_a,
            overflow: &overflow_a,
            pv: &pv_a,
        },
        window: 40,
        cutoff: 1,
        dropoff: 20,
    };
    let b = RpsWordParams {
        rows: RpsRows {
            rows: &rows_b,
            length: rows_b.len(),
            past: None,
        },
        table: RpsTable {
            cells: &cells_b,
            overflow: &overflow_b,
            pv: &pv_b,
        },
        ..a
    };
    let query_a = [1u8; 40];
    let query_b: Vec<u8> = (0..60).map(|i| if i % 7 == 6 { 20 } else { 3 }).collect();
    let hits = |scratch: &mut RpsWordScratch, params: &RpsWordParams<'_>, rows: usize, q: &[u8]| {
        scratch.prepare_for(rows, rows as i32, 40, q.len());
        let mut hits = Vec::new();
        word_finder(params, scratch, q, &mut hits);
        hits
    };
    let fresh_b = hits(
        &mut RpsWordScratch::new(Default::default()),
        &b,
        900,
        &query_b,
    );
    assert!(!fresh_b.is_empty());
    let mut scratch = RpsWordScratch::new(Default::default());
    assert!(!hits(&mut scratch, &a, 984, &query_a).is_empty());
    assert_eq!(scratch.diag.array.len(), 1024);
    let kept = scratch.diag.array.as_ptr();
    assert_eq!(hits(&mut scratch, &b, 900, &query_b), fresh_b);
    assert_eq!(scratch.diag.array.as_ptr(), kept, "the table is kept");
    // A volume whose table is larger than the bucket ring uses the bucket-sliding ring
    // `(RPS_BUCKET_SIZE + longest).next_power_of_two()` (4,096 instead of 8,192).
    scratch.prepare_for(5035, 5035, 40, 60);
    assert_eq!(scratch.diag.array.len(), 4096);
    assert_eq!(scratch.diag.offset, 40);
    assert!(scratch.diag.array.iter().all(|&entry| entry == 0));
}

// PS12m: the word finders read the hits in passes of buckets; the initial HSPs (and their
// order) are the same for one pass of every bucket and for passes of 1, 2 and 64 buckets, the
// two-hit and the one-hit finder, the diagonals of a long query sharing the table's entries.
#[test]
fn the_initial_hsps_do_not_depend_on_the_groups_of_buckets() {
    use crate::seqdb::rps::RPS_BACKBONE_SIZE;
    let rows_count = 9_000usize;
    let mut rows = vec![[-4i32; 28]; rows_count];
    for row in rows.iter_mut().take(rows_count - 1) {
        row[1] = 5;
    }
    rows[rows_count - 1] = [-32767; 28];
    let mut cells = vec![[0i32; 4]; RPS_BACKBONE_SIZE + 1];
    let mut pv = vec![0u64; RPS_BACKBONE_SIZE / 64];
    let w = (1 << 10) | (1 << 5) | 1;
    // The word A A A hits every 37th row, in ascending order (five buckets, one overflow run).
    let entries: Vec<i32> = (2..rows_count as i32 - 1).step_by(37).collect();
    cells[w] = [entries.len() as i32, entries[0], 0, 0];
    pv[w / 64] |= 1 << (w % 64);
    let overflow = entries[1..].to_vec();
    let mut query = Vec::new();
    for run in [30usize, 12, 55, 9, 80] {
        query.extend(std::iter::repeat_n(1u8, run));
        query.extend([20u8, 20]);
    }
    let hits_of = |window: i32, group_buckets: usize| {
        let params = RpsWordParams {
            rows: RpsRows {
                rows: &rows,
                length: rows.len(),
                past: None,
            },
            table: RpsTable {
                cells: &cells,
                overflow: &overflow,
                pv: &pv,
            },
            window,
            cutoff: 1,
            dropoff: 20,
        };
        let mut scratch = RpsWordScratch {
            diag: RpsDiagTable::new(),
            buckets: RpsBuckets::with_group(rows_count, Default::default(), group_buckets),
        };
        scratch.diag.reset_for(rows_count as i32, window);
        let mut hits = Vec::new();
        word_finder(&params, &mut scratch, &query, &mut hits);
        hits
    };
    for window in [40, 0] {
        let whole = hits_of(window, usize::MAX / 2);
        assert!(whole.len() > 10, "window {window}: {} HSPs", whole.len());
        for group_buckets in [1, 2, 64] {
            assert_eq!(
                hits_of(window, group_buckets),
                whole,
                "{window} {group_buckets}"
            );
        }
    }
}

// Bucket-sliding ring parity: for both two-hit and one-hit word finders, across multiple groups,
// empty-bucket gaps, multiple queries reusing one worker scratch, and short/long subjects where
// `rows + subject.len() <= diag_array_length`, the bucket-sliding ring produces the exact same
// `RpsInitHsp` list as NCBI's full diagonal table (after `sort_init_hsps_by_score`, and in fact
// before sorting too when a group completes in one pass).
#[test]
fn bucket_sliding_ring_matches_full_table_for_two_hit_and_one_hit() {
    use crate::seqdb::rps::RPS_BACKBONE_SIZE;
    // 140,000 rows (> 64 buckets, so 2 groups of 64 buckets; NCBI table length = 262,144,
    // bucket ring = 4,096 for subjects <= 2,048 aa).
    let rows_count = 140_000usize;
    let mut rows = vec![[-4i32; 28]; rows_count];
    for row in rows.iter_mut().take(rows_count - 1) {
        row[1] = 5;
        row[2] = 4;
    }
    rows[rows_count - 1] = [-32767; 28];
    let mut cells = vec![[0i32; 4]; RPS_BACKBONE_SIZE + 1];
    let mut pv = vec![0u64; RPS_BACKBONE_SIZE / 64];
    let mut overflow = Vec::new();
    // Place hits across buckets 0, 1, 2 (adjacent), skip buckets 3..10 (empty gap > ring slots),
    // hit buckets 10, 11 (adjacent across gap), and cross group boundary at bucket 63 -> 64 -> 68.
    // Step 11 (< window 40) creates dense two-hit pairs on the same diagonals across bucket bounds.
    let mut entries_a = Vec::new();
    for &b in &[0usize, 1, 2, 10, 11, 63, 64, 68] {
        let base = b * RPS_BUCKET_SIZE;
        for step in (0..RPS_BUCKET_SIZE).step_by(11) {
            let pos = base + step;
            if pos >= 2 && pos + 1 < rows_count {
                entries_a.push(pos as i32);
            }
        }
    }
    let w1 = (1 << 10) | (1 << 5) | 1;
    cells[w1] = [entries_a.len() as i32, entries_a[0], 0, 0];
    overflow.extend_from_slice(&entries_a[1..]);
    pv[w1 / 64] |= 1 << (w1 % 64);
    // Inline word (<= 3 hits) spanning group 0 and group 1.
    let w2 = (2 << 10) | (2 << 5) | 2;
    cells[w2] = [3, 2045, 2055, (64 * RPS_BUCKET_SIZE + 12) as i32];
    pv[w2 / 64] |= 1 << (w2 % 64);

    let mut q1 = Vec::new();
    for run in [25usize, 9, 45, 11, 60] {
        q1.extend(std::iter::repeat_n(1u8, run));
        q1.extend([2u8, 2, 2, 20, 20]);
    }
    let mut q2 = vec![20u8; 2_500];
    for &pos in &[0usize, 500, 1_500, 2_040, 2_450] {
        q2[pos..pos + 30].fill(1);
    }

    for window in [40i32, 0] {
        let params = RpsWordParams {
            rows: RpsRows {
                rows: &rows,
                length: rows.len(),
                past: None,
            },
            table: RpsTable {
                cells: &cells,
                overflow: &overflow,
                pv: &pv,
            },
            window,
            cutoff: 1,
            dropoff: 20,
        };
        let ncbi_len = diag_array_length(rows_count as i32, window);
        let mut full_scratch = RpsWordScratch {
            diag: RpsDiagTable::new(),
            buckets: RpsBuckets::new(rows_count, Default::default()),
        };
        full_scratch.diag.reset_for(ncbi_len, window);

        let mut ring_scratch = RpsWordScratch::new(Default::default());
        ring_scratch.prepare_for(rows_count, rows_count as i32, window, q2.len());
        assert_eq!(ring_scratch.diag.array.len(), 8192);
        assert!(ring_scratch.diag.array.len() < ncbi_len as usize);

        for query in [&q1[..], &q2[..], &q1[..]] {
            let mut full_hits = Vec::new();
            word_finder(&params, &mut full_scratch, query, &mut full_hits);
            let mut ring_hits = Vec::new();
            word_finder(&params, &mut ring_scratch, query, &mut ring_hits);
            assert!(!full_hits.is_empty());
            assert_eq!(ring_hits, full_hits, "window={window} len={}", query.len());

            // Also force the multi-pass group fallback (`max_hits = 50`) for both group 0 and
            // group 1 and verify sorted init_hsps match the full table bit-for-bit.
            let mut mp_scratch = RpsWordScratch::new(Default::default());
            mp_scratch.prepare_for(rows_count, rows_count as i32, window, query.len());
            let mut mp_hits = Vec::new();
            let ring_len = bucket_ring_length(query.len()) as usize;
            if window > 0 {
                run_sliding::<_, true>(
                    &params,
                    &params.rows,
                    &mut mp_scratch,
                    query,
                    ring_len,
                    50,
                    &mut mp_hits,
                );
            } else {
                run_sliding::<_, false>(
                    &params,
                    &params.rows,
                    &mut mp_scratch,
                    query,
                    ring_len,
                    50,
                    &mut mp_hits,
                );
            }
            let sort_key = |h: &RpsInitHsp| {
                (
                    std::cmp::Reverse(h.ungapped.score),
                    h.ungapped.s_start,
                    std::cmp::Reverse(h.ungapped.length),
                    h.ungapped.q_start,
                )
            };
            full_hits.sort_by_key(sort_key);
            mp_hits.sort_by_key(sort_key);
            assert_eq!(
                mp_hits,
                full_hits,
                "multipass window={window} len={}",
                query.len()
            );
        }
    }
}
