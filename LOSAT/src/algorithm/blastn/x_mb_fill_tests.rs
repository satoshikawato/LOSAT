use super::super::compute_mb_pv_params;
use super::super::x_mb_lookup_par::tests::{random_case, repeat_case, Lcg};
use super::*;

#[derive(PartialEq, Debug)]
struct Built {
    hashtable: Vec<u32>,
    next_pos: Vec<u32>,
    pv_array: Vec<u32>,
    helper_array: Vec<u32>,
    words: Vec<(u32, u32)>,
    total_positions: usize,
    ambiguous_skipped: usize,
}

/// The arrays a fill makes; `h24`: in 3-byte elements (LOSAT_X_MB24), returned as 4-byte values.
fn build(words: &Words<'_>, h24: bool, fill: impl FnOnce(&mut XMbTables<'_>)) -> Built {
    let table_size = words.kmer_mask as usize + 1;
    let entries: usize = words.queries_blastna.iter().map(|q| q.len()).sum();
    let (pv_size, pv_array_bts) =
        compute_mb_pv_params(table_size, entries, false, words.lut_word_length);
    let next_len = words
        .queries_blastna
        .iter()
        .zip(words.query_offsets)
        .map(|(q, &o)| o.max(0) as usize + q.len())
        .max()
        .unwrap_or(0)
        + 1;
    let width = |n: usize| if h24 { x_mb_narrow::words_for(n) } else { n };
    let mut b = Built {
        hashtable: vec![0; width(table_size)],
        next_pos: vec![0; width(next_len)],
        pv_array: vec![0; pv_size],
        helper_array: vec![0; (table_size >> K_COMPRESSION_FACTOR_SHIFT).max(1)],
        words: Vec::new(),
        total_positions: 0,
        ambiguous_skipped: 0,
    };
    fill(&mut XMbTables {
        hashtable: &mut b.hashtable,
        next_pos: &mut b.next_pos,
        pv_array: &mut b.pv_array,
        pv_array_bts,
        helper_array: &mut b.helper_array,
        words: &mut b.words,
        total_positions: &mut b.total_positions,
        ambiguous_skipped: &mut b.ambiguous_skipped,
        x_h24: h24,
    });
    if h24 {
        b.hashtable = x_mb_narrow::widen(&b.hashtable, table_size);
        b.next_pos = x_mb_narrow::widen(&b.next_pos, next_len);
    }
    b
}

/// Both queues (prefetch and batch) and both cell orders build the arrays and counters of
/// `x_fill_mb_queued` on random queries (ambiguity codes, masks, several queries,
/// `lut_word_length` 4..12, runs shorter and longer than the queue), on single queries with
/// fewer words than the queue holds (`lut_word_length` 6, 7, 10, 12) and on long repeats (long
/// chains, helper counts); newest first, so do the instances for 3-byte elements (LOSAT_X_MB24)
/// of both and of the queued loop.
#[test]
fn the_run_loop_builds_the_queued_loops_arrays() {
    let mut rng = Lcg(0x4D42_F111);
    for case in 0..60 {
        let (queries, offsets, masks, word_length, lut_word_length) = if case % 6 == 4 {
            // fewer words than the queue holds, at the table widths `random_case` leaves out
            let lut = [6usize, 7, 10, 12][(case / 6) % 4];
            let len = lut + rng.below(60) as usize;
            let seq: Vec<u8> = (0..len)
                .map(|_| {
                    if rng.below(30) == 0 {
                        14
                    } else {
                        rng.below(4) as u8
                    }
                })
                .collect();
            (vec![seq], vec![0], vec![Vec::new()], lut, lut)
        } else if case % 6 == 5 {
            let lut = [11usize, 12][case % 2];
            (
                vec![repeat_case(&mut rng, 9000)],
                vec![0],
                vec![Vec::new()],
                lut + 16,
                lut,
            )
        } else {
            random_case(&mut rng)
        };
        let refs: Vec<&[u8]> = queries.iter().map(|q| q.as_slice()).collect();
        for ascending_cells in [false, true] {
            let words = Words {
                queries_blastna: &refs,
                query_offsets: &offsets,
                query_masks: &masks,
                word_length,
                lut_word_length,
                kmer_mask: (1u64 << (2 * lut_word_length)) - 1,
                ascending_cells,
            };
            let reference = build(&words, false, |t| {
                x_fill_mb_queued::<false>(
                    x_mb_queue(),
                    words.queries_blastna,
                    words.query_offsets,
                    words.query_masks,
                    word_length,
                    lut_word_length,
                    words.kmer_mask,
                    None,
                    0,
                    ascending_cells,
                    t,
                )
            });
            assert!(reference.total_positions > 0 || queries.iter().all(|q| q.len() < 80));
            for (batch, h24) in [(false, false), (true, false), (false, true), (true, true)] {
                if h24 && ascending_cells {
                    continue;
                }
                let got = build(&words, h24, |t| {
                    assert!(in_range(&words, t));
                    // SAFETY: checked above.
                    let (n, a) = unsafe {
                        match (ascending_cells, batch, h24) {
                            (true, false, _) => fill_runs::<true, false, false>(&words, t),
                            (false, false, false) => fill_runs::<false, false, false>(&words, t),
                            (true, true, _) => fill_runs::<true, true, false>(&words, t),
                            (false, true, false) => fill_runs::<false, true, false>(&words, t),
                            (false, false, true) => fill_runs::<false, false, true>(&words, t),
                            (false, true, true) => fill_runs::<false, true, true>(&words, t),
                        }
                    };
                    *t.total_positions += n;
                    *t.ambiguous_skipped += a;
                });
                assert!(
                    got == reference,
                    "case {case} ascending {ascending_cells} batch {batch} h24 {h24}"
                );
            }
            if !ascending_cells {
                for queue in [XMbQueue::Delay, XMbQueue::Batch(64)] {
                    let got = build(&words, true, |t| {
                        x_fill_mb_queued::<true>(
                            queue,
                            words.queries_blastna,
                            words.query_offsets,
                            words.query_masks,
                            word_length,
                            lut_word_length,
                            words.kmer_mask,
                            None,
                            0,
                            false,
                            t,
                        )
                    });
                    assert!(got == reference, "case {case} queued h24");
                }
            }
        }
    }
}
