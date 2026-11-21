use super::*;
use crate::algorithm::blastp::encoding::EncodedProtein;
use crate::utils::glibc_fma_libm::GlibcLibm;

const ECUTOFF: f64 = 0.002;

fn hsp(oid: u32, e_value: f64) -> BlastpHsp {
    BlastpHsp {
        num: 1,
        identity: 50.0,
        length: 10,
        mismatch: 5,
        gapopen: 0,
        q_start: 1,
        q_end: 10,
        s_start: 1,
        s_end: 10,
        gapped_q_start: 0,
        gapped_s_start: 0,
        e_value,
        bit_score: 20.0,
        num_ident: 5,
        query_context: 0,
        query_frame: 0,
        subject_frame: 0,
        query_length: 10,
        q_idx: 0,
        s_idx: oid,
        raw_score: 40,
        gap_info: None,
        num_positives: 6,
        comp_adjustment_method: 2,
    }
}

/// A list of subject `oid` with best E-value `best_evalue`; without HSPs when `empty`.
fn list(oid: u32, best_evalue: f64, empty: bool) -> BlastpHspList {
    BlastpHspList {
        oid,
        query_index: 0,
        hsps: if empty {
            Vec::new()
        } else {
            vec![hsp(oid, best_evalue)]
        },
        best_evalue,
    }
}

fn hit_list(lists: Vec<BlastpHspList>) -> Option<BlastpHitList> {
    Some(BlastpHitList {
        hsplist_count: lists.len(),
        hsplist_max: 1050,
        worst_evalue: 0.0,
        low_score: 0,
        heapified: false,
        hsplist_array: lists,
        hsplist_current: 0,
        num_hits: 0,
    })
}

fn redone(oid: u32, best_evalue: f64) -> Result<BlastpPostprocessResult> {
    Ok(BlastpPostprocessResult {
        hsp_list: Some(list(oid, best_evalue, false)),
        best_score: 40,
        best_evalue,
    })
}

fn heaps_of(count: usize, heap_threshold: usize) -> Vec<BlastCompoHeap> {
    (0..count)
        .map(|_| BlastCompoHeap::new(heap_threshold, ECUTOFF))
        .collect()
}

/// Each query alone in its batch.
fn alone(count: usize) -> Vec<Range<usize>> {
    (0..count).map(|q_idx| q_idx..q_idx + 1).collect()
}

fn at(q_idx: usize, position: usize) -> XMatch {
    XMatch { q_idx, position }
}

/// The window from `start` of at most `jobs` jobs, without a cap on jobs at risk.
fn window_from(
    lists: &[Option<BlastpHitList>],
    heaps: &[BlastCompoHeap],
    batches: &[Range<usize>],
    end_above: &[Option<f64>],
    start: XMatch,
    jobs: usize,
) -> XWindow {
    let caps = vec![usize::MAX; lists.len()];
    x_window(lists, heaps, batches, end_above, &caps, start, jobs)
}

impl XWindow {
    fn at_risk_jobs(&self) -> usize {
        self.at_risk.iter().map(|&(_, jobs)| jobs).sum()
    }
}

fn replayed(skipped: usize, wasted: usize, dropped: &[usize]) -> XReplayed {
    XReplayed {
        skipped,
        wasted,
        dropped: dropped.to_vec(),
    }
}

/// All matches of `lists` in stream order.
fn stream(lists: &[Option<BlastpHitList>]) -> Vec<XMatch> {
    let window = window_from(
        lists,
        &heaps_of(lists.len(), 1),
        &alone(lists.len()),
        &vec![None; lists.len()],
        at(0, 0),
        usize::MAX,
    );
    assert_eq!(window.next, None);
    window.matches
}

#[test]
fn jobs_hold_up_to_four_matches_of_one_query_and_windows_end_at_their_job_budget() {
    let lists = vec![
        hit_list(vec![
            list(1, 1e-9, false),
            list(2, 1e-8, true),
            list(3, 1e-7, false),
            list(7, 1e-7, false),
            list(8, 1e-7, false),
            list(9, 1e-6, false),
            list(10, 1e-6, false),
        ]),
        None,
        hit_list(vec![list(4, 1e-9, false)]),
        hit_list(vec![list(5, 1e-9, false), list(6, 1e-8, false)]),
        None,
    ];
    assert_eq!(x_matches(&lists, &[]), 9);
    let (heaps, batches, none) = (heaps_of(5, 1), alone(5), vec![None; 5]);
    let window = |start, jobs| window_from(&lists, &heaps, &batches, &none, start, jobs);
    // The list without HSPs (query 0, position 1) is no match.
    let all = window(at(0, 0), usize::MAX);
    assert_eq!(
        all.matches,
        [
            at(0, 0),
            at(0, 2),
            at(0, 3),
            at(0, 4),
            at(0, 5),
            at(0, 6),
            at(2, 0),
            at(3, 0),
            at(3, 1)
        ]
    );
    assert_eq!(all.jobs, [0..4, 4..6, 6..7, 7..9]);
    assert_eq!((all.next, all.ended), (None, vec![false; 9]));
    // One job: the window ends inside query 0, and the next one starts at its fifth match
    // with the same job as above.
    let first = window(at(0, 0), 1);
    assert_eq!((first.jobs, first.next), (vec![0..4], Some(at(0, 5))));
    let second = window(at(0, 5), 1);
    assert_eq!(second.matches, [at(0, 5), at(0, 6)]);
    assert_eq!(second.next, Some(at(2, 0)));
    // The queries without matches after the last match end the stream.
    let last = window(at(2, 0), 2);
    assert_eq!((last.jobs, last.next), (vec![0..1, 1..3], None));
    assert!(window(at(4, 0), 1).jobs.is_empty());
}

/// Two queries of one batch, one subject each to keep (`-max_target_seqs 1`).
fn two_query_batch() -> (
    Vec<Option<BlastpHitList>>,
    Vec<BlastCompoHeap>,
    Vec<Range<usize>>,
) {
    let lists = vec![
        hit_list(vec![list(10, 1e-20, false)]),
        hit_list(vec![
            list(20, 1e-20, false),
            list(21, 1.0, false),
            list(22, 0.005, false),
        ]),
    ];
    (lists, heaps_of(2, 1), vec![0..2, 0..2])
}

#[test]
fn the_replay_drops_what_the_early_termination_test_ends_and_reports_the_next_error() {
    let (lists, mut heaps, batches) = two_query_batch();
    let matches = stream(&lists);
    assert_eq!(matches.len(), 4);
    // The match of best E 1.0 comes after both heaps are filled to the cutoff: the test
    // ends it, so its error is dropped. The next one (0.005 <= 5 * 0.002) is redone, and its
    // error is the one reported.
    let results = vec![
        Some(redone(10, 1e-20)),
        Some(redone(20, 1e-20)),
        Some(Err(anyhow::anyhow!("ended"))),
        Some(Err(anyhow::anyhow!("reported"))),
    ];
    let error = x_replay(&matches, results, &lists, &mut heaps, &batches)
        .expect_err("the redone match's error");
    assert_eq!(error.to_string(), "reported");
    assert_eq!((heaps[0].n, heaps[1].n), (1, 1));
}

#[test]
fn the_replay_inserts_in_stream_order_and_counts_ended_results() {
    let (lists, heaps, batches) = two_query_batch();
    let matches = stream(&lists);
    // The third match is ended with a result (computed) or without one (not computed).
    for third in [Some(Err(anyhow::anyhow!("ended"))), None] {
        let computed = third.is_some();
        let mut heaps = heaps.clone();
        let results = vec![
            Some(redone(10, 1e-20)),
            Some(redone(20, 1e-20)),
            third,
            Some(redone(22, 0.005)),
        ];
        let ended =
            x_replay(&matches, results, &lists, &mut heaps, &batches).expect("no error is redone");
        assert_eq!(
            ended,
            if computed {
                replayed(0, 1, &[1])
            } else {
                replayed(1, 0, &[])
            }
        );
        // Query 1 keeps its best subject: 22 (E 0.005, above the cutoff) does not replace 20.
        let hit_lists = fill_results_from_compo_heaps(1, &mut heaps);
        let oids: Vec<Vec<u32>> = hit_lists
            .iter()
            .map(|hit_list| {
                let hit_list = hit_list.as_ref().expect("a hit list per query");
                hit_list.hsplist_array[..hit_list.hsplist_count]
                    .iter()
                    .map(|list| list.oid)
                    .collect()
            })
            .collect();
        assert_eq!(oids, [vec![10], vec![20]]);
    }
}

#[test]
#[should_panic(expected = "ended before its window, not ended at its replay")]
fn a_match_without_a_result_is_ended_at_its_replay() {
    let (lists, mut heaps, batches) = two_query_batch();
    let matches = stream(&lists);
    // The first match cannot be ended (no heap is filled): it needs its result.
    let results = vec![None, None, None, None];
    let _ = x_replay(&matches, results, &lists, &mut heaps, &batches);
}

#[test]
fn a_query_alone_in_its_batch_reads_only_its_own_heap() {
    // Query 0 has no match, so in a batch with query 1 the test is never true for query 1;
    // alone in its batch, query 1's second match is ended once its own heap is filled.
    let lists = vec![
        None,
        hit_list(vec![list(20, 1e-20, false), list(21, 1.0, false)]),
    ];
    for (batches, ended) in [
        (vec![0..2, 0..2], replayed(0, 0, &[])),
        (vec![0..1, 1..2], replayed(0, 1, &[1])),
    ] {
        let mut heaps = heaps_of(2, 1);
        let matches = stream(&lists);
        let results = vec![Some(redone(20, 1e-20)), Some(redone(21, 1.0))];
        assert_eq!(
            x_replay(&matches, results, &lists, &mut heaps, &batches).expect("no error"),
            ended,
            "{batches:?}"
        );
    }
}

/// The pre-check of a window rests on this: for a fixed E-value, once the early termination
/// test is true on a batch's heaps it stays true after any later `would_insert` / `insert`
/// (the replay's steps) on them.
#[test]
fn the_early_termination_test_stays_true_after_later_inserts() {
    // A fixed xorshift sequence.
    let mut state = 0x9e37_79b9_7f4a_7c15_u64;
    let mut next = move |bound: usize| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state % bound as u64) as usize
    };
    // E-values around the cutoff (0.002), its stretch (0.01) and the epsilon of the heap's
    // comparison (1e-180).
    let evalues = [
        0.0,
        1e-200,
        1e-181,
        1e-5,
        ECUTOFF,
        0.004,
        0.01,
        0.010_000_1,
        0.05,
        1.0,
        10.0,
    ];
    let mut true_then_tested = 0;
    for heap_threshold in [1, 2, 3, 5] {
        for _ in 0..200 {
            let mut heaps = heaps_of(3, heap_threshold);
            let mut ended = [false; 11];
            for _ in 0..40 {
                let q_idx = next(3);
                // Inserts within the cutoff twice as often, so that heaps fill.
                let evalue = evalues[if next(3) == 0 {
                    next(evalues.len())
                } else {
                    next(5)
                }];
                let (score, oid) = (next(4) as i32, next(6) as u32);
                if heaps[q_idx].would_insert(evalue, score, oid) {
                    heaps[q_idx].insert(list(oid, evalue, false), evalue, score, oid);
                }
                for (index, &tested) in evalues.iter().enumerate() {
                    let now = blast_compo_early_termination(tested, &heaps);
                    assert!(
                        now || !ended[index],
                        "E {tested}: true before an insert of E {evalue}, false after"
                    );
                    true_then_tested += usize::from(ended[index]);
                    ended[index] = now;
                }
            }
        }
    }
    assert!(true_then_tested > 1000, "{true_then_tested}");
}

#[test]
fn matches_may_be_ended_for_a_query_whose_batch_heaps_can_all_be_filled() {
    // Batch 0..3: queries 0 and 2 have two matches, query 1 one; query 3 has a list
    // without HSPs only, query 4 none.
    let lists = vec![
        hit_list(vec![list(1, 1e-9, false), list(2, 1.0, false)]),
        hit_list(vec![list(1, 1e-9, false)]),
        hit_list(vec![
            list(1, 1e-9, false),
            list(2, 1.0, false),
            list(3, 1.0, true),
        ]),
        hit_list(vec![list(1, 1.0, true)]),
        None,
    ];
    let batches = vec![0..3, 0..3, 0..3, 3..4, 4..5];
    let stretch = Some(COMPO_HEAP_EVALUE_STRETCH * ECUTOFF);
    // One subject to keep: only the last query of the batch.
    assert_eq!(
        x_end_above(&x_query_matches(&lists, &[]), &heaps_of(5, 1), &batches),
        [None, None, stretch, None, None]
    );
    // Two: query 1 cannot fill its heap.
    assert_eq!(
        x_end_above(&x_query_matches(&lists, &[]), &heaps_of(5, 2), &batches),
        [None; 5]
    );
    // A heap after a query that is filled already counts as such.
    let mut heaps = heaps_of(5, 1);
    heaps[2].insert(list(1, 1e-9, false), 1e-9, 40, 1);
    assert_eq!(
        x_end_above(&x_query_matches(&lists, &[]), &heaps, &batches),
        [None, stretch, stretch, None, None]
    );
}

/// The driver with results made from the preliminary lists (`redone` of each list's best E)
/// instead of the redo; returns the counts and how many matches were computed.
fn drive(
    lists: &mut [Option<BlastpHitList>],
    heaps: &mut [BlastCompoHeap],
    batches: &[Range<usize>],
    limits: XWindowLimits,
    window_seen: impl FnMut(&[Option<BlastpHitList>], &XWindow),
) -> (XCounts, usize) {
    drive_expanding(lists, heaps, batches, limits, None, window_seen)
}

/// `drive` with the compact lists of `expand` expanded ahead of the windows.
fn drive_expanding(
    lists: &mut [Option<BlastpHitList>],
    heaps: &mut [BlastCompoHeap],
    batches: &[Range<usize>],
    limits: XWindowLimits,
    expand: Option<&mut XExpandAhead<'_>>,
    mut window_seen: impl FnMut(&[Option<BlastpHitList>], &XWindow),
) -> (XCounts, usize) {
    let mut computed = 0;
    let counts = x_redo_windows(
        lists,
        heaps,
        batches,
        limits,
        false,
        expand,
        &mut |lists, _, window| {
            window_seen(lists, window);
            window
                .matches
                .iter()
                .zip(&window.ended)
                .map(|(&matched, &ended)| {
                    (!ended).then(|| {
                        computed += 1;
                        let local_match = x_local_match(lists, matched);
                        redone(local_match.oid, local_match.best_evalue)
                    })
                })
                .collect()
        },
    )
    .expect("no error");
    (counts, computed)
}

#[test]
fn a_query_alone_in_its_batch_computes_one_window_of_jobs_at_risk_once_its_heap_fills() {
    // `-max_target_seqs 1`: the first match fills the heap, and the test ends every later
    // one (E 1.0 > 5 x 0.002). The one-thread loop redoes only the first.
    const MATCHES: u32 = 2000;
    const THREADS: usize = 4;
    let mut lists = vec![hit_list(
        (0..MATCHES)
            .map(|oid| list(oid, if oid == 0 { 1e-20 } else { 1.0 }, false))
            .collect(),
    )];
    let mut heaps = heaps_of(1, 1);
    let limits = XWindowLimits {
        jobs: X_WINDOW_JOBS_PER_THREAD * THREADS,
        at_risk_jobs: THREADS,
    };
    let (counts, computed) = drive(&mut lists, &mut heaps, &alone(1), limits, |_, _| {});
    // The first window holds the jobs at risk that the cap allows (4 x 4 matches; only the
    // first match fits in the empty heap, and it shares its job with three at risk); the test
    // ends 15 of them at the replay, which returns the cap to four, and every later match
    // before its window, which leaves no job at risk.
    assert_eq!(computed, X_JOB_MATCHES * THREADS);
    assert_eq!(counts.wasted, X_JOB_MATCHES * THREADS - 1);
    assert_eq!(counts.skipped, MATCHES as usize - X_JOB_MATCHES * THREADS);
    assert_eq!(counts.skipped + counts.wasted, MATCHES as usize - 1);
    // The later windows hold their whole job budget (1,024 matches).
    assert_eq!((counts.windows, counts.split), (3, 2));
    assert_eq!((counts.matches, counts.released), (MATCHES as usize, 1));
    assert_eq!(heaps[0].n, 1);
}

/// The windows of a search where the heaps can be filled but are not: the cap on jobs at risk
/// grows while the test ends nothing.
#[test]
fn the_cap_on_jobs_at_risk_grows_while_the_test_ends_nothing() {
    // `-max_target_seqs 5`, a query alone in its batch with 4,000 matches whose results stay
    // above the cutoff (E 1.0): the heap keeps five lists but is never filled to its cutoff,
    // so the test ends nothing, though every match after the fifth may be ended.
    const MATCHES: u32 = 4000;
    const THREADS: usize = 4;
    let mut lists = vec![hit_list(
        (0..MATCHES).map(|oid| list(oid, 1.0, false)).collect(),
    )];
    let mut heaps = heaps_of(1, 5);
    let limits = XWindowLimits {
        jobs: X_WINDOW_JOBS_PER_THREAD * THREADS,
        at_risk_jobs: THREADS,
    };
    let mut at_risk = Vec::new();
    let (counts, computed) = drive(&mut lists, &mut heaps, &alone(1), limits, |_, window| {
        at_risk.push(window.at_risk_jobs())
    });
    // The first window: the first job (its matches fit in the empty heap) and four at risk;
    // then the cap doubles with each window up to the job budget (256): 1,000 jobs in nine
    // windows, where a cap fixed at four takes 250.
    assert_eq!(at_risk, [4, 8, 16, 32, 64, 128, 256, 256, 235]);
    assert_eq!((counts.windows, counts.split), (9, 8));
    assert_eq!(
        (counts.skipped, counts.wasted, computed),
        (0, 0, MATCHES as usize)
    );
    assert_eq!(heaps[0].n, 5);
}

/// Batches whose heaps never fill before one that fills: each query's cap on jobs at risk
/// grows with its own work only, so the batch that fills drops at most its own cap's results.
#[test]
fn a_batch_that_fills_drops_at_most_its_own_cap_after_batches_that_never_fill() {
    const THREADS: usize = 4;
    // `-max_target_seqs 1`, each query alone in its batch: queries 0 and 1 keep a list above
    // the cutoff (E 1.0), so their heaps are never filled; query 2's first match fills its
    // heap, and the test ends the others.
    let weak = |count: u32| hit_list((0..count).map(|oid| list(oid, 1.0, false)).collect());
    let mut lists = vec![
        weak(200),
        weak(200),
        hit_list(
            (0..400)
                .map(|oid| list(oid, if oid == 0 { 1e-20 } else { 1.0 }, false))
                .collect(),
        ),
    ];
    let mut heaps = heaps_of(3, 1);
    let limits = XWindowLimits {
        jobs: X_WINDOW_JOBS_PER_THREAD * THREADS,
        at_risk_jobs: THREADS,
    };
    let mut at_risk = Vec::new();
    let (counts, computed) = drive(&mut lists, &mut heaps, &alone(3), limits, |_, window| {
        at_risk.push(window.at_risk.clone())
    });
    // The caps of queries 0 and 1 double over their own windows (50 jobs each; the first
    // match of each fits in the empty heap), and a window holds jobs at risk of two queries
    // where the first one's end within its cap. Query 2 starts at four: its first window
    // computes 16 matches, the test drops 15 of them, and the later 384 are ended before
    // their window.
    assert_eq!(
        at_risk,
        [
            vec![(0, 4)],
            vec![(0, 8)],
            vec![(0, 16)],
            vec![(0, 22), (1, 4)],
            vec![(1, 8)],
            vec![(1, 16)],
            vec![(1, 22), (2, 4)],
            vec![],
        ]
    );
    assert_eq!((counts.wasted, counts.skipped, computed), (15, 384, 416));
    assert_eq!((counts.windows, counts.split), (8, 7));
}

#[test]
fn a_replayed_query_with_an_unfilled_heap_and_room_in_the_query_heap_leave_jobs_not_at_risk() {
    // A replayed query of the batch: query 0's one match is replayed when the window starts at
    // query 1 (eight matches, `-max_target_seqs 1`).
    let lists = vec![
        hit_list(vec![list(1, 1.0, false)]),
        hit_list((0..8).map(|oid| list(oid, 1.0, false)).collect()),
    ];
    let batches = vec![0..2, 0..2];
    let end_above = x_end_above(&x_query_matches(&lists, &[]), &heaps_of(2, 1), &batches);
    assert_eq!(end_above, [None, Some(COMPO_HEAP_EVALUE_STRETCH * ECUTOFF)]);
    let at_risk = |heaps: &[BlastCompoHeap], start| {
        window_from(&lists, heaps, &batches, &end_above, start, usize::MAX).at_risk_jobs()
    };
    // Query 0 kept a list above the cutoff: its heap is not filled and stays so, so no
    // match of query 1 can be ended.
    let mut heaps = heaps_of(2, 1);
    heaps[0].insert(list(1, 1.0, false), 1.0, 40, 1);
    assert_eq!(at_risk(&heaps, at(1, 0)), 0);
    // Within the cutoff: both jobs of query 1 are at risk (all its matches but the first).
    let mut heaps = heaps_of(2, 1);
    heaps[0].insert(list(1, 1e-9, false), 1e-9, 40, 1);
    assert_eq!(at_risk(&heaps, at(1, 0)), 2);
    // A window that starts at query 0 cannot count on its heap yet.
    assert_eq!(at_risk(&heaps_of(2, 1), at(0, 0)), 2);

    // Room in the query's heap: a query alone with 12 matches, `-max_target_seqs 8`, and a
    // window from its fifth match. The heap gains at most one list per match, so with `kept`
    // lists the window's first 8 - `kept` matches cannot be ended.
    let lists = vec![hit_list((0..12).map(|oid| list(oid, 1.0, false)).collect())];
    let end_above = x_end_above(&x_query_matches(&lists, &[]), &heaps_of(1, 8), &alone(1));
    for (kept, jobs_at_risk) in [(0, 0), (4, 1)] {
        let mut heaps = heaps_of(1, 8);
        for oid in 0..kept {
            heaps[0].insert(list(oid, 1.0, false), 1.0, 40, oid);
        }
        let window = window_from(&lists, &heaps, &alone(1), &end_above, at(0, 4), usize::MAX);
        assert_eq!(window.at_risk_jobs(), jobs_at_risk, "{kept} kept");
    }
}

#[test]
fn a_query_that_spans_windows_keeps_its_list_until_its_last_match_is_replayed() {
    let matches = |count: u32| hit_list((0..count).map(|oid| list(oid, 1e-9, false)).collect());
    let mut lists = vec![
        matches(2),
        None,
        matches(9),
        hit_list(vec![list(0, 1e-9, true)]),
        matches(1),
    ];
    let mut heaps = heaps_of(5, 1);
    let mut present = Vec::new();
    let (counts, computed) = drive(
        &mut lists,
        &mut heaps,
        &alone(5),
        XWindowLimits {
            jobs: 1,
            at_risk_jobs: 1,
        },
        |lists, window| {
            present.push((
                window.matches[0],
                lists.iter().map(Option::is_some).collect::<Vec<_>>(),
            ))
        },
    );
    // One job per window: query 2 spans three. A query's list is there until the window
    // that holds its last match is replayed.
    let (yes, no) = (true, false);
    assert_eq!(
        present,
        [
            (at(0, 0), vec![yes, no, yes, yes, yes]),
            (at(2, 0), vec![no, no, yes, yes, yes]),
            (at(2, 4), vec![no, no, yes, yes, yes]),
            (at(2, 8), vec![no, no, yes, yes, yes]),
            (at(4, 0), vec![no, no, no, no, yes]),
        ]
    );
    assert_eq!((counts.windows, counts.split, computed), (5, 2, 12));
    assert_eq!(counts.released, 4);
    assert!(lists.iter().all(Option::is_none));
    // Every match is within the cutoff: each heap keeps all of its query's.
    assert_eq!(
        heaps.iter().map(|heap| heap.n).collect::<Vec<_>>(),
        [2, 0, 9, 0, 1]
    );
}

/// The compact list of a preliminary hit list: the same fields,
/// subjects, best E-values and HSPs per list.
fn compact_of(list: &Option<BlastpHitList>) -> Option<BlastpHitList<XPrelimHsp>> {
    list.as_ref().map(|list| BlastpHitList {
        hsplist_count: list.hsplist_count,
        hsplist_max: list.hsplist_max,
        worst_evalue: list.worst_evalue,
        low_score: list.low_score,
        heapified: list.heapified,
        hsplist_array: list
            .hsplist_array
            .iter()
            .map(|hsp_list| BlastpHspList {
                oid: hsp_list.oid,
                query_index: hsp_list.query_index,
                hsps: hsp_list
                    .hsps
                    .iter()
                    .map(|hsp| XPrelimHsp {
                        hsp: BlastpPreliminaryHsp {
                            query_start: hsp.q_start as i32 - 1,
                            query_end: hsp.q_end as i32,
                            subject_start: hsp.s_start as i32 - 1,
                            subject_end: hsp.s_end as i32,
                            gapped_query_start: 0,
                            gapped_subject_start: 0,
                            raw_score: hsp.raw_score,
                            query_context: 0,
                            query_frame: 0,
                            subject_frame: 0,
                            query_length: hsp.query_length,
                            q_idx: 0,
                            s_idx: hsp_list.oid,
                        },
                        e_value: hsp.e_value,
                    })
                    .collect(),
                best_evalue: hsp_list.best_evalue,
            })
            .collect(),
        hsplist_current: list.hsplist_current,
        num_hits: list.num_hits,
    })
}

/// A window as the tests record it.
type Seen = (
    Vec<XMatch>,
    Vec<bool>,
    Vec<Range<usize>>,
    Option<XMatch>,
    Vec<(usize, usize)>,
);

fn seen(window: &XWindow) -> Seen {
    (
        window.matches.clone(),
        window.ended.clone(),
        window.jobs.clone(),
        window.next,
        window.at_risk.clone(),
    )
}

/// The driver on `lists` expanded before the windows and on their compact lists expanded ahead
/// of the windows: the same windows, counts, matches computed and heaps. Returns, per window of
/// the second run, which lists were there, and its counts.
fn same_windows_expanding_ahead(
    lists: &[Option<BlastpHitList>],
    heap_threshold: usize,
    batches: &[Range<usize>],
    limits: XWindowLimits,
) -> (Vec<Vec<bool>>, XCounts) {
    let count = lists.len();
    let mut full_lists = lists.to_vec();
    let mut full_heaps = heaps_of(count, heap_threshold);
    let mut full_windows = Vec::new();
    let (full_counts, full_computed) = drive(
        &mut full_lists,
        &mut full_heaps,
        batches,
        limits,
        |_, window| full_windows.push(seen(window)),
    );
    let mut compact: Vec<_> = lists.iter().map(compact_of).collect();
    let mut ahead_lists: Vec<Option<BlastpHitList>> = vec![None; count];
    let mut expand = XExpandAhead::new(&ahead_lists, &mut compact);
    let mut heaps = heaps_of(count, heap_threshold);
    let mut windows = Vec::new();
    let mut present = Vec::new();
    let (counts, computed) = drive_expanding(
        &mut ahead_lists,
        &mut heaps,
        batches,
        limits,
        Some(&mut expand),
        |lists, window| {
            windows.push(seen(window));
            present.push(lists.iter().map(Option::is_some).collect());
        },
    );
    assert_eq!(windows, full_windows);
    assert_eq!(
        (
            counts.windows,
            counts.split,
            counts.matches,
            counts.skipped,
            counts.wasted,
            counts.released,
            computed
        ),
        (
            full_counts.windows,
            full_counts.split,
            full_counts.matches,
            full_counts.skipped,
            full_counts.wasted,
            full_counts.released,
            full_computed
        )
    );
    let ns = |heaps: &[BlastCompoHeap]| heaps.iter().map(|heap| heap.n).collect::<Vec<_>>();
    assert_eq!(ns(&heaps), ns(&full_heaps));
    drop(expand);
    assert!(compact.iter().all(Option::is_none), "every list expanded");
    assert!(
        ahead_lists.iter().all(Option::is_none),
        "every list released"
    );
    (present, counts)
}

/// Before each window the compact lists are expanded in query order until the queries from
/// the window's start hold more jobs than a window, and a query's list is released after the
/// window that holds its last match, so the lists there are those of about one window.
#[test]
fn windows_expand_the_compact_lists_as_they_reach_them_and_choose_the_same_windows() {
    let matches = |count: u32| hit_list((0..count).map(|oid| list(oid, 1e-9, false)).collect());
    // Jobs per query: 1, 0, 3, 0 (a list without HSPs), 1, 2, 2.
    let lists = vec![
        matches(2),
        None,
        matches(9),
        hit_list(vec![list(0, 1e-9, true)]),
        matches(1),
        matches(6),
        matches(5),
    ];
    // Heaps that keep every list: no job is at risk, and each window holds two jobs.
    let limits = XWindowLimits {
        jobs: 2,
        at_risk_jobs: 1,
    };
    let (present, _) = same_windows_expanding_ahead(&lists, 100, &alone(7), limits);
    let (yes, no) = (true, false);
    assert_eq!(
        present,
        [
            // Queries 0-2 (four jobs) for the window of query 0 and query 2's first job.
            vec![yes, no, yes, no, no, no, no],
            // Query 2's last two jobs, with queries 3 and 4 (one job) ahead.
            vec![no, no, yes, yes, yes, no, no],
            vec![no, no, no, no, yes, yes, no],
            vec![no, no, no, no, no, yes, yes],
            vec![no, no, no, no, no, no, yes],
        ]
    );
}

/// The same windows where the early termination test ends matches and the caps on jobs at
/// risk end windows: queries 0 and 1 never fill their heaps, query 2's first match fills its
/// heap (`-max_target_seqs 1`), and with one batch of three the test also reads the heaps of
/// the queries before.
#[test]
fn windows_expanding_ahead_end_and_drop_the_matches_of_the_full_expansion() {
    const THREADS: usize = 4;
    let weak = |count: u32| hit_list((0..count).map(|oid| list(oid, 1.0, false)).collect());
    let strong = hit_list(
        (0..400)
            .map(|oid| list(oid, if oid == 0 { 1e-20 } else { 1.0 }, false))
            .collect(),
    );
    let lists = vec![weak(200), weak(200), strong];
    // Windows of 256 jobs (all 200 expanded before the first) and of 8 jobs, where each query's
    // list is expanded when the windows reach it and the caps on jobs at risk end windows.
    for (batches, limits) in [
        (
            alone(3),
            XWindowLimits {
                jobs: X_WINDOW_JOBS_PER_THREAD * THREADS,
                at_risk_jobs: THREADS,
            },
        ),
        (
            alone(3),
            XWindowLimits {
                jobs: 8,
                at_risk_jobs: 1,
            },
        ),
        (
            vec![0..3; 3],
            XWindowLimits {
                jobs: 8,
                at_risk_jobs: 1,
            },
        ),
    ] {
        let (present, counts) = same_windows_expanding_ahead(&lists, 1, &batches, limits);
        if limits.jobs == 8 && batches[2] == (2..3) {
            // Query 0's list alone at the first window; query 2's matches ended before
            // their window and at their replay.
            assert_eq!(present[0], [true, false, false]);
            assert!(counts.skipped > 0 && counts.wasted > 0, "{counts:?}");
        }
    }
}

#[test]
fn results_differ_by_any_field_or_by_error() {
    let a = redone(1, 1e-10);
    assert_eq!(
        x_result_difference(&a, &redone(1, 1e-10).map_err(|e| e.to_string())),
        None
    );
    let mut b = redone(1, 1e-10).expect("a result");
    b.hsp_list.as_mut().expect("a list").hsps[0].gapped_s_start = 1;
    assert!(x_result_difference(&a, &Ok(b)).is_some());
    assert!(x_result_difference(&a, &Err("failed".to_string())).is_some());
    assert_eq!(
        x_result_difference(&Err(anyhow::anyhow!("failed")), &Err("failed".to_string())),
        None
    );
}

/// A fixed xorshift sequence.
struct XorShift(u64);

impl XorShift {
    fn below(&mut self, bound: usize) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 % bound as u64) as usize
    }
}

#[test]
fn the_lpt_order_is_a_permutation_largest_first_with_ties_in_index_order() {
    assert!(x_lpt_order(&[]).is_empty());
    assert_eq!(x_lpt_order(&[3, 9, 3, 0, 9, 4]), [1, 4, 5, 0, 2, 3]);
    let mut rng = XorShift(0x2545_f491_4f6c_dd1d);
    for len in [1, 2, 7, 64, 300] {
        // Few distinct costs (many ties), and many.
        for distinct in [1, 3, 1000] {
            let costs: Vec<u64> = (0..len).map(|_| rng.below(distinct) as u64).collect();
            let order = x_lpt_order(&costs);
            let mut jobs = order.clone();
            jobs.sort_unstable();
            assert_eq!(jobs, (0..len).collect::<Vec<_>>(), "a permutation");
            for pair in order.windows(2) {
                let (a, b) = (pair[0], pair[1]);
                assert!(
                    costs[a] > costs[b] || (costs[a] == costs[b] && a < b),
                    "job {a} (cost {}) before job {b} (cost {})",
                    costs[a],
                    costs[b]
                );
            }
        }
    }
}

/// The scope runs every job once, with at most one scratch per task, and returns the results
/// by job index; also when the jobs run nested parallel work, where a thread that waits in it
/// may run another task's loop (no task waits for another one).
#[test]
fn the_lpt_scope_runs_every_job_once_and_returns_the_results_by_job_index() {
    for threads in [2, 3, 7] {
        crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
            for jobs in [1_usize, 2, 5, 40] {
                let costs: Vec<u64> = (0..jobs as u64).map(|job| (job * 7919) % 13).collect();
                let order = x_lpt_order(&costs);
                let tasks = threads.min(jobs);
                let (scratches, runs) = (AtomicUsize::new(0), AtomicUsize::new(0));
                let results = x_compute_window_lpt(
                    pool,
                    tasks,
                    &order,
                    || {
                        scratches.fetch_add(1, Ordering::Relaxed);
                        Vec::new()
                    },
                    |scratch: &mut Vec<usize>, job| {
                        runs.fetch_add(1, Ordering::Relaxed);
                        scratch.push(job);
                        let nested: u64 = (0..64_u64).into_par_iter().map(|x| x * job as u64).sum();
                        (job, nested, rayon::current_thread_index().is_some())
                    },
                );
                let expected: Vec<(usize, u64, bool)> = (0..jobs)
                    .map(|job| (job, 2016 * job as u64, true))
                    .collect();
                assert_eq!(results, expected, "threads {threads} jobs {jobs}");
                assert_eq!(runs.load(Ordering::Relaxed), jobs);
                let made = scratches.load(Ordering::Relaxed);
                assert!(
                    (1..=tasks).contains(&made),
                    "{made} scratches for {tasks} tasks"
                );
            }
            Ok(())
        })
        .expect("a pool");
    }
}

/// A job that panics: the panic comes out of the scope (the other tasks finish first, so
/// nothing hangs) and the same pool computes the next call correctly.
#[test]
fn a_panicking_job_propagates_and_the_pool_stays_usable() {
    for threads in [2, 3, 7] {
        crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
            let jobs = 40_usize;
            let order: Vec<usize> = (0..jobs).collect();
            let run = |panicking: bool| {
                x_compute_window_lpt(
                    pool,
                    threads.min(jobs),
                    &order,
                    || (),
                    |_: &mut (), job| {
                        assert!(!(panicking && (job == 5 || job == 30)), "job {job} panics");
                        job * 3
                    },
                )
            };
            let caught = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| run(true)));
            assert!(caught.is_err(), "threads {threads}: the panic was lost");
            let expected: Vec<usize> = (0..jobs).map(|job| job * 3).collect();
            assert_eq!(run(false), expected, "threads {threads}");
            Ok(())
        })
        .expect("a pool");
    }
}

const AMINO_ACIDS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

/// A redo's inputs made up from related proteins: random queries; subjects that hold a
/// mutated piece of one query between random flanks; per query a preliminary hit list with a
/// list for each related subject (one HSP over the piece, its gapped start in the middle), a
/// few lists of unrelated subjects, and one list without HSPs; one query without a list.
struct XRedoFixture {
    contexts: Vec<QueryContext>,
    search_spaces: Vec<SearchSpace>,
    subjects: Vec<EncodedProtein>,
    lists: Vec<Option<BlastpHitList>>,
    batches: Vec<Range<usize>>,
    evalue_adjust_contexts: Vec<crate::algorithm::blastp::kappa::EvalueAdjustContext>,
    gapped_params: KarlinParams,
    gapped_gumbel: BlastGumbelBlk,
    matrix_info: crate::core::composition_adjustment::adjust_scores::BlastMatrixInfo,
}

/// A list of subject `oid` for query `q_idx` (of `query_length` letters) with one HSP of
/// `length` letters from query offset `q0` and subject offset `s0` (0-based).
fn aligned_list(
    q_idx: usize,
    oid: usize,
    (q0, s0, length): (usize, usize, usize),
    query_length: usize,
    e_value: f64,
) -> BlastpHspList {
    let mut hsp = hsp(oid as u32, e_value);
    hsp.q_start = q0 + 1;
    hsp.q_end = q0 + length;
    hsp.s_start = s0 + 1;
    hsp.s_end = s0 + length;
    hsp.gapped_q_start = (q0 + length / 2) as i32;
    hsp.gapped_s_start = (s0 + length / 2) as i32;
    hsp.length = length;
    hsp.raw_score = 2 * length as i32;
    hsp.query_length = query_length;
    hsp.q_idx = q_idx as u32;
    BlastpHspList {
        oid: oid as u32,
        query_index: q_idx as u32,
        hsps: vec![hsp],
        best_evalue: e_value,
    }
}

fn redo_fixture() -> XRedoFixture {
    let spec = crate::config::ProteinScoringSpec::default();
    let gapped_params = crate::stats::tables::lookup_protein_params(&spec);
    let mut rng = XorShift(0x9e37_79b9_7f4a_7c15);
    let random = |rng: &mut XorShift, len: usize| -> Vec<u8> {
        (0..len).map(|_| AMINO_ACIDS[rng.below(20)]).collect()
    };
    let queries: Vec<Vec<u8>> = (0..6)
        .map(|_| {
            let len = 120 + rng.below(240);
            random(&mut rng, len)
        })
        .collect();
    let mut texts = Vec::new();
    let mut related = Vec::new();
    for oid in 0..48 {
        let q_idx = rng.below(queries.len());
        let query = &queries[q_idx];
        let length = 50 + rng.below(query.len() - 60);
        let q0 = rng.below(query.len() - length + 1);
        let left = rng.below(80);
        let mut subject = random(&mut rng, left);
        for &residue in &query[q0..q0 + length] {
            // About a quarter of the letters changed.
            subject.push(if rng.below(4) == 0 {
                AMINO_ACIDS[rng.below(20)]
            } else {
                residue
            });
        }
        let right = rng.below(80);
        subject.extend(random(&mut rng, right));
        related.push((q_idx, oid, (q0, left, length)));
        texts.push(subject);
    }
    let mut per_query: Vec<Vec<BlastpHspList>> = vec![Vec::new(); queries.len()];
    for &(q_idx, oid, piece) in &related {
        let e_value = 10_f64.powi(-(3 + rng.below(40) as i32));
        per_query[q_idx].push(aligned_list(
            q_idx,
            oid,
            piece,
            queries[q_idx].len(),
            e_value,
        ));
    }
    for (q_idx, lists) in per_query.iter_mut().enumerate() {
        let unrelated: Vec<usize> = related
            .iter()
            .filter(|&&(other, ..)| other != q_idx)
            .map(|&(_, oid, _)| oid)
            .collect();
        for &oid in unrelated.iter().take(3) {
            let length = 30;
            let q0 = rng.below(queries[q_idx].len() - length);
            let s0 = rng.below(texts[oid].len() - length);
            lists.push(aligned_list(
                q_idx,
                oid,
                (q0, s0, length),
                queries[q_idx].len(),
                0.5,
            ));
        }
        lists.push(BlastpHspList {
            oid: unrelated[3] as u32,
            query_index: q_idx as u32,
            hsps: Vec::new(),
            best_evalue: 1.0,
        });
        lists.sort_by(|a, b| a.best_evalue.total_cmp(&b.best_evalue));
    }
    let mut lists: Vec<Option<BlastpHitList>> = per_query.into_iter().map(hit_list).collect();
    lists[4] = None;
    let contexts: Vec<QueryContext> = queries
        .iter()
        .enumerate()
        .map(|(q_idx, query)| {
            let encoded = super::super::encoding::encode_protein_sequence(query);
            QueryContext {
                q_idx: q_idx as u32,
                frame: 0,
                aa_seq: encoded.aa_seq,
                aa_seq_nomask: None,
                aa_len: encoded.aa_len,
                orig_len: encoded.aa_len,
                frame_base: 0,
                is_valid: true,
                karlin_params: gapped_params,
            }
        })
        .collect();
    let subjects: Vec<EncodedProtein> = texts
        .iter()
        .map(|text| super::super::encoding::encode_protein_sequence(text))
        .collect();
    let db_letters: usize = subjects.iter().map(|subject| subject.aa_len).sum();
    let search_spaces: Vec<SearchSpace> = contexts
        .iter()
        .map(|ctx| {
            let length_adjustment = 20_i64;
            let effective_query_len = (ctx.aa_len as i64 - length_adjustment) as f64;
            let effective_db_len =
                (db_letters as i64 - subjects.len() as i64 * length_adjustment) as f64;
            SearchSpace {
                effective_space: effective_query_len * effective_db_len,
                length_adjustment,
            }
        })
        .collect();
    let ideal =
        crate::stats::compute_blosum62_ideal_karlin_params().expect("BLOSUM62 ideal parameters");
    let batches = vec![0..3, 0..3, 0..3, 3..6, 3..6, 3..6];
    let evalue_adjust_contexts = crate::algorithm::blastp::kappa::evalue_adjust_contexts(
        &contexts,
        &search_spaces,
        &batches,
        0,
        None,
    )
    .expect("no refusal without unified P-values");
    XRedoFixture {
        contexts,
        search_spaces,
        subjects,
        lists,
        batches,
        evalue_adjust_contexts,
        gapped_params,
        gapped_gumbel: crate::stats::spouge::lookup_protein_gumbel_params(&spec, db_letters as i64)
            .expect("BLOSUM62 11/1 Gumbel parameters"),
        matrix_info: crate::core::composition_adjustment::adjust_scores::build_matrix_info(
            ScoringMatrix::Blosum62,
            ideal.lambda / 32.0,
        )
        .expect("BLOSUM62 matrix information"),
    }
}

/// The redo parameters of blastp's defaults (composition matrix adjustment, scaled by 32), as
/// `run_resolved_in_pool` builds them.
fn redo_params(fixture: &XRedoFixture, q_idx: usize) -> BlastRedoAlignParams {
    use crate::core::composition_adjustment::redo_alignment::{
        BlastCompoAdjustMode, BlastCompoGappingParams,
    };
    const SCALE: f64 = 32.0;
    let scaled_lambda = fixture.gapped_params.lambda / SCALE;
    BlastRedoAlignParams {
        matrix_info: fixture.matrix_info.clone(),
        gapping_params: BlastCompoGappingParams {
            gap_open: 11 * SCALE as i32,
            gap_extend: SCALE as i32,
            x_dropoff: (25.0 * std::f64::consts::LN_2 / scaled_lambda) as i32,
            context: std::cell::Cell::new(None),
        },
        compo_adjust_mode: BlastCompoAdjustMode::CompositionMatrixAdjust,
        alphsize: crate::utils::matrix::BLASTAA_SIZE as i32,
        composition_test_index: 0,
        log_k: fixture.gapped_params.k.glibc_ln(),
        score_divisor: SCALE,
        smith_waterman: false,
        near_identical_cutoff: 1.74 * std::f64::consts::LN_2 / scaled_lambda,
        position_based: false,
        ccat_query_length: fixture.contexts[q_idx].aa_len as i32,
        query_is_translated: false,
        subject_is_translated: false,
        cutoff_score: 1,
        do_link_hsps: false,
    }
}

fn query_raw(ctx: &QueryContext) -> &[u8] {
    &ctx.aa_seq[1..ctx.aa_seq.len() - 1]
}

/// `x_compute_window` with the switch gives the results of the `par_iter` path, match for
/// match (error or not, scores, E-value bits, every HSP field), on windows from random starts
/// with random job budgets and matches ended at random, at 1, 2, 3 and 7 threads, with the
/// search's SEG store, and when the parameters of one query fail.
#[test]
fn windows_computed_with_jobs_largest_first_equal_the_par_iter_windows() {
    let fixture = redo_fixture();
    let count = fixture.lists.len();
    let heaps = heaps_of(count, 500);
    let all = stream(&fixture.lists);
    let mut rng = XorShift(0x5851_f42d_4c95_7f2d);
    let windows: Vec<XWindow> = (0..10)
        .map(|index| {
            let start = all[rng.below(all.len())];
            let mut window = window_from(
                &fixture.lists,
                &heaps,
                &fixture.batches,
                &vec![None; count],
                start,
                2 + rng.below(10),
            );
            // None ended, about one in four, about one in two.
            let share = index % 3;
            for ended in &mut window.ended {
                *ended = share > 0 && rng.below(4) < share;
            }
            window
        })
        .collect();
    let failing = all[all.len() / 2].q_idx;
    let (mut compared, mut with_list, mut errors, mut lpt_windows) = (0, 0, 0, 0);
    for threads in [1, 2, 3, 7] {
        crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
            for fail in [None, Some(failing)] {
                let build = |q_idx: usize| -> Result<BlastRedoAlignParams> {
                    if Some(q_idx) == fail {
                        anyhow::bail!("no parameters for query {q_idx}");
                    }
                    Ok(redo_params(&fixture, q_idx))
                };
                let search = XRedoParSearch {
                    pool,
                    threads,
                    contexts: &fixture.contexts,
                    query_nomask: query_raw,
                    search_spaces: &fixture.search_spaces,
                    subjects: crate::algorithm::blastp::x_redo_subjects::XRedoSubjects::Encoded(
                        &fixture.subjects,
                    ),
                    query_batches: &fixture.batches,
                    evalue_adjust_contexts: &fixture.evalue_adjust_contexts,
                    matrix: ScoringMatrix::Blosum62,
                    gapped_params: &fixture.gapped_params,
                    gapped_gumbel: Some(&fixture.gapped_gumbel),
                    evalue: 10.0,
                    build_redo_align_params: &build,
                    timing: false,
                    engine_positives: false,
                };
                let shared = || Arc::new(XSharedSubjectRanges::new());
                for window in &windows {
                    let (reference, ran) =
                        x_compute_window(&search, &fixture.lists, &heaps, window, &shared(), false);
                    assert!(!ran, "the par_iter path");
                    let (lpt, ran) =
                        x_compute_window(&search, &fixture.lists, &heaps, window, &shared(), true);
                    let acts = threads > 1
                        && window.jobs.len() > 1
                        && !window.ended.iter().all(|&ended| ended);
                    assert_eq!(ran, acts, "threads {threads}: {window:?}");
                    lpt_windows += usize::from(ran);
                    assert_eq!(lpt.len(), reference.len());
                    for (position, (computed, reference)) in lpt.iter().zip(reference).enumerate() {
                        let at = (threads, fail, window.matches[position]);
                        match (computed, reference) {
                            (None, None) => {}
                            (Some(computed), Some(reference)) => {
                                with_list += usize::from(
                                    reference
                                        .as_ref()
                                        .is_ok_and(|result| result.hsp_list.is_some()),
                                );
                                errors += usize::from(reference.is_err());
                                let reference = reference.map_err(|error| format!("{error:#}"));
                                assert_eq!(
                                    x_result_difference(computed, &reference),
                                    None,
                                    "{at:?}"
                                );
                                compared += 1;
                            }
                            (computed, reference) => panic!(
                                "{at:?}: computed {} / {}",
                                computed.is_some(),
                                reference.is_some()
                            ),
                        }
                    }
                }
            }
            Ok(())
        })
        .expect("a pool");
    }
    // Real work on every side: matches computed, alignments kept, the failing query's errors,
    // and windows the switch computed.
    assert!(
        compared > 200 && with_list > 50 && errors > 0 && lpt_windows > 30,
        "compared {compared} with a list {with_list} errors {errors} lpt windows {lpt_windows}"
    );
}

/// EXPERIMENT (LOSAT_X_BPSEGWIN): windows computed with the search's SEG store emptied after
/// each one (`XSharedSubjectRanges::x_clear`) give the results of the same windows computed
/// with one store kept for the whole search, match for match (error or not, scores, E-value
/// bits, every HSP field), on windows from random starts with random job budgets and matches
/// ended at random, at 1, 2 and 7 threads, and when the parameters of one query fail; the
/// later windows mask ranges again that an earlier window had masked.
#[test]
fn windows_with_the_seg_store_emptied_after_each_equal_windows_with_one_store() {
    let fixture = redo_fixture();
    let count = fixture.lists.len();
    let heaps = heaps_of(count, 500);
    let all = stream(&fixture.lists);
    let mut rng = XorShift(0x2545_f491_4f6c_dd1d);
    let windows: Vec<XWindow> = (0..10)
        .map(|index| {
            let start = all[rng.below(all.len())];
            let mut window = window_from(
                &fixture.lists,
                &heaps,
                &fixture.batches,
                &vec![None; count],
                start,
                2 + rng.below(10),
            );
            // None ended, about one in four, about one in two.
            let share = index % 3;
            for ended in &mut window.ended {
                *ended = share > 0 && rng.below(4) < share;
            }
            window
        })
        .collect();
    let failing = all[all.len() / 2].q_idx;
    let (mut compared, mut with_list, mut errors, mut masked_again) = (0, 0, 0, 0);
    for threads in [1, 2, 7] {
        crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
            for fail in [None, Some(failing)] {
                let build = |q_idx: usize| -> Result<BlastRedoAlignParams> {
                    if Some(q_idx) == fail {
                        anyhow::bail!("no parameters for query {q_idx}");
                    }
                    Ok(redo_params(&fixture, q_idx))
                };
                let search = XRedoParSearch {
                    pool,
                    threads,
                    contexts: &fixture.contexts,
                    query_nomask: query_raw,
                    search_spaces: &fixture.search_spaces,
                    subjects: crate::algorithm::blastp::x_redo_subjects::XRedoSubjects::Encoded(
                        &fixture.subjects,
                    ),
                    query_batches: &fixture.batches,
                    evalue_adjust_contexts: &fixture.evalue_adjust_contexts,
                    matrix: ScoringMatrix::Blosum62,
                    gapped_params: &fixture.gapped_params,
                    gapped_gumbel: Some(&fixture.gapped_gumbel),
                    evalue: 10.0,
                    build_redo_align_params: &build,
                    timing: false,
                    engine_positives: false,
                };
                let kept = Arc::new(XSharedSubjectRanges::new());
                let emptied = Arc::new(XSharedSubjectRanges::new());
                let mut window_entries = 0;
                for window in &windows {
                    let (reference, _) =
                        x_compute_window(&search, &fixture.lists, &heaps, window, &kept, false);
                    let (computed, _) =
                        x_compute_window(&search, &fixture.lists, &heaps, window, &emptied, false);
                    window_entries += emptied.x_len();
                    emptied.x_clear();
                    assert_eq!(emptied.x_len(), 0);
                    assert_eq!(computed.len(), reference.len());
                    for (position, (computed, reference)) in
                        computed.iter().zip(reference).enumerate()
                    {
                        let at = (threads, fail, window.matches[position]);
                        match (computed, reference) {
                            (None, None) => {}
                            (Some(computed), Some(reference)) => {
                                with_list += usize::from(
                                    reference
                                        .as_ref()
                                        .is_ok_and(|result| result.hsp_list.is_some()),
                                );
                                errors += usize::from(reference.is_err());
                                let reference = reference.map_err(|error| format!("{error:#}"));
                                assert_eq!(
                                    x_result_difference(computed, &reference),
                                    None,
                                    "{at:?}"
                                );
                                compared += 1;
                            }
                            (computed, reference) => panic!(
                                "{at:?}: computed {} / {}",
                                computed.is_some(),
                                reference.is_some()
                            ),
                        }
                    }
                }
                // The entries the emptied store made beyond the kept store's: ranges masked
                // again by a later window.
                masked_again += window_entries - kept.x_len();
            }
            Ok(())
        })
        .expect("a pool");
    }
    assert!(
        compared > 150 && with_list > 40 && errors > 0 && masked_again > 0,
        "compared {compared} with a list {with_list} errors {errors} masked again {masked_again}"
    );
}
