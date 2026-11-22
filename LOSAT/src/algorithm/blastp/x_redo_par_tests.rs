use super::*;
use crate::algorithm::blastp::encoding::EncodedProtein;
use crate::core::composition_adjustment::smith_waterman_match::SwRedoEndless;
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

/// A computed match without a deferred heap test.
fn done(result: Result<BlastpPostprocessResult>) -> Option<XComputed> {
    Some(result.into())
}

/// A computed match whose redo failed with `message`.
fn failed(message: &str) -> Option<XComputed> {
    done(Err(anyhow::anyhow!("{message}")))
}

/// The text of an error of the redo: the batch that the error of a batch's search names
/// (`HeaderReadError`, as `redo_run::sw_batch_failed` and `redo_run::sw_endless` make it), and the
/// message.
fn error_text(result: &Result<()>) -> Option<String> {
    result.as_ref().err().map(|error| {
        match error.downcast_ref::<crate::algorithm::common::database::HeaderReadError>() {
            Some(batch) => format!("batch {}: {:#}", batch.query, batch.error),
            None => format!("{error:#}"),
        }
    })
}

/// The text of the error of a redo that does not end, in the batch that starts at `batch_start`.
fn endless_in(batch_start: usize) -> Option<String> {
    error_text(&Err(sw_endless(batch_start)))
}

/// The text of the error of the batch that starts at `batch_start` and fails at its end, its first
/// query's heap holding a list (`holds`) or none.
fn batch_failed_text(batch_start: usize, holds: bool) -> Option<String> {
    let mut first_query = BlastCompoHeap::new(1, ECUTOFF);
    if holds {
        first_query.insert(list(0, 1e-9, false), 1e-9, 40, 0);
    }
    error_text(&Err(sw_batch_failed(batch_start, &first_query)))
}

/// The replay of the X-drop redo's window.
fn replay(
    matches: &[XMatch],
    results: XWindowResults,
    lists: &[Option<BlastpHitList>],
    heaps: &mut [BlastCompoHeap],
    batches: &[Range<usize>],
) -> Result<XReplayed> {
    x_replay::<false>(
        matches,
        results,
        lists,
        heaps,
        batches,
        &mut XSwBatch::new(&[]),
        &mut None,
    )
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
    x_window::<false>(lists, heaps, batches, end_above, &caps, None, start, jobs)
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
        turned: 0,
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
        done(redone(10, 1e-20)),
        done(redone(20, 1e-20)),
        failed("ended"),
        failed("reported"),
    ];
    let error = replay(&matches, results, &lists, &mut heaps, &batches)
        .expect_err("the redone match's error");
    assert_eq!(error.to_string(), "reported");
    assert_eq!((heaps[0].n, heaps[1].n), (1, 1));
}

#[test]
fn the_replay_inserts_in_stream_order_and_counts_ended_results() {
    let (lists, heaps, batches) = two_query_batch();
    let matches = stream(&lists);
    // The third match is ended with a result (computed) or without one (not computed).
    for third in [failed("ended"), None] {
        let computed = third.is_some();
        let mut heaps = heaps.clone();
        let results = vec![
            done(redone(10, 1e-20)),
            done(redone(20, 1e-20)),
            third,
            done(redone(22, 0.005)),
        ];
        let ended =
            replay(&matches, results, &lists, &mut heaps, &batches).expect("no error is redone");
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
    let _ = replay(&matches, results, &lists, &mut heaps, &batches);
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
        let results = vec![done(redone(20, 1e-20)), done(redone(21, 1.0))];
        assert_eq!(
            replay(&matches, results, &lists, &mut heaps, &batches).expect("no error"),
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
    let counts = x_redo_windows::<false>(
        lists,
        heaps,
        batches,
        limits,
        false,
        expand,
        &[],
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
                        redone(local_match.oid, local_match.best_evalue).into()
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
                    sw_excluded: &[],
                };
                let shared = || Arc::new(XSharedSubjectRanges::new());
                for window in &windows {
                    let (reference, ran) = x_compute_window::<_, false>(
                        &search,
                        &fixture.lists,
                        &heaps,
                        window,
                        &shared(),
                        false,
                    );
                    assert!(!ran, "the par_iter path");
                    let (lpt, ran) = x_compute_window::<_, false>(
                        &search,
                        &fixture.lists,
                        &heaps,
                        window,
                        &shared(),
                        true,
                    );
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
                                let reference = reference.result;
                                with_list += usize::from(
                                    reference
                                        .as_ref()
                                        .is_ok_and(|result| result.hsp_list.is_some()),
                                );
                                errors += usize::from(reference.is_err());
                                let reference = reference.map_err(|error| format!("{error:#}"));
                                assert_eq!(
                                    x_result_difference(&computed.result, &reference),
                                    None,
                                    "{at:?}"
                                );
                                assert!(computed.deferred.is_none(), "{at:?}");
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
                    sw_excluded: &[],
                };
                let kept = Arc::new(XSharedSubjectRanges::new());
                let emptied = Arc::new(XSharedSubjectRanges::new());
                let mut window_entries = 0;
                for window in &windows {
                    let (reference, _) = x_compute_window::<_, false>(
                        &search,
                        &fixture.lists,
                        &heaps,
                        window,
                        &kept,
                        false,
                    );
                    let (computed, _) = x_compute_window::<_, false>(
                        &search,
                        &fixture.lists,
                        &heaps,
                        window,
                        &emptied,
                        false,
                    );
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
                                let reference = reference.result;
                                with_list += usize::from(
                                    reference
                                        .as_ref()
                                        .is_ok_and(|result| result.hsp_list.is_some()),
                                );
                                errors += usize::from(reference.is_err());
                                let reference = reference.map_err(|error| format!("{error:#}"));
                                assert_eq!(
                                    x_result_difference(&computed.result, &reference),
                                    None,
                                    "{at:?}"
                                );
                                assert!(computed.deferred.is_none(), "{at:?}");
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

// ---- The Smith-Waterman redo in windows (`-use_sw_tback`, item X4-mt). ----

/// What the Smith-Waterman redo of a synthetic match gives: its first round's E-value and
/// score (its heap test, when the test is asked) and the outcome after a TRUE test.
#[derive(Clone, Copy, Debug, PartialEq)]
enum SwKind {
    /// The first round's score is not significant: no test, no alignment.
    NotAsked,
    /// Alignments whose best E-value and score are the first round's.
    Kept,
    /// No alignment left after the redo (all above `-evalue`).
    Purged,
    /// A failed redo (`SwRedoFailed`).
    Failed,
    /// A redo that does not end (`SwRedoEndless`).
    Endless,
}

/// A synthetic match of the Smith-Waterman redo: the E-value and score of its first round, which
/// its heap test asks, and those of its result after a TRUE answer (its best alignment's, which
/// the heap's insert test asks), drawn apart as the redo's rounds make them.
#[derive(Clone, Copy, Debug)]
struct SwCase {
    evalue: f64,
    score: i32,
    final_evalue: f64,
    final_score: i32,
    kind: SwKind,
}

fn no_alignment() -> Result<BlastpPostprocessResult> {
    Ok(BlastpPostprocessResult {
        hsp_list: None,
        best_score: 0,
        best_evalue: f64::MAX,
    })
}

/// A job's computation of a synthetic match: the heap test answered as `test` says (as
/// `x_redo_match_sw` answers it), the outcome of its case after a TRUE answer.
fn sw_job_match(test: XSwTest<'_>, oid: u32, case: SwCase) -> XComputed {
    let mut deferred = None;
    let answer = |deferred: &mut Option<(f64, i32)>| match test {
        XSwTest::Heap(heap, before) => {
            x_sw_heap_test(heap, before, (case.evalue, case.score, oid), deferred)
        }
        XSwTest::Continued(asked) => {
            assert_eq!(asked, (case.evalue, case.score), "the continued test");
            true
        }
    };
    let result = if case.kind == SwKind::NotAsked || !answer(&mut deferred) {
        no_alignment()
    } else {
        match case.kind {
            SwKind::Kept => Ok(BlastpPostprocessResult {
                hsp_list: Some(list(oid, case.final_evalue, false)),
                best_score: case.final_score,
                best_evalue: case.final_evalue,
            }),
            SwKind::Purged => no_alignment(),
            SwKind::Failed => Err(SwRedoFailed.into()),
            SwKind::Endless => Err(SwRedoEndless.into()),
            SwKind::NotAsked => unreachable!(),
        }
    };
    XComputed {
        result,
        deferred,
        first_round_s: 0.0,
    }
}

/// The windows of the Smith-Waterman redo on synthetic matches (`cases[q][position]`), each
/// job's matches computed by `sw_job_match`; returns the redo's result and the matches
/// computed.
fn sw_drive(
    lists: &mut [Option<BlastpHitList>],
    cases: &[Vec<SwCase>],
    heaps: &mut [BlastCompoHeap],
    batches: &[Range<usize>],
    limits: XWindowLimits,
) -> (Result<XCounts>, usize) {
    sw_drive_excluded(lists, cases, heaps, batches, limits, &[], |_| {})
}

/// `sw_drive` with the marks of the lists that the fetch excludes (`excluded`, by query and
/// place), and each window given to `window_seen` before its jobs run.
fn sw_drive_excluded(
    lists: &mut [Option<BlastpHitList>],
    cases: &[Vec<SwCase>],
    heaps: &mut [BlastCompoHeap],
    batches: &[Range<usize>],
    limits: XWindowLimits,
    excluded: &[FetchExcludedMatch],
    mut window_seen: impl FnMut(&XWindow),
) -> (Result<XCounts>, usize) {
    sw_drive_expanding(
        lists,
        cases,
        heaps,
        batches,
        limits,
        None,
        excluded,
        |_, window| window_seen(window),
    )
}

/// `sw_drive_excluded` with the compact lists of `expand` expanded ahead of the windows, each
/// window given to `window_seen` with the lists as they are while its jobs run.
#[allow(clippy::too_many_arguments)]
fn sw_drive_expanding(
    lists: &mut [Option<BlastpHitList>],
    cases: &[Vec<SwCase>],
    heaps: &mut [BlastCompoHeap],
    batches: &[Range<usize>],
    limits: XWindowLimits,
    expand: Option<&mut XExpandAhead<'_>>,
    excluded: &[FetchExcludedMatch],
    mut window_seen: impl FnMut(&[Option<BlastpHitList>], &XWindow),
) -> (Result<XCounts>, usize) {
    let mut computed = 0;
    let counts = x_redo_windows::<true>(
        lists,
        heaps,
        batches,
        limits,
        false,
        expand,
        excluded,
        &mut |lists, heaps, window| {
            window_seen(lists, window);
            window
                .jobs
                .iter()
                .flat_map(|job| job.clone())
                .map(|index| {
                    let matched = window.matches[index];
                    (!window.ended[index]).then(|| {
                        computed += 1;
                        let first = window
                            .matches
                            .partition_point(|other| other.q_idx < matched.q_idx);
                        let test =
                            match window.continued.binary_search_by_key(&index, |&(at, _)| at) {
                                Ok(found) => XSwTest::Continued(window.continued[found].1),
                                Err(_) => XSwTest::Heap(&heaps[matched.q_idx], index - first),
                            };
                        sw_job_match(
                            test,
                            x_local_match(lists, matched).oid,
                            cases[matched.q_idx][matched.position],
                        )
                    })
                })
                .collect()
        },
    );
    (counts, computed)
}

/// NCBI's one-thread Smith-Waterman redo of the synthetic matches (the steps of
/// `redo_run::one_thread_redo`): the early termination test, the heap test (a FALSE ends the
/// match with status 0), the outcome, the insert, the batch's status and its check at the
/// batch's last query.
fn sw_one_thread(
    lists: &[Option<BlastpHitList>],
    cases: &[Vec<SwCase>],
    heaps: &mut [BlastCompoHeap],
    batches: &[Range<usize>],
) -> Result<()> {
    sw_one_thread_excluded(lists, cases, heaps, batches, &[])
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3540
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3560-3575
/// `sw_one_thread` with the lists that the fetch excludes (`excluded`, by query and place) in
/// the stream as NCBI's loop reaches them: each one in its place, its early termination test on
/// the heaps as they are there, and, where the test does not hold, the status back to 0. This is
/// NCBI's order itself, not the reading after each query's loop that `one_thread_redo` and the
/// windows share (`redo_run::query_fetch_excluded_reached`).
fn sw_one_thread_excluded(
    lists: &[Option<BlastpHitList>],
    cases: &[Vec<SwCase>],
    heaps: &mut [BlastCompoHeap],
    batches: &[Range<usize>],
    excluded: &[FetchExcludedMatch],
) -> Result<()> {
    let mut failed = false;
    for q_idx in 0..lists.len() {
        let batch = batches[q_idx].clone();
        if batch.start == q_idx {
            failed = false;
        }
        let mut marks = excluded
            .iter()
            .filter(|mark| mark.query == q_idx)
            .peekable();
        // Whether NCBI's redo reaches an excluded list there (its early termination test does
        // not hold), so that the excluded fetch sets the status back to 0.
        let reached = |mark: &FetchExcludedMatch, heaps: &[BlastCompoHeap]| {
            !blast_compo_early_termination(mark.best_evalue, &heaps[batch.clone()])
        };
        if let Some(hit_list) = lists[q_idx].as_ref() {
            for (position, local_match) in hit_list.hsplist_array[..hit_list.hsplist_count]
                .iter()
                .enumerate()
            {
                // The excluded lists before this one (`position` kept lists before each).
                while let Some(mark) = marks.next_if(|mark| mark.position <= position) {
                    failed &= !reached(mark, heaps);
                }
                if local_match.hsps.is_empty()
                    || blast_compo_early_termination(local_match.best_evalue, &heaps[batch.clone()])
                {
                    continue;
                }
                let (case, oid) = (cases[q_idx][position], local_match.oid);
                if case.kind == SwKind::NotAsked
                    || !heaps[q_idx].would_insert(case.evalue, case.score, oid)
                {
                    failed = false;
                    continue;
                }
                match case.kind {
                    SwKind::Kept => {
                        failed = false;
                        let (evalue, score) = (case.final_evalue, case.final_score);
                        if heaps[q_idx].would_insert(evalue, score, oid) {
                            heaps[q_idx].insert(list(oid, evalue, false), evalue, score, oid);
                        }
                    }
                    SwKind::Purged => failed = false,
                    SwKind::Failed => failed = true,
                    SwKind::Endless => return Err(sw_endless(batch.start)),
                    SwKind::NotAsked => unreachable!(),
                }
            }
        }
        // The excluded lists after the query's last list.
        for mark in marks {
            failed &= !reached(mark, heaps);
        }
        if failed && q_idx + 1 == batch.end {
            return Err(sw_batch_failed(batch.start, &heaps[batch.start]));
        }
    }
    Ok(())
}

/// Each heap as its pops give it (the subjects, worst first), with its count and worst E.
fn heap_contents(heaps: &[BlastCompoHeap]) -> Vec<(usize, u64, Vec<u32>)> {
    heaps
        .iter()
        .map(|heap| {
            let mut heap = heap.clone();
            let (n, worst) = (heap.n, heap.worst_evalue.to_bits());
            let mut oids = Vec::new();
            while let Some(list) = heap.pop() {
                oids.push(list.oid);
            }
            (n, worst, oids)
        })
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:251-276
#[test]
fn the_job_answers_the_heap_test_where_it_is_certain_and_defers_the_rest() {
    let mut heap = BlastCompoHeap::new(2, ECUTOFF);
    let mut deferred = None;
    let mut test = |heap: &BlastCompoHeap, before, evalue, score, oid| {
        deferred = None;
        let answer = x_sw_heap_test(heap, before, (evalue, score, oid), &mut deferred);
        (answer, deferred)
    };
    // Room for the matches before: TRUE for certain, whatever the E-value.
    assert_eq!(test(&heap, 1, 5.0, 30, 9), (true, None));
    // Within the cutoff: TRUE for certain, however full the heap may be.
    assert_eq!(test(&heap, 7, ECUTOFF, 30, 9), (true, None));
    // The heap may be full by then and the empty heap's test says TRUE: deferred (the job's
    // redo stops at the test, as on FALSE).
    assert_eq!(test(&heap, 2, 5.0, 30, 9), (false, Some((5.0, 30))));
    heap.insert(list(1, 0.5, false), 0.5, 40, 1);
    heap.insert(list(2, 0.01, false), 0.01, 50, 2);
    // Full above the cutoff: better than the root (E 0.5) is deferred, worse is FALSE.
    assert_eq!(test(&heap, 0, 0.1, 30, 9), (false, Some((0.1, 30))));
    assert_eq!(test(&heap, 0, 0.6, 30, 9), (false, None));
    // A tie of E-values: the higher score is better, the lower is not.
    assert_eq!(test(&heap, 0, 0.5, 41, 9), (false, Some((0.5, 41))));
    assert_eq!(test(&heap, 0, 0.5, 39, 9), (false, None));
    // Without side effect: the array was not turned into a heap.
    assert!(heap.heap_array.is_none());
}

#[test]
fn smith_waterman_jobs_are_at_risk_only_as_x_drop_jobs_are() {
    // One query, six matches, two lists to keep, no match the early termination test may
    // end. The Smith-Waterman jobs' tests of the third to sixth matches are deferred (the heap
    // may hold two by then), which wastes nothing: neither redo has a job at risk, and a cap
    // of one job at risk does not end the window.
    let lists = vec![hit_list((0..6).map(|oid| list(oid, 1e-9, false)).collect())];
    let heaps = heaps_of(1, 2);
    let none = [None];
    for caps in [vec![usize::MAX], vec![1]] {
        let xdrop = x_window::<false>(&lists, &heaps, &alone(1), &none, &caps, None, at(0, 0), 8);
        let sw = x_window::<true>(&lists, &heaps, &alone(1), &none, &caps, None, at(0, 0), 8);
        assert_eq!((sw.jobs.len(), sw.at_risk_jobs()), (2, 0));
        assert_eq!(sw, xdrop);
    }
}

/// A match whose first round and result have the E-value `evalue` and the score 40.
fn sw_case(evalue: f64, kind: SwKind) -> SwCase {
    SwCase {
        evalue,
        score: 40,
        final_evalue: evalue,
        final_score: 40,
        kind,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1463-1476
/// A heap test that the job defers and the replay answers FALSE ends the match as NCBI ends it
/// (no alignment, status 0): the job's redo stopped at the test, so the failed or endless redo
/// after it is never computed; where the test is TRUE for certain, the same redos fail the
/// batch and end the search.
#[test]
fn a_deferred_heap_test_that_turns_false_masks_a_failed_and_an_endless_redo() {
    let lists = || vec![hit_list((0..3).map(|oid| list(oid, 1e-9, false)).collect())];
    let cases = vec![vec![
        sw_case(1e-20, SwKind::Kept),
        sw_case(0.5, SwKind::Failed),
        sw_case(0.5, SwKind::Endless),
    ]];
    let one_window = XWindowLimits {
        jobs: 8,
        at_risk_jobs: 8,
    };
    for (heap_threshold, expected) in [(1, None), (3, endless_in(0))] {
        let mut heaps = heaps_of(1, heap_threshold);
        let (counts, computed) = sw_drive(&mut lists(), &cases, &mut heaps, &alone(1), one_window);
        assert_eq!(computed, 3);
        let mut model = heaps_of(1, heap_threshold);
        let reference = sw_one_thread(&lists(), &cases, &mut model, &alone(1));
        assert_eq!(
            error_text(&reference),
            expected,
            "threshold {heap_threshold}"
        );
        match counts {
            Ok(counts) => {
                assert_eq!(expected, None);
                // Both later tests were deferred and turned FALSE.
                assert_eq!((counts.turned, counts.continued), (Some(2), 0));
            }
            Err(error) => assert_eq!(error_text(&Err(error)), expected),
        }
        assert_eq!(heap_contents(&heaps), heap_contents(&model));
    }
}

/// A heap test that the job defers and the replay answers TRUE stops the replay at its match:
/// the next window's first job finishes the match (a continuation) and the carried matches
/// after it are replayed with their results, before the jobs that follow them (which count
/// them as the query's matches before theirs). Every match is computed once, and each
/// continued one once more.
#[test]
fn a_deferred_heap_test_that_turns_true_is_finished_by_a_continuation_in_the_next_window() {
    // One query, one list to keep. The first match's redo asks no test (no insert), so the
    // second's deferred test is TRUE; the third is better than the second (TRUE); the others
    // are worse (FALSE), the last a failed redo.
    let lists = || {
        vec![hit_list(
            (0..10).map(|oid| list(oid, 1e-9, false)).collect(),
        )]
    };
    let mut query = vec![
        sw_case(0.5, SwKind::NotAsked),
        sw_case(0.5, SwKind::Kept),
        sw_case(0.4, SwKind::Kept),
    ];
    query.extend((3..9).map(|_| sw_case(0.6, SwKind::Kept)));
    query.push(sw_case(0.7, SwKind::Failed));
    let cases = vec![query];
    // One job per window: the replay stops at the second match, then (in the second window,
    // the carry and the next job) at the third; the third window holds the carried matches only
    // (six: two jobs, two budgets), whose five deferred tests turn FALSE, and the last job's
    // matches are computed after the better list is kept and are FALSE for certain. Two and
    // eight: every test after the third match's turns FALSE.
    for (jobs, turned) in [(1, 5), (2, 7), (8, 7)] {
        let limits = XWindowLimits {
            jobs,
            at_risk_jobs: jobs,
        };
        let mut heaps = heaps_of(1, 1);
        let (counts, computed) = sw_drive(&mut lists(), &cases, &mut heaps, &alone(1), limits);
        let counts = counts.expect("no error");
        let mut model = heaps_of(1, 1);
        sw_one_thread(&lists(), &cases, &mut model, &alone(1)).expect("no error");
        assert_eq!(heap_contents(&heaps), heap_contents(&model), "jobs {jobs}");
        assert_eq!(
            (counts.turned, counts.continued, computed),
            (Some(turned), 2, 12),
            "jobs {jobs}"
        );
    }
    // The deferred TRUE at a window's last match (nothing carried after it) and at the
    // stream's last match (the last window holds the continuation only).
    let lists = || {
        vec![
            hit_list(vec![list(0, 1e-9, false), list(1, 1e-9, false)]),
            hit_list(vec![list(2, 1e-9, false), list(3, 1e-9, false)]),
        ]
    };
    let pair = vec![sw_case(0.5, SwKind::NotAsked), sw_case(0.5, SwKind::Kept)];
    let cases = vec![pair.clone(), pair];
    for jobs in [1, 8] {
        let limits = XWindowLimits {
            jobs,
            at_risk_jobs: jobs,
        };
        let mut heaps = heaps_of(2, 1);
        let (counts, computed) = sw_drive(&mut lists(), &cases, &mut heaps, &alone(2), limits);
        let counts = counts.expect("no error");
        let mut model = heaps_of(2, 1);
        sw_one_thread(&lists(), &cases, &mut model, &alone(2)).expect("no error");
        assert_eq!(heap_contents(&heaps), heap_contents(&model), "jobs {jobs}");
        assert_eq!(
            (counts.turned, counts.continued, computed),
            (Some(0), 2, 6),
            "jobs {jobs}"
        );
    }
    // Two queries, each a match that inserts its list and a better one whose deferred test
    // is TRUE. At the stop in the first query, the second query's carried matches are
    // replayed ahead on a copy of its heap: its test is TRUE there, so both continuations run
    // in the second window.
    let cases = vec![vec![sw_case(0.5, SwKind::Kept), sw_case(0.4, SwKind::Kept)]; 2];
    let limits = XWindowLimits {
        jobs: 8,
        at_risk_jobs: 8,
    };
    let mut heaps = heaps_of(2, 1);
    let (counts, computed) = sw_drive(&mut lists(), &cases, &mut heaps, &alone(2), limits);
    let counts = counts.expect("no error");
    let mut model = heaps_of(2, 1);
    sw_one_thread(&lists(), &cases, &mut model, &alone(2)).expect("no error");
    assert_eq!(heap_contents(&heaps), heap_contents(&model));
    assert_eq!(
        (counts.turned, counts.continued, counts.windows, computed),
        (Some(0), 2, 2, 6)
    );
}

/// At a stop, a carried deferred test that `x_sw_continued` leaves to the replay because the
/// early termination test may end its match still counts as one that may insert a list, so the
/// bound continues no later match of the query once the heap may be filled before it (review 2 of
/// X4-mt, N-5: the bound reads only the heap's count). A continuation's result is replayed without
/// its heap test, so a later match continued by mistake would insert a list that NCBI's loop does
/// not keep.
#[test]
fn a_carried_test_that_the_early_termination_test_may_end_counts_toward_the_bound() {
    // One query, two lists to keep. The first two matches ask no test (no insert); the third's
    // deferred test is TRUE, so the replay stops there; the fourth's preliminary E-value lets the
    // early termination test end it once the heap is filled, which it is not yet when its
    // deferred test is TRUE (the second stop); the fifth's deferred test turns FALSE after the
    // third's and fourth's lists are kept, and its result after a TRUE would be kept before them.
    let lists = || {
        vec![hit_list(vec![
            list(0, 1e-9, false),
            list(1, 1e-9, false),
            list(2, 1e-9, false),
            list(3, 0.5, false),
            list(4, 1e-9, false),
        ])]
    };
    let kept = |final_evalue| SwCase {
        evalue: 0.5,
        score: 40,
        final_evalue,
        final_score: 40,
        kind: SwKind::Kept,
    };
    let cases = vec![vec![
        sw_case(0.5, SwKind::NotAsked),
        sw_case(0.5, SwKind::NotAsked),
        kept(1e-5),
        kept(1e-5),
        kept(1e-30),
    ]];
    let one_window = XWindowLimits {
        jobs: 8,
        at_risk_jobs: 8,
    };
    let mut heaps = heaps_of(1, 2);
    let (counts, computed) = sw_drive(&mut lists(), &cases, &mut heaps, &alone(1), one_window);
    let counts = counts.expect("no error");
    let mut model = heaps_of(1, 2);
    sw_one_thread(&lists(), &cases, &mut model, &alone(1)).expect("no error");
    assert_eq!(heap_contents(&heaps), heap_contents(&model));
    assert_eq!(heaps[0].n, 2, "the third's and fourth's lists");
    // Two stops, each continued in the next window, nothing by the bound; the fifth's test
    // turned FALSE. Five matches in the first window and two continuations.
    assert_eq!(
        (
            counts.turned,
            counts.continued,
            counts.continued_by.stop,
            counts.continued_by.bound,
            counts.windows,
            computed
        ),
        (Some(1), 2, 2, 0, 3, 7)
    );
}

/// A continuation is replayed before the end of its batch: its redo sets the batch's status
/// (clearing a failure before it, or failing the batch) or ends the search (an endless redo),
/// before any match of a later batch.
#[test]
fn a_continuation_sets_its_batch_status_before_the_batch_ends() {
    let one_job = XWindowLimits {
        jobs: 1,
        at_risk_jobs: 1,
    };
    let one_window = XWindowLimits {
        jobs: 64,
        at_risk_jobs: 64,
    };
    let lists = || {
        vec![
            hit_list(vec![list(0, 1e-9, false), list(1, 1e-9, false)]),
            hit_list(vec![list(2, 1e-9, false)]),
        ]
    };
    let later = vec![sw_case(1e-20, SwKind::Endless)];
    for (first, second, expected) in [
        // A failed redo, then the continuation's redo: status 0, the later batch ends it.
        (
            sw_case(1e-20, SwKind::Failed),
            sw_case(0.5, SwKind::Kept),
            endless_in(1),
        ),
        // The continuation's redo fails: the batch fails at its end.
        (
            sw_case(0.5, SwKind::NotAsked),
            sw_case(0.5, SwKind::Failed),
            batch_failed_text(0, false),
        ),
        // The continuation's redo never ends.
        (
            sw_case(0.5, SwKind::NotAsked),
            sw_case(0.5, SwKind::Endless),
            endless_in(0),
        ),
    ] {
        let cases = vec![vec![first, second], later.clone()];
        for limits in [one_job, one_window] {
            let mut heaps = heaps_of(2, 1);
            let (result, computed) = sw_drive(&mut lists(), &cases, &mut heaps, &alone(2), limits);
            let mut model = heaps_of(2, 1);
            let reference = sw_one_thread(&lists(), &cases, &mut model, &alone(2));
            let result = result.map(|_| ());
            assert_eq!(error_text(&reference), expected);
            assert_eq!(error_text(&result), error_text(&reference), "{limits:?}");
            assert_eq!(heap_contents(&heaps), heap_contents(&model), "{limits:?}");
            // The second match's job and its continuation.
            assert!(computed >= 3, "{limits:?}");
        }
    }
}

/// The batch's status crosses windows: a failed redo at the end of a batch fails it when the
/// replay passes the batch's end (in a later window, or after the last one); a later match of
/// the batch that the redo reaches clears it, one that the early termination test ends does
/// not; and a failed batch comes before an endless redo of a later batch computed ahead.
#[test]
fn the_batch_status_crosses_windows_and_fails_the_batch_at_its_end() {
    let one_job = XWindowLimits {
        jobs: 1,
        at_risk_jobs: 1,
    };
    let one_window = XWindowLimits {
        jobs: 64,
        at_risk_jobs: 64,
    };
    // Batch 0..2: query 0 ends with a failed redo; query 1's matches are ended by the test
    // (its heap and query 0's filled to the cutoff, E 1.0) or reached (E 1e-9).
    let run = |second: f64, limits: XWindowLimits, batches: Vec<Range<usize>>| {
        let lists = vec![
            hit_list(vec![list(1, 1e-9, false), list(2, 1e-9, false)]),
            hit_list(vec![list(3, 1e-9, false), list(4, second, false)]),
            hit_list(vec![list(5, 1e-9, false)]),
        ];
        let cases = vec![
            vec![sw_case(1e-20, SwKind::Kept), sw_case(1e-20, SwKind::Failed)],
            vec![sw_case(1e-20, SwKind::Kept), sw_case(1e-20, SwKind::Kept)],
            vec![sw_case(1e-20, SwKind::Endless)],
        ];
        let mut heaps = heaps_of(3, 1);
        let (result, _) = sw_drive(&mut lists.clone(), &cases, &mut heaps, &batches, limits);
        let mut model = heaps_of(3, 1);
        let reference = sw_one_thread(&lists, &cases, &mut model, &batches);
        let result = result.map(|_| ());
        assert_eq!(error_text(&result), error_text(&reference));
        assert_eq!(heap_contents(&heaps), heap_contents(&model));
        error_text(&result)
    };
    let failed = batch_failed_text(0, true);
    let endless = endless_in(2);
    for limits in [one_job, one_window] {
        // Query 1's first match is reached and clears the status; batch 2's redo never ends.
        assert_eq!(run(1e-9, limits, vec![0..2, 0..2, 2..3]), endless);
        // Query 1 alone in a batch: query 0's failed status fails batch 0 before batch 1.
        assert_eq!(run(1e-9, limits, vec![0..1, 1..2, 2..3]), failed);
        // All three in one batch: the endless redo of query 2 ends it.
        assert_eq!(run(1e-9, limits, vec![0..3; 3]), endless_in(0));
    }
    // A match ended by the test leaves the status that the match before it set: query 1's
    // third match (E 1.0) comes after both heaps are filled to the cutoff.
    let lists = vec![
        hit_list(vec![list(1, 1e-9, false), list(2, 1e-9, false)]),
        hit_list(vec![
            list(3, 1e-9, false),
            list(4, 1e-9, false),
            list(5, 1.0, false),
        ]),
    ];
    let cases = vec![
        vec![sw_case(1e-20, SwKind::Kept), sw_case(1e-20, SwKind::Kept)],
        vec![
            sw_case(1e-20, SwKind::Kept),
            sw_case(1e-20, SwKind::Failed),
            sw_case(1e-20, SwKind::Kept),
        ],
    ];
    for limits in [one_job, one_window] {
        let mut heaps = heaps_of(2, 1);
        let (result, _) = sw_drive(
            &mut lists.clone(),
            &cases,
            &mut heaps,
            &[0..2, 0..2],
            limits,
        );
        let error = error_text(&result.map(|_| ()));
        // Query 0 holds two lists: the batch's failure where NCBI's clean-up crashes.
        assert_eq!(error, batch_failed_text(0, true));
        let mut model = heaps_of(2, 1);
        let reference = sw_one_thread(&lists, &cases, &mut model, &[0..2, 0..2]);
        assert_eq!(error_text(&reference), error);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3540
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3560-3575
/// A list that the fetch excludes (line RG X4-excl) after a failed redo clears the batch's
/// failure unless the early termination test ends it there; one before the failed redo does
/// not count. The replay reads a query's marks when it passes the query: across a window's
/// border (one job per window), at the stream's end, and for a query of the batch without a
/// match, before the batch ends; a query of another batch does not count.
#[test]
fn an_excluded_list_after_a_failed_redo_clears_the_failure_unless_the_redo_ends_early_there() {
    let one_job = XWindowLimits {
        jobs: 1,
        at_risk_jobs: 1,
    };
    let one_window = XWindowLimits {
        jobs: 64,
        at_risk_jobs: 64,
    };
    let mark = |query, position, best_evalue| FetchExcludedMatch {
        query,
        position,
        best_evalue,
    };
    let run = |lists: &[Option<BlastpHitList>],
               cases: &[Vec<SwCase>],
               batches: &[Range<usize>],
               marks: &[FetchExcludedMatch],
               limits: XWindowLimits| {
        let mut heaps = heaps_of(lists.len(), 1);
        let (result, _) = sw_drive_excluded(
            &mut lists.to_vec(),
            cases,
            &mut heaps,
            batches,
            limits,
            marks,
            |_| {},
        );
        let mut model = heaps_of(lists.len(), 1);
        let reference = sw_one_thread_excluded(lists, cases, &mut model, batches, marks);
        let result = result.map(|_| ());
        let at = format!("{marks:?} {limits:?}");
        assert_eq!(error_text(&result), error_text(&reference), "{at}");
        assert_eq!(heap_contents(&heaps), heap_contents(&model), "{at}");
        error_text(&result).is_some()
    };
    // Query 0 alone: three kept redos fill its heap to the cutoff, a failed redo at place 3
    // (the last match of the first job), and a list at place 4 (E 1.0) that the test ends.
    let lists = vec![hit_list(vec![
        list(1, 1e-9, false),
        list(2, 1e-9, false),
        list(3, 1e-9, false),
        list(4, 1e-9, false),
        list(5, 1.0, false),
    ])];
    let kept = sw_case(1e-20, SwKind::Kept);
    let cases = vec![vec![kept, kept, kept, sw_case(1e-20, SwKind::Failed), kept]];
    for (marks, fails) in [
        (vec![], true),
        // Before the failed redo: it sets the status again.
        (vec![mark(0, 3, 1e-9)], true),
        // After it, between the windows of one job: the test does not end it (E 1e-9).
        (vec![mark(0, 4, 1e-9)], false),
        // After it, but the test ends it (E 1.0, above 5 times the cutoff of the filled heap).
        (vec![mark(0, 4, 1.0)], true),
        // After the last list: read at the stream's end.
        (vec![mark(0, 5, 1e-9)], false),
        (vec![mark(0, 4, 1.0), mark(0, 5, 0.001)], false),
        (vec![mark(0, 2, 1e-9), mark(0, 5, 1.0)], true),
    ] {
        for limits in [one_job, one_window] {
            assert_eq!(
                run(&lists, &cases, &alone(1), &marks, limits),
                fails,
                "{marks:?} {limits:?}"
            );
        }
    }
    // Batch 0..2: query 0's failed redo, then query 1, whose lists the fetch all excludes (no
    // match); query 2 in the next batch. Query 1's heap is empty, so the test ends none of its
    // lists: the batch's failure is cleared before it ends, but not by a list of another batch.
    let lists = vec![
        hit_list(vec![list(1, 1e-9, false), list(2, 1e-9, false)]),
        hit_list(Vec::new()),
        hit_list(vec![list(3, 1e-9, false)]),
    ];
    let cases = vec![
        vec![kept, sw_case(1e-20, SwKind::Failed)],
        Vec::new(),
        vec![kept],
    ];
    for (batches, marks, fails) in [
        (vec![0..2, 0..2, 2..3], vec![], true),
        (vec![0..2, 0..2, 2..3], vec![mark(1, 0, 1.0)], false),
        (vec![0..1, 1..2, 2..3], vec![mark(1, 0, 1.0)], true),
        (vec![0..1, 1..3, 1..3], vec![mark(1, 0, 1.0)], true),
    ] {
        for limits in [one_job, one_window] {
            assert_eq!(
                run(&lists, &cases, &batches, &marks, limits),
                fails,
                "{batches:?} {marks:?} {limits:?}"
            );
        }
    }
}

/// The windows of the Smith-Waterman redo give the one-thread loop's heaps and its error (the
/// failed batch, the endless redo, or none) on random searches: queries in batches, matches
/// whose first-round tests are certain, deferred and turned, failed and endless redos, early
/// termination, `-max_target_seqs` 1 to 500, windows of one job up to a few, small caps on
/// jobs at risk; in every other search, lists that the fetch excludes at random places, which
/// NCBI's loop reaches in its order (`sw_one_thread_excluded`). A match's result has its own
/// E-value and score, drawn apart from its first round's (the test's): the look-ahead of
/// `x_sw_continued` must insert with the one and test with the other. Each way a carried match
/// is continued (the stop, the exact look-ahead, the bound) is reached many times.
#[test]
fn smith_waterman_windows_give_the_one_thread_heaps_and_error_on_random_searches() {
    let mut rng = XorShift(0x9e37_79b9_7f4a_7c15);
    // The excluded lists' own sequence, so that the searches are drawn as without them.
    let mut mark_rng = XorShift(0x2545_f491_4f6c_dd1d);
    // The results' E-values and scores, their own sequence too.
    let mut final_rng = XorShift(0x5851_f42d_4c95_7f2d);
    let mut changed_by_excluded = 0;
    let prelim = [1e-30, 1e-5, 0.001, 0.005, 0.02, 1.0];
    let first_round = [1e-30, 1e-5, ECUTOFF, 0.004, 0.01, 0.5, 2.0];
    let (mut turned, mut continued, mut ended) = (0, 0, 0);
    let mut continued_by = XContinuedBy::default();
    let (mut batch_failures, mut endless, mut clean) = (0, 0, 0);
    for run in 0..1500 {
        let queries = 1 + rng.below(8);
        let mut batches = Vec::new();
        while batches.len() < queries {
            let start = batches.len();
            let end = (start + 1 + rng.below(3)).min(queries);
            batches.extend((start..end).map(|_| start..end));
        }
        let heap_threshold = [1, 2, 3, 5, 500][rng.below(5)];
        let mut lists = Vec::new();
        let mut cases = Vec::new();
        for _ in 0..queries {
            if rng.below(6) == 0 {
                lists.push(None);
                cases.push(Vec::new());
                continue;
            }
            let count = rng.below(25);
            let mut evalues: Vec<f64> = (0..count)
                .map(|_| prelim[rng.below(prelim.len())])
                .collect();
            evalues.sort_by(f64::total_cmp);
            lists.push(hit_list(
                evalues
                    .iter()
                    .enumerate()
                    .map(|(oid, &evalue)| list(oid as u32, evalue, rng.below(10) == 0))
                    .collect(),
            ));
            cases.push(
                (0..count)
                    .map(|_| {
                        let evalue = first_round[rng.below(first_round.len())];
                        let score = 30 + rng.below(4) as i32;
                        let kind = match rng.below(100) {
                            0..=9 => SwKind::NotAsked,
                            10..=19 => SwKind::Purged,
                            20..=29 => SwKind::Failed,
                            30 if run % 4 == 0 => SwKind::Endless,
                            _ => SwKind::Kept,
                        };
                        // The result: a first-round E-value, or three times or a seventh of the
                        // test's, and the test's score within three.
                        let final_evalue = match final_rng.below(first_round.len() + 2) {
                            pick if pick < first_round.len() => first_round[pick],
                            pick if pick == first_round.len() => evalue * 3.0,
                            _ => evalue / 7.0,
                        };
                        let final_score = score - 3 + final_rng.below(7) as i32;
                        SwCase {
                            evalue,
                            score,
                            final_evalue,
                            final_score,
                            kind,
                        }
                    })
                    .collect(),
            );
        }
        let limits = XWindowLimits {
            jobs: 1 + rng.below(6),
            at_risk_jobs: 1 + rng.below(3),
        };
        let mut excluded = Vec::new();
        if run % 2 == 1 {
            for (q_idx, hit_list) in lists.iter().enumerate() {
                let Some(hit_list) = hit_list else {
                    continue;
                };
                let mut places: Vec<usize> = (0..mark_rng.below(3))
                    .map(|_| mark_rng.below(hit_list.hsplist_count + 1))
                    .collect();
                places.sort_unstable();
                excluded.extend(places.into_iter().map(|position| FetchExcludedMatch {
                    query: q_idx,
                    position,
                    best_evalue: prelim[mark_rng.below(prelim.len())],
                }));
            }
        }
        let mut heaps = heaps_of(queries, heap_threshold);
        let (result, _) = sw_drive_excluded(
            &mut lists.clone(),
            &cases,
            &mut heaps,
            &batches,
            limits,
            &excluded,
            |_| {},
        );
        let mut model = heaps_of(queries, heap_threshold);
        let reference = sw_one_thread_excluded(&lists, &cases, &mut model, &batches, &excluded);
        // The searches whose outcome the excluded lists change (18 of the 1,500).
        if !excluded.is_empty() {
            let mut unmarked = heaps_of(queries, heap_threshold);
            let without = sw_one_thread(&lists, &cases, &mut unmarked, &batches);
            changed_by_excluded += usize::from(error_text(&without) != error_text(&reference));
        }
        let at = format!("run {run}: {limits:?} threshold {heap_threshold} {batches:?}");
        if let Ok(counts) = &result {
            turned += counts.turned.expect("the Smith-Waterman counts");
            continued += counts.continued;
            continued_by.add(counts.continued_by);
            ended += counts.skipped + counts.wasted;
            let by = counts.continued_by;
            assert_eq!(by.stop + by.ahead + by.bound, counts.continued);
        }
        let result = result.map(|_| ());
        assert_eq!(error_text(&result), error_text(&reference), "{at}");
        assert_eq!(heap_contents(&heaps), heap_contents(&model), "{at}");
        match error_text(&reference) {
            None => clean += 1,
            Some(text) if text.contains("does not end") => endless += 1,
            Some(_) => batch_failures += 1,
        }
    }
    assert!(
        turned > 100
            && continued > 100
            && continued_by.stop > 100
            && continued_by.ahead > 20
            && continued_by.bound > 20
            && ended > 100
            && batch_failures > 50
            && endless > 20
            && clean > 300
            && changed_by_excluded > 15,
        "turned {turned} continued {continued} {continued_by:?} ended {ended} failed batches \
             {batch_failures} endless {endless} clean {clean} changed by excluded lists \
             {changed_by_excluded}"
    );
}

/// A query whose replay stops at each of its matches after the first: a list kept above the
/// cutoff (E 0.9, preliminary E 0.5) and then `stops` matches each better than the last
/// (E 0.8 / (k + 1)), so that with `-max_target_seqs 1` each one's deferred heap test is TRUE at
/// its replay.
fn stopping_query(stops: u32) -> (Option<BlastpHitList>, Vec<SwCase>) {
    let cases = (0..=stops)
        .map(|position| {
            let evalue = if position == 0 {
                0.9
            } else {
                0.8 / f64::from(position + 1)
            };
            SwCase {
                evalue,
                score: 30,
                final_evalue: evalue,
                final_score: 30,
                kind: SwKind::Kept,
            }
        })
        .collect();
    (
        hit_list((0..=stops).map(|oid| list(oid, 0.5, false)).collect()),
        cases,
    )
}

/// The carry-growth order of review 1 of X4-mt (S1), at a budget of 16 jobs: query 0, alone in its
/// batch with `-max_target_seqs 1`, stops the replay at each of its 40 matches after the first;
/// 2,000 queries of four matches within the cutoff follow (one job each).
fn carry_growth_order() -> (Vec<Option<BlastpHitList>>, Vec<Vec<SwCase>>, XWindowLimits) {
    const QUERIES: usize = 1 + 2000;
    let (first, first_cases) = stopping_query(40);
    let mut lists = vec![first];
    let mut cases = vec![first_cases];
    for _ in 1..QUERIES {
        lists.push(hit_list((0..4).map(|oid| list(oid, 1e-9, false)).collect()));
        cases.push(vec![sw_case(1e-9, SwKind::Kept); 4]);
    }
    let limits = XWindowLimits {
        jobs: 16,
        at_risk_jobs: 8,
    };
    (lists, cases, limits)
}

/// The carry of a stopped replay is bounded (review 1 of X4-mt, S1), in the carry-growth order. A
/// window takes at most `jobs` new jobs (the first clause of the budget: review 2's mutant S1b,
/// which drops it, takes 32 in the first window), and with the carried matches (in jobs of four)
/// at most `2 x jobs`, so it holds at most `2 x X_JOB_MATCHES x jobs` matches with the carried
/// ones, and a carry at most one fewer; where the carried matches did not count, the carry grew by
/// about `X_JOB_MATCHES x (jobs - 1)` matches per stopped window (2,329 matches here, with every
/// carried match's result and its query's list held).
#[test]
fn a_carry_counts_toward_the_window_budget_and_stays_within_one_window() {
    let (lists, cases, limits) = carry_growth_order();
    let batches = alone(lists.len());
    let mut heaps = heaps_of(lists.len(), 1);
    let mut held = Vec::new();
    let (counts, _) = sw_drive_excluded(
        &mut lists.clone(),
        &cases,
        &mut heaps,
        &batches,
        limits,
        &[],
        |window| held.push((window.carried, window.matches.len(), window.stream_jobs())),
    );
    let counts = counts.expect("no error");
    let mut model = heaps_of(lists.len(), 1);
    sw_one_thread(&lists, &cases, &mut model, &batches).expect("no error");
    assert_eq!(heap_contents(&heaps), heap_contents(&model));
    let bound = 2 * X_JOB_MATCHES * limits.jobs;
    assert_eq!(
        held.iter().map(|&(_, _, new_jobs)| new_jobs).max(),
        Some(limits.jobs)
    );
    assert_eq!(
        held.iter().map(|&(_, matches, _)| matches).max(),
        Some(bound)
    );
    assert_eq!(
        held.iter().map(|&(carried, _, _)| carried).max(),
        Some(bound - 1)
    );
    assert_eq!(counts.carried_max, bound - 1);
    assert_eq!(
        (counts.continued, counts.continued_by),
        (
            40,
            XContinuedBy {
                stop: 40,
                ahead: 0,
                bound: 0
            }
        )
    );
}

/// The expanded lists that a carry holds (review 2 of X4-mt, S-3): the carry-growth order with its
/// compact lists expanded ahead of the windows (BPREDOWINEXPAND). After each window its reads and
/// the next match lie within the expanded lists (`XExpandAhead::check`, which `x_redo_windows`
/// asserts; review 2's mutant S1b, a window of up to two budgets of new jobs, reads past them);
/// the heaps are the one-thread loop's, and a window takes at most one budget of new jobs. While a
/// window's jobs run, the lists held besides query 0's (each list whole) are, without a carry, those
/// of one budget and the next job read ahead (`X_JOB_MATCHES x (jobs + 1)` matches), and with a
/// carry also those of the queries of the carried matches, up to the carry's bound: three budgets
/// here (31 queries carried behind the window's start, with query 0's last three matches, and 17
/// read ahead).
#[test]
fn a_carry_holds_the_expanded_lists_of_about_three_budgets() {
    let (lists, cases, limits) = carry_growth_order();
    let count = lists.len();
    let batches = alone(count);
    let mut compact: Vec<_> = lists.iter().map(compact_of).collect();
    let mut ahead_lists: Vec<Option<BlastpHitList>> = vec![None; count];
    let mut expand = XExpandAhead::new(&ahead_lists, &mut compact);
    let mut heaps = heaps_of(count, 1);
    // Per window: whether it has a carry, its new jobs, and the matches of the lists held after
    // query 0's.
    let mut seen = Vec::new();
    let (counts, _) = sw_drive_expanding(
        &mut ahead_lists,
        &cases,
        &mut heaps,
        &batches,
        limits,
        Some(&mut expand),
        &[],
        |lists, window| {
            let held: usize = lists[1..].iter().flatten().map(x_list_matches).sum();
            seen.push((window.carried > 0, window.stream_jobs(), held));
        },
    );
    let counts = counts.expect("no error");
    let mut model = heaps_of(count, 1);
    sw_one_thread(&lists, &cases, &mut model, &batches).expect("no error");
    assert_eq!(heap_contents(&heaps), heap_contents(&model));
    assert_eq!(counts.continued, 40);
    drop(expand);
    assert!(compact.iter().all(Option::is_none), "every list expanded");
    assert!(
        ahead_lists.iter().all(Option::is_none),
        "every list released"
    );
    assert_eq!(
        seen.iter().map(|&(_, new_jobs, _)| new_jobs).max(),
        Some(limits.jobs)
    );
    let most_held = |carry: bool| {
        seen.iter()
            .filter(|&&(carried, _, _)| carried == carry)
            .map(|&(_, _, held)| held)
            .max()
    };
    assert_eq!(most_held(false), Some(X_JOB_MATCHES * (limits.jobs + 1)));
    assert_eq!(most_held(true), Some(3 * X_JOB_MATCHES * limits.jobs));
}

/// The results at risk that stopped replays carry on are bounded per stopped window (review 1 of
/// X4-mt, S2, as restated in X4mt-9). Query 0 (alone in its batch, `-max_target_seqs 1`) stops the
/// replay `stops` times; queries 1 and 2 form a batch whose heaps fill to the cutoff at query 2's
/// first match, after which the early termination test ends query 2's 399 others. Query 2's cap
/// does not grow while the replays stop before it, so each stopped window takes at most that cap
/// of new jobs at risk of query 2 and carries their results on, and the one window that drops
/// query 2's results drops at most `4 x cap x (stops + 1)` of them (`16 x (stops + 1) - 1` here:
/// query 2's first match fills the heaps) and at most the results a window holds with the carried
/// ones (`2 x 4 x 32`; 251 after 20 stops). Where a stopped window's jobs at risk grew the cap,
/// 1, 2 and 5 stops dropped 47, 111 and 251 results.
#[test]
fn a_stopped_replay_drops_at_most_a_cap_of_results_at_risk_per_stopped_window() {
    let limits = XWindowLimits {
        jobs: 32,
        at_risk_jobs: 4,
    };
    for (stops, wasted) in [(0, 15), (1, 31), (2, 47), (5, 95), (20, 251)] {
        let (first, first_cases) = stopping_query(stops as u32);
        let lists = vec![
            first,
            hit_list(vec![list(0, 1e-9, false)]),
            hit_list(
                (0..400)
                    .map(|oid| list(oid, if oid == 0 { 1e-5 } else { 1.0 }, false))
                    .collect(),
            ),
        ];
        let cases = vec![
            first_cases,
            vec![sw_case(1e-9, SwKind::Kept)],
            vec![sw_case(1e-9, SwKind::Kept); 400],
        ];
        let batches = [0..1, 1..3, 1..3];
        let mut heaps = heaps_of(3, 1);
        let (counts, computed) = sw_drive(&mut lists.clone(), &cases, &mut heaps, &batches, limits);
        let counts = counts.expect("no error");
        let mut model = heaps_of(3, 1);
        sw_one_thread(&lists, &cases, &mut model, &batches).expect("no error");
        assert_eq!(
            heap_contents(&heaps),
            heap_contents(&model),
            "{stops} stops"
        );
        let bound = (X_JOB_MATCHES * limits.at_risk_jobs * (stops + 1))
            .min(2 * X_JOB_MATCHES * limits.jobs);
        assert!(
            counts.wasted <= bound,
            "{stops} stops: {} results dropped, bound {bound}",
            counts.wasted
        );
        assert_eq!(
            (counts.wasted, counts.skipped, counts.continued),
            (wasted, 399 - wasted, stops),
            "{stops} stops"
        );
        // Query 2's first match and its matches at risk, query 1's match, query 0's matches and
        // their continuations.
        assert_eq!(
            computed,
            (1 + wasted) + 1 + (stops + 1) + stops,
            "{stops} stops"
        );
    }
}

/// The Smith-Waterman jobs on real matches (`x_redo_match_sw`): at 1, 2, 3 and 7 threads, with
/// and without the jobs largest first, each computed match equals the one-thread loop's call
/// (`postprocess_preliminary_hits_smith_waterman` with fresh scratch) whose heap test the
/// query's heap at the window's start answers, and the job defers exactly the tests that this
/// heap cannot answer for certain.
#[test]
fn smith_waterman_jobs_equal_the_one_thread_call_on_the_window_start_heaps() {
    let fixture = redo_fixture();
    let count = fixture.lists.len();
    // Two sets of heaps that keep two lists. The first: by query, empty, one list, and full
    // above the cutoff (worst E 0.5). The second: full within the cutoff (worst E 1e-25), but
    // for each batch's last query, which stays empty so that no batch ends early: there every
    // test of an E-value above the cutoff is FALSE for certain (the fixture asks such tests of
    // the first query of each batch).
    let heap_sets: Vec<Vec<BlastCompoHeap>> = [false, true]
        .into_iter()
        .map(|within| {
            let mut heaps = heaps_of(count, 2);
            for (q_idx, heap) in heaps.iter_mut().enumerate() {
                let kept: &[(f64, i32)] = if !within {
                    &[(1e-30, 300), (0.5, 60)][..q_idx % 3]
                } else if q_idx + 1 == fixture.batches[q_idx].end {
                    &[]
                } else {
                    &[(1e-30, 300), (1e-25, 250)]
                };
                for (index, &(evalue, score)) in kept.iter().enumerate() {
                    let oid = 1000 + index as u32;
                    heap.insert(list(oid, evalue, false), evalue, score, oid);
                }
            }
            heaps
        })
        .collect();
    let build = |q_idx: usize| -> Result<BlastRedoAlignParams> {
        Ok(BlastRedoAlignParams {
            smith_waterman: true,
            ..redo_params(&fixture, q_idx)
        })
    };
    let (mut compared, mut with_list, mut deferred_tests, mut certain_false) = (0, 0, 0, 0);
    let mut continuations = 0;
    let stream_end = XMatch {
        q_idx: count,
        position: 0,
    };
    for heaps in &heap_sets {
        let all = stream(&fixture.lists);
        let mut rng = XorShift(0x2545_f491_4f6c_dd1d);
        let caps = vec![usize::MAX; count];
        let windows: Vec<XWindow> = (0..8)
            .map(|index| {
                let start = all[rng.below(all.len())];
                let mut window = x_window::<true>(
                    &fixture.lists,
                    heaps,
                    &fixture.batches,
                    &vec![None; count],
                    &caps,
                    None,
                    start,
                    2 + rng.below(10),
                );
                for ended in &mut window.ended {
                    *ended = index % 2 == 1 && rng.below(4) == 0;
                }
                window
            })
            .collect();
        // Windows after a replay that stopped at a match of each window above whose test is
        // asked (built at the first thread count, with the test's E-value and score).
        let mut carried_windows: Vec<XWindow> = Vec::new();
        for threads in [1, 2, 3, 7] {
            crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
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
                    sw_excluded: &[],
                };
                // The one-thread loop's call with fresh scratch, its heap test answered by
                // `answer`; and the test's E-value and score, if it was asked.
                let one_thread = |matched: XMatch, answer: &dyn Fn(f64, i32) -> bool| {
                    let q_idx = matched.q_idx;
                    let local_match = x_local_match(&fixture.lists, matched);
                    let hits: Vec<BlastpPreliminaryHsp> = local_match
                        .hsps
                        .iter()
                        .map(preliminary_hit_from_local_match_hsp)
                        .collect();
                    let mut scratch = XScratch::new(
                        BlastpKappaSubjectRangeCache::new(None),
                        ScoringMatrix::Blosum62,
                    );
                    let mut sw_scratch = SwMatchScratch::new();
                    let subject = search.subjects.read(local_match.oid).expect("a subject");
                    let mut asked = None;
                    let result = postprocess_preliminary_hits_smith_waterman(
                        &hits,
                        &x_query_workspace(&search, q_idx),
                        &mut scratch.composition_workspace,
                        query_raw(&fixture.contexts[q_idx]),
                        subject.residues(),
                        ScoringMatrix::Blosum62,
                        &fixture.gapped_params,
                        Some(&fixture.gapped_gumbel),
                        &build(q_idx).expect("parameters"),
                        10.0,
                        &mut scratch.gap_scratch,
                        &mut scratch.range_cache,
                        false,
                        &mut sw_scratch,
                        &mut |evalue, score| {
                            asked = Some((evalue, score));
                            answer(evalue, score)
                        },
                    );
                    (result.map_err(|error| format!("{error:#}")), asked)
                };
                if carried_windows.is_empty() {
                    for window in &windows {
                        let Some(index) = (0..window.matches.len()).find(|&index| {
                            !window.ended[index]
                                && one_thread(window.matches[index], &|_: f64, _: i32| false)
                                    .1
                                    .is_some()
                        }) else {
                            continue;
                        };
                        let matched = window.matches[index];
                        let asked = one_thread(matched, &|_: f64, _: i32| false)
                            .1
                            .expect("asked");
                        let stopped = XComputed {
                            result: no_alignment(),
                            deferred: Some(asked),
                            first_round_s: 0.0,
                        };
                        let carry = XCarry {
                            matches: window.matches[index..].to_vec(),
                            results: std::iter::once(Some(stopped))
                                .chain(window.ended[index + 1..].iter().map(|&ended| {
                                    (!ended).then(|| XComputed::from(no_alignment()))
                                }))
                                .collect(),
                        };
                        carried_windows.push(x_window::<true>(
                            &fixture.lists,
                            heaps,
                            &fixture.batches,
                            &vec![None; count],
                            &caps,
                            Some(&carry),
                            window.next.unwrap_or(stream_end),
                            carry.matches.len().div_ceil(X_JOB_MATCHES) + window.jobs.len() + 2,
                        ));
                    }
                }
                for window in windows.iter().chain(&carried_windows) {
                    for lpt in [false, true] {
                        let shared = Arc::new(XSharedSubjectRanges::new());
                        let (computed, _) = x_compute_window::<_, true>(
                            &search,
                            &fixture.lists,
                            heaps,
                            window,
                            &shared,
                            lpt,
                        );
                        let in_jobs: Vec<usize> =
                            window.jobs.iter().flat_map(|job| job.clone()).collect();
                        assert_eq!(computed.len(), in_jobs.len());
                        assert_eq!(
                            in_jobs.len() + window.carried,
                            window.matches.len() + window.continued.len(),
                            "the carried matches but the continued ones are in no job"
                        );
                        for (index, computed) in in_jobs.into_iter().zip(&computed) {
                            let matched = window.matches[index];
                            let at = (threads, lpt, matched);
                            let Some(computed) = computed else {
                                assert!(window.ended[index], "{at:?}");
                                continue;
                            };
                            let q_idx = matched.q_idx;
                            let oid = x_local_match(&fixture.lists, matched).oid;
                            let heap = &heaps[q_idx];
                            let before =
                                index - window.matches.partition_point(|other| other.q_idx < q_idx);
                            // A test the job cannot answer for certain on the window's start.
                            let deferrable = |evalue: f64| {
                                evalue > heap.ecutoff && heap.n + before >= heap.heap_threshold
                            };
                            let continued = window
                                .continued
                                .iter()
                                .find(|&&(at, _)| at == index)
                                .map(|&(_, test)| test);
                            // The one-thread call answered as the job answers: TRUE for a
                            // continuation, FALSE for a deferred test, else the heap at the
                            // window's start.
                            let (reference, asked) =
                                one_thread(matched, &|evalue: f64, score: i32| {
                                    continued.is_some()
                                        || (!deferrable(evalue)
                                            && heap.clone().would_insert(evalue, score, oid))
                                });
                            with_list += usize::from(
                                reference
                                    .as_ref()
                                    .is_ok_and(|result| result.hsp_list.is_some()),
                            );
                            assert_eq!(
                                x_result_difference(&computed.result, &reference),
                                None,
                                "{at:?}"
                            );
                            if let Some(test) = continued {
                                assert_eq!(asked, Some(test), "{at:?}");
                                assert_eq!(computed.deferred, None, "{at:?}");
                                continuations += 1;
                                continue;
                            }
                            let expected = asked.filter(|&(evalue, score)| {
                                deferrable(evalue) && heap.would_insert_view(evalue, score, oid)
                            });
                            assert_eq!(computed.deferred, expected, "{at:?}");
                            deferred_tests += usize::from(expected.is_some());
                            certain_false += usize::from(asked.is_some_and(|(evalue, score)| {
                                !heap.would_insert_view(evalue, score, oid)
                            }));
                            compared += 1;
                        }
                    }
                }
                Ok(())
            })
            .expect("a pool");
        }
    }
    assert!(
        compared > 150
            && with_list > 30
            && deferred_tests > 10
            && certain_false > 10
            && continuations > 50,
        "compared {compared} with a list {with_list} deferred {deferred_tests} certain FALSE \
             {certain_false} continuations {continuations}"
    );
}
