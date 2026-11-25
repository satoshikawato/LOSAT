use super::*;
use crate::algorithm::blastp::hsp::BlastpHsp;

fn make_hsp(oid: u32, e_value: f64) -> BlastpHsp {
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
        comp_adjustment_method: 0,
    }
}

fn make_list(oid: u32, query_index: u32, best_evalue: f64) -> BlastpHspList {
    BlastpHspList {
        oid,
        query_index,
        hsps: vec![make_hsp(oid, best_evalue)],
        best_evalue,
    }
}

fn make_hit_list(lists: Vec<BlastpHspList>) -> Option<BlastpHitList> {
    Some(BlastpHitList {
        hsplist_count: lists.len(),
        hsplist_max: 500,
        worst_evalue: 0.0,
        low_score: 0,
        heapified: false,
        hsplist_array: lists,
        hsplist_current: 0,
        num_hits: 0,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:91-96,158-204,678-686
#[test]
fn collect_batch_events_sorts_by_ordinal_asc_then_query_desc() {
    // Ranks 0, 1 map to DB ordinals 10, 30; excluded mark has ordinal 20 and ordinal 10 on q1.
    let lists = vec![
        make_hit_list(vec![make_list(1, 0, 1e-5), make_list(0, 0, 1e-2)]),
        make_hit_list(vec![make_list(1, 1, 1e-4)]),
        make_hit_list(vec![make_list(0, 2, 1e-9)]),
    ];
    let ordinals = [10u32, 30u32];
    let excluded = vec![
        FetchExcludedMatch {
            query: 1,
            position: 0,
            ordinal: 10,
            best_evalue: 0.5,
        },
        FetchExcludedMatch {
            query: 0,
            position: 1,
            ordinal: 20,
            best_evalue: 0.01,
        },
        // Query outside batch 0..2 is ignored.
        FetchExcludedMatch {
            query: 2,
            position: 0,
            ordinal: 5,
            best_evalue: 0.001,
        },
    ];
    let events = collect_batch_events(0..2, &lists, Some(&ordinals), &excluded);
    let keys: Vec<(u32, usize, bool)> = events
        .iter()
        .map(|e| {
            (
                e.ordinal(),
                e.q_idx(),
                matches!(e, Cbs0Event::Excluded { .. }),
            )
        })
        .collect();
    assert_eq!(
        keys,
        vec![
            (10, 1, true),  // ordinal 10, query 1 (excluded) before query 0
            (10, 0, false), // ordinal 10, query 0
            (20, 0, true),  // ordinal 20, query 0 (excluded)
            (30, 1, false), // ordinal 30, query 1 before query 0
            (30, 0, false), // ordinal 30, query 0
        ]
    );
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
mod par_tests {
    use super::super::par::{drive_oid_windows, OidComputed, OidLimits};
    use super::*;
    use crate::algorithm::blastp::kappa::BlastpPostprocessResult;
    use crate::algorithm::blastp::x_redo_par::x_sw_heap_test;
    use crate::core::composition_adjustment::smith_waterman_match::{
        SwBatchStatus, SwRedoEndless, SwRedoFailed,
    };

    struct XorShift(u64);
    impl XorShift {
        fn next(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }
        fn below(&mut self, n: usize) -> usize {
            (self.next() as usize) % n
        }
    }

    #[derive(Clone, Copy)]
    enum Kind {
        Kept,
        Purged,
        Failed,
        Endless,
        NotAsked,
    }

    #[derive(Clone, Copy)]
    struct Case {
        first_evalue: f64,
        first_score: i32,
        final_evalue: f64,
        final_score: i32,
        kind: Kind,
    }

    fn eval_case(
        oid: u32,
        q_idx: usize,
        case: Case,
        inserts: bool,
    ) -> Result<BlastpPostprocessResult> {
        if !inserts {
            return Ok(BlastpPostprocessResult {
                hsp_list: None,
                best_score: 0,
                best_evalue: f64::MAX,
            });
        }
        match case.kind {
            Kind::Kept => Ok(BlastpPostprocessResult {
                hsp_list: Some(make_list(oid, q_idx as u32, case.final_evalue)),
                best_score: case.final_score,
                best_evalue: case.final_evalue,
            }),
            Kind::Purged | Kind::NotAsked => Ok(BlastpPostprocessResult {
                hsp_list: None,
                best_score: 0,
                best_evalue: f64::MAX,
            }),
            Kind::Failed => Err(SwRedoFailed.into()),
            Kind::Endless => Err(SwRedoEndless.into()),
        }
    }

    fn serial_oid_reference(
        batch: Range<usize>,
        events: &[Cbs0Event],
        lists: &[Option<BlastpHitList>],
        cases: &[Vec<Case>],
        heaps: &mut [BlastCompoHeap],
    ) -> Result<()> {
        let mut sw_status = SwBatchStatus::default();
        let notes = FetchNotes::none();
        for &event in events {
            match event {
                Cbs0Event::Excluded { best_evalue, .. } => {
                    let term = blast_compo_early_termination(best_evalue, &heaps[batch.clone()]);
                    sw_status.reaches_redo(term, true);
                }
                Cbs0Event::Kept {
                    ordinal,
                    q_idx,
                    position,
                } => {
                    let local_match = get_m(lists, q_idx, position);
                    if blast_compo_early_termination(local_match.best_evalue, &heaps[batch.clone()])
                    {
                        continue;
                    }
                    let case = cases[q_idx][position];
                    let inserts = match case.kind {
                        Kind::NotAsked => true,
                        _ => heaps[q_idx].would_insert(
                            case.first_evalue,
                            case.first_score,
                            local_match.oid,
                        ),
                    };
                    let res = eval_case(local_match.oid, q_idx, case, inserts);
                    let Some(post) = sw_status.redone(res, || {
                        sw_oid_stop_error(&notes, batch.clone(), heaps, Some((ordinal, q_idx)))
                    })?
                    else {
                        continue;
                    };
                    if let Some(hsp_list) = post.hsp_list {
                        if heaps[q_idx].would_insert(
                            post.best_evalue,
                            post.best_score,
                            local_match.oid,
                        ) {
                            heaps[q_idx].insert(
                                hsp_list,
                                post.best_evalue,
                                post.best_score,
                                local_match.oid,
                            );
                        }
                    }
                }
            }
        }
        if sw_status.failed() {
            return Err(sw_oid_stop_error(&notes, batch, heaps, None));
        }
        Ok(())
    }

    fn heap_snapshot(heaps: &[BlastCompoHeap]) -> Vec<Vec<(u32, u64, i32)>> {
        heaps
            .iter()
            .map(|h| {
                let mut v: Vec<(u32, u64, i32)> = h
                    .records()
                    .iter()
                    .map(|r| (r.subject_index, r.best_evalue.to_bits(), r.best_score))
                    .collect();
                v.sort_unstable();
                v
            })
            .collect()
    }

    fn err_summary(r: &Result<()>) -> Option<String> {
        r.as_ref().err().map(|e| {
            match e.downcast_ref::<crate::blastinput::app::ErrorAfterBatchReports>() {
                Some(b) => format!("batch {:?}: {:#}", b.query, b.error),
                None => format!("{e:#}"),
            }
        })
    }

    // Proves that parallel OID-ordered windows match one-thread OID order across random multi-query
    // batches with deferred/turned/continued heap tests, early termination, excluded fetches,
    // failed batches, and endless matches.
    #[test]
    fn oid_windows_match_serial_reference_on_random_batches() {
        const ECUTOFF: f64 = 0.002;
        let mut rng = XorShift(0x7a12_3456_9abc_def1);
        let prelim_e = [1e-30, 1e-5, 0.001, 0.005, 0.02, 1.0];
        let first_e = [1e-30, 1e-5, ECUTOFF, 0.004, 0.01, 0.5, 2.0];
        let (mut total_turned, mut total_continued, mut total_ended) = (0, 0, 0);
        let notes = FetchNotes::none();

        for run in 0..800 {
            let num_queries = 1 + rng.below(5);
            let num_subjects = 2 + rng.below(18);
            let heap_threshold = [1, 2, 3, 5, 50][rng.below(5)];
            let mut lists = Vec::with_capacity(num_queries);
            let mut cases = Vec::with_capacity(num_queries);
            let mut excluded = Vec::new();
            for q_idx in 0..num_queries {
                let mut q_lists = Vec::new();
                let mut q_cases = Vec::new();
                for s in 0..num_subjects {
                    let roll = rng.below(10);
                    if roll < 3 {
                        continue;
                    }
                    let pe = prelim_e[rng.below(prelim_e.len())];
                    if roll == 3 && run % 2 == 1 {
                        excluded.push(FetchExcludedMatch {
                            query: q_idx,
                            position: q_lists.len(),
                            ordinal: s as u32,
                            best_evalue: pe,
                        });
                        continue;
                    }
                    q_lists.push(make_list(s as u32, q_idx as u32, pe));
                    let fe = first_e[rng.below(first_e.len())];
                    let fs = 30 + rng.below(5) as i32;
                    let kind = match rng.below(100) {
                        0..=9 => Kind::NotAsked,
                        10..=19 => Kind::Purged,
                        20..=29 => Kind::Failed,
                        30 if run % 4 == 0 => Kind::Endless,
                        _ => Kind::Kept,
                    };
                    q_cases.push(Case {
                        first_evalue: fe,
                        first_score: fs,
                        final_evalue: first_e[rng.below(first_e.len())],
                        final_score: fs - 2 + rng.below(5) as i32,
                        kind,
                    });
                }
                lists.push(make_hit_list(q_lists));
                cases.push(q_cases);
            }
            let batch = 0..num_queries;
            let events = collect_batch_events(batch.clone(), &lists, None, &excluded);
            let make_heaps = || {
                (0..num_queries)
                    .map(|_| BlastCompoHeap::new(heap_threshold, ECUTOFF))
                    .collect::<Vec<_>>()
            };
            let mut ref_heaps = make_heaps();
            let ref_res =
                serial_oid_reference(batch.clone(), &events, &lists, &cases, &mut ref_heaps);

            let limits = OidLimits {
                jobs: 1 + rng.below(5),
                at_risk_matches: 1 + rng.below(6),
            };
            let mut win_heaps = make_heaps();
            let win_res = drive_oid_windows(
                batch.clone(),
                &events,
                &lists,
                &mut win_heaps,
                limits,
                false,
                &notes,
                &mut |heaps, window, _carried| {
                    let mut out = Vec::new();
                    for job in &window.jobs {
                        for &idx in job {
                            let Cbs0Event::Kept {
                                q_idx, position, ..
                            } = window.events[idx]
                            else {
                                continue;
                            };
                            let local_match = get_m(&lists, q_idx, position);
                            let case = cases[q_idx][position];
                            let is_cont = window
                                .continued
                                .binary_search_by_key(&idx, |&(at, _)| at)
                                .is_ok();
                            let mut deferred = None;
                            let inserts = if is_cont || matches!(case.kind, Kind::NotAsked) {
                                true
                            } else {
                                x_sw_heap_test(
                                    &heaps[q_idx],
                                    window.before_same_query[idx],
                                    (case.first_evalue, case.first_score, local_match.oid),
                                    &mut deferred,
                                )
                            };
                            let result = eval_case(local_match.oid, q_idx, case, inserts);
                            out.push((
                                idx,
                                OidComputed {
                                    result,
                                    deferred,
                                    resume: None,
                                },
                            ));
                        }
                    }
                    out
                },
            );
            if let Ok(counts) = &win_res {
                total_turned += counts.turned;
                total_continued += counts.continued;
                total_ended += counts.skipped + counts.wasted;
            }
            let win_unit = win_res.map(|_| ());
            assert_eq!(
                err_summary(&win_unit),
                err_summary(&ref_res),
                "run {run} error mismatch"
            );
            if ref_res.is_ok() {
                assert_eq!(
                    heap_snapshot(&win_heaps),
                    heap_snapshot(&ref_heaps),
                    "run {run} heap mismatch"
                );
            }
        }
        assert!(total_turned > 0, "expected deferred tests turning FALSE");
        assert!(total_continued > 0, "expected continued matches");
        assert!(total_ended > 0, "expected early-terminated matches");
    }
}
