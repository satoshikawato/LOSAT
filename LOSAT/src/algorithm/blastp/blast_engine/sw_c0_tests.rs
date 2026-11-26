//! Unit tests for BLASTP's `-comp_based_stats 0 -use_sw_tback` OID-order stream and windowed
//! parallel replay (`sw_c0`).

use super::*;
use crate::algorithm::blastp::hsp::{fill_results_from_compo_heaps, BlastpHsp};
use crate::core::composition_adjustment::smith_waterman_match::{SwRedoEndless, SwRedoFailed};

fn dummy_hsp(oid: u32, q_idx: usize, score: i32, evalue: f64) -> BlastpHsp {
    BlastpHsp {
        num: 1,
        identity: 50.0,
        length: 20,
        mismatch: 10,
        gapopen: 0,
        q_start: 1,
        q_end: 20,
        s_start: 1,
        s_end: 20,
        gapped_q_start: 0,
        gapped_s_start: 0,
        e_value: evalue,
        bit_score: score as f64,
        num_ident: 10,
        query_context: 0,
        query_frame: 0,
        subject_frame: 0,
        query_length: 40,
        q_idx: q_idx as u32,
        s_idx: oid,
        raw_score: score,
        gap_info: None,
        num_positives: 12,
        comp_adjustment_method: 0,
    }
}

fn make_hsp_list(oid: u32, q_idx: usize, best_evalue: f64, empty: bool) -> BlastpHspList {
    BlastpHspList {
        oid,
        query_index: q_idx as u32,
        hsps: if empty {
            Vec::new()
        } else {
            vec![dummy_hsp(oid, q_idx, 50, best_evalue)]
        },
        best_evalue,
    }
}

fn make_hit_list(items: Vec<BlastpHspList>) -> Option<BlastpHitList> {
    let count = items.len();
    Some(BlastpHitList {
        hsplist_array: items,
        hsplist_count: count,
        hsplist_max: count.max(1),
        hsplist_current: 0,
        worst_evalue: f64::MAX,
        low_score: 0,
        heapified: false,
        num_hits: 0,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:91-96,158-207,678-686
#[test]
fn batch_stream_orders_by_database_ordinal_asc_and_query_desc_with_excluded_marks() {
    // Ranks 0, 1, 2 map to database ordinals 10, 30, 50; excluded marks sit at ordinals 20 and 50.
    let listed_ordinals = [10u32, 30, 50];
    let lists = vec![
        make_hit_list(vec![
            make_hsp_list(2, 0, 1e-10, false), // ordinal 50
            make_hsp_list(0, 0, 1e-5, false),  // ordinal 10
        ]),
        make_hit_list(vec![
            make_hsp_list(1, 1, 1e-8, false), // ordinal 30
            make_hsp_list(0, 1, 1e-2, false), // ordinal 10
            make_hsp_list(2, 1, 1e-1, true),  // empty: skipped
        ]),
        make_hit_list(vec![make_hsp_list(0, 2, 1e-3, false)]), // ordinal 10
    ];
    let excluded = vec![
        FetchExcludedMatch {
            query: 0,
            position: 1,
            ordinal: 20,
            best_evalue: 0.5,
        },
        FetchExcludedMatch {
            query: 2,
            position: 0,
            ordinal: 20,
            best_evalue: 0.25,
        },
    ];
    let stream = build_batch_stream(0, &(0..3), &lists, &excluded, Some(&listed_ordinals));
    let summary: Vec<(u32, usize, bool)> = stream
        .iter()
        .map(|e| {
            (
                e.ordinal,
                e.q_idx,
                matches!(e.kind, SwC0Kind::Excluded { .. }),
            )
        })
        .collect();
    assert_eq!(
        summary,
        vec![
            (10, 2, false),
            (10, 1, false),
            (10, 0, false),
            (20, 2, true),
            (20, 0, true),
            (30, 1, false),
            (50, 0, false),
        ]
    );
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
mod par_tests {
    use super::*;
    use crate::algorithm::blastp::x_redo_par::x_sw_heap_test;

    const ECUTOFF: f64 = 0.002;

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
        fn below(&mut self, bound: usize) -> usize {
            (self.next() as usize) % bound
        }
    }

    #[derive(Clone, Copy)]
    enum Kind {
        NotAsked,
        Purged,
        Kept,
        Failed,
        Endless,
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
        case: Case,
        oid: u32,
        q_idx: usize,
        answer: &mut dyn FnMut(f64, i32) -> SwHeapAnswer,
    ) -> Result<BlastpPostprocessResult> {
        match case.kind {
            Kind::NotAsked => Ok(BlastpPostprocessResult {
                hsp_list: None,
                best_score: 0,
                best_evalue: f64::MAX,
            }),
            _ => match answer(case.first_evalue, case.first_score) {
                SwHeapAnswer::Known(false) | SwHeapAnswer::Deferred => {
                    Ok(BlastpPostprocessResult {
                        hsp_list: None,
                        best_score: 0,
                        best_evalue: f64::MAX,
                    })
                }
                SwHeapAnswer::Known(true) => match case.kind {
                    Kind::NotAsked => unreachable!(),
                    Kind::Purged => Ok(BlastpPostprocessResult {
                        hsp_list: None,
                        best_score: case.final_score,
                        best_evalue: case.final_evalue,
                    }),
                    Kind::Kept => Ok(BlastpPostprocessResult {
                        hsp_list: Some(BlastpHspList {
                            oid,
                            query_index: q_idx as u32,
                            hsps: vec![dummy_hsp(oid, q_idx, case.final_score, case.final_evalue)],
                            best_evalue: case.final_evalue,
                        }),
                        best_score: case.final_score,
                        best_evalue: case.final_evalue,
                    }),
                    Kind::Failed => Err(SwRedoFailed.into()),
                    Kind::Endless => Err(SwRedoEndless.into()),
                },
            },
        }
    }

    fn one_thread_model(
        lists: &[Option<BlastpHitList>],
        cases: &[Vec<Case>],
        heaps: &mut [BlastCompoHeap],
        batches: &[Range<usize>],
        excluded: &[FetchExcludedMatch],
        listed_ordinals: Option<&[u32]>,
    ) -> Result<()> {
        let notes = FetchNotes::none();
        for (batch_idx, batch) in batches.iter().enumerate() {
            let stream = build_batch_stream(batch_idx, batch, lists, excluded, listed_ordinals);
            let mut status = SwBatchStatus::default();
            for entry in stream {
                match entry.kind {
                    SwC0Kind::Excluded { best_evalue } => {
                        let early =
                            blast_compo_early_termination(best_evalue, &heaps[batch.clone()]);
                        let _ = status.reaches_redo(early, true);
                    }
                    SwC0Kind::Kept { position } => {
                        let local_match = local_match_of(lists, entry.q_idx, position);
                        if blast_compo_early_termination(
                            local_match.best_evalue,
                            &heaps[batch.clone()],
                        ) {
                            continue;
                        }
                        let s_idx = local_match.oid;
                        let case = cases[entry.q_idx][position];
                        let redone = eval_case(case, s_idx, entry.q_idx, &mut |ev, sc| {
                            SwHeapAnswer::Known(heaps[entry.q_idx].would_insert(ev, sc, s_idx))
                        });
                        let Some(post) = status.redone(redone, || {
                            sw_endless_after_oid(
                                &notes,
                                batch.clone(),
                                heaps,
                                entry.ordinal,
                                entry.q_idx,
                            )
                        })?
                        else {
                            continue;
                        };
                        if let Some(hsp_list) = post.hsp_list {
                            if heaps[entry.q_idx].would_insert(
                                post.best_evalue,
                                post.best_score,
                                s_idx,
                            ) {
                                heaps[entry.q_idx].insert(
                                    hsp_list,
                                    post.best_evalue,
                                    post.best_score,
                                    s_idx,
                                );
                            }
                        }
                    }
                }
            }
            if status.failed() {
                return Err(sw_batch_failed_after_oid(&notes, batch.clone(), heaps));
            }
        }
        Ok(())
    }

    fn heap_snapshot(heaps: &[BlastCompoHeap], limit: usize) -> Vec<Vec<(u32, i32, u64)>> {
        let mut copy = heaps.to_vec();
        fill_results_from_compo_heaps(limit, &mut copy)
            .into_iter()
            .map(|list| {
                list.map_or_else(Vec::new, |hl| {
                    hl.hsplist_array
                        .iter()
                        .take(hl.hsplist_count)
                        .map(|m| (m.oid, m.hsps[0].raw_score, m.best_evalue.to_bits()))
                        .collect()
                })
            })
            .collect()
    }

    fn err_str(res: &Result<()>) -> Option<String> {
        res.as_ref().err().map(|e| format!("{e:#}"))
    }

    // Randomized stress test: compares `sw_c0_run_windows` against `one_thread_model` across 1,500
    // random multi-batch searches with interleaved OID order, non-monotone preliminary E-values,
    // excluded marks, speculative deferred heap tests (turned vs kept), failures, endless redos,
    // and tiny window/at-risk limits.
    #[test]
    fn c0_windows_match_one_thread_model_on_random_interleaved_searches() {
        let mut rng = XorShift(0x9e37_79b9_7f4a_7c15);
        let mut final_rng = XorShift(0x5851_f42d_4c95_7f2d);
        let prelim = [1e-30, 1e-5, 0.001, 0.005, 0.02, 1.0];
        let first_round = [1e-30, 1e-5, ECUTOFF, 0.004, 0.01, 0.5, 2.0];
        let (mut turned, mut continued, mut ended) = (0usize, 0usize, 0usize);
        let (mut batch_failures, mut endless, mut clean) = (0usize, 0usize, 0usize);
        let notes = FetchNotes::none();

        for run in 0..1500 {
            let queries = 1 + rng.below(7);
            let mut batches = Vec::new();
            let mut q = 0;
            while q < queries {
                let end = (q + 1 + rng.below(3)).min(queries);
                batches.push(q..end);
                q = end;
            }
            let heap_threshold = [1, 2, 3, 5, 20][rng.below(5)];
            // Map ranks 0..30 to scattered ordinals so excluded ordinals interleave between them.
            let ordinals: Vec<u32> = (0..30).map(|i| (i * 2) as u32).collect();
            let mut lists = Vec::new();
            let mut cases = Vec::new();
            for q_idx in 0..queries {
                if rng.below(6) == 0 {
                    lists.push(None);
                    cases.push(Vec::new());
                    continue;
                }
                let mut chosen_oids: Vec<u32> = (0..20)
                    .filter(|_| rng.below(2) == 0)
                    .map(|i| i as u32)
                    .collect();
                // Unsorted preliminary E-values per OID (and query's list sorted by E-value,
                // while stream re-sorts by OID).
                let mut pairs: Vec<(u32, f64, bool)> = chosen_oids
                    .drain(..)
                    .map(|oid| (oid, prelim[rng.below(prelim.len())], rng.below(10) == 0))
                    .collect();
                pairs.sort_by(|a, b| a.1.total_cmp(&b.1));
                let mut q_cases = Vec::with_capacity(pairs.len());
                let hsp_lists: Vec<BlastpHspList> = pairs
                    .iter()
                    .map(|&(oid, ev, empty)| {
                        let first_evalue = first_round[rng.below(first_round.len())];
                        let first_score = 30 + rng.below(4) as i32;
                        let kind = match rng.below(100) {
                            0..=9 => Kind::NotAsked,
                            10..=19 => Kind::Purged,
                            20..=29 => Kind::Failed,
                            30 if run % 4 == 0 => Kind::Endless,
                            _ => Kind::Kept,
                        };
                        let final_evalue = match final_rng.below(first_round.len() + 2) {
                            pick if pick < first_round.len() => first_round[pick],
                            pick if pick == first_round.len() => first_evalue * 3.0,
                            _ => first_evalue / 7.0,
                        };
                        let final_score = first_score - 3 + final_rng.below(7) as i32;
                        q_cases.push(Case {
                            first_evalue,
                            first_score,
                            final_evalue,
                            final_score,
                            kind,
                        });
                        make_hsp_list(oid, q_idx, ev, empty)
                    })
                    .collect();
                lists.push(make_hit_list(hsp_lists));
                cases.push(q_cases);
            }
            let mut excluded = Vec::new();
            if run % 2 == 1 {
                for q_idx in 0..queries {
                    let count = rng.below(3);
                    for _ in 0..count {
                        // Odd ordinals sit strictly between kept ordinals (even).
                        let ord = (rng.below(20) * 2 + 1) as u32;
                        excluded.push(FetchExcludedMatch {
                            query: q_idx,
                            position: 0,
                            ordinal: ord,
                            best_evalue: prelim[rng.below(prelim.len())],
                        });
                    }
                }
            }
            let win_m = 1 + rng.below(12);
            let limits = par::TestWindowLimits {
                matches: win_m,
                fill_matches: (1 + rng.below(6)).min(win_m),
                small_heap_matches: (1 + rng.below(10)).min(win_m),
                at_risk: 1 + rng.below(4),
            };
            let mut heaps: Vec<BlastCompoHeap> = (0..queries)
                .map(|_| BlastCompoHeap::new(heap_threshold, ECUTOFF))
                .collect();
            let mut lists_copy = lists.clone();
            let mut compact: Vec<Option<BlastpHitList<XPrelimHsp>>> = vec![None; queries];
            let res = par::test_run_windows(
                &mut lists_copy,
                &mut compact,
                &mut heaps,
                &batches,
                &excluded,
                Some(&ordinals),
                &notes,
                limits,
                false,
                &mut |lists, heaps, window: &par::TestWindow| {
                    window
                        .jobs
                        .iter()
                        .map(|&at| {
                            let entry = window.entries[at];
                            let SwC0Kind::Kept { position } = entry.kind else {
                                unreachable!()
                            };
                            let local_match = local_match_of(lists, entry.q_idx, position);
                            let case = cases[entry.q_idx][position];
                            let mut deferred = None;
                            let result =
                                eval_case(case, local_match.oid, entry.q_idx, &mut |ev, sc| {
                                    let ins = x_sw_heap_test(
                                        &heaps[entry.q_idx],
                                        window.before[at],
                                        (ev, sc, local_match.oid),
                                        &mut deferred,
                                    );
                                    SwHeapAnswer::Known(if deferred.is_some() { true } else { ins })
                                });
                            par::TestComputed { result, deferred }
                        })
                        .collect()
                },
            );
            let mut model_heaps: Vec<BlastCompoHeap> = (0..queries)
                .map(|_| BlastCompoHeap::new(heap_threshold, ECUTOFF))
                .collect();
            let expected = one_thread_model(
                &lists,
                &cases,
                &mut model_heaps,
                &batches,
                &excluded,
                Some(&ordinals),
            );
            let win_unit = res.as_ref().map(|_| ());
            assert_eq!(
                err_str(&win_unit.map_err(|e| anyhow::anyhow!("{e:#}"))),
                err_str(&expected),
                "run {run} limits {limits:?}"
            );
            assert_eq!(
                heap_snapshot(&heaps, heap_threshold),
                heap_snapshot(&model_heaps, heap_threshold),
                "run {run} heaps differ"
            );
            match res {
                Ok(counts) => {
                    clean += 1;
                    turned += counts.turned;
                    continued += counts.continued;
                    ended += counts.skipped + counts.wasted;
                }
                Err(err) if format!("{err:#}").contains("does not end") => endless += 1,
                Err(_) => batch_failures += 1,
            }
        }
        assert!(
            clean > 200
                && batch_failures > 100
                && endless > 20
                && turned > 100
                && continued > 100
                && ended > 50,
            "clean={clean} fail={batch_failures} endless={endless} turned={turned} cont={continued} ended={ended}"
        );
    }
}
