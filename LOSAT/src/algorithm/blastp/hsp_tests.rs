use super::*;

#[test]
fn an_hsp_of_length_zero_sorts_after_a_longer_one_that_starts_alike() {
    // Both start at subject offset 162; the longer one ends later and so comes first.
    assert_eq!(subject_offset_end(163, 162), (162, 162));
    assert_eq!(subject_offset_end(163, 169), (162, 169));
    assert_eq!(subject_offset_end(1, 0), (0, 0));
    let empty = make_hit(16, 26, 25, 163, 162);
    let longer = make_hit(16, 26, 32, 163, 169);
    assert_eq!(score_compare_hits(&longer, &empty), Ordering::Less);
    assert_eq!(score_compare_hits(&empty, &longer), Ordering::Greater);
}

fn make_hit(score: i32, q_start: usize, q_end: usize, s_start: usize, s_end: usize) -> Hit {
    Hit {
        identity: 0.0,
        length: q_end.saturating_sub(q_start).saturating_add(1),
        mismatch: 0,
        gapopen: 0,
        q_start,
        q_end,
        s_start,
        s_end,
        e_value: 0.0,
        bit_score: 0.0,
        num_ident: 0,
        query_frame: 0,
        query_length: 0,
        q_idx: 0,
        s_idx: 0,
        raw_score: score,
        // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-143
        // BlastHSP stores query/subject BlastSeg offsets used by HSP comparators.
        sort_query_offset: 0,
        sort_query_end: 0,
        sort_subject_offset: 0,
        sort_subject_end: 0,
        has_sort_offsets: false,
        gap_info: None,
        num_positives: 0,
    }
}

fn make_blastp_hsp(
    score: i32,
    q_start: usize,
    q_end: usize,
    s_start: usize,
    s_end: usize,
) -> BlastpHsp {
    BlastpHsp::from_hit(make_hit(score, q_start, q_end, s_start, s_end))
}

#[test]
fn test_score_compare_hits_prefers_higher_score_then_offsets() {
    let left = make_hit(100, 10, 30, 20, 40);
    let right = make_hit(90, 5, 25, 10, 30);
    assert_eq!(score_compare_hits(&left, &right), Ordering::Less);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1355-1377,1379-1381
#[test]
fn test_sort_hsplist_by_score_matches_ncbi_tie_breakers() {
    let mut list = BlastpHspList {
        oid: 0,
        query_index: 0,
        hsps: vec![
            make_blastp_hsp(100, 20, 40, 30, 50),
            make_blastp_hsp(120, 15, 35, 20, 40),
            make_blastp_hsp(100, 10, 25, 20, 35),
            make_blastp_hsp(100, 20, 30, 20, 30),
        ],
        best_evalue: 0.0,
    };

    sort_hsplist_by_score(&mut list);

    let ordered: Vec<(i32, usize, usize, usize, usize)> = list
        .hsps
        .iter()
        .map(|hsp| {
            (
                hsp.raw_score,
                hsp.s_start,
                hsp.s_end,
                hsp.q_start,
                hsp.q_end,
            )
        })
        .collect();
    assert_eq!(
        ordered,
        vec![
            (120, 20, 40, 15, 35),
            (100, 20, 35, 10, 25),
            (100, 20, 30, 20, 30),
            (100, 30, 50, 20, 40),
        ]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1414-1455
#[test]
fn test_sort_hsplist_by_evalue_uses_score_compare_as_tie_breaker() {
    let mut list = BlastpHspList {
        oid: 0,
        query_index: 0,
        hsps: vec![
            make_blastp_hsp(90, 20, 40, 30, 50),
            make_blastp_hsp(110, 15, 35, 20, 40),
            make_blastp_hsp(110, 10, 25, 20, 35),
        ],
        best_evalue: 0.0,
    };
    list.hsps[0].e_value = 1e-20;
    list.hsps[1].e_value = 1e-40;
    list.hsps[2].e_value = 1e-40;

    sort_hsplist_by_evalue(&mut list);

    let ordered: Vec<(f64, i32, usize, usize)> = list
        .hsps
        .iter()
        .map(|hsp| (hsp.e_value, hsp.raw_score, hsp.s_start, hsp.q_start))
        .collect();
    assert_eq!(
        ordered,
        vec![
            (1e-40, 110, 20, 15),
            (1e-40, 110, 20, 10),
            (1e-20, 90, 30, 20),
        ]
    );
}

/// One HSP list of the final hit list with the given (score, q_start, q_end, s_start,
/// s_end, E-value) HSPs in the given order.
fn hsp_list_of(oid: u32, hsps: &[(i32, usize, usize, usize, usize, f64)]) -> BlastpHspList {
    BlastpHspList {
        oid,
        query_index: 0,
        hsps: hsps
            .iter()
            .map(|&(score, q_start, q_end, s_start, s_end, e_value)| {
                let mut hsp = make_blastp_hsp(score, q_start, q_end, s_start, s_end);
                hsp.e_value = e_value;
                hsp
            })
            .collect(),
        best_evalue: 0.0,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1199-1203 (cited at
// `PairwiseHit::seqalign_num_ident`): each collected HSP carries its Seq-align's `num_ident`,
// the engine's identities, which `EngineCounts` keeps where the display recount replaced
// them; its keys count the queries from the run's first, also for a window of queries.
#[test]
fn collected_hits_carry_the_seqalign_identities() {
    let hit_lists: Vec<Option<BlastpHitList>> = (0..2)
        .map(|_| {
            let mut list = hsp_list_of(
                0,
                &[(100, 1, 50, 1, 50, 1e-20), (90, 60, 99, 60, 99, 1e-10)],
            );
            for (hsp, identities) in list.hsps.iter_mut().zip([40, 30]) {
                hsp.num_ident = identities;
            }
            let mut hit_list = BlastpHitList::new(10);
            hit_list.update(list);
            Some(hit_list)
        })
        .collect();
    let mut engine_counts = super::super::seqalign::EngineCounts::default();
    engine_counts.save((1, 0, 1), (28, 35));
    let counts = |first: usize, end: usize| -> Vec<(usize, usize)> {
        collect_hits_and_methods_from_hit_lists(&hit_lists[first..end], first, &engine_counts)
            .iter()
            .map(|(hit, _, seqalign_num_ident, _)| (hit.num_ident, *seqalign_num_ident))
            .collect()
    };
    assert_eq!(counts(0, 1), [(40, 40), (30, 30)]);
    assert_eq!(counts(1, 2), [(40, 40), (30, 28)]);
    assert_eq!(counts(0, 2), [(40, 40), (30, 30), (40, 40), (30, 28)]);
}

fn hsp_order(list: &BlastpHspList) -> Vec<(i32, usize, usize)> {
    list.hsps
        .iter()
        .map(|hsp| (hsp.raw_score, hsp.q_start, hsp.s_start))
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1574-1577
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1392-1395
// The two lists of the OrthoFinder-style benchmark (2026-10-13) as the redo leaves them,
// after s_HSPListNormalizeScores: sorted by the scaled scores (43971 before 43959, 48512
// before 48503), which both round to the same score, with E-values below 1.0e-180. NCBI
// prints them in ScoreCompareHSPs order, the smaller subject offset first.
#[test]
fn test_hit_list_to_seqalign_order_sorts_rounded_score_ties_below_evalue_epsilon() {
    let mut hit_list = BlastpHitList::new(10);
    hit_list.update(hsp_list_of(
        784,
        &[
            (1378, 173, 783, 65, 657, 5.0e-182),
            (1374, 92, 674, 96, 657, 9.825702395207e-182),
            (1374, 229, 797, 65, 643, 1.118498368325e-181),
            (1269, 81, 618, 113, 657, 6.24e-166),
        ],
    ));
    hit_list.update(hsp_list_of(
        1211,
        &[
            (1516, 604, 1120, 700, 1233, 6.002427993935e-191),
            (1516, 1, 514, 1, 534, 6.600259301753e-191),
        ],
    ));

    hit_list_to_seqalign_order(&mut hit_list);

    assert_eq!(hit_list.hsplist_count, 2);
    assert_eq!(
        hsp_order(&hit_list.hsplist_array[0]),
        vec![
            (1378, 173, 65),
            (1374, 229, 65),
            (1374, 92, 96),
            (1269, 81, 113)
        ]
    );
    assert_eq!(
        hsp_order(&hit_list.hsplist_array[1]),
        vec![(1516, 1, 1), (1516, 604, 700)]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1431-1434
// E-values of 1.0e-180 or more are compared as they are: the HSP with the higher scaled
// score has the smaller E-value and stays first although its subject offset is larger. An
// HSP of score 0 is in no Seq-align, and a list with only such HSPs is left out.
#[test]
fn test_hit_list_to_seqalign_order_keeps_evalue_order_and_drops_score_zero() {
    let mut hit_list = BlastpHitList::new(10);
    hit_list.update(hsp_list_of(
        0,
        &[
            (100, 10, 60, 40, 90, 1.0e-20),
            (100, 1, 50, 1, 50, 2.0e-20),
            (0, 1, 5, 1, 5, 9.0),
        ],
    ));
    hit_list.update(hsp_list_of(1, &[(0, 1, 5, 1, 5, 9.0)]));
    hit_list.update(hsp_list_of(2, &[(80, 1, 40, 1, 40, 1.0e-10)]));

    hit_list_to_seqalign_order(&mut hit_list);

    assert_eq!(hit_list.hsplist_count, 2);
    assert_eq!(
        hsp_order(&hit_list.hsplist_array[0]),
        vec![(100, 10, 40), (100, 1, 1)]
    );
    assert_eq!(hit_list.hsplist_array[1].oid, 2);
}

#[test]
fn test_purge_hits_with_common_endpoints_keeps_longer_match() {
    let short = make_hit(100, 1, 50, 1, 50);
    let long = make_hit(100, 1, 80, 1, 80);
    let purged = purge_hits_with_common_endpoints(vec![short, long]);
    assert_eq!(purged.len(), 1);
    assert_eq!(purged[0].q_end, 80);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2475-2534
#[test]
fn test_purge_hits_with_common_endpoints_ex_keeps_single_hsp_list_order_without_regrouping() {
    let left = make_hit(110, 5, 40, 50, 85);
    let right = make_hit(100, 1, 30, 10, 39);

    let (purged, extra_start) = purge_hits_with_common_endpoints_ex(vec![left, right], false);

    assert_eq!(purged.len(), 2);
    assert_eq!(extra_start, 2);
    assert_eq!(purged[0].raw_score, 100);
    assert_eq!(purged[1].raw_score, 110);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478-2527
#[test]
fn test_purge_hits_with_common_endpoints_no_duplicate_allocates_no_removal_marks() {
    let hits = vec![
        make_hit(100, 1, 10, 1, 10),
        make_hit(90, 20, 30, 20, 30),
        make_hit(80, 40, 50, 40, 50),
    ];

    let (purged, extra_start, stats) = purge_hits_for_subject_ex_with_stats(hits, true);

    assert_eq!(extra_start, 3);
    assert_eq!(purged.len(), 3);
    assert_eq!(
        purged
            .iter()
            .map(|hit| (hit.q_start, hit.q_end, hit.s_start, hit.s_end))
            .collect::<Vec<_>>(),
        vec![(1, 10, 1, 10), (20, 30, 20, 30), (40, 50, 40, 50)]
    );
    assert_eq!(stats.start.removed_count, 0);
    assert_eq!(stats.end.removed_count, 0);
    assert!(!stats.start.allocated_removal_marks);
    assert!(!stats.end.allocated_removal_marks);
    assert!(!stats.start.compacted);
    assert!(!stats.end.compacted);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2480-2500
#[test]
fn test_purge_hits_with_common_endpoints_sparse_duplicate_uses_index_marks() {
    let best = make_hit(100, 1, 80, 1, 80);
    let duplicate = make_hit(100, 1, 50, 1, 50);
    let distinct = make_hit(90, 100, 120, 100, 120);

    let (purged, extra_start, stats) =
        purge_hits_for_subject_ex_with_stats(vec![duplicate, distinct, best], true);

    assert_eq!(extra_start, 2);
    assert_eq!(purged.len(), 2);
    assert_eq!(purged[0].q_end, 80);
    assert_eq!(purged[1].q_start, 100);
    assert_eq!(stats.start.removed_count, 1);
    assert!(stats.start.allocated_removal_marks);
    assert!(!stats.start.used_bitmap);
    assert!(stats.start.compacted);
    assert_eq!(stats.end.removed_count, 0);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2488-2500
#[test]
fn test_purge_hits_with_common_endpoints_dense_duplicate_switches_to_bitmap() {
    let mut hits = Vec::new();
    for extra_len in 0..512 {
        hits.push(make_hit(100, 1, 20 + extra_len, 1, 20 + extra_len));
    }

    let (purged, extra_start, stats) = purge_hits_for_subject_ex_with_stats(hits, true);

    assert_eq!(extra_start, 1);
    assert_eq!(purged.len(), 1);
    assert_eq!(purged[0].q_end, 531);
    assert_eq!(stats.start.removed_count, 511);
    assert!(stats.start.allocated_removal_marks);
    assert!(stats.start.used_bitmap);
    assert!(stats.start.compacted);
}

#[test]
fn test_hit_list_update_and_prune_keep_best_subjects() {
    let mut hit_list = BlastpHitList::new(2);

    let mut first = BlastpHspList {
        oid: 0,
        query_index: 0,
        hsps: vec![BlastpHsp::from_hit(make_hit(50, 1, 10, 1, 10))],
        best_evalue: 0.0,
    };
    let mut second = BlastpHspList {
        oid: 1,
        query_index: 0,
        hsps: vec![BlastpHsp::from_hit(make_hit(60, 1, 10, 2, 11))],
        best_evalue: 0.0,
    };
    let mut third = BlastpHspList {
        oid: 2,
        query_index: 0,
        hsps: vec![BlastpHsp::from_hit(make_hit(40, 1, 10, 3, 12))],
        best_evalue: 0.0,
    };

    first.hsps[0].e_value = 1e-10;
    second.hsps[0].e_value = 1e-20;
    third.hsps[0].e_value = 1e-5;

    hit_list.update(first);
    hit_list.update(second);
    hit_list.update(third);
    hit_list.sort_by_evalue();
    hit_list.prune_by_size(2);

    assert_eq!(hit_list.hsplist_count, 2);
    assert_eq!(hit_list.hsplist_array[0].oid, 1);
    assert_eq!(hit_list.hsplist_array[1].oid, 0);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1976-2005
#[test]
fn test_reap_hsplist_by_evalue_preserves_ncbi_compaction_order() {
    let mut hsp_list = BlastpHspList {
        oid: 0,
        query_index: 0,
        hsps: vec![
            BlastpHsp::from_hit(make_hit(80, 1, 10, 1, 10)),
            BlastpHsp::from_hit(make_hit(70, 2, 11, 2, 11)),
            BlastpHsp::from_hit(make_hit(60, 3, 12, 3, 12)),
        ],
        best_evalue: 0.0,
    };
    hsp_list.hsps[0].e_value = 1e-20;
    hsp_list.hsps[1].e_value = 1e-2;
    hsp_list.hsps[2].e_value = 1e-10;

    reap_hsplist_by_evalue(&mut hsp_list, 1e-5);

    assert_eq!(hsp_list.hsps.len(), 2);
    assert_eq!(hsp_list.hsps[0].raw_score, 80);
    assert_eq!(hsp_list.hsps[1].raw_score, 60);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3243-3297
#[test]
fn test_hit_list_update_replaces_equal_worst_with_newer_oid() {
    let mut hit_list = BlastpHitList::new(1);

    let mut first = BlastpHspList {
        oid: 1,
        query_index: 0,
        hsps: vec![BlastpHsp::from_hit(make_hit(50, 1, 10, 1, 10))],
        best_evalue: 0.0,
    };
    let mut second = BlastpHspList {
        oid: 2,
        query_index: 0,
        hsps: vec![BlastpHsp::from_hit(make_hit(50, 1, 10, 2, 11))],
        best_evalue: 0.0,
    };

    first.hsps[0].e_value = 1e-20;
    second.hsps[0].e_value = 1e-20;

    hit_list.update(first);
    hit_list.update(second);
    hit_list.sort_by_evalue();

    assert_eq!(hit_list.hsplist_count, 1);
    assert_eq!(hit_list.hsplist_array[0].oid, 2);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:257-284
#[test]
fn test_compo_heap_would_insert_uses_ncbi_score_and_subject_tie_breakers() {
    let mut heap = BlastCompoHeap::new(1, 1e-50);
    let mut first = BlastpHspList {
        oid: 1,
        query_index: 0,
        hsps: vec![BlastpHsp::from_hit(make_hit(50, 1, 10, 1, 10))],
        best_evalue: 1e-20,
    };
    first.hsps[0].e_value = 1e-20;
    heap.insert(first, 1e-20, 50, 1);

    assert!(heap.would_insert(1e-20, 50, 2));
    assert!(!heap.would_insert(1e-20, 50, 0));
}

/// A fixed xorshift sequence below `bound`.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn xorshift_below(state: &mut u64, bound: usize) -> usize {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    (*state % bound as u64) as usize
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:251-276
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:329-400
/// `would_insert_view` answers as `would_insert` on random heaps of both forms (an array that
/// reached its threshold before any insert or test turned it into a heap, and a heap), ties of
/// keys included; and a FALSE stays FALSE after any later inserts, those the redo makes (after
/// a TRUE of `would_insert`) and any other. The windows of the parallel redo stop a
/// Smith-Waterman match at its first round on a FALSE of the heap as a window's start left it
/// (`x_redo_par.rs`), and rest on both.
#[test]
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn test_compo_heap_would_insert_view_equals_would_insert_and_a_false_stays_false() {
    const ECUTOFF: f64 = 0.002;
    // Around the cutoff, its stretch and the heap's comparison of E-values.
    let evalues = [0.0, 1e-200, 1e-30, 1e-5, ECUTOFF, 0.004, 0.01, 1.0, 10.0];
    let mut state = 0x2545_f491_4f6c_dd1d_u64;
    let list = |oid: u32, evalue: f64| {
        let mut list = BlastpHspList {
            oid,
            query_index: 0,
            hsps: vec![BlastpHsp::from_hit(make_hit(60, 1, 10, 1, 10))],
            best_evalue: evalue,
        };
        list.hsps[0].e_value = evalue;
        list
    };
    let (mut compared, mut array_form, mut false_kept) = (0, 0, 0);
    for heap_threshold in [1, 2, 3, 5, 8] {
        for _ in 0..100 {
            let mut heap = BlastCompoHeap::new(heap_threshold, ECUTOFF);
            let mut refused: Vec<(f64, i32, u32)> = Vec::new();
            for _ in 0..30 {
                for _ in 0..4 {
                    let evalue = evalues[xorshift_below(&mut state, evalues.len())];
                    let score = xorshift_below(&mut state, 5) as i32;
                    let oid = xorshift_below(&mut state, 12) as u32;
                    let view = heap.would_insert_view(evalue, score, oid);
                    assert_eq!(
                        view,
                        heap.clone().would_insert(evalue, score, oid),
                        "{evalue} {score} {oid} on {heap:?}"
                    );
                    compared += 1;
                    array_form +=
                        usize::from(heap.heap_array.is_none() && heap.n >= heap_threshold);
                    if !view {
                        refused.push((evalue, score, oid));
                    }
                }
                // Inserts within the cutoff more often, so that heaps fill to it.
                let evalue = evalues[xorshift_below(&mut state, 6)];
                let score = xorshift_below(&mut state, 5) as i32;
                let oid = xorshift_below(&mut state, 12) as u32;
                if xorshift_below(&mut state, 2) == 0 || heap.would_insert(evalue, score, oid) {
                    heap.insert(list(oid, evalue), evalue, score, oid);
                }
                for &(evalue, score, oid) in &refused {
                    assert!(
                        !heap.would_insert_view(evalue, score, oid)
                            && !heap.clone().would_insert(evalue, score, oid),
                        "{evalue} {score} {oid}: FALSE before, TRUE on {heap:?}"
                    );
                    false_kept += 1;
                }
            }
        }
    }
    assert!(
        compared > 50_000 && array_form > 500 && false_kept > 50_000,
        "compared {compared} array form {array_form} FALSE kept {false_kept}"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:391-414
#[test]
fn test_compo_heap_pop_returns_worst_record_first() {
    let mut heap = BlastCompoHeap::new(2, 1e-50);

    let mut better = BlastpHspList {
        oid: 7,
        query_index: 0,
        hsps: vec![BlastpHsp::from_hit(make_hit(70, 1, 10, 1, 10))],
        best_evalue: 1e-20,
    };
    let mut worse = BlastpHspList {
        oid: 3,
        query_index: 0,
        hsps: vec![BlastpHsp::from_hit(make_hit(60, 1, 10, 2, 11))],
        best_evalue: 1e-10,
    };
    better.hsps[0].e_value = 1e-20;
    worse.hsps[0].e_value = 1e-10;

    heap.insert(better, 1e-20, 70, 7);
    heap.insert(worse, 1e-10, 60, 3);

    let first = heap.pop().expect("worst record present");
    let second = heap.pop().expect("best record present");

    assert_eq!(first.oid, 3);
    assert_eq!(second.oid, 7);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1560-1580
// The test reads every heap it is given (NCBI: those of the queries of one query batch):
// while the heap of one query is as BlastCompo_HeapInitialize left it, no match of another
// query is ended, however full that query's own heap is. The query-parallel redo of
// blast_engine.rs relies on this.
#[test]
fn test_compo_early_termination_reads_every_heap_of_the_batch() {
    let ecutoff = 0.002;
    let significant = |oid: u32, evalue: f64| {
        let mut list = BlastpHspList {
            oid,
            query_index: 0,
            hsps: vec![BlastpHsp::from_hit(make_hit(60, 1, 10, 1, 10))],
            best_evalue: evalue,
        };
        list.hsps[0].e_value = evalue;
        list
    };
    let mut heaps = vec![
        BlastCompoHeap::new(2, ecutoff),
        BlastCompoHeap::new(2, ecutoff),
    ];
    assert!(!heaps[1].filled_to_cutoff());
    heaps[0].insert(significant(1, 1e-20), 1e-20, 60, 1);
    heaps[0].insert(significant(2, 1e-10), 1e-10, 60, 2);
    assert!(heaps[0].filled_to_cutoff());

    let weak = COMPO_HEAP_EVALUE_STRETCH * ecutoff * 2.0;
    // The first query's heap alone ends a weak match ...
    assert!(blast_compo_early_termination(weak, &heaps[..1]));
    // ... but not while the second query's heap is untouched or only partly filled.
    assert!(!blast_compo_early_termination(weak, &heaps));
    heaps[1].insert(significant(3, 1e-20), 1e-20, 60, 3);
    assert!(!blast_compo_early_termination(weak, &heaps));
    heaps[1].insert(significant(4, 1e-10), 1e-10, 60, 4);
    assert!(blast_compo_early_termination(weak, &heaps));
    // A match within EVALUE_STRETCH of the cutoff is redone whatever the heaps hold.
    assert!(!blast_compo_early_termination(
        COMPO_HEAP_EVALUE_STRETCH * ecutoff,
        &heaps
    ));
}

// No NCBI counterpart: the heap of a redo that ends no match early (unified P-values).
#[test]
fn a_heap_without_early_end_is_never_filled_to_its_cutoff() {
    let ecutoff = 0.002;
    let significant = |oid: u32, evalue: f64| {
        let mut list = BlastpHspList {
            oid,
            query_index: 0,
            hsps: vec![BlastpHsp::from_hit(make_hit(60, 1, 10, 1, 10))],
            best_evalue: evalue,
        };
        list.hsps[0].e_value = evalue;
        list
    };
    let mut heaps = vec![
        BlastCompoHeap::new(2, ecutoff),
        BlastCompoHeap::new(2, ecutoff).without_early_end(),
    ];
    for heap in &mut heaps {
        heap.insert(significant(1, 1e-20), 1e-20, 60, 1);
        heap.insert(significant(2, 1e-10), 1e-10, 60, 2);
    }
    let weak = COMPO_HEAP_EVALUE_STRETCH * ecutoff * 2.0;
    assert!(heaps[0].filled_to_cutoff() && heaps[0].may_fill(0));
    assert!(blast_compo_early_termination(weak, &heaps[..1]));
    // The same matches in a heap without an early end: it keeps them and ends nothing.
    assert_eq!(heaps[1].n, 2);
    assert!(!heaps[1].filled_to_cutoff() && !heaps[1].may_fill(5));
    assert!(!blast_compo_early_termination(weak, &heaps[1..]));
    assert!(!blast_compo_early_termination(weak, &heaps));
    // A match within the cutoff would still be inserted, as in any heap.
    assert!(heaps[1].would_insert(1e-30, 70, 3));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2494-2514
#[test]
fn test_fill_results_from_compo_heaps_builds_hitlist_via_hitlist_update() {
    let mut heaps = vec![BlastCompoHeap::new(2, 1e-50)];

    let mut first = BlastpHspList {
        oid: 1,
        query_index: 0,
        hsps: vec![BlastpHsp::from_hit(make_hit(50, 1, 10, 1, 10))],
        best_evalue: 1e-10,
    };
    let mut second = BlastpHspList {
        oid: 2,
        query_index: 0,
        hsps: vec![BlastpHsp::from_hit(make_hit(60, 1, 10, 2, 11))],
        best_evalue: 1e-20,
    };
    first.hsps[0].e_value = 1e-10;
    second.hsps[0].e_value = 1e-20;

    heaps[0].insert(first, 1e-10, 50, 1);
    heaps[0].insert(second, 1e-20, 60, 2);

    let mut hit_lists = fill_results_from_compo_heaps(2, &mut heaps);
    let hit_list = hit_lists[0].as_mut().expect("query hit list initialized");
    hit_list.sort_by_evalue();

    assert_eq!(hit_list.hsplist_count, 2);
    assert_eq!(hit_list.hsplist_array[0].oid, 2);
    assert_eq!(hit_list.hsplist_array[1].oid, 1);
}

// PHI-BLAST's occurrence (`BlastpHsp::phi_occurrence`, item PS8a) takes the padding of the struct:
// an HSP of every other search keeps its size (x86-64 and the other 64-bit targets).
#[cfg(target_pointer_width = "64")]
#[test]
fn the_phi_occurrence_keeps_the_size_of_an_hsp() {
    assert_eq!(std::mem::size_of::<BlastpHsp>(), 168);
}
