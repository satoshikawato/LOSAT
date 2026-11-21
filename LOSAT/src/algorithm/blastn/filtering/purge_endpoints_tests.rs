use super::*;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/gapinfo.h:45-54
fn edit_script_consumption(edit_ops: &[GapEditOp]) -> (usize, usize) {
    edit_ops.iter().fold(
        (0usize, 0usize),
        |(query_consumed, subject_consumed), op| match op {
            GapEditOp::Sub(count) => (
                query_consumed + *count as usize,
                subject_consumed + *count as usize,
            ),
            GapEditOp::Del(count) => (query_consumed, subject_consumed + *count as usize),
            GapEditOp::Ins(count) => (query_consumed + *count as usize, subject_consumed),
        },
    )
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-148
fn opaque_equal_endpoint_hsps() -> (BlastnHsp, BlastnHsp) {
    let g_ops = vec![
        GapEditOp::Sub(10),
        GapEditOp::Ins(1),
        GapEditOp::Sub(1),
        GapEditOp::Del(1),
        GapEditOp::Sub(10),
    ];
    let s_ops = vec![GapEditOp::Sub(22)];
    let opaque_query = [1u8; 22];
    let opaque_subject = [1u8; 22];
    let g_stats = stats_from_edit_ops(&opaque_query, &opaque_subject, 0, 0, &g_ops);
    let s_stats = stats_from_edit_ops(&opaque_query, &opaque_subject, 0, 0, &s_ops);

    let make_hsp = |gap_info: Vec<GapEditOp>,
                    (matches, mismatches, gap_opens, _gap_letters, align_length): (
        usize,
        usize,
        usize,
        usize,
        usize,
    )| BlastnHsp {
        identity: matches as f64 * 100.0 / align_length as f64,
        length: align_length,
        mismatch: mismatches,
        gapopen: gap_opens,
        q_start: 6,
        q_end: 27,
        s_start: 12,
        s_end: 33,
        e_value: 1.0e-9,
        bit_score: 31.0,
        num_ident: matches,
        query_frame: 1,
        query_length: 64,
        q_idx: 2,
        s_idx: 3,
        raw_score: 42,
        internal_q_offset_0: 5,
        internal_q_end_0: 27,
        internal_s_offset_0: 11,
        internal_s_end_0: 33,
        internal_query_context_offset: 0,
        gap_info: Some(gap_info),
        num_positives: matches,
    };

    (make_hsp(g_ops, g_stats), make_hsp(s_ops, s_stats))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2480-2535
// Repeated removals in both passes must keep the active order and the
// reverse removal order of live trimmed tail entries around NULL slots.
#[test]
fn common_endpoint_shifts_preserve_trimmed_tail_and_active_boundary() {
    let make = |start: usize, end: usize, score: i32| {
        let (_, mut h) = opaque_equal_endpoint_hsps();
        h.query_frame = 1;
        h.query_length = 128;
        h.q_idx = 0;
        h.s_idx = 0;
        h.raw_score = score;
        h.internal_q_offset_0 = start as i32;
        h.internal_s_offset_0 = start;
        h.internal_q_end_0 = end as i32;
        h.internal_s_end_0 = end;
        h.q_start = start + 1;
        h.s_start = start + 1;
        h.q_end = end;
        h.s_end = end;
        h.gap_info = Some(vec![GapEditOp::Sub((end - start) as u32)]);
        h
    };
    let input = vec![
        make(40, 50, 80),
        make(0, 20, 90),
        make(60, 70, 50),
        make(0, 8, 80),
        make(30, 50, 100),
        make(0, 10, 100),
        make(20, 50, 90),
    ];
    for purge in [false, true] {
        let (result, active) = purge_hsps_for_subject_ex(input.clone(), purge);
        assert_eq!(active, 3);
        let spans: Vec<_> = result
            .iter()
            .map(|h| {
                assert_eq!(h.internal_q_offset_0 as usize, h.internal_s_offset_0);
                assert_eq!(h.internal_q_end_0 as usize, h.internal_s_end_0);
                assert_eq!(
                    h.gap_info,
                    Some(vec![GapEditOp::Sub(
                        (h.internal_q_end_0 - h.internal_q_offset_0) as u32
                    )])
                );
                (h.internal_q_offset_0, h.internal_q_end_0)
            })
            .collect();
        let mut expected = vec![(0, 10), (30, 50), (60, 70)];
        if !purge {
            expected.extend([(20, 30), (10, 20)]);
        }
        assert_eq!(spans, expected);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478,2504,2268-2387
#[test]
fn handle_sort_keeps_payload_addresses_equal_order_and_tail() {
    let (a, b) = opaque_equal_endpoint_hsps();
    let storage = vec![Some(a.clone()), Some(b.clone()), Some(a)];
    let before: Vec<_> = storage
        .iter()
        .map(|h| {
            let h = h.as_ref().unwrap();
            (h as *const BlastnHsp, format!("{h:?}"))
        })
        .collect();
    let mut handles = vec![
        NonZeroUsize::new(2),
        None,
        NonZeroUsize::new(1),
        NonZeroUsize::new(3),
    ];
    for compare in [
        compare_query_offset_hsps_for_common_endpoint as fn(&BlastnHsp, &BlastnHsp) -> Ordering,
        compare_query_end_hsps_for_common_endpoint,
    ] {
        qsort_hsp_array_prefix_by(&storage, &mut handles, 3, compare);
        assert_eq!(
            handles,
            vec![
                NonZeroUsize::new(2),
                NonZeroUsize::new(1),
                None,
                NonZeroUsize::new(3)
            ]
        );
        for (slot, (address, payload)) in storage.iter().zip(&before) {
            let h = slot.as_ref().unwrap();
            assert_eq!(h as *const BlastnHsp, *address);
            assert_eq!(format!("{h:?}"), *payload);
        }
    }
    assert_eq!(
        std::mem::size_of::<HspHandle>(),
        std::mem::size_of::<usize>()
    );
}

// No NCBI counterpart: `x_in_handle_order` gives the HSPs of the live handles in their order,
// for every order and number of live handles of a short list, in the input's allocation.
#[test]
fn handle_order_is_applied_in_the_input_allocation() {
    let (a, _) = opaque_equal_endpoint_hsps();
    let n = 5usize;
    let mut perms: Vec<Vec<usize>> = vec![Vec::new()];
    for _ in 0..n {
        let mut next = Vec::new();
        for p in &perms {
            for s in (0..n).filter(|s| !p.contains(s)) {
                let mut q = p.clone();
                q.push(s);
                next.push(q);
            }
        }
        perms = next;
    }
    assert_eq!(perms.len(), 120);
    for perm in &perms {
        for live in 0..=n {
            let hsps: Vec<BlastnHsp> = (0..n)
                .map(|i| {
                    let mut h = a.clone();
                    h.raw_score = i as i32;
                    h
                })
                .collect();
            let address = hsps.as_ptr() as usize;
            let storage: Vec<Option<BlastnHsp>> = hsps.into_iter().map(Some).collect();
            assert_eq!(storage.as_ptr() as usize, address);
            // the live handles, with a removed (None) entry after the first one
            let mut handles: Vec<HspHandle> = perm[..live]
                .iter()
                .map(|&s| NonZeroUsize::new(s + 1))
                .collect();
            handles.insert(live.min(1), None);
            let out = x_in_handle_order(storage, &handles);
            let scores: Vec<usize> = out.iter().map(|h| h.raw_score as usize).collect();
            assert_eq!(scores, perm[..live].to_vec());
            if live > 0 {
                assert_eq!(out.as_ptr() as usize, address);
            }
        }
    }
}

// No NCBI counterpart: a deleted HSP's empty slot (`None`, as the purge leaves it) that a swap
// moves forward is dropped with the slots past the live ones.
#[test]
fn handle_order_drops_an_empty_slot_moved_forward() {
    let (a, _) = opaque_equal_endpoint_hsps();
    let with_score = |score: i32| {
        let mut h = a.clone();
        h.raw_score = score;
        Some(h)
    };
    let storage = vec![None, with_score(1), with_score(2)];
    let handles = vec![NonZeroUsize::new(3), None, NonZeroUsize::new(2)];
    let out = x_in_handle_order(storage, &handles);
    assert_eq!(
        out.iter().map(|h| h.raw_score).collect::<Vec<_>>(),
        vec![2, 1]
    );
}

// No NCBI counterpart: a live handle to an emptied slot is an error, as in the copy
// `x_in_handle_order` replaced (`take` of an empty slot).
#[test]
#[should_panic(expected = "live output handle requires a unique HSP")]
fn handle_order_rejects_a_handle_to_an_empty_slot() {
    let (_, b) = opaque_equal_endpoint_hsps();
    let storage = vec![None, Some(b)];
    let handles = vec![NonZeroUsize::new(1)];
    let _ = x_in_handle_order(storage, &handles);
}

// No NCBI counterpart: a repeated live handle is an error, as in the copy `x_in_handle_order`
// replaced (`take` of an empty slot).
#[test]
#[should_panic(expected = "live output handle requires a unique HSP")]
fn handle_order_rejects_a_repeated_handle() {
    let (a, b) = opaque_equal_endpoint_hsps();
    let storage = vec![Some(a), Some(b)];
    let handles = vec![NonZeroUsize::new(2), NonZeroUsize::new(2)];
    let _ = x_in_handle_order(storage, &handles);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2490-2499,2516-2525
// Also complements unit_tests/api/blasthits_unit_test.cpp:1163-1204.
#[test]
fn final_active_slot_can_be_trimmed_or_deleted_in_either_pass() {
    let make = |start: usize, end: usize, score: i32| {
        let (_, mut h) = opaque_equal_endpoint_hsps();
        h.query_frame = 1;
        h.query_length = 128;
        h.raw_score = score;
        h.internal_q_offset_0 = start as i32;
        h.internal_s_offset_0 = start;
        h.internal_q_end_0 = end as i32;
        h.internal_s_end_0 = end;
        h.q_start = start + 1;
        h.s_start = start + 1;
        h.q_end = end;
        h.s_end = end;
        h.gap_info = Some(vec![GapEditOp::Sub((end - start) as u32)]);
        h
    };
    for (keeper, removed, trimmed) in [
        (make(0, 10, 100), make(0, 20, 90), Some((10, 20))),
        (make(0, 10, 100), make(0, 8, 90), None),
        (make(10, 20, 100), make(0, 20, 90), Some((0, 10))),
        (make(10, 20, 100), make(15, 20, 90), None),
    ] {
        for purge in [false, true] {
            let (result, active) =
                purge_hsps_for_subject_ex(vec![removed.clone(), keeper.clone()], purge);
            assert_eq!(active, 1);
            assert_eq!(format!("{:?}", result[0]), format!("{keeper:?}"));
            let expected_tail = if purge { None } else { trimmed };
            assert_eq!(result.len(), 1 + usize::from(expected_tail.is_some()));
            if let Some((start, end)) = expected_tail {
                let tail = &result[1];
                assert_eq!(
                    (
                        tail.internal_q_offset_0,
                        tail.internal_q_end_0,
                        tail.internal_s_offset_0,
                        tail.internal_s_end_0
                    ),
                    (start as i32, end as i32, start, end)
                );
                assert_eq!(
                    tail.gap_info,
                    Some(vec![GapEditOp::Sub((end - start) as u32)])
                );
            }
        }
    }
}

// Test-only post-sort injection seam for the NCBI first-survivor loop.
// It intentionally omits qsort so tests can supply either permitted order
// without depending on the host libc's ordering of comparator-equal items.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478-2500
fn purge_after_injected_common_start_sort(post_sort_order: Vec<BlastnHsp>) -> Vec<BlastnHsp> {
    let mut hsp_storage: Vec<_> = post_sort_order.into_iter().map(Some).collect();
    let mut hsp_array: Vec<HspHandle> = (0..hsp_storage.len())
        .map(|index| NonZeroUsize::new(index + 1))
        .collect();
    let mut hsp_count = hsp_array.len();
    let mut i = 0usize;
    while i < hsp_count {
        let j = 1usize;
        while i + j < hsp_count && same_common_start(&hsp_storage, &hsp_array, i, i + j) {
            hsp_count -= 1;
            let removed = hsp_array[i + j].take().unwrap();
            drop(hsp_storage[removed.get() - 1].take().unwrap());
            hsp_array.copy_within(i + j + 1..hsp_count + 1, i + j);
            hsp_array[hsp_count] = None;
        }
        i += j;
    }
    hsp_array
        .into_iter()
        .take(hsp_count)
        .flatten()
        .map(|index| hsp_storage[index.get() - 1].take().unwrap())
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2268-2387
// Both common-endpoint comparators return zero after context, endpoints,
// and score tie; neither comparator examines GapEditScript.
#[test]
fn equivalent_scripts_keep_endpoints_and_score() {
    let (g_path, s_path) = opaque_equal_endpoint_hsps();
    let g_ops = g_path
        .gap_info
        .as_deref()
        .expect("G path has an edit script");
    let s_ops = s_path
        .gap_info
        .as_deref()
        .expect("S path has an edit script");

    assert_eq!(edit_script_consumption(g_ops), (22, 22));
    assert_eq!(edit_script_consumption(s_ops), (22, 22));
    assert_eq!(g_path.internal_q_offset_0, s_path.internal_q_offset_0);
    assert_eq!(g_path.internal_q_end_0, s_path.internal_q_end_0);
    assert_eq!(g_path.internal_s_offset_0, s_path.internal_s_offset_0);
    assert_eq!(g_path.internal_s_end_0, s_path.internal_s_end_0);
    assert_eq!(g_path.raw_score, s_path.raw_score);
    assert_ne!(g_ops, s_ops);
    assert_eq!(
        (
            g_path.num_ident,
            g_path.length,
            g_path.mismatch,
            g_path.gapopen
        ),
        (21, 23, 0, 2)
    );
    assert_eq!(
        (
            s_path.num_ident,
            s_path.length,
            s_path.mismatch,
            s_path.gapopen
        ),
        (22, 22, 0, 0)
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2285-2387
#[test]
fn common_endpoint_comparator_ignores_edit_script() {
    let (g_path, s_path) = opaque_equal_endpoint_hsps();

    assert_ne!(g_path.gap_info, s_path.gap_info);
    assert_eq!(
        compare_query_offset_hsps_for_common_endpoint(&g_path, &s_path),
        Ordering::Equal
    );
    assert_eq!(
        compare_query_offset_hsps_for_common_endpoint(&s_path, &g_path),
        Ordering::Equal
    );
    assert_eq!(
        compare_query_end_hsps_for_common_endpoint(&g_path, &s_path),
        Ordering::Equal
    );
    assert_eq!(
        compare_query_end_hsps_for_common_endpoint(&s_path, &g_path),
        Ordering::Equal
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478-2500
// The common-start purge keeps hsp_array[i] and frees hsp_array[i+j].
#[test]
fn source_compatible_first_survivor_accepts_either_equal_order() {
    let (g_path, s_path) = opaque_equal_endpoint_hsps();
    let g_script = g_path.gap_info.clone();
    let s_script = s_path.gap_info.clone();

    let g_then_s = purge_after_injected_common_start_sort(vec![g_path.clone(), s_path.clone()]);
    let s_then_g = purge_after_injected_common_start_sort(vec![s_path, g_path]);

    assert_eq!(g_then_s.len(), 1);
    assert_eq!(s_then_g.len(), 1);
    assert_eq!(g_then_s[0].gap_info, g_script);
    assert_eq!(s_then_g[0].gap_info, s_script);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2488-2498
// In purge mode the non-survivor is freed; the keeper's GapEditScript is
// neither trimmed nor combined with the removed HSP.
#[test]
fn source_compatible_purge_never_synthesizes_edit_script() {
    let (g_path, s_path) = opaque_equal_endpoint_hsps();
    let exact_inputs = [
        (
            g_path.gap_info.clone(),
            g_path.num_ident,
            g_path.length,
            g_path.mismatch,
            g_path.gapopen,
        ),
        (
            s_path.gap_info.clone(),
            s_path.num_ident,
            s_path.length,
            s_path.mismatch,
            s_path.gapopen,
        ),
    ];

    let survivors = purge_after_injected_common_start_sort(vec![g_path.clone(), s_path.clone()])
        .into_iter()
        .chain(purge_after_injected_common_start_sort(vec![s_path, g_path]))
        .collect::<Vec<_>>();
    assert_eq!(survivors.len(), 2);

    for survivor in survivors {
        let observed = (
            survivor.gap_info,
            survivor.num_ident,
            survivor.length,
            survivor.mismatch,
            survivor.gapopen,
        );
        assert!(exact_inputs.contains(&observed));
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2285-2318
#[test]
fn test_purge_common_start_keeps_longer_end_on_score_tie() {
    fn make_hit(q_end: usize, s_end: usize) -> BlastnHsp {
        // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-148
        BlastnHsp {
            identity: 100.0,
            length: q_end,
            mismatch: 0,
            gapopen: 0,
            q_start: 1,
            q_end,
            s_start: 1,
            s_end,
            e_value: 0.0,
            bit_score: 0.0,
            num_ident: q_end,
            query_frame: 1,
            query_length: 100,
            q_idx: 0,
            s_idx: 0,
            raw_score: 100,
            internal_q_offset_0: 0,
            internal_q_end_0: q_end as i32,
            internal_s_offset_0: 0,
            internal_s_end_0: s_end,
            internal_query_context_offset: 0,
            gap_info: None,
            num_positives: q_end,
        }
    }

    let short = make_hit(50, 50);
    let long = make_hit(80, 80);

    let (purged, _) = purge_hsps_with_common_endpoints_ex(vec![short, long], true);
    assert_eq!(purged.len(), 1);
    assert_eq!(purged[0].q_end, 80);
    assert_eq!(purged[0].s_end, 80);
}
