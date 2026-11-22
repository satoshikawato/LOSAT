use super::*;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:608-633:
// if(sum0>=max0) { max0=sum0; best[0]=H; }
// if(sum1>=max1) { max1=sum1; best[1]=H; }
// Reconciled T01 list scan, retained verbatim as the tree's reference.
fn reference_maxima(
    pool_hsp_links: &[HspLink],
    active_head: usize,
    cutoffs: [i32; 2],
    ignore_small_gaps: bool,
) -> ([Option<usize>; 2], [i64; 2]) {
    let mut best = [None; 2];
    let mut best_sum = [-i64::from(cutoffs[0]), -i64::from(cutoffs[1])];
    if !ignore_small_gaps {
        let mut cur = active_head;
        while cur != SENTINEL_IDX {
            // NCBI line 613: if(sum0>=max0)
            if pool_hsp_links[cur].sum[0] >= best_sum[0] {
                best_sum[0] = pool_hsp_links[cur].sum[0];
                best[0] = Some(cur);
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:618-623:
            // if(sum1>=max1) { max1=sum1; best[1]=H; }
            // Both independent maxima are complete before path validation.
            if pool_hsp_links[cur].sum[1] >= best_sum[1] {
                best_sum[1] = pool_hsp_links[cur].sum[1];
                best[1] = Some(cur);
            }
            cur = pool_hsp_links[cur].next_active;
        }
    } else {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:624-633:
        //    if(sum>=maxscore) { maxscore=sum; best[1]=H; }
        let mut cur = active_head;
        while cur != SENTINEL_IDX {
            // NCBI line 618: if(sum1>=max1)
            if pool_hsp_links[cur].sum[1] >= best_sum[1] {
                best_sum[1] = pool_hsp_links[cur].sum[1];
                best[1] = Some(cur);
            }
            cur = pool_hsp_links[cur].next_active;
        }
    }

    (best, best_sum)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:610-633,755-758,975-978:
// for (H=hp_start->next; H!=NULL; H=H->next) { /* stored sums */ }
// Finite state fixture with stable leaves and a separately maintained list.
fn tree_links(scores: &[[i64; 2]], active: &[bool]) -> (Vec<HspLink>, usize) {
    let mut links: Vec<_> = scores
        .iter()
        .map(|&sum| HspLink {
            score: 0,
            ctx_idx: 0,
            q_off_trim: 0,
            s_off_trim: 0,
            q_end_trim: 0,
            s_end_trim: 0,
            sum,
            xsum: [0.0; 2],
            num: [1; 2],
            link: [SENTINEL_IDX; 2],
            changed: false,
            linked_to: 0,
            start_of_chain: false,
            linked_set: false,
            next_active: SENTINEL_IDX,
            prev_active: SENTINEL_IDX,
        })
        .collect();
    let mut head = SENTINEL_IDX;
    let mut previous = SENTINEL_IDX;
    for (index, &enabled) in active.iter().enumerate() {
        if enabled {
            if previous == SENTINEL_IDX {
                head = index;
            } else {
                links[previous].next_active = index;
            }
            links[index].prev_active = previous;
            previous = index;
        }
    }
    (links, head)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:608-633:
// if(sum0>=max0) { max0=sum0; best[0]=H; }
// if(sum1>=max1) { max1=sum1; best[1]=H; }
#[test]
fn dual_maximum_boundaries_and_channel_ties() {
    let cases: Vec<(Vec<[i64; 2]>, Vec<bool>, [Option<usize>; 2], [i64; 2])> = vec![
        (vec![], vec![], [None, None], [-10, -20]),
        (vec![[4, 5]], vec![true], [Some(0), Some(0)], [4, 5]),
        (vec![[-10, -20]], vec![true], [Some(0), Some(0)], [-10, -20]),
        (vec![[-11, -21]], vec![true], [None, None], [-10, -20]),
        (
            vec![[9, 1], [2, 8]],
            vec![true, true],
            [Some(0), Some(1)],
            [9, 8],
        ),
        (
            vec![[9, 9], [9, 9], [9, 9]],
            vec![true; 3],
            [Some(2), Some(2)],
            [9, 9],
        ),
        (
            vec![[99, 99], [7, 8], [99, 99]],
            vec![false, true, false],
            [Some(1), Some(1)],
            [7, 8],
        ),
        (vec![[99, 99]], vec![false], [None, None], [-10, -20]),
    ];
    for (scores, active, expected_best, expected_sum) in cases {
        let (links, head) = tree_links(&scores, &active);
        let mut tree = DualMaximumTree::new(scores.len());
        tree.rebuild(&links, head);
        assert_eq!(tree.maxima([10, 20], false), (expected_best, expected_sum));
        for ignore in [false, true] {
            assert_eq!(
                tree.maxima([10, 20], ignore),
                reference_maxima(&links, head, [10, 20], ignore)
            );
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:907-908,936-937,965-978:
// best[0]->hsp_link.sum[0] += (best[0]->hsp_link.num[0])*cutoff[0];
#[test]
fn dual_maximum_surviving_addback_and_multi_remove() {
    let mut scores = vec![[9, 1], [2, 12], [3, 11], [1, 4], [0, 2]];
    let mut active = vec![true; 5];
    let (links, head) = tree_links(&scores, &active);
    let mut tree = DualMaximumTree::new(5);
    tree.rebuild(&links, head);
    assert_eq!(tree.capacity, 8);
    // Small winner survives selection/removal of the independent large chain.
    scores[0][0] += 10;
    tree.update(0, DualMaximum::active(0, scores[0]));
    for index in [1, 2, 4] {
        active[index] = false;
        tree.update(index, DualMaximum::INACTIVE);
        let (links, head) = tree_links(&scores, &active);
        assert_eq!(
            tree.maxima([10, 20], false),
            reference_maxima(&links, head, [10, 20], false)
        );
    }
    assert_eq!(tree.maxima([10, 20], false), ([Some(0), Some(3)], [19, 4]));
    // Explicit increase, decrease, nonmaximum and maximum removals.
    for (index, sum) in [(3, [30, 40]), (3, [-100, -100]), (0, [-10, -20])] {
        scores[index] = sum;
        tree.update(index, DualMaximum::active(index, sum));
        let (links, head) = tree_links(&scores, &active);
        assert_eq!(
            tree.maxima([10, 20], false),
            reference_maxima(&links, head, [10, 20], false)
        );
    }
    for index in [3, 0] {
        tree.update(index, DualMaximum::INACTIVE);
    }
    assert_eq!(tree.maxima([10, 20], false), ([None, None], [-10, -20]));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:659,755-758,869-873:
// if(!use_current_max) { ... H->hsp_link.sum[index] = new_sum; ... }
#[test]
fn dual_maximum_rebuild_replaces_unrepaired_state() {
    let mut tree = DualMaximumTree::new(5);
    for (scores, active) in [
        (
            vec![[99, 1], [2, 80], [70, 70], [8, 9], [6, 7]],
            vec![true; 5],
        ),
        (
            vec![[-40, -40], [3, 4], [90, 90], [11, 22], [33, 12]],
            vec![true, false, false, true, true],
        ),
        (vec![[2, 2]; 5], vec![false, true, true, false, false]),
        (vec![[100, 100]; 5], vec![false; 5]),
    ] {
        // Change only reference state: no point updates before the rebuild.
        let (links, head) = tree_links(&scores, &active);
        tree.rebuild(&links, head);
        for ignore in [false, true] {
            assert_eq!(
                tree.maxima([10, 20], ignore),
                reference_maxima(&links, head, [10, 20], ignore)
            );
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:659,755-758,869-873,975-978:
// if(!use_current_max) { ... H->hsp_link.sum[index] = new_sum; ... }
// Deterministic finite transitions, including rebuilt arbitrary stored state.
#[test]
fn dual_maximum_deterministic_state_transitions() {
    let mut seed = 0x0911_2026_u32;
    for size in [0, 1, 2, 3, 7, 8, 9, 31, 65] {
        let mut scores = vec![[0; 2]; size];
        let mut active = vec![true; size];
        let mut tree = DualMaximumTree::new(size);
        let allocation = tree.nodes.as_ptr();
        for round in 0..160 {
            for index in 0..size {
                seed = seed.wrapping_mul(1664525).wrapping_add(1013904223);
                if round % 11 == 0 || seed.is_multiple_of(7) {
                    scores[index] = [(seed % 101) as i64 - 50, ((seed >> 8) % 101) as i64 - 50];
                    active[index] = seed & 8 != 0;
                    tree.update(
                        index,
                        if active[index] {
                            DualMaximum::active(index, scores[index])
                        } else {
                            DualMaximum::INACTIVE
                        },
                    );
                }
            }
            let (links, head) = tree_links(&scores, &active);
            if round % 11 == 0 {
                tree.rebuild(&links, head);
            }
            for cutoffs in [[0, 0], [10, 20], [-10, -20], [i32::MAX, i32::MAX]] {
                for ignore in [false, true] {
                    assert_eq!(
                        tree.maxima(cutoffs, ignore),
                        reference_maxima(&links, head, cutoffs, ignore),
                        "size={size} round={round}"
                    );
                }
            }
            assert_eq!(allocation, tree.nodes.as_ptr());
        }
    }
}

/// Helper to create a mock UngappedHit with specified coordinates
pub(super) fn mock_hit(
    q_aa_start: usize,
    q_aa_end: usize,
    s_aa_start: usize,
    s_aa_end: usize,
    id: i32,
) -> UngappedHit {
    UngappedHit {
        q_idx: 0,
        s_idx: 0,
        ctx_idx: 0,
        s_f_idx: 0,
        q_frame: 1,
        s_frame: 1,
        q_aa_start: u32::try_from(q_aa_start).unwrap(),
        q_aa_end: u32::try_from(q_aa_end).unwrap(),
        s_aa_start: u32::try_from(s_aa_start).unwrap(),
        s_aa_end: u32::try_from(s_aa_end).unwrap(),
        q_seed_off: u32::try_from(q_aa_start).unwrap(),
        s_seed_off: u32::try_from(s_aa_start).unwrap(),
        q_orig_len: 1000,
        s_orig_len: 1000,
        raw_score: id, // Use raw_score as ID for verification
        e_value: 1e-10,
        num_ident: 0, // Mock value for tests
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1379-1381
        hsp_list_order: u32::try_from(id).unwrap(),
        ordering_method: 0,
        linked_set: false,
        start_of_chain: false,
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:63-65
        link_id: 0,
        chain_next_link_id: None,
        hsp_link_num: 0,
        num: 0,
    }
}

#[test]
fn context_groups_are_query_and_strand_whatever_the_frames_before() {
    // A 3-nt query has only frames +1 and -1 in LOSAT (`ctx_idx` 0 and 1); the next
    // query's frames are `ctx_idx` 2..7. NCBI's groups (`context / 3`, six contexts per
    // query) are 2 * query + strand.
    let frames = [1i8, 2, 3, -1, -2, -3];
    for (offset, &frame) in frames.iter().enumerate() {
        let mut hit = mock_hit(0, 10, 0, 10, 1);
        hit.q_idx = 1;
        hit.q_frame = frame;
        hit.ctx_idx = u32::try_from(2 + offset).unwrap();
        assert_eq!(
            translated_context_group(&hit),
            if frame > 0 { 2 } else { 3 },
            "{frame}"
        );
    }
}

#[test]
fn test_replay_ncbi_output_list_inserts_chain_members_after_head() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1018-1020
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1046-1056
    let mut member = mock_hit(1, 10, 1, 10, 10);
    member.link_id = 0;
    member.linked_set = true;
    member.start_of_chain = false;

    let mut head = mock_hit(20, 30, 20, 30, 20);
    head.link_id = 1;
    head.linked_set = true;
    head.start_of_chain = true;
    head.chain_next_link_id = Some(member.link_id);

    let mut single = mock_hit(40, 50, 40, 50, 30);
    single.link_id = 2;
    single.start_of_chain = true;

    let replayed = replay_ncbi_output_list(vec![&mut member, &mut head, &mut single]);
    let replayed_ids: Vec<usize> = replayed.iter().map(|hit| hit.link_id as usize).collect();
    assert_eq!(replayed_ids, vec![1, 0, 2]);
}

#[test]
fn replay_gives_every_chain_hsp_the_head_count() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1044-1061
    let mut member = mock_hit(1, 10, 1, 10, 10);
    member.link_id = 0;
    member.linked_set = true;
    member.num = 1;
    let mut head = mock_hit(20, 30, 20, 30, 20);
    head.link_id = 1;
    head.linked_set = true;
    head.start_of_chain = true;
    head.chain_next_link_id = Some(member.link_id);
    head.hsp_link_num = 2;
    head.num = 1;
    let mut single = mock_hit(40, 50, 40, 50, 30);
    single.link_id = 2;
    single.start_of_chain = true;
    single.hsp_link_num = 7;
    single.num = 1;

    let replayed = replay_ncbi_output_list(vec![&mut member, &mut head, &mut single]);
    let nums: Vec<(usize, i32)> = replayed
        .iter()
        .map(|hit| (hit.link_id as usize, hit.num))
        .collect();
    // The chain (head, member) has 2 HSPs; an HSP without a link keeps 1.
    assert_eq!(nums, vec![(1, 2), (0, 2), (2, 1)]);
}

#[test]
fn test_link_sort_tie_preserves_current_hsp_list_order_for_equal_hits() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:483-486
    let mut later = mock_hit(100, 150, 200, 250, 20);
    later.hsp_list_order = 2;
    let mut earlier = mock_hit(100, 150, 200, 250, 10);
    earlier.hsp_list_order = 1;
    let mut hits = vec![later, earlier];

    sort_hsps_by_ncbi_link_order(&mut hits, rev_compare_hsps_tbx);

    let ordered: Vec<usize> = hits.iter().map(|hit| hit.hsp_list_order as usize).collect();
    assert_eq!(ordered, vec![2, 1]);
}

/// Verify LOSAT's HSP sort order matches NCBI's s_RevCompareHSPsTbx
///
/// NCBI comparison logic (c++/src/algo/blast/core/link_hsps.c:359-375):
///
/// C qsort: compare(a,b) > 0 means "a comes after b"
/// So returning 1 when `h1` is below `h2` means smaller values come AFTER = DESCENDING order
#[test]
fn test_hsp_sort_order_matches_ncbi() {
    // Test data: hits with varying coordinates
    // ID = raw_score for identification after sort
    //
    // Expected order after NCBI-style sort (all descending):
    // 1. q_off=100, q_end=150, s_off=200, s_end=250 (ID=1) - largest q_off/q_end/s_off/s_end
    // 2. q_off=100, q_end=150, s_off=150, s_end=200 (ID=2) - same q_*, smaller s_off
    // 3. q_off=100, q_end=140, s_off=200, s_end=250 (ID=3) - same q_off, smaller q_end
    // 4. q_off=50,  q_end=100, s_off=300, s_end=350 (ID=4) - smallest q_off
    let mut hits = [
        mock_hit(50, 100, 300, 350, 4),  // Should be 4th (smallest q_off)
        mock_hit(100, 140, 200, 250, 3), // Should be 3rd (smaller q_end)
        mock_hit(100, 150, 200, 250, 1), // Should be 1st (largest all)
        mock_hit(100, 150, 150, 200, 2), // Should be 2nd (smaller s_off)
    ];

    // NCBI parity: Use frame-relative coordinates (not concatenated absolute)
    #[inline]
    fn frame_relative_coords(hit: &UngappedHit) -> (i32, i32, i32, i32) {
        (
            hit.q_aa_start as i32,
            hit.q_aa_end as i32,
            hit.s_aa_start as i32,
            hit.s_aa_end as i32,
        )
    }

    // LOSAT sort: all fields use b.cmp(&a) = descending
    hits.sort_by(|a, b| {
        let (aqo, aqe, aso, ase) = frame_relative_coords(a);
        let (bqo, bqe, bso, bse) = frame_relative_coords(b);
        bqo.cmp(&aqo)
            .then(bqe.cmp(&aqe))
            .then(bso.cmp(&aso))
            .then(bse.cmp(&ase))
    });

    // Verify order by checking raw_score (used as ID)
    let sorted_ids: Vec<i32> = hits.iter().map(|h| h.raw_score).collect();
    let expected_ids = vec![1, 2, 3, 4];

    println!("LOSAT sorted order (by ID): {:?}", sorted_ids);
    println!("Expected NCBI order:        {:?}", expected_ids);

    assert_eq!(
        sorted_ids, expected_ids,
        "HSP sort order mismatch! LOSAT={:?}, NCBI expected={:?}",
        sorted_ids, expected_ids
    );
}

/// Test that ties are handled correctly (stable sort behavior)
#[test]
fn test_hsp_sort_identical_coords() {
    // Two hits with identical coordinates but different IDs
    let mut hits = [
        mock_hit(100, 150, 200, 250, 1),
        mock_hit(100, 150, 200, 250, 2),
    ];

    // NCBI parity: Use frame-relative coordinates (not concatenated absolute)
    #[inline]
    fn frame_relative_coords(hit: &UngappedHit) -> (i32, i32, i32, i32) {
        (
            hit.q_aa_start as i32,
            hit.q_aa_end as i32,
            hit.s_aa_start as i32,
            hit.s_aa_end as i32,
        )
    }

    hits.sort_by(|a, b| {
        let (aqo, aqe, aso, ase) = frame_relative_coords(a);
        let (bqo, bqe, bso, bse) = frame_relative_coords(b);
        bqo.cmp(&aqo)
            .then(bqe.cmp(&aqe))
            .then(bso.cmp(&aso))
            .then(bse.cmp(&ase))
    });

    // Sort should be stable - original order preserved for equal elements
    let sorted_ids: Vec<i32> = hits.iter().map(|h| h.raw_score).collect();
    // Rust's sort_by is stable, so [1, 2] should remain [1, 2]
    assert_eq!(
        sorted_ids,
        vec![1, 2],
        "Stable sort should preserve order for equal elements"
    );
}

/// Verify NCBI comparison semantics: all fields use DESCENDING order
/// This documents the correct interpretation of s_RevCompareHSPsTbx
#[test]
fn test_ncbi_comparison_semantics() {
    // NCBI s_RevCompareHSPsTbx (c++/src/algo/blast/core/link_hsps.c:367-370):
    //
    // In C qsort: positive return means first arg comes AFTER second arg
    // So: h1.offset < h2.offset → h1 after h2 → smaller comes after → DESCENDING
    //
    // The status document incorrectly stated this was ASCENDING.
    // This test documents the correct interpretation.

    // Test: same query coords, different subject offsets
    // NCBI descending: larger s_off should come first
    let mut hits = [
        mock_hit(100, 150, 50, 100, 1),  // s_off=50 (smaller)
        mock_hit(100, 150, 200, 250, 2), // s_off=200 (larger)
    ];

    // NCBI parity: Use frame-relative coordinates (not concatenated absolute)
    #[inline]
    fn frame_relative_coords(hit: &UngappedHit) -> (i32, i32, i32, i32) {
        (
            hit.q_aa_start as i32,
            hit.q_aa_end as i32,
            hit.s_aa_start as i32,
            hit.s_aa_end as i32,
        )
    }

    hits.sort_by(|a, b| {
        let (aqo, aqe, aso, ase) = frame_relative_coords(a);
        let (bqo, bqe, bso, bse) = frame_relative_coords(b);
        bqo.cmp(&aqo)
            .then(bqe.cmp(&aqe))
            .then(bso.cmp(&aso)) // LOSAT: b.cmp(&a) = DESCENDING
            .then(bse.cmp(&ase))
    });

    let sorted_ids: Vec<i32> = hits.iter().map(|h| h.raw_score).collect();

    // DESCENDING: larger s_off (ID=2) should come first
    assert_eq!(
        sorted_ids,
        vec![2, 1],
        "Subject offset should be sorted DESCENDING (larger first). Got: {:?}",
        sorted_ids
    );

    println!("✓ Confirmed: LOSAT uses DESCENDING order for subject.offset (matches NCBI)");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:478-486,990-994
#[test]
fn test_final_link_sorts_preserve_payloads_addresses_and_all_ties() {
    use crate::algorithm::tblastx::ungapped_hit_sort::sort_ungapped_hits_by_index_replay;
    let mut hits: Vec<_> = (0..96)
        .map(|i| {
            let mut hit = mock_hit(
                (i / 3 % 8) * 20,
                (i / 3 % 8) * 20 + 15,
                (i / 3 % 5) * 30,
                (i / 3 % 5) * 30 + 15,
                i as i32,
            );
            hit.ctx_idx = u32::try_from(i / 24).unwrap();
            hit.s_frame = if i / 12 % 2 == 0 { -1 } else { 1 };
            hit.link_id = u32::try_from(i).unwrap();
            hit.e_value = f64::from_bits(0x3ff0000000000000 + i as u64);
            hit
        })
        .collect();
    sort_hsps_by_ncbi_link_order(&mut hits, rev_compare_hsps_tbx);
    let mut expected = hits.clone();
    let mut addresses = vec![std::ptr::null(); hits.len()];
    for hit in &hits {
        addresses[hit.link_id as usize] = hit as *const UngappedHit;
    }
    let mut references: Vec<_> = hits.iter_mut().collect();
    let comparators: [fn(&UngappedHit, &UngappedHit) -> Ordering; 2] =
        [rev_compare_hsps_transl, fwd_compare_hsps_transl];
    for compare in comparators {
        sort_ungapped_hits_by_index_replay(&mut expected, compare);
        sort_linked_hsp_references(&mut references, compare);
        for (actual, original) in references.iter().zip(&expected) {
            assert_eq!(
                &**actual as *const UngappedHit,
                addresses[actual.link_id as usize]
            );
            assert_eq!(format!("{:?}", **actual), format!("{original:?}"));
            assert_eq!(actual.e_value.to_bits(), original.e_value.to_bits());
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1018-1020,1046-1056,1080-1085
// Preserve the existing Rust replay's singleton and repeated-index cases.
#[test]
fn test_link_reference_replay_singleton_and_repeated_indices() {
    assert!(replay_ncbi_output_list(Vec::new()).is_empty());
    let mut member = mock_hit(1, 10, 1, 10, 10);
    member.link_id = 7;
    member.linked_set = true;
    member.start_of_chain = false;
    let singleton = replay_ncbi_output_list(vec![&mut member]);
    assert_eq!(singleton.len(), 1);
    assert_eq!(format!("{:?}", singleton[0]), format!("{member:?}"));

    let mut head_a = mock_hit(20, 30, 20, 30, 20);
    head_a.link_id = 1;
    head_a.start_of_chain = true;
    head_a.linked_set = true;
    head_a.chain_next_link_id = Some(7);
    let mut head_b = head_a.clone();
    head_b.link_id = 2;
    let replayed = replay_ncbi_output_list(vec![&mut member, &mut head_a, &mut head_b]);
    assert_eq!(
        replayed.iter().map(|h| h.link_id).collect::<Vec<_>>(),
        [1, 7, 2, 7]
    );
    assert_eq!(format!("{:?}", replayed[1]), format!("{:?}", replayed[3]));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558,990-994,1080-1085
#[test]
fn test_final_link_replay_parallel_matches_sequential_complete_payloads() {
    let params = KarlinParams::default();
    let contexts: Vec<_> = (0..12)
        .map(|index| QueryContext {
            q_idx: (index / 6) as u32,
            frame: [1, 2, 3, -1, -2, -3][index % 6],
            aa_seq: Vec::new(),
            aa_seq_nomask: None,
            aa_len: 2000,
            orig_len: 6000,
            frame_base: (index * 2001) as i32,
            is_valid: true,
            karlin_params: params,
        })
        .collect();
    let linking = LinkingParams {
        subject_len_nucl: 6000,
        gap_decay_rate: 0.5,
        cutoffs: LinkHspCutoffs {
            cutoff_small_gap: 20,
            cutoff_big_gap: 20,
            gap_prob: 0.5,
        },
    };
    let hits: Vec<_> = (0..96)
        .map(|index| {
            let position = (index % 12) * 30;
            let mut hit = mock_hit(
                position,
                position + 20,
                position,
                position + 20,
                40 + (index % 4) as i32,
            );
            hit.ctx_idx = u32::try_from(index / 8).unwrap();
            hit.q_idx = (hit.ctx_idx as usize / 6) as u32;
            hit.q_frame = contexts[hit.ctx_idx as usize].frame;
            hit.s_frame = if index / 4 % 2 == 0 { -1 } else { 1 };
            hit.q_orig_len = 6000;
            hit.s_orig_len = 6000;
            hit.hsp_list_order = u32::try_from(index).unwrap();
            hit
        })
        .collect();
    for count in [0, 1, hits.len()] {
        let run = |allow_parallel| {
            apply_sum_stats_even_gap_linking_with_parallel(
                hits[..count].to_vec(),
                &params,
                &linking,
                &contexts,
                &[0; 6],
                &[0; 12],
                &[4_000_000; 12],
                allow_parallel,
            )
        };
        let sequential = run(false);
        let parallel = run(true);
        assert_eq!(sequential.len(), count);
        assert_eq!(parallel.len(), sequential.len());
        for (a, b) in sequential.iter().zip(&parallel) {
            assert_eq!(format!("{a:?}"), format!("{b:?}"));
            assert_eq!(a.e_value.to_bits(), b.e_value.to_bits());
        }
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:827-863,876-885
// for (H2_index=H_index-1; H2_index>1;) { ... if(b0) H2_index=next_larger; }
// The end sentinels 0/1 cannot select a predecessor, even with tracing enabled.
#[test]
fn large_gap_scan_preserves_initial_state_and_zero_one_jumps() {
    let initial = (5, 3, -0.0f64, 9usize);
    let empty = scan_large_gap_predecessors(
        2,
        &[],
        &[],
        &[],
        10,
        10,
        false,
        initial.0,
        initial.1,
        initial.2,
        initial.3,
    );
    assert_eq!(
        (empty.0, empty.1, empty.2.to_bits(), empty.3),
        (initial.0, initial.1, initial.2.to_bits(), initial.3)
    );
    for jump in [0, 1] {
        let helper = LhHelper {
            hsp_idx: SENTINEL_IDX,
            q_off_trim: 20,
            s_off_trim: 20,
            sum: [0, 5],
            next_larger: jump,
            maxsum1: 0,
        };
        let helpers = vec![helper; 3];
        for trace in [false, true] {
            let result = scan_large_gap_predecessors(
                3,
                &helpers,
                &[],
                &[],
                10,
                10,
                trace,
                initial.0,
                initial.1,
                initial.2,
                initial.3,
            );
            assert_eq!(
                (result.0, result.1, result.2.to_bits(), result.3),
                (initial.0, initial.1, initial.2.to_bits(), initial.3)
            );
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:838-859
// if (!(b0|b1|b2)) { H2=H2_helper->ptr; H_hsp_xsum=H2->hsp_link.xsum[index]; }
// Strict ties retain the first visited HSP, using the helper's mapped index.
#[test]
fn large_gap_scan_preserves_strict_ties_mapping_and_xsum_bits() {
    let (mut links, _) = tree_links(&[[0, 10], [0, 10]], &[true, true]);
    links[0].num[1] = 7;
    links[1].num[1] = 8;
    links[0].xsum[1] = f64::from_bits(1.0f64.to_bits() + 1);
    links[1].xsum[1] = f64::from_bits(2.0f64.to_bits() + 1);
    let hits = vec![mock_hit(0, 1, 0, 1, 0), mock_hit(2, 3, 2, 3, 1)];
    let helper = LhHelper {
        hsp_idx: 1,
        q_off_trim: 20,
        s_off_trim: 20,
        sum: [0, 10],
        next_larger: 0,
        maxsum1: 0,
    };
    let mut helpers = vec![helper; 4];
    helpers[3].hsp_idx = 0;
    helpers[3].next_larger = 2;
    for (query_offset, subject_offset, expected) in [(20, 20, 0), (10, 20, 1), (20, 10, 1)] {
        helpers[3].q_off_trim = query_offset;
        helpers[3].s_off_trim = subject_offset;
        for trace in [false, true] {
            let result = scan_large_gap_predecessors(
                4,
                &helpers,
                &links,
                &hits,
                10,
                10,
                trace,
                0,
                0,
                0.0,
                SENTINEL_IDX,
            );
            assert_eq!(
                (result.0, result.1, result.2.to_bits(), result.3),
                (
                    links[expected].sum[1],
                    links[expected].num[1],
                    links[expected].xsum[1].to_bits(),
                    expected
                )
            );
        }
    }
}

// -----------------------------------------------------------------------
// Incremental kernel (`linking_incr.rs`, the default)
// -----------------------------------------------------------------------

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:414-419,589
// These tests compare the incremental kernel with the port of s_BlastEvenGapLinkHSPs on the
// same groups. The reference results are the literal port above, not hand-written numbers.
use crate::algorithm::tblastx::sum_stats_linking::linking_incr::{
    link_hsp_group_incr_with, IncrOptions, IncrStats,
};

pub(super) fn fast_test_cutoffs() -> LinkHspCutoffs {
    LinkHspCutoffs {
        cutoff_small_gap: 41,
        cutoff_big_gap: 44,
        gap_prob: 0.5,
    }
}

pub(super) fn fast_test_contexts() -> (KarlinParams, Vec<QueryContext>, Vec<f64>) {
    let params = KarlinParams {
        lambda: 0.3176,
        k: 0.134,
        h: 0.401,
        alpha: 0.7916,
        beta: -3.2,
    };
    let contexts: Vec<_> = (0..6)
        .map(|index| QueryContext {
            q_idx: 0,
            frame: [1, 2, 3, -1, -2, -3][index],
            aa_seq: Vec::new(),
            aa_seq_nomask: None,
            aa_len: 100_000,
            orig_len: 300_000,
            frame_base: (index * 100_001) as i32,
            is_valid: true,
            karlin_params: params,
        })
        .collect();
    let log_k: Vec<f64> = contexts
        .iter()
        .map(|c| c.karlin_params.k.glibc_ln())
        .collect();
    (params, contexts, log_k)
}

/// HSPs of a fixed test: (q_start, q_end, s_start, s_end, score).
fn fast_test_hits(spec: &[(usize, usize, usize, usize, i32)]) -> Vec<UngappedHit> {
    spec.iter()
        .enumerate()
        .map(|(index, &(q0, q1, s0, s1, score))| {
            let mut hit = mock_hit(q0, q1, s0, s1, score);
            hit.hsp_list_order = u32::try_from(index).unwrap();
            hit.q_orig_len = 300_000;
            hit.s_orig_len = 300_000;
            hit
        })
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:414-419,974,1049
// The helper runs the port of s_BlastEvenGapLinkHSPs and the incremental kernel on the same
// sorted group and compares every field of the result (E-values bit for bit).
/// Links `hits` (any order) with the literal port and with the incremental
/// kernel and requires the same HSPs in the same order with every field
/// equal. Returns that result and what the incremental kernel did.
fn link_with_both_kernels(
    mut hits: Vec<UngappedHit>,
    cutoffs: &LinkHspCutoffs,
    options: IncrOptions,
    label: &str,
) -> (Vec<UngappedHit>, IncrStats) {
    let (params, contexts, log_k) = fast_test_contexts();
    sort_hsps_by_ncbi_link_order(&mut hits, rev_compare_hsps_tbx);
    for (link_id, hit) in hits.iter_mut().enumerate() {
        hit.link_id = u32::try_from(link_id).unwrap();
        hit.chain_next_link_id = None;
        hit.num = 1;
    }
    let mut helpers = Vec::new();
    let mut links = Vec::new();
    let expected = link_hsp_group_ncbi(
        hits.clone(),
        &params,
        cutoffs,
        0.5,
        false,
        300_000,
        &contexts,
        &[0; 6],
        &[49; 6],
        &[9_372_428_362; 6],
        &log_k,
        &mut helpers,
        &mut links,
    );
    let (actual, stats) = link_hsp_group_incr_with(
        hits,
        cutoffs,
        0.5,
        300_000,
        &contexts,
        &[49; 6],
        &[9_372_428_362; 6],
        &log_k,
        options,
    );
    assert_eq!(expected.len(), actual.len(), "{label}");
    for (a, b) in expected.iter().zip(&actual) {
        assert_eq!(a.e_value.to_bits(), b.e_value.to_bits(), "{label}");
        assert_eq!(format!("{a:?}"), format!("{b:?}"), "{label}");
    }
    (actual, stats)
}

/// Every pass checked against a pass from scratch, with the u64 key and
/// with the u128 key.
const VERIFIED: [IncrOptions; 2] = [
    IncrOptions {
        verify: true,
        wide_keys: false,
    },
    IncrOptions {
        verify: true,
        wide_keys: true,
    },
];

/// The HSP of a fixed test, found by its query and subject start.
fn linked_hit(hits: &[UngappedHit], q_aa_start: usize, s_aa_start: usize) -> &UngappedHit {
    hits.iter()
        .find(|h| h.q_aa_start as usize == q_aa_start && h.s_aa_start as usize == s_aa_start)
        .expect("HSP of the fixed test")
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:827-861
// The NCBI scan visits only HSPs before H in list order, so J is never a candidate of H. The
// test makes the Fenwick tree of the index-1 sweep answer with J, which forces the fallback scan.
// H spans 4 residues and J starts one residue before it. J follows H in
// list order, so NCBI never offers J to H; but J's trimmed start (99 + 5)
// lies beyond H's trimmed end (104 - 1) on both sequences, so the tree of
// the sweep answers H's large-gap search with J. The kernel must then
// decide by NCBI's scan.
#[test]
fn incremental_kernel_scans_when_the_tree_answers_with_a_later_hsp() {
    let cutoffs = fast_test_cutoffs();
    for options in VERIFIED {
        // Nothing precedes H: it stays alone.
        let (linked, stats) = link_with_both_kernels(
            fast_test_hits(&[(100, 104, 200, 204, 50), (99, 129, 260, 290, 120)]),
            &cutoffs,
            options,
            "4-residue HSP without a predecessor",
        );
        assert_eq!(stats.fallbacks, 1, "{stats:?}");
        assert!(stats.verified_passes > 0, "{stats:?}");
        assert_eq!(stats.wide_groups, u64::from(options.wide_keys), "{stats:?}");
        let h = linked_hit(&linked, 100, 200);
        let j = linked_hit(&linked, 99, 260);
        assert!(!h.linked_set && h.chain_next_link_id.is_none(), "{h:?}");
        assert!(!j.linked_set && j.chain_next_link_id.is_none(), "{j:?}");

        // K precedes H and is its only candidate. Its sum (60 - 44) is
        // smaller than J's (120 - 44), so the tree still answers with J:
        // H must link to K.
        let (linked, stats) = link_with_both_kernels(
            fast_test_hits(&[
                (100, 104, 200, 204, 50),
                (99, 129, 260, 290, 120),
                (140, 170, 240, 270, 60),
            ]),
            &cutoffs,
            options,
            "4-residue HSP with one predecessor",
        );
        assert_eq!(stats.fallbacks, 1, "{stats:?}");
        let h = linked_hit(&linked, 100, 200);
        let j = linked_hit(&linked, 99, 260);
        let k = linked_hit(&linked, 140, 240);
        assert!(h.linked_set && h.start_of_chain, "{h:?}");
        assert_eq!(h.chain_next_link_id, Some(k.link_id));
        assert!(!j.linked_set, "{j:?}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:635-649,781-795
// This test reaches the recomputation pass after a removal (use_current_max = 0) and checks
// that the HSPs the kernel does not compute again, and the choices it keeps without a search,
// equal the ones NCBI computes by searching.
#[test]
fn incremental_kernel_keeps_unchanged_choices_in_a_later_pass() {
    let cutoffs = fast_test_cutoffs();
    // Removing the best chain a1 -> a2 leaves e, whose recorded chain ran
    // through a2, as the stale maximum, so NCBI recomputes
    // (link_hsps.c:639-649). Only e selected a removed HSP: the pass
    // computes e again (both methods) and leaves b2 ("no HSP") and
    // b1 -> b2 as they are, without visiting them.
    let spec = [
        (1000, 1030, 1000, 1030, 100), // a2
        (960, 990, 960, 990, 100),     // a1 -> a2
        (955, 985, 962, 992, 90),      // e  -> a2
        (3040, 3070, 140, 170, 60),    // b2
        (3000, 3030, 100, 130, 60),    // b1 -> b2
    ];
    for options in VERIFIED {
        let (linked, stats) =
            link_with_both_kernels(fast_test_hits(&spec), &cutoffs, options, "unchanged HSPs");
        assert_eq!((stats.rounds, stats.passes), (3, 2), "{stats:?}");
        assert_eq!(stats.removed, 2, "{stats:?}");
        // Five HSPs in the first pass, e in the second.
        assert_eq!((stats.visited0, stats.searched0), (6, 6), "{stats:?}");
        assert_eq!((stats.visited1, stats.searched1), (6, 6), "{stats:?}");
        assert_eq!((stats.kept0, stats.kept1), (0, 0), "{stats:?}");
        // The second sweep reads the three remaining HSPs.
        assert_eq!((stats.swept, stats.fallbacks), (3, 0), "{stats:?}");
        assert_eq!(stats.verified_passes, 2, "{stats:?}");
        let a2 = linked_hit(&linked, 1000, 1000);
        let a1 = linked_hit(&linked, 960, 960);
        let e = linked_hit(&linked, 955, 962);
        let b2 = linked_hit(&linked, 3040, 140);
        let b1 = linked_hit(&linked, 3000, 100);
        assert_eq!(a1.chain_next_link_id, Some(a2.link_id));
        assert_eq!(b1.chain_next_link_id, Some(b2.link_id));
        assert!(!e.linked_set && e.chain_next_link_id.is_none(), "{e:?}");
    }

    // P selects Q1 over Q2 (equal sums, Q1 nearer) under both methods, and
    // H selects P. X, whose only candidate is Q1, gives the best chain
    // X -> Q1; its removal leaves H -> P -> Q1 as the stale maximum. The
    // second pass computes P again: it selects Q2 at the same sum, so its
    // link changes and H is visited, but P's sum did not change and H
    // keeps P without a search under both methods.
    let spec = [
        (1110, 1140, 1110, 1140, 60),  // Q2
        (1100, 1130, 1140, 1170, 60),  // Q1
        (1050, 1080, 1075, 1105, 60),  // P  -> Q1, then Q2
        (1040, 1070, 1095, 1125, 100), // X  -> Q1
        (1025, 1055, 1020, 1050, 60),  // H  -> P
    ];
    for options in VERIFIED {
        let (linked, stats) =
            link_with_both_kernels(fast_test_hits(&spec), &cutoffs, options, "kept choices");
        assert_eq!(stats.passes, 2, "{stats:?}");
        assert_eq!((stats.kept0, stats.kept1), (1, 1), "{stats:?}");
        assert_eq!(stats.verified_passes, 2, "{stats:?}");
        let q2 = linked_hit(&linked, 1110, 1110);
        let q1 = linked_hit(&linked, 1100, 1140);
        let p = linked_hit(&linked, 1050, 1075);
        let x = linked_hit(&linked, 1040, 1095);
        let h = linked_hit(&linked, 1025, 1020);
        assert_eq!(x.chain_next_link_id, Some(q1.link_id));
        assert!(h.linked_set && h.start_of_chain, "{h:?}");
        assert_eq!(h.chain_next_link_id, Some(p.link_id));
        assert_eq!(p.chain_next_link_id, Some(q2.link_id));
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:738-745
// NCBI takes the candidate with the larger sum (p), not the nearer one (k). The test checks that
// a choice stays correct while another candidate of H loses its chain.
// H has two small-gap candidates, p (sum 69) and k (sum 5 + 59 through
// a2). The first pass selects p. Removing a1 -> a2 forces a second pass
// in which k searches again and drops to 5, while p is untouched: H is
// not computed again and keeps p, and the pass from scratch confirms it.
#[test]
fn incremental_kernel_keeps_a_choice_while_another_candidate_loses_its_chain() {
    let cutoffs = fast_test_cutoffs();
    let spec = [
        (1000, 1030, 1000, 1030, 100), // a2
        (975, 1010, 915, 950, 110),    // p
        (960, 990, 960, 990, 100),     // a1 -> a2
        (955, 985, 962, 992, 90),      // e  -> a2, the stale maximum
        (945, 975, 930, 960, 46),      // k  -> a2
        (910, 940, 880, 910, 43),      // H  -> p (k is nearer but smaller)
    ];
    for options in VERIFIED {
        let (linked, stats) = link_with_both_kernels(
            fast_test_hits(&spec),
            &cutoffs,
            options,
            "kept choice, changed rival",
        );
        assert_eq!(stats.passes, 2, "{stats:?}");
        // Six HSPs in the first pass, e and k (not H) in the second.
        assert_eq!((stats.visited0, stats.searched0), (8, 8), "{stats:?}");
        assert_eq!(stats.verified_passes, 2, "{stats:?}");
        let h = linked_hit(&linked, 910, 880);
        let p = linked_hit(&linked, 975, 915);
        let k = linked_hit(&linked, 945, 930);
        assert!(h.linked_set && h.start_of_chain, "{h:?}");
        assert_eq!(h.chain_next_link_id, Some(p.link_id));
        assert!(!k.linked_set, "{k:?}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:414-419,589
// The reference for the comparison is the port of s_BlastEvenGapLinkHSPs.
// Differential test of the incremental kernel against the literal port
// on random HSP groups built to reach the rare paths: equal sums, equal
// coordinates, HSPs of four to six residues (whose candidates can follow
// them in list order), dense clusters and long chains, with and without
// small gaps, half of them with the u128 key. Every pass of the smaller
// groups is also checked against a pass from scratch, and the cases must
// reach every branch of the kernel.
#[test]
fn incremental_kernel_matches_ncbi_kernel_on_random_groups() {
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n.max(1)
        }
    }

    let (_, contexts, _) = fast_test_contexts();
    let cutoff_sets = [
        fast_test_cutoffs(),
        LinkHspCutoffs {
            cutoff_small_gap: 20,
            cutoff_big_gap: 20,
            gap_prob: 0.5,
        },
        LinkHspCutoffs {
            cutoff_small_gap: 0,
            cutoff_big_gap: 30,
            gap_prob: 0.0,
        },
        LinkHspCutoffs {
            cutoff_small_gap: 30,
            cutoff_big_gap: 25,
            gap_prob: 0.5,
        },
    ];

    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let mut total_hsps = 0usize;
    let mut total_chained = 0usize;
    let mut total = IncrStats::default();
    let cases: usize = std::env::var("LOSAT_FUZZ_CASES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(3000);
    for case in 0..cases {
        let cutoffs = &cutoff_sets[case % cutoff_sets.len()];
        let style = (case / cutoff_sets.len()) % 6;
        let n = 1 + rng.below(if case % 50 == 0 { 900 } else { 140 }) as usize;
        let span: u64 = match style {
            0 => 3000,
            1 => 300,
            2 => 60,
            3 => 1200,
            4 => 150,
            _ => 20_000,
        };
        let base_score = cutoffs.cutoff_big_gap.max(cutoffs.cutoff_small_gap);
        let mut hits: Vec<UngappedHit> = Vec::with_capacity(n);
        for index in 0..n {
            let short = style == 4 || rng.below(9) == 0;
            let len = if short {
                4 + rng.below(3) as usize
            } else {
                7 + rng.below(if style == 3 { 200 } else { 45 }) as usize
            };
            let (q0, s0) = if style == 3 || (style == 5 && rng.below(3) > 0) {
                // near one diagonal: long chains
                let q = rng.below(span) as usize;
                let jitter = rng.below(25) as usize;
                (q, q + jitter)
            } else {
                (rng.below(span) as usize, rng.below(span) as usize)
            };
            // few distinct scores: many equal sums
            let score = base_score - 3
                + if style == 2 {
                    rng.below(6) as i32
                } else {
                    rng.below(40) as i32
                };
            let mut hit = mock_hit(q0, q0 + len, s0, s0 + len, score.max(1));
            hit.ctx_idx = u32::try_from(rng.below(3)).unwrap();
            hit.q_frame = contexts[hit.ctx_idx as usize].frame;
            hit.s_frame = 1 + rng.below(3) as i8;
            hit.q_orig_len = 300_000;
            hit.s_orig_len = 300_000;
            hit.hsp_list_order = u32::try_from(index).unwrap();
            if index > 0 && rng.below(12) == 0 {
                // the coordinates of an earlier HSP, in another frame
                let other = rng.below(index as u64) as usize;
                hit.q_aa_start = hits[other].q_aa_start;
                hit.q_aa_end = hits[other].q_aa_end;
                hit.s_aa_start = hits[other].s_aa_start;
                hit.s_aa_end = hits[other].s_aa_end;
            }
            hits.push(hit);
        }

        let label = format!("case {case} style {style} n {n}");
        // The pass from scratch of `verify` costs what the literal port
        // costs, so only the smaller groups are verified.
        let options = IncrOptions {
            verify: n <= 140,
            wide_keys: case % 2 == 1,
        };
        let (linked, stats) = link_with_both_kernels(hits, cutoffs, options, &label);
        total.add(&stats);
        total_hsps += linked.len();
        total_chained += linked.iter().filter(|h| h.linked_set).count();
    }
    // The cases must exercise chains, not only single HSPs ...
    assert!(
        total_chained * 10 > total_hsps,
        "{total_chained} of {total_hsps}"
    );
    // ... and every branch of the incremental kernel.
    if cases >= 3000 {
        let t = &total;
        assert!(t.kept0 > 0 && t.searched0 > 0, "{t:?}");
        assert!(t.kept1 > 0 && t.searched1 > 0, "{t:?}");
        assert!(t.same0 > 0 && t.same1 > 0, "{t:?}");
        assert!(t.removed > 0 && t.swept > 0 && t.fallbacks > 0, "{t:?}");
        assert!(t.verified_passes > 0 && t.wide_groups > 0, "{t:?}");
        eprintln!(
            "incremental_kernel_matches_ncbi_kernel_on_random_groups: {total_hsps} HSPs, {t:?}"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:753,868,907-908
// Groups whose sums leave Int4, where NCBI's C would wrap. LOSAT adds the sums in i64 (Owner
// decision 2026-10-10; linking_incr.rs, "Sums"), and the reference is the port of
// s_BlastEvenGapLinkHSPs, which also adds in i64.
/// Differential test on groups whose sums leave the Int4 range: the
/// literal port against the incremental kernel with either key (every
/// pass checked against a pass from scratch). Requires the same result
/// from all three, and both kinds of group (with and without a pass value
/// outside Int4).
#[test]
fn kernels_add_sums_beyond_int4() {
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n.max(1)
        }
    }

    // A small Lambda keeps xsum, and so the E-values that choose the
    // ordering method, in the usual range for raw scores near 2^30.
    let params = KarlinParams {
        lambda: 1.0e-8,
        k: 0.134,
        h: 0.401,
        alpha: 0.7916,
        beta: -3.2,
    };
    let (_, mut contexts, _) = fast_test_contexts();
    for context in contexts.iter_mut() {
        context.karlin_params = params;
    }
    let log_k: Vec<f64> = contexts
        .iter()
        .map(|c| c.karlin_params.k.glibc_ln())
        .collect();
    let link_ncbi = |hits: Vec<UngappedHit>, cutoffs: &LinkHspCutoffs| {
        let mut helpers = Vec::new();
        let mut links = Vec::new();
        link_hsp_group_ncbi(
            hits,
            &params,
            cutoffs,
            0.5,
            false,
            300_000,
            &contexts,
            &[0; 6],
            &[49; 6],
            &[9_372_428_362; 6],
            &log_k,
            &mut helpers,
            &mut links,
        )
    };
    let link_incr = |hits: Vec<UngappedHit>, cutoffs: &LinkHspCutoffs, wide_keys| {
        link_hsp_group_incr_with(
            hits,
            cutoffs,
            0.5,
            300_000,
            &contexts,
            &[49; 6],
            &[9_372_428_362; 6],
            &log_k,
            IncrOptions {
                verify: true,
                wide_keys,
            },
        )
    };
    let assert_same = |a: &[UngappedHit], b: &[UngappedHit], label: &str| {
        assert_eq!(a.len(), b.len(), "{label}");
        for (x, y) in a.iter().zip(b) {
            assert_eq!(x.e_value.to_bits(), y.e_value.to_bits(), "{label}");
            assert_eq!(format!("{x:?}"), format!("{y:?}"), "{label}");
        }
    };

    // (cutoff[0], cutoff[1], gap_prob, scores from, to); cutoff[0] == 0 ignores small gaps.
    let settings: [(i32, i32, f64, i64, i64); 6] = [
        (41, 44, 0.5, 300_000_000, 1_100_000_000),
        (41, 44, 0.5, 45, 900_000_000),
        (600_000_000, 700_000_000, 0.5, 700_000_001, 1_700_000_000),
        (0, 500_000_000, 0.0, 500_000_001, 1_500_000_000),
        (41, -300_000_000, 0.5, 200_000_000, 900_000_000),
        (900_000_000, 950_000_000, 0.5, 950_000_001, 1_300_000_000),
    ];
    let mut rng = Rng(0xD1B5_4A32_D192_ED03);
    let cases: usize = std::env::var("LOSAT_FUZZ_CASES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(20_000);
    // Groups without / with a pass value outside Int4.
    let mut cases_of = [0usize; 2];
    for case in 0..cases {
        let (c0, c1, gap_prob, lo, hi) = settings[case % settings.len()];
        let cutoffs = LinkHspCutoffs {
            cutoff_small_gap: c0,
            cutoff_big_gap: c1,
            gap_prob,
        };
        let ignore_small_gaps = cutoffs.ignore_small_gaps();
        let n = 1 + rng.below(if case % 40 == 0 { 120 } else { 14 }) as usize;
        let span = [60u64, 200, 1200][case % 3];
        let mut hits: Vec<UngappedHit> = Vec::with_capacity(n);
        for index in 0..n {
            let len = if rng.below(9) == 0 {
                4 + rng.below(3) as usize
            } else {
                7 + rng.below(40) as usize
            };
            let q = rng.below(span) as usize;
            let s = if rng.below(3) > 0 {
                q + rng.below(20) as usize
            } else {
                rng.below(span) as usize
            };
            let score = (lo + rng.below((hi - lo + 1) as u64) as i64) as i32;
            let mut hit = mock_hit(q, q + len, s, s + len, score);
            hit.ctx_idx = u32::try_from(rng.below(3)).unwrap();
            hit.q_frame = contexts[hit.ctx_idx as usize].frame;
            hit.s_frame = 1 + rng.below(3) as i8;
            hit.q_orig_len = 300_000;
            hit.s_orig_len = 300_000;
            hit.hsp_list_order = u32::try_from(index).unwrap();
            hits.push(hit);
        }
        sort_hsps_by_ncbi_link_order(&mut hits, rev_compare_hsps_tbx);
        for (link_id, hit) in hits.iter_mut().enumerate() {
            hit.link_id = u32::try_from(link_id).unwrap();
            hit.chain_next_link_id = None;
            hit.num = 1;
        }
        let label = format!("case {case} n {n} cutoffs {c0} {c1} ignore {ignore_small_gaps}");

        let expected = link_ncbi(hits.clone(), &cutoffs);
        let (incr, stats_incr) = link_incr(hits.clone(), &cutoffs, false);
        assert_same(&expected, &incr, &format!("{label}: incremental kernel"));
        let (wide, _) = link_incr(hits, &cutoffs, true);
        assert_same(
            &expected,
            &wide,
            &format!("{label}: incremental kernel, u128 key"),
        );
        cases_of[usize::from(stats_incr.max_pass_sum > i64::from(i32::MAX))] += 1;
    }
    eprintln!("kernels_add_sums_beyond_int4: {cases} cases; (inside, beyond) Int4 {cases_of:?}");
    assert!(cases_of.iter().all(|&c| c > 0), "{cases_of:?}");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:753,907-908
// The keys of the incremental kernel (linking_incr.rs, "Sums"): a group whose pass values
// exceed the u64 key is linked with the u128 key, and an add-back beyond the u64 key is
// clamped; both give the result of the port of s_BlastEvenGapLinkHSPs.
/// One chain of 300 HSPs (a) whose pass values exceed 2^39 and (b) whose
/// pass values fit the u64 key but whose add-back does not.
#[test]
fn incremental_kernel_links_groups_beyond_its_u64_key() {
    let chain = |score: i32| -> Vec<UngappedHit> {
        fast_test_hits(
            &(0..300)
                .map(|i| (10 * i, 10 * i + 8, 10 * i, 10 * i + 8, score))
                .collect::<Vec<_>>(),
        )
    };
    let plain = IncrOptions {
        verify: true,
        wide_keys: false,
    };
    // (a) pass values up to 300 * (2e9 - 41) > 2^39: the u128 key.
    let cutoffs = fast_test_cutoffs();
    let (_, stats) = link_with_both_kernels(chain(2_000_000_000), &cutoffs, plain, "(a)");
    assert_eq!(stats.wide_groups, 1, "{stats:?}");
    assert!(stats.max_pass_sum > 1 << 39, "{stats:?}");
    // (b) pass values up to 300 * 1e8 with cutoffs near 2e9: the u64 key; the add-back of
    // 300 * 2e9 goes beyond it and is clamped.
    let cutoffs = LinkHspCutoffs {
        cutoff_small_gap: 2_000_000_000,
        cutoff_big_gap: 2_050_000_000,
        gap_prob: 0.5,
    };
    let (_, stats) = link_with_both_kernels(chain(2_100_000_000), &cutoffs, plain, "(b)");
    assert_eq!(stats.wide_groups, 0, "{stats:?}");
    assert!(stats.clamped > 0, "{stats:?}");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:756,871,918,931
// `H_hsp_num` is an Int4 stored into the Int2 count: the count of a set of 32,770 HSPs wraps
// to -32766, which is not "> 1", and the sum statistics give a count below 1 the E-value 0
// (s_BlastSumP, blast_stat.c:4390-4391: it returns 0 when r is below 1). Both kernels keep
// NCBI's count (the literal port overflowed an i16 here in builds with overflow checks).
#[test]
fn a_set_of_more_than_32767_hsps_takes_the_wrapped_int2_count() {
    let chain = |n: usize| -> Vec<UngappedHit> {
        fast_test_hits(
            &(0..n)
                .map(|i| (3 * i, 3 * i + 2, 3 * i, 3 * i + 2, 60))
                .collect::<Vec<_>>(),
        )
    };
    let cutoffs = fast_test_cutoffs();
    let options = IncrOptions {
        verify: true,
        wide_keys: false,
    };
    // The count of the set's first HSP (`hsp_link.num`), then every HSP's after the output
    // list is hooked up (link_hsps.c:1044-1061 copies the first one's count to the set).
    let counts = |n: usize| {
        let (mut linked, _) = link_with_both_kernels(chain(n), &cutoffs, options, "chain");
        assert_eq!(linked.len(), n);
        let heads: Vec<i16> = linked
            .iter()
            .filter(|h| h.start_of_chain)
            .map(|h| h.hsp_link_num)
            .collect();
        assert!(linked.iter().all(|h| h.linked_set), "{:?}", &linked[..2]);
        let evalues: Vec<f64> = linked.iter().map(|h| h.e_value).collect();
        let nums: Vec<i32> = replay_ncbi_output_list(linked.iter_mut().collect())
            .iter()
            .map(|h| h.num)
            .collect();
        (heads, evalues, nums)
    };
    let (heads, evalues, nums) = counts(32_770);
    assert_eq!(heads, [-32766]);
    assert!(evalues.iter().all(|&e| e == 0.0), "{:?}", &evalues[..2]);
    assert!(nums.iter().all(|&n| n == -32766), "{:?}", &nums[..2]);
    // One HSP fewer than the wrap keeps its count.
    let (heads, _, nums) = counts(32_767);
    assert_eq!(heads, [32_767]);
    assert!(nums.iter().all(|&n| n == 32_767), "{:?}", &nums[..2]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:495,939-952
// With large gaps only, the E-value of a set is not divided by `1 - gap_prob`: the same
// HSPs give the same E-values whatever gap_prob is.
#[test]
fn large_gaps_only_do_not_divide_by_the_gap_probability() {
    let hits = || {
        fast_test_hits(&[
            (100, 130, 200, 230, 120),
            (300, 330, 400, 430, 110),
            (500, 530, 600, 630, 100),
        ])
    };
    let large_only = |gap_prob: f64| LinkHspCutoffs {
        cutoff_small_gap: 0,
        cutoff_big_gap: 30,
        gap_prob,
    };
    for options in VERIFIED {
        let (with_zero, _) =
            link_with_both_kernels(hits(), &large_only(0.0), options, "gap_prob 0");
        let (with_half, _) =
            link_with_both_kernels(hits(), &large_only(0.5), options, "gap_prob 0.5");
        // One set of the three HSPs, with an E-value below the INT4_MAX ceiling.
        let heads: Vec<i16> = with_zero
            .iter()
            .filter(|h| h.start_of_chain)
            .map(|h| h.hsp_link_num)
            .collect();
        assert_eq!(heads, [3], "{with_zero:?}");
        assert!(with_zero.iter().all(|h| h.e_value < 1.0), "{with_zero:?}");
        for (a, b) in with_zero.iter().zip(&with_half) {
            assert_eq!(a.e_value.to_bits(), b.e_value.to_bits(), "{a:?} {b:?}");
        }
    }
}
