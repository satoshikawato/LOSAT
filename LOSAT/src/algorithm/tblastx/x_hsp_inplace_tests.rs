use super::super::tests::{fast_test_contexts, fast_test_cutoffs, mock_hit};
use super::*;

/// xorshift64, for reproducible random lists (as in the linking fuzz test).
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
    /// A random permutation of `0..n` (Fisher-Yates).
    fn shuffled(&mut self, n: usize) -> Vec<usize> {
        let mut values: Vec<usize> = (0..n).collect();
        for index in (1..n).rev() {
            let other = self.below(index as u64 + 1) as usize;
            values.swap(index, other);
        }
        values
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:483-486
// No computation: the test calls the reference linking body past the dispatch. That body
// sorts the list once with the index-replay sort, which is the in-place sort here as well.
/// The body of `apply_sum_stats_even_gap_linking_with_parallel`. The linking dispatch takes
/// the flag unless the list is empty.
#[allow(clippy::too_many_arguments)]
fn reference_linking(
    hits: Vec<UngappedHit>,
    params: &KarlinParams,
    linking_params: &LinkingParams,
    query_contexts: &[QueryContext],
    subject_frame_bases: &[i32],
    length_adj_per_context: &[i64],
    eff_searchsp_per_context: &[i64],
    allow_parallel: bool,
) -> Vec<UngappedHit> {
    let n = hits.len();
    LINKING_REFERENCE_CALL.with(|flag| flag.set(true));
    let linked = apply_sum_stats_even_gap_linking_with_parallel(
        hits,
        params,
        linking_params,
        query_contexts,
        subject_frame_bases,
        length_adj_per_context,
        eff_searchsp_per_context,
        allow_parallel,
    );
    let linking_left = LINKING_REFERENCE_CALL.with(|flag| flag.replace(false));
    assert_eq!(
        linking_left,
        n == 0,
        "the dispatch of the reference linking call takes the flag when it is reached"
    );
    linked
}

/// The same HSPs in the same order, every field, E-values bit for bit.
fn assert_same_list(label: &str, actual: &[UngappedHit], expected: &[UngappedHit]) {
    assert_eq!(actual.len(), expected.len(), "{label}: length");
    for (index, (a, b)) in actual.iter().zip(expected).enumerate() {
        assert_eq!(
            a.e_value.to_bits(),
            b.e_value.to_bits(),
            "{label}: e_value of HSP {index}"
        );
        assert_eq!(format!("{a:?}"), format!("{b:?}"), "{label}: HSP {index}");
    }
    // The comparison of `x_assert_same_hits` must agree with the one above.
    x_assert_same_hits(label, actual, expected);
}

const FRAMES: [i8; 6] = [1, 2, 3, -1, -2, -3];

/// A random HSP of one of `queries` queries: few distinct coordinates and scores when
/// `span` is small (many comparator-equal HSPs), every record distinct in `hsp_list_order`,
/// `link_id` and the bits of `e_value`.
fn random_hit(rng: &mut Rng, index: usize, span: u64, queries: usize) -> UngappedHit {
    let len = 4 + rng.below(40) as usize;
    let q0 = rng.below(span) as usize;
    let s0 = if rng.below(3) == 0 {
        rng.below(span) as usize
    } else {
        q0 + rng.below(25) as usize
    };
    let mut hit = mock_hit(q0, q0 + len, s0, s0 + len, 38 + rng.below(40) as i32);
    hit.ctx_idx = u32::try_from(rng.below(6 * queries as u64)).unwrap();
    hit.q_idx = (hit.ctx_idx as usize / 6) as u32;
    hit.q_frame = FRAMES[hit.ctx_idx as usize % 6];
    hit.s_frame = FRAMES[rng.below(6) as usize];
    hit.s_f_idx = u8::try_from(rng.below(6)).unwrap();
    hit.q_orig_len = 300_000;
    hit.s_orig_len = 300_000;
    hit.hsp_list_order = u32::try_from(index).unwrap();
    hit.link_id = u32::try_from(index).unwrap();
    hit.e_value = f64::from_bits(0x3ff0_0000_0000_0000 + index as u64);
    hit
}

fn random_hits(rng: &mut Rng, n: usize, span: u64, queries: usize) -> Vec<UngappedHit> {
    (0..n)
        .map(|index| random_hit(rng, index, span, queries))
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1380-1381
// The arrangement the reference builds from a copy (`hits[dst] = original[src].clone()`)
// against the same arrangement made by swaps: n = 0, 1, 2, identity, reversal, one cycle
// through every position (both directions), pair swaps and random permutations.
#[test]
fn permutation_in_place_equals_copy_and_replay() {
    fn check(source: Vec<usize>, label: &str) {
        let items: Vec<(usize, String)> = (0..source.len())
            .map(|index| (index, format!("record {index}")))
            .collect();
        let expected: Vec<(usize, String)> = source.iter().map(|&src| items[src].clone()).collect();
        let mut actual = items.clone();
        let mut work = source.clone();
        apply_permutation_in_place(&mut actual, &mut work);
        assert_eq!(actual, expected, "{label}: {source:?}");
        assert!(
            work.iter().enumerate().all(|(index, &src)| index == src),
            "{label}: the permutation is left as the identity"
        );
    }

    check(Vec::new(), "n = 0");
    check(vec![0], "n = 1");
    check(vec![0, 1], "n = 2 identity");
    check(vec![1, 0], "n = 2 swap");
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    for n in [3usize, 4, 7, 64, 1000, 4097] {
        check((0..n).collect(), "identity");
        check((0..n).rev().collect(), "reversal");
        check((0..n).map(|i| (i + 1) % n).collect(), "one cycle, forward");
        check(
            (0..n).map(|i| (i + n - 1) % n).collect(),
            "one cycle, backward",
        );
        check(
            (0..n)
                .map(|i| if (i ^ 1) < n { i ^ 1 } else { i })
                .collect(),
            "pair swaps",
        );
        for round in 0..20 {
            check(rng.shuffled(n), &format!("random {round}"));
        }
    }
}

/// The end of the reference linking call: its two final sorts and its replay.
fn reference_replay(mut results: Vec<UngappedHit>) -> Vec<UngappedHit> {
    let mut references: Vec<&mut UngappedHit> = results.iter_mut().collect();
    sort_linked_hsp_references(&mut references, rev_compare_hsps_transl);
    sort_linked_hsp_references(&mut references, fwd_compare_hsps_transl);
    replay_ncbi_output_list(references)
}

/// Random HSPs linked into random chains the way the kernels leave them: unique `link_id`s
/// in an order unrelated to the list order, heads with `start_of_chain` and the chain's
/// count, members with `linked_set`, single HSPs with neither link nor `linked_set`.
fn chained_hits(rng: &mut Rng, n: usize, max_chain: u64) -> Vec<UngappedHit> {
    let mut hits = random_hits(rng, n, 40, 2);
    let ids = rng.shuffled(n);
    for (hit, &id) in hits.iter_mut().zip(&ids) {
        hit.link_id = u32::try_from(id).unwrap();
        hit.num = 1;
        hit.start_of_chain = rng.below(2) == 0;
    }
    let positions = rng.shuffled(n);
    let mut next = 0usize;
    while next < n {
        let len = (1 + rng.below(max_chain) as usize).min(n - next);
        let chain = &positions[next..next + len];
        next += len;
        if len == 1 {
            continue;
        }
        for (k, &position) in chain.iter().enumerate() {
            let next_id = chain.get(k + 1).map(|&p| hits[p].link_id);
            let hit = &mut hits[position];
            hit.linked_set = true;
            hit.start_of_chain = k == 0;
            hit.chain_next_link_id = next_id;
            hit.hsp_link_num = if k == 0 { len as i16 } else { 1 };
            hit.ordering_method = (len % 2) as u8;
        }
    }
    hits
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1018-1020,1044-1061,1080-1082
// Lists whose chains hold every HSP exactly once (what linking produces): the in-place
// replay must give the reference's list, in the list's own allocation.
#[test]
fn replay_in_place_equals_the_reference_replay_on_chains() {
    let mut rng = Rng(0x2545_F491_4F6C_DD1D);
    let mut chained = 0usize;
    for case in 0..400 {
        let n = match case % 8 {
            0 => 2,
            1 => 3,
            2 => 2000,
            _ => 2 + rng.below(300) as usize,
        };
        let max_chain = [1u64, 2, 6, 40][case % 4];
        let hits = chained_hits(&mut rng, n, max_chain);
        let expected = reference_replay(hits.clone());
        let allocation = hits.as_ptr();
        let (actual, in_place) = replay_in_place_counted(hits, &mut LinkLap::off(), false);
        let label = format!("case {case} n {n} max_chain {max_chain}");
        assert!(in_place, "{label}: every HSP is replayed once");
        assert_eq!(actual.as_ptr(), allocation, "{label}: no second list");
        assert_same_list(&label, &actual, &expected);
        chained += actual.iter().filter(|hit| hit.num > 1).count();
    }
    assert!(chained > 10_000, "{chained}");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1018-1020,1046-1056
// Link fields NCBI's linking does not produce, which the reference replay still answers: a
// chain member no head reaches (left out), a member two heads reach (replayed twice), both
// at once (as many entries as HSPs, but not every HSP), a loop, a link to an unknown HSP and
// a repeated `link_id`. The replay must then return the reference's list (its clones).
#[test]
fn replay_falls_back_when_the_replay_is_not_a_permutation() {
    fn base(n: usize) -> Vec<UngappedHit> {
        (0..n)
            .map(|index| {
                let at = 10 + 30 * index;
                let mut hit = mock_hit(at, at + 20, at, at + 20, 50 + index as i32);
                hit.link_id = u32::try_from(index).unwrap();
                hit.hsp_list_order = u32::try_from(index).unwrap();
                hit.num = 1;
                hit.start_of_chain = true;
                hit
            })
            .collect()
    }
    fn head(hit: &mut UngappedHit, next: usize, count: i16) {
        hit.linked_set = true;
        hit.start_of_chain = true;
        hit.chain_next_link_id = Some(u32::try_from(next).unwrap());
        hit.hsp_link_num = count;
    }
    fn member(hit: &mut UngappedHit) {
        hit.linked_set = true;
        hit.start_of_chain = false;
    }
    fn check(hits: Vec<UngappedHit>, expected_len: usize, expect_in_place: bool, label: &str) {
        let expected = reference_replay(hits.clone());
        assert_eq!(expected.len(), expected_len, "{label}: reference length");
        let (actual, in_place) = replay_in_place_counted(hits, &mut LinkLap::off(), false);
        assert_eq!(in_place, expect_in_place, "{label}");
        assert_same_list(label, &actual, &expected);
    }

    // A member no head reaches is left out.
    let mut hits = base(5);
    member(&mut hits[2]);
    check(hits, 4, false, "orphan member");

    // A member two heads reach is replayed twice, with the count of the later head.
    let mut hits = base(5);
    head(&mut hits[0], 3, 2);
    head(&mut hits[4], 3, 7);
    member(&mut hits[3]);
    check(hits, 6, false, "shared member");

    // Both: as many entries as HSPs, one HSP twice and one missing.
    let mut hits = base(5);
    head(&mut hits[0], 3, 2);
    head(&mut hits[4], 3, 7);
    member(&mut hits[3]);
    member(&mut hits[1]);
    check(hits, 5, false, "shared member and orphan member");

    // A loop: the walk stops after more hops than HSPs (the two single HSPs, then head,
    // member, head, member, head).
    let mut hits = base(4);
    head(&mut hits[1], 2, 2);
    member(&mut hits[2]);
    hits[2].chain_next_link_id = Some(1);
    check(hits, 2 + 5, false, "loop");

    // A link to an HSP that is not in the list ends the chain: still every HSP once.
    let mut hits = base(4);
    head(&mut hits[1], 99, 2);
    check(hits, 4, true, "link to an unknown HSP");

    // A repeated link_id: the later HSP of the sorted order answers for the id.
    let mut hits = base(5);
    head(&mut hits[0], 3, 2);
    member(&mut hits[3]);
    member(&mut hits[4]);
    hits[4].link_id = 3;
    check(hits, 4, false, "repeated link_id");
}

/// Two queries of six contexts each, with the Karlin block of the linking tests.
fn two_query_contexts() -> (KarlinParams, Vec<QueryContext>) {
    let (params, contexts, _) = fast_test_contexts();
    let mut all = contexts.clone();
    for (index, context) in contexts.iter().enumerate() {
        let mut second = context.clone();
        second.q_idx = 1;
        second.frame_base = ((6 + index) * 100_001) as i32;
        all.push(second);
    }
    (params, all)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:414-419,483-486,990-994
// The whole linking call: this path against `apply_sum_stats_even_gap_linking_with_parallel`
// on random HSP lists of one and two queries (up to eight groups, chains, comparator-equal
// HSPs), serially and on the pool, with the cutoff sets of the kernel fuzz test.
#[test]
fn linking_in_place_equals_the_reference_call() {
    let (params, contexts) = two_query_contexts();
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
    ];
    let cases: usize = std::env::var("LOSAT_FUZZ_CASES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(400);
    let mut rng = Rng(0xA076_1D64_78BD_642F);
    let mut total_hsps = 0usize;
    let mut total_chained = 0usize;
    let mut most_groups = 0usize;
    for case in 0..cases {
        let n = match case % 10 {
            0 => case / 10 % 4,
            1 => 2500,
            _ => 1 + rng.below(400) as usize,
        };
        let queries = 1 + case % 2;
        let span = [60u64, 400, 3000, 20_000][case / 2 % 4];
        let mut hits = random_hits(&mut rng, n, span, queries);
        if case % 3 == 0 {
            // A list that arrives with spare capacity, as the search's lists do.
            hits.reserve(n + 7);
        }
        let linking = LinkingParams {
            subject_len_nucl: 300_000,
            gap_decay_rate: 0.5,
            cutoffs: cutoff_sets[case % cutoff_sets.len()],
        };
        for allow_parallel in [false, true] {
            let label = format!(
                "case {case} n {n} queries {queries} span {span} parallel {allow_parallel}"
            );
            let expected = reference_linking(
                hits.clone(),
                &params,
                &linking,
                &contexts,
                &[0; 6],
                &[49; 12],
                &[9_372_428_362; 12],
                allow_parallel,
            );
            let actual = link_in_place(
                hits.clone(),
                &params,
                &linking,
                &contexts,
                &[0; 6],
                &[49; 12],
                &[9_372_428_362; 12],
                allow_parallel,
                rev_compare_hsps_tbx,
                x_link_slices::link_modes,
            );
            assert_same_list(&label, &actual, &expected);
            assert_eq!(actual.len(), n, "{label}: every HSP is in the result");
            if allow_parallel {
                total_hsps += actual.len();
                total_chained += actual.iter().filter(|hit| hit.linked_set).count();
            }
        }
        let mut sorted = hits;
        sort_in_place(
            &mut sorted,
            rev_compare_hsps_tbx,
            x_link_sort::Site::First(false),
        );
        most_groups =
            most_groups.max(group_starts(&sorted).len() + usize::from(!sorted.is_empty()));
    }
    if cases >= 400 {
        // The lists must exercise chains and every group of two queries.
        assert!(
            total_chained * 10 > total_hsps,
            "{total_chained} of {total_hsps}"
        );
        assert_eq!(most_groups, 8);
    }
}

/// An HSP of the one-hit word finder of length 0, -1 or -2 that starts at query offset 0 (so
/// it comes after every other HSP of its group in the link sort), with NCBI's `Int4` ends
/// (below 0 when the subject start is small) and a score of 1 to `max_score`.
fn random_inverted_hit(rng: &mut Rng, index: usize, queries: usize, max_score: u64) -> UngappedHit {
    let mut hit = random_hit(rng, index, 60, queries);
    let length = -(rng.below(3) as i32);
    let s_start = rng.below(6) as u32;
    hit.q_aa_start = 0;
    hit.s_aa_start = s_start;
    hit.q_seed_off = 1;
    hit.s_seed_off = s_start + 1;
    hit.set_int4_ends(length, s_start as i32 + length);
    hit.raw_score = 1 + rng.below(max_score) as i32;
    hit
}

fn one_hit_linking(
    hits: Vec<UngappedHit>,
    linking: &LinkingParams,
    parallel: bool,
) -> Vec<UngappedHit> {
    let (params, contexts) = two_query_contexts();
    link_one_hit_with_int4_ends(
        hits,
        &params,
        linking,
        &contexts,
        &[0; 6],
        &[49; 12],
        &[9_372_428_362; 12],
        parallel,
    )
}

fn default_linking(
    hits: Vec<UngappedHit>,
    linking: &LinkingParams,
    parallel: bool,
) -> Vec<UngappedHit> {
    let (params, contexts) = two_query_contexts();
    link_in_place(
        hits,
        &params,
        linking,
        &contexts,
        &[0; 6],
        &[49; 12],
        &[9_372_428_362; 12],
        parallel,
        rev_compare_hsps_tbx,
        x_link_slices::link_modes,
    )
}

fn test_cutoff_sets() -> [LinkHspCutoffs; 3] {
    [
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
    ]
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:483-486
// The `Int4` comparator orders an end below 0 below every other end; the `u32` one above.
#[test]
fn int4_link_comparator_reads_negative_ends_as_int4() {
    let positive = mock_hit(0, 5, 3, 8, 50);
    let mut inverted = mock_hit(0, 0, 1, 1, 7);
    inverted.set_int4_ends(-1, 0);
    assert_eq!(
        rev_compare_hsps_tbx_int4(&positive, &inverted),
        Ordering::Less
    );
    assert_eq!(
        rev_compare_hsps_tbx(&positive, &inverted),
        Ordering::Greater
    );
    let other = mock_hit(2, 9, 4, 11, 50);
    assert_eq!(
        rev_compare_hsps_tbx_int4(&positive, &other),
        rev_compare_hsps_tbx(&positive, &other)
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:414-419
// On lists of positive HSPs the one-hit form of the linking call (`Int4` comparator, literal
// port for every group) is `link_in_place` (incremental kernel), serially and on the pool.
#[test]
fn int4_linking_of_positive_lists_is_link_in_place() {
    let mut rng = Rng(0x5851_F42D_4C95_7F2D);
    let cutoff_sets = test_cutoff_sets();
    for case in 0..120 {
        let n = 1 + rng.below(160) as usize;
        let queries = 1 + case % 2;
        let span = [60u64, 400, 3000][case % 3];
        let hits = random_hits(&mut rng, n, span, queries);
        let linking = LinkingParams {
            subject_len_nucl: 300_000,
            gap_decay_rate: 0.5,
            cutoffs: cutoff_sets[case % cutoff_sets.len()],
        };
        for parallel in [false, true] {
            let label = format!("case {case} n {n} parallel {parallel}");
            let expected = default_linking(hits.clone(), &linking, parallel);
            let actual = link_in_place_int4(
                hits.clone(),
                &two_query_contexts().0,
                &linking,
                &two_query_contexts().1,
                &[0; 6],
                &[49; 12],
                &[9_372_428_362; 12],
                parallel,
            );
            assert_same_list(&label, &actual, &expected);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:687,699 and :452-453
// HSPs of length 0 or less whose score is at most the linking cutoffs stay alone in the link
// (no chain takes them or is taken by them); placed after the other HSPs of their group, they
// change no value of the others. The linking call keeps them (the re-evaluation's fetch of the
// subject sees them, `subject_link::drop_unread_int4_inverted` drops them after it); linked with
// them and dropped, the list is the list linked without them: the same HSPs, E-values and order
// (the input sorted by score, as the search gives it).
#[test]
fn int4_linking_keeps_hsps_of_length_zero_or_less_that_change_nothing() {
    let mut rng = Rng(0x2545_F491_4F6C_DD1D);
    let cutoff_sets = test_cutoff_sets();
    for case in 0..150 {
        let n = rng.below(120) as usize;
        let queries = 1 + case % 2;
        let span = [60u64, 400, 3000][case % 3];
        let mut base = random_hits(&mut rng, n, span, queries);
        crate::algorithm::tblastx::blast_gapalign::sort_ungapped_hits_by_score_ncbi(&mut base);
        let inverted: Vec<UngappedHit> = (0..1 + rng.below(12) as usize)
            .map(|k| random_inverted_hit(&mut rng, n + k, queries, 20))
            .collect();
        let inverted_count = inverted.len();
        let mut with_inverted = base.clone();
        with_inverted.extend(inverted);
        let linking = LinkingParams {
            subject_len_nucl: 300_000,
            gap_decay_rate: 0.5,
            cutoffs: cutoff_sets[case % cutoff_sets.len()],
        };
        for parallel in [false, true] {
            let label = format!("case {case} n {n} parallel {parallel}");
            let mut actual = one_hit_linking(with_inverted.clone(), &linking, parallel);
            let kept = actual.iter().filter(|hit| hit.is_int4_inverted()).count();
            assert_eq!(kept, inverted_count, "{label}");
            actual.retain(|hit| !hit.is_int4_inverted());
            // The HSPs of length 0 or less score below the others, so they come last in the
            // score sort and every other HSP keeps its list position (`link_id`).
            let expected = default_linking(base.clone(), &linking, parallel);
            assert_same_list(&label, &actual, &expected);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:838-893
// HSPs of length 0 or less with any score, anywhere in the list: every HSP leaves the linking
// call once (those of length 0 or less too, for the re-evaluation's fetch of the subject),
// serially and on the pool alike.
#[test]
fn int4_linking_lets_every_hsp_out_once() {
    let mut rng = Rng(0x1405_7B7E_F767_814F);
    let cutoff_sets = test_cutoff_sets();
    let mut linked_with_inverted = 0usize;
    for case in 0..150 {
        let n = rng.below(100) as usize;
        let queries = 1 + case % 2;
        let span = [60u64, 400][case % 2];
        let mut hits = random_hits(&mut rng, n, span, queries);
        let inverted_count = 1 + rng.below(20) as usize;
        for k in 0..inverted_count {
            let mut hit = random_inverted_hit(&mut rng, n + k, queries, 90);
            let q = rng.below(span) as u32;
            let length = hit.q_end_int4();
            hit.q_aa_start = q;
            hit.s_aa_start = q + rng.below(5) as u32;
            let s = hit.s_aa_start as i32;
            hit.set_int4_ends(q as i32 + length, s + length);
            hits.push(hit);
        }
        crate::algorithm::tblastx::blast_gapalign::sort_ungapped_hits_by_score_int4(&mut hits);
        let linking = LinkingParams {
            subject_len_nucl: 300_000,
            gap_decay_rate: 0.5,
            cutoffs: cutoff_sets[case % cutoff_sets.len()],
        };
        let mut serial = one_hit_linking(hits.clone(), &linking, false);
        let pooled = one_hit_linking(hits.clone(), &linking, true);
        assert_same_list(&format!("case {case}: pool"), &pooled, &serial);
        let mut orders: Vec<u32> = serial.iter().map(|hit| hit.hsp_list_order).collect();
        orders.sort_unstable();
        let all = (n + inverted_count) as u32;
        assert_eq!(orders, (0..all).collect::<Vec<_>>(), "case {case}");
        serial.retain(|hit| !hit.is_int4_inverted());
        assert_eq!(serial.len(), n, "case {case}");
        let without: Vec<UngappedHit> = hits
            .iter()
            .filter(|hit| !hit.is_int4_inverted())
            .cloned()
            .collect();
        let alone = default_linking(without, &linking, false);
        if alone
            .iter()
            .zip(serial.iter())
            .any(|(a, b)| a.e_value.to_bits() != b.e_value.to_bits() || a.num != b.num)
        {
            linked_with_inverted += 1;
        }
    }
    // Some lists link an HSP of length 0 or less into a chain of others, which changes their
    // E-values or counts.
    assert!(linked_with_inverted > 0);
}
