use super::super::tests::{fast_test_contexts, fast_test_cutoffs, mock_hit};
use super::super::{
    apply_sum_stats_even_gap_linking_with_parallel, fwd_compare_hsps_transl,
    link_hsp_group_incremental, replay_ncbi_output_list, rev_compare_hsps_tbx,
    rev_compare_hsps_transl, sort_hsps_by_ncbi_link_order, sort_linked_hsp_references,
    translated_context_group,
};
use super::*;
use crate::utils::glibc_fma_libm::GlibcLibm;

/// xorshift64, for reproducible random lists (as in `x_hsp_inplace.rs`).
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

const FRAMES: [i8; 6] = [1, 2, 3, -1, -2, -3];

/// A random HSP of one of `queries` queries (the generator of `x_hsp_inplace.rs`'s tests):
/// few distinct coordinates and scores when `span` is small, every record distinct in
/// `hsp_list_order` and in the bits of `e_value`.
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

/// The list as the group step of a linking call receives it: link ids, no chain link and
/// `num = 1` set as the call does, then sorted by `s_RevCompareHSPsTbx`.
fn sorted(mut hits: Vec<UngappedHit>) -> Vec<UngappedHit> {
    for (link_id, hit) in hits.iter_mut().enumerate() {
        hit.link_id = u32::try_from(link_id).unwrap();
        hit.chain_next_link_id = None;
        hit.num = 1;
    }
    sort_hsps_by_ncbi_link_order(&mut hits, rev_compare_hsps_tbx);
    hits
}

/// The grouping loop of `apply_sum_stats_even_gap_linking_with_parallel`, statement for
/// statement.
fn push_grouping(hits: Vec<UngappedHit>) -> Vec<Vec<UngappedHit>> {
    let mut frame_groups: Vec<Vec<UngappedHit>> = Vec::new();
    let mut current_group: Vec<UngappedHit> = Vec::new();

    for hit in hits {
        let ctx_group = translated_context_group(&hit);
        let s_sign = hit.s_frame.signum();

        if let Some(prev) = current_group.last() {
            let prev_ctx_group = translated_context_group(prev);
            let prev_s_sign = prev.s_frame.signum();

            if prev_ctx_group != ctx_group || prev_s_sign != s_sign {
                // Frame switch - start new group
                if !current_group.is_empty() {
                    frame_groups.push(std::mem::take(&mut current_group));
                }
            }
        }
        current_group.push(hit);
    }
    if !current_group.is_empty() {
        frame_groups.push(current_group);
    }
    frame_groups
}

/// The group step of the reference body (`linking.rs`): the push grouping, each group
/// through the per-group function the body chooses with one scratch sized as there, and the
/// linked groups concatenated in group order.
fn reference_groups(sorted: &[UngappedHit], call: &Call<'_>, use_legacy: bool) -> Vec<UngappedHit> {
    let total_hits = sorted.len();
    let mut pools = BufferPools {
        lh_helpers: Vec::with_capacity((total_hits + 5).max(1024)),
        hsp_links: Vec::with_capacity(total_hits),
    };
    let mut ordered = Vec::new();
    for group_hits in push_grouping(sorted.to_vec()) {
        let processed = if use_legacy {
            link_hsp_group_ncbi(
                group_hits,
                call.params,
                call.cutoffs,
                call.linking_params.gap_decay_rate,
                false,
                call.linking_params.subject_len_nucl,
                call.query_contexts,
                call.subject_frame_bases,
                call.length_adj_per_context,
                call.eff_searchsp_per_context,
                call.log_k_by_ctx,
                &mut pools.lh_helpers,
                &mut pools.hsp_links,
            )
        } else {
            link_hsp_group_incremental(
                group_hits,
                call.params,
                call.cutoffs,
                call.linking_params.gap_decay_rate,
                call.linking_params.subject_len_nucl,
                call.query_contexts,
                call.subject_frame_bases,
                call.length_adj_per_context,
                call.eff_searchsp_per_context,
                call.log_k_by_ctx,
                &mut pools.lh_helpers,
                &mut pools.hsp_links,
            )
        };
        ordered.extend(processed);
    }
    ordered
}

/// The end of the linking call: its two final sorts and its replay (the reference's).
fn reference_replay(mut results: Vec<UngappedHit>) -> Vec<UngappedHit> {
    if results.len() <= 1 {
        return results;
    }
    let mut references: Vec<&mut UngappedHit> = results.iter_mut().collect();
    sort_linked_hsp_references(&mut references, rev_compare_hsps_transl);
    sort_linked_hsp_references(&mut references, fwd_compare_hsps_transl);
    replay_ncbi_output_list(references)
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

// No NCBI counterpart: `body` once on the test thread and once inside pools of 1, 2 and 8
// threads (a build without the parallel feature has no pool and runs it on the test thread
// only), as in `x_link_sort.rs`'s tests. The flag tells `body` whether it runs in a pool.
fn on_and_off_pools(body: impl Fn(&str, bool) + Sync) {
    body("no pool", false);
    #[cfg(feature = "parallel")]
    for threads in [1usize, 2, 8] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("test pool");
        pool.install(|| body(&format!("pool of {threads}"), true));
    }
}

/// Does this build link groups on the pool?
const PARALLEL_BUILD: bool = cfg!(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
));

/// The cutoff sets of the kernel fuzz test (as in `x_hsp_inplace.rs`'s tests).
fn cutoff_sets() -> [LinkHspCutoffs; 3] {
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:522-531
// Module doc, item 1: the slices are the reference's groups (the push loop), the same records
// in the same order, and each one is a view of the list at its run's position: no HSP, one
// HSP, one group, sorted lists of one and three queries, and groups of one HSP each.
#[test]
fn the_slices_are_the_groups_of_the_reference() {
    fn check(mut hits: Vec<UngappedHit>, label: &str) -> usize {
        let expected = push_grouping(hits.clone());
        let starts = group_starts(&hits);
        let base = hits.as_ptr() as usize;
        let size = std::mem::size_of::<UngappedHit>();
        let groups = group_slices(&mut hits, &starts);
        assert_eq!(groups.len(), expected.len(), "{label}: number of groups");
        let mut position = 0usize;
        for (index, (slice, group)) in groups.iter().zip(&expected).enumerate() {
            assert_eq!(
                slice.as_ptr() as usize,
                base + position * size,
                "{label}: group {index} lies in the list at its run's position"
            );
            position += slice.len();
            assert_same_list(&format!("{label}: group {index}"), slice, group);
        }
        groups.len()
    }

    let mut rng = Rng(0xD1B5_4A32_D192_ED03);
    assert_eq!(check(Vec::new(), "no HSP"), 0);
    assert_eq!(check(random_hits(&mut rng, 1, 50, 1), "one HSP"), 1);
    let mut one_group = random_hits(&mut rng, 300, 50, 1);
    for hit in one_group.iter_mut() {
        hit.q_frame = 2;
        hit.s_frame = -1;
    }
    assert_eq!(check(one_group, "one group"), 1);
    for queries in [1usize, 3] {
        for n in [2usize, 5, 64, 3000] {
            let hits = sorted(random_hits(&mut rng, n, 60, queries));
            let groups = check(hits, &format!("sorted, {queries} queries, n {n}"));
            assert!(groups <= 4 * queries);
            if n == 3000 {
                assert_eq!(groups, 4 * queries);
            }
        }
    }
    let mut singles = random_hits(&mut rng, 97, 50, 1);
    for (index, hit) in singles.iter_mut().enumerate() {
        hit.q_idx = (index / 2) as u32;
        hit.q_frame = 1;
        hit.s_frame = if index % 2 == 0 { 1 } else { -1 };
    }
    assert_eq!(check(singles, "groups of one HSP"), 97);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558,1080-1082
// Module doc, items 1-5: on random lists of one and two queries (up to eight groups, chains,
// comparator-equal HSPs) the slices give the reference's group step, every field, with the
// incremental kernel and with the literal port, serially and on pools of 1, 2 and 8 threads,
// in the list's own allocation; after the final sorts and the replay the result is the whole
// linking call's.
#[test]
fn the_slices_give_the_reference_list_serially_and_on_pools() {
    let (params, contexts) = two_query_contexts();
    let log_k: Vec<f64> = contexts
        .iter()
        .map(|c| c.karlin_params.k.glibc_ln())
        .collect();
    let cutoffs = cutoff_sets();
    let cases: usize = std::env::var("LOSAT_FUZZ_CASES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(120);
    let mut rng = Rng(0xA076_1D64_78BD_642F);
    let (mut most_groups, mut chained, mut total) = (0usize, 0usize, 0usize);
    for case in 0..cases {
        let n = match case % 10 {
            0 => case / 10 % 4,
            1 => 2500,
            _ => 1 + rng.below(400) as usize,
        };
        let queries = 1 + case % 2;
        let span = [60u64, 400, 3000, 20_000][case / 2 % 4];
        let hits = random_hits(&mut rng, n, span, queries);
        let linking = LinkingParams {
            subject_len_nucl: 300_000,
            gap_decay_rate: 0.5,
            cutoffs: cutoffs[case % cutoffs.len()],
        };
        let list = sorted(hits.clone());
        let expected_groups = push_grouping(list.clone()).len();
        most_groups = most_groups.max(expected_groups);
        for use_legacy in [false, true] {
            let call = |allow_parallel| Call {
                params: &params,
                cutoffs: &linking.cutoffs,
                linking_params: &linking,
                query_contexts: &contexts,
                subject_frame_bases: &[0; 6],
                length_adj_per_context: &[49; 12],
                eff_searchsp_per_context: &[9_372_428_362; 12],
                log_k_by_ctx: &log_k,
                allow_parallel,
            };
            let expected = reference_groups(&list, &call(false), use_legacy);
            if !use_legacy {
                total += expected.len();
                chained += expected.iter().filter(|hit| hit.linked_set).count();
            }
            on_and_off_pools(|where_, _| {
                for allow_parallel in [false, true] {
                    let label = format!(
                        "case {case} n {n} queries {queries} span {span} legacy {use_legacy} {where_} parallel {allow_parallel}"
                    );
                    let input = list.clone();
                    let allocation = input.as_ptr();
                    let before = PARALLEL_LINKINGS.load(AtomicOrdering::Relaxed);
                    let (actual, groups) = link_groups_as(
                        input,
                        &call(allow_parallel),
                        &mut LinkLap::off(),
                        use_legacy,
                        use_legacy,
                    );
                    assert_eq!(actual.as_ptr(), allocation, "{label}: no second list");
                    assert_eq!(groups, expected_groups, "{label}");
                    assert_same_list(&label, &actual, &expected);
                    if PARALLEL_BUILD && allow_parallel && groups > 1 {
                        assert!(
                            PARALLEL_LINKINGS.load(AtomicOrdering::Relaxed) > before,
                            "{label}: the groups ran on the pool"
                        );
                    }
                }
            });
        }
        // The whole linking call: the reference body against the slices
        // followed by the reference's final sorts and replay.
        for allow_parallel in [false, true] {
            let call = Call {
                params: &params,
                cutoffs: &linking.cutoffs,
                linking_params: &linking,
                query_contexts: &contexts,
                subject_frame_bases: &[0; 6],
                length_adj_per_context: &[49; 12],
                eff_searchsp_per_context: &[9_372_428_362; 12],
                log_k_by_ctx: &log_k,
                allow_parallel,
            };
            let whole_call = apply_sum_stats_even_gap_linking_with_parallel(
                hits.clone(),
                &params,
                &linking,
                &contexts,
                &[0; 6],
                &[49; 12],
                &[9_372_428_362; 12],
                allow_parallel,
            );
            let (linked, _) =
                link_groups_as(list.clone(), &call, &mut LinkLap::off(), false, false);
            assert_same_list(
                &format!("case {case} whole call, parallel {allow_parallel}"),
                &reference_replay(linked),
                &whole_call,
            );
        }
    }
    if cases >= 120 {
        // The lists must exercise chains and every group of two queries.
        assert!(chained * 10 > total, "{chained} of {total}");
        assert_eq!(most_groups, 8);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
// The edges: an empty list (no group, no kernel), one HSP, one group of many HSPs with the
// pool allowed, and a list that arrives with spare capacity (given back, as the reference's cut
// does).
#[test]
fn empty_single_and_one_group_lists() {
    let (params, contexts) = two_query_contexts();
    let log_k: Vec<f64> = contexts
        .iter()
        .map(|c| c.karlin_params.k.glibc_ln())
        .collect();
    let linking = LinkingParams {
        subject_len_nucl: 300_000,
        gap_decay_rate: 0.5,
        cutoffs: fast_test_cutoffs(),
    };
    let call = Call {
        params: &params,
        cutoffs: &linking.cutoffs,
        linking_params: &linking,
        query_contexts: &contexts,
        subject_frame_bases: &[0; 6],
        length_adj_per_context: &[49; 12],
        eff_searchsp_per_context: &[9_372_428_362; 12],
        log_k_by_ctx: &log_k,
        allow_parallel: true,
    };
    let (linked, groups) = link_groups_as(Vec::new(), &call, &mut LinkLap::off(), false, false);
    assert!(linked.is_empty());
    assert_eq!(groups, 0);

    let mut rng = Rng(0x5851_F42D_4C95_7F2D);
    let one = sorted(random_hits(&mut rng, 1, 50, 1));
    let (linked, groups) = link_groups_as(one.clone(), &call, &mut LinkLap::off(), false, false);
    assert_eq!(groups, 1);
    assert_same_list("one HSP", &linked, &reference_groups(&one, &call, false));

    let mut one_group = random_hits(&mut rng, 500, 300, 1);
    for hit in one_group.iter_mut() {
        hit.q_frame = -3;
        hit.ctx_idx = 5;
        hit.s_frame = 2;
    }
    let one_group = sorted(one_group);
    let (linked, groups) =
        link_groups_as(one_group.clone(), &call, &mut LinkLap::off(), false, false);
    assert_eq!(groups, 1);
    assert_same_list(
        "one group",
        &linked,
        &reference_groups(&one_group, &call, false),
    );
    assert!(
        linked.iter().any(|hit| hit.linked_set),
        "the group has chains"
    );

    let mut spare = sorted(random_hits(&mut rng, 300, 400, 2));
    let expected = reference_groups(&spare, &call, false);
    spare.reserve(1000);
    assert!(spare.capacity() > spare.len());
    let (linked, _) = link_groups_as(spare, &call, &mut LinkLap::off(), false, false);
    assert_same_list("spare capacity", &linked, &expected);
}
