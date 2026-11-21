use super::super::rev_compare_hsps_tbx;
use super::super::tests::mock_hit;
use super::*;
use crate::algorithm::tblastx::blast_gapalign::score_compare_ungapped_hits_ncbi;

type Compare = fn(&UngappedHit, &UngappedHit) -> Ordering;

/// The four comparators the index sorts receive (module doc, item 3).
const COMPARATORS: [(&str, Compare); 4] = [
    ("ScoreCompareHSPs", score_compare_ungapped_hits_ncbi),
    ("s_RevCompareHSPsTbx", rev_compare_hsps_tbx),
    ("s_RevCompareHSPsTransl", rev_compare_hsps_transl),
    ("s_FwdCompareHSPsTransl", fwd_compare_hsps_transl),
];

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

/// `n` random HSPs of two queries. A small `span` gives few distinct offsets and scores, so
/// many HSPs that a comparator calls equal; every record differs in `hsp_list_order`.
fn random_hits(rng: &mut Rng, n: usize, span: u64) -> Vec<UngappedHit> {
    (0..n)
        .map(|index| {
            let q0 = rng.below(span) as usize;
            let s0 = rng.below(span) as usize;
            let q_len = 4 + rng.below(3) as usize;
            let s_len = 4 + rng.below(3) as usize;
            let score = 38 + rng.below(span.min(40)) as i32;
            let mut hit = mock_hit(q0, q0 + q_len, s0, s0 + s_len, score);
            hit.ctx_idx = u32::try_from(rng.below(12)).unwrap();
            hit.q_idx = (hit.ctx_idx as usize / 6) as u32;
            hit.q_frame = FRAMES[hit.ctx_idx as usize % 6];
            hit.s_frame = FRAMES[rng.below(6) as usize];
            hit.hsp_list_order = u32::try_from(index).unwrap();
            hit
        })
        .collect()
}

/// The lists the sorts are tried on: random, already in the comparator's order, and in the
/// reverse of it (whole runs that descend, with equal neighbours inside them).
fn arrangements(rng: &mut Rng, n: usize, span: u64, compare: Compare) -> Vec<Vec<UngappedHit>> {
    let random = random_hits(rng, n, span);
    let mut ascending = random.clone();
    ascending.sort_by(compare);
    let mut descending = ascending.clone();
    descending.reverse();
    vec![random, ascending, descending]
}

// No NCBI counterpart: `body` once on the test thread and once inside pools of 1, 2 and 8
// threads (a serial build has no pool and runs it on the test thread only). The flag tells
// `body` whether it runs on a pool worker.
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

/// Does this build have the parallel sort?
const PARALLEL_BUILD: bool = cfg!(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
));

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:359-366
// Module doc, "Terms" and item 3: each comparator, and the composite of the final sort, is a
// total preorder on random records (reversal, and transitivity of `Less` and of `Equal`).
#[test]
fn the_comparators_are_total_preorders() {
    let composite: Compare =
        |a, b| fwd_compare_hsps_transl(a, b).then_with(|| rev_compare_hsps_transl(a, b));
    let mut all = COMPARATORS.to_vec();
    all.push(("composite of the final sort", composite));
    let mut rng = Rng(0x243F_6A88_85A3_08D3);
    for (name, compare) in all {
        let hits = random_hits(&mut rng, 60, 2);
        let mut equal_pairs = 0usize;
        for a in &hits {
            assert_eq!(compare(a, a), Ordering::Equal, "{name}");
            for b in &hits {
                let ab = compare(a, b);
                assert_eq!(compare(b, a), ab.reverse(), "{name}");
                equal_pairs += usize::from(ab == Ordering::Equal);
                for c in &hits {
                    let bc = compare(b, c);
                    if ab != Ordering::Greater && bc != Ordering::Greater {
                        let expected = if ab == Ordering::Equal && bc == Ordering::Equal {
                            Ordering::Equal
                        } else {
                            Ordering::Less
                        };
                        assert_eq!(compare(a, c), expected, "{name}");
                    }
                }
            }
        }
        // Different records that the comparator calls equal must be among them.
        assert!(equal_pairs > hits.len(), "{name}: {equal_pairs}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:485-486
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1380-1381
// Module doc, item 2: this path's index order equals today's for each comparator, on and off
// pools of 1, 2 and 8 threads, for lengths around rayon's insertion-sort limit (20) and chunk
// length (2,000) with every list sent to the parallel sort, and for lengths below, at and
// above `PARALLEL_MIN` with the real threshold. Random lists, lists already in order, and
// lists in the reverse order (rayon reverses a chunk that descends strictly).
#[test]
fn the_index_order_is_todays_at_any_number_of_threads() {
    let short: [usize; 11] = [0, 1, 2, 3, 20, 21, 1999, 2000, 2001, 3999, 4001];
    let mut rng = Rng(0x0123_4567_89AB_CDEF);
    for (name, compare) in COMPARATORS {
        // (threshold, list, today's index order)
        let mut cases: Vec<(usize, Vec<UngappedHit>, Vec<usize>)> = Vec::new();
        let mut add = |parallel_min: usize, hits: Vec<UngappedHit>| {
            let mut today: Vec<usize> = (0..hits.len()).collect();
            sort_today(&mut today, &hits, compare);
            cases.push((parallel_min, hits, today));
        };
        for n in short {
            for span in [3u64, 100_000] {
                for hits in arrangements(&mut rng, n, span, compare) {
                    add(0, hits);
                }
            }
        }
        add(PARALLEL_MIN, random_hits(&mut rng, PARALLEL_MIN - 1, 60));
        for hits in arrangements(&mut rng, PARALLEL_MIN, 60, compare) {
            add(PARALLEL_MIN, hits);
        }
        add(PARALLEL_MIN, random_hits(&mut rng, 70_001, 60));
        let equal_neighbours: usize = cases
            .iter()
            .map(|(_, hits, today)| {
                today
                    .windows(2)
                    .filter(|pair| compare(&hits[pair[0]], &hits[pair[1]]) == Ordering::Equal)
                    .count()
            })
            .sum();
        // The lists must hold comparator-equal HSPs, whose order only stability decides.
        assert!(equal_neighbours > 10_000, "{name}: {equal_neighbours}");
        on_and_off_pools(|where_, in_pool| {
            for (case, (parallel_min, hits, today)) in cases.iter().enumerate() {
                let n = hits.len();
                // Not allowed to use the pool: today's statement, once is enough.
                for on_pool in [false, true] {
                    if !on_pool && in_pool {
                        continue;
                    }
                    let label = format!("{name} case {case} n {n} {where_} on_pool {on_pool}");
                    let mut new: Vec<usize> = (0..n).collect();
                    let parallel = sort_stable(&mut new, hits, compare, on_pool, *parallel_min);
                    assert_eq!(&new, today, "{label}");
                    assert_eq!(
                        parallel,
                        PARALLEL_BUILD && on_pool && n >= *parallel_min,
                        "{label}"
                    );
                }
            }
        });
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:485-486
// (S) of the module doc does not need the indices to start as `0..n`: the same equality from
// a shuffled index vector, where the stable order is "comparator, then position in the
// vector".
#[test]
fn the_index_order_is_todays_from_any_starting_order() {
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    for (name, compare) in COMPARATORS {
        for n in [50usize, 2500, 12_000] {
            let hits = random_hits(&mut rng, n, 6);
            let mut start: Vec<usize> = (0..n).collect();
            for index in (1..n).rev() {
                start.swap(index, rng.below(index as u64 + 1) as usize);
            }
            let mut today = start.clone();
            sort_today(&mut today, &hits, compare);
            on_and_off_pools(|where_, _| {
                let mut new = start.clone();
                sort_stable(&mut new, &hits, compare, true, 0);
                assert_eq!(new, today, "{name} n {n} {where_}");
            });
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:990-994
// Module doc, item 1: one stable sort with "second, ties by first" equals the two stable
// sorts in sequence. The pair of the linking call, and pairs of comparators whose ties do
// not coincide; the lists hold pairs of HSPs on which the two comparators disagree, on which
// only the second ties, and on which both tie.
#[test]
fn one_composite_sort_equals_the_two_sorts_in_sequence() {
    let pairs: [(&str, Compare, Compare); 5] = [
        (
            "the linking call: s_RevCompareHSPsTransl, then s_FwdCompareHSPsTransl",
            rev_compare_hsps_transl,
            fwd_compare_hsps_transl,
        ),
        (
            "ScoreCompareHSPs, then s_FwdCompareHSPsTransl",
            score_compare_ungapped_hits_ncbi,
            fwd_compare_hsps_transl,
        ),
        (
            "s_RevCompareHSPsTbx, then ScoreCompareHSPs",
            rev_compare_hsps_tbx,
            score_compare_ungapped_hits_ncbi,
        ),
        (
            "s_RevCompareHSPsTbx, then s_FwdCompareHSPsTransl",
            rev_compare_hsps_tbx,
            fwd_compare_hsps_transl,
        ),
        (
            "ScoreCompareHSPs, then s_RevCompareHSPsTbx",
            score_compare_ungapped_hits_ncbi,
            rev_compare_hsps_tbx,
        ),
    ];
    let mut rng = Rng(0xD1B5_4A32_D192_ED03);
    for (pair_index, (name, first, second)) in pairs.into_iter().enumerate() {
        let mut disagree = 0usize;
        let mut second_ties_only = 0usize;
        let mut both_tie = 0usize;
        let mut cases: Vec<(usize, Vec<UngappedHit>)> = Vec::new();
        for (n, span, parallel_min) in [
            (0usize, 2u64, 0usize),
            (1, 2, 0),
            (2, 2, 0),
            (300, 2, 0),
            (300, 40, 0),
            (5000, 2, 0),
            (5000, 100_000, 0),
            (20_000, 1, 0),
            (PARALLEL_MIN - 1, 3, PARALLEL_MIN),
            (PARALLEL_MIN + 3, 3, PARALLEL_MIN),
        ] {
            let hits = random_hits(&mut rng, n, span);
            for pair in hits.windows(2) {
                let (f, s) = (first(&pair[0], &pair[1]), second(&pair[0], &pair[1]));
                disagree += usize::from(f != Ordering::Equal && f == s.reverse());
                second_ties_only += usize::from(s == Ordering::Equal && f != Ordering::Equal);
                both_tie += usize::from(s == Ordering::Equal && f == Ordering::Equal);
            }
            cases.push((parallel_min, hits));
        }
        assert!(disagree > 100, "{name}: {disagree}");
        assert!(both_tie > 10, "{name}: {both_tie}");
        if pair_index == 0 {
            // The pair of the linking call reads the same keys: a tie of one is a tie of both.
            assert_eq!(second_ties_only, 0, "{name}");
        } else {
            assert!(second_ties_only > 10, "{name}: {second_ties_only}");
        }
        on_and_off_pools(|where_, _| {
            for (parallel_min, hits) in &cases {
                let n = hits.len();
                let mut today: Vec<usize> = (0..n).collect();
                sort_today(&mut today, hits, first);
                sort_today(&mut today, hits, second);
                for on_pool in [false, true] {
                    let mut order: Vec<usize> = (0..n).collect();
                    final_sort_as(&mut order, hits, first, second, on_pool, *parallel_min);
                    assert_eq!(order, today, "{name} n {n} {where_} on_pool {on_pool}");
                }
            }
        });
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-149
// Module doc, "The pool and `-num_threads`": a linking call's sorts follow its flag; a
// dispatch sort follows the calling thread, which is not a rayon worker off a pool (the
// thread of a `-num_threads 1` search) and is one inside `ThreadPool::install`, where the
// current pool has the threads it was built with.
// EXPERIMENT (LOSAT_X_TBXBOARD): the board's chunk sorts and merges give today's order for every
// comparator, list arrangement and chunk count (lists below and above the chunk size), on pools
// of 2 and 8 threads.
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
#[test]
fn the_board_sort_is_todays() {
    let mut rng = Rng(0x1319_8A2E_0370_7344);
    for threads in [2usize, 8] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("test pool");
        for (name, compare) in COMPARATORS {
            for n in [0usize, 1, 7, 8191, 8192, 30_001, 70_000] {
                for span in [2u64, 50] {
                    for hits in arrangements(&mut rng, n, span, compare) {
                        let mut expected: Vec<usize> = (0..n).collect();
                        sort_today(&mut expected, &hits, compare);
                        let mut actual: Vec<usize> = (0..n).collect();
                        pool.install(|| {
                            let board = crate::algorithm::tblastx::x_work_board::XBoard::default();
                            x_board_sort(&board, &mut actual, &hits, &compare);
                        });
                        assert_eq!(
                            actual, expected,
                            "{name}, n {n}, span {span}, {threads} threads"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn the_pool_decision_follows_the_flag_or_the_calling_thread() {
    on_and_off_pools(|where_, in_pool| {
        for site in [Site::First(false), Site::Final(false)] {
            assert!(!pool_allowed(site), "{site:?} {where_}");
        }
        for site in [Site::First(true), Site::Final(true)] {
            assert_eq!(pool_allowed(site), PARALLEL_BUILD, "{site:?} {where_}");
        }
        assert_eq!(
            pool_allowed(Site::Dispatch),
            PARALLEL_BUILD && in_pool,
            "{where_}"
        );
    });
    #[cfg(feature = "parallel")]
    {
        assert_eq!(rayon::current_thread_index(), None);
        for threads in [1usize, 2, 8] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("test pool");
            let (index, size) =
                pool.install(|| (rayon::current_thread_index(), rayon::current_num_threads()));
            assert!(index.is_some_and(|index| index < threads));
            assert_eq!(size, threads);
        }
    }
}
