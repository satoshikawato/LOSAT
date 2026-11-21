use super::*;

fn key(begin: i32, end: i32, score: i32) -> CullKey {
    CullKey {
        begin,
        end,
        score,
        subject_offset: i64::from(begin),
    }
}

/// The payload of a test HSP is its score and its position in the input.
type Payload = (i32, usize);

fn forest(limit: i64) -> CullingForest<Payload> {
    CullingForest::new(limit)
}

fn by_score(a: &Payload, b: &Payload) -> Ordering {
    b.0.cmp(&a.0)
}

fn kept(forest: CullingForest<Payload>) -> Vec<(u32, Vec<i32>)> {
    forest
        .finalize(|_| 0, by_score)
        .into_iter()
        .flat_map(|(_, lists)| lists)
        .map(|(oid, list)| (oid, list.into_iter().map(|h| h.0).collect()))
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:338-340: the preliminary writer adds
// 3 to a limit above 1. NCBI's Int4 sum overflows near INT4_MAX; LOSAT's does not.
#[test]
fn the_preliminary_limit_adds_three_above_one_without_overflow() {
    assert_eq!(CullingForest::<()>::prelim_max(0), 0);
    assert_eq!(CullingForest::<()>::prelim_max(1), 1);
    assert_eq!(CullingForest::<()>::prelim_max(2), 5);
    assert_eq!(CullingForest::<()>::prelim_max(2147483644), 2147483647);
    assert_eq!(CullingForest::<()>::prelim_max(2147483645), 2147483648);
    assert_eq!(CullingForest::<()>::prelim_max(i32::MAX), 2147483650);
}

// NCBI frees a culled HSP at once (`s_HSPFree`): slots, entries and nodes are reused, so
// a stream of ever better overlapping HSPs, each culling its predecessor, and a stream of
// dominated HSPs both keep the arenas at a few items, and the survivors are unchanged.
#[test]
fn freed_slots_entries_and_nodes_are_reused() {
    let mut f = forest(1);
    for i in 0..1000 {
        f.save(0, 300, 0, key(10, 40, 10 + i), (10 + i, i as usize));
    }
    for i in 0..1000 {
        f.save(0, 300, 1, key(10, 40, 5), (5, i as usize));
    }
    assert!(f.hsps.len() <= 3 && f.entries.len() <= 3 && f.nodes.len() <= 3);
    assert_eq!(kept(f), vec![(0, vec![1009])]);
    // Disjoint HSPs fork the root; a long HSP then culls them all and frees the nodes.
    let mut f = forest(1);
    for i in 0..100 {
        f.save(0, 1000, 0, key(i * 5, i * 5 + 4, 10), (10, i as usize));
    }
    let nodes = f.nodes.len();
    f.save(0, 1000, 0, key(0, 600, 5000), (5000, 100));
    for i in 0..100 {
        f.save(0, 1000, 0, key(i * 5, i * 5 + 4, 10), (10, i as usize));
    }
    assert_eq!(f.hsps.len(), 101);
    assert!(f.nodes.len() <= nodes + 3);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:136-163: an HSP never dominates
// itself (`r != y`), so a limit of 1 keeps the best HSP of overlapping ones, whichever
// comes first.
#[test]
fn a_limit_of_one_keeps_the_dominating_hsp() {
    let mut f = forest(1);
    f.save(0, 300, 0, key(10, 40, 50), (50, 0));
    f.save(0, 300, 0, key(12, 40, 30), (30, 1));
    assert_eq!(kept(f), vec![(0, vec![50])]);
    let mut f = forest(1);
    f.save(0, 300, 0, key(12, 40, 30), (30, 0));
    f.save(0, 300, 0, key(10, 40, 50), (50, 1));
    assert_eq!(kept(f), vec![(0, vec![50])]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:79-120: with equal score and range
// the lower OID dominates; another tree has its own HSPs.
#[test]
fn ties_are_broken_by_oid_and_trees_do_not_interact() {
    let mut f = forest(1);
    f.save(0, 300, 1, key(10, 40, 50), (50, 0));
    f.save(0, 300, 0, key(10, 40, 50), (50, 1));
    f.save(1, 300, 0, key(10, 40, 50), (50, 2));
    let all = f.finalize(|cid| cid, by_score);
    assert_eq!(all.len(), 2);
    assert_eq!(all[0].1, vec![(0, vec![(50, 1)])]);
    assert_eq!(all[1].1, vec![(0, vec![(50, 2)])]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:100-112: with the same OID the HSP
// with the lower subject offset dominates, and equal keys dominate each other, so the
// HSP saved later loses.
#[test]
fn same_oid_ties_follow_the_subject_offset_and_the_insertion_order() {
    let mut f = forest(1);
    let mut late = key(10, 40, 50);
    late.subject_offset = 90;
    f.save(0, 300, 0, late, (50, 0));
    f.save(0, 300, 0, key(10, 40, 50), (50, 1));
    let lists = f.finalize(|_| 0, by_score);
    assert_eq!(lists[0].1, vec![(0, vec![(50, 1)])]);
    let mut f = forest(1);
    f.save(0, 300, 0, key(10, 40, 50), (50, 0));
    f.save(0, 300, 0, key(10, 40, 50), (50, 1));
    let lists = f.finalize(|_| 0, by_score);
    assert_eq!(lists[0].1, vec![(0, vec![(50, 0)])]);
}

// The merit does not wrap (NCBI's Int4 does: there 2147483646 and 2147483647 act as a
// limit of 1 in the preliminary stage): the limits at the top of the range keep both
// HSPs, as any limit above the number of dominating HSPs does.
#[test]
fn merits_near_int4_max_do_not_wrap() {
    for limit in [2147483645, 2147483646, i32::MAX] {
        let mut f = forest(CullingForest::<()>::prelim_max(limit));
        f.save(0, 300, 0, key(10, 40, 50), (50, 0));
        f.save(0, 300, 0, key(12, 40, 30), (30, 1));
        assert_eq!(kept(f), vec![(0, vec![50, 30])]);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:258-299,453-458: a leaf with 20
// HSPs forks; HSPs that do not overlap are all kept.
#[test]
fn many_disjoint_hsps_fork_and_survive() {
    let mut f = forest(1);
    for i in 0..40 {
        f.save(0, 3000, 0, key(i * 20, i * 20 + 10, 30), (30, i as usize));
    }
    assert_eq!(kept(f)[0].1.len(), 40);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:396-423,516-581: the rip takes a
// node's list, then the left subtree, then the right subtree; the lists of a query are
// made in the order of the first HSP of each subject.
#[test]
fn final_layout_follows_the_rip_order() {
    let mut f = forest(5);
    // The first HSP spans the middle (150) and stays in the root; twenty more fork
    // the root into a left and a right child.
    f.save(0, 300, 7, key(140, 160, 40), (40, 0));
    for i in 0..10 {
        f.save(
            0,
            300,
            3,
            key(i * 10, i * 10 + 5, 20),
            (20, 100 + i as usize),
        );
    }
    for i in 0..10 {
        f.save(
            0,
            300,
            5,
            key(200 + i * 10, 205 + i * 10, 21),
            (21, 200 + i as usize),
        );
    }
    let lists = f.finalize(|_| 0, by_score);
    let oids: Vec<u32> = lists[0].1.iter().map(|(oid, _)| *oid).collect();
    assert_eq!(oids, vec![7, 3, 5]);
    assert_eq!(lists[0].1[1].1.len(), 10);
    assert_eq!(lists[0].1[2].1.len(), 10);
}

// Queries come in order, each with its own lists and the lists of an OID found again.
#[test]
fn queries_come_in_order_with_their_own_lists() {
    let mut f = forest(10);
    f.save(3, 100, 4, key(0, 10, 9), (9, 0));
    f.save(1, 100, 4, key(0, 10, 8), (8, 1));
    f.save(2, 100, 2, key(50, 60, 7), (7, 2));
    f.save(2, 100, 4, key(70, 80, 6), (6, 3));
    let out = f.finalize(|cid| cid / 2, by_score);
    assert_eq!(out.len(), 2);
    assert_eq!(out[0].0, 0);
    assert_eq!(out[0].1, vec![(4, vec![(8, 1)])]);
    assert_eq!(out[1].0, 1);
    // Tree 2 rips its newest HSP first (subject 4), then tree 3 adds to subject 4's list,
    // which is sorted again.
    assert_eq!(out[1].1, vec![(4, vec![(9, 0), (6, 3)]), (2, vec![(7, 2)])]);
}
