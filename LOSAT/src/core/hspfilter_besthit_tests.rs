use super::*;

fn key(begin: i32, len: i32, score: i32, evalue: f64) -> BestHitKey {
    BestHitKey {
        begin,
        len,
        score,
        evalue,
    }
}

fn params(overhang: f64, score_edge: f64) -> BestHitParams {
    BestHitParams {
        overhang,
        score_edge,
    }
}

/// The ids of the kept HSPs of each query, per subject.
fn kept(lists: BestHitLists<u32>) -> Vec<(usize, Vec<(u32, Vec<u32>)>)> {
    lists.finalize(|a, b| a.cmp(b))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.hpp:1354-1393:
// the option not given takes its default.
#[test]
fn either_option_turns_the_filter_on() {
    assert_eq!(BestHitParams::from_options(None, None), None);
    assert_eq!(
        BestHitParams::from_options(Some(0.25), None),
        Some(params(0.25, 0.1))
    );
    assert_eq!(
        BestHitParams::from_options(None, Some(0.3)),
        Some(params(0.1, 0.3))
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:271-288: a
// contained HSP with a worse e-value and a lower score per letter (by more than the score
// edge) is dropped; a contained HSP with a higher score per letter is kept.
#[test]
fn a_contained_worse_hsp_is_dropped() {
    let mut lists = BestHitLists::writer(params(0.1, 0.1), 10);
    // 100 letters scoring 200 (2 per letter).
    lists.save(0, 0, key(100, 100, 200, 1e-50), 1);
    // Inside it, 50 letters scoring 80: 2 is above 1.6 / 0.9, dropped.
    lists.save(0, 1, key(120, 50, 80, 1e-20), 2);
    // Inside it, 50 letters scoring 95: 2 is not above 1.9 / 0.9, kept (and it does not
    // remove the first, which does not lie inside it).
    lists.save(0, 2, key(120, 50, 95, 1e-20), 3);
    assert_eq!(kept(lists), vec![(0, vec![(0, vec![1]), (2, vec![3])])]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:290-323: a new
// better HSP removes the HSPs it contains, within its overhang (the pipe's) but not past it.
#[test]
fn a_better_hsp_removes_the_ones_it_contains() {
    let mut lists = BestHitLists::pipe(params(0.1, 0.1), 10);
    lists.save(0, 0, key(105, 40, 40, 1e-5), 1);
    lists.save(0, 1, key(95, 40, 40, 1e-5), 2);
    // 100..200 with an overhang of 10: holds 95..135 but not 85..125.
    lists.save(0, 2, key(85, 40, 40, 1e-5), 3);
    lists.save(0, 3, key(100, 100, 300, 1e-60), 4);
    assert_eq!(kept(lists), vec![(0, vec![(2, vec![3]), (3, vec![4])])]);
    // The writer has no overhang: 95..135 stays.
    let mut lists = BestHitLists::writer(params(0.1, 0.1), 10);
    lists.save(0, 1, key(95, 40, 40, 1e-5), 2);
    lists.save(0, 0, key(105, 40, 40, 1e-5), 1);
    lists.save(0, 3, key(100, 100, 300, 1e-60), 4);
    assert_eq!(kept(lists), vec![(0, vec![(1, vec![2]), (3, vec![4])])]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:291-296: an
// overhang whose `Int4` overflows is `INT_MIN` (the oracle's conversion), and the new HSP
// removes nothing; a smaller overhang removes the HSP it holds.
#[test]
fn an_overhang_past_int4_removes_nothing() {
    let run = |overhang: f64| {
        let mut lists = BestHitLists::pipe(params(overhang, 0.1), 10);
        lists.save(0, 0, key(1000, 100, 100, 1e-5), 1);
        lists.save(0, 1, key(900, 300, 900, 1e-50), 2);
        kept(lists)
    };
    // 2 * 300 * 0.4999999 / 2e-7 is about 1.5e9: an `Int4`.
    assert_eq!(run(0.4999999), vec![(0, vec![(1, vec![2])])]);
    // 2 * 300 * 0.49999999999 / 2e-11 is about 1.5e13: past it.
    assert_eq!(
        run(0.49999999999),
        vec![(0, vec![(1, vec![2]), (0, vec![1])])]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:194-231: the
// subjects come in the order of their first HSP in the list (by start), each list sorted.
#[test]
fn the_lists_follow_the_starts_and_queries_apart() {
    let mut lists = BestHitLists::writer(params(0.1, 0.1), 10);
    lists.save(1, 7, key(500, 10, 30, 1e-3), 9);
    lists.save(0, 4, key(300, 10, 30, 1e-3), 6);
    lists.save(0, 5, key(100, 10, 30, 1e-3), 5);
    lists.save(0, 4, key(50, 10, 30, 1e-3), 8);
    assert_eq!(
        kept(lists),
        vec![
            (0, vec![(4, vec![6, 8]), (5, vec![5])]),
            (1, vec![(7, vec![9])])
        ]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:88-186,340-345:
// past twice `prelim_hitlist_size` HSPs the list is rebuilt without the overhangs, so an
// HSP that the stored overhang of another held is no longer dropped.
#[test]
fn a_rebuild_drops_the_overhangs() {
    let run = |prelim_hitlist_size: usize| {
        // The pipe's overhang of 0.4: 100..200 is stored as 60..240.
        let mut lists = BestHitLists::pipe(params(0.4, 0.1), prelim_hitlist_size);
        lists.save(0, 0, key(100, 100, 300, 1e-10), 1);
        lists.save(0, 1, key(2000, 10, 10, 1e-3), 2);
        lists.save(0, 2, key(3000, 10, 10, 1e-3), 3);
        // 70..220 at 1 per letter: within 60..240, which scores 3 per letter with a better
        // e-value, so it is dropped; within 100..200 (after a rebuild) it is not, and it
        // does not remove 100..200 (whose e-value is better).
        lists.save(0, 3, key(70, 150, 150, 1e-5), 4);
        kept(lists)
    };
    assert_eq!(
        run(10),
        vec![(0, vec![(0, vec![1]), (1, vec![2]), (2, vec![3])])]
    );
    // A `prelim_hitlist_size` of 1 rebuilds at the third HSP.
    assert_eq!(
        run(1),
        vec![(
            0,
            vec![(3, vec![4]), (0, vec![1]), (1, vec![2]), (2, vec![3])]
        )]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:141-161: 2 plus the
// plus strand's three frames, or the minus strand's when the first plus frame is empty.
#[test]
fn a_translated_query_length_is_its_frames_plus_two() {
    // 100 letters: frames of 33, 33 and 32 residues on each strand.
    assert_eq!(translated_query_length(&[33, 33, 32, 33, 33, 32]), 100);
    assert_eq!(translated_query_length(&[0, 0, 0, 33, 33, 32]), 100);
    assert_eq!(translated_query_length(&[33, 33, 32, 0, 0, 0]), 100);
    // 2 letters: no frame has a residue.
    assert_eq!(translated_query_length(&[0; 6]), 2);
}
