use super::*;

fn align(score: i32, q: (i32, i32), s: (i32, i32)) -> Align {
    Align {
        score,
        rule: EMatrixAdjustRule::CompoScaleOldMatrix,
        query_start: q.0,
        query_end: q.1,
        match_start: s.0,
        match_end: s.1,
        frame: 0,
        edit_script: Vec::new(),
    }
}

// `s_WithDistinctEnds`: a new alignment that shares an end with a better one is dropped; one
// that beats it replaces every alignment that shares an end and goes to the head.
#[test]
fn distinct_ends_keep_the_best_of_a_shared_end() {
    let mut list = vec![align(100, (0, 50), (0, 50)), align(80, (60, 90), (60, 90))];
    with_distinct_ends::<XDropRedo>(align(90, (0, 40), (0, 40)), &mut list, true);
    assert_eq!(list.len(), 2);
    with_distinct_ends::<XDropRedo>(align(120, (10, 90), (10, 90)), &mut list, true);
    assert_eq!(
        list.iter().map(|a| a.score).collect::<Vec<_>>(),
        vec![120, 100]
    );
    // With different adjustments a similar end (contained, same diagonal) counts.
    with_distinct_ends::<XDropRedo>(align(110, (5, 45), (5, 45)), &mut list, false);
    assert_eq!(list.len(), 2);
}

// `s_HitlistReapContained` keeps the first of two nested HSPs unless the later scores more.
#[test]
fn contained_hsps_are_reaped() {
    let mut redone = vec![
        Redone {
            align: align(100, (0, 50), (0, 50)),
            e_value: 0.0,
            num: 1,
        },
        Redone {
            align: align(90, (10, 40), (10, 40)),
            e_value: 0.0,
            num: 1,
        },
        Redone {
            align: align(80, (60, 70), (60, 70)),
            e_value: 0.0,
            num: 1,
        },
    ];
    reap_contained::<XDropRedo>(&mut redone);
    assert_eq!(
        redone.iter().map(|r| r.align.score).collect::<Vec<_>>(),
        vec![100, 80]
    );
}

// `s_preliminaryTestNearIdentical`: the subject span must reach min(query length, 50) and
// the score per aligned column the cutoff.
#[test]
fn the_preliminary_near_identical_test_is_ncbis() {
    let a = align(1000, (0, 60), (0, 60));
    assert!(preliminary_test_near_identical(100, &a, 16.0));
    assert!(!preliminary_test_near_identical(100, &a, 17.0));
    let short = align(1000, (0, 20), (0, 20));
    assert!(!preliminary_test_near_identical(100, &short, 1.0));
    assert!(preliminary_test_near_identical(20, &short, 1.0));
}
