use super::*;

#[test]
fn the_plan_hands_out_the_largest_queries_first() {
    let order = plan(&[3, 0, 3, 3, 1, 2], Stage::Prelim).expect("no query above a quarter");
    assert_eq!(order, vec![0, 2, 3, 5, 4]);
}

#[test]
fn the_owner_loop_stays_for_one_query_or_a_dominant_one() {
    assert!(plan(&[0, 7, 0], Stage::Traceback).is_none());
    assert!(plan(&[], Stage::Traceback).is_none());
    // 5 of 16 is above a quarter, 4 of 16 is not.
    assert!(plan(&[5, 4, 4, 3], Stage::Traceback).is_none());
    assert!(plan(&[4, 4, 4, 4], Stage::Traceback).is_some());
}

#[test]
fn group_in_place_keeps_each_querys_order() {
    let mut items = [(2u32, 'a'), (0, 'b'), (2, 'c'), (1, 'd'), (0, 'e')];
    let counts = counts(items.iter().map(|item| item.0), 4);
    assert_eq!(counts, vec![2, 1, 2, 0]);
    group_in_place(&mut items, |item| item.0, &counts);
    assert_eq!(items, [(0, 'b'), (0, 'e'), (1, 'd'), (2, 'a'), (2, 'c')]);
    assert_eq!(starts(&counts), vec![0, 2, 3, 5, 5]);
    // Longer cycles: a reversed list of 3 queries x 50 items, each query's items numbered.
    let mut items: Vec<(u32, u32)> = (0..150u32).rev().map(|i| (i % 3, i)).collect();
    let three = super::counts(items.iter().map(|item| item.0), 3);
    group_in_place(&mut items, |item| item.0, &three);
    let mut want: Vec<(u32, u32)> = (0..150u32).rev().map(|i| (i % 3, i)).collect();
    want.sort_by_key(|item| item.0);
    assert_eq!(items, want);
}

#[test]
fn in_query_order_reads_slices() {
    assert!(in_query_order([0u32, 0, 1, 3, 3].into_iter()));
    assert!(in_query_order(std::iter::empty()));
    assert!(!in_query_order([0u32, 2, 1].into_iter()));
}

#[test]
fn the_rule_names_empty_stages_and_dominant_queries() {
    assert_eq!(choose(&[0, 0, 0]), Choice::Empty);
    assert_eq!(choose(&[]), Choice::Empty);
    assert_eq!(choose(&[0, 7, 0]), Choice::Owner);
    assert_eq!(choose(&[5, 4, 4, 3]), Choice::Owner);
    assert_eq!(choose(&[4, 4, 4, 4]), Choice::Queries);
}

#[test]
fn sort_stable_keeps_ties_in_order() {
    let items: Vec<(u32, u32)> = (0..2_000u32).map(|i| ((i * 7) % 13, i)).collect();
    let mut want = items.clone();
    want.sort_by_key(|item| item.0);
    let mut serial = items.clone();
    sort_stable(&SearchPool::x_serial(), &mut serial, |a, b| a.0.cmp(&b.0));
    assert_eq!(serial, want);
    #[cfg(feature = "parallel")]
    crate::utils::threading::with_search_pool(4, "blastn", |pool| {
        let mut pooled = items.clone();
        sort_stable(pool, &mut pooled, |a, b| a.0.cmp(&b.0));
        assert_eq!(pooled, want);
        Ok(())
    })
    .unwrap();
}

/// On a pool of four threads every query of the order runs exactly once and its result lands
/// at its index, also with more threads than queries.
#[cfg(feature = "parallel")]
#[test]
fn run_on_a_pool_runs_each_query_once() {
    use std::sync::atomic::AtomicUsize;
    crate::utils::threading::with_search_pool(4, "blastn", |pool| {
        for order in [(0..40u32).rev().collect::<Vec<_>>(), vec![3u32, 1]] {
            let calls: Vec<AtomicUsize> = (0..41).map(|_| AtomicUsize::new(0)).collect();
            let results = run(
                pool,
                Stage::Prelim,
                &order,
                41,
                || (),
                |_, query| {
                    calls[query as usize].fetch_add(1, Ordering::SeqCst);
                    query * 3
                },
            );
            for (query, result) in results.iter().enumerate() {
                let ran = order.contains(&(query as u32));
                assert_eq!(*result, ran.then_some(query as u32 * 3), "query {query}");
                assert_eq!(calls[query].load(Ordering::SeqCst), usize::from(ran));
            }
        }
        Ok(())
    })
    .unwrap();
}

#[test]
fn run_gives_every_query_its_own_result() {
    let order = [4u32, 1, 2];
    let serial = run(
        &SearchPool::x_serial(),
        Stage::Traceback,
        &order,
        6,
        || 0usize,
        |calls, query| {
            *calls += 1;
            query * 10
        },
    );
    assert_eq!(serial, vec![None, Some(10), Some(20), None, Some(40), None]);
}
