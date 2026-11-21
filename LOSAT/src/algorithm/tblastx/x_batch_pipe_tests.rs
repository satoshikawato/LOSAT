use super::*;
use std::sync::atomic::AtomicUsize;

fn pool(threads: usize) -> rayon::ThreadPool {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("a test pool")
}

/// A fake batch search: `subjects` "subjects" on the pool, each ending its seed stage with
/// `range.start + 1` HSPs, then a value made from the range; it records how many searches run
/// at once.
fn fake_search<'a>(
    running: &'a AtomicUsize,
    most: &'a AtomicUsize,
    subjects: usize,
    fail_at: Option<usize>,
) -> impl Fn(std::ops::Range<usize>, Option<&XSeedEnd<'_>>) -> Result<Vec<usize>> + Sync + 'a {
    move |range, seed_end| {
        let now = running.fetch_add(1, Ordering::SeqCst) + 1;
        most.fetch_max(now, Ordering::SeqCst);
        use rayon::prelude::*;
        (0..subjects).into_par_iter().for_each(|s| {
            let mut guard = XSeedGuard::new(seed_end);
            let mut x = 0u64;
            for k in 0..20_000u64 {
                x = x.wrapping_mul(31).wrapping_add(k ^ s as u64);
            }
            std::hint::black_box(x);
            guard.end(range.start + 1);
            // the "chain" after the seed stage
            for k in 0..20_000u64 {
                x = x.wrapping_mul(17).wrapping_add(k);
            }
            std::hint::black_box(x);
        });
        running.fetch_sub(1, Ordering::SeqCst);
        if fail_at == Some(range.start) {
            anyhow::bail!("batch {} failed", range.start);
        }
        Ok(range.collect())
    }
}

fn run_loop(
    threads: usize,
    batches: usize,
    limit: u64,
    fail_at: Option<usize>,
) -> (Vec<Result<Vec<usize>>>, usize, u64) {
    let plan: Vec<_> = (0..batches).map(|b| b * 2..b * 2 + 2).collect();
    let running = AtomicUsize::new(0);
    let most = AtomicUsize::new(0);
    let search = fake_search(&running, &most, 3, fail_at);
    let pipe = XPipe::with_limit(plan.clone(), 3, &search, limit).expect("two batches or more");
    let results = pool(threads).install(|| {
        rayon::in_place_scope(|scope| {
            let _closer = pipe.closer();
            let mut out = Vec::new();
            for range in plan.iter().cloned() {
                let result = pipe.take(scope, range);
                let failed = result.is_err();
                out.push(result);
                if failed {
                    break;
                }
            }
            out
        })
    });
    (
        results,
        most.load(Ordering::SeqCst),
        pipe.counts.early.load(Ordering::SeqCst),
    )
}

#[test]
fn the_loop_gets_every_batch_in_order_at_any_thread_count() {
    for threads in [1, 2, 4, 8, 16] {
        let (results, most, early) = run_loop(threads, 6, u64::MAX, None);
        let values: Vec<Vec<usize>> = results.into_iter().map(|r| r.expect("no error")).collect();
        let expected: Vec<Vec<usize>> = (0..6).map(|b| vec![b * 2, b * 2 + 1]).collect();
        assert_eq!(values, expected, "threads {threads}");
        assert!(most <= 6, "threads {threads}");
        assert_eq!(
            early, 5,
            "every batch after the first starts at the seed end before it"
        );
    }
}

#[test]
fn a_limit_below_one_batch_keeps_one_search_at_a_time() {
    for threads in [2, 8] {
        let (results, most, early) = run_loop(threads, 4, 0, None);
        assert_eq!(results.len(), 4);
        assert!(results.iter().all(Result::is_ok));
        assert_eq!(early, 0);
        assert_eq!(most, 1, "threads {threads}");
    }
}

#[test]
fn an_error_stops_the_early_starts_and_is_taken_in_order() {
    let (results, _, _) = run_loop(4, 6, u64::MAX, Some(4));
    assert_eq!(results.len(), 3);
    assert!(results[0].is_ok() && results[1].is_ok());
    assert_eq!(
        format!("{:#}", results[2].as_ref().unwrap_err()),
        "batch 4 failed"
    );
}

#[test]
fn a_key_the_plan_does_not_hold_is_searched_on_the_loop() {
    let plan = vec![0..2, 2..4, 4..6];
    let running = AtomicUsize::new(0);
    let most = AtomicUsize::new(0);
    let search = fake_search(&running, &most, 2, None);
    let pipe = XPipe::with_limit(plan, 2, &search, u64::MAX).expect("three batches");
    let values = pool(4).install(|| {
        rayon::in_place_scope(|scope| {
            vec![
                pipe.take(scope, 0..2).expect("planned"),
                pipe.take(scope, 2..3).expect("not planned"),
            ]
        })
    });
    assert_eq!(values, vec![vec![0, 1], vec![2]]);
    assert_eq!(pipe.counts.key_mismatch.load(Ordering::SeqCst), 1);
}

#[test]
fn a_chain_links_on_the_pool_in_the_last_batch_or_over_the_limit() {
    let noop = |_: u64| {};
    let held = XHeld::with_limit(10);
    {
        let seed_end = XSeedEnd::new(3, &noop, &held, false);
        let mut a = XSeedGuard::new(Some(&seed_end));
        let mut b = XSeedGuard::new(Some(&seed_end));
        assert!(!a.end(3), "room left, a later batch may run");
        assert!(b.end(8), "11 HSPs held, over the limit");
        assert_eq!(held.now.load(Ordering::SeqCst), 11);
        drop(a);
        assert_eq!(held.now.load(Ordering::SeqCst), 8);
        // A subject that returns before its seed stage ends holds nothing.
        drop(XSeedGuard::new(Some(&seed_end)));
        drop(b);
    }
    assert_eq!(held.now.load(Ordering::SeqCst), 0);
    assert_eq!(held.most.load(Ordering::SeqCst), 11);
    let last = XSeedEnd::new(1, &noop, &held, true);
    assert!(XSeedGuard::new(Some(&last)).end(1), "the plan's last batch");
    assert!(XSeedGuard::new(None).end(1), "no pipe");
}

#[test]
fn a_panic_in_an_early_batch_is_raised_on_the_loop() {
    for threads in [1, 4, 8, 16] {
        let plan = vec![0..2, 2..4, 4..6];
        let running = AtomicUsize::new(0);
        let most = AtomicUsize::new(0);
        let inner = fake_search(&running, &most, 3, None);
        let search = |range: std::ops::Range<usize>, seed_end: Option<&XSeedEnd<'_>>| {
            let start = range.start;
            let result = inner(range, seed_end);
            assert!(start != 2, "batch {start} panicked");
            result
        };
        let pipe = XPipe::with_limit(plan.clone(), 3, &search, u64::MAX).expect("three batches");
        let caught = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            pool(threads).install(|| {
                rayon::in_place_scope(|scope| {
                    let _closer = pipe.closer();
                    for range in plan.iter().cloned() {
                        let _ = pipe.take(scope, range);
                    }
                })
            })
        }));
        let panic = caught.expect_err("the loop raises the batch's panic");
        let message = panic.downcast_ref::<String>().cloned().unwrap_or_default();
        assert!(
            message.contains("batch 2 panicked"),
            "threads {threads}: {message}"
        );
        assert_eq!(
            pipe.counts.taken.load(Ordering::SeqCst),
            2,
            "threads {threads}"
        );
        // Batch 1 may have started batch 2 before it panicked.
        assert!(
            pipe.counts.early.load(Ordering::SeqCst) >= 1,
            "threads {threads}"
        );
    }
}

#[test]
fn after_the_loop_ends_no_batch_starts_early() {
    let plan = vec![0..2, 2..4, 4..6];
    let running = AtomicUsize::new(0);
    let most = AtomicUsize::new(0);
    let search = fake_search(&running, &most, 2, None);
    let pipe = XPipe::with_limit(plan, 2, &search, u64::MAX).expect("three batches");
    pool(2).install(|| {
        rayon::in_place_scope(|scope| {
            drop(pipe.closer());
            pipe.seed_ended(scope, 0);
        })
    });
    assert_eq!(pipe.counts.early.load(Ordering::SeqCst), 0);
    assert!(!pipe.started[1].load(Ordering::SeqCst));
}

#[test]
fn a_plan_of_one_batch_makes_no_pipe() {
    let running = AtomicUsize::new(0);
    let most = AtomicUsize::new(0);
    let search = fake_search(&running, &most, 2, None);
    assert!(XPipe::with_limit(vec![0..1], 2, &search, u64::MAX).is_none());
}
