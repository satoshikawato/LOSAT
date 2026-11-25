use super::*;
use std::sync::atomic::AtomicBool;
use std::time::{Duration, Instant};

fn work(i: usize) -> u64 {
    let rounds = 500 + (i * 7_919) % 5_000;
    (0..rounds as u64).fold(i as u64, |acc, x| {
        acc.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(x)
    })
}

fn pool(threads: usize) -> rayon::ThreadPool {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap()
}

// Every item runs once and the values come back in index order: outside a pool and in pools of
// 1, 2, 4 and 8 threads, for sets of none, one, a few and many items.
#[test]
fn every_item_runs_once_in_index_order() {
    for n in [0, 1, 2, 3, 17, 200] {
        let want: Vec<u64> = (0..n).map(work).collect();
        let runs: Vec<AtomicUsize> = (0..n).map(|_| AtomicUsize::new(0)).collect();
        let got = x_map(n, |i| {
            runs[i].fetch_add(1, Ordering::Relaxed);
            work(i)
        });
        assert_eq!(got, want, "outside a pool, n {n}");
        assert!(runs.iter().all(|r| r.load(Ordering::Relaxed) == 1));
        for threads in [1, 2, 4, 8] {
            let runs: Vec<AtomicUsize> = (0..n).map(|_| AtomicUsize::new(0)).collect();
            let got = pool(threads).install(|| {
                x_map(n, |i| {
                    runs[i].fetch_add(1, Ordering::Relaxed);
                    work(i)
                })
            });
            assert_eq!(got, want, "{threads} threads, n {n}");
            assert!(runs.iter().all(|r| r.load(Ordering::Relaxed) == 1));
        }
    }
}

// Item 1 is taken by the other thread of a pool of two and runs for a while; meanwhile the owner,
// which has run item 0, has a job of the pool in its own queue (pushed by item 0). A thread that
// waited in the pool (`join`, `par_iter`) would run that job before it returned; the owner blocks
// instead, and the job runs on the other thread once item 1 has ended (or after the owner returned).
#[test]
fn the_owner_runs_no_pool_job_while_it_waits() {
    let pool = pool(2);
    pool.install(|| {
        let owner = rayon::current_thread_index();
        let in_map = Arc::new(AtomicBool::new(false));
        let ran_on_owner_in_map = Arc::new(AtomicBool::new(false));
        let item1_started = AtomicBool::new(false);
        let job_spawned = AtomicBool::new(false);
        in_map.store(true, Ordering::SeqCst);
        let got = x_map(2, |i| {
            if i == 0 {
                // Wait (up to 10 s) for the other thread to take item 1, then push a job of the
                // pool on this thread's queue.
                let start = Instant::now();
                while !item1_started.load(Ordering::SeqCst)
                    && start.elapsed() < Duration::from_secs(10)
                {
                    std::thread::yield_now();
                }
                let (in_map, ran) = (Arc::clone(&in_map), Arc::clone(&ran_on_owner_in_map));
                rayon::spawn(move || {
                    if rayon::current_thread_index() == owner && in_map.load(Ordering::SeqCst) {
                        ran.store(true, Ordering::SeqCst);
                    }
                });
                job_spawned.store(true, Ordering::SeqCst);
            } else {
                assert_ne!(rayon::current_thread_index(), owner);
                item1_started.store(true, Ordering::SeqCst);
                std::thread::sleep(Duration::from_millis(100));
            }
            i
        });
        in_map.store(false, Ordering::SeqCst);
        assert_eq!(got, vec![0, 1]);
        assert!(item1_started.load(Ordering::SeqCst) && job_spawned.load(Ordering::SeqCst));
        assert!(!ran_on_owner_in_map.load(Ordering::SeqCst));
    });
}

// A panic of an item is raised by the owner, whether the owner ran the item (item 0) or a helper
// job did (item 1, taken by the other thread while the owner runs item 0); the owner waits for
// the items the helpers took before it raises it.
#[test]
fn a_panicking_item_is_raised_by_the_owner() {
    let pool = pool(2);
    pool.install(|| {
        let result =
            std::panic::catch_unwind(|| x_map(4, |i| if i == 0 { panic!("item 0") } else { i }));
        let panic = result.expect_err("the owner's item panicked");
        assert_eq!(panic.downcast_ref::<&str>(), Some(&"item 0"));
        let owner = rayon::current_thread_index();
        let item1_started = AtomicBool::new(false);
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            x_map(2, |i| {
                if i == 0 {
                    let start = Instant::now();
                    while !item1_started.load(Ordering::SeqCst)
                        && start.elapsed() < Duration::from_secs(10)
                    {
                        std::thread::yield_now();
                    }
                    0
                } else {
                    item1_started.store(true, Ordering::SeqCst);
                    if rayon::current_thread_index() != owner {
                        panic!("item 1 on a helper");
                    }
                    1
                }
            })
        }));
        let panic = result.expect_err("the helper's item panicked");
        assert_eq!(panic.downcast_ref::<&str>(), Some(&"item 1 on a helper"));
    });
}

// An item's own set runs on the item's thread, one item after the other (an item makes no pool
// call), and gives the same values.
#[test]
fn an_items_own_set_runs_on_its_thread() {
    let got = pool(4).install(|| {
        x_map(8, |i| {
            let here = rayon::current_thread_index();
            let inner = x_map(16, |j| {
                assert_eq!(rayon::current_thread_index(), here);
                assert!(in_item());
                work(i * 16 + j)
            });
            inner.into_iter().fold(0u64, u64::wrapping_add)
        })
    });
    let want: Vec<u64> = (0..8)
        .map(|i| {
            (0..16)
                .map(|j| work(i * 16 + j))
                .fold(0u64, u64::wrapping_add)
        })
        .collect();
    assert_eq!(got, want);
    assert!(!in_item());
}
