use super::{search_in_order, search_in_order_bounded, WaitingBound};
use std::sync::atomic::{AtomicUsize, Ordering};

fn pool(threads: usize) -> rayon::ThreadPool {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("a test pool")
}

// Uneven work per item (a few slow items), so that results arrive out of order.
fn work(i: usize) -> u64 {
    let rounds = if i % 37 == 5 {
        20_000
    } else {
        50 + (i * 7919) % 300
    };
    (0..rounds as u64).fold(i as u64, |acc, x| acc.wrapping_mul(31).wrapping_add(x))
}

#[test]
fn hands_over_every_item_once_in_item_order() {
    for (threads, n, window) in [
        (1, 10, 4),
        (4, 0, 8),
        (4, 1, 8),
        (4, 500, 1),
        (4, 500, 16),
        (8, 1000, 64),
        (3, 7, 100),
    ] {
        for assist in [false, true] {
            let handed = pool(threads).install(|| {
                let mut handed: Vec<(usize, u64)> = Vec::new();
                search_in_order(
                    n,
                    window,
                    assist,
                    &|| 0usize,
                    &|calls: &mut usize, i| {
                        *calls += 1;
                        work(i)
                    },
                    &mut |i, value| handed.push((i, value)),
                );
                handed
            });
            let want: Vec<(usize, u64)> = (0..n).map(|i| (i, work(i))).collect();
            assert_eq!(
                handed, want,
                "threads {threads} n {n} window {window} assist {assist}"
            );
        }
    }
}

// A search with a nested parallel stage (as the speculative traceback has): a thread that
// waits for its nested work must not pick up a second loop of this pass.
#[test]
fn a_nested_parallel_stage_does_not_block_the_pass() {
    use rayon::prelude::*;
    let threads = 4;
    let nested = |i: usize| -> u64 {
        (0..64usize)
            .into_par_iter()
            .map(|x| work(i + x))
            .reduce(|| 0, u64::wrapping_add)
    };
    let want: Vec<(usize, u64)> = (0..400).map(|i| (i, nested(i))).collect();
    for assist in [false, true] {
        let handed = pool(threads).install(|| {
            let mut handed: Vec<(usize, u64)> = Vec::new();
            search_in_order(
                400,
                8,
                assist,
                &|| (),
                &|_: &mut (), i| nested(i),
                &mut |i, value| handed.push((i, value)),
            );
            handed
        });
        assert_eq!(handed, want, "assist {assist}");
    }
}

// A thread of the pool that is busy when the pass starts holds back every item until it has
// started its loop (`wait_for_every_loop`): before then it could take a nested job of an item
// and start its loop inside the wait for that job.
#[test]
fn no_item_is_searched_before_every_thread_has_started_its_loop() {
    use std::sync::atomic::AtomicBool;
    use std::sync::Arc;
    use std::time::Duration;
    let pool = pool(4);
    let busy = Arc::new(AtomicBool::new(false));
    let release = Arc::new(AtomicBool::new(false));
    {
        let (busy, release) = (Arc::clone(&busy), Arc::clone(&release));
        pool.spawn(move || {
            busy.store(true, Ordering::SeqCst);
            while !release.load(Ordering::SeqCst) {
                std::hint::spin_loop();
            }
        });
    }
    while !busy.load(Ordering::SeqCst) {
        std::thread::yield_now();
    }
    let releaser = {
        let release = Arc::clone(&release);
        std::thread::spawn(move || {
            std::thread::sleep(Duration::from_millis(200));
            release.store(true, Ordering::SeqCst);
        })
    };
    // items searched while the busy thread had not been released
    let early = AtomicUsize::new(0);
    let handed = pool.install(|| {
        let mut handed = Vec::new();
        search_in_order(
            32,
            4,
            true,
            &|| (),
            &|_: &mut (), i| {
                if !release.load(Ordering::SeqCst) {
                    early.fetch_add(1, Ordering::SeqCst);
                }
                work(i)
            },
            &mut |i, _| handed.push(i),
        );
        handed
    });
    releaser.join().unwrap();
    assert_eq!(handed, (0..32).collect::<Vec<_>>());
    assert_eq!(
        early.load(Ordering::SeqCst),
        0,
        "items were searched before the busy thread started its loop"
    );
}

// With `assist`, the threads that wait for the window run the chunks that a slow item's
// search hands out: a pass whose first item is long and parallel (the long subject of a `-db`
// search) ends with the other threads' help, in item order.
//
// The first chunk (run inline by the thread that searches item 0, while the other chunks wait
// in its deque) holds until another thread has run a chunk, so that the test does not depend
// on the other threads being scheduled before item 0 ends on a loaded host (CI ran it at ~50 %
// failures). Without the assist the other threads never take a chunk, and the deadline ends
// the wait so that the assertion fails instead of the test hanging.
#[test]
fn waiting_threads_run_a_long_items_nested_jobs_with_assist() {
    use rayon::prelude::*;
    use std::collections::HashSet;
    use std::sync::Mutex;
    use std::time::{Duration, Instant};
    let threads = 4;
    let runners: Mutex<HashSet<usize>> = Mutex::new(HashSet::new());
    let search = |_: &mut (), i: usize| -> u64 {
        if i == 0 {
            (0..256usize)
                .into_par_iter()
                .map(|x| {
                    runners
                        .lock()
                        .unwrap()
                        .insert(rayon::current_thread_index().unwrap());
                    if x == 0 {
                        let deadline = Instant::now() + Duration::from_secs(30);
                        while runners.lock().unwrap().len() < 2 && Instant::now() < deadline {
                            std::thread::yield_now();
                        }
                    }
                    (0..20_000u64).fold(x as u64, |a, y| a.wrapping_mul(31).wrapping_add(y))
                })
                .reduce(|| 0, u64::wrapping_add)
        } else {
            work(i)
        }
    };
    let handed = pool(threads).install(|| {
        let mut handed: Vec<usize> = Vec::new();
        search_in_order(64, 8, true, &|| (), &search, &mut |i, _| handed.push(i));
        handed
    });
    assert_eq!(handed, (0..64).collect::<Vec<_>>());
    assert!(
        runners.lock().unwrap().len() > 1,
        "the long item ran on one thread"
    );
}

#[test]
fn claims_stay_within_the_window() {
    let threads = 6;
    let window = 12;
    let handed = AtomicUsize::new(0);
    let worst = AtomicUsize::new(0);
    pool(threads).install(|| {
        search_in_order(
            2000,
            window,
            false,
            &|| (),
            &|_: &mut (), i| {
                let ahead = i - handed.load(Ordering::SeqCst).min(i);
                worst.fetch_max(ahead, Ordering::SeqCst);
                work(i)
            },
            &mut |i, _| handed.store(i + 1, Ordering::SeqCst),
        );
    });
    // The owner counts an item as handed over just before `hand_over` runs.
    assert!(worst.load(Ordering::SeqCst) <= window);
}

#[test]
fn makes_one_state_per_thread_at_most() {
    let threads = 4;
    let states = AtomicUsize::new(0);
    pool(threads).install(|| {
        search_in_order(
            300,
            32,
            false,
            &|| {
                states.fetch_add(1, Ordering::SeqCst);
            },
            &|_: &mut (), i| work(i),
            &mut |_, _| {},
        );
    });
    let made = states.load(Ordering::SeqCst);
    assert!((1..=threads).contains(&made), "{made} states");
}

#[test]
fn a_failed_search_stops_the_pass() {
    for failing in [0usize, 3, 299] {
        let outcome = std::panic::catch_unwind(|| {
            pool(4).install(|| {
                search_in_order(
                    300,
                    16,
                    false,
                    &|| (),
                    &|_: &mut (), i| {
                        assert!(i != failing, "item {i} fails");
                        work(i)
                    },
                    &mut |_, _| {},
                );
            })
        });
        assert!(outcome.is_err(), "item {failing}");
    }
}

// While the bytes are above the bound, no item past `floor` is claimed until the owner hands
// the first ones over; below it, the claims go on to the whole window. Item 0 holds the start
// of the window until the other threads have gone as far as they may.
#[test]
fn the_bytes_narrow_the_window_to_its_floor() {
    const N: usize = 64;
    const WINDOW: usize = 32;
    const FLOOR: usize = 8;
    for (waiting, last_while_held) in [(1usize, FLOOR - 1), (0, WINDOW - 1)] {
        let waiting_bytes = AtomicUsize::new(waiting);
        let highest = AtomicUsize::new(0);
        let mut handed = Vec::new();
        pool(4).install(|| {
            search_in_order_bounded(
                N,
                WINDOW,
                true,
                Some(WaitingBound {
                    floor: FLOOR,
                    waiting_bytes: &waiting_bytes,
                    max_bytes: 0,
                }),
                &|| (),
                &|_: &mut (), i| {
                    highest.fetch_max(i, Ordering::SeqCst);
                    if i == 0 {
                        let deadline =
                            std::time::Instant::now() + std::time::Duration::from_secs(30);
                        while highest.load(Ordering::SeqCst) < last_while_held
                            && std::time::Instant::now() < deadline
                        {
                            std::thread::yield_now();
                        }
                        std::thread::sleep(std::time::Duration::from_millis(50));
                        assert_eq!(highest.load(Ordering::SeqCst), last_while_held);
                    }
                    i
                },
                &mut |i, value| {
                    assert_eq!(i, value);
                    handed.push(i);
                },
            );
        });
        assert_eq!(handed, (0..N).collect::<Vec<_>>(), "waiting {waiting}");
    }
}

#[test]
fn runs_on_one_thread_outside_a_pool() {
    let mut handed = Vec::new();
    search_in_order(
        20,
        4,
        false,
        &|| (),
        &|_: &mut (), i| i * 2,
        &mut |i, value| handed.push((i, value)),
    );
    assert_eq!(handed, (0..20).map(|i| (i, i * 2)).collect::<Vec<_>>());
}

// B17 probe (run by hand with `--ignored`): many passes whose first item hands out nested jobs
// at once, on more threads than CPUs. A pool thread that steals one of those jobs before it
// has started its loop runs the loop nested inside the stolen job's join wait (rayon-core 1.13
// `take_local_job` pops the thread's own broadcast deque there), above unfinished work of item
// 0. A pass that makes no progress for LOSAT_B17_STALL_S seconds ends the process (exit 3).
#[test]
#[ignore]
fn probe_b17_many_passes_with_a_nested_first_item() {
    use rayon::prelude::*;
    use std::time::{Duration, Instant};
    let var = |name: &str, default: usize| {
        std::env::var(name)
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(default)
    };
    let threads = var("LOSAT_B17_THREADS", 32);
    let passes = var("LOSAT_B17_PASSES", 100_000);
    let stall = var("LOSAT_B17_STALL_S", 20) as u64;
    let assist = var("LOSAT_B17_ASSIST", 1) == 1;
    let chunks = var("LOSAT_B17_CHUNKS", 256);
    let done = AtomicUsize::new(0);
    let start = Instant::now();
    let light = |x: usize| (0..200u64).fold(x as u64, |a, y| a.wrapping_mul(31).wrapping_add(y));
    std::thread::scope(|s| {
        s.spawn(|| {
            let mut last = 0;
            let mut since = Instant::now();
            loop {
                std::thread::sleep(Duration::from_millis(200));
                let now = done.load(Ordering::Relaxed);
                if now >= passes {
                    return;
                }
                if now != last {
                    last = now;
                    since = Instant::now();
                } else if since.elapsed() > Duration::from_secs(stall) {
                    eprintln!(
                        "B17 HANG at pass {now} of {passes} after {:.1} s (threads {threads} assist {assist})",
                        start.elapsed().as_secs_f64()
                    );
                    std::process::exit(3);
                }
            }
        });
        pool(threads).install(|| {
            for p in 0..passes {
                search_in_order(
                    16,
                    2,
                    assist,
                    &|| (),
                    &|_: &mut (), i: usize| -> u64 {
                        if i == 0 {
                            (0..chunks)
                                .into_par_iter()
                                .map(light)
                                .reduce(|| 0, u64::wrapping_add)
                        } else {
                            light(i)
                        }
                    },
                    &mut |_, _| {},
                );
                done.store(p + 1, Ordering::Relaxed);
            }
        });
    });
    eprintln!(
        "B17 {passes} passes in {:.3} s (threads {threads} assist {assist})",
        start.elapsed().as_secs_f64()
    );
}
