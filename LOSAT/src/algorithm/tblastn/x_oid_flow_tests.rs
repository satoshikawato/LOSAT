use super::*;
use std::sync::atomic::{AtomicUsize, Ordering};

/// The waits the tests run: the flow's condition only, and with the pool's jobs (100 us).
const HELPS: [Option<Duration>; 2] = [None, Some(Duration::from_micros(100))];

/// Runs a flow as `run_preliminary_subjects` does: the owner takes back the lists in order
/// inside a scope of the pool whose other threads run `work`. Returns the subjects in the
/// order the owner received them and the peak of waiting lists.
fn run_flow(
    threads: usize,
    chunks: Vec<Range<usize>>,
    window: usize,
    help: Option<Duration>,
    search: &(impl Fn(usize) -> usize + Sync),
) -> (Vec<usize>, usize) {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    let flow = XOidFlow::new(chunks, window, threads - 1).x_with_help(help);
    let received = pool.install(|| {
        rayon::scope(|scope| {
            for _ in 1..threads {
                scope.spawn(|_| flow.work(search));
            }
            let _close = XCloseOnDrop(&flow);
            let mut received = Vec::new();
            while let Some(lists) = flow.next_in_order(search).unwrap() {
                received.extend(lists);
            }
            received
        })
    });
    (received, flow.peak_waiting())
}

// No NCBI counterpart: the owner receives every subject once and in OID order whatever order
// the chunks finish in, and never more than `window` finished lists wait.
#[test]
fn the_owner_takes_back_the_lists_in_oid_order_within_the_window() {
    // Chunks of 1 to 3 subjects; the work of a subject varies so that later chunks often
    // finish first.
    let mut chunks = Vec::new();
    let mut start = 0;
    for size in [1, 3, 1, 2, 2, 1, 3, 1, 1, 2, 3, 1, 1, 1, 2, 1, 3, 2, 1, 1] {
        chunks.push(start..start + size);
        start += size;
    }
    let subjects = start;
    let search = |subject: usize| {
        let spin = (subject * 7919) % 13;
        std::thread::sleep(std::time::Duration::from_micros(200 * spin as u64));
        subject
    };
    for help in HELPS {
        for (threads, window) in [(1, 1), (2, 1), (2, 2), (3, 4), (4, 2), (8, 4), (8, 32)] {
            let (received, peak) = run_flow(threads, chunks.clone(), window, help, &search);
            assert_eq!(
                received,
                (0..subjects).collect::<Vec<_>>(),
                "{threads} threads, window {window}, help {help:?}"
            );
            assert!(
                peak <= window,
                "{threads} threads, window {window}, help {help:?}: peak {peak}"
            );
        }
    }
}

// No NCBI counterpart: searches that run parts of their work on the pool (as the units of a
// long subject beyond a chunk's budget) still give every subject once and in OID order, with
// the smallest windows too: a thread waiting for a part runs other jobs of the pool, and the
// threads start together, so that no loop of the flow is left in the pool to start there.
#[test]
fn searches_with_parts_on_the_pool_keep_the_order() {
    use rayon::prelude::*;
    let chunks: Vec<Range<usize>> = (0..24).map(|i| i..i + 1).collect();
    let search = |subject: usize| {
        let parts: usize = (0..6)
            .into_par_iter()
            .with_max_len(1)
            .map(|part| {
                std::thread::sleep(std::time::Duration::from_micros(
                    100 * ((subject + part) % 5) as u64,
                ));
                part
            })
            .sum();
        assert_eq!(parts, 15);
        subject
    };
    for help in HELPS {
        for (threads, window) in [(2, 1), (3, 1), (4, 1), (4, 2), (8, 1), (8, 4)] {
            let (received, peak) = run_flow(threads, chunks.clone(), window, help, &search);
            assert_eq!(
                received,
                (0..24).collect::<Vec<_>>(),
                "{threads} threads, window {window}, help {help:?}"
            );
            assert!(
                peak <= window,
                "{threads} threads, window {window}, help {help:?}: peak {peak}"
            );
        }
    }
}

// No NCBI counterpart (LOSAT_X_TBNFLOWHELP): the owner took a short subject while a long one
// runs its parts on the pool, as on a genome database, and every chunk it could take is taken.
// With the help the owner runs some of those parts while it waits, without it none; the
// subjects come back once and in OID order either way. On two threads the other thread takes
// chunk 0: it completes the start and takes the next chunk under the same lock.
#[test]
fn the_waiting_owner_runs_the_parts_of_a_long_subject() {
    use rayon::prelude::*;
    let owner = AtomicUsize::new(usize::MAX);
    let owner_parts = AtomicUsize::new(0);
    let parts_run = AtomicUsize::new(0);
    let wait_for_owner = std::sync::atomic::AtomicBool::new(false);
    // Subject 0 is long (64 parts of 1 ms on the pool); the others are short. With the help
    // the other thread's parts wait (10 s at most) until the owner has run one, so that a
    // slow owner cannot find them all run.
    let search = |subject: usize| {
        if subject == 0 {
            (0..64).into_par_iter().with_max_len(1).for_each(|_| {
                std::thread::sleep(std::time::Duration::from_millis(1));
                parts_run.fetch_add(1, Ordering::SeqCst);
                if rayon::current_thread_index() == Some(owner.load(Ordering::SeqCst)) {
                    owner_parts.fetch_add(1, Ordering::SeqCst);
                } else if wait_for_owner.load(Ordering::SeqCst) {
                    let start = Instant::now();
                    while owner_parts.load(Ordering::SeqCst) == 0
                        && start.elapsed() < Duration::from_secs(10)
                    {
                        std::thread::sleep(std::time::Duration::from_millis(1));
                    }
                }
            });
        }
        subject
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(2)
        .build()
        .unwrap();
    for help in HELPS {
        owner_parts.store(0, Ordering::SeqCst);
        parts_run.store(0, Ordering::SeqCst);
        wait_for_owner.store(help.is_some(), Ordering::SeqCst);
        let flow = XOidFlow::new((0..12).map(|i| i..i + 1).collect(), 2, 1).x_with_help(help);
        let received = pool.install(|| {
            rayon::scope(|scope| {
                owner.store(rayon::current_thread_index().unwrap(), Ordering::SeqCst);
                scope.spawn(|_| flow.work(&search));
                let _close = XCloseOnDrop(&flow);
                let mut received = Vec::new();
                while let Some(lists) = flow.next_in_order(&search).unwrap() {
                    received.extend(lists);
                }
                received
            })
        });
        assert_eq!(received, (0..12).collect::<Vec<_>>(), "help {help:?}");
        assert_eq!(parts_run.load(Ordering::SeqCst), 64, "help {help:?}");
        let by_owner = owner_parts.load(Ordering::SeqCst);
        if help.is_some() {
            assert!(by_owner > 0, "the owner ran no part");
        } else {
            assert_eq!(by_owner, 0);
        }
    }
}

// LOSAT_X_TBNLOOPFREE: every thread of the pool, and the owner after its reduction, ends its
// loop over the chunks without the buffers it kept between subjects (here the kept frames).
#[test]
fn a_thread_frees_its_kept_frames_when_its_loop_ends() {
    use crate::algorithm::tblastn::x_subject_translate::{x_kept_frames, x_reuse_frame};
    // Nothing to check with LOSAT_X_TBNLOOPFREE=0, as in the gate's arm of the switches off.
    if !x_loop_free_on() {
        return;
    }
    let threads = 4;
    let search = |subject: usize| {
        x_reuse_frame(vec![0; 64]);
        assert!(x_kept_frames() > 0);
        subject
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    let flow = XOidFlow::new((0..64).map(|i| i..i + 1).collect(), 8, threads - 1);
    let received = pool.install(|| {
        let received = rayon::scope(|scope| {
            for _ in 1..threads {
                scope.spawn(|_| flow.work(&search));
            }
            let _close = XCloseOnDrop(&flow);
            let mut received = Vec::new();
            while let Some(lists) = flow.next_in_order(&search).unwrap() {
                received.extend(lists);
            }
            received
        });
        x_loop_free();
        received
    });
    assert_eq!(received, (0..64).collect::<Vec<_>>());
    assert_eq!(pool.broadcast(|_| x_kept_frames()), vec![0; threads]);
}

// No NCBI counterpart: when the owner stops early (an error in the reduction) the other
// threads stop taking chunks and the scope ends.
#[test]
fn closing_stops_the_threads() {
    let searched = AtomicUsize::new(0);
    let search = |subject: usize| {
        searched.fetch_add(1, Ordering::SeqCst);
        std::thread::sleep(std::time::Duration::from_millis(1));
        subject
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap();
    let flow = XOidFlow::new((0..1000).map(|i| i..i + 1).collect(), 8, 3);
    pool.install(|| {
        rayon::scope(|scope| {
            for _ in 1..4 {
                scope.spawn(|_| flow.work(&search));
            }
            let _close = XCloseOnDrop(&flow);
            for _ in 0..3 {
                flow.next_in_order(&search).unwrap().unwrap();
            }
        })
    });
    // The chunks taken before the close lie within the window after the three taken back.
    assert!(searched.load(Ordering::SeqCst) <= 3 + 8);
}

// No NCBI counterpart: a panic in a search, on a thread of the pool or on the owner, ends the
// scope with that panic instead of leaving the owner or the threads waiting (a thread's panic
// fails the owner's wait; the owner's panic closes the flow).
#[test]
fn a_panicking_search_ends_the_scope_with_its_panic() {
    let search = |subject: usize| {
        std::thread::sleep(std::time::Duration::from_millis(2));
        if subject == 5 {
            panic!("test panic in subject 5");
        }
        subject
    };
    for threads in [1, 2, 4] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let flow = XOidFlow::new((0..20).map(|i| i..i + 1).collect(), 4, threads - 1);
        let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            pool.install(|| {
                rayon::scope(|scope| {
                    for _ in 1..threads {
                        scope.spawn(|_| flow.work(&search));
                    }
                    let _close = XCloseOnDrop(&flow);
                    while flow.next_in_order(&search)?.is_some() {}
                    anyhow::Ok(())
                })
            })
        }));
        let panic = outcome.expect_err("the panic reaches the caller");
        assert_eq!(
            panic.downcast_ref::<&str>().copied(),
            Some("test panic in subject 5"),
            "{threads} threads"
        );
    }
}

// Probe (run by hand with `--ignored`, after line B's B17 probe in `x_bnstream.rs`): many
// searches whose first subject hands out parts on the pool at once, on more threads than
// CPUs (run it under `taskset -c 0,1`), with the waits that run the pool's jobs. A search that
// makes no progress for LOSAT_B11C_STALL_S seconds ends the process (exit 3).
#[test]
#[ignore]
fn probe_b11c_many_flows_with_a_long_first_subject() {
    use rayon::prelude::*;
    let var = |name: &str, default: usize| {
        std::env::var(name)
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(default)
    };
    let threads = var("LOSAT_B11C_THREADS", 32);
    let passes = var("LOSAT_B11C_PASSES", 20_000);
    let stall = var("LOSAT_B11C_STALL_S", 20) as u64;
    let help = match var("LOSAT_B11C_HELP_US", 100) {
        0 => None,
        us => Some(Duration::from_micros(us as u64)),
    };
    let done = AtomicUsize::new(0);
    let start = Instant::now();
    let light = |x: usize| (0..200u64).fold(x as u64, |a, y| a.wrapping_mul(31).wrapping_add(y));
    let search = |subject: usize| -> usize {
        if subject == 0 {
            (0..256)
                .into_par_iter()
                .map(light)
                .reduce(|| 0, u64::wrapping_add) as usize
                % 7
        } else {
            light(subject) as usize % 7
        }
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
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
                        "B11C HANG at pass {now} of {passes} after {:.1} s (threads {threads} help {help:?})",
                        start.elapsed().as_secs_f64()
                    );
                    std::process::exit(3);
                }
            }
        });
        pool.install(|| {
            for p in 0..passes {
                let flow = XOidFlow::new((0..16).map(|i| i..i + 1).collect(), 2, threads - 1)
                    .x_with_help(help);
                let received = rayon::scope(|scope| {
                    for _ in 1..threads {
                        scope.spawn(|_| flow.work(&search));
                    }
                    let _close = XCloseOnDrop(&flow);
                    let mut received = 0;
                    while let Some(lists) = flow.next_in_order(&search).unwrap() {
                        received += lists.len();
                    }
                    received
                });
                assert_eq!(received, 16);
                done.store(p + 1, Ordering::Relaxed);
            }
        });
    });
    eprintln!(
        "B11C {passes} passes in {:.3} s (threads {threads} help {help:?})",
        start.elapsed().as_secs_f64()
    );
}
