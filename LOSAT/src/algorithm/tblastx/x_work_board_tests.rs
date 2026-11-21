use super::*;
use std::sync::atomic::AtomicU64;

fn work(s: usize, i: usize) -> u64 {
    let rounds = 2_000 + ((s * 7_919 + i * 104_729) % 20_000);
    (0..rounds as u64).fold((s * 1_000 + i) as u64, |acc, x| {
        acc.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(x)
    })
}

// Subjects that put their units on the board from the threads of a broadcast (the TBXNEXT
// loop), the other threads helping, and each subject linking a few groups on the board after its
// units: every unit and group runs once, each subject gets its own results in unit order, and
// the run ends at 1, 2, 4 and 8 threads with fewer, as many and more subjects than threads.
#[test]
fn every_item_runs_once_and_every_owner_gets_its_results() {
    for threads in [1usize, 2, 4, 8] {
        for subjects in [1usize, 2, 3, threads, threads + 3] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("test pool");
            let board = XBoard::default();
            let batch = board.batch(subjects);
            let next = AtomicUsize::new(0);
            let ran = AtomicU64::new(0);
            let checked = AtomicUsize::new(0);
            pool.install(|| {
                rayon::broadcast(|_| {
                    loop {
                        let s = next.fetch_add(1, Ordering::Relaxed);
                        if s >= subjects {
                            break;
                        }
                        let _entered = enter(&board, &batch);
                        if s % 5 == 4 {
                            // A subject too short to search: no units.
                            batch.close_without_units(s);
                            checked.fetch_add(1, Ordering::Relaxed);
                            continue;
                        }
                        let units = 3 + (s * 11) % 40;
                        let out: Vec<OnceLock<u64>> = (0..units).map(|_| OnceLock::new()).collect();
                        board.run_units(&batch, s, units, &|i, taker: &mut XTaker| {
                            taker.scratch.get_or_insert_with(|| Box::new(0u64));
                            ran.fetch_add(1, Ordering::Relaxed);
                            out[i].set(work(s, i)).expect("a unit runs once");
                        });
                        for (i, cell) in out.iter().enumerate() {
                            assert_eq!(cell.get(), Some(&work(s, i)), "subject {s} unit {i}");
                        }
                        let groups: Vec<OnceLock<u64>> = (0..4).map(|_| OnceLock::new()).collect();
                        board.run_groups(4, &|g, _| {
                            groups[g].set(work(s + 100, g)).expect("a group runs once");
                        });
                        for (g, cell) in groups.iter().enumerate() {
                            assert_eq!(cell.get(), Some(&work(s + 100, g)));
                        }
                        checked.fetch_add(1, Ordering::Relaxed);
                    }
                    board.help(&batch, &AtomicBool::new(false));
                });
            });
            board.end_batch(&batch);
            let expected: u64 = (0..subjects)
                .filter(|s| s % 5 != 4)
                .map(|s| (3 + (s * 11) % 40) as u64)
                .sum();
            assert_eq!(
                ran.load(Ordering::Relaxed),
                expected,
                "{threads} threads, {subjects} subjects"
            );
            assert_eq!(checked.load(Ordering::Relaxed), subjects);
        }
    }
}

// `drained` reads each slot once per batch (its cursor never goes back) and is true only when
// every slot is closed and taken; a slot closed without units does not hold the unit order's
// hint (`lowest`) back.
#[test]
fn drained_and_lowest_pass_closed_slots_once() {
    let board = XBoard::default();
    let batch = board.batch(4);
    batch.close_without_units(0);
    batch.close_without_units(2);
    assert!(!batch.drained());
    assert_eq!(batch.drained_below.load(Ordering::Relaxed), 1);
    let ran = AtomicUsize::new(0);
    board.run_units(&batch, 1, 3, &|_, _: &mut XTaker| {
        ran.fetch_add(1, Ordering::Relaxed);
    });
    assert_eq!(ran.load(Ordering::Relaxed), 3);
    assert!(
        batch.slots[1].lock().unwrap().is_none(),
        "the slot is cleared"
    );
    assert!(!batch.drained());
    assert_eq!(batch.drained_below.load(Ordering::Relaxed), 3);
    assert!(!batch.run_one(&mut XTaker::default()));
    assert_eq!(batch.lowest.load(Ordering::Relaxed), 3);
    batch.close_without_units(3);
    assert!(batch.drained());
    assert_eq!(batch.drained_below.load(Ordering::Relaxed), 4);
}

// An item that panics on another thread counts as run, so its owner returns, and the panic
// reaches the caller (a hang here would stop the test at the watchdog).
#[test]
fn a_panicking_item_on_another_thread_does_not_hang_its_owner() {
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .expect("test pool");
        let board = XBoard::default();
        let other_took = AtomicBool::new(false);
        let owner_done = AtomicBool::new(false);
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            pool.install(|| {
                rayon::broadcast(|ctx| {
                    if ctx.index() == 0 {
                        board.run_set(
                            4,
                            &|_, _| {
                                if rayon::current_thread_index() != Some(0) {
                                    other_took.store(true, Ordering::SeqCst);
                                    panic!("an item panics on the other thread");
                                }
                                let start = std::time::Instant::now();
                                while !other_took.load(Ordering::SeqCst)
                                    && start.elapsed() < std::time::Duration::from_secs(5)
                                {
                                    std::hint::spin_loop();
                                }
                            },
                            false,
                        );
                        owner_done.store(true, Ordering::SeqCst);
                    } else {
                        let mut taker = XTaker::default();
                        while !owner_done.load(Ordering::SeqCst) {
                            if !board.take_one(&mut taker) {
                                std::thread::yield_now();
                            }
                        }
                    }
                });
            });
        }));
        let _ = tx.send((result.is_err(), other_took.load(Ordering::SeqCst)));
    });
    let (panicked, other_took) = rx
        .recv_timeout(std::time::Duration::from_secs(60))
        .expect("the owner returned");
    if other_took {
        assert!(panicked, "the item's panic reaches the caller");
    }
}

// An owner whose own item panics while another thread runs an item of its set waits for that
// item (it borrows the owner's frame), then unwinds: an item that ends below the set's length
// notifies nobody, so a wait without a timeout would hang here.
#[test]
fn an_owner_that_unwinds_waits_for_the_item_on_another_thread() {
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .expect("test pool");
        let board = XBoard::default();
        let other_running = AtomicBool::new(false);
        let other_ended = AtomicBool::new(false);
        let owner_done = AtomicBool::new(false);
        let ended_at_unwind = AtomicBool::new(false);
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            pool.install(|| {
                rayon::broadcast(|ctx| {
                    if ctx.index() == 0 {
                        struct Done<'a>(&'a AtomicBool, &'a AtomicBool, &'a AtomicBool);
                        impl Drop for Done<'_> {
                            fn drop(&mut self) {
                                self.2
                                    .store(self.1.load(Ordering::SeqCst), Ordering::SeqCst);
                                self.0.store(true, Ordering::SeqCst);
                            }
                        }
                        let _done = Done(&owner_done, &other_ended, &ended_at_unwind);
                        board.run_set(
                            4,
                            &|_, _| {
                                if rayon::current_thread_index() != Some(0) {
                                    other_running.store(true, Ordering::SeqCst);
                                    std::thread::sleep(std::time::Duration::from_millis(50));
                                    other_ended.store(true, Ordering::SeqCst);
                                    return;
                                }
                                let start = std::time::Instant::now();
                                while !other_running.load(Ordering::SeqCst)
                                    && start.elapsed() < std::time::Duration::from_secs(5)
                                {
                                    std::hint::spin_loop();
                                }
                                panic!("the owner's own item panics");
                            },
                            false,
                        );
                    } else {
                        let mut taker = XTaker::default();
                        while !owner_done.load(Ordering::SeqCst) {
                            if !other_running.load(Ordering::SeqCst) {
                                let _ = board.take_one(&mut taker);
                            }
                            std::thread::yield_now();
                        }
                    }
                });
            });
        }));
        let _ = tx.send((
            result.is_err(),
            other_running.load(Ordering::SeqCst),
            ended_at_unwind.load(Ordering::SeqCst),
        ));
    });
    let (panicked, other_ran, ended_at_unwind) = rx
        .recv_timeout(std::time::Duration::from_secs(60))
        .expect("the owner returned");
    assert!(panicked, "the owner's panic reaches the caller");
    if other_ran {
        assert!(
            ended_at_unwind,
            "the owner waited for the other thread's item"
        );
    }
}

// Subject after subject puts a group set on the board inside a subject loop, whose threads take
// no pool job: one helper job waits in the queues, not one per set, and no set stays listed.
#[test]
fn group_sets_leave_one_helper_job_in_the_queues() {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(2)
        .build()
        .expect("test pool");
    let board = XBoard::default();
    let checked = AtomicBool::new(false);
    pool.install(|| {
        rayon::broadcast(|ctx| {
            if ctx.index() == 0 {
                for _ in 0..1_000 {
                    board.run_groups(6, &|_, _| {});
                }
                let jobs = Arc::strong_count(&board.groups) - 1;
                let sets = board.groups.sets.lock().unwrap().len();
                checked.store(true, Ordering::SeqCst);
                assert_eq!((jobs, sets), (1, 0));
            } else {
                while !checked.load(Ordering::SeqCst) {
                    std::thread::yield_now();
                }
            }
        });
    });
}

// A thread that helps stays, asleep, until a subject puts its units on the board, and takes one
// (the owner's first unit waits for it, up to a bound); a stopped loop wakes it without a ring.
#[test]
fn a_helper_sleeps_until_the_units_come() {
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let board = XBoard::default();
        let batch = board.batch(2);
        let helper_took = AtomicBool::new(false);
        let helper = std::sync::Mutex::new(None);
        std::thread::scope(|scope| {
            scope.spawn(|| {
                *helper.lock().unwrap() = Some(std::thread::current().id());
                board.help(&batch, &AtomicBool::new(false));
            });
            while helper.lock().unwrap().is_none() {
                std::thread::yield_now();
            }
            std::thread::sleep(std::time::Duration::from_millis(20));
            let helper_id = helper.lock().unwrap().expect("set above");
            board.run_units(&batch, 0, 8, &|_, _: &mut XTaker| {
                if std::thread::current().id() == helper_id {
                    helper_took.store(true, Ordering::SeqCst);
                    return;
                }
                let start = std::time::Instant::now();
                while !helper_took.load(Ordering::SeqCst)
                    && start.elapsed() < std::time::Duration::from_secs(20)
                {
                    std::thread::yield_now();
                }
            });
            batch.close_without_units(1);
        });
        let stopped_board = XBoard::default();
        let stopped_batch = stopped_board.batch(1);
        let stopped = AtomicBool::new(false);
        std::thread::scope(|scope| {
            scope.spawn(|| stopped_board.help(&stopped_batch, &stopped));
            std::thread::sleep(std::time::Duration::from_millis(20));
            stopped.store(true, Ordering::SeqCst);
        });
        let _ = tx.send(helper_took.load(Ordering::SeqCst));
    });
    let helper_took = rx
        .recv_timeout(std::time::Duration::from_secs(60))
        .expect("the helpers returned");
    assert!(helper_took, "the helper stayed for the units");
}
