use super::for_each_subject;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

fn pool(threads: usize) -> rayon::ThreadPool {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("a test pool")
}

// A few slow subjects first (a genome's chromosomes), then many fast ones.
fn work(i: usize) -> u64 {
    let rounds = if i < 4 {
        2_000_000
    } else {
        200 + (i * 7919) % 500
    };
    (0..rounds as u64).fold(i as u64, |acc, x| acc.wrapping_mul(31).wrapping_add(x))
}

#[test]
fn searches_every_subject_once_with_one_state_per_thread() {
    for (threads, n) in [(1, 10), (4, 0), (4, 1), (4, 3), (4, 500), (8, 1000), (3, 7)] {
        let inits = AtomicUsize::new(0);
        let seen = Mutex::new(vec![0usize; n]);
        let threads_used = Mutex::new(std::collections::BTreeSet::new());
        pool(threads).install(|| {
            for_each_subject(
                n,
                &|| {
                    inits.fetch_add(1, Ordering::Relaxed);
                    Vec::<usize>::new()
                },
                &|mine: &mut Vec<usize>, i| {
                    std::hint::black_box(work(i));
                    mine.push(i);
                    seen.lock().unwrap()[i] += 1;
                    threads_used
                        .lock()
                        .unwrap()
                        .insert(rayon::current_thread_index());
                },
                None,
            );
        });
        let seen = seen.into_inner().unwrap();
        assert!(
            seen.iter().all(|&c| c == 1),
            "threads {threads}, n {n}: {seen:?}"
        );
        let inits = inits.into_inner();
        assert!(
            inits <= threads.min(n),
            "threads {threads}, n {n}: {inits} states"
        );
        assert_eq!(inits, threads_used.into_inner().unwrap().len());
    }
}

#[test]
fn the_first_subjects_are_searched_at_once_on_different_threads() {
    // Each of the first four subjects waits until all four have started: only threads that
    // take the next subject when they are free get there (one thread searching the four in
    // order would wait for ever; the wait gives up after a minute and the test fails).
    let threads = 4;
    let started = Mutex::new(0usize);
    let all_started = std::sync::Condvar::new();
    let owner = Mutex::new(vec![None; 64]);
    pool(threads).install(|| {
        for_each_subject(
            64,
            &|| (),
            &|_: &mut (), i| {
                if i < threads {
                    let mut count = started.lock().unwrap();
                    *count += 1;
                    all_started.notify_all();
                    let (count, timeout) = all_started
                        .wait_timeout_while(count, std::time::Duration::from_secs(60), |c| {
                            *c < threads
                        })
                        .unwrap();
                    assert!(
                        !timeout.timed_out(),
                        "{} of the first subjects started",
                        *count
                    );
                }
                owner.lock().unwrap()[i] = rayon::current_thread_index();
            },
            None,
        );
    });
    let owner = owner.into_inner().unwrap();
    let first: std::collections::BTreeSet<_> = owner[..threads].iter().collect();
    assert_eq!(first.len(), threads, "the first subjects ran on {owner:?}");
}

#[test]
fn a_panic_reaches_the_caller_and_stops_the_other_threads() {
    let searched = AtomicUsize::new(0);
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        pool(4).install(|| {
            for_each_subject(
                100_000,
                &|| (),
                &|_: &mut (), i| {
                    searched.fetch_add(1, Ordering::Relaxed);
                    if i == 10 {
                        panic!("subject 10");
                    }
                    // About 50 µs a subject: the other threads take subjects until the panicking
                    // one has printed its message and unwound to `take_subjects`, and with
                    // subjects of about 1 µs (`work(i + 4)`) they took all of them on a Windows
                    // runner at opt-level 3 (2026-10-23). All of them now take seconds.
                    std::hint::black_box(
                        (0..50_000u64)
                            .fold(i as u64, |acc, x| acc.wrapping_mul(31).wrapping_add(x)),
                    );
                },
                None,
            );
        })
    }));
    assert!(result.is_err());
    assert!(searched.into_inner() < 100_000);
}

#[test]
fn outside_a_pool_the_subjects_run_in_order_on_this_thread() {
    let order = Mutex::new(Vec::new());
    for_each_subject(
        5,
        &|| (),
        &|_: &mut (), i| order.lock().unwrap().push(i),
        None,
    );
    assert_eq!(order.into_inner().unwrap(), vec![0, 1, 2, 3, 4]);
}
