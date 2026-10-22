//! EXPERIMENT (LOSAT_X_TBXNEXT): the subjects of a TBLASTX search handed to the threads of the
//! pool one at a time, each thread taking the next subject when it is free, with one worker state
//! per thread.
//!
//! This file is a child module of `blast_engine/run_impl.rs` (declared there with `#[path]`, as
//! `x_seed_release.rs` is).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1413
//! ```c
//!     itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
//!
//!     /* iterate over all subject sequences */
//!     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//!            != BLAST_SEQSRC_EOF) {
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:519-526
//! ```c++
//!     /* If internal iterator is uninitialized/invalid, retrieve the next chunk
//!        from the BlastSeqSrc */
//!     if (itr->current_pos == UINT4_MAX) {
//!         status = s_SeqDbGetNextChunk(seqdb_handle, itr);
//!         if (status == BLAST_SEQSRC_ERROR || status == BLAST_SEQSRC_EOF) {
//!             return status;
//!         }
//!     }
//! ```
//! Each NCBI search thread runs this loop with its own iterator over one shared database: a thread
//! that has searched its subjects asks for the next ones, so a long subject keeps one thread busy
//! while the others go on with the rest. Without the switch LOSAT gives the subject indices to
//! rayon as one indexed range (`(0..n).into_par_iter().for_each_init`): rayon cuts the range into
//! pieces before the search and a thread searches its piece in order, so the long subjects that
//! a genome lists first (the chromosomes) fall into one piece and one thread searches them one
//! after the other (3 transcripts against the fruit fly genome: 11.4 s at 8 threads, 11.6 s at
//! one). With the switch each thread of the pool takes the next index from one counter when it
//! is free, so the long subjects go to different threads. The threaded Wasm module
//! (`wasm-threads`) runs the same loop: its subject loop in `run_impl.rs` keeps each subject's
//! HSPs with the subject index, without the switch in rayon's pieces (`map_init`), with it in one
//! list the threads append to, and orders the lists by that index after the loop.
//!
//! What a thread does with a subject is unchanged: the same `process_subject` on a worker state
//! made by the same `init`. Each thread makes its state once, when it takes its first subject (a
//! thread that takes none makes none), where rayon made one per piece. The search's results do
//! not depend on which thread searched which subject or in which order: each subject's HSPs are
//! sent as one list (`WorkerState::tx`) and `report::final_hit_order` orders them by query and
//! subject index, as for rayon's pieces, whose order of completion is not fixed either. The
//! tables of the seed release (`x_seed_release.rs`) count the subjects whose stages have ended, not their order.
//!
//! The threads run their loops as one broadcast job each, not as jobs a thread may steal, so a
//! thread that waits inside a search (the unit path of PAIRPAR runs a subject's units on
//! the pool) runs only stolen units, which end on their own. A loop never waits: it takes indices
//! until the counter passes the last one. A panic in a search stops the other loops before their
//! next subject, and `rayon::broadcast` raises it on the calling thread.

use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::OnceLock;

// No NCBI counterpart: the variable is read once per process; it chooses how the subject loop is
// scheduled, not a value NCBI computes.
/// LOSAT_X_TBXNEXT: a strict-set switch (on unless it is `0`, or `LOSAT_X_DEFAULT=0` and it is
/// not set; `utils/xswitch.rs`).
pub(super) fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBXNEXT"))
}

// No NCBI counterpart for the counter itself (scheduling only): NCBI's iterators share the
// database's chunk state (module text); here one counter over the subject indices plays that part.
/// Runs `search` for every subject index `0..n`, once each, on the threads of the pool this thread
/// runs in: each thread takes the next index when it is free, with its own state from `init`, made
/// when it takes its first index. Call it from a thread of the search pool (`SearchPool::install`);
/// outside a pool, or in a pool of one thread, the indices are searched in order on this thread.
pub(super) fn for_each_subject<S>(
    n: usize,
    init: &(impl Fn() -> S + Sync),
    search: &(impl Fn(&mut S, usize) + Sync),
) {
    if n == 0 {
        return;
    }
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if rayon::current_thread_index().is_some() && rayon::current_num_threads() > 1 {
        let next = AtomicUsize::new(0);
        let stopped = AtomicBool::new(false);
        rayon::broadcast(|_| take_subjects(n, &next, &stopped, init, search));
        return;
    }
    let next = AtomicUsize::new(0);
    let stopped = AtomicBool::new(false);
    take_subjects(n, &next, &stopped, init, search);
}

/// The loop of one thread: the next index from `next` until it passes `n - 1` or another thread's
/// search panicked (`stopped`).
fn take_subjects<S>(
    n: usize,
    next: &AtomicUsize,
    stopped: &AtomicBool,
    init: &(impl Fn() -> S + Sync),
    search: &(impl Fn(&mut S, usize) + Sync),
) {
    struct StopOnPanic<'a>(&'a AtomicBool);
    impl Drop for StopOnPanic<'_> {
        fn drop(&mut self) {
            if std::thread::panicking() {
                self.0.store(true, Ordering::Relaxed);
            }
        }
    }
    let _stop = StopOnPanic(stopped);
    let mut state: Option<S> = None;
    while !stopped.load(Ordering::Relaxed) {
        let index = next.fetch_add(1, Ordering::Relaxed);
        if index >= n {
            break;
        }
        search(state.get_or_insert_with(init), index);
    }
}

#[cfg(all(
    test,
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
mod tests {
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
            for_each_subject(64, &|| (), &|_: &mut (), i| {
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
            });
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
                for_each_subject(100_000, &|| (), &|_: &mut (), i| {
                    searched.fetch_add(1, Ordering::Relaxed);
                    if i == 10 {
                        panic!("subject 10");
                    }
                    std::hint::black_box(work(i + 4));
                });
            })
        }));
        assert!(result.is_err());
        assert!(searched.into_inner() < 100_000);
    }

    #[test]
    fn outside_a_pool_the_subjects_run_in_order_on_this_thread() {
        let order = Mutex::new(Vec::new());
        for_each_subject(5, &|| (), &|_: &mut (), i| order.lock().unwrap().push(i));
        assert_eq!(order.into_inner().unwrap(), vec![0, 1, 2, 3, 4]);
    }
}
