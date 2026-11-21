//! EXPERIMENT (LOSAT_X_TBXNEXT): the subjects of a TBLASTX search handed to the threads of the
//! pool one at a time, each thread taking the next subject when it is free, with one worker state
//! per thread.
//!
//! This file is a child module of `blast_engine/run_impl.rs` (declared there with `#[path]`, as
//! `x_seed_release.rs` is).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1413
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:519-526
//! Each NCBI search thread runs this loop with its own iterator over one shared database: a thread
//! that has searched its subjects asks for the next ones, so a long subject keeps one thread busy
//! while the others go on with the rest. Without the switch LOSAT gives the subject indices to
//! rayon as one indexed range (`(0..n).into_par_iter().for_each_init`): rayon cuts the range into
//! pieces before the search and a thread searches its piece in order, so the long subjects that
//! a genome lists first (the chromosomes) fall into one piece and one thread searches them one
//! after the other (3 transcripts against the fruit fly genome: 11.4 s at 8 threads, 11.6 s at
//! one). With the switch each thread of the pool takes the next index from one counter when it
//! is free, so the long subjects go to different threads. The threaded Wasm module
//! (`wasm-threads`) runs the same loop, with the collector the subject loop keeps in place of
//! the writer thread (`hit_collector.rs`).
//!
//! What a thread does with a subject is unchanged: the same `process_subject` on a worker state
//! made by the same `init`. Each thread makes its state once, when it takes its first subject (a
//! thread that takes none makes none), where rayon made one per piece. The search's results do
//! not depend on which thread searched which subject or in which order: each subject's HSPs are
//! sent as one message with its index (`WorkerState::tx`) and the batch's collector takes them in
//! subject order (`hit_collector.rs`), as for rayon's pieces, whose order of completion is not
//! fixed either. The tables of the seed release (`x_seed_release.rs`) count the subjects whose
//! stages have ended, not their order.
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
/// LOSAT_X_TBXNEXT: a strict-set switch (on unless it is `0`; `utils/xswitch.rs`).
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
    // EXPERIMENT (LOSAT_X_TBXBOARD): what a thread does once no subject is left to take (it helps
    // with the batch's units, tblastx/x_work_board.rs `XBoard::help`), given the loop's stop flag.
    help: Option<&(dyn Fn(&AtomicBool) + Sync)>,
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
        rayon::broadcast(|_| {
            take_subjects(n, &next, &stopped, init, search);
            if let Some(help) = help {
                if !stopped.load(Ordering::Relaxed) {
                    help(&stopped);
                }
            }
        });
        return;
    }
    let _ = help;
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
#[path = "x_subject_next_tests.rs"]
mod tests;
