//! Graduated 2026-10-22 (was `LOSAT_X_TBNOIDNEXT`): the chunks of subjects of a TBLASTN search
//! handed to the threads of the pool one at a time, without a wait between entries, and their lists
//! taken back in OID order.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:307-344
//! ```c++
//!     CSeqDBLockHold locked(m_Atlas);
//! ...
//!    	m_Atlas.Lock(locked);
//!     if (! state_obj) {
//!         state_obj = & m_NextChunkOID;
//!     }
//! ...
//!     begin_chunk = * state_obj;
//!
//!     // fill the cache for all sequence in mmaped slice
//!     if (m_NumThreads) {
//!         SSeqResBuffer * buffer = m_CachedSeqs[cacheID];
//!         x_FillSeqBuffer(buffer, begin_chunk);
//!         end_chunk = begin_chunk + static_cast<int>(buffer->results.size());
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1413
//! ```c
//!     /* iterate over all subject sequences */
//!     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//!            != BLAST_SEQSRC_EOF) {
//! ```
//! Each NCBI search thread asks the database for the next chunk, under its lock, as soon as it has
//! searched the subjects of the last one; no thread waits for another. LOSAT's pool earlier
//! (`LOSAT_X_TBNOIDNEXT=0`, graduated 2026-10-22) took a few chunks per thread and waited until
//! all of them were searched before it handed over the next ones, so the threads with short
//! chunks waited for the longest one. Here the threads take the next chunk when they are free,
//! as NCBI's do.
//!
//! LOSAT keeps its OID-ordered reduction, the order of the one-thread search: the thread that
//! owns the search takes back the lists chunk by chunk in OID order and enters them into the
//! collector as before, and it searches a chunk itself while the next one in order is not ready.
//! The lists of chunks finished ahead of the oldest one not yet taken back wait; at most `window`
//! chunks are taken or waiting beyond it, and a thread waits before a new chunk when that window
//! is full (no NCBI counterpart: NCBI's threads write their lists into the shared HSP stream at
//! once; this bound only keeps the memory of the waiting lists bounded). The chunks and the
//! search of each subject are those of the chunks (`x_oid_chunk_bases`), so the reports are the
//! same.
//!
//! The threads start together: each thread's loop and the owner wait until every thread of the pool
//! has entered its loop before the first chunk is searched (no NCBI counterpart). A search may run
//! the units of its subject on the pool (a subject longer than a chunk's budget), and a thread that
//! waits for a unit taken by another thread runs other jobs of the pool meanwhile; were a thread's
//! loop still queued then, it would start there, on top of the chunk being searched, and could wait
//! for the window while that chunk is the one the owner waits for. Once every loop runs,
//! only units are left for a waiting thread to take.

use std::collections::BTreeMap;
use std::ops::Range;
use std::sync::{Condvar, Mutex, MutexGuard, PoisonError};

use anyhow::{bail, Result};

/// The chunks of one search, the window, and what the threads have done with them.
pub(super) struct XOidFlow<T> {
    chunks: Vec<Range<usize>>,
    window: usize,
    /// The threads that run `work`, besides the owner.
    workers: usize,
    state: Mutex<XFlowState<T>>,
    changed: Condvar,
}

struct XFlowState<T> {
    /// The next chunk a thread takes.
    next: usize,
    /// The chunks taken back by the owner of the search (all of those before this index).
    taken_back: usize,
    /// The lists of finished chunks the owner has not taken back yet.
    finished: BTreeMap<usize, Vec<T>>,
    /// Set when the owner stops taking lists back (the end, an error, a panic).
    closed: bool,
    /// Set when a thread stopped in the middle of a chunk (a panic).
    failed: bool,
    /// The most finished lists that waited at once (for the tests and the diagnostic line).
    peak_waiting: usize,
    /// The threads that have entered `work`.
    started: usize,
}

impl<T: Send> XOidFlow<T> {
    /// A flow over `chunks` (consecutive ranges of subject indices, in OID order) with at most
    /// `window` chunks (at least one) taken or waiting beyond the oldest one not taken back, for
    /// the owner and `workers` threads that run `work`.
    pub(super) fn new(chunks: Vec<Range<usize>>, window: usize, workers: usize) -> Self {
        XOidFlow {
            chunks,
            window: window.max(1),
            workers,
            state: Mutex::new(XFlowState {
                next: 0,
                taken_back: 0,
                finished: BTreeMap::new(),
                closed: false,
                failed: false,
                peak_waiting: 0,
                started: 0,
            }),
            changed: Condvar::new(),
        }
    }

    fn lock(&self) -> MutexGuard<'_, XFlowState<T>> {
        self.state.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// Waits until every worker has entered `work` (or the flow is closed or failed), so that no
    /// loop of the flow is still queued in the pool when a search starts (see the module text).
    fn start_together<'s>(
        &self,
        mut state: MutexGuard<'s, XFlowState<T>>,
    ) -> MutexGuard<'s, XFlowState<T>> {
        while state.started < self.workers && !state.closed && !state.failed {
            state = self
                .changed
                .wait(state)
                .unwrap_or_else(PoisonError::into_inner);
        }
        state
    }

    /// The next chunk, if one is left and it lies within the window; marks it taken.
    fn take(&self, state: &mut XFlowState<T>) -> Option<usize> {
        let index = state.next;
        if index < self.chunks.len() && index < state.taken_back.saturating_add(self.window) {
            state.next += 1;
            Some(index)
        } else {
            None
        }
    }

    /// Searches chunk `index` (each subject in order, outside the lock) and leaves its lists for
    /// the owner. A panic of `search` marks the flow failed, so that the owner stops waiting.
    fn search_chunk(&self, index: usize, search: &(impl Fn(usize) -> T + Sync)) {
        struct Failed<'a, T>(&'a XOidFlow<T>);
        impl<T> Drop for Failed<'_, T> {
            fn drop(&mut self) {
                if std::thread::panicking() {
                    let mut state = self.0.state.lock().unwrap_or_else(PoisonError::into_inner);
                    state.failed = true;
                    drop(state);
                    self.0.changed.notify_all();
                }
            }
        }
        let _failed = Failed(self);
        let lists: Vec<T> = self.chunks[index].clone().map(search).collect();
        let mut state = self.lock();
        state.finished.insert(index, lists);
        state.peak_waiting = state.peak_waiting.max(state.finished.len());
        drop(state);
        self.changed.notify_all();
    }

    /// The loop of a thread of the pool: take the next chunk when there is one within the window,
    /// search it, and wait while the window is full, until no chunk is left or the owner closes.
    pub(super) fn work(&self, search: &(impl Fn(usize) -> T + Sync)) {
        let mut state = self.lock();
        state.started += 1;
        self.changed.notify_all();
        state = self.start_together(state);
        loop {
            if state.closed || state.failed || state.next >= self.chunks.len() {
                return;
            }
            if let Some(index) = self.take(&mut state) {
                drop(state);
                self.search_chunk(index, search);
                state = self.lock();
            } else {
                state = self
                    .changed
                    .wait(state)
                    .unwrap_or_else(PoisonError::into_inner);
            }
        }
    }

    /// For the owner of the search: the lists of the next chunk in OID order, or `None` after the
    /// last one. While that chunk is not finished the owner searches the next chunk within the
    /// window itself, and waits only when every chunk it could take is taken.
    pub(super) fn next_in_order(
        &self,
        search: &(impl Fn(usize) -> T + Sync),
    ) -> Result<Option<Vec<T>>> {
        let mut state = self.start_together(self.lock());
        loop {
            let wanted = state.taken_back;
            if wanted >= self.chunks.len() {
                return Ok(None);
            }
            if let Some(lists) = state.finished.remove(&wanted) {
                state.taken_back += 1;
                drop(state);
                self.changed.notify_all();
                return Ok(Some(lists));
            }
            if state.failed {
                bail!(
                    "a thread of the TBLASTN subject search stopped before its chunk was searched"
                );
            }
            if let Some(index) = self.take(&mut state) {
                drop(state);
                self.search_chunk(index, search);
                state = self.lock();
            } else {
                state = self
                    .changed
                    .wait(state)
                    .unwrap_or_else(PoisonError::into_inner);
            }
        }
    }

    /// Stops the threads: no chunk is taken after this (the owner's end, error or panic).
    pub(super) fn close(&self) {
        self.lock().closed = true;
        self.changed.notify_all();
    }

    /// The number of chunks.
    pub(super) fn chunks(&self) -> usize {
        self.chunks.len()
    }

    /// The most finished lists that waited for the owner at once.
    pub(super) fn peak_waiting(&self) -> usize {
        self.lock().peak_waiting
    }
}

/// Closes the flow when dropped, so that the threads of the pool stop whether the owner of the
/// search returns, fails or panics.
pub(super) struct XCloseOnDrop<'a, T: Send>(pub(super) &'a XOidFlow<T>);

impl<T: Send> Drop for XCloseOnDrop<'_, T> {
    fn drop(&mut self) {
        self.0.close();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    /// Runs a flow as `run_preliminary_subjects` does: the owner takes back the lists in order
    /// inside a scope of the pool whose other threads run `work`. Returns the subjects in the
    /// order the owner received them and the peak of waiting lists.
    fn run_flow(
        threads: usize,
        chunks: Vec<Range<usize>>,
        window: usize,
        search: &(impl Fn(usize) -> usize + Sync),
    ) -> (Vec<usize>, usize) {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let flow = XOidFlow::new(chunks, window, threads - 1);
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
        for (threads, window) in [(1, 1), (2, 1), (2, 2), (3, 4), (4, 2), (8, 4), (8, 32)] {
            let (received, peak) = run_flow(threads, chunks.clone(), window, &search);
            assert_eq!(
                received,
                (0..subjects).collect::<Vec<_>>(),
                "{threads} threads, window {window}"
            );
            assert!(
                peak <= window,
                "{threads} threads, window {window}: peak {peak}"
            );
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
        for (threads, window) in [(2, 1), (3, 1), (4, 1), (4, 2), (8, 1), (8, 4)] {
            let (received, peak) = run_flow(threads, chunks.clone(), window, &search);
            assert_eq!(
                received,
                (0..24).collect::<Vec<_>>(),
                "{threads} threads, window {window}"
            );
            assert!(
                peak <= window,
                "{threads} threads, window {window}: peak {peak}"
            );
        }
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
}
