//! Graduated 2026-10-22 (was `LOSAT_X_TBNOIDNEXT`): the chunks of subjects of a TBLASTN search
//! handed to the threads of the pool one at a time, without a wait between entries, and their lists
//! taken back in OID order.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:307-344
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1413
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
//! only units and the linking of a subject's query contexts (`stage_d_linking/x_link_ctx.rs`) are
//! left for a waiting thread to take.
//!
//! `LOSAT_X_TBNFLOWHELP` (on unless `=0`): a thread that waits for the flow
//! after the start (the owner for the next chunk in order, a thread of the pool for the window)
//! runs the pool's jobs meanwhile instead of sleeping on the flow's condition (no NCBI
//! counterpart: NCBI's threads do not wait for one another). Those jobs are units of the subjects
//! being searched and the linking of their query contexts, which do not wait for the flow; a
//! unit's result goes back to its subject by unit number and a context's linked values by input
//! position, so the reports are the same. With `=0` the owner sleeps while a long subject's
//! units wait in the pool, when it took a short subject first and every chunk it could take is
//! taken.

use std::collections::BTreeMap;
use std::ops::Range;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Condvar, Mutex, MutexGuard, OnceLock, PoisonError};
use std::time::{Duration, Instant};

use anyhow::{bail, Result};

/// The chunks of one search, the window, and what the threads have done with them.
pub(super) struct XOidFlow<T> {
    chunks: Vec<Range<usize>>,
    window: usize,
    /// The threads that run `work`, besides the owner.
    workers: usize,
    /// `LOSAT_X_TBNFLOWHELP`: when set, a waiting thread runs the pool's jobs, and waits at most
    /// this long for a change when there is none (`wait`).
    help: Option<Duration>,
    /// Diagnostic (`LOSAT_X_STATS`): the time the owner and the other threads waited, and the time
    /// and number of the pool's jobs they ran while waiting (`stats_line`).
    stats: Option<XWaitStats>,
    state: Mutex<XFlowState<T>>,
    changed: Condvar,
}

/// Diagnostic (`LOSAT_X_STATS`): nanoseconds and counts of the waits of the flow; index 0 is the
/// owner, 1 the other threads.
#[derive(Default)]
struct XWaitStats {
    idle_ns: [AtomicU64; 2],
    help_ns: [AtomicU64; 2],
    jobs: [AtomicU64; 2],
}

// No NCBI counterpart (memory only; no value changes): reads the LOSAT_X_TBNLOOPFREE variable once.
/// LOSAT_X_TBNLOOPFREE, in the strict set (on unless `=0`): whether a thread frees the buffers it
/// keeps between subjects when its loop over the chunks ends (`x_loop_free`).
pub(super) fn x_loop_free_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNLOOPFREE"))
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1412-1413,1604-1605,1693-1694
/// NCBI's search thread frees its word-finder state (diagonal table, initial hits, offset pairs)
/// when its loop over the subjects ends, and its caller then the gapped alignment memory; a
/// subject's translation is freed after the subject (blast_engine.c:627-631). Here a thread's loop
/// over the chunks of the flow (or the owner's reduction) ends: the buffers the thread kept from
/// one subject to the next go then (`search_gapped::XUnitBuffers`, the frames kept for the next
/// plan of `x_subject_translate`, the letter buffer of `x_subject_store`). Kept, they stayed in
/// the thread's heap through the query batch's redo or traceback, under its blocks, until the next
/// batch's loop. A thread that searches a subject after it makes new ones, as for its first
/// subject (a unit on the pool always starts so): a new diagonal table's entries all lie beyond
/// the two-hit window, as a kept table's entries do after NCBI's offset step
/// (blast_extend.c:166-171), so no reported value changes. The presence bytes of `x_scan_byte`
/// (32 KiB a table, which every thread reads) stay: the chunks of a split query batch are searched
/// in flows of their own with the same table, and building them again in each flow cost more than
/// they hold (81 KB of subjects at 32 threads: 1.2 % more instructions and 1,600 more waits).
/// `LOSAT_X_TBNLOOPFREE=0` keeps them all.
pub(super) fn x_loop_free() {
    if !x_loop_free_on() {
        return;
    }
    super::search_gapped::x_free_unit_buffers();
    super::x_subject_translate::x_free_kept_frames();
    super::x_subject_store::x_free_letters();
}

/// `LOSAT_X_TBNFLOWHELP` (strict set; `=0` off): how long a waiting thread that found no job
/// waits for a change of the flow before it looks for a job again, 100 us unless the variable
/// gives a number of microseconds from 2 to 1,000,000.
pub(super) fn x_flow_help() -> Option<Duration> {
    static HELP: OnceLock<Option<Duration>> = OnceLock::new();
    *HELP.get_or_init(|| {
        crate::utils::xswitch::value("LOSAT_X_TBNFLOWHELP", 2..=1_000_000, 100)
            .map(|us| Duration::from_micros(us as u64))
    })
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
    /// Counts the changes above that a waiting thread may wait for (each one notifies), so that
    /// a thread that let go of the lock to look for the pool's jobs sees one it missed.
    changes: u64,
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
            help: x_flow_help(),
            stats: std::env::var_os("LOSAT_X_STATS").map(|_| XWaitStats::default()),
            state: Mutex::new(XFlowState {
                next: 0,
                taken_back: 0,
                finished: BTreeMap::new(),
                closed: false,
                failed: false,
                peak_waiting: 0,
                started: 0,
                changes: 0,
            }),
            changed: Condvar::new(),
        }
    }

    /// The same flow with `help` in place of the switch's (for the tests).
    #[cfg(test)]
    fn x_with_help(mut self, help: Option<Duration>) -> Self {
        self.help = help;
        self
    }

    fn lock(&self) -> MutexGuard<'_, XFlowState<T>> {
        self.state.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// Waits for a change of the flow after the start (`owner`: the owner of the search, else a
    /// thread in `work`; `idle`: the waits of the caller since it last ran a job or took a chunk).
    /// With `LOSAT_X_TBNFLOWHELP` a thread of the pool runs one pending job of the pool and
    /// returns, so that its caller looks at the flow again. When there is none it returns at once
    /// if the flow changed meanwhile, and otherwise waits for a change for at most `help` times
    /// 2^`idle` (up to 16 times `help`), after which it looks for a job again. Every loop of the
    /// flow runs by then (`start_together`), so the job is a unit of a subject being searched or
    /// the linking of one of its query contexts, and never a loop.
    fn wait<'s>(
        &'s self,
        state: MutexGuard<'s, XFlowState<T>>,
        owner: bool,
        idle: &mut u32,
    ) -> MutexGuard<'s, XFlowState<T>> {
        let start = self.stats.as_ref().map(|_| Instant::now());
        let who = usize::from(!owner);
        let (state, helped) = match self.help {
            Some(help) if rayon::current_thread_index().is_some() => {
                let seen = state.changes;
                drop(state);
                if let Some(rayon::Yield::Executed) = rayon::yield_now() {
                    *idle = 0;
                    (self.lock(), true)
                } else {
                    let state = self.lock();
                    if state.changes != seen {
                        (state, false)
                    } else {
                        let timeout = help * (1 << (*idle).min(4));
                        *idle += 1;
                        let state = self
                            .changed
                            .wait_timeout(state, timeout)
                            .unwrap_or_else(PoisonError::into_inner)
                            .0;
                        (state, false)
                    }
                }
            }
            _ => (
                self.changed
                    .wait(state)
                    .unwrap_or_else(PoisonError::into_inner),
                false,
            ),
        };
        if let (Some(stats), Some(start)) = (self.stats.as_ref(), start) {
            let ns = u64::try_from(start.elapsed().as_nanos()).unwrap_or(u64::MAX);
            if helped {
                stats.help_ns[who].fetch_add(ns, Ordering::Relaxed);
                stats.jobs[who].fetch_add(1, Ordering::Relaxed);
            } else {
                stats.idle_ns[who].fetch_add(ns, Ordering::Relaxed);
            }
        }
        state
    }

    /// Diagnostic (`LOSAT_X_STATS`): the waits of the flow, in milliseconds, for stderr (one line
    /// per flow, that is per query batch searched by the flow).
    pub(super) fn stats_line(&self) -> Option<String> {
        let stats = self.stats.as_ref()?;
        let ms = |v: &AtomicU64| v.load(Ordering::Relaxed) as f64 / 1e6;
        Some(format!(
            "[X_STATS] TBNFLOWHELP help={} owner_idle_ms={:.1} owner_help_ms={:.1} owner_jobs={} workers_idle_ms={:.1} workers_help_ms={:.1} workers_jobs={}",
            self.help.map_or(0, |d| d.as_micros()),
            ms(&stats.idle_ns[0]),
            ms(&stats.help_ns[0]),
            stats.jobs[0].load(Ordering::Relaxed),
            ms(&stats.idle_ns[1]),
            ms(&stats.help_ns[1]),
            stats.jobs[1].load(Ordering::Relaxed),
        ))
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
                    state.changes += 1;
                    drop(state);
                    self.0.changed.notify_all();
                }
            }
        }
        let _failed = Failed(self);
        let plan = crate::algorithm::common::x_db_stream::x_plan_chunk();
        let lists: Vec<T> = self.chunks[index].clone().map(search).collect();
        drop(plan);
        let mut state = self.lock();
        state.finished.insert(index, lists);
        state.peak_waiting = state.peak_waiting.max(state.finished.len());
        state.changes += 1;
        drop(state);
        self.changed.notify_all();
    }

    /// The loop of a thread of the pool: take the next chunk when there is one within the window,
    /// search it, and wait while the window is full, until no chunk is left or the owner closes.
    pub(super) fn work(&self, search: &(impl Fn(usize) -> T + Sync)) {
        let mut state = self.lock();
        state.started += 1;
        state.changes += 1;
        self.changed.notify_all();
        state = self.start_together(state);
        let mut idle = 0;
        loop {
            if state.closed || state.failed || state.next >= self.chunks.len() {
                drop(state);
                x_loop_free();
                return;
            }
            if let Some(index) = self.take(&mut state) {
                drop(state);
                self.search_chunk(index, search);
                idle = 0;
                state = self.lock();
            } else {
                state = self.wait(state, false, &mut idle);
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
        let mut idle = 0;
        loop {
            let wanted = state.taken_back;
            if wanted >= self.chunks.len() {
                return Ok(None);
            }
            if let Some(lists) = state.finished.remove(&wanted) {
                state.taken_back += 1;
                state.changes += 1;
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
                idle = 0;
                state = self.lock();
            } else {
                state = self.wait(state, true, &mut idle);
            }
        }
    }

    /// Stops the threads: no chunk is taken after this (the owner's end, error or panic).
    pub(super) fn close(&self) {
        let mut state = self.lock();
        state.closed = true;
        state.changes += 1;
        drop(state);
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
#[path = "x_oid_flow_tests.rs"]
mod tests;
