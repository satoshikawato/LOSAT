//! Diagnostic (LOSAT_TIMING): where the time of a TBLASTX linking call goes.
//!
//! No NCBI counterpart: NCBI BLAST+ has no such timer. The module reads a monotonic clock and
//! adds differences to counters; it computes no value NCBI computes and no value the search
//! reads. With `LOSAT_TIMING` not set it never reads the clock.
//!
//! # What is measured
//!
//! `[TIMING] x_subject` (`x_search_timing.rs`) gives each of a subject's two linking calls one
//! number, `link1` and `link2`. On large pairs those two are a third to two thirds of an
//! eight-thread search and gain much less from the threads than the seed stage does. Two lines
//! are added, both per batch, after the `x_subject` line, one per linking call:
//!
//! - `[TIMING] x_link1: reset .. sort .. groups .. kernels .. assemble .. final_sorts .. replay ..
//!   sum .. (calls=N hsps=N groups=N largest_group=N kernel_sum ..s kernel_max ..s)`
//! - `[TIMING] x_link2:` the same for the second call.
//!
//! The stages are the steps of `x_hsp_inplace::link_in_place`, the linking call of the default
//! path, in the order they run (`LinkStage`). One clock per call (`LinkLap`) is read at every
//! stage boundary on the thread that made the call, and the time since the previous reading goes
//! to the stage that just ended; `sum` is the time from the call's start to its last reading, so
//! the stages of a call add up to it exactly. What follows the last reading (the return and the
//! drop of the call's locals) is in `link1` / `link2` of `x_subject` only.
//! `kernels` is the wall time of the per-group linking, on the pool when the call may use it;
//! `kernel_sum` adds the times of the per-group calls on whatever threads ran them and
//! `kernel_max` is the longest of them, the least `kernels` can be with any number of threads.
//! `hsps` is the length of the list, `groups` the number of groups and `largest_group` the number
//! of HSPs in the largest one. With several subjects every value is a sum over the calls, except
//! `largest_group` and `kernel_max`, which are the largest over the calls (a sum of maxima over
//! the subjects is no bound on any one call).
//!
//! A linking call that takes the reference body instead (`LOSAT_DIAGNOSTICS`,
//! `LOSAT_DEBUG_CUTOFFS`) is not timed here: its line
//! shows `calls=0`, and `link1` / `link2` of `x_subject` still give its whole time.
//!
//! # How a call learns that it is timed
//!
//! The linking function has callers that are not a search's (unit tests), so its arguments stay as
//! they are. The batch function sets a flag of the calling thread immediately before each of its
//! two calls (`XLinkTimer::arm`), naming the call and the batch's counters; `LinkLap::start`, the
//! first statement of `link_in_place`, takes the flag; the batch function clears it after the call
//! (`XLinkTimer::disarm`) in case the call took another body. Nothing between `arm` and
//! `LinkLap::start` hands the thread to other work. The pool can give the thread another subject
//! while it waits for its own groups, but only after the flag was taken, and that subject's calls
//! arm, take and clear the flag before the thread comes back, so each call adds to the counters it
//! was armed with.

use std::cell::Cell;
use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};
use std::sync::Arc;
use std::time::Instant;

// No NCBI counterpart: names the two linking calls of a subject for the timer.
/// Which of a subject's two linking calls.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum LinkCall {
    /// On every HSP the seed stage saved (`XStage::Link1`).
    First,
    /// On the re-evaluated list (`XStage::Link2`).
    Second,
}

const CALLS: [(LinkCall, &str); 2] = [(LinkCall::First, "x_link1"), (LinkCall::Second, "x_link2")];

// No NCBI counterpart: names of the intervals of `link_in_place` for the timer.
/// The stages of one linking call in `x_hsp_inplace::link_in_place`, in the order they run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum LinkStage {
    /// The link ids and `num = 1` of every HSP, the per-context `ln K`.
    Reset,
    /// The sort of the list by reverse position (`sort_in_place`, `rev_compare_hsps_tbx`).
    Sort,
    /// The cut of the sorted list into its groups as sub-slices (`group_starts`, `group_slices`).
    Groups,
    /// The per-group linking, serial or on the pool, with the results back in group order.
    Kernels,
    /// The assembly: empty, the linked groups already lie in group order inside the one list.
    Assemble,
    /// The single final sort of the list's indices (`x_link_sort::final_sort`, the order of
    /// `rev_compare_hsps_transl` then `fwd_compare_hsps_transl`).
    FinalSorts,
    /// The walk of the chains in the sorted order and the move of the records into it.
    Replay,
}

const STAGE_COUNT: usize = 7;

/// The stages in the order of their discriminants, with the names the lines print.
const STAGES: [(LinkStage, &str); STAGE_COUNT] = [
    (LinkStage::Reset, "reset"),
    (LinkStage::Sort, "sort"),
    (LinkStage::Groups, "groups"),
    (LinkStage::Kernels, "kernels"),
    (LinkStage::Assemble, "assemble"),
    (LinkStage::FinalSorts, "final_sorts"),
    (LinkStage::Replay, "replay"),
];

// No NCBI counterpart: counters of the timer; the search never reads them.
/// The sums over the timed calls of one kind.
struct CallTimes {
    stage_ns: [AtomicU64; STAGE_COUNT],
    call_ns: AtomicU64,
    calls: AtomicU64,
    hsps: AtomicU64,
    groups: AtomicU64,
    largest_group: AtomicU64,
    kernel_sum_ns: AtomicU64,
    kernel_max_ns: AtomicU64,
}

impl CallTimes {
    fn new() -> Self {
        Self {
            stage_ns: std::array::from_fn(|_| AtomicU64::new(0)),
            call_ns: AtomicU64::new(0),
            calls: AtomicU64::new(0),
            hsps: AtomicU64::new(0),
            groups: AtomicU64::new(0),
            largest_group: AtomicU64::new(0),
            kernel_sum_ns: AtomicU64::new(0),
            kernel_max_ns: AtomicU64::new(0),
        }
    }

    fn line(&self, name: &str) -> String {
        let load = |counter: &AtomicU64| counter.load(AtomicOrdering::Relaxed);
        let seconds = |counter: &AtomicU64| load(counter) as f64 / 1e9;
        let mut line = format!("[TIMING] {name}:");
        for (stage, stage_name) in STAGES {
            line.push_str(&format!(
                " {stage_name} {:.3}s",
                seconds(&self.stage_ns[stage as usize])
            ));
        }
        line.push_str(&format!(
            " sum {:.3}s (calls={} hsps={} groups={} largest_group={} kernel_sum {:.3}s kernel_max {:.3}s)",
            seconds(&self.call_ns),
            load(&self.calls),
            load(&self.hsps),
            load(&self.groups),
            load(&self.largest_group),
            seconds(&self.kernel_sum_ns),
            seconds(&self.kernel_max_ns)
        ));
        line
    }
}

// No NCBI counterpart: counters of the timer; the search never reads them.
/// The linking times of one batch search. Shared by the threads that run its subjects.
pub(super) struct XLinkTimes {
    calls: [CallTimes; 2],
}

impl XLinkTimes {
    fn new() -> Self {
        Self {
            calls: [CallTimes::new(), CallTimes::new()],
        }
    }

    fn lines(&self) -> [String; 2] {
        CALLS.map(|(call, name)| self.calls[call as usize].line(name))
    }
}

// No NCBI counterpart: the flag of the module header ("How a call learns that it is timed").
thread_local! {
    static ARMED: Cell<Option<(LinkCall, Arc<XLinkTimes>)>> = const { Cell::new(None) };
}

// No NCBI counterpart: the batch function's side of the timer.
/// The linking timer of one batch search: its counters when `LOSAT_TIMING` is set, nothing
/// otherwise.
pub(super) struct XLinkTimer {
    times: Option<Arc<XLinkTimes>>,
}

impl XLinkTimer {
    // No NCBI counterpart: `enabled` is whether LOSAT_TIMING is set (read by the batch function).
    /// Counters at zero. With `enabled` false no method does anything.
    pub(super) fn new(enabled: bool) -> Self {
        Self {
            times: enabled.then(|| Arc::new(XLinkTimes::new())),
        }
    }

    // No NCBI counterpart: sets the calling thread's flag for the linking call that follows.
    /// The next linking call of this thread is `call` of this batch.
    #[inline]
    pub(super) fn arm(&self, call: LinkCall) {
        if let Some(times) = &self.times {
            ARMED.with(|armed| armed.set(Some((call, Arc::clone(times)))));
        }
    }

    // No NCBI counterpart: clears a flag that the linking call did not take.
    /// The linking call that was armed has returned.
    #[inline]
    pub(super) fn disarm(&self) {
        if self.times.is_some() {
            ARMED.with(|armed| armed.set(None));
        }
    }

    // No NCBI counterpart: the two lines of the module header.
    /// The two lines, as text. `None` when the timer is off.
    pub(super) fn lines(&self) -> Option<[String; 2]> {
        self.times.as_deref().map(XLinkTimes::lines)
    }
}

struct Running {
    call: LinkCall,
    times: Arc<XLinkTimes>,
    start: Instant,
    last: Instant,
}

// No NCBI counterpart: the running clock of the timer for one linking call.
/// The clock of one linking call: every `lap` gives the time since the previous one (or since
/// the start) to a stage, and the drop counts the call and its time up to the last `lap`.
pub(super) struct LinkLap {
    running: Option<Running>,
}

impl LinkLap {
    // No NCBI counterpart: takes the calling thread's flag and, if it was set, reads the clock.
    /// The clock of the linking call that starts here. It must be the first local of the
    /// function, so that it is dropped after everything else the function owns.
    #[inline]
    pub(super) fn start() -> Self {
        let running = ARMED.with(Cell::take).map(|(call, times)| {
            let start = Instant::now();
            Running {
                call,
                times,
                start,
                last: start,
            }
        });
        Self { running }
    }

    // No NCBI counterpart: a clock that is never read, for callers that are not a timed call (the
    // tests).
    /// A clock that measures nothing.
    #[cfg(test)]
    pub(super) fn off() -> Self {
        Self { running: None }
    }

    fn counters(&self) -> Option<&CallTimes> {
        self.running
            .as_ref()
            .map(|running| &running.times.calls[running.call as usize])
    }

    // No NCBI counterpart: one clock reading at a stage boundary.
    /// The stage `stage` ends here.
    #[inline]
    pub(super) fn lap(&mut self, stage: LinkStage) {
        if let Some(running) = &mut self.running {
            let now = Instant::now();
            running.times.calls[running.call as usize].stage_ns[stage as usize].fetch_add(
                now.duration_since(running.last).as_nanos() as u64,
                AtomicOrdering::Relaxed,
            );
            running.last = now;
        }
    }

    // No NCBI counterpart: the sizes the line prints next to the times.
    /// The call's list has `hsps` records in the groups of `group_sizes`.
    #[inline]
    pub(super) fn sizes(&self, hsps: usize, group_sizes: impl Iterator<Item = usize>) {
        if let Some(counters) = self.counters() {
            let (mut groups, mut largest) = (0u64, 0u64);
            for size in group_sizes {
                groups += 1;
                largest = largest.max(size as u64);
            }
            counters
                .hsps
                .fetch_add(hsps as u64, AtomicOrdering::Relaxed);
            counters.groups.fetch_add(groups, AtomicOrdering::Relaxed);
            counters
                .largest_group
                .fetch_max(largest, AtomicOrdering::Relaxed);
        }
    }

    // No NCBI counterpart: the clock of the per-group calls of this linking call.
    /// The timer of this call's per-group linking, for the threads that run the groups.
    #[inline]
    pub(super) fn kernels(&self) -> KernelTimes {
        KernelTimes {
            enabled: self.running.is_some(),
            sum_ns: AtomicU64::new(0),
            max_ns: AtomicU64::new(0),
        }
    }

    // No NCBI counterpart: adds the per-group times of this call to the batch's counters.
    /// The per-group linking of this call is over.
    #[inline]
    pub(super) fn add_kernels(&self, kernels: &KernelTimes) {
        if let Some(counters) = self.counters() {
            counters.kernel_sum_ns.fetch_add(
                kernels.sum_ns.load(AtomicOrdering::Relaxed),
                AtomicOrdering::Relaxed,
            );
            counters.kernel_max_ns.fetch_max(
                kernels.max_ns.load(AtomicOrdering::Relaxed),
                AtomicOrdering::Relaxed,
            );
        }
    }
}

impl Drop for LinkLap {
    // No NCBI counterpart: counts a timed linking call; the clock is not read again.
    fn drop(&mut self) {
        if let Some(running) = &self.running {
            let counters = &running.times.calls[running.call as usize];
            counters.call_ns.fetch_add(
                running.last.duration_since(running.start).as_nanos() as u64,
                AtomicOrdering::Relaxed,
            );
            counters.calls.fetch_add(1, AtomicOrdering::Relaxed);
        }
    }
}

// No NCBI counterpart: the clock of the per-group linking calls of one linking call.
/// The times of the per-group calls of one linking call, shared by the threads that run them.
pub(super) struct KernelTimes {
    enabled: bool,
    sum_ns: AtomicU64,
    max_ns: AtomicU64,
}

impl KernelTimes {
    // No NCBI counterpart: two clock readings around one per-group call.
    /// The result of `kernel`, with its time added when the linking call is timed.
    #[inline]
    pub(super) fn time<T>(&self, kernel: impl FnOnce() -> T) -> T {
        if !self.enabled {
            return kernel();
        }
        let start = Instant::now();
        let linked = kernel();
        let elapsed = start.elapsed().as_nanos() as u64;
        self.sum_ns.fetch_add(elapsed, AtomicOrdering::Relaxed);
        self.max_ns.fetch_max(elapsed, AtomicOrdering::Relaxed);
        linked
    }
}

#[cfg(test)]
#[path = "x_link_timing_tests.rs"]
mod tests;
