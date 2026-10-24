//! Diagnostic (LOSAT_TIMING): where the time of one TBLASTX batch search goes, so that every
//! part of `[TIMING] search_total` has a number.
//!
//! This file is a child module of `blast_engine/run_impl.rs` (declared there with `#[path]`, as
//! `x_seed_release.rs` is).
//!
//! No NCBI counterpart: NCBI BLAST+ has no such timer. The module reads a monotonic clock and
//! adds differences to counters; it computes no value NCBI computes and no value the search
//! reads. With `LOSAT_TIMING` not set it never reads the clock.
//!
//! # What is measured
//!
//! The batch function (`run_impl.rs`) measured `search_total` (from the end of the set-up of the
//! search's parameters to the collected hit list) and, inside it, only single calls: the scan and
//! the extension of the reference frame loop, the re-evaluation of each HSP, the second linking
//! call, and the unit path's own units and merge (`[TIMING] x_pairpar`). On a 0.66 Mb by 0.61 Mb
//! pair the time between those was 1.10 s at one thread and 0.50 s at eight, with no number.
//!
//! Two lines are added, both per batch, after the lines that were there:
//!
//! - `[TIMING] x_search: setup .. subjects .. collect ..` cuts `search_total` into three wall
//!   times that add up to it: before the first subject (cutoffs, the unit path's preparation),
//!   the subject loop, and after it (the collector thread's join, the final e-value cut).
//! - `[TIMING] x_subject: <stage> .. ... sum .. (subjects=N)` cuts the time of the per-subject
//!   closure (`process_subject`) into the stages of `XStage`, in the order they run. One clock
//!   per subject (`XLap`) is read at every stage boundary, and the time since the previous
//!   reading goes to the stage that just ended, so the stages of a subject add up to its whole
//!   time (`sum`) with nothing left over. The clock is the closure's first local, so it is
//!   dropped last: the last stage, `end`, includes the release of the subject's translation and
//!   of whatever else the closure still owns. With several subjects the stage times are sums
//!   over the subjects, on whatever threads ran them, so `sum` can exceed `subjects` of the
//!   first line; with one subject the two are the same interval but for the loop's own set-up.
//!
//! `units` and `merge` are the calls whose insides `[TIMING] x_pairpar` already times, and
//! `link2` is the call `[TIMING] sum_stats_linking` times; `reevaluate` is the whole list pass,
//! of which `[TIMING] reevaluate` times the per-HSP calls. The old lines are unchanged.

use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};
use std::time::{Duration, Instant};

// No NCBI counterpart: names of the intervals of `process_subject` for the timer.
/// The stages of one subject in `process_subject`, in the order they run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum XStage {
    /// Entering the seed stage and clearing the worker's diagonal table.
    Enter,
    /// The six-frame translation of the subject (and of its compressed bases, if it has
    /// ambiguity codes).
    Translate,
    /// The subject's linking cutoffs and the unit path's view of the subject.
    Prepare,
    /// The unit path's units (`XPairPar::run_subject_units`).
    Units,
    /// The end of the unit path's stage (frees the group tables at the batch's last subject).
    UnitsRelease,
    /// The merge of the unit lists (`XPairPar::merge_subject`).
    Merge,
    /// The reference frame loop (scan, two-hit extension, per-frame append and sort).
    FrameLoop,
    /// The end of the lookup stage (frees the lookup table and the
    /// worker's diagonal table at the batch's last subject).
    SeedRelease,
    /// The first linking call, on every HSP the seed stage saved.
    Link1,
    /// The score-order check (and sort) and the e-value cut after the first linking call.
    Reap,
    /// The re-evaluation of the list and the score-order check (and sort) after it.
    Reevaluate,
    /// The second linking call, on the re-evaluated list.
    Link2,
    /// The conversion of the linked HSPs into output hits.
    Convert,
    /// Handing the hits to the collector.
    Emit,
    /// After the last stage boundary that was reached: the drop of what the closure owns (the
    /// translation, lists not consumed), or all of a subject that is too short to search.
    End,
}

const STAGE_COUNT: usize = 15;

/// The stages in the order of their discriminants, with the names the line prints.
const STAGES: [(XStage, &str); STAGE_COUNT] = [
    (XStage::Enter, "enter"),
    (XStage::Translate, "translate"),
    (XStage::Prepare, "prepare"),
    (XStage::Units, "units"),
    (XStage::UnitsRelease, "units_release"),
    (XStage::Merge, "merge"),
    (XStage::FrameLoop, "frame_loop"),
    (XStage::SeedRelease, "seed_release"),
    (XStage::Link1, "link1"),
    (XStage::Reap, "reap"),
    (XStage::Reevaluate, "reevaluate"),
    (XStage::Link2, "link2"),
    (XStage::Convert, "convert"),
    (XStage::Emit, "emit"),
    (XStage::End, "end"),
];

// No NCBI counterpart: counters of the timer; the search never reads them.
/// The stage times of one batch search. Shared by the threads that run its subjects.
pub(super) struct XSearchTimes {
    enabled: bool,
    stage_ns: [AtomicU64; STAGE_COUNT],
    subject_ns: AtomicU64,
    subjects: AtomicU64,
}

impl XSearchTimes {
    // No NCBI counterpart: `enabled` is whether LOSAT_TIMING is set (read by the batch function).
    /// Counters at zero. With `enabled` false no method reads the clock.
    pub(super) fn new(enabled: bool) -> Self {
        Self {
            enabled,
            stage_ns: std::array::from_fn(|_| AtomicU64::new(0)),
            subject_ns: AtomicU64::new(0),
            subjects: AtomicU64::new(0),
        }
    }

    // No NCBI counterpart: a clock reading for the first line's three intervals.
    /// The clock now, or `None` when the timer is off.
    pub(super) fn now(&self) -> Option<Instant> {
        self.enabled.then(Instant::now)
    }

    // No NCBI counterpart: starts the clock of one call of `process_subject`.
    /// The clock of one subject. It must be the first local of the closure, so that it is
    /// dropped after everything else the closure owns.
    pub(super) fn subject(&self) -> XLap<'_> {
        let start = self.now();
        XLap {
            times: self,
            start,
            last: start,
        }
    }

    fn stage(&self, stage: XStage) -> u64 {
        self.stage_ns[stage as usize].load(AtomicOrdering::Relaxed)
    }

    // No NCBI counterpart: the two lines of the module header.
    /// The two lines, as text. `search_start` is where `search_total` starts, `subjects_start`
    /// and `subjects_end` are `now()` before and after the subject loop, and `search_total` is
    /// the value the `[TIMING] search_total` line prints. `None` when the timer is off.
    pub(super) fn lines(
        &self,
        search_start: Instant,
        subjects_start: Option<Instant>,
        subjects_end: Option<Instant>,
        search_total: Duration,
    ) -> Option<[String; 2]> {
        let (subjects_start, subjects_end) = (subjects_start?, subjects_end?);
        let setup = subjects_start.duration_since(search_start);
        let subjects = subjects_end.duration_since(subjects_start);
        let collect = search_total.saturating_sub(subjects_end.duration_since(search_start));
        let search = format!(
            "[TIMING] x_search: setup {:.3}s subjects {:.3}s collect {:.3}s",
            setup.as_secs_f64(),
            subjects.as_secs_f64(),
            collect.as_secs_f64()
        );
        let mut subject = String::from("[TIMING] x_subject:");
        for (stage, name) in STAGES {
            subject.push_str(&format!(" {name} {:.3}s", self.stage(stage) as f64 / 1e9));
        }
        subject.push_str(&format!(
            " sum {:.3}s (subjects={})",
            self.subject_ns.load(AtomicOrdering::Relaxed) as f64 / 1e9,
            self.subjects.load(AtomicOrdering::Relaxed)
        ));
        Some([search, subject])
    }
}

// No NCBI counterpart: the running clock of the timer for one subject.
/// The clock of one subject: every `lap` gives the time since the previous one (or since the
/// start) to a stage, and the drop gives the rest to `XStage::End`.
pub(super) struct XLap<'a> {
    times: &'a XSearchTimes,
    start: Option<Instant>,
    last: Option<Instant>,
}

impl XLap<'_> {
    // No NCBI counterpart: one clock reading at a stage boundary.
    /// The stage `stage` ends here.
    #[inline]
    pub(super) fn lap(&mut self, stage: XStage) {
        if let Some(last) = self.last {
            let now = Instant::now();
            self.times.stage_ns[stage as usize].fetch_add(
                now.duration_since(last).as_nanos() as u64,
                AtomicOrdering::Relaxed,
            );
            self.last = Some(now);
        }
    }
}

impl Drop for XLap<'_> {
    // No NCBI counterpart: the last clock reading of a subject.
    fn drop(&mut self) {
        if let (Some(start), Some(last)) = (self.start, self.last) {
            let now = Instant::now();
            self.times.stage_ns[XStage::End as usize].fetch_add(
                now.duration_since(last).as_nanos() as u64,
                AtomicOrdering::Relaxed,
            );
            self.times.subject_ns.fetch_add(
                now.duration_since(start).as_nanos() as u64,
                AtomicOrdering::Relaxed,
            );
            self.times.subjects.fetch_add(1, AtomicOrdering::Relaxed);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn stage_sum(times: &XSearchTimes) -> u64 {
        STAGES.iter().map(|(stage, _)| times.stage(*stage)).sum()
    }

    #[test]
    fn the_stage_table_follows_the_discriminants() {
        for (index, (stage, _)) in STAGES.iter().enumerate() {
            assert_eq!(*stage as usize, index, "{stage:?}");
        }
        let mut names: Vec<&str> = STAGES.iter().map(|(_, name)| *name).collect();
        names.sort_unstable();
        names.dedup();
        assert_eq!(names.len(), STAGE_COUNT);
    }

    #[test]
    fn with_the_timer_off_nothing_is_read_or_counted() {
        let times = XSearchTimes::new(false);
        assert!(times.now().is_none());
        {
            let mut lap = times.subject();
            assert!(lap.start.is_none() && lap.last.is_none());
            lap.lap(XStage::Units);
            lap.lap(XStage::Link1);
        }
        assert_eq!(stage_sum(&times), 0);
        assert_eq!(times.subject_ns.load(AtomicOrdering::Relaxed), 0);
        assert_eq!(times.subjects.load(AtomicOrdering::Relaxed), 0);
        assert!(times
            .lines(Instant::now(), None, None, Duration::ZERO)
            .is_none());
    }

    #[test]
    fn the_stages_of_a_subject_add_up_to_its_time() {
        let times = XSearchTimes::new(true);
        for _ in 0..3 {
            let mut lap = times.subject();
            std::thread::sleep(Duration::from_millis(2));
            lap.lap(XStage::Translate);
            std::thread::sleep(Duration::from_millis(3));
            lap.lap(XStage::Link1);
            std::thread::sleep(Duration::from_millis(1));
        }
        let total = times.subject_ns.load(AtomicOrdering::Relaxed);
        assert_eq!(stage_sum(&times), total);
        assert_eq!(times.subjects.load(AtomicOrdering::Relaxed), 3);
        assert!(times.stage(XStage::Translate) >= 6_000_000);
        assert!(times.stage(XStage::Link1) >= 9_000_000);
        assert!(times.stage(XStage::End) >= 3_000_000);
        assert_eq!(times.stage(XStage::Units), 0);
    }

    #[test]
    fn a_subject_that_returns_early_is_all_in_the_last_stage() {
        let times = XSearchTimes::new(true);
        drop(times.subject());
        assert_eq!(times.subjects.load(AtomicOrdering::Relaxed), 1);
        assert_eq!(
            times.stage(XStage::End),
            times.subject_ns.load(AtomicOrdering::Relaxed)
        );
    }

    #[test]
    fn the_first_line_adds_up_to_search_total() {
        let times = XSearchTimes::new(true);
        let search_start = Instant::now();
        let subjects_start = search_start + Duration::from_millis(250);
        let subjects_end = subjects_start + Duration::from_millis(1500);
        let [search, subject] = times
            .lines(
                search_start,
                Some(subjects_start),
                Some(subjects_end),
                Duration::from_millis(2000),
            )
            .expect("timer on");
        assert_eq!(
            search,
            "[TIMING] x_search: setup 0.250s subjects 1.500s collect 0.250s"
        );
        assert!(subject.starts_with("[TIMING] x_subject: enter 0.000s translate 0.000s "));
        assert!(subject.ends_with(" end 0.000s sum 0.000s (subjects=0)"));
        assert_eq!(subject.matches(" 0.000s").count(), STAGE_COUNT + 1);
    }

    #[test]
    fn the_counters_are_shared_between_the_threads_of_a_search() {
        fn shared<S: Send + Sync>() {}
        shared::<XSearchTimes>();
    }
}
