//! Diagnostic (LOSAT_TIMING): the fields every program prints in one line, so that the time
//! outside the search has a number in every program and the programs can be compared.
//!
//! `main` takes the clock at its entry (`start`) and prints at its end (`print`):
//!
//! `[TIMING] common: since_start ..s outside_search ..s read_fasta ..s encode ..s translate ..s
//! seg ..s build_lookup ..s search ..s output ..s (search_sum ..s)`
//!
//! - `since_start` is the wall time from the entry of `main` to the line (the whole run, output
//!   written). It is printed apart from the programs' own `total` lines, which start later.
//! - `search` is the wall time during which at least one interval marked `Search` ran (the union
//!   of the intervals, so batches searched side by side count once), `search_sum` adds them up,
//!   and `outside_search` = `since_start` - `search`: the fixed costs of the run.
//! - The other fields add up the wall times of the intervals a program marks with `clock` and
//!   `add` (or `add_duration`), on the thread that runs them: reading the FASTA query and
//!   subject files (`read_fasta`), turning the letters into the search's codes or packed form
//!   (`encode`), translating (`translate`), the low-complexity filter (`seg`: SEG, or DUST for
//!   BLASTN), building the lookup tables (`build_lookup`), the search (`search`) and writing the
//!   report (`output`). A program marks the intervals its driver runs one after the other; work
//!   that runs inside the search on its threads (the per-subject translation of TBLASTX, for
//!   example) is part of `search` and of the program's own `[TIMING]` lines. A field a program
//!   does not mark stays 0. Intervals marked on several threads at once add up their wall times.
//! - `add_duration` takes an interval that ends at the call (an existing timer's `elapsed()`).
//!
//! No NCBI counterpart: a timer; it computes no value NCBI computes and none the search reads.
//! With `LOSAT_TIMING` not set it never reads the clock.

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Mutex, OnceLock};
use std::time::{Duration, Instant};

/// The intervals of the common line, in the order it prints them.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Phase {
    ReadFasta,
    Encode,
    Translate,
    Seg,
    BuildLookup,
    Search,
    Output,
}

const PHASES: [(Phase, &str); 7] = [
    (Phase::ReadFasta, "read_fasta"),
    (Phase::Encode, "encode"),
    (Phase::Translate, "translate"),
    (Phase::Seg, "seg"),
    (Phase::BuildLookup, "build_lookup"),
    (Phase::Search, "search"),
    (Phase::Output, "output"),
];

// No NCBI counterpart: the nanoseconds of each phase (index = `Phase as usize`).
static NS: [AtomicU64; 7] = [
    AtomicU64::new(0),
    AtomicU64::new(0),
    AtomicU64::new(0),
    AtomicU64::new(0),
    AtomicU64::new(0),
    AtomicU64::new(0),
    AtomicU64::new(0),
];

// No NCBI counterpart: the clock at the entry of `main` (only with LOSAT_TIMING).
static START: OnceLock<Instant> = OnceLock::new();

// No NCBI counterpart: the intervals marked `Search` (start, end), for their union.
static SEARCH_SPANS: Mutex<Vec<(Instant, Instant)>> = Mutex::new(Vec::new());

/// Whether `LOSAT_TIMING` is set (read once).
pub fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_TIMING").is_some())
}

/// Takes the clock of `since_start` (the entry of `main`; nothing without `LOSAT_TIMING`).
pub fn start() {
    if on() {
        let _ = START.set(Instant::now());
    }
}

/// PROBE (v4-A s17, LOSAT_X_S17STAMP with LOSAT_TIMING): prints `XSTAMP name us` (microseconds since
/// the entry of `main`).
pub fn x_stamp(name: &str) {
    static ON: OnceLock<bool> = OnceLock::new();
    if !*ON.get_or_init(|| std::env::var_os("LOSAT_X_S17STAMP").is_some()) {
        return;
    }
    if let Some(t) = START.get() {
        eprintln!("XSTAMP {name} {}", t.elapsed().as_micros());
    }
}

/// PROBE (v4-A s17): microseconds since the entry of `main` (0 without LOSAT_TIMING).
pub fn x_since_us() -> u64 {
    START.get().map_or(0, |t| t.elapsed().as_micros() as u64)
}

/// The clock at the start of an interval: `Some` only with `LOSAT_TIMING`.
#[inline]
pub fn clock() -> Option<Instant> {
    if on() {
        Some(Instant::now())
    } else {
        None
    }
}

/// Adds the time since `t0` (from `clock`) to `phase`; nothing when `t0` is `None`.
#[inline]
pub fn add(phase: Phase, t0: Option<Instant>) {
    if let Some(t0) = t0 {
        span(phase, t0, Instant::now());
    }
}

/// Adds an interval a program measured itself (an existing timer), which ends now, to `phase`.
#[inline]
pub fn add_duration(phase: Phase, d: Duration) {
    let end = Instant::now();
    span(phase, end.checked_sub(d).unwrap_or(end), end);
}

fn span(phase: Phase, start: Instant, end: Instant) {
    NS[phase as usize].fetch_add(
        end.saturating_duration_since(start).as_nanos() as u64,
        Ordering::Relaxed,
    );
    if phase == Phase::Search {
        SEARCH_SPANS
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .push((start, end));
    }
}

/// The total length of the union of the intervals.
fn union(spans: &mut [(Instant, Instant)]) -> Duration {
    spans.sort_unstable();
    let mut total = Duration::ZERO;
    let mut current: Option<(Instant, Instant)> = None;
    for &(start, end) in spans.iter() {
        match current {
            Some((s, e)) if start <= e => current = Some((s, e.max(end))),
            _ => {
                if let Some((s, e)) = current {
                    total += e.saturating_duration_since(s);
                }
                current = Some((start, end));
            }
        }
    }
    if let Some((s, e)) = current {
        total += e.saturating_duration_since(s);
    }
    total
}

/// The common line (`since` is the time since the entry of `main`, `None` when `start` was not
/// called; `search` the union of the `Search` intervals).
fn line(since: Option<Duration>, search: Duration) -> String {
    let since = since.unwrap_or(Duration::ZERO);
    let mut s = format!(
        "[TIMING] common: since_start {:.3}s outside_search {:.3}s",
        since.as_secs_f64(),
        since.saturating_sub(search).as_secs_f64()
    );
    for (phase, name) in PHASES {
        let secs = if phase == Phase::Search {
            search.as_secs_f64()
        } else {
            NS[phase as usize].load(Ordering::Relaxed) as f64 / 1e9
        };
        s.push_str(&format!(" {name} {secs:.3}s"));
    }
    s.push_str(&format!(
        " (search_sum {:.3}s)",
        NS[Phase::Search as usize].load(Ordering::Relaxed) as f64 / 1e9
    ));
    s
}

/// Prints the common line to stderr (only with `LOSAT_TIMING`); `main` calls it last.
pub fn print() {
    if !on() {
        return;
    }
    let since = START.get().map(|t| t.elapsed());
    let search = union(&mut SEARCH_SPANS.lock().unwrap_or_else(|e| e.into_inner()));
    eprintln!("{}", line(since, search));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_line_names_every_phase_in_order() {
        let l = line(
            Some(Duration::from_millis(1500)),
            Duration::from_millis(1000),
        );
        assert!(
            l.starts_with("[TIMING] common: since_start 1.500s outside_search 0.500s read_fasta ")
        );
        assert!(l.contains(" search 1.000s output "), "{l}");
        let names: Vec<&str> = l.split_whitespace().skip(2).step_by(2).collect::<Vec<_>>();
        assert_eq!(
            names,
            [
                "since_start",
                "outside_search",
                "read_fasta",
                "encode",
                "translate",
                "seg",
                "build_lookup",
                "search",
                "output",
                "(search_sum"
            ]
        );
    }

    // Batches searched side by side count once in `search`, and a gap between intervals is
    // outside it.
    #[test]
    fn the_search_is_the_union_of_its_intervals() {
        let t = Instant::now();
        let at = |ms: u64| t + Duration::from_millis(ms);
        let mut spans = vec![
            (at(10), at(40)),
            (at(0), at(20)),
            (at(30), at(35)),
            (at(50), at(60)),
            (at(60), at(70)),
        ];
        assert_eq!(union(&mut spans), Duration::from_millis(60));
        assert_eq!(union(&mut []), Duration::ZERO);
    }

    #[test]
    fn durations_add_up_per_phase() {
        let before = NS[Phase::Seg as usize].load(Ordering::Relaxed);
        add_duration(Phase::Seg, Duration::from_nanos(7));
        add_duration(Phase::Seg, Duration::from_nanos(5));
        assert!(NS[Phase::Seg as usize].load(Ordering::Relaxed) >= before + 12);
        add(Phase::Seg, None);
    }
}
