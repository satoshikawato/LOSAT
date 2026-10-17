//! Diagnostic (LOSAT_TIMING): the fields every program prints in one line, so that the time
//! outside the search has a number in every program and the programs can be compared.
//!
//! `main` takes the clock at its entry (`start`) and prints at its end (`print`):
//!
//! `[TIMING] common: since_start ..s read_fasta ..s encode ..s translate ..s seg ..s
//! build_lookup ..s search ..s output ..s`
//!
//! - `since_start` is the wall time from the entry of `main` to the line (the whole run, output
//!   written). It is printed apart from the programs' own `total` lines, which start later.
//! - The other fields add up the wall times of the intervals a program marks with `clock` and
//!   `add` (or `add_duration`), on the thread that runs them: reading the FASTA query and
//!   subject files (`read_fasta`), turning the letters into the search's codes or packed form
//!   (`encode`), translating (`translate`), the low-complexity filter (`seg`: SEG, or DUST for
//!   BLASTN), building the lookup tables (`build_lookup`), the search (`search`) and writing the
//!   report (`output`). A program marks the intervals its driver runs one after the other; work
//!   that runs inside the search on its threads (the per-subject translation of TBLASTX, for
//!   example) is part of `search` and of the program's own `[TIMING]` lines. A field a program
//!   does not mark stays 0. Intervals marked on several threads at once add up their wall times.
//!
//! No NCBI counterpart: a timer; it computes no value NCBI computes and none the search reads.
//! With `LOSAT_TIMING` not set it never reads the clock.

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;
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
        add_duration(phase, t0.elapsed());
    }
}

/// Adds a duration a program measured itself (an existing timer) to `phase`.
#[inline]
pub fn add_duration(phase: Phase, d: Duration) {
    NS[phase as usize].fetch_add(d.as_nanos() as u64, Ordering::Relaxed);
}

/// The common line, or `None` without `LOSAT_TIMING` (`since` is the time since the entry of
/// `main`, `None` when `start` was not called).
fn line(since: Option<Duration>) -> String {
    let mut s = format!(
        "[TIMING] common: since_start {:.3}s",
        since.map_or(0.0, |d| d.as_secs_f64())
    );
    for (phase, name) in PHASES {
        let ns = NS[phase as usize].load(Ordering::Relaxed);
        s.push_str(&format!(" {name} {:.3}s", ns as f64 / 1e9));
    }
    s
}

/// Prints the common line to stderr (only with `LOSAT_TIMING`); `main` calls it last.
pub fn print() {
    if !on() {
        return;
    }
    eprintln!("{}", line(START.get().map(|t| t.elapsed())));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_line_names_every_phase_in_order() {
        let l = line(Some(Duration::from_millis(1500)));
        assert!(l.starts_with("[TIMING] common: since_start 1.500s read_fasta "));
        let names: Vec<&str> = l.split_whitespace().skip(2).step_by(2).collect::<Vec<_>>();
        assert_eq!(
            names,
            [
                "since_start",
                "read_fasta",
                "encode",
                "translate",
                "seg",
                "build_lookup",
                "search",
                "output"
            ]
        );
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
