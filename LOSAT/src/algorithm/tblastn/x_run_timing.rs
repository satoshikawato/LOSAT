//! Diagnostic (LOSAT_TIMING): where the wall time of one TBLASTN run goes, and how many threads
//! worked in each part of it, so that the parts that do not run in parallel have numbers.
//!
//! No NCBI counterpart: NCBI BLAST+ has no such timer. The module reads a monotonic clock and the
//! process's CPU time and prints them on standard error; it computes no value NCBI computes and no
//! value the search reads. Without `LOSAT_TIMING` it reads neither.
//!
//! # What is measured
//!
//! The command line's run (`args.rs` `run`) starts one timer ([`start`]). The thread that started
//! it calls [`lap`] at the end of each step, which adds the wall time and the CPU time of the
//! process since the previous lap to the step's part; laps from other threads (query batches
//! searched side by side) are ignored, and so are laps without a timer (the library entry and the
//! tests). Parts of the same name (one per query batch or query set) are added. When the run ends
//! (the guard is dropped, after the subjects are freed) one line is printed:
//!
//! `[TIMING] tbn_run: threads=N total=..s <part> ..s busy .. ...`
//!
//! `busy` is CPU / wall of the part: the number of threads that worked in it on average (the CPU
//! time of all the threads of the process, ended ones too, from `/proc/self/stat` at 10 ms
//! resolution; Linux only, `-` elsewhere). The parts add up to `total`; the time before the timer
//! starts (the options) and after the line (the process exit) is outside.
//!
//! The subject loop (`run_preliminary_subjects`) and the composition redo split their part in
//! two ([`lap_split`]). The subject loop: `prelim_search`, the time the calling thread waits for
//! the pool or searches subjects itself (one thread; with `LOSAT_X_TBNOIDNEXT` only its own
//! searches, its waits for a chunk stay in the rest), and `prelim_reduce`, the rest (the
//! OID-ordered collection of the subjects' lists). The redo: `redo_pool`, the waits for the pool,
//! and `redo_serial`, the rest (the early-termination tests, the heap updates, the matches redone
//! in order). `busy` of a part the calling thread spent searching itself is `-`. The subject loop
//! also prints, per query set, [`XSubjectLoop`]'s line:
//!
//! `[TIMING] tbn_subjects: subjects=.. threads=.. entries=.. pool=..s pool_busy=.. work=..s work_busy=.. own=..s work_per_subject=..us entry_ms=..`
//!
//! `work` is the sum over the subjects of the wall time of each subject's search (on whatever
//! thread ran it). The search runs on a thread of the pool (`with_search_pool`), so the calling
//! thread searches subjects of its own entries too; `work_busy` = (work - own) / pool is the number
//! of threads that were inside a subject on average during the entries, the calling thread
//! included. Against `pool_busy` it separates time a worker spends between subjects (the pool's own
//! scheduling, waiting at the end of an entry) from the work; against the one-thread run's `work`
//! it shows how much a subject's search slows down when many threads search at once. `own` is the
//! calling thread's searches outside the entries (the one-thread loop, the flow of
//! `LOSAT_X_TBNOIDNEXT`). `entry_ms` lists the wall time of each of the first 40 entries.

use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Mutex, OnceLock};
use std::time::{Duration, Instant};

// No NCBI counterpart: reads LOSAT_TIMING once; it does not change any value NCBI computes.
/// Whether `LOSAT_TIMING` is set.
pub(super) fn enabled() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_TIMING").is_some())
}

// No NCBI counterpart: the CPU time of the process for the timer only.
/// The CPU time (user + system) of every thread of the process, ended ones too; `None` where
/// `/proc/self/stat` is not there.
pub(super) fn process_cpu() -> Option<Duration> {
    #[cfg(target_os = "linux")]
    {
        cpu_of_stat(&std::fs::read_to_string("/proc/self/stat").ok()?)
    }
    #[cfg(not(target_os = "linux"))]
    {
        None
    }
}

/// utime + stime of a `/proc/<pid>/stat` line.
#[cfg_attr(not(target_os = "linux"), allow(dead_code))]
fn cpu_of_stat(stat: &str) -> Option<Duration> {
    // The second field is the command name in parentheses, which may hold spaces and
    // parentheses; utime and stime (fields 14 and 15, proc(5)) are the 12th and 13th after the
    // state, field 3.
    let rest = stat.get(stat.rfind(')')? + 2..)?;
    let mut fields = rest.split(' ');
    let utime: u64 = fields.nth(11)?.parse().ok()?;
    let stime: u64 = fields.next()?.parse().ok()?;
    // In clock ticks of USER_HZ, which is 100 on the targets LOSAT builds for.
    Some(Duration::from_millis(
        utime.saturating_add(stime).saturating_mul(10),
    ))
}

// No NCBI counterpart: one named part of the run for the timer.
struct XPart {
    name: &'static str,
    wall: Duration,
    cpu: Option<Duration>,
}

// No NCBI counterpart: the state of the run's timer.
struct XRun {
    thread: std::thread::ThreadId,
    threads: usize,
    start: Instant,
    last: Instant,
    last_cpu: Option<Duration>,
    parts: Vec<XPart>,
}

impl XRun {
    fn add(&mut self, name: &'static str, wall: Duration, cpu: Option<Duration>) {
        match self.parts.iter_mut().find(|part| part.name == name) {
            Some(part) => {
                part.wall += wall;
                part.cpu = part.cpu.zip(cpu).map(|(a, b)| a + b);
            }
            None => self.parts.push(XPart { name, wall, cpu }),
        }
    }
}

static ACTIVE: AtomicBool = AtomicBool::new(false);
static RUN: Mutex<Option<XRun>> = Mutex::new(None);

// No NCBI counterpart: ends the run's timer and prints its line.
/// The run's timer; dropping it adds the last part (`drop`) and prints the line.
pub(super) struct XRunGuard(());

/// Starts the run's timer when `LOSAT_TIMING` is set; the guard is the first local of the run, so
/// that it is dropped after everything the run owns.
pub(super) fn start(threads: usize) -> Option<XRunGuard> {
    if !enabled() {
        return None;
    }
    let now = Instant::now();
    let mut run = RUN.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
    *run = Some(XRun {
        thread: std::thread::current().id(),
        threads,
        start: now,
        last: now,
        last_cpu: process_cpu(),
        parts: Vec::new(),
    });
    ACTIVE.store(true, Ordering::Release);
    Some(XRunGuard(()))
}

// No NCBI counterpart: the run's timer, as `lap` and `lap_split` share it.
fn with_run(work: impl FnOnce(&mut XRun, Instant, Option<Duration>)) {
    if !ACTIVE.load(Ordering::Acquire) {
        return;
    }
    let mut run = RUN.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
    let Some(run) = run.as_mut() else {
        return;
    };
    if run.thread != std::thread::current().id() {
        return;
    }
    let now = Instant::now();
    let cpu = process_cpu();
    work(run, now, cpu);
    run.last = now;
    run.last_cpu = cpu;
}

/// Adds the time since the previous lap to the part `name`.
pub(super) fn lap(name: &'static str) {
    with_run(|run, now, cpu| {
        let wall = now - run.last;
        let cpu = cpu
            .zip(run.last_cpu)
            .map(|(now, last)| now.saturating_sub(last));
        run.add(name, wall, cpu);
    });
}

/// Adds the time since the previous lap to two parts: `inner` gets `inner_wall` / `inner_cpu`
/// (measured by the caller inside the interval), `rest` gets the remainder.
pub(super) fn lap_split(
    rest: &'static str,
    inner: &'static str,
    inner_wall: Duration,
    inner_cpu: Option<Duration>,
) {
    with_run(|run, now, cpu| {
        let wall = now - run.last;
        let cpu = cpu
            .zip(run.last_cpu)
            .map(|(now, last)| now.saturating_sub(last));
        run.add(inner, inner_wall.min(wall), inner_cpu);
        run.add(
            rest,
            wall.saturating_sub(inner_wall),
            cpu.zip(inner_cpu)
                .map(|(all, inner)| all.saturating_sub(inner)),
        );
    });
}

fn busy(wall: Duration, cpu: Option<Duration>) -> String {
    match cpu {
        Some(cpu) if !wall.is_zero() => format!("{:.1}", cpu.as_secs_f64() / wall.as_secs_f64()),
        _ => "-".to_string(),
    }
}

impl Drop for XRunGuard {
    fn drop(&mut self) {
        lap("drop");
        ACTIVE.store(false, Ordering::Release);
        let mut run = RUN.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
        let Some(run) = run.take() else {
            return;
        };
        let mut line = format!(
            "[TIMING] tbn_run: threads={} total={:.3}s",
            run.threads,
            (run.last - run.start).as_secs_f64()
        );
        for part in &run.parts {
            line.push_str(&format!(
                " {} {:.3}s busy {}",
                part.name,
                part.wall.as_secs_f64(),
                busy(part.wall, part.cpu)
            ));
        }
        eprintln!("{line}");
    }
}

// No NCBI counterpart: the waits of the calling thread for the pool, for the timer.
/// The time the calling thread waits for the pool (the subject loop's entries, the redo's
/// batches); nothing is read without `LOSAT_TIMING`. Atomics only, so that the closures that run
/// on the pool can share the value with the calling thread.
pub(super) struct XPoolWaits {
    on: bool,
    calls: AtomicU64,
    wall_nanos: AtomicU64,
    cpu_nanos: AtomicU64,
    cpu_missing: AtomicBool,
    /// The wall time of each of the first entries, in milliseconds.
    first: Mutex<Vec<f64>>,
}

/// How many entries [`XSubjectLoop`]'s line lists one by one.
const LISTED_ENTRIES: usize = 40;

fn nanos(duration: Duration) -> u64 {
    u64::try_from(duration.as_nanos()).unwrap_or(u64::MAX)
}

impl XPoolWaits {
    pub(super) fn new() -> Self {
        Self {
            on: enabled(),
            calls: AtomicU64::new(0),
            wall_nanos: AtomicU64::new(0),
            cpu_nanos: AtomicU64::new(0),
            cpu_missing: AtomicBool::new(false),
            first: Mutex::new(Vec::new()),
        }
    }

    /// One entry into the pool, timed when the timer is on.
    pub(super) fn entry<R>(&self, work: impl FnOnce() -> R) -> R {
        if !self.on {
            return work();
        }
        let start = Instant::now();
        let start_cpu = process_cpu();
        let result = work();
        let wall = start.elapsed();
        self.calls.fetch_add(1, Ordering::Relaxed);
        self.wall_nanos.fetch_add(nanos(wall), Ordering::Relaxed);
        let mut first = self
            .first
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        if first.len() < LISTED_ENTRIES {
            first.push(wall.as_secs_f64() * 1e3);
        }
        drop(first);
        match process_cpu().zip(start_cpu) {
            Some((now, then)) => {
                self.cpu_nanos
                    .fetch_add(nanos(now.saturating_sub(then)), Ordering::Relaxed);
            }
            None => self.cpu_missing.store(true, Ordering::Relaxed),
        }
        result
    }

    fn wall(&self) -> Duration {
        Duration::from_nanos(self.wall_nanos.load(Ordering::Relaxed))
    }

    fn cpu(&self) -> Option<Duration> {
        (!self.cpu_missing.load(Ordering::Relaxed))
            .then(|| Duration::from_nanos(self.cpu_nanos.load(Ordering::Relaxed)))
    }

    /// Adds the laps `inner` (the waits) and `rest` (the remainder since the previous lap).
    pub(super) fn finish(&self, rest: &'static str, inner: &'static str) {
        if self.on {
            lap_split(rest, inner, self.wall(), self.cpu());
        }
    }
}

// No NCBI counterpart: the subject loop's numbers for the timer.
/// The subject loop's waits for the pool and the subjects' own search time; nothing is read
/// without `LOSAT_TIMING`.
pub(super) struct XSubjectLoop {
    pool: XPoolWaits,
    owner: std::thread::ThreadId,
    /// Set while the calling thread is in an entry: the search runs on a pool thread
    /// (`with_search_pool`), so the calling thread searches subjects of its own entries too.
    in_entry: AtomicBool,
    work_nanos: AtomicU64,
    /// The calling thread's searches outside the entries (the one-thread loop, the flow).
    own_nanos: AtomicU64,
}

impl XSubjectLoop {
    pub(super) fn new() -> Self {
        Self {
            pool: XPoolWaits::new(),
            owner: std::thread::current().id(),
            in_entry: AtomicBool::new(false),
            work_nanos: AtomicU64::new(0),
            own_nanos: AtomicU64::new(0),
        }
    }

    /// One entry into the pool (`pool.install`), timed when the timer is on.
    pub(super) fn entry<R>(&self, work: impl FnOnce() -> R) -> R {
        self.in_entry.store(true, Ordering::Relaxed);
        let result = self.pool.entry(work);
        self.in_entry.store(false, Ordering::Relaxed);
        result
    }

    /// One subject's search, timed when the timer is on (from any thread).
    pub(super) fn subject<R>(&self, work: impl FnOnce() -> R) -> R {
        if !self.pool.on {
            return work();
        }
        let start = Instant::now();
        let result = work();
        let spent = nanos(start.elapsed());
        self.work_nanos.fetch_add(spent, Ordering::Relaxed);
        if std::thread::current().id() == self.owner && !self.in_entry.load(Ordering::Relaxed) {
            self.own_nanos.fetch_add(spent, Ordering::Relaxed);
        }
        result
    }

    /// Prints the line of the loop and adds its laps (`prelim_search`, `prelim_reduce`).
    pub(super) fn finish(&self, subjects: usize, threads: usize) {
        let pool = &self.pool;
        if !pool.on {
            return;
        }
        let work = Duration::from_nanos(self.work_nanos.load(Ordering::Relaxed));
        let own = Duration::from_nanos(self.own_nanos.load(Ordering::Relaxed));
        let per_subject = if subjects == 0 {
            0.0
        } else {
            work.as_secs_f64() * 1e6 / subjects as f64
        };
        let first: Vec<String> = pool
            .first
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .iter()
            .map(|ms| format!("{ms:.1}"))
            .collect();
        eprintln!(
            "[TIMING] tbn_subjects: subjects={} threads={} entries={} pool={:.3}s pool_busy={} work={:.3}s work_busy={} own={:.3}s work_per_subject={:.1}us entry_ms={}",
            subjects,
            threads,
            pool.calls.load(Ordering::Relaxed),
            pool.wall().as_secs_f64(),
            busy(pool.wall(), pool.cpu()),
            work.as_secs_f64(),
            busy(pool.wall(), Some(work.saturating_sub(own))),
            own.as_secs_f64(),
            per_subject,
            first.join(",")
        );
        lap_split(
            "prelim_reduce",
            "prelim_search",
            pool.wall() + own,
            pool.cpu().filter(|_| own.is_zero()),
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_cpu_time_is_utime_plus_stime_after_a_command_name_with_spaces() {
        // proc(5): pid (comm) state ppid pgrp session tty_nr tpgid flags minflt cminflt majflt
        // cmajflt utime stime ...
        let line = "4242 (LOSAT (x) y) R 1 2 3 4 5 6 7 8 9 10 250 37 0 0 20 0 33 0 1 2 3\n";
        assert_eq!(cpu_of_stat(line), Some(Duration::from_millis(2_870)));
        assert_eq!(cpu_of_stat("4242 (LOSAT) R 1 2"), None);
        #[cfg(target_os = "linux")]
        assert!(process_cpu().is_some());
    }

    #[test]
    fn parts_of_the_same_name_are_added_in_the_order_they_first_ran() {
        let now = Instant::now();
        let mut run = XRun {
            thread: std::thread::current().id(),
            threads: 1,
            start: now,
            last: now,
            last_cpu: None,
            parts: Vec::new(),
        };
        let ms = Duration::from_millis;
        run.add("setup", ms(2), Some(ms(1)));
        run.add("search", ms(10), Some(ms(70)));
        run.add("setup", ms(3), Some(ms(2)));
        run.add("search", ms(5), None);
        let parts: Vec<_> = run.parts.iter().map(|p| (p.name, p.wall, p.cpu)).collect();
        assert_eq!(
            parts,
            vec![("setup", ms(5), Some(ms(3))), ("search", ms(15), None)]
        );
        assert_eq!(busy(ms(10), Some(ms(70))), "7.0");
        assert_eq!(busy(ms(10), None), "-");
    }
}
