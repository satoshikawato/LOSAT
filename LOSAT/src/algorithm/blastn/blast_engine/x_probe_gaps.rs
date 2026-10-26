//! PROBE (line ST s2, branch probe/st-bnahead, never a PR): clocks of the steps of blastn's batch
//! loop. With `LOSAT_PROBE_BNGAPS` set, each `mark` keeps (batch, step, wall ns, process CPU ns)
//! and `dump` writes them to stderr as `[PROBE_BNGAPS] batch step wall_ns cpu_ns`; a step's
//! CPU over its wall time is the number of threads busy in it. No NCBI counterpart.

use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::{Mutex, OnceLock};
use std::time::Instant;

type Event = (u32, &'static str, u64, u64);

static ON: OnceLock<bool> = OnceLock::new();
static T0: OnceLock<Instant> = OnceLock::new();
static BATCH: AtomicU32 = AtomicU32::new(0);
static EVENTS: Mutex<Vec<Event>> = Mutex::new(Vec::new());

pub(crate) fn on() -> bool {
    *ON.get_or_init(|| std::env::var_os("LOSAT_PROBE_BNGAPS").is_some())
}

#[cfg(target_os = "linux")]
fn cpu_ns() -> u64 {
    #[repr(C)]
    struct Timespec {
        sec: i64,
        nsec: i64,
    }
    extern "C" {
        fn clock_gettime(clock: i32, time: *mut Timespec) -> i32;
    }
    // CLOCK_PROCESS_CPUTIME_ID: every thread of the process.
    let mut time = Timespec { sec: 0, nsec: 0 };
    unsafe { clock_gettime(2, &mut time) };
    time.sec as u64 * 1_000_000_000 + time.nsec as u64
}

#[cfg(not(target_os = "linux"))]
fn cpu_ns() -> u64 {
    0
}

pub(crate) fn next_batch() {
    if on() {
        BATCH.fetch_add(1, Ordering::Relaxed);
    }
}

pub(crate) fn mark(step: &'static str) {
    if !on() {
        return;
    }
    let wall = T0.get_or_init(Instant::now).elapsed().as_nanos() as u64;
    let cpu = cpu_ns();
    let batch = BATCH.load(Ordering::Relaxed);
    EVENTS
        .lock()
        .unwrap_or_else(|e| e.into_inner())
        .push((batch, step, wall, cpu));
}

pub(crate) fn dump() {
    if !on() {
        return;
    }
    let events = std::mem::take(&mut *EVENTS.lock().unwrap_or_else(|e| e.into_inner()));
    let mut text = String::new();
    for (batch, step, wall, cpu) in events {
        text.push_str(&format!("[PROBE_BNGAPS] {batch} {step} {wall} {cpu}\n"));
    }
    eprint!("{text}");
}
