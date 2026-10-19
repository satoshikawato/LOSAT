//! EXPERIMENT (R7_PROBE=gt; line R's R7 probe, `v4-R/DESIGN_R7.md`; not for merge): the time of
//! each greedy traceback call of blastn (`greedy_gapped_alignment_with_traceback`), by the
//! alignment length it returns, in power-of-two buckets. Printed at exit as `[R7] gt ...` lines.
//! No value changes: the call site only reads a clock around the call when the probe is on.

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;
use std::time::Instant;

const BUCKETS: usize = 24;
static CALLS: [AtomicU64; BUCKETS] = [const { AtomicU64::new(0) }; BUCKETS];
static NANOS: [AtomicU64; BUCKETS] = [const { AtomicU64::new(0) }; BUCKETS];
static LETTERS: [AtomicU64; BUCKETS] = [const { AtomicU64::new(0) }; BUCKETS];
static PRECOMPUTED: AtomicU64 = AtomicU64::new(0);
static NONE: AtomicU64 = AtomicU64::new(0);
static MAX_NS: AtomicU64 = AtomicU64::new(0);
static MAX_LEN: AtomicU64 = AtomicU64::new(0);

fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| {
        std::env::var("R7_PROBE")
            .map(|v| v.split(',').any(|t| t.trim() == "gt"))
            .unwrap_or(false)
    })
}

/// The start of a call (None when the probe is off).
#[inline(never)]
pub(crate) fn start() -> Option<Instant> {
    on().then(Instant::now)
}

/// The end of a call that returned an alignment of `len` columns (`None`: no alignment);
/// `precomputed`: the call was made ahead (LOSAT_X_GREEDYSPEC) and not here.
#[inline(never)]
pub(crate) fn stop(t0: Instant, len: Option<usize>, precomputed: bool) {
    let ns = t0.elapsed().as_nanos() as u64;
    if precomputed {
        PRECOMPUTED.fetch_add(1, Ordering::Relaxed);
    }
    let Some(len) = len else {
        NONE.fetch_add(1, Ordering::Relaxed);
        return;
    };
    let b = (usize::BITS - len.max(1).leading_zeros()) as usize - 1;
    let b = b.min(BUCKETS - 1);
    CALLS[b].fetch_add(1, Ordering::Relaxed);
    NANOS[b].fetch_add(ns, Ordering::Relaxed);
    LETTERS[b].fetch_add(len as u64, Ordering::Relaxed);
    if ns > MAX_NS.fetch_max(ns, Ordering::Relaxed) {
        MAX_LEN.store(len as u64, Ordering::Relaxed);
    }
}

/// Prints the buckets (called at exit).
pub(crate) fn print() {
    if !on() {
        return;
    }
    let r = Ordering::Relaxed;
    for b in 0..BUCKETS {
        let c = CALLS[b].load(r);
        if c == 0 {
            continue;
        }
        eprintln!(
            "[R7] gt len={}-{} calls={} s={:.6} letters={}",
            1u64 << b,
            (1u64 << (b + 1)) - 1,
            c,
            NANOS[b].load(r) as f64 / 1e9,
            LETTERS[b].load(r)
        );
    }
    eprintln!(
        "[R7] gt none={} precomputed={} max_call_s={:.6} max_call_len={}",
        NONE.load(r),
        PRECOMPUTED.load(r),
        MAX_NS.load(r) as f64 / 1e9,
        MAX_LEN.load(r)
    );
}
