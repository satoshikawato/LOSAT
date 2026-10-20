//! EXPERIMENT (`xstats` feature): global work counters, printed by `main`
//! when LOSAT_X_STATS is set.
//!
//! No NCBI counterpart: these counters count calls, rows and cells of the ported code and print
//! them on request; they do not change any value NCBI computes.

use std::sync::atomic::{AtomicU64, Ordering};

macro_rules! counters {
    ($($name:ident),* $(,)?) => {
        $(pub static $name: AtomicU64 = AtomicU64::new(0);)*
        pub fn print() {
            if std::env::var_os("LOSAT_X_STATS").is_none() {
                return;
            }
            $(eprintln!("[X_STATS] {}={}", stringify!($name), $name.load(Ordering::Relaxed));)*
            #[cfg(feature = "xstats")]
            print_f1();
        }
    };
}

counters!(
    DP_SO_CALLS,
    DP_SO_ROWS,
    DP_SO_CELLS,
    DP_TB_CALLS,
    DP_TB_ROWS,
    DP_TB_CELLS,
    DP_RS_CALLS,
    DP_RS_ROWS,
    DP_RS_CELLS,
    NT_SO_CALLS,
    NT_SO_ROWS,
    NT_SO_CELLS,
    NT_TB_CALLS,
    NT_TB_ROWS,
    NT_TB_CELLS,
    REDO_CALLS,
    COMPO_ADJ_CALLS,
    NEWTON_ITERS,
    SEG_CALLS,
    SEG_RESIDUES,
    DP_SO_NS,
    DP_TB_NS,
    NT_SO_NS,
    NT_TB_NS,
    DP_FAST_CALLS,
    DP_FALLBACK_CALLS,
    REDO_NS,
    COMPO_NS,
    AHEAD_PRELIM_EVALUATED,
    AHEAD_PRELIM_TAKEN,
    AHEAD_TRACE_EVALUATED,
    AHEAD_TRACE_TAKEN,
    SEG_SHARE_CHECKED,
);

#[inline(always)]
pub fn add(counter: &AtomicU64, n: u64) {
    #[cfg(feature = "xstats")]
    counter.fetch_add(n, Ordering::Relaxed);
    #[cfg(not(feature = "xstats"))]
    let _ = (counter, n);
}

/// Wall-clock accumulation (only with the `xstats` feature).
#[inline(always)]
pub fn now() -> Option<std::time::Instant> {
    #[cfg(feature = "xstats")]
    return Some(std::time::Instant::now());
    #[cfg(not(feature = "xstats"))]
    None
}

#[inline(always)]
pub fn add_ns(counter: &AtomicU64, t0: Option<std::time::Instant>) {
    if let Some(t0) = t0 {
        counter.fetch_add(t0.elapsed().as_nanos() as u64, Ordering::Relaxed);
    }
}

// EXPERIMENT (perf-of F1, `xstats` feature): the work ledger of BLASTP (`F1_*`). The statics exist in
// every build (no code); the adds, the thread-local accumulators and the print are `xstats` only.
macro_rules! f1_counters {
    ($($name:ident),* $(,)?) => {
        $(pub static $name: AtomicU64 = AtomicU64::new(0);)*
        #[cfg(feature = "xstats")]
        fn print_f1() {
            $(eprintln!("[X_STATS] {}={}", stringify!($name), $name.load(Ordering::Relaxed));)*
        }
    };
}

f1_counters!(
    F1_SCAN_POS,
    F1_SCAN_CELLS,
    F1_HITS,
    F1_HIT_MASKED,
    F1_HIT_FLAGRESET,
    F1_HIT_TOOFAR,
    F1_HIT_OVERLAP,
    F1_HIT_CTXREJECT,
    F1_HIT_EXTEND,
    F1_UNG_CALLS,
    F1_UNG_LEFT_CELLS,
    F1_UNG_RIGHT_CALLS,
    F1_UNG_RIGHT_CELLS,
    F1_UNG_SAVED,
    F1_PG_INIT,
    F1_PG_CONTAINED,
    F1_PG_EXT,
    F1_PG_RS_CALLS,
    F1_PG_RS_ROWS,
    F1_PG_RS_CELLS,
    F1_PG_RS_REDO,
    F1_PG_SO_CALLS,
    F1_PG_SO_SIMD_CELLS,
    F1_TB_CALLS,
    F1_TB_SIMD_CELLS,
    F1_RD_PAIRS,
    F1_RD_HSPS,
);

// Slots of the per-thread accumulators (flushed once per subject by `f1_flush`): the counters
// that are added once per hit or per ungapped extension.
pub const F1T_SCAN_CELLS: usize = 0;
pub const F1T_HITS: usize = 1;
pub const F1T_MASKED: usize = 2;
pub const F1T_FLAGRESET: usize = 3;
pub const F1T_TOOFAR: usize = 4;
pub const F1T_OVERLAP: usize = 5;
pub const F1T_CTXREJECT: usize = 6;
pub const F1T_EXTEND: usize = 7;
pub const F1T_UNG_CALLS: usize = 8;
pub const F1T_UNG_LEFT_CELLS: usize = 9;
pub const F1T_UNG_RIGHT_CALLS: usize = 10;
pub const F1T_UNG_RIGHT_CELLS: usize = 11;
pub const F1T_UNG_SAVED: usize = 12;
#[cfg(feature = "xstats")]
const F1T_SLOTS: usize = 13;

#[cfg(feature = "xstats")]
thread_local! {
    static F1_TL: [std::cell::Cell<u64>; F1T_SLOTS] = const { [const { std::cell::Cell::new(0) }; F1T_SLOTS] };
}

/// Adds `n` to a per-thread slot (`F1T_*`).
#[inline(always)]
pub fn f1_tl(slot: usize, n: u64) {
    #[cfg(feature = "xstats")]
    F1_TL.with(|t| t[slot].set(t[slot].get() + n));
    #[cfg(not(feature = "xstats"))]
    let _ = (slot, n);
}

/// Moves the per-thread slots into the global counters.
#[inline(always)]
pub fn f1_flush() {
    #[cfg(feature = "xstats")]
    {
        let targets: [&AtomicU64; F1T_SLOTS] = [
            &F1_SCAN_CELLS,
            &F1_HITS,
            &F1_HIT_MASKED,
            &F1_HIT_FLAGRESET,
            &F1_HIT_TOOFAR,
            &F1_HIT_OVERLAP,
            &F1_HIT_CTXREJECT,
            &F1_HIT_EXTEND,
            &F1_UNG_CALLS,
            &F1_UNG_LEFT_CELLS,
            &F1_UNG_RIGHT_CALLS,
            &F1_UNG_RIGHT_CELLS,
            &F1_UNG_SAVED,
        ];
        F1_TL.with(|t| {
            for (slot, target) in t.iter().zip(targets) {
                let v = slot.replace(0);
                if v != 0 {
                    target.fetch_add(v, Ordering::Relaxed);
                }
            }
        });
    }
}
