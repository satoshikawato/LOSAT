//! EXPERIMENT (`xstats` feature): global work counters, printed by `main`
//! when LOSAT_X_STATS is set.

use std::sync::atomic::{AtomicU64, Ordering};

macro_rules! counters {
    ($($name:ident),* $(,)?) => {
        $(pub static $name: AtomicU64 = AtomicU64::new(0);)*
        pub fn print() {
            if std::env::var_os("LOSAT_X_STATS").is_none() {
                return;
            }
            $(eprintln!("[X_STATS] {}={}", stringify!($name), $name.load(Ordering::Relaxed));)*
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
