//! TEMPORARY (plan v4 B21 measurement, probe branch only, never for dev): a timeline of leaf work
//! by thread (`LOSAT_XTLINE=<file>`), as the B5 patch had.
use std::sync::{Mutex, OnceLock};
use std::time::Instant;

static BASE: OnceLock<Instant> = OnceLock::new();
type Ev = (&'static str, u32, i64, u64, u64);
static EV: Mutex<Vec<Ev>> = Mutex::new(Vec::new());

pub fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_XTLINE").is_some())
}

pub fn base() -> Instant {
    *BASE.get_or_init(Instant::now)
}

#[inline]
pub fn rec(kind: &'static str, id: u32, t0: Instant, t1: Instant) {
    if !on() {
        return;
    }
    let b = base();
    let th = rayon::current_thread_index()
        .map(|i| i as i64)
        .unwrap_or(-1);
    let a = t0.saturating_duration_since(b).as_nanos() as u64;
    let z = t1.saturating_duration_since(b).as_nanos() as u64;
    EV.lock()
        .unwrap_or_else(|e| e.into_inner())
        .push((kind, id, th, a, z));
}

/// Kernel time buckets (ns): 0 pass setup (removals), 1 index 0, 2 index 1, 3 tree update, 4 rounds.
static NS: [std::sync::atomic::AtomicU64; 5] = [const { std::sync::atomic::AtomicU64::new(0) }; 5];

#[inline]
pub fn add_ns(bucket: usize, ns: u64) {
    NS[bucket].fetch_add(ns, std::sync::atomic::Ordering::Relaxed);
}

pub fn dump() {
    if !on() {
        return;
    }
    let ns: Vec<f64> = NS
        .iter()
        .map(|a| a.load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9)
        .collect();
    eprintln!(
        "[XTLINE] kernel pass_setup {:.3}s index0 {:.3}s index1 {:.3}s tree {:.3}s rounds_total {:.3}s",
        ns[0], ns[1], ns[2], ns[3], ns[4]
    );
    let path = std::env::var("LOSAT_XTLINE").unwrap_or_default();
    let ev = EV.lock().unwrap_or_else(|e| e.into_inner());
    let mut s = String::new();
    for (k, id, th, a, z) in ev.iter() {
        s.push_str(&format!("{k}\t{id}\t{th}\t{a}\t{z}\n"));
    }
    let _ = std::fs::write(path, s);
}
