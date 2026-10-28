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
    *BASE.get_or_init(|| {
        let b = Instant::now();
        if on() {
            start_rss(b);
        }
        b
    })
}

/// RssAnon samples every `LOSAT_XTLINE_RSS` ms (s22 RSS probe): kind "rss", id = RssAnon KiB
/// (statm resident - shared), thread -2.
fn start_rss(b: Instant) {
    let ms: u64 = match std::env::var("LOSAT_XTLINE_RSS")
        .ok()
        .and_then(|v| v.parse().ok())
    {
        Some(v) if v > 0 => v,
        _ => return,
    };
    std::thread::spawn(move || loop {
        if let Ok(s) = std::fs::read_to_string("/proc/self/statm") {
            let mut it = s.split_ascii_whitespace().skip(1);
            let res: u64 = it.next().and_then(|x| x.parse().ok()).unwrap_or(0);
            let sh: u64 = it.next().and_then(|x| x.parse().ok()).unwrap_or(0);
            let t = b.elapsed().as_nanos() as u64;
            let kib = (res.saturating_sub(sh) * 4) as u32;
            EV.lock()
                .unwrap_or_else(|e| e.into_inner())
                .push(("rss", kib, -2, t, t));
        }
        std::thread::sleep(std::time::Duration::from_millis(ms));
    });
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

pub fn dump() {
    if !on() {
        return;
    }
    let path = std::env::var("LOSAT_XTLINE").unwrap_or_default();
    let ev = EV.lock().unwrap_or_else(|e| e.into_inner());
    let mut s = String::new();
    for (k, id, th, a, z) in ev.iter() {
        s.push_str(&format!("{k}\t{id}\t{th}\t{a}\t{z}\n"));
    }
    let _ = std::fs::write(path, s);
}

/// PROBE (b26 timeline): records `kind` from `t0` to the drop.
pub(crate) struct XRecOnDrop(
    pub(crate) &'static str,
    pub(crate) u32,
    pub(crate) std::time::Instant,
);

impl Drop for XRecOnDrop {
    fn drop(&mut self) {
        rec(self.0, self.1, self.2, std::time::Instant::now());
    }
}
