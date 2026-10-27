//! PROBE (line A s16, not for a PR): a timeline of the batch flow's units under `LOSAT_X_BXPROBE=1`,
//! printed to stderr at the end of the flow as `XPROBE unit kind thread ns`.
use std::cell::Cell;
use std::sync::{Mutex, OnceLock};
use std::time::Instant;

pub(crate) const READ: u8 = 0;
pub(crate) const PREP0: u8 = 1;
pub(crate) const PREP1: u8 = 2;
pub(crate) const PLACE: u8 = 3;
pub(crate) const BEGIN: u8 = 4;
pub(crate) const QSETUP: u8 = 5;
pub(crate) const BUILT: u8 = 6;
pub(crate) const CHUNKS: u8 = 7;
pub(crate) const ENDED: u8 = 8;
pub(crate) const DONE: u8 = 9;
pub(crate) const WRITTEN: u8 = 10;
pub(crate) const BLOCK0: u8 = 11;
pub(crate) const BLOCK1: u8 = 12;
pub(crate) const INPUT_END: u8 = 13;

fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_BXPROBE").is_some_and(|v| v == "1"))
}

fn base() -> Instant {
    static BASE: OnceLock<Instant> = OnceLock::new();
    *BASE.get_or_init(Instant::now)
}

static EV: Mutex<Vec<(u32, u8, u32, u64)>> = Mutex::new(Vec::new());

thread_local! {
    static UNIT: Cell<u32> = const { Cell::new(u32::MAX) };
}

pub(crate) fn start() {
    if on() {
        let _ = base();
    }
}

pub(crate) fn mark(unit: usize, kind: u8) {
    if on() {
        let ns = base().elapsed().as_nanos() as u64;
        let thread = rayon::current_thread_index().map_or(u32::MAX, |t| t as u32);
        EV.lock().unwrap().push((unit as u32, kind, thread, ns));
    }
}

pub(crate) fn mark_here(kind: u8) {
    if on() {
        let unit = UNIT.with(Cell::get);
        if unit != u32::MAX {
            mark(unit as usize, kind);
        }
    }
}

pub(crate) fn enter(unit: usize) -> u32 {
    UNIT.with(|u| u.replace(unit as u32))
}

pub(crate) fn leave(previous: u32) {
    UNIT.with(|u| u.set(previous));
}

pub(crate) fn dump() {
    if on() {
        use std::io::Write;
        let ev = std::mem::take(&mut *EV.lock().unwrap());
        let mut out = String::new();
        for (unit, kind, thread, ns) in ev {
            out.push_str(&format!("XPROBE {unit} {kind} {thread} {ns}\n"));
        }
        let _ = std::io::stderr().write_all(out.as_bytes());
    }
}
