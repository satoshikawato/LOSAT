//! PROBE (v4-B10 B11e, never merged): with `LOSAT_X_B11EPROBE=1`, tblastn prints `[B11E]` lines
//! with the thread and the start and length of every preliminary subject, every (frame, chunk)
//! unit, every redo item and redo chunk, and the phases of each batch, so that a per-thread
//! timeline can be drawn. Counting only; no value changes.
use std::cell::{Cell, RefCell};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;
use std::time::Instant;

static ON: OnceLock<bool> = OnceLock::new();
static START: OnceLock<Instant> = OnceLock::new();
/// The batch being searched (set at its start).
pub(crate) static BATCH: AtomicUsize = AtomicUsize::new(0);
static NEXT_BATCH: AtomicUsize = AtomicUsize::new(0);

pub(crate) fn on() -> bool {
    *ON.get_or_init(|| {
        let on = std::env::var("LOSAT_X_B11EPROBE").is_ok_and(|v| v == "1");
        if on {
            START.get_or_init(Instant::now);
        }
        on
    })
}

/// Seconds since the probe was first asked.
pub(crate) fn now() -> f64 {
    START.get_or_init(Instant::now).elapsed().as_secs_f64()
}

pub(crate) fn thread() -> i64 {
    rayon::current_thread_index().map_or(-1, |t| t as i64)
}

pub(crate) fn ns(since: Instant) -> u64 {
    since.elapsed().as_nanos() as u64
}

pub(crate) fn ms(ns: u64) -> f64 {
    ns as f64 / 1e6
}

/// PROBE (B11e): `LOSAT_X_B11ELPTMAX` in place of `X_REDO_LPT_MAX_THREADS` (scheduling only).
pub(crate) fn lpt_max(default: usize) -> usize {
    static MAX: OnceLock<usize> = OnceLock::new();
    *MAX.get_or_init(|| {
        std::env::var("LOSAT_X_B11ELPTMAX")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(default)
    })
}

pub(crate) fn batch() -> usize {
    BATCH.load(Ordering::Relaxed)
}

/// A new batch: its index, and a phase line.
pub(crate) fn batch_start(queries: usize, subjects: usize) {
    let b = NEXT_BATCH.fetch_add(1, Ordering::SeqCst);
    BATCH.store(b, Ordering::SeqCst);
    eprintln!(
        "[B11E] batch b={b} phase=start t={:.4} queries={queries} subjects={subjects}",
        now()
    );
}

pub(crate) fn phase(name: &str) {
    eprintln!("[B11E] batch b={} phase={name} t={:.4}", batch(), now());
}

/// The parts of the last `x_preliminary_protein_hsps_by_unit` call of this thread.
#[derive(Default, Clone, Copy, Debug)]
pub(crate) struct PrelimParts {
    pub plan_ns: u64,
    pub units_ns: u64,
    pub merge_ns: u64,
    pub units: usize,
    pub parallel: bool,
}

/// The parts of the unit running on this thread.
#[derive(Default, Clone, Copy, Debug)]
pub(crate) struct UnitParts {
    pub scan_ns: u64,
    pub gapped_ns: u64,
    pub purge_ns: u64,
    pub hits: usize,
    pub hsps: usize,
}

/// The parts of the last `link_preliminary_hsps` call of this thread.
#[derive(Default, Clone, Debug)]
pub(crate) struct LinkParts {
    pub n: usize,
    pub prep_ns: u64,
    pub sort_ns: u64,
    pub window_ns: u64,
    pub loop_ns: u64,
    pub final_ns: u64,
    pub heads: u64,
    /// Per context: (HSPs, heads, loop ns).
    pub by_context: Vec<(usize, u64, u64)>,
}

impl LinkParts {
    /// "ctx:hsps:heads:ms" of the three contexts with the most loop time.
    pub(crate) fn top(&self) -> String {
        let mut order: Vec<usize> = (0..self.by_context.len())
            .filter(|&c| self.by_context[c].0 > 0)
            .collect();
        order.sort_by(|&a, &b| self.by_context[b].2.cmp(&self.by_context[a].2));
        order
            .iter()
            .take(3)
            .map(|&c| {
                let (n, heads, ns) = self.by_context[c];
                format!("{c}:{n}:{heads}:{:.1}", ms(ns))
            })
            .collect::<Vec<_>>()
            .join(",")
    }

    pub(crate) fn contexts(&self) -> usize {
        self.by_context.iter().filter(|c| c.0 > 0).count()
    }
}

thread_local! {
    pub(crate) static PRELIM: Cell<PrelimParts> = Cell::new(PrelimParts::default());
    pub(crate) static UNIT: Cell<UnitParts> = Cell::new(UnitParts::default());
    pub(crate) static LINK: RefCell<LinkParts> = RefCell::new(LinkParts::default());
    /// The subject whose preliminary stage this thread runs (set by `subject_core`).
    pub(crate) static OID: Cell<usize> = const { Cell::new(usize::MAX) };
}

pub(crate) fn take_link() -> LinkParts {
    LINK.with(|r| std::mem::take(&mut *r.borrow_mut()))
}

/// Prints one `[B11E] redo` line when dropped (the redo item's end, whatever it returns).
pub(crate) struct ItemGuard {
    context: usize,
    oid: usize,
    hsps: usize,
    t0: f64,
    start: Instant,
}

impl ItemGuard {
    pub(crate) fn new(context: usize, oid: usize, hsps: usize) -> Self {
        Self {
            context,
            oid,
            hsps,
            t0: now(),
            start: Instant::now(),
        }
    }
}

impl Drop for ItemGuard {
    fn drop(&mut self) {
        eprintln!(
            "[B11E] redo b={} ctx={} oid={} thr={} t0={:.4} dt={:.2} n={}",
            batch(),
            self.context,
            self.oid,
            thread(),
            self.t0,
            ms(ns(self.start)),
            self.hsps
        );
    }
}
