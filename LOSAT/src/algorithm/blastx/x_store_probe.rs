//! EXPERIMENT (LOSAT_X_BXSTOREPROBE=1): a measurement of how BLASTX `-db` uses the stored
//! database sequences (`common::database::StoredProteins`), for the premise check of a
//! block-wise store (v4-D12, D13). Not shipped code: every function here does nothing unless the
//! variable is `1`, and the callers keep their calls out of the per-subject loops or behind the
//! `on()` check. Lines start with `XBXSTOREPROBE` and go to standard error.
//!
//! No NCBI counterpart: it counts reads of LOSAT's own store and changes no value NCBI computes.

use std::collections::HashSet;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering::Relaxed};
use std::sync::{Mutex, OnceLock};

/// `LOSAT_X_BXSTOREPROBE=1`, read once.
pub(crate) fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var("LOSAT_X_BXSTOREPROBE").is_ok_and(|value| value.trim() == "1"))
}

/// The process's RssAnon in KiB (/proc/self/status), 0 where it cannot be read.
pub(crate) fn rss_anon_kib() -> u64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|status| {
            status.lines().find_map(|line| {
                line.strip_prefix("RssAnon:")
                    .and_then(|rest| rest.trim().trim_end_matches("kB").trim().parse().ok())
            })
        })
        .unwrap_or(0)
}

/// A count of things in flight and the largest it reached.
pub(crate) struct Gauge {
    current: AtomicUsize,
    max: AtomicUsize,
}

impl Gauge {
    const fn new() -> Self {
        Self {
            current: AtomicUsize::new(0),
            max: AtomicUsize::new(0),
        }
    }

    fn enter(&'static self) -> InFlight {
        let now = self.current.fetch_add(1, Relaxed) + 1;
        self.max.fetch_max(now, Relaxed);
        InFlight(self)
    }

    fn max(&self) -> usize {
        self.max.load(Relaxed)
    }
}

/// One thing in flight; leaves its gauge when dropped.
pub(crate) struct InFlight(&'static Gauge);

impl Drop for InFlight {
    fn drop(&mut self) {
        self.0.current.fetch_sub(1, Relaxed);
    }
}

static BATCHES: Gauge = Gauge::new();
static CHUNKS: Gauge = Gauge::new();
static TASKS: Gauge = Gauge::new();
static FLOW_RUNS: AtomicU64 = AtomicU64::new(0);
static FLOW_TASKS: AtomicU64 = AtomicU64::new(0);
static FLOW_K: AtomicUsize = AtomicUsize::new(0);
static FLOW_HELPERS: AtomicUsize = AtomicUsize::new(0);
static BATCH_COUNT: AtomicU64 = AtomicU64::new(0);

/// A query chunk being scanned (`None` when the probe is off).
pub(crate) fn enter_chunk() -> Option<InFlight> {
    on().then(|| CHUNKS.enter())
}

/// A subject task of the subject flow being computed (`None` when the probe is off).
pub(crate) fn enter_task() -> Option<InFlight> {
    on().then(|| TASKS.enter())
}

/// One subject flow (`x_subject_flow::run`): `tasks` tasks of `k` subjects, `helpers` helper
/// threads.
pub(crate) fn flow_run(k: usize, tasks: usize, helpers: usize) {
    if on() {
        FLOW_RUNS.fetch_add(1, Relaxed);
        FLOW_TASKS.fetch_add(tasks as u64, Relaxed);
        FLOW_K.store(k, Relaxed);
        FLOW_HELPERS.fetch_max(helpers, Relaxed);
    }
}

/// The stages after the preliminary scan that read a subject's residues.
#[derive(Clone, Copy)]
pub(crate) enum Stage {
    /// The ungapped path's evaluation (`Runtime::subject`, `-ungapped`).
    Link = 0,
    /// The composition redo (speculated in parallel and serial).
    Redo = 1,
    /// The ordinary traceback (`-comp_based_stats 0`).
    Traceback = 2,
}

const STAGE_NAMES: [&str; 3] = ["link", "redo", "traceback"];

/// One query batch's counts (one `search_core` batch range).
pub(crate) struct Batch {
    queries: usize,
    q0: String,
    passes: AtomicUsize,
    chunks: AtomicUsize,
    pool_threads: AtomicUsize,
    prelim: Mutex<HashSet<u32>>,
    prelim_bytes: AtomicU64,
    oids: [Mutex<HashSet<u32>>; 3],
    bytes: [AtomicU64; 3],
    calls: [AtomicU64; 3],
    _in_flight: InFlight,
}

impl Batch {
    /// A batch of `queries` queries whose first is `q0`.
    pub(crate) fn begin(q0: &str, queries: usize) -> Box<Self> {
        Box::new(Self {
            queries,
            q0: q0.split_whitespace().next().unwrap_or("").to_string(),
            passes: AtomicUsize::new(0),
            chunks: AtomicUsize::new(0),
            pool_threads: AtomicUsize::new(0),
            prelim: Mutex::new(HashSet::new()),
            prelim_bytes: AtomicU64::new(0),
            oids: Default::default(),
            bytes: Default::default(),
            calls: Default::default(),
            _in_flight: BATCHES.enter(),
        })
    }

    /// The scan's counts: complete passes over the subjects, chunks, and the pool's threads (0
    /// for a serial search).
    pub(crate) fn scanned(&self, passes: usize, chunks: usize, pool_threads: usize) {
        self.passes.store(passes, Relaxed);
        self.chunks.store(chunks, Relaxed);
        self.pool_threads.store(pool_threads, Relaxed);
    }

    /// Subject `oid` (`length` residues) has preliminary HSPs.
    pub(crate) fn prelim(&self, oid: usize, length: usize) {
        if self.prelim.lock().unwrap().insert(oid as u32) {
            self.prelim_bytes.fetch_add(length as u64, Relaxed);
        }
    }

    /// Stage `stage` reads the residues of subject `oid` (`length` residues).
    pub(crate) fn read(&self, stage: Stage, oid: usize, length: usize) {
        let s = stage as usize;
        self.calls[s].fetch_add(1, Relaxed);
        if self.oids[s].lock().unwrap().insert(oid as u32) {
            self.bytes[s].fetch_add(length as u64, Relaxed);
        }
    }

    /// Writes the batch's line.
    pub(crate) fn end(self: Box<Self>) {
        let index = BATCH_COUNT.fetch_add(1, Relaxed);
        let mut line = format!(
            "XBXSTOREPROBE batch n={index} q0={} queries={} chunks={} passes={} pool_threads={} prelim_oids={} prelim_bytes={}",
            self.q0,
            self.queries,
            self.chunks.load(Relaxed),
            self.passes.load(Relaxed),
            self.pool_threads.load(Relaxed),
            self.prelim.lock().unwrap().len(),
            self.prelim_bytes.load(Relaxed),
        );
        let prelim = self.prelim.lock().unwrap();
        for (s, name) in STAGE_NAMES.iter().enumerate() {
            let oids = self.oids[s].lock().unwrap();
            let outside = oids.iter().filter(|oid| !prelim.contains(oid)).count();
            line.push_str(&format!(
                " {name}_oids={} {name}_bytes={} {name}_calls={} {name}_outside_prelim={outside}",
                oids.len(),
                self.bytes[s].load(Relaxed),
                self.calls[s].load(Relaxed),
            ));
        }
        line.push_str(&format!(" rss_anon_kib={}", rss_anon_kib()));
        eprintln!("{line}");
    }
}

/// The report stage of one query batch (`ReportedSubjects::for_batch`): `oids` reported
/// subjects whose residues (`bytes`) the report records copy from the store.
pub(crate) fn report(q0: &str, queries: usize, oids: usize, bytes: u64) {
    eprintln!(
        "XBXSTOREPROBE report q0={} queries={queries} report_oids={oids} report_bytes={bytes} rss_anon_kib={}",
        q0.split_whitespace().next().unwrap_or(""),
        rss_anon_kib()
    );
}

/// The run's concurrency line.
pub(crate) fn run_summary(threads: usize) {
    eprintln!(
        "XBXSTOREPROBE run threads={threads} batches={} max_batches_in_flight={} max_chunks_in_flight={} max_subject_tasks_in_flight={} flow_runs={} flow_tasks={} flow_k={} flow_helpers={} rss_anon_kib={}",
        BATCH_COUNT.load(Relaxed),
        BATCHES.max(),
        CHUNKS.max(),
        TASKS.max(),
        FLOW_RUNS.load(Relaxed),
        FLOW_TASKS.load(Relaxed),
        FLOW_K.load(Relaxed),
        FLOW_HELPERS.load(Relaxed),
        rss_anon_kib()
    );
}
