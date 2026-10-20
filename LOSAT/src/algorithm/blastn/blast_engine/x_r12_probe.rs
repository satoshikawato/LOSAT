//! EXPERIMENT (exp/r12-fold, v4-R R12; evidence only, never merged): probes of the two-hit fold of
//! dc-megablast (`x_pair_seed::disc::fold_two_hit`), chosen once from `R12_PROBE` (comma list):
//! - `cnt`: the fold runs as a copy that counts its work (seeds, the `last_hit` skips, the cells
//!   each retrieve and each insert walks, the insert's outcome, the window test, the ungapped
//!   extensions and the saved hits, the stage tables' occupancy) and prints one `[R12]` line at
//!   the end. Same decisions and the same table: the counting walks are `retrieve` / `insert`
//!   with a counter.
//! - `copy`: the fold runs as the same copy without the counters (the control of the arms that
//!   change the copy).
//!
//! Unset or empty = the dev code paths.
use super::{DiagHashTable, DIAGHASH_NUM_BUCKETS};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

#[derive(Clone, Copy, Debug, Default)]
pub(super) struct Probe {
    pub(super) cnt: bool,
    pub(super) copy: bool,
}

impl Probe {
    #[inline(always)]
    pub(super) fn on(&self) -> bool {
        self.cnt || self.copy
    }
}

pub(super) fn probe() -> Probe {
    static P: OnceLock<Probe> = OnceLock::new();
    *P.get_or_init(|| {
        let mut p = Probe::default();
        if let Ok(v) = std::env::var("R12_PROBE") {
            for t in v.split(',').map(str::trim) {
                match t {
                    "" => {}
                    "cnt" => p.cnt = true,
                    "copy" => p.copy = true,
                    _ => panic!("R12_PROBE: unknown item {t}"),
                }
            }
            eprintln!("[R12_PROBE] {p:?}");
        }
        p
    })
}

/// Walk-length histogram bins: 0, 1, ..., 15, 16+.
pub(super) const BINS: usize = 17;

/// A fold call's counters (summed into the process counters at the end of the call).
#[derive(Default)]
pub(super) struct Counts {
    pub(super) seeds: u64,
    pub(super) skipped: u64,
    pub(super) lv_skip: u64,
    pub(super) ret_found: u64,
    pub(super) ret_cells: u64,
    pub(super) ret_hist: [u64; BINS],
    pub(super) first: u64,
    pub(super) second: u64,
    pub(super) ext: u64,
    pub(super) saved: u64,
    pub(super) ins: u64,
    pub(super) ins_cells: u64,
    pub(super) ins_hist: [u64; BINS],
    pub(super) ins_match: u64,
    pub(super) ins_take: u64,
    pub(super) ins_append: u64,
    pub(super) ctx_change: u64,
}

const N: usize = 17 + 2 * BINS;
static TOTAL: [AtomicU64; N] = [const { AtomicU64::new(0) }; N];
static STAGES: AtomicU64 = AtomicU64::new(0);
static OCC_SUM: AtomicU64 = AtomicU64::new(0);
static OCC_MAX: AtomicU64 = AtomicU64::new(0);
static CAP_MAX: AtomicU64 = AtomicU64::new(0);

impl Counts {
    pub(super) fn flush(&self) {
        let v = [
            self.seeds,
            self.skipped,
            self.lv_skip,
            self.ret_found,
            self.ret_cells,
            self.first,
            self.second,
            self.ext,
            self.saved,
            self.ins,
            self.ins_cells,
            self.ins_match,
            self.ins_take,
            self.ins_append,
            self.ctx_change,
            0,
            0,
        ];
        for (i, x) in v.iter().enumerate() {
            if *x != 0 {
                TOTAL[i].fetch_add(*x, Ordering::Relaxed);
            }
        }
        for b in 0..BINS {
            if self.ret_hist[b] != 0 {
                TOTAL[17 + b].fetch_add(self.ret_hist[b], Ordering::Relaxed);
            }
            if self.ins_hist[b] != 0 {
                TOTAL[17 + BINS + b].fetch_add(self.ins_hist[b], Ordering::Relaxed);
            }
        }
    }
}

/// The stage table at the end of a stage (after the merge).
pub(super) fn stage_end(hash: &DiagHashTable) {
    STAGES.fetch_add(1, Ordering::Relaxed);
    OCC_SUM.fetch_add(u64::from(hash.occupancy), Ordering::Relaxed);
    OCC_MAX.fetch_max(u64::from(hash.occupancy), Ordering::Relaxed);
    CAP_MAX.fetch_max(u64::from(hash.capacity), Ordering::Relaxed);
}

pub(super) fn report() {
    let p = probe();
    if !p.cnt {
        return;
    }
    let t: Vec<u64> = TOTAL.iter().map(|a| a.load(Ordering::Relaxed)).collect();
    let names = [
        "seeds",
        "skipped",
        "lv_skip",
        "ret_found",
        "ret_cells",
        "first",
        "second",
        "ext",
        "saved",
        "ins",
        "ins_cells",
        "ins_match",
        "ins_take",
        "ins_append",
        "ctx_change",
    ];
    let mut s = String::new();
    for (i, n) in names.iter().enumerate() {
        s.push_str(&format!(" {n}={}", t[i]));
    }
    let hist = |o: usize| {
        (0..BINS)
            .map(|b| t[o + b].to_string())
            .collect::<Vec<_>>()
            .join(",")
    };
    eprintln!(
        "[R12] cnt{s} ret_hist={} ins_hist={} stages={} occ_sum={} occ_max={} cap_max={} cell_bytes={} buckets={}",
        hist(17),
        hist(17 + BINS),
        STAGES.load(Ordering::Relaxed),
        OCC_SUM.load(Ordering::Relaxed),
        OCC_MAX.load(Ordering::Relaxed),
        CAP_MAX.load(Ordering::Relaxed),
        std::mem::size_of::<super::DiagHashCell>(),
        DIAGHASH_NUM_BUCKETS
    );
}

impl DiagHashTable {
    /// `retrieve` that also counts the cells it walks.
    #[inline(always)]
    pub(super) fn x_r12_retrieve(&self, diag: i32, c: &mut Counts) -> Option<(i32, i32, bool)> {
        let bucket =
            ((diag as u32).wrapping_mul(0x9E370001) % DIAGHASH_NUM_BUCKETS as u32) as usize;
        let mut index = self.backbone[bucket];
        let mut k = 0usize;
        let mut r = None;
        while index != 0 {
            let cell = &self.chain[index as usize];
            k += 1;
            if cell.diag == diag {
                r = Some((cell.level, cell.hit_len, cell.hit_saved));
                break;
            }
            index = cell.next;
        }
        c.ret_cells += k as u64;
        c.ret_hist[k.min(BINS - 1)] += 1;
        if r.is_some() {
            c.ret_found += 1;
        }
        r
    }

    /// The walk of `insert` (read only): counts the cells and the outcome the insert will have.
    #[inline(always)]
    pub(super) fn x_r12_insert_count(
        &self,
        diag: i32,
        s_off: i32,
        window_size: i32,
        c: &mut Counts,
    ) {
        let bucket =
            ((diag as u32).wrapping_mul(0x9E370001) % DIAGHASH_NUM_BUCKETS as u32) as usize;
        let mut index = self.backbone[bucket];
        let mut k = 0usize;
        let mut outcome = 2u8;
        while index != 0 {
            let cell = &self.chain[index as usize];
            k += 1;
            if cell.diag == diag {
                outcome = 0;
                break;
            }
            if s_off - cell.level > window_size {
                outcome = 1;
                break;
            }
            index = cell.next;
        }
        c.ins += 1;
        c.ins_cells += k as u64;
        c.ins_hist[k.min(BINS - 1)] += 1;
        match outcome {
            0 => c.ins_match += 1,
            1 => c.ins_take += 1,
            _ => c.ins_append += 1,
        }
    }
}

impl DiagHashTable {
    /// The fold's retrieve (`retrieve_below` when the table has level bounds, else `retrieve`,
    /// without the shadow) with the counters: a skipped walk counts in `lv_skip`.
    #[inline(always)]
    pub(super) fn x_r12_ret(
        &self,
        diag: i32,
        x_lv: bool,
        seeds: &super::x_diag_lv::LvSeeds,
        s_off_pos: usize,
        c: &mut Counts,
    ) -> Option<(i32, i32, bool)> {
        if x_lv {
            let bucket =
                ((diag as u32).wrapping_mul(0x9E370001) % DIAGHASH_NUM_BUCKETS as u32) as usize;
            let bound = seeds.bound(s_off_pos);
            if matches!(self.max_level.get(bucket), Some(&m) if m < bound) {
                c.lv_skip += 1;
                return None;
            }
        }
        self.x_r12_retrieve(diag, c)
    }
}
