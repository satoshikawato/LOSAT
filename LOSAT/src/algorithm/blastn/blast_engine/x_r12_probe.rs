//! EXPERIMENT (exp/r12-fold, v4-R R12; evidence only, never merged): probes of the two-hit fold of
//! dc-megablast (`x_pair_seed::disc::fold_two_hit`), chosen once from `R12_PROBE` (comma list):
//! - `cnt`: the fold runs as a copy that counts its work (seeds, the `last_hit` skips, the cells
//!   each retrieve and each insert walks, the insert's outcome, the window test, the ungapped
//!   extensions and the saved hits, the stage tables' occupancy) and prints one `[R12]` line at
//!   the end. Same decisions and the same table: the counting walks are `retrieve` / `insert`
//!   with a counter.
//! - `copy`: the fold runs as the same copy without the counters (the control of the arms that
//!   change the copy).
//! - `tow` (with `copy`): when `word_length == lut_word_length` (dc-megablast) the copy does not call
//!   `type_of_word`, whose first statement returns `(1, 0, q_off, s_off)` then.
//! - `pfd` (with `copy`): before seed `k`, the copy prefetches the head cell of the chain of seed
//!   `k + 16`'s bucket (read only).
//! - `lv2` (with `copy` or `cnt`; needs `LOSAT_X_DIAGLV`'s level bounds): a side table per fold table
//!   and stage (`Side`) keeps, per bucket, the latest insert's `(diag, level, hit_saved)` and an upper
//!   bound `other_max` of the levels of all the bucket's other cells. When the DIAGLV bound `b` of the
//!   seed is above `other_max`, every cell but the latest insert's lies below `b` and makes the
//!   decisions of no cell (DIAGLV's (ii)); the latest insert's cell is the first cell of its diagonal
//!   in the chain (an insert stops at the first cell of its diagonal or the first stale cell and
//!   writes its diagonal there, and no insert of this bucket came after it), so the retrieve's walk
//!   would return it for that diagonal. The copy then answers without the walk: that cell's values
//!   for its diagonal, no cell for any other. At a stage's start (and for a table not seen in this
//!   stage) `other_max` is the bucket's DIAGLV bound `max_level` and no latest insert is known;
//!   every insert of the copy moves the previous latest level into `other_max`. Same decisions, same
//!   table (the inserts are unchanged).
//! - `insi` (with `copy`): the copy's inserts run an `#[inline(always)]` copy of
//!   `DiagHashTable::insert` (`x_r12_insert_inl`, the same statements) instead of the call.
//!
//! Unset or empty = the dev code paths.
use super::{DiagHashTable, DIAGHASH_NUM_BUCKETS};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

#[derive(Clone, Copy, Debug, Default)]
pub(super) struct Probe {
    pub(super) cnt: bool,
    pub(super) copy: bool,
    pub(super) tow: bool,
    pub(super) pfd: bool,
    pub(super) lv2: bool,
    pub(super) insi: bool,
}

impl Probe {
    #[inline(always)]
    pub(super) fn on(&self) -> bool {
        self.cnt || self.copy
    }

    /// The variant bits of the copy: 1 `tow`, 2 `pfd`, 4 `lv2`, 8 `insi`.
    pub(super) fn var(&self) -> u8 {
        u8::from(self.tow)
            | (u8::from(self.pfd) << 1)
            | (u8::from(self.lv2) << 2)
            | (u8::from(self.insi) << 3)
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
                    "tow" => p.tow = true,
                    "pfd" => p.pfd = true,
                    "lv2" => p.lv2 = true,
                    "insi" => p.insi = true,
                    _ => panic!("R12_PROBE: unknown item {t}"),
                }
            }
            assert!(
                !(p.cnt && (p.tow || p.pfd)),
                "R12_PROBE: cnt runs the plain copy (no tow / pfd)"
            );
            assert!(
                !(p.cnt && p.insi),
                "R12_PROBE: cnt runs the plain copy (no insi)"
            );
            if (p.tow || p.pfd || p.lv2 || p.insi) && !p.cnt {
                p.copy = true;
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
    pub(super) lv2_hit: u64,
    pub(super) lv2_found: u64,
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
            self.lv2_hit,
            self.lv2_found,
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
        "lv2_hit",
        "lv2_found",
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

/// Seeds ahead of the `pfd` prefetch.
pub(super) const PFD: usize = 16;

impl DiagHashTable {
    /// `pfd`: prefetch the head cell of the chain of `diag`'s bucket (reads the backbone only).
    #[inline(always)]
    pub(super) fn x_r12_prefetch(&self, diag: i32) {
        let bucket =
            ((diag as u32).wrapping_mul(0x9E370001) % DIAGHASH_NUM_BUCKETS as u32) as usize;
        let index = self.backbone[bucket] as usize;
        if index != 0 && index < self.chain.len() {
            super::x_r4_probe::prefetch(self.chain.as_ptr().wrapping_add(index));
        }
    }
}

/// `lv2`: the side table of one fold table in one stage (see the module comment).
pub(super) struct Side {
    stage: u64,
    on: bool,
    last_diag: Vec<i32>,
    last_level: Vec<i32>,
    last_saved: Vec<bool>,
    last_valid: Vec<bool>,
    other_max: Vec<i32>,
}

static SIDES: std::sync::Mutex<Vec<(usize, Box<Side>)>> = std::sync::Mutex::new(Vec::new());
static STAGE_ID: AtomicU64 = AtomicU64::new(1);

/// A new stage id (one per call of the seed stage).
pub(super) fn next_stage() -> u64 {
    STAGE_ID.fetch_add(1, Ordering::Relaxed)
}

/// The side table of `table` in stage `stage` (reset when the stage differs). The caller holds the
/// table exclusively (the partition's lock, or the owner's table), so the side table is its own.
pub(super) fn side_for(table: &DiagHashTable, stage: u64) -> *mut Side {
    let key = table as *const DiagHashTable as usize;
    let mut v = SIDES
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    let i = match v.iter().position(|(k, _)| *k == key) {
        Some(i) => i,
        None => {
            v.push((
                key,
                Box::new(Side {
                    stage: 0,
                    on: false,
                    last_diag: vec![0; DIAGHASH_NUM_BUCKETS],
                    last_level: vec![0; DIAGHASH_NUM_BUCKETS],
                    last_saved: vec![false; DIAGHASH_NUM_BUCKETS],
                    last_valid: vec![false; DIAGHASH_NUM_BUCKETS],
                    other_max: vec![i32::MAX; DIAGHASH_NUM_BUCKETS],
                }),
            ));
            v.len() - 1
        }
    };
    let side = &mut v[i].1;
    if side.stage != stage {
        side.stage = stage;
        side.on = table.max_level.len() == DIAGHASH_NUM_BUCKETS;
        side.last_valid.fill(false);
        if side.on {
            side.other_max.copy_from_slice(&table.max_level);
        } else {
            side.other_max.fill(i32::MAX);
        }
    }
    &mut **side as *mut Side
}

impl Side {
    /// The fold's retrieve with the side table: DIAGLV's skip, then the memo of the latest insert,
    /// else the walk.
    #[inline(always)]
    pub(super) fn retrieve(
        &self,
        table: &DiagHashTable,
        diag: i32,
        seeds: &super::x_diag_lv::LvSeeds,
        s_off_pos: usize,
        c: &mut Counts,
        cnt: bool,
    ) -> Option<(i32, i32, bool)> {
        let b = ((diag as u32).wrapping_mul(0x9E370001) % DIAGHASH_NUM_BUCKETS as u32) as usize;
        let bound = seeds.bound(s_off_pos);
        if !self.on {
            return if cnt {
                table.x_r12_retrieve(diag, c)
            } else {
                table.retrieve(diag)
            };
        }
        if table.max_level[b] < bound {
            if cnt {
                c.lv_skip += 1;
            }
            return None;
        }
        if self.last_valid[b] && self.other_max[b] < bound {
            let found = self.last_diag[b] == diag;
            if cnt {
                c.lv2_hit += 1;
                if found {
                    c.lv2_found += 1;
                }
                // `cnt` checks the memo against the walk: the same cell when found, else no cell or a
                // cell below the bound.
                match table.retrieve(diag) {
                    Some((level, _, saved)) if found => assert!(
                        level == self.last_level[b] && saved == self.last_saved[b],
                        "R12 lv2: memo differs from the walk"
                    ),
                    Some((level, _, _)) => {
                        assert!(level < bound, "R12 lv2: a live cell the memo misses")
                    }
                    None => assert!(!found, "R12 lv2: memo found a cell the walk misses"),
                }
            }
            return found.then(|| (self.last_level[b], 0, self.last_saved[b]));
        }
        if cnt {
            table.x_r12_retrieve(diag, c)
        } else {
            table.retrieve(diag)
        }
    }

    /// Every insert of the copy: `(diag, level, hit_saved)` becomes the bucket's latest insert.
    #[inline(always)]
    pub(super) fn note(&mut self, diag: i32, level: i32, hit_saved: bool) {
        let b = ((diag as u32).wrapping_mul(0x9E370001) % DIAGHASH_NUM_BUCKETS as u32) as usize;
        if self.last_valid[b] && self.other_max[b] < self.last_level[b] {
            self.other_max[b] = self.last_level[b];
        }
        self.last_diag[b] = diag;
        self.last_level[b] = level;
        self.last_saved[b] = hit_saved;
        self.last_valid[b] = true;
    }
}

impl DiagHashTable {
    /// `insi`: `insert`'s statements, inlined into the copy.
    #[allow(clippy::too_many_arguments)]
    #[inline(always)]
    pub(super) fn x_r12_insert_inl(
        &mut self,
        diag: i32,
        level: i32,
        hit_len: i32,
        hit_saved: bool,
        s_off: i32,
        window_size: i32,
    ) -> bool {
        let bucket =
            ((diag as u32).wrapping_mul(0x9E370001) % DIAGHASH_NUM_BUCKETS as u32) as usize;
        if let Some(max_level) = self.max_level.get_mut(bucket) {
            if *max_level < level {
                *max_level = level;
            }
        }
        let mut index = self.backbone[bucket];
        while index != 0 {
            let cell = &mut self.chain[index as usize];
            if cell.diag == diag {
                cell.level = level;
                cell.hit_len = hit_len;
                cell.hit_saved = hit_saved;
                return true;
            }
            if s_off - cell.level > window_size {
                cell.diag = diag;
                cell.level = level;
                cell.hit_len = hit_len;
                cell.hit_saved = hit_saved;
                return true;
            }
            index = cell.next;
        }
        if self.occupancy == self.capacity {
            self.capacity = self.capacity.saturating_mul(2);
            self.chain
                .resize(self.capacity as usize, super::DiagHashCell::default());
        }
        let cell_index = self.occupancy;
        let cell = &mut self.chain[cell_index as usize];
        cell.diag = diag;
        cell.level = level;
        cell.hit_len = hit_len;
        cell.hit_saved = hit_saved;
        cell.next = self.backbone[bucket];
        self.backbone[bucket] = cell_index;
        self.occupancy = self.occupancy.saturating_add(1);
        true
    }
}
