//! EXPERIMENT (LOSAT_X_PAIRPARFOLD; a strict-set switch since 2026-10-17): the
//! diagonal-hash fold of the strips paths of PAIRPAR (blastn's one-hit fold, the
//! parent's `fold`) and `LOSAT_X_PAIRPARDISC` (dc-megablast's two-hit fold,
//! `disc::fold_two_hit`) split by diagonal-hash bucket over the threads of the search pool,
//! instead of done by the owner thread alone. No NCBI counterpart (scheduling only).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:361-381
//! (`s_BlastDiagHashRetrieve`)
//! ```c
//!     Uint4 bucket = ((Uint4) diag * 0x9E370001) % DIAGHASH_NUM_BUCKETS;
//!     Uint4 index = table->backbone[bucket];
//!
//!     while (index) {
//!         if (table->chain[index].diag == diag) {
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:396-449
//! (`s_BlastDiagHashInsert`)
//! ```c
//!     Uint4 bucket = ((Uint4) diag * 0x9E370001) % DIAGHASH_NUM_BUCKETS;
//!     Uint4 index = table->backbone[bucket];
//! ...
//!         if (table->chain[index].diag == diag) {
//! ...
//!             if (s_off - table->chain[index].level > window_size) {
//! ...
//!     if (table->occupancy == table->capacity) {
//!         table->capacity *= 2;
//!         table->chain =
//!             realloc(table->chain, table->capacity * sizeof(DiagHashCell));
//! ...
//!     cell = table->chain + table->occupancy;
//! ...
//!     cell->next = table->backbone[bucket];
//!     table->backbone[bucket] = table->occupancy;
//!     table->occupancy++;
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:825-841,853-858,939-941
//! (`s_BlastnDiagHashExtendInitialHit`)
//! ```c
//!     Int4 Delta = MIN(word_params->options->scan_range, window_size - word_length);
//! ...
//!     rc = s_BlastDiagHashRetrieve(hash_table, diag, &last_hit, &s_l, &hit_saved);
//! ...
//!     if (s_off_pos < last_hit) return 0;
//!
//!     if (two_hits && (hit_saved || s_end_pos > last_hit + window_size )) {
//! ...
//!         if (word_type == 1) {
//!             /* try off-diagonal */
//! ...
//!             if (Delta < 0) Delta = 0;
//!             for (delta = 1; delta <= Delta; ++delta) {
//! ...
//!     s_BlastDiagHashInsert(hash_table, diag, s_end_pos,
//!                           (hit_ready) ? 0 : s_end_pos - s_off_pos,
//!                           hit_ready, s_off_pos, window_size + Delta + 1);
//! ```
//!
//! How the stage runs (`seed_stage_part`): the participants are the owner (the thread that
//! searches the chunk) and `threads - 1` jobs of an `in_place_scope`. Each participant claims
//! the next strip of the scan (`plan_strips`; the parent's `scan_strip` or
//! `disc::scan_strip_disc`) when it lies within `ahead` strips of the slowest partition, cuts
//! the strip's seeds by partition (`bucket % parts`, `bucket` NCBI's hash of the seed's
//! diagonal; a stable counting sort, so each partition keeps the seeds in scan order) and folds
//! partition p of the next strip of any partition p it can lock (`try_lock`), into a table of
//! that partition (`PartState::hash`). The owner alone can scan every strip and fold every
//! partition, so no participant waits for a job that may never start. A helper job that starts
//! on a participant's thread (a thread that waits inside a nested parallel call of its own scan
//! or fold) returns at once; a helper on any other thread that finds nothing to do leaves once
//! the whole stage has made no progress (no strip scanned, no partition folded) for
//! `HELPER_STALL` (the owner never leaves). A thread that waits inside a piece of a
//! participant's nested call that it stole may run a helper job above the work it interrupted;
//! that helper waits for the stage while the stage waits for that work (the pattern that hung
//! LOSAT_X_BNSTREAM, #218; test `a_nested_parallel_call_in_every_participant_does_not_hang`)
//! until the stall exit releases it. At the end the hits are put in the order of their seeds'
//! keys (strip, position in the strip's scan order) and the changed buckets of every partition
//! are written into the stage's table with the reference's cell numbers (`merge`).
//!
//! Every stage copies the stage's table into every partition's table (`rebuild`), so the
//! partition tables never carry anything from one stage to the next; only their allocations
//! (and those of the strip and merge buffers) stay in the subject scratch (`PartScratch`, a
//! field of `run.rs` `SubjectScratch`).
//!
//! Why the output is the same:
//! - A seed's fold reads and writes the table only through its own diagonal's bucket: the
//!   retrieve walks `backbone[bucket]`'s chain, the insert walks the same chain and changes a
//!   cell of it (the cell of the diagonal, or a stale cell it takes over) or appends a cell and
//!   links it at that bucket's head, and `LOSAT_X_DIAGLV` reads and raises `max_level[bucket]`.
//!   The append also reads and writes `occupancy` and `capacity`, the only state shared by the
//!   buckets (below). `s_TypeOfWord` (the lookup table, the query masks, the subject) and the
//!   ungapped extension (the query, the subject, the matrix) read only the inputs, and with
//!   `scan_range == 0` (in the scope of both strips paths) `Delta` is 0 and the off-diagonal
//!   loop, the only reader of other diagonals, runs no iteration. The hit list is only appended
//!   to.
//! - So the chain a seed reads is changed only by seeds of its bucket, all of which are in its
//!   partition, and a partition folds its seeds in the reference order (the strips in order, a
//!   strip's seeds in the scan order). Since partition p's table starts with the stage's chains
//!   of p's buckets, every seed makes the reference's decisions (the `last_hit` test, the
//!   window test, `s_TypeOfWord`, the extension and the cutoff) and the reference's changes to
//!   its chain, by induction over the reference order.
//! - The reference numbers an appended cell `occupancy` at the time of the append, so the
//!   appends of the stage get the numbers occ0, occ0 + 1, ... in the reference order of the
//!   appending seeds, which is the order of their keys (one insert per seed). The merge gives
//!   the append of rank r among all partitions' appends the number occ0 + r and maps the
//!   `next` links and the heads that point to appended cells; a cell's contents and its place in
//!   its chain are those of its partition. The capacity: the reference doubles it before an
//!   append when `occupancy == capacity`, so from c0 at the stage's start it ends at the
//!   smallest c0 * 2^k that holds the final occupancy, which the merge sets. The bounds of
//!   `LOSAT_X_DIAGLV` are per bucket. The table is therefore the reference's cell for cell.
//! - The reference pushes a seed's hit when the seed is folded and each seed pushes at most
//!   one, so the hits in key order are the reference's list element for element.
//! - No floating-point value is computed in the fold, so no order of floating-point
//!   operations changes.

use super::super::{DiagHashCell, DiagHashTable, UngappedHit, DIAGHASH_NUM_BUCKETS};
use super::{FoldCounts, Seed};
use std::cell::Cell;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, OnceLock, PoisonError};

/// No NCBI counterpart: reads LOSAT_X_PAIRPARFOLD once (a strict-set switch: on unless it is `0`).
/// It does not change any value NCBI computes.
fn switch_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_PAIRPARFOLD"))
}

/// No NCBI counterpart: the number of partitions of a stage on `threads` threads (one per
/// thread, at most one per bucket), or `None` when the switch is off or the pool has one
/// thread (the owner folds).
pub(super) fn partitions(threads: usize) -> Option<usize> {
    if switch_on() {
        // A stage on a pool has two threads or more; one thread only in the PAIRPARDISCSHADOW and
        // PAIRPARNASHADOW runs.
        crate::utils::xreach::x_reach!("LOSAT_X_PAIRPARFOLD", threads >= 2);
    }
    (switch_on() && threads >= 2).then(|| threads.min(DIAGHASH_NUM_BUCKETS))
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:368
/// ```c
///     Uint4 bucket = ((Uint4) diag * 0x9E370001) % DIAGHASH_NUM_BUCKETS;
/// ```
#[inline(always)]
fn bucket_of_diag(diag: i32) -> usize {
    ((diag as u32).wrapping_mul(0x9E370001) % DIAGHASH_NUM_BUCKETS as u32) as usize
}

/// The bucket of a seed's diagonal, with the diagonal as the folds compute it (`diag = s_off -
/// q_off`, as `i32`).
#[inline(always)]
fn bucket_of(seed: &Seed) -> usize {
    bucket_of_diag((seed.s as isize - seed.q as isize) as i32)
}

/// What a fold call of one partition writes besides the table and the hit list: the keys of
/// its hits and of its inserts that appended a cell, the buckets its inserts changed, and its
/// time in ungapped extensions (LOSAT_TIMING, flushed once per stage). A fold without
/// partitions (`PT` false) writes nothing here.
#[derive(Default)]
pub(super) struct PartOut {
    /// The strip being folded (the high half of the keys).
    pub(super) strip: u64,
    pub(super) hit_keys: Vec<u64>,
    appends: Vec<(u64, u32)>,
    touched: Vec<bool>,
    dirty: Vec<u16>,
    pub(super) ungapped_ns: u64,
    pub(super) ungapped_calls: u64,
}

impl PartOut {
    /// The key of seed `k` of a fold call: (strip, the seed's position in the strip's scan
    /// order).
    #[inline(always)]
    pub(super) fn key(&self, pos: &[u32], k: usize) -> u64 {
        (self.strip << 32) | u64::from(pos[k])
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:939-941
/// ```c
///     s_BlastDiagHashInsert(hash_table, diag, s_end_pos,
///                           (hit_ready) ? 0 : s_end_pos - s_off_pos,
///                           hit_ready, s_off_pos, window_size + Delta + 1);
/// ```
/// A fold's insert (`DiagHashTable::insert_memo`; with `LEAN`, LOSAT_X_FOLDLEAN, the same statements
/// inlined, `insert_inl`). Without `MEMO` the table has no memo (LOSAT_X_DIAGLV2) and the insert
/// does not test for one: `insert_plain`, the statements before that switch (with `LEAN` inlined,
/// `insert_inl::<false>`). With partitions (`PT`) it also records
/// the diagonal's bucket and, when the insert appended a cell, the seed's key with the cell's
/// number in the partition's table.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
pub(super) fn insert<const PT: bool, const LEAN: bool, const MEMO: bool>(
    hash: &mut DiagHashTable,
    out: &mut PartOut,
    key: u64,
    diag: i32,
    level: i32,
    hit_len: i32,
    hit_saved: bool,
    s_off: i32,
    window_size: i32,
) {
    if PT {
        let occupancy = hash.occupancy;
        if LEAN {
            hash.insert_inl::<MEMO>(diag, level, hit_len, hit_saved, s_off, window_size);
        } else if MEMO {
            hash.insert_memo(diag, level, hit_len, hit_saved, s_off, window_size);
        } else {
            hash.insert_plain(diag, level, hit_len, hit_saved, s_off, window_size);
        }
        let bucket = bucket_of_diag(diag);
        if !out.touched[bucket] {
            out.touched[bucket] = true;
            out.dirty.push(bucket as u16);
        }
        if hash.occupancy != occupancy {
            out.appends.push((key, occupancy));
        }
    } else if LEAN {
        hash.insert_inl::<MEMO>(diag, level, hit_len, hit_saved, s_off, window_size);
    } else if MEMO {
        hash.insert_memo(diag, level, hit_len, hit_saved, s_off, window_size);
    } else {
        hash.insert_plain(diag, level, hit_len, hit_saved, s_off, window_size);
    }
}

/// No NCBI counterpart (scheduling only): a partition's state during a stage, on cache lines of
/// its own. The participants poll every partition's lock (`try_lock`) while one of them folds
/// into it; without the alignment the lock word and the first fields shared a line with the
/// fields the folding thread writes, or with the neighbouring partition's, as the size of
/// `PartState` placed them (DiagHashTable's memo of LOSAT_X_DIAGLV2 made DQ n16 x1.028 slower
/// with both switches off: memory-ordering machine clears x1.37, same instructions).
#[repr(align(128))]
struct OwnLines<T>(T);

/// The fold state of one partition (its allocations kept in the subject scratch between
/// stages; its contents set by `rebuild` at the start of every stage).
pub(super) struct PartState {
    /// The partition's table: a copy of the stage's table at the start of the stage, of which
    /// the partition reads and changes only the chains of its own buckets.
    pub(super) hash: DiagHashTable,
    pub(super) hits: Vec<UngappedHit>,
    pub(super) out: PartOut,
    pub(super) counts: FoldCounts,
    /// the time in ungapped extensions of the folds of this partition (LOSAT_TIMING)
    pub(super) scan_pause_ns: u64,
    fold_ns: u64,
    /// The next strip this partition folds.
    next: usize,
    /// The numbers of the partition's appended cells in the stage's table.
    map: Vec<u32>,
}

impl PartState {
    fn new() -> Self {
        PartState {
            hash: DiagHashTable {
                num_buckets: 0,
                occupancy: 0,
                capacity: 0,
                backbone: Vec::new(),
                chain: Vec::new(),
                offset: 0,
                window: 0,
                max_level: Vec::new(),
                x_lv_shadow: false,
                x_last: Vec::new(),
                #[cfg(debug_assertions)]
                x_last_fresh: true,
            },
            hits: Vec::new(),
            out: PartOut {
                touched: vec![false; DIAGHASH_NUM_BUCKETS],
                ..PartOut::default()
            },
            counts: FoldCounts::default(),
            scan_pause_ns: 0,
            fold_ns: 0,
            next: 0,
            map: Vec::new(),
        }
    }

    /// The per-stage fields back to the start of a stage.
    fn reset(&mut self) {
        self.hits.clear();
        self.out.strip = 0;
        self.out.hit_keys.clear();
        self.out.appends.clear();
        for &bucket in &self.out.dirty {
            self.out.touched[bucket as usize] = false;
        }
        self.out.dirty.clear();
        self.out.ungapped_ns = 0;
        self.out.ungapped_calls = 0;
        self.counts = FoldCounts::default();
        self.scan_pause_ns = 0;
        self.fold_ns = 0;
        self.next = 0;
        self.map.clear();
    }
}

/// The seeds of one strip cut by partition: partition p's seeds are `seeds[offs[p]..offs[p +
/// 1]]`, in scan order, with their positions in the strip's scan order in `pos`.
#[derive(Default)]
struct StripBuf {
    seeds: Vec<Seed>,
    pos: Vec<u32>,
    offs: Vec<u32>,
}

/// No NCBI counterpart (scheduling only): the allocations of the partitioned fold's partition
/// tables and buffers, kept in the subject scratch (`run.rs` `SubjectScratch`) from one stage
/// to the next. No content is carried over: every stage copies the stage's table into the
/// partition tables (`rebuild`).
#[derive(Default)]
pub(in super::super) struct PartScratch {
    parts: Vec<PartState>,
    spare: Vec<StripBuf>,
    raw: Vec<Vec<Seed>>,
    order: Vec<(u64, u32, u32)>,
    apps: Vec<(u64, u32, u32)>,
}

/// How a strips stage folds when the switch is on: the scratch and the partition count.
pub(super) struct FoldPart<'a> {
    pub(super) scratch: &'a mut PartScratch,
    pub(super) parts: usize,
}

/// What one partitioned stage did.
#[derive(Default)]
pub(super) struct PartStats {
    pub(super) counts: FoldCounts,
    /// the owner's time in ungapped extensions (in the partitions it folded): the strips paths
    /// take it off the stage's wall time for `[TIMING] scan_lookup`, which so stays a figure of
    /// the owner's wall time as without the switch
    pub(super) scan_pause_ns: u64,
    /// the time in folds summed over the participants (CPU time, not wall time)
    pub(super) fold_ns: u64,
    pub(super) owner_scan_ns: u64,
    pub(super) wait_ns: u64,
    pub(super) ungapped_ns: u64,
    pub(super) ungapped_calls: u64,
    pub(super) merge_ns: u64,
    /// helpers that left a stalled stage before its end
    pub(super) helpers_left: u64,
}

thread_local! {
    // No NCBI counterpart: whether this thread runs a participant of a partitioned stage.
    static IN_STAGE: Cell<bool> = const { Cell::new(false) };
}

/// No NCBI counterpart: marks this thread as a participant while it lives.
struct InStage(bool);

impl InStage {
    fn enter() -> Self {
        InStage(IN_STAGE.with(|flag| flag.replace(true)))
    }
}

impl Drop for InStage {
    fn drop(&mut self) {
        IN_STAGE.with(|flag| flag.set(self.0));
    }
}

/// No NCBI counterpart: marks the stage as failed if a participant unwinds, so that the owner
/// stops (and panics) and the helpers stop.
struct Failing<'a> {
    failed: &'a AtomicBool,
    done: bool,
}

impl Drop for Failing<'_> {
    fn drop(&mut self) {
        if !self.done {
            self.failed.store(true, Ordering::Release);
        }
    }
}

/// No NCBI counterpart (scheduling only): how long a helper that finds nothing to do waits for
/// the stage to make progress before it leaves. A helper waits while a partition that lags
/// holds back the window of strips; as long as the stage makes progress, the helper stays.
/// The stage stops making progress for this long in three cases:
/// - the work it needs runs below a waiting helper on that helper's own thread: a thread that
///   is not a participant stole a piece of a participant's nested parallel call and then,
///   while waiting inside it, a helper job (on a participant's own thread the guard `IN_STAGE`
///   makes the job return at once). Only the stall exit ends this wait, so the helper has to
///   leave;
/// - the OS deschedules the owner (or the participant that holds the lagging partition) on a
///   loaded machine;
/// - one long fold unit (a partition of a strip with much work) holds back the window.
///
/// In the last two cases the idle helpers leave for good, and the rest of the stage runs on
/// fewer threads. The output does not change (any participant can stop at any time and the
/// owner finishes the stage), but the stage is slower.
const HELPER_STALL: std::time::Duration = std::time::Duration::from_millis(10);

/// No NCBI counterpart: a participant's wait when it found nothing to do. It spins a little,
/// then yields; a native helper that keeps finding nothing sleeps 50 us between looks. On
/// wasm32 it only yields: `std::thread::sleep` panics under the browser's WASI shim.
#[inline]
fn pause(idle: u32, owner: bool) {
    if idle < 64 {
        std::hint::spin_loop();
        return;
    }
    #[cfg(not(target_arch = "wasm32"))]
    if !owner && idle >= 4096 {
        std::thread::sleep(std::time::Duration::from_micros(50));
        return;
    }
    let _ = owner;
    std::thread::yield_now();
}

fn lock<T>(m: &Mutex<T>) -> std::sync::MutexGuard<'_, T> {
    m.lock().unwrap_or_else(PoisonError::into_inner)
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684
/// ```c
///     while(s_DetermineScanningOffsets(subject, word_length, lut_word_length, scan_range)) {
///
///         hitsfound = scansub(lookup_wrap, subject, offset_pairs, max_hits, &scan_range[1]);
/// ...
///         hits_extended += extend(offset_pairs, hitsfound, word_params,
///                                 lookup_wrap, query, subject, matrix,
///                                 query_info, ewp, init_hitlist, scan_range[2] + lut_word_length);
///     }
/// ```
/// The seed stage of `n` strips on `threads` threads (this one and `threads - 1` jobs of an
/// `in_place_scope`) with the fold cut into `parts` partitions (module comment): `scan` makes
/// the seeds of a strip, `fold` folds a partition's seeds of one strip (`PartOut::strip` set)
/// into the partition's state; the hits go to `hits` and the table changes to `hash`, both as
/// the reference fold makes them.
#[allow(clippy::too_many_arguments)]
pub(super) fn seed_stage_part<FS, FF>(
    n: usize,
    threads: usize,
    ahead: usize,
    hash: &mut DiagHashTable,
    hits: &mut Vec<UngappedHit>,
    fold_part: FoldPart<'_>,
    scan: &FS,
    fold: &FF,
    timed: bool,
) -> PartStats
where
    FS: Fn(usize, &mut Vec<Seed>) + Sync,
    FF: Fn(&mut PartState, &[Seed], &[u32]) + Sync,
{
    let FoldPart { scratch, parts } = fold_part;
    assert!(
        (1..=DIAGHASH_NUM_BUCKETS).contains(&parts),
        "LOSAT_X_PAIRPARFOLD: {parts} partitions"
    );
    let ahead = ahead.max(1);
    rebuild(scratch, hash, parts);
    let occ0 = hash.occupancy;

    let states: Vec<Mutex<OwnLines<PartState>>> = scratch
        .parts
        .drain(..)
        .map(|st| Mutex::new(OwnLines(st)))
        .collect();
    let spare = Mutex::new(std::mem::take(&mut scratch.spare));
    let raw_spare = Mutex::new(std::mem::take(&mut scratch.raw));
    let part_next: Vec<AtomicUsize> = (0..parts).map(|_| AtomicUsize::new(0)).collect();
    let slots: Vec<Mutex<Option<Arc<StripBuf>>>> = (0..n).map(|_| Mutex::new(None)).collect();
    let ready: Vec<AtomicBool> = (0..n).map(|_| AtomicBool::new(false)).collect();
    let done: Vec<AtomicUsize> = (0..n).map(|_| AtomicUsize::new(0)).collect();
    let cursor = AtomicUsize::new(0);
    let failed = AtomicBool::new(false);
    // strips scanned plus partitions of strips folded, and the helpers that left a stalled stage
    let progress = AtomicUsize::new(0);
    let left = AtomicUsize::new(0);

    // A strip's seeds, cut by partition with a stable counting sort.
    let produce = |i: usize, raw: &mut Vec<Seed>, fill: &mut Vec<u32>| {
        raw.clear();
        scan(i, raw);
        let mut buf = lock(&spare).pop().unwrap_or_default();
        let StripBuf { seeds, pos, offs } = &mut buf;
        offs.clear();
        offs.resize(parts + 1, 0);
        for seed in raw.iter() {
            offs[bucket_of(seed) % parts + 1] += 1;
        }
        for p in 0..parts {
            offs[p + 1] += offs[p];
        }
        fill.clear();
        fill.extend_from_slice(&offs[..parts]);
        let zero = Seed {
            context: 0,
            q: 0,
            s: 0,
            s_range: 0,
        };
        seeds.clear();
        seeds.resize(raw.len(), zero);
        pos.clear();
        pos.resize(raw.len(), 0);
        for (k, seed) in raw.iter().enumerate() {
            let p = bucket_of(seed) % parts;
            let at = fill[p] as usize;
            fill[p] += 1;
            seeds[at] = *seed;
            pos[at] = k as u32;
        }
        *lock(&slots[i]) = Some(Arc::new(buf));
        ready[i].store(true, Ordering::Release);
    };

    // (the owner's time in scans, in waits and in the ungapped extensions of its folds)
    let participant = |me: usize, owner: bool| -> (u64, u64, u64) {
        // A helper job that this thread steals while it waits inside a participant (a nested
        // parallel call) does not join: it could wait for the partition that participant holds.
        if !owner && IN_STAGE.with(Cell::get) {
            return (0, 0, 0);
        }
        let _in_stage = InStage::enter();
        let mut failing = Failing {
            failed: &failed,
            done: false,
        };
        let mut raw = lock(&raw_spare).pop().unwrap_or_default();
        let mut fill: Vec<u32> = Vec::with_capacity(parts);
        let (mut scan_ns, mut wait_ns, mut pause_ns) = (0u64, 0u64, 0u64);
        let mut wait_start: Option<std::time::Instant> = None;
        let mut idle = 0u32;
        let mut stall: Option<(std::time::Instant, usize)> = None;
        loop {
            if failed.load(Ordering::Acquire) {
                if owner {
                    panic!("LOSAT_X_PAIRPARFOLD: a participant of the partitioned fold failed");
                }
                failing.done = true;
                return (scan_ns, wait_ns, pause_ns);
            }
            let mut did = false;
            for k in 0..parts {
                let p = (me + k) % parts;
                let i = part_next[p].load(Ordering::Acquire);
                if i >= n || !ready[i].load(Ordering::Acquire) {
                    continue;
                }
                let Ok(mut guard) = states[p].try_lock() else {
                    continue;
                };
                let st = &mut guard.0;
                let i = st.next;
                if i >= n || !ready[i].load(Ordering::Acquire) {
                    continue;
                }
                let strip = lock(&slots[i]).clone().expect("a ready strip");
                let (a, b) = (strip.offs[p] as usize, strip.offs[p + 1] as usize);
                if a < b {
                    st.out.strip = i as u64;
                    let t0 = timed.then(std::time::Instant::now);
                    let pause0 = st.scan_pause_ns;
                    fold(st, &strip.seeds[a..b], &strip.pos[a..b]);
                    if let Some(t0) = t0 {
                        st.fold_ns += t0.elapsed().as_nanos() as u64;
                    }
                    if owner {
                        pause_ns += st.scan_pause_ns - pause0;
                    }
                }
                drop(strip);
                st.next = i + 1;
                part_next[p].store(i + 1, Ordering::Release);
                drop(guard);
                progress.fetch_add(1, Ordering::Relaxed);
                if done[i].fetch_add(1, Ordering::AcqRel) + 1 == parts {
                    if let Some(buf) = lock(&slots[i]).take() {
                        if let Ok(buf) = Arc::try_unwrap(buf) {
                            lock(&spare).push(buf);
                        }
                    }
                }
                did = true;
            }
            let min_next = part_next
                .iter()
                .map(|next| next.load(Ordering::Acquire))
                .min()
                .unwrap_or(n);
            if min_next >= n {
                failing.done = true;
                if let Some(t0) = wait_start {
                    wait_ns += t0.elapsed().as_nanos() as u64;
                }
                lock(&raw_spare).push(raw);
                return (scan_ns, wait_ns, pause_ns);
            }
            let c = cursor.load(Ordering::Relaxed);
            if c < n && c < min_next + ahead {
                let i = cursor.fetch_add(1, Ordering::Relaxed);
                if i < n {
                    let t0 = (timed && owner).then(std::time::Instant::now);
                    produce(i, &mut raw, &mut fill);
                    progress.fetch_add(1, Ordering::Relaxed);
                    if let Some(t0) = t0 {
                        scan_ns += t0.elapsed().as_nanos() as u64;
                    }
                    did = true;
                }
            }
            if did {
                idle = 0;
                stall = None;
                if let Some(t0) = wait_start.take() {
                    wait_ns += t0.elapsed().as_nanos() as u64;
                }
            } else {
                if timed && owner && wait_start.is_none() {
                    wait_start = Some(std::time::Instant::now());
                }
                idle = idle.saturating_add(1);
                if !owner && idle.is_multiple_of(64) {
                    // the owner can finish the stage alone
                    let now = progress.load(Ordering::Relaxed);
                    match stall {
                        Some((since, seen)) if seen == now => {
                            if since.elapsed() >= HELPER_STALL {
                                left.fetch_add(1, Ordering::Relaxed);
                                failing.done = true;
                                lock(&raw_spare).push(raw);
                                return (scan_ns, wait_ns, pause_ns);
                            }
                        }
                        _ => stall = Some((std::time::Instant::now(), now)),
                    }
                }
                pause(idle, owner);
            }
        }
    };

    let (owner_scan_ns, wait_ns, owner_pause_ns) = run_participants(threads, &participant);

    let mut states: Vec<PartState> = states
        .into_iter()
        .map(|m| m.into_inner().unwrap_or_else(PoisonError::into_inner).0)
        .collect();
    let t_merge = timed.then(std::time::Instant::now);
    let mut stats = PartStats {
        scan_pause_ns: owner_pause_ns,
        owner_scan_ns,
        wait_ns,
        helpers_left: left.into_inner() as u64,
        ..PartStats::default()
    };
    for st in &states {
        stats.counts.seeds += st.counts.seeds;
        stats.counts.skipped += st.counts.skipped;
        stats.fold_ns += st.fold_ns;
        stats.ungapped_ns += st.out.ungapped_ns;
        stats.ungapped_calls += st.out.ungapped_calls;
    }
    merge(
        hash,
        hits,
        &mut states,
        occ0,
        &mut scratch.order,
        &mut scratch.apps,
    );
    if let Some(t0) = t_merge {
        stats.merge_ns = t0.elapsed().as_nanos() as u64;
    }
    // (the allocations only: the next stage's `rebuild` sets every partition's contents)
    scratch.parts = states;
    scratch.spare = spare.into_inner().unwrap_or_else(PoisonError::into_inner);
    scratch.raw = raw_spare
        .into_inner()
        .unwrap_or_else(PoisonError::into_inner);
    stats
}

/// No NCBI counterpart (scheduling only): `participant(0, true)` on this thread (the owner)
/// and `participant(h, false)` for `h` in `1..threads` as jobs of an `in_place_scope` of its
/// pool; outside a pool, only the owner. Returns the owner's result.
fn run_participants(
    threads: usize,
    participant: &(impl Fn(usize, bool) -> (u64, u64, u64) + Sync),
) -> (u64, u64, u64) {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if threads > 1 && rayon::current_thread_index().is_some() {
        return rayon::in_place_scope(|scope| {
            for h in 1..threads {
                scope.spawn(move |_| {
                    participant(h, false);
                });
            }
            participant(0, true)
        });
    }
    let _ = threads;
    participant(0, true)
}

/// No NCBI counterpart: every partition's table made a copy of the stage's table, and its
/// per-stage fields reset (at the start of every stage; only the allocations are reused).
fn rebuild(scratch: &mut PartScratch, hash: &DiagHashTable, parts: usize) {
    scratch.parts.truncate(parts);
    while scratch.parts.len() < parts {
        scratch.parts.push(PartState::new());
    }
    for st in &mut scratch.parts {
        st.reset();
        let t = &mut st.hash;
        t.num_buckets = hash.num_buckets;
        t.occupancy = hash.occupancy;
        t.backbone.clone_from(&hash.backbone);
        t.chain.clear();
        t.chain
            .extend_from_slice(&hash.chain[..hash.occupancy as usize]);
        t.chain
            .resize(hash.capacity as usize, DiagHashCell::default());
        t.capacity = hash.capacity;
        t.offset = hash.offset;
        t.window = hash.window;
        t.max_level.clone_from(&hash.max_level);
        t.x_lv_shadow = hash.x_lv_shadow;
        t.x_last.clone_from(&hash.x_last);
        #[cfg(debug_assertions)]
        {
            t.x_last_fresh = hash.x_last_fresh;
        }
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:428-447
/// ```c
///     if (table->occupancy == table->capacity) {
///         table->capacity *= 2;
///         table->chain =
///             realloc(table->chain, table->capacity * sizeof(DiagHashCell));
/// ...
///     cell = table->chain + table->occupancy;
///     cell->diag = diag;
///     cell->level = level;
///     cell->hit_len = len;
///     cell->hit_saved = hit_saved;
///     cell->next = table->backbone[bucket];
///     table->backbone[bucket] = table->occupancy;
///     table->occupancy++;
/// ```
/// The end of a partitioned stage (module comment): the hits in key order after `hits`; the
/// appended cells numbered `occ0 + rank` of their keys; the chains of the changed buckets of
/// every partition written into `hash` with those numbers (each cell at its number, its link
/// and the bucket's head mapped; the partition tables themselves are left as they are); the
/// capacity of the reference's doublings.
fn merge(
    hash: &mut DiagHashTable,
    hits: &mut Vec<UngappedHit>,
    states: &mut [PartState],
    occ0: u32,
    order: &mut Vec<(u64, u32, u32)>,
    apps: &mut Vec<(u64, u32, u32)>,
) {
    let parts = states.len();
    order.clear();
    for (p, st) in states.iter().enumerate() {
        debug_assert_eq!(st.hits.len(), st.out.hit_keys.len());
        order.extend(
            st.out
                .hit_keys
                .iter()
                .enumerate()
                .map(|(j, &key)| (key, p as u32, j as u32)),
        );
    }
    order.sort_unstable();
    hits.reserve(order.len());
    for &(_, p, j) in order.iter() {
        hits.push(states[p as usize].hits[j as usize]);
    }

    apps.clear();
    for (p, st) in states.iter().enumerate() {
        apps.extend(
            st.out
                .appends
                .iter()
                .map(|&(key, local)| (key, p as u32, local)),
        );
    }
    apps.sort_unstable();
    for st in states.iter_mut() {
        st.map.clear();
        st.map.resize(st.out.appends.len(), 0);
    }
    for (rank, &(_, p, local)) in apps.iter().enumerate() {
        // a partition numbers its appends occ0, occ0 + 1, ... in its key order
        states[p as usize].map[(local - occ0) as usize] = occ0 + rank as u32;
    }
    let occ_final = occ0 + apps.len() as u32;
    let mut capacity = hash.capacity;
    while capacity < occ_final {
        capacity = capacity.saturating_mul(2);
    }
    if capacity as usize > hash.chain.len() {
        hash.chain
            .resize(capacity as usize, DiagHashCell::default());
    }
    hash.capacity = capacity;
    hash.occupancy = occ_final;
    let levels = hash.x_lv_on();
    for (p, st) in states.iter().enumerate() {
        let PartState {
            hash: t, out, map, ..
        } = st;
        // a cell's number in the stage: the cells below occ0 keep theirs
        let remap = |index: u32| -> u32 {
            if index < occ0 {
                index
            } else {
                map[(index - occ0) as usize]
            }
        };
        for &bucket in &out.dirty {
            let b = bucket as usize;
            debug_assert_eq!(b % parts, p);
            let mut local = t.backbone[b];
            hash.backbone[b] = remap(local);
            while local != 0 {
                let mut cell = t.chain[local as usize];
                let next = cell.next;
                cell.next = remap(next);
                hash.chain[remap(local) as usize] = cell;
                local = next;
            }
            if levels {
                hash.max_level[b] = t.max_level[b];
            }
            // EXPERIMENT (LOSAT_X_DIAGLV2): the bucket's memo is its partition's.
            if let Some(last) = hash.x_last.get_mut(b) {
                *last = t.x_last[b];
            }
        }
        #[cfg(debug_assertions)]
        {
            hash.x_last_fresh &= t.x_last_fresh;
        }
    }
}

// No NCBI counterpart: the counters of the summary fields, per fold (0: one hit, 1: two hits,
// 2: one hit of a direct or Na table, the child module `na`).
pub(super) const ONE_HIT: usize = 0;
pub(super) const TWO_HITS: usize = 1;
pub(super) const NA_ONE_HIT: usize = 2;
static MERGES: [AtomicU64; 3] = [const { AtomicU64::new(0) }; 3];
static PARTS: [AtomicU64; 3] = [const { AtomicU64::new(0) }; 3];
static MERGE_NS: [AtomicU64; 3] = [const { AtomicU64::new(0) }; 3];
static LEFT: [AtomicU64; 3] = [const { AtomicU64::new(0) }; 3];

/// No NCBI counterpart: counts one partitioned stage of `kind` (once per stage).
pub(super) fn record(kind: usize, parts: usize, stats: &PartStats) {
    MERGES[kind].fetch_add(1, Ordering::Relaxed);
    PARTS[kind].fetch_max(parts as u64, Ordering::Relaxed);
    MERGE_NS[kind].fetch_add(stats.merge_ns, Ordering::Relaxed);
    LEFT[kind].fetch_add(stats.helpers_left, Ordering::Relaxed);
}

/// No NCBI counterpart: the fields the shadow summary lines of `kind` end with when the switch
/// is on: `fold=part`, the partition count of the stages and the number of partitioned stages
/// merged. Without the switch the lines are those of the strips paths alone.
pub(super) fn shadow_fields(kind: usize) -> String {
    if !switch_on() {
        return String::new();
    }
    let load = |c: &AtomicU64| c.load(Ordering::Relaxed);
    format!(
        " fold=part partitions={} merges={}",
        load(&PARTS[kind]),
        load(&MERGES[kind])
    )
}

/// No NCBI counterpart: the fields the LOSAT_TIMING lines of `kind` end with when the switch is
/// on: the partition count of the stages, the number of partitioned stages merged, the time in
/// merges and the number of helpers that left a stalled stage. Without the switch the lines are
/// those of the strips paths alone.
pub(super) fn timing_fields(kind: usize) -> String {
    if !switch_on() {
        return String::new();
    }
    let load = |c: &AtomicU64| c.load(Ordering::Relaxed);
    format!(
        " partitions={} merges={} merge={:.3}s helpers_left={}",
        load(&PARTS[kind]),
        load(&MERGES[kind]),
        load(&MERGE_NS[kind]) as f64 / 1e9,
        load(&LEFT[kind])
    )
}

// (every test runs stages in a pool of its own)
#[cfg(all(
    test,
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
mod tests {
    use super::super::super::{advance_diag_table_offset, x_diag_lv::LvSeeds};
    use super::super::copy_hash;
    use super::*;

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 33
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n.max(1)
        }
    }

    type HashState = (
        Vec<u32>,
        Vec<(i32, i32, bool, i32, u32)>,
        u32,
        u32,
        i32,
        Vec<i32>,
        Vec<super::super::super::x_diag_lv::LvLast>,
    );

    /// Everything of a table a later seed or the shadow can read: the heads, the cells up to
    /// the occupancy, the occupancy, the capacity, the offset, the level bounds and the memo
    /// (LOSAT_X_DIAGLV2).
    fn hash_state(h: &DiagHashTable) -> HashState {
        assert_eq!(
            h.chain.len(),
            h.capacity as usize,
            "capacity and chain length"
        );
        (
            h.backbone.clone(),
            h.chain[..h.occupancy as usize]
                .iter()
                .map(|c| (c.diag, c.level, c.hit_saved, c.hit_len, c.next))
                .collect(),
            h.occupancy,
            h.capacity,
            h.offset,
            h.max_level.clone(),
            h.x_last.clone(),
        )
    }

    // A seed outcome of the one-hit fold, as `fold` takes it: the retrieve, the `last_hit`
    // test, then an insert of a seed that is extended below the cutoff (`hit_len` =
    // `s_end_pos - s_off_pos`, not saved) or above it (a hit; the ungapped end, saved).
    #[derive(Clone, Copy, Debug)]
    struct Op {
        diag: i32,
        s_off_pos: i32,
        s_end_pos: i32,
        pass: bool,
        ungapped_end: i32,
    }

    /// The seed of op `index`: its diagonal as `s - q`.
    fn seed_of(index: usize, op: &Op) -> Seed {
        const Q0: i64 = 1 << 20;
        Seed {
            context: index as u32,
            q: Q0 as u32,
            s: (Q0 + i64::from(op.diag)) as u32,
            s_range: 0,
        }
    }

    fn last_hit<const MEMO: bool>(hash: &DiagHashTable, op: &Op) -> i32 {
        let found = if hash.x_lv_on() {
            hash.retrieve_below::<MEMO>(op.diag, &LvSeeds::new(false, 11, 0), op.s_off_pos as usize)
        } else {
            hash.retrieve(op.diag)
        };
        found.map_or(0, |(level, _, _)| level)
    }

    /// The reference: op after op on one table; the hits are op indices. With `MEMO` it reads
    /// and keeps the table's memo (LOSAT_X_DIAGLV2) as the folds of a strips stage do, without it
    /// neither, as the reference fold.
    fn reference<const MEMO: bool>(
        hash: &mut DiagHashTable,
        ops: &[Op],
        first: usize,
        window: i32,
    ) -> Vec<usize> {
        let insert = if MEMO {
            DiagHashTable::insert_memo
        } else {
            DiagHashTable::insert_plain
        };
        let mut hits = Vec::new();
        for (k, op) in ops.iter().enumerate() {
            if op.s_off_pos < last_hit::<MEMO>(hash, op) {
                continue;
            }
            if op.pass {
                hits.push(first + k);
                insert(
                    hash,
                    op.diag,
                    op.ungapped_end,
                    0,
                    true,
                    op.s_off_pos,
                    window,
                );
            } else {
                insert(
                    hash,
                    op.diag,
                    op.s_end_pos,
                    op.s_end_pos - op.s_off_pos,
                    false,
                    op.s_off_pos,
                    window,
                );
            }
        }
        hits
    }

    fn hit_of(index: usize) -> UngappedHit {
        UngappedHit {
            context_idx: 0,
            query_idx: 0,
            query_frame: 1,
            query_context_offset: 0,
            seed_q_off: index,
            seed_s_off: 0,
            qs: 0,
            qe: 0,
            ss: 0,
            se: 0,
            score: 0,
        }
    }

    /// The same ops as a partitioned stage: strips cut at `cuts`, `parts` partitions on
    /// `threads` threads (in a pool of that many). `nested` runs a parallel iterator inside
    /// every scan and fold call. With `owner_waits` a fold on the owner's thread first sleeps
    /// 100 us while no helper has folded yet in this stage, so that the helpers get to fold on
    /// any machine (even on one CPU, where the owner could otherwise fold the whole stage
    /// before a helper job starts). Returns the hits (op indices), the stage's counters and
    /// the number of threads that folded.
    #[allow(clippy::too_many_arguments)]
    fn partitioned(
        pool: &rayon::ThreadPool,
        hash: &mut DiagHashTable,
        scratch: &mut PartScratch,
        ops: &[Op],
        first: usize,
        cuts: &[usize],
        parts: usize,
        ahead: usize,
        window: i32,
        nested: bool,
        owner_waits: bool,
    ) -> (Vec<usize>, PartStats, usize) {
        let threads = pool.current_num_threads();
        let folders = Mutex::new(std::collections::BTreeSet::new());
        let helper_folded = AtomicBool::new(false);
        pool.install(|| {
            let owner = rayon::current_thread_index();
            let scan = |i: usize, out: &mut Vec<Seed>| {
                if nested {
                    use rayon::prelude::*;
                    let sum: u64 = (0..64u64).into_par_iter().map(|v| v * v).sum();
                    assert_eq!(sum, 85344);
                }
                for k in cuts[i]..cuts[i + 1] {
                    out.push(seed_of(first + k, &ops[k]));
                }
            };
            let fold = |st: &mut PartState, seeds: &[Seed], pos: &[u32]| {
                if nested {
                    use rayon::prelude::*;
                    let sum: u64 = (0..64u64).into_par_iter().map(|v| v + 1).sum();
                    assert_eq!(sum, 2080);
                }
                let here = rayon::current_thread_index();
                lock(&folders).insert(here);
                if owner_waits {
                    if here != owner {
                        helper_folded.store(true, Ordering::Relaxed);
                    } else if !helper_folded.load(Ordering::Relaxed) {
                        std::thread::sleep(std::time::Duration::from_micros(100));
                    }
                }
                for (k, seed) in seeds.iter().enumerate() {
                    let index = seed.context as usize;
                    let op = &ops[index - first];
                    let key = st.out.key(pos, k);
                    st.counts.seeds += 1;
                    if op.s_off_pos < last_hit::<true>(&st.hash, op) {
                        st.counts.skipped += 1;
                        continue;
                    }
                    if op.pass {
                        st.hits.push(hit_of(index));
                        st.out.hit_keys.push(key);
                        insert::<true, false, true>(
                            &mut st.hash,
                            &mut st.out,
                            key,
                            op.diag,
                            op.ungapped_end,
                            0,
                            true,
                            op.s_off_pos,
                            window,
                        );
                    } else {
                        insert::<true, false, true>(
                            &mut st.hash,
                            &mut st.out,
                            key,
                            op.diag,
                            op.s_end_pos,
                            op.s_end_pos - op.s_off_pos,
                            false,
                            op.s_off_pos,
                            window,
                        );
                    }
                }
            };
            let mut hits = Vec::new();
            let stats = seed_stage_part(
                cuts.len() - 1,
                threads,
                ahead,
                hash,
                &mut hits,
                FoldPart { scratch, parts },
                &scan,
                &fold,
                false,
            );
            let folders = lock(&folders).len();
            (hits.iter().map(|h| h.seed_q_off).collect(), stats, folders)
        })
    }

    /// Random ops: diagonals of a few buckets' worth (d + 512 k share a bucket) and many
    /// others, positions from `s` on, a few per position.
    fn random_ops(rng: &mut Lcg, n: usize, s: &mut i32, diagonals: usize) -> Vec<Op> {
        let base = rng.below(1000) as i32 - 500;
        let diags: Vec<i32> = (0..diagonals)
            .map(|k| {
                if k < 16 {
                    base + 512 * k as i32
                } else {
                    rng.below(40_000) as i32 - 20_000
                }
            })
            .collect();
        (0..n)
            .map(|_| {
                *s += rng.below(4) as i32;
                let s_end = *s + 11 + rng.below(4) as i32;
                Op {
                    diag: diags[rng.below(diags.len() as u64) as usize],
                    s_off_pos: *s,
                    s_end_pos: s_end,
                    pass: rng.below(3) == 0,
                    ungapped_end: s_end + rng.below(200) as i32,
                }
            })
            .collect()
    }

    fn random_cuts(rng: &mut Lcg, n: usize) -> Vec<usize> {
        let mut cuts = vec![0usize];
        while *cuts.last().unwrap() < n {
            // some empty strips
            let step = if rng.below(10) == 0 {
                0
            } else {
                1 + rng.below(80) as usize
            };
            cuts.push((cuts.last().unwrap() + step).min(n));
        }
        cuts
    }

    fn advance(hash: &mut DiagHashTable, subject_len: usize) {
        let mut diag_table_offset = 0;
        let (mut levels, mut lens) = (Vec::new(), Vec::new());
        let window = hash.window;
        advance_diag_table_offset(
            &mut diag_table_offset,
            window,
            subject_len,
            false,
            &mut levels,
            &mut lens,
            hash,
        );
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:396-449
    // ```c
    //     while (index) {
    //         if (table->chain[index].diag == diag) {
    // ...
    //             if (s_off - table->chain[index].level > window_size) {
    // ...
    //     if (table->occupancy == table->capacity) {
    //         table->capacity *= 2;
    // ```
    // `seed_stage_part` against the reference on random op streams, for 2, 3, 7 and 16
    // partitions on 1..8 threads, with and without level bounds, all stages of a case on one
    // scratch: a stage on a new table that grows past 256 cells (many diagonals within the
    // insert window), a second stage after the offset moved on (the next subject or chunk), a
    // third after the clear of the offset (occupancy 1), a chunk the reference folds, a stage
    // after it, a stage of no strips, and a stage with another partition count. Every stage
    // gives the reference's hits and table cell for cell (heads, cells, occupancy, capacity,
    // bounds). In every other case on 2+ threads the owner sleeps in its folds until a helper
    // has folded (`partitioned`), so that these cases check the equality on folds by several
    // threads on any machine, which the test asserts for each of them.
    #[test]
    fn seed_stage_part_equals_the_reference_fold() {
        let mut rng = Lcg(0x9A27_F01D);
        let mut case = 0usize;
        let mut waited_cases = 0usize;
        for parts in [2usize, 3, 7, 16] {
            for threads in 1..=8usize {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                for levels in [false, true] {
                    case += 1;
                    let window = [12, 40, 5000][case % 3];
                    // the bounds (LOSAT_X_DIAGLV) and the memo (LOSAT_X_DIAGLV2) as `levels` says,
                    // whatever the environment says
                    let mut reference_hash = DiagHashTable::new(16);
                    reference_hash.max_level = Vec::new();
                    reference_hash.x_last = Vec::new();
                    if levels {
                        reference_hash.x_lv_start(false);
                        reference_hash.x_lv2_start();
                    }
                    let mut hash = copy_hash(&reference_hash);
                    let mut scratch = PartScratch::default();
                    let mut s = 50i32;
                    let mut first = 0usize;
                    let mut most_folders = 0usize;
                    // one of the two cases of each (parts, threads), alternating with threads
                    let owner_waits =
                        threads >= 2 && (threads + usize::from(levels)).is_multiple_of(2);
                    for stage in 0..7 {
                        let at = format!("case {case} parts {parts} threads {threads} levels {levels} stage {stage}");
                        match stage {
                            1 | 4 => {
                                advance(&mut reference_hash, 3000);
                                advance(&mut hash, 3000);
                            }
                            2 => {
                                // the clear: offset >= INT4_MAX / 4
                                reference_hash.offset = i32::MAX / 4;
                                hash.offset = i32::MAX / 4;
                                advance(&mut reference_hash, 3000);
                                advance(&mut hash, 3000);
                                assert_eq!(hash.occupancy, 1);
                                s = 50;
                            }
                            _ => {}
                        }
                        let n = if stage == 5 {
                            0
                        } else {
                            500 + rng.below(2500) as usize
                        };
                        let diagonals = if stage == 0 {
                            1200
                        } else {
                            60 + rng.below(600) as usize
                        };
                        let ops = random_ops(&mut rng, n, &mut s, diagonals);
                        // as the search: a strips stage makes the memo (LOSAT_X_DIAGLV2) again at
                        // its start (`seed_stage`) and keeps it; the reference fold (stage 3)
                        // neither reads nor keeps it
                        let expected = if stage == 3 {
                            reference::<false>(&mut reference_hash, &ops, first, window)
                        } else {
                            if !reference_hash.x_last.is_empty() {
                                reference_hash.x_lv2_reset();
                            }
                            reference::<true>(&mut reference_hash, &ops, first, window)
                        };
                        if stage == 3 {
                            // a chunk folded by the reference (the partition tables do not see it)
                            let _ = reference::<false>(&mut hash, &ops, first, window);
                        } else {
                            if !hash.x_last.is_empty() {
                                hash.x_lv2_reset();
                            }
                            let stage_parts = if stage == 6 { parts + 1 } else { parts };
                            let cuts = if n == 0 {
                                vec![0]
                            } else {
                                random_cuts(&mut rng, n)
                            };
                            let ahead = 1 + rng.below(2 * threads as u64) as usize;
                            let (got, stats, folders) = partitioned(
                                &pool,
                                &mut hash,
                                &mut scratch,
                                &ops,
                                first,
                                &cuts,
                                stage_parts,
                                ahead,
                                window,
                                false,
                                owner_waits,
                            );
                            most_folders = most_folders.max(folders);
                            assert_eq!(got, expected, "{at}: hits");
                            assert_eq!(stats.counts.seeds, n as u64, "{at}");
                        }
                        assert_eq!(
                            hash_state(&hash),
                            hash_state(&reference_hash),
                            "{at}: table"
                        );
                        if stage == 0 {
                            assert!(
                                hash.occupancy > 256,
                                "{at}: the table did not grow ({})",
                                hash.occupancy
                            );
                        }
                        first += n;
                    }
                    if owner_waits {
                        // the helpers folded too
                        waited_cases += 1;
                        assert!(
                            most_folders >= 2,
                            "case {case} parts {parts} threads {threads} levels {levels}: the owner waited in its folds, but no stage was folded on more than one thread"
                        );
                    }
                }
            }
        }
        assert_eq!(waited_cases, 28);
    }

    // A parallel iterator inside every scan and fold call of every participant (the pattern
    // that hung LOSAT_X_BNSTREAM, #218): a thread waiting inside it (the participant, or a thread
    // that stole a piece of the call) may steal a helper job of the stage, which must not wait
    // for the stage while the stage waits for the call below it (the helper returns at once on a
    // participant's thread, and leaves a stage that has stalled for `HELPER_STALL`). The stages
    // end and give the reference's hits and table.
    #[test]
    fn a_nested_parallel_call_in_every_participant_does_not_hang() {
        let mut rng = Lcg(0x4E57_ED00);
        for threads in 2..=8usize {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            let mut reference_hash = DiagHashTable::new(16);
            let mut hash = copy_hash(&reference_hash);
            let mut scratch = PartScratch::default();
            let mut s = 50i32;
            let mut first = 0usize;
            for stage in 0..3 {
                let n = 2000 + rng.below(2000) as usize;
                let ops = random_ops(&mut rng, n, &mut s, 400);
                let expected = reference::<true>(&mut reference_hash, &ops, first, 40);
                let cuts = random_cuts(&mut rng, n);
                let (got, _, _) = partitioned(
                    &pool,
                    &mut hash,
                    &mut scratch,
                    &ops,
                    first,
                    &cuts,
                    threads,
                    2 * threads,
                    40,
                    true,
                    false,
                );
                assert_eq!(got, expected, "threads {threads} stage {stage}");
                assert_eq!(
                    hash_state(&hash),
                    hash_state(&reference_hash),
                    "threads {threads} stage {stage}"
                );
                first += n;
            }
        }
    }

    // A participant that panics (in a fold on the owner, in a fold on a helper, in a scan) makes
    // the stage panic instead of hang and leaves the table untouched, and the next stage on the
    // same scratch gives the reference's result.
    #[test]
    fn a_panicking_participant_ends_the_stage() {
        let mut rng = Lcg(0x9A41_C000);
        for threads in [2usize, 4, 8] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            for (case, place) in ["owner fold", "helper fold", "scan"]
                .into_iter()
                .enumerate()
            {
                let mut reference_hash = DiagHashTable::new(16);
                let mut hash = copy_hash(&reference_hash);
                let mut scratch = PartScratch::default();
                let mut s = 50i32;
                let n = 3000;
                let ops = random_ops(&mut rng, n, &mut s, 300);
                let expected = reference::<true>(&mut reference_hash, &ops, 0, 40);
                let cuts = random_cuts(&mut rng, n);
                let strips = cuts.len() - 1;
                let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    pool.install(|| {
                        let owner = rayon::current_thread_index();
                        let scan = |i: usize, out: &mut Vec<Seed>| {
                            if place == "scan" && i == strips / 2 {
                                panic!("a test scan panics");
                            }
                            for k in cuts[i]..cuts[i + 1] {
                                out.push(seed_of(k, &ops[k]));
                            }
                        };
                        let fold = |st: &mut PartState, seeds: &[Seed], _pos: &[u32]| {
                            let on_owner = rayon::current_thread_index() == owner;
                            match place {
                                "owner fold" if on_owner && st.out.strip >= 2 => {
                                    panic!("a test fold panics on the owner")
                                }
                                "helper fold" if !on_owner => {
                                    panic!("a test fold panics on a helper")
                                }
                                // the owner leaves the folds to the helpers, or they to it
                                "helper fold" | "owner fold" => {
                                    std::thread::sleep(std::time::Duration::from_millis(1))
                                }
                                _ => {}
                            }
                            st.counts.seeds += seeds.len() as u64;
                        };
                        let mut hits = Vec::new();
                        seed_stage_part(
                            strips,
                            threads,
                            2,
                            &mut hash,
                            &mut hits,
                            FoldPart {
                                scratch: &mut scratch,
                                parts: threads,
                            },
                            &scan,
                            &fold,
                            false,
                        )
                        .counts
                        .seeds
                    })
                }));
                assert!(
                    result.is_err(),
                    "threads {threads} case {case} ({place}): no panic"
                );
                // the table is the one before the stage: a merge never ran
                let fresh = DiagHashTable::new(16);
                assert_eq!(hash_state(&hash), hash_state(&fresh), "{place}");
                let (got, _, _) = partitioned(
                    &pool,
                    &mut hash,
                    &mut scratch,
                    &ops,
                    0,
                    &cuts,
                    threads,
                    2,
                    40,
                    false,
                    false,
                );
                assert_eq!(got, expected, "threads {threads} {place}");
                assert_eq!(
                    hash_state(&hash),
                    hash_state(&reference_hash),
                    "threads {threads} {place}"
                );
            }
        }
    }

    // Every stage starts from a copy of the stage's table: a change of the table between two
    // partitioned stages on one scratch that keeps the occupancy (here the level of every cell
    // raised, which no fold does) is seen by the next stage, which gives the reference's result
    // on the changed table and not the one on the table before the change.
    #[test]
    fn a_table_changed_between_stages_is_seen() {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap();
        let mut rng = Lcg(0xBAD0_7AB1);
        let mut reference_hash = DiagHashTable::new(16);
        let mut hash = copy_hash(&reference_hash);
        let mut scratch = PartScratch::default();
        let mut s = 50i32;
        let ops = random_ops(&mut rng, 2000, &mut s, 300);
        let expected = reference::<true>(&mut reference_hash, &ops, 0, 40);
        let cuts = random_cuts(&mut rng, ops.len());
        let (got, _, _) = partitioned(
            &pool,
            &mut hash,
            &mut scratch,
            &ops,
            0,
            &cuts,
            2,
            2,
            40,
            false,
            false,
        );
        assert_eq!(got, expected, "first stage: hits");
        assert_eq!(
            hash_state(&hash),
            hash_state(&reference_hash),
            "first stage: table"
        );
        // the same diagonals again, further on
        let more: Vec<Op> = ops
            .iter()
            .map(|op| Op {
                s_off_pos: op.s_off_pos + 1000,
                s_end_pos: op.s_end_pos + 1000,
                ungapped_end: op.ungapped_end + 1000,
                ..*op
            })
            .collect();
        let unchanged = reference::<true>(&mut copy_hash(&reference_hash), &more, 2000, 40);
        let occupancy = hash.occupancy;
        for h in [&mut hash, &mut reference_hash] {
            for cell in &mut h.chain[1..occupancy as usize] {
                cell.level = i32::MAX / 2;
            }
            // the bounds (LOSAT_X_DIAGLV) and the memo (LOSAT_X_DIAGLV2) describe the cells the
            // inserts wrote; made again from the changed cells
            if h.x_lv_on() {
                let shadow = h.x_lv_shadow;
                let memo = !h.x_last.is_empty();
                h.x_lv_start(shadow);
                if memo {
                    h.x_lv2_start();
                }
            }
        }
        let expected = reference::<true>(&mut reference_hash, &more, 2000, 40);
        assert_ne!(expected, unchanged, "the change does not matter");
        let cuts = random_cuts(&mut rng, more.len());
        let (got, _, _) = partitioned(
            &pool,
            &mut hash,
            &mut scratch,
            &more,
            2000,
            &cuts,
            2,
            2,
            40,
            false,
            false,
        );
        assert_eq!(got, expected, "second stage: hits");
        assert_eq!(
            hash_state(&hash),
            hash_state(&reference_hash),
            "second stage: table"
        );
    }
}
