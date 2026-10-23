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
//! the whole stage has made no progress (no strip scanned, no partition folded; with
//! LOSAT_X_PAIRPARHOME no strip handed out, no strip folded by any partition) for
//! `HELPER_STALL` (the owner never leaves). A thread that waits inside a piece of a
//! participant's nested call that it stole may run a helper job above the work it interrupted;
//! that helper waits for the stage while the stage waits for that work (the pattern that hung
//! LOSAT_X_BNSTREAM, #218; test `a_nested_parallel_call_in_every_participant_does_not_hang`)
//! until the stall exit releases it. At the end the hits are put in the order of their seeds'
//! keys (strip, position in the strip's scan order) and the changed buckets of every partition
//! are written into the stage's table with the reference's cell numbers (`merge`).
//! With LOSAT_X_PAIRPARHOME (on by default) the participants share the folds otherwise: each
//! partition has a home participant that folds it run by run of ready strips, and the others
//! fold it only in the cases `seed_stage_home` lists; the argument below needs only that each
//! partition folds its strips in order with one folder at a time, which both ways keep.
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
use std::sync::{Arc, Mutex, OnceLock, PoisonError, RwLock};

/// No NCBI counterpart: reads LOSAT_X_PAIRPARFOLD once (a strict-set switch: on unless it is `0`,
/// or `LOSAT_X_DEFAULT=0` and it is not set). It does not change any value NCBI computes.
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

/// No NCBI counterpart: reads LOSAT_X_PAIRPARHOME once (a strict-set switch: on unless it is `0`,
/// or `LOSAT_X_DEFAULT=0` and it is not set): the participants of a partitioned stage fold
/// their home partitions (`seed_stage_home`) instead of any partition they can lock. It does not
/// change any value NCBI computes.
fn home_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_PAIRPARHOME"))
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
/// A fold's insert (`DiagHashTable::insert`). With partitions (`PT`) it also records the
/// diagonal's bucket and, when the insert appended a cell, the seed's key with the cell's
/// number in the partition's table.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
pub(super) fn insert<const PT: bool>(
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
        hash.insert(diag, level, hit_len, hit_saved, s_off, window_size);
        let bucket = bucket_of_diag(diag);
        if !out.touched[bucket] {
            out.touched[bucket] = true;
            out.dirty.push(bucket as u16);
        }
        if hash.occupancy != occupancy {
            out.appends.push((key, occupancy));
        }
    } else {
        hash.insert(diag, level, hit_len, hit_saved, s_off, window_size);
    }
}

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
    /// LOSAT_X_PAIRPARHOME: the runs folded (locks of a partition that went through one strip
    /// or more, empty partitions of strips included), and the strips of partitions folded by a
    /// participant whose home the partition is not (a partition without a participant at home,
    /// the slowest partition, any partition once every strip was handed out)
    pub(super) runs: u64,
    pub(super) away: u64,
    /// the partitions of strips of the stage (strips x partitions)
    pub(super) units: u64,
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
/// the stage to make progress (a strip scanned or a partition's strip folded; with
/// LOSAT_X_PAIRPARHOME a strip handed out or a partition's strip folded) before it leaves. A helper waits while a partition that lags
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
    // EXPERIMENT (LOSAT_X_PAIRPARHOME): home partitions, folded run by run of ready strips.
    let home = home_on();
    if home {
        // every partitioned stage takes the switch's path
        crate::utils::xreach::x_reach!("LOSAT_X_PAIRPARHOME", true);
    }
    seed_stage_in(
        home, n, threads, ahead, hash, hits, fold_part, scan, fold, timed,
    )
}

/// No NCBI counterpart (scheduling only): `seed_stage_part` with the sharing of the folds given
/// (`home`: LOSAT_X_PAIRPARHOME, `seed_stage_home`; else any participant folds any partition it
/// can lock).
#[allow(clippy::too_many_arguments)]
fn seed_stage_in<FS, FF>(
    home: bool,
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
    if home {
        return seed_stage_home(
            n, threads, ahead, hash, hits, scratch, parts, occ0, scan, fold, timed,
        );
    }

    let states: Vec<Mutex<PartState>> = scratch.parts.drain(..).map(Mutex::new).collect();
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
                let Ok(mut st) = states[p].try_lock() else {
                    continue;
                };
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
                    fold(&mut st, &strip.seeds[a..b], &strip.pos[a..b]);
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
                drop(st);
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
        .map(|m| m.into_inner().unwrap_or_else(PoisonError::into_inner))
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

/// No NCBI counterpart (scheduling only, LOSAT_X_PAIRPARHOME): a value on cache lines of its own
/// (128 bytes: a line and the adjacent line the prefetcher pairs with it).
#[derive(Default)]
#[repr(align(128))]
struct OwnLines<T>(T);

/// No NCBI counterpart (scheduling only, LOSAT_X_PAIRPARHOME): marks the home partitions of a
/// participant (`p % threads == me`) as having their participant in the stage while it lives,
/// so that a helper that leaves a stalled stage (or unwinds) hands them to every participant.
/// (The flags change twice per participant and stage, so they share lines; the partitions'
/// next strips, which change with every strip folded, have lines of their own.)
struct AtHome<'a> {
    homes: &'a [AtomicBool],
    me: usize,
    threads: usize,
}

impl<'a> AtHome<'a> {
    fn enter(homes: &'a [AtomicBool], me: usize, threads: usize) -> Self {
        for home in homes.iter().skip(me).step_by(threads) {
            home.store(true, Ordering::Relaxed);
        }
        AtHome { homes, me, threads }
    }
}

impl Drop for AtHome<'_> {
    fn drop(&mut self) {
        for home in self.homes.iter().skip(self.me).step_by(self.threads) {
            home.store(false, Ordering::Relaxed);
        }
    }
}

fn read<T>(m: &RwLock<T>) -> std::sync::RwLockReadGuard<'_, T> {
    m.read().unwrap_or_else(PoisonError::into_inner)
}

fn write<T>(m: &RwLock<T>) -> std::sync::RwLockWriteGuard<'_, T> {
    m.write().unwrap_or_else(PoisonError::into_inner)
}

/// No NCBI counterpart (scheduling only, LOSAT_X_PAIRPARHOME): why a participant folds a
/// partition, which says how far its run goes (`seed_stage_home`).
#[derive(Clone, Copy, PartialEq, Eq)]
enum Run {
    /// its home: every ready strip
    Home,
    /// a partition without a participant at home in the stage: until its participant comes
    Orphan,
    /// the slowest partition, by a participant with nothing else to do: one strip
    Slowest,
    /// any partition once every strip is handed out: every ready strip
    Tail,
}

/// No NCBI counterpart (scheduling only, LOSAT_X_PAIRPARHOME): the stage of `seed_stage_part`
/// after `rebuild`, with the folds shared by home partitions. Partition p is the home of
/// participant `p % threads`, which folds p's strips as they become ready, every ready strip
/// from p's next one on in one lock (a run); it scans a strip within the window between runs.
/// Another participant folds p only (`Run`) while p has no participant at home in the stage
/// (its helper job has not started, started on a participant's thread and returned, or left a
/// stalled stage), until that participant comes; one strip when p is the slowest partition
/// and the participant has nothing else to do (the window is closed); and every ready strip
/// when the participant has nothing else to do once every strip is handed out. Partition p keeps its strips' order and one folder at a time
/// (its lock) as in the stage without the switch, so the output is the same (module comment);
/// only who folds which run changes. A strip's buffer is given back once every partition has
/// folded it (the slowest partition's next strip is past it; `freed`), without a count per
/// strip; the strips are read under a read lock instead of a shared pointer per fold unit; the
/// helpers' stall test sums the cursor and the partitions' next strips instead of a shared
/// count per fold unit; and each partition's next strip and state are on lines of their own.
/// The owner alone can still scan every strip and fold every partition.
#[allow(clippy::too_many_arguments)]
fn seed_stage_home<FS, FF>(
    n: usize,
    threads: usize,
    ahead: usize,
    hash: &mut DiagHashTable,
    hits: &mut Vec<UngappedHit>,
    scratch: &mut PartScratch,
    parts: usize,
    occ0: u32,
    scan: &FS,
    fold: &FF,
    timed: bool,
) -> PartStats
where
    FS: Fn(usize, &mut Vec<Seed>) + Sync,
    FF: Fn(&mut PartState, &[Seed], &[u32]) + Sync,
{
    let threads = threads.max(1);
    let states: Vec<Mutex<OwnLines<PartState>>> = scratch
        .parts
        .drain(..)
        .map(|st| Mutex::new(OwnLines(st)))
        .collect();
    let spare = Mutex::new(std::mem::take(&mut scratch.spare));
    let raw_spare = Mutex::new(std::mem::take(&mut scratch.raw));
    // per partition: the next strip it folds (set under its lock by the participant that folded
    // the strip before) and whether a participant whose home it is takes part in the stage
    let lanes: Vec<OwnLines<AtomicUsize>> = (0..parts).map(|_| OwnLines::default()).collect();
    let homes: Vec<AtomicBool> = (0..parts).map(|_| AtomicBool::new(false)).collect();
    let slots: Vec<RwLock<Option<StripBuf>>> = (0..n).map(|_| RwLock::new(None)).collect();
    let ready: Vec<AtomicBool> = (0..n).map(|_| AtomicBool::new(false)).collect();
    let cursor = AtomicUsize::new(0);
    // the strips below it have been given back (every partition folded them)
    let freed = AtomicUsize::new(0);
    let failed = AtomicBool::new(false);
    let left = AtomicUsize::new(0);
    let runs = AtomicU64::new(0);
    let away = AtomicU64::new(0);

    // The slowest partition's next strip (every partition has folded the strips below it); the
    // strips below it go back to the spare buffers, by the participant that moves `freed`.
    let slowest = || -> usize {
        let min_next = lanes
            .iter()
            .map(|lane| lane.0.load(Ordering::Acquire))
            .min()
            .unwrap_or(n);
        // (a load first: idle participants call this in every turn, and a write would take the
        // line from all of them; `freed` moves only under the spare lock, and the buffers move
        // into the spare list under it, so that a producer that waits for the lock finds them;
        // lock order: the spare list, then a slot)
        if freed.load(Ordering::Relaxed) < min_next {
            let mut spare = lock(&spare);
            let f = freed.load(Ordering::Relaxed);
            if f < min_next {
                freed.store(min_next, Ordering::Relaxed);
                for slot in &slots[f..min_next] {
                    if let Some(buf) = write(slot).take() {
                        spare.push(buf);
                    }
                }
            }
        }
        min_next
    };

    // A strip's seeds, cut by partition with a stable counting sort.
    let produce = |i: usize, raw: &mut Vec<Seed>, fill: &mut Vec<u32>| {
        raw.clear();
        scan(i, raw);
        // (no spare buffer: first the strips every partition has passed, which wait for the
        // next participant that reads the slowest partition)
        let mut buf = lock(&spare).pop();
        if buf.is_none() {
            slowest();
            buf = lock(&spare).pop();
        }
        let mut buf = buf.unwrap_or_default();
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
        *write(&slots[i]) = Some(buf);
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
        let _at_home = AtHome::enter(&homes, me, threads);
        let mut raw = lock(&raw_spare).pop().unwrap_or_default();
        let mut fill: Vec<u32> = Vec::with_capacity(parts);
        let (mut scan_ns, mut wait_ns, mut pause_ns) = (0u64, 0u64, 0u64);
        let (mut my_runs, mut my_away) = (0u64, 0u64);
        let mut wait_start: Option<std::time::Instant> = None;
        let mut idle = 0u32;
        let mut stall: Option<(std::time::Instant, usize)> = None;
        // the slowest partition's next strip as this participant last read it (it only grows)
        let mut seen = 0usize;

        // Folds partition p's ready strips from its next one on, in order, if it can lock p, and
        // returns how many: every ready strip (`Run::Home`, `Run::Tail`), while p has no
        // participant at home (`Run::Orphan`), or one strip (`Run::Slowest`: the window opens and
        // p stays with its home).
        let mut run = |p: usize, kind: Run, pause_ns: &mut u64| -> u64 {
            let lane = &lanes[p].0;
            let i = lane.load(Ordering::Acquire);
            if i >= n || !ready[i].load(Ordering::Acquire) {
                return 0;
            }
            let Ok(mut guard) = states[p].try_lock() else {
                return 0;
            };
            let st = &mut guard.0;
            let mut i = st.next;
            let mut folded = 0u64;
            while i < n && ready[i].load(Ordering::Acquire) {
                let slot = read(&slots[i]);
                let strip = slot.as_ref().expect("a ready strip");
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
                        *pause_ns += st.scan_pause_ns - pause0;
                    }
                }
                drop(slot);
                i += 1;
                st.next = i;
                lane.store(i, Ordering::Release);
                folded += 1;
                match kind {
                    Run::Home | Run::Tail => {}
                    Run::Orphan if !homes[p].load(Ordering::Relaxed) => {}
                    Run::Orphan | Run::Slowest => break,
                }
            }
            if folded > 0 {
                my_runs += 1;
                if kind != Run::Home {
                    my_away += folded;
                }
            }
            folded
        };

        loop {
            if failed.load(Ordering::Acquire) {
                if owner {
                    panic!("LOSAT_X_PAIRPARFOLD: a participant of the partitioned fold failed");
                }
                failing.done = true;
                break;
            }
            let mut did = false;
            // its home partitions, and those without a participant at home in the stage
            for k in 0..parts {
                let p = (me + k) % parts;
                let kind = if p % threads == me {
                    Run::Home
                } else if !homes[p].load(Ordering::Relaxed) {
                    Run::Orphan
                } else {
                    continue;
                };
                did |= run(p, kind, &mut pause_ns) > 0;
            }
            let c = cursor.load(Ordering::Relaxed);
            if c >= n || c >= seen + ahead || !did {
                seen = slowest();
                if seen >= n {
                    failing.done = true;
                    break;
                }
            }
            if c < n && c < seen + ahead {
                let i = cursor.fetch_add(1, Ordering::Relaxed);
                if i < n {
                    let t0 = (timed && owner).then(std::time::Instant::now);
                    produce(i, &mut raw, &mut fill);
                    if let Some(t0) = t0 {
                        scan_ns += t0.elapsed().as_nanos() as u64;
                    }
                    did = true;
                }
            }
            if !did {
                // nothing of its own: a strip of the slowest partition, or once every strip is
                // handed out every ready strip of any partition
                let tail = cursor.load(Ordering::Relaxed) >= n;
                for k in 0..parts {
                    let p = (me + k) % parts;
                    if p % threads == me {
                        continue;
                    }
                    let kind = if tail { Run::Tail } else { Run::Slowest };
                    if (tail || lanes[p].0.load(Ordering::Acquire) == seen)
                        && run(p, kind, &mut pause_ns) > 0
                    {
                        did = true;
                        break;
                    }
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
                    // the owner can finish the stage alone; the stage's progress: strips handed
                    // out plus strips folded by each partition
                    let now = cursor.load(Ordering::Relaxed)
                        + lanes
                            .iter()
                            .map(|lane| lane.0.load(Ordering::Relaxed))
                            .sum::<usize>();
                    match stall {
                        Some((since, seen_progress)) if seen_progress == now => {
                            if since.elapsed() >= HELPER_STALL {
                                left.fetch_add(1, Ordering::Relaxed);
                                failing.done = true;
                                break;
                            }
                        }
                        _ => stall = Some((std::time::Instant::now(), now)),
                    }
                }
                pause(idle, owner);
            }
        }
        if let Some(t0) = wait_start {
            wait_ns += t0.elapsed().as_nanos() as u64;
        }
        runs.fetch_add(my_runs, Ordering::Relaxed);
        away.fetch_add(my_away, Ordering::Relaxed);
        lock(&raw_spare).push(raw);
        (scan_ns, wait_ns, pause_ns)
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
        runs: runs.into_inner(),
        away: away.into_inner(),
        units: (n * parts) as u64,
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
    // (the allocations only: the next stage's `rebuild` sets every partition's contents; the
    // participant that saw every partition past the last strip gave every strip back)
    debug_assert!(slots.iter().all(|slot| read(slot).is_none()));
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
static RUNS: [AtomicU64; 3] = [const { AtomicU64::new(0) }; 3];
static AWAY: [AtomicU64; 3] = [const { AtomicU64::new(0) }; 3];
static UNITS: [AtomicU64; 3] = [const { AtomicU64::new(0) }; 3];

/// No NCBI counterpart: counts one partitioned stage of `kind` (once per stage).
pub(super) fn record(kind: usize, parts: usize, stats: &PartStats) {
    MERGES[kind].fetch_add(1, Ordering::Relaxed);
    PARTS[kind].fetch_max(parts as u64, Ordering::Relaxed);
    MERGE_NS[kind].fetch_add(stats.merge_ns, Ordering::Relaxed);
    LEFT[kind].fetch_add(stats.helpers_left, Ordering::Relaxed);
    RUNS[kind].fetch_add(stats.runs, Ordering::Relaxed);
    AWAY[kind].fetch_add(stats.away, Ordering::Relaxed);
    UNITS[kind].fetch_add(stats.units, Ordering::Relaxed);
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
/// merges and the number of helpers that left a stalled stage; with LOSAT_X_PAIRPARHOME also the
/// runs folded and the partitions of strips folded away from home (`PartStats::runs`, `away`), and
/// the partitions of strips of the stages (`units`: strips x partitions). Without the
/// switch the lines are those of the strips paths alone.
pub(super) fn timing_fields(kind: usize) -> String {
    if !switch_on() {
        return String::new();
    }
    let load = |c: &AtomicU64| c.load(Ordering::Relaxed);
    let home = if home_on() {
        format!(
            " runs={} away={} units={}",
            load(&RUNS[kind]),
            load(&AWAY[kind]),
            load(&UNITS[kind])
        )
    } else {
        String::new()
    };
    format!(
        " partitions={} merges={} merge={:.3}s helpers_left={}{home}",
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
    );

    /// Everything of a table a later seed or the shadow can read: the heads, the cells up to
    /// the occupancy, the occupancy, the capacity, the offset and the level bounds.
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

    fn last_hit(hash: &DiagHashTable, op: &Op) -> i32 {
        let found = if hash.x_lv_on() {
            hash.retrieve_below(op.diag, &LvSeeds::new(false, 11, 0), op.s_off_pos as usize)
        } else {
            hash.retrieve(op.diag)
        };
        found.map_or(0, |(level, _, _)| level)
    }

    /// The reference: op after op on one table; the hits are op indices.
    fn reference(hash: &mut DiagHashTable, ops: &[Op], first: usize, window: i32) -> Vec<usize> {
        let mut hits = Vec::new();
        for (k, op) in ops.iter().enumerate() {
            if op.s_off_pos < last_hit(hash, op) {
                continue;
            }
            if op.pass {
                hits.push(first + k);
                hash.insert(op.diag, op.ungapped_end, 0, true, op.s_off_pos, window);
            } else {
                hash.insert(
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
    /// every scan and fold call. With `owner_waits` a fold on the owner's thread first waits
    /// until a helper has folded in this stage (100 us at a time, 200 ms at most in the stage),
    /// so that the helpers get to fold on any machine (even on one CPU, or on a loaded one,
    /// where the owner could otherwise fold the whole stage before a helper job starts). `home`: LOSAT_X_PAIRPARHOME's way of sharing the folds.
    /// Returns the hits (op indices), the stage's counters and the number of threads that
    /// folded.
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
        home: bool,
    ) -> (Vec<usize>, PartStats, usize) {
        let threads = pool.current_num_threads();
        partitioned_on(
            pool,
            threads,
            None,
            hash,
            scratch,
            ops,
            first,
            cuts,
            parts,
            ahead,
            window,
            nested,
            owner_waits,
            home,
        )
    }

    /// How `partitioned_on` slows a participant down.
    #[derive(Clone, Copy, Debug)]
    enum Slow {
        /// every scan call on this pool thread sleeps 1 ms first
        Scans(usize),
        /// the first scan on the owner's thread sleeps 60 ms (well past `HELPER_STALL`)
        OwnerFirstScan,
    }

    /// `partitioned` with the stage's thread count given (`threads` participants whatever the
    /// pool's size) and a participant slowed down as `slow` says.
    #[allow(clippy::too_many_arguments)]
    fn partitioned_on(
        pool: &rayon::ThreadPool,
        threads: usize,
        slow: Option<Slow>,
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
        home: bool,
    ) -> (Vec<usize>, PartStats, usize) {
        let folders = Mutex::new(std::collections::BTreeSet::new());
        let helper_folded = AtomicBool::new(false);
        let owner_waited = AtomicUsize::new(0);
        let owner_scanned = AtomicBool::new(false);
        pool.install(|| {
            let owner = rayon::current_thread_index();
            let slow_scans = match slow {
                Some(Slow::Scans(t)) => Some(t),
                _ => None,
            };
            let scan = |i: usize, out: &mut Vec<Seed>| {
                if slow_scans.is_some() && rayon::current_thread_index() == slow_scans {
                    std::thread::sleep(std::time::Duration::from_millis(1));
                }
                if matches!(slow, Some(Slow::OwnerFirstScan))
                    && rayon::current_thread_index() == owner
                    && !owner_scanned.swap(true, Ordering::Relaxed)
                {
                    std::thread::sleep(std::time::Duration::from_millis(60));
                }
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
                    } else {
                        // (200 ms at most in a stage: on a loaded machine a helper may start late)
                        while !helper_folded.load(Ordering::Relaxed)
                            && owner_waited.load(Ordering::Relaxed) < 2000
                        {
                            std::thread::sleep(std::time::Duration::from_micros(100));
                            owner_waited.fetch_add(1, Ordering::Relaxed);
                        }
                    }
                }
                for (k, seed) in seeds.iter().enumerate() {
                    let index = seed.context as usize;
                    let op = &ops[index - first];
                    let key = st.out.key(pos, k);
                    st.counts.seeds += 1;
                    if op.s_off_pos < last_hit(&st.hash, op) {
                        st.counts.skipped += 1;
                        continue;
                    }
                    if op.pass {
                        st.hits.push(hit_of(index));
                        st.out.hit_keys.push(key);
                        insert::<true>(
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
                        insert::<true>(
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
            let stats = seed_stage_in(
                home,
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
                    let mut reference_hash = DiagHashTable::new(16);
                    if levels {
                        reference_hash.x_lv_start(false);
                    }
                    let mut hash = copy_hash(&reference_hash);
                    let mut scratch = PartScratch::default();
                    let mut s = 50i32;
                    let mut first = 0usize;
                    let mut most_folders = 0usize;
                    let mut most_ahead = 0usize;
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
                        let expected = reference(&mut reference_hash, &ops, first, window);
                        if stage == 3 {
                            // a chunk folded by the reference (the partition tables do not see it)
                            let _ = reference(&mut hash, &ops, first, window);
                        } else {
                            let stage_parts = if stage == 6 { parts + 1 } else { parts };
                            let cuts = if n == 0 {
                                vec![0]
                            } else {
                                random_cuts(&mut rng, n)
                            };
                            let ahead = 1 + rng.below(2 * threads as u64) as usize;
                            most_ahead = most_ahead.max(ahead);
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
                                // both ways of sharing the folds (LOSAT_X_PAIRPARHOME), the
                                // stages of a case alternating
                                (case + stage) % 2 == 1,
                            );
                            most_folders = most_folders.max(folders);
                            assert_eq!(got, expected, "{at}: hits");
                            assert_eq!(stats.counts.seeds, n as u64, "{at}");
                            // the strips' buffers stay within the window: the strips claimed
                            // past the slowest partition (the window and the claims that passed
                            // its test at the same time; a strip given back never to the spare
                            // buffers would make one buffer per strip)
                            assert!(
                                scratch.spare.len() <= most_ahead + threads,
                                "{at}: {} strip buffers, window {most_ahead}",
                                scratch.spare.len()
                            );
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
                let expected = reference(&mut reference_hash, &ops, first, 40);
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
                    // both ways of sharing the folds (LOSAT_X_PAIRPARHOME): stage 1 without
                    stage != 1,
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
            for (case, (place, home)) in ["owner fold", "helper fold", "scan"]
                .into_iter()
                .flat_map(|place| [(place, false), (place, true)])
                .enumerate()
            {
                let mut reference_hash = DiagHashTable::new(16);
                let mut hash = copy_hash(&reference_hash);
                let mut scratch = PartScratch::default();
                let mut s = 50i32;
                let n = 3000;
                let ops = random_ops(&mut rng, n, &mut s, 300);
                let expected = reference(&mut reference_hash, &ops, 0, 40);
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
                        seed_stage_in(
                            home,
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
                    "threads {threads} case {case} ({place}, home {home}): no panic"
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
                    !home,
                );
                assert_eq!(got, expected, "threads {threads} {place} home {home}");
                assert_eq!(
                    hash_state(&hash),
                    hash_state(&reference_hash),
                    "threads {threads} {place}"
                );
            }
        }
    }

    // LOSAT_X_PAIRPARHOME: a home participant whose scans are slow (1 ms each on one pool
    // thread; its partition lags while it scans, and it does not hold its partition's lock
    // then) has its partition folded by the others (`Run::Slowest`, and at the tail
    // `Run::Tail`), and a stage on more participants than the pool has threads (the helper jobs
    // start only after the owner finished the stage alone) has its partitions folded without a
    // participant at home (`Run::Orphan`). Both give the reference's hits and table. (Who folds
    // what depends on the timing: the slow cases get up to five stages to fold a strip away
    // from home.)
    #[test]
    fn home_partitions_are_folded_away_from_home() {
        let mut rng = Lcg(0x4053_A11D);
        for (pool_threads, threads, slow) in [(4usize, 4usize, true), (8, 8, true), (1, 4, false)] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(pool_threads)
                .build()
                .unwrap();
            let owner = pool.install(rayon::current_thread_index);
            let slow = slow.then(|| Slow::Scans((owner.unwrap_or(0) + 1) % pool_threads));
            let mut reference_hash = DiagHashTable::new(16);
            let mut hash = copy_hash(&reference_hash);
            let mut scratch = PartScratch::default();
            let mut s = 50i32;
            let mut first = 0usize;
            let mut away = false;
            for stage in 0..if slow.is_some() { 5 } else { 2 } {
                let at =
                    format!("pool {pool_threads} threads {threads} slow {slow:?} stage {stage}");
                let n = 3000;
                let ops = random_ops(&mut rng, n, &mut s, 400);
                let expected = reference(&mut reference_hash, &ops, first, 40);
                let cuts: Vec<usize> = (0..=60).map(|k| k * n / 60).collect();
                let (got, stats, _) = partitioned_on(
                    &pool,
                    threads,
                    slow,
                    &mut hash,
                    &mut scratch,
                    &ops,
                    first,
                    &cuts,
                    threads,
                    2,
                    40,
                    false,
                    false,
                    true,
                );
                assert_eq!(got, expected, "{at}: hits");
                assert_eq!(
                    hash_state(&hash),
                    hash_state(&reference_hash),
                    "{at}: table"
                );
                if pool_threads == 1 {
                    // the owner alone, every partition but its own without a participant at home
                    assert_eq!(stats.away, 60 * (threads as u64 - 1), "{at}");
                }
                first += n;
                away |= stats.away > 0;
                if away && slow.is_some() {
                    break;
                }
            }
            assert!(
                away,
                "pool {pool_threads} threads {threads} slow {slow:?}: no strip folded away from home"
            );
        }
    }

    // A stage that stalls for longer than `HELPER_STALL` (the owner's first scan sleeps 60 ms,
    // and no partition can fold past that strip) makes the idle helpers leave
    // (`helpers_left`); the owner finishes the stage alone, with LOSAT_X_PAIRPARHOME folding the
    // helpers' home partitions away from home (they have no participant at home any more).
    // Both ways give the reference's hits and table. (On a loaded machine the helper jobs may
    // start only after the owner's sleep, and then nothing stalls, so each way gets up to five
    // stages to see a helper leave.)
    #[test]
    fn a_stalled_stage_lets_the_helpers_leave() {
        let mut rng = Lcg(0x57A1_1ED0);
        for threads in [2usize, 4] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            for home in [false, true] {
                let mut reference_hash = DiagHashTable::new(16);
                let mut hash = copy_hash(&reference_hash);
                let mut scratch = PartScratch::default();
                let mut s = 50i32;
                let mut first = 0usize;
                let mut left = false;
                for attempt in 0..5 {
                    let at = format!("threads {threads} home {home} attempt {attempt}");
                    let n = 2000;
                    let ops = random_ops(&mut rng, n, &mut s, 300);
                    let expected = reference(&mut reference_hash, &ops, first, 40);
                    let cuts: Vec<usize> = (0..=40).map(|k| k * n / 40).collect();
                    let (got, stats, _) = partitioned_on(
                        &pool,
                        threads,
                        Some(Slow::OwnerFirstScan),
                        &mut hash,
                        &mut scratch,
                        &ops,
                        first,
                        &cuts,
                        threads,
                        2 * threads,
                        40,
                        false,
                        false,
                        home,
                    );
                    assert_eq!(got, expected, "{at}: hits");
                    assert_eq!(
                        hash_state(&hash),
                        hash_state(&reference_hash),
                        "{at}: table"
                    );
                    first += n;
                    if stats.helpers_left > 0 {
                        if home {
                            assert!(stats.away > 0, "{at}: no strip folded away from home");
                        }
                        left = true;
                        break;
                    }
                }
                assert!(
                    left,
                    "threads {threads} home {home}: no helper left a stalled stage in five stages"
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
        let expected = reference(&mut reference_hash, &ops, 0, 40);
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
            true,
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
        let unchanged = reference(&mut copy_hash(&reference_hash), &more, 2000, 40);
        let occupancy = hash.occupancy;
        for h in [&mut hash, &mut reference_hash] {
            for cell in &mut h.chain[1..occupancy as usize] {
                cell.level = i32::MAX / 2;
            }
        }
        let expected = reference(&mut reference_hash, &more, 2000, 40);
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
