//! EXPERIMENT (LOSAT_X_PAIRPARFOLD; a strict-set switch since 2026-10-17): the
//! diagonal-hash fold of the strips paths of PAIRPAR (blastn's one-hit fold, the
//! parent's `fold`) and `LOSAT_X_PAIRPARDISC` (dc-megablast's two-hit fold,
//! `disc::fold_two_hit`) split by diagonal-hash bucket over the threads of the search pool,
//! instead of done by the owner thread alone. No NCBI counterpart (scheduling only).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:361-381
//! (`s_BlastDiagHashRetrieve`)
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:396-449
//! (`s_BlastDiagHashInsert`)
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:825-841,853-858,939-941
//! (`s_BlastnDiagHashExtendInitialHit`)
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
//! EXPERIMENT (`seed_stage_part_lists`: LOSAT_X_DCSEEDLEAN for dc-megablast's stages,
//! LOSAT_X_BNSEEDLISTS for those of `x_pair_seed.rs` and `na`): the home loop's strips are
//! scanned straight into the partitions' lists (`PartSink`: each seed in the partition of the
//! cut, numbered in scan order, so the lists are those of the stable cut), and the folds do not
//! record the buckets they change, so the merge writes back every bucket of each partition (an
//! unchanged one has the stage's chain, its cells below the stage's occupancy, so writing it
//! back changes nothing).
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
use super::{FoldCounts, Resume, Seed, SeedOut};

// EXPERIMENT (LOSAT_X_SEEDPIECE): the home stage of `seed_stage_part_lists` in seed pieces of a
// fixed size (`x_pair_seed_piece.rs`).
#[path = "x_pair_seed_piece.rs"]
mod piece;
pub(super) use piece::PieceSink;
use std::cell::Cell;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, OnceLock, PoisonError, RwLock};

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
        // A stage on a pool has two threads or more; one thread in the dc-megablast stages
        // without a pool (LOSAT_X_BNDCONE).
        crate::utils::xreach::x_reach!("LOSAT_X_PAIRPARFOLD", threads >= 2);
    }
    (switch_on() && threads >= 2).then(|| threads.min(DIAGHASH_NUM_BUCKETS))
}

/// No NCBI counterpart: reads LOSAT_X_PAIRPARHOME once (a strict-set switch: on unless it is
/// `0`): the participants of a partitioned stage fold their home partitions (`seed_stage_home`)
/// instead of any partition they can lock. It does not change any value NCBI computes.
fn home_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_PAIRPARHOME"))
}

/// No NCBI counterpart: reads LOSAT_X_PAIRPARMERGE once (a strict-set switch: on unless it is
/// `0`): the merge of a partitioned stage puts the stage's hits in key order strip by strip, on a
/// pool in ranges of strips (`hits_by_strips`), instead of with one sort of all of them. It does
/// not change any value NCBI computes.
fn merge_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_PAIRPARMERGE"))
}

/// No NCBI counterpart (scheduling only, LOSAT_X_PAIRPARMERGE): the ranges of strips per thread of
/// the pool in `hits_by_strips` (more ranges than threads, so that a range with many hits does not
/// hold the others back).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
const MERGE_RANGES_PER_THREAD: usize = 2;

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:368
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
    /// EXPERIMENT (LOSAT_X_DCSEEDLEAN, LOSAT_X_BNSEEDLISTS): the stage's folds do not record the
    /// buckets they change
    /// (`touched`, `dirty`; `insert_t` without `TRACK`), and the merge writes back every bucket of
    /// the partition.
    all_dirty: bool,
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
/// A fold's insert (`DiagHashTable::insert_memo`; with `LEAN`, LOSAT_X_FOLDLEAN, the same statements
/// inlined, `insert_inl`). Without `MEMO` the table has no memo (LOSAT_X_DIAGLV2) and the insert
/// does not test for one: `insert_plain`, the statements before that switch (with `LEAN` inlined,
/// `insert_inl::<false>`). With partitions (`PT`) it also records
/// the diagonal's bucket and, when the insert appended a cell, the seed's key with the cell's
/// number in the partition's table.
#[cfg(all(
    test,
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
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
    insert_t::<PT, LEAN, MEMO, true>(
        hash,
        out,
        key,
        diag,
        level,
        hit_len,
        hit_saved,
        s_off,
        window_size,
    );
}

/// `insert`; EXPERIMENT (LOSAT_X_DCSEEDLEAN, LOSAT_X_BNSEEDLISTS): without `TRACK` a partition's
/// insert does not record the diagonal's bucket (`PartOut::all_dirty`).
#[allow(clippy::too_many_arguments)]
#[inline(always)]
pub(super) fn insert_t<const PT: bool, const LEAN: bool, const MEMO: bool, const TRACK: bool>(
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
        // a stage whose folds do not record their buckets merges every bucket
        debug_assert!(
            TRACK || out.all_dirty,
            "an untracked fold outside its stage"
        );
        let occupancy = hash.occupancy;
        if LEAN {
            hash.insert_inl::<MEMO>(diag, level, hit_len, hit_saved, s_off, window_size);
        } else if MEMO {
            hash.insert_memo(diag, level, hit_len, hit_saved, s_off, window_size);
        } else {
            hash.insert_plain(diag, level, hit_len, hit_saved, s_off, window_size);
        }
        if TRACK {
            let bucket = bucket_of_diag(diag);
            if !out.touched[bucket] {
                out.touched[bucket] = true;
                out.dirty.push(bucket as u16);
            }
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
/// with both switches off: memory-ordering machine clears x1.37, same instructions). With
/// LOSAT_X_PAIRPARHOME (`seed_stage_home`) also each partition's next strip, which its folder
/// writes after every strip and the others read.
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
                x_last: Vec::new(),
                #[cfg(debug_assertions)]
                x_last_fresh: true,
                x_map: super::super::XDiagMap::default(),
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
/// 1]]`, in scan order, with their positions in the strip's scan order in `pos`; or, where the
/// scan filled the partitions' lists (`split`, EXPERIMENT LOSAT_X_DCSEEDLEAN /
/// LOSAT_X_BNSEEDLISTS), `lists[p]`.
#[derive(Default)]
struct StripBuf {
    seeds: Vec<Seed>,
    pos: Vec<u32>,
    offs: Vec<u32>,
    lists: Vec<PartSeeds>,
    split: bool,
}

impl StripBuf {
    /// Partition p's seeds of the strip and their positions in the strip's scan order.
    #[inline(always)]
    fn part(&self, p: usize) -> (&[Seed], &[u32]) {
        if self.split {
            let list = &self.lists[p];
            (&list.seeds, &list.pos)
        } else {
            let (a, b) = (self.offs[p] as usize, self.offs[p + 1] as usize);
            (&self.seeds[a..b], &self.pos[a..b])
        }
    }
}

/// No NCBI counterpart (scheduling only, LOSAT_X_DCSEEDLEAN, LOSAT_X_BNSEEDLISTS): one
/// partition's seeds of a strip, in scan order, with their positions in the strip's scan order.
#[derive(Default)]
struct PartSeeds {
    seeds: Vec<Seed>,
    pos: Vec<u32>,
}

/// No NCBI counterpart (scheduling only, LOSAT_X_DCSEEDLEAN, LOSAT_X_BNSEEDLISTS): the receiver
/// of a strip's seeds that puts each seed in its partition's list (`table`: the partition of
/// each bucket, `bucket % parts`), numbering them in scan order, so the strip needs no cut: the
/// lists are those of the stable counting sort of `seed_stage_home`'s `produce`.
pub(super) struct PartSink<'a> {
    table: &'a [u16; DIAGHASH_NUM_BUCKETS],
    lists: &'a mut [PartSeeds],
    k: u32,
}

impl SeedOut for PartSink<'_> {
    #[inline(always)]
    fn push(&mut self, seed: Seed) {
        let p = self.table[bucket_of(&seed)] as usize;
        debug_assert!(p < self.lists.len());
        // SAFETY: the table holds `bucket % parts` and there are `parts` lists (`produce`).
        let list = unsafe { self.lists.get_unchecked_mut(p) };
        list.seeds.push(seed);
        list.pos.push(self.k);
        self.k += 1;
    }
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
    /// EXPERIMENT (LOSAT_X_SEEDPIECE): the chains of the pool's slots (allocated when first
    /// used) and the pool's seeds (`piece::Store`).
    arenas: Vec<piece::Links>,
    store: piece::Store,
    /// (tests) LOSAT_X_SEEDPIECE's stage as a test asks for it: `None` as the switch says,
    /// `Some(0)` the strips' lists, `Some(k)` pieces of `k` seeds.
    #[cfg(test)]
    pub(super) x_pieces: Option<usize>,
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
    /// the slowest partition, any partition once every strip was handed out; EXPERIMENT
    /// LOSAT_X_SEEDPIECE: the pieces of partitions so folded)
    pub(super) runs: u64,
    pub(super) away: u64,
    /// the partitions of strips of the stage (strips x partitions; EXPERIMENT LOSAT_X_SEEDPIECE:
    /// the pieces made, empty ones included, x partitions)
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

/// No NCBI counterpart (scheduling only): the turns without work in a row that a participant
/// spends spinning (`pause`) before it yields or sleeps; a helper watches the stage for a stall
/// (`Stall`) from the turn after them on.
const SPIN_TURNS: u32 = 64;

/// No NCBI counterpart (scheduling only): a helper's watch for a stalled stage, one form for the
/// three stages (`seed_stage_part`, `seed_stage_home`, `piece::seed_stage_pieces`). On every turn
/// without work after the `SPIN_TURNS` spinning ones (each such turn then yields or sleeps in
/// `pause`) the helper reads the stage's progress, and it leaves once the progress has not
/// changed for `HELPER_STALL` by the clock. It looks on every such turn, not on every 64th: a
/// turn that yields can take up to a scheduler quantum on a loaded machine (macOS's
/// `sched_yield`, libpthread's `swtch_pri` with 0, lowers the thread's priority for up to
/// `thread_depress_time`, one quantum of 10 ms: XNU osfmk/kern/syscall_subr.c `swtch_pri`,
/// osfmk/kern/sched_prim.c `sched_timeshare_timebase_init`), so 64 turns can take longer than
/// a stall, the next look finds the stage moving again, and a stall of many times
/// `HELPER_STALL` goes unseen. A look costs a clock read and the loads of the progress, beside
/// the yield or the 50 us sleep of its turn; turns with work do not look.
#[derive(Default)]
struct Stall(Option<(std::time::Instant, u64)>);

impl Stall {
    /// The participant did some work: the stage moves, so the watch starts again.
    #[inline]
    fn reset(&mut self) {
        self.0 = None;
    }

    /// After `idle` turns without work in a row: whether the stage, whose progress `progress`
    /// reads (only on the turns that look), has made no progress for `HELPER_STALL`.
    #[inline]
    fn leaves(&mut self, idle: u32, progress: impl FnOnce() -> u64) -> bool {
        if idle < SPIN_TURNS {
            return false;
        }
        let now = progress();
        match self.0 {
            Some((since, seen)) if seen == now => since.elapsed() >= HELPER_STALL,
            _ => {
                self.0 = Some((std::time::Instant::now(), now));
                false
            }
        }
    }
}

/// No NCBI counterpart: a participant's wait when it found nothing to do. It spins a little,
/// then yields; a native helper that keeps finding nothing sleeps 50 us between looks. On
/// wasm32 it only yields: `std::thread::sleep` panics under the browser's WASI shim.
#[inline]
fn pause(idle: u32, owner: bool) {
    if idle < SPIN_TURNS {
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
        home,
        n,
        threads,
        ahead,
        hash,
        hits,
        fold_part,
        scan,
        None::<&NoLists>,
        fold,
        timed,
    )
}

/// No NCBI counterpart (scheduling only): a stage without a scan into the partitions' lists.
type NoLists = fn(usize, &mut PartSink<'_>);

/// No NCBI counterpart (scheduling only): EXPERIMENT (LOSAT_X_DCSEEDLEAN, LOSAT_X_BNSEEDLISTS)
/// `seed_stage_part` of
/// a caller whose `fold` does not record the buckets it changes (`insert_t` without `TRACK`):
/// the merge writes back every bucket of each partition (`PartOut::all_dirty`; a bucket no fold
/// changed has the stage's chain, cells below the stage's occupancy, so writing it back changes
/// nothing). With LOSAT_X_PAIRPARHOME the strips are scanned by `scan_lists` into their
/// partitions' lists (`PartSink`) instead of by `scan` into one list that `produce` cuts.
/// EXPERIMENT (LOSAT_X_SEEDPIECE): with LOSAT_X_PAIRPARHOME the strips are scanned by
/// `scan_piece` in pieces of a fixed size (`piece::seed_stage_pieces`) instead.
#[allow(clippy::too_many_arguments)]
pub(super) fn seed_stage_part_lists<FS, FL, FP, FF>(
    n: usize,
    threads: usize,
    ahead: usize,
    hash: &mut DiagHashTable,
    hits: &mut Vec<UngappedHit>,
    fold_part: FoldPart<'_>,
    scan: &FS,
    scan_lists: &FL,
    scan_piece: &FP,
    fold: &FF,
    timed: bool,
) -> PartStats
where
    FS: Fn(usize, &mut Vec<Seed>) + Sync,
    FL: Fn(usize, &mut PartSink<'_>) + Sync,
    FP: Fn(usize, Option<Resume>, &mut PieceSink<'_>) -> Option<Resume> + Sync,
    FF: Fn(&mut PartState, &[Seed], &[u32]) + Sync,
{
    let home = home_on();
    if home {
        crate::utils::xreach::x_reach!("LOSAT_X_PAIRPARHOME", true);
    }
    #[cfg(test)]
    let pieces = match fold_part.scratch.x_pieces {
        None => piece::seeds(),
        Some(0) => None,
        k => k,
    };
    #[cfg(not(test))]
    let pieces = piece::seeds();
    if let (true, Some(piece_seeds)) = (home, pieces) {
        let FoldPart { scratch, parts } = fold_part;
        assert!(
            (1..=DIAGHASH_NUM_BUCKETS).contains(&parts),
            "LOSAT_X_PAIRPARFOLD: {parts} partitions"
        );
        rebuild(scratch, hash, parts, true);
        let occ0 = hash.occupancy;
        return piece::seed_stage_pieces(
            n,
            threads,
            ahead.max(1),
            hash,
            hits,
            scratch,
            parts,
            occ0,
            piece_seeds,
            scan_piece,
            fold,
            timed,
        );
    }
    seed_stage_in(
        home,
        n,
        threads,
        ahead,
        hash,
        hits,
        fold_part,
        scan,
        Some(scan_lists),
        fold,
        timed,
    )
}

/// No NCBI counterpart (scheduling only): `seed_stage_part` with the sharing of the folds given
/// (`home`: LOSAT_X_PAIRPARHOME, `seed_stage_home`; else any participant folds any partition it
/// can lock). EXPERIMENT (LOSAT_X_DCSEEDLEAN, LOSAT_X_BNSEEDLISTS): with `scan_lists`
/// (`seed_stage_part_lists`) the folds do not record their buckets, and the home stage scans
/// into the partitions' lists.
#[allow(clippy::too_many_arguments)]
fn seed_stage_in<FS, FL, FF>(
    home: bool,
    n: usize,
    threads: usize,
    ahead: usize,
    hash: &mut DiagHashTable,
    hits: &mut Vec<UngappedHit>,
    fold_part: FoldPart<'_>,
    scan: &FS,
    scan_lists: Option<&FL>,
    fold: &FF,
    timed: bool,
) -> PartStats
where
    FS: Fn(usize, &mut Vec<Seed>) + Sync,
    FL: Fn(usize, &mut PartSink<'_>) + Sync,
    FF: Fn(&mut PartState, &[Seed], &[u32]) + Sync,
{
    let FoldPart { scratch, parts } = fold_part;
    assert!(
        (1..=DIAGHASH_NUM_BUCKETS).contains(&parts),
        "LOSAT_X_PAIRPARFOLD: {parts} partitions"
    );
    let ahead = ahead.max(1);
    rebuild(scratch, hash, parts, scan_lists.is_some());
    let occ0 = hash.occupancy;
    if home {
        return seed_stage_home(
            n, threads, ahead, hash, hits, scratch, parts, occ0, scan, scan_lists, fold, timed,
        );
    }

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
        buf.split = false;
        let StripBuf {
            seeds, pos, offs, ..
        } = &mut buf;
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
        let mut stall = Stall::default();
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
                stall.reset();
                if let Some(t0) = wait_start.take() {
                    wait_ns += t0.elapsed().as_nanos() as u64;
                }
            } else {
                if timed && owner && wait_start.is_none() {
                    wait_start = Some(std::time::Instant::now());
                }
                idle = idle.saturating_add(1);
                // (the owner can finish the stage alone)
                if !owner && stall.leaves(idle, || progress.load(Ordering::Relaxed) as u64) {
                    left.fetch_add(1, Ordering::Relaxed);
                    failing.done = true;
                    lock(&raw_spare).push(raw);
                    return (scan_ns, wait_ns, pause_ns);
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
        n,
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
/// The owner alone can still scan every strip and fold every partition. EXPERIMENT
/// (LOSAT_X_DCSEEDLEAN, LOSAT_X_BNSEEDLISTS): with `scan_lists` a strip is scanned into its
/// partitions' lists (`PartSink`) instead of by `scan` into one list cut by partition (the same
/// seeds of each partition in the same order, with the same positions).
#[allow(clippy::too_many_arguments)]
fn seed_stage_home<FS, FL, FF>(
    n: usize,
    threads: usize,
    ahead: usize,
    hash: &mut DiagHashTable,
    hits: &mut Vec<UngappedHit>,
    scratch: &mut PartScratch,
    parts: usize,
    occ0: u32,
    scan: &FS,
    scan_lists: Option<&FL>,
    fold: &FF,
    timed: bool,
) -> PartStats
where
    FS: Fn(usize, &mut Vec<Seed>) + Sync,
    FL: Fn(usize, &mut PartSink<'_>) + Sync,
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
    let lanes: Vec<OwnLines<AtomicUsize>> =
        (0..parts).map(|_| OwnLines(AtomicUsize::new(0))).collect();
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

    // EXPERIMENT (LOSAT_X_DCSEEDLEAN, LOSAT_X_BNSEEDLISTS): the partition of each bucket, for
    // `PartSink`.
    let table: Option<[u16; DIAGHASH_NUM_BUCKETS]> =
        scan_lists.map(|_| std::array::from_fn(|b| (b % parts) as u16));
    // A strip's seeds, cut by partition with a stable counting sort (or scanned into the
    // partitions' lists).
    let produce = |i: usize, raw: &mut Vec<Seed>, fill: &mut Vec<u32>| {
        if scan_lists.is_none() {
            raw.clear();
            scan(i, raw);
        }
        // (no spare buffer: first the strips every partition has passed, which wait for the
        // next participant that reads the slowest partition)
        let mut buf = lock(&spare).pop();
        if buf.is_none() {
            slowest();
            buf = lock(&spare).pop();
        }
        let mut buf = buf.unwrap_or_default();
        if let (Some(scan_lists), Some(table)) = (scan_lists, &table) {
            buf.lists.resize_with(parts, PartSeeds::default);
            for list in &mut buf.lists {
                list.seeds.clear();
                list.pos.clear();
            }
            scan_lists(
                i,
                &mut PartSink {
                    table,
                    lists: &mut buf.lists,
                    k: 0,
                },
            );
            buf.split = true;
            *write(&slots[i]) = Some(buf);
            ready[i].store(true, Ordering::Release);
            return;
        }
        buf.split = false;
        let StripBuf {
            seeds, pos, offs, ..
        } = &mut buf;
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
        let mut stall = Stall::default();
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
                let (seeds, pos) = strip.part(p);
                if !seeds.is_empty() {
                    st.out.strip = i as u64;
                    let t0 = timed.then(std::time::Instant::now);
                    let pause0 = st.scan_pause_ns;
                    fold(st, seeds, pos);
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
                stall.reset();
                if let Some(t0) = wait_start.take() {
                    wait_ns += t0.elapsed().as_nanos() as u64;
                }
            } else {
                if timed && owner && wait_start.is_none() {
                    wait_start = Some(std::time::Instant::now());
                }
                idle = idle.saturating_add(1);
                // (the owner can finish the stage alone; the stage's progress: strips handed out
                // plus strips folded by each partition)
                if !owner
                    && stall.leaves(idle, || {
                        (cursor.load(Ordering::Relaxed)
                            + lanes
                                .iter()
                                .map(|lane| lane.0.load(Ordering::Relaxed))
                                .sum::<usize>()) as u64
                    })
                {
                    left.fetch_add(1, Ordering::Relaxed);
                    failing.done = true;
                    break;
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
        n,
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
fn rebuild(scratch: &mut PartScratch, hash: &DiagHashTable, parts: usize, all_dirty: bool) {
    scratch.parts.truncate(parts);
    while scratch.parts.len() < parts {
        scratch.parts.push(PartState::new());
    }
    for st in &mut scratch.parts {
        st.reset();
        st.out.all_dirty = all_dirty;
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
        t.x_last.clone_from(&hash.x_last);
        #[cfg(debug_assertions)]
        {
            t.x_last_fresh = hash.x_last_fresh;
        }
    }
}

/// No NCBI counterpart (scheduling only, LOSAT_X_PAIRPARMERGE): the stage's hits appended to
/// `hits` in the order of their keys, as `merge`'s sort puts them, without that sort: the keys
/// of a partition ascend (it folds the strips in order, and a strip's seeds of the partition in
/// their scan order, the counting sort of `produce` being stable; one hit per seed at most), and
/// a key is (strip, position in the strip's scan order), so the hits of strip i are one run per
/// partition, after those of the strips before i. Strip by strip, the runs are gathered (in
/// `runs`) and the few hits sorted (`hits_of_strips`). On a pool of two threads or more, the
/// strips of a stage of two strips or more are cut into ranges that the pool orders side by side,
/// each range writing the slice of the output after the hits of the ranges before it; otherwise
/// this thread goes through the strips. Every stage takes this path, whatever its number of hits:
/// the stages of a search are either large (the blastn and blastn-short folds, tens to hundreds
/// of thousands of hits, where the ranges on the pool take a third to a fifth of the sort's time)
/// or small (the two-hit fold of dc-megablast and megablast's, at most a few thousand hits,
/// where the pool's start costs tens of microseconds per stage, a ten-thousandth of the search).
/// Returns `false`, having written nothing, where the switch is off (`merge` sorts the hits then).
fn hits_by_strips(
    hits: &mut Vec<UngappedHit>,
    states: &[PartState],
    strips: usize,
    runs: &mut Vec<(u64, u32, u32)>,
) -> bool {
    if !merge_on() {
        return false;
    }
    // every stage takes the switch's path (once per stage)
    crate::utils::xreach::x_reach!("LOSAT_X_PAIRPARMERGE", true);
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    {
        // (the pool's size only on a pool: off one, `current_num_threads` would start rayon's
        // global pool)
        let threads = if rayon::current_thread_index().is_some() {
            rayon::current_num_threads()
        } else {
            1
        };
        if threads >= 2 && strips >= 2 {
            hits_by_strips_on_pool(hits, states, strips, threads * MERGE_RANGES_PER_THREAD);
            return true;
        }
    }
    hits_by_strips_here(hits, states, strips, runs);
    true
}

/// No NCBI counterpart (LOSAT_X_PAIRPARMERGE): `hits_by_strips` on this thread, the strips in
/// order from each partition's first hit.
fn hits_by_strips_here(
    hits: &mut Vec<UngappedHit>,
    states: &[PartState],
    strips: usize,
    runs: &mut Vec<(u64, u32, u32)>,
) {
    let mut cur = vec![0usize; states.len()];
    hits.reserve(states.iter().map(|st| st.hits.len()).sum());
    hits_of_strips(states, 0..strips, &mut cur, runs, |hit| hits.push(hit));
    // every key is of one of the stage's strips
    assert!(cur.iter().zip(states).all(|(&c, st)| c == st.hits.len()));
}

/// No NCBI counterpart (LOSAT_X_PAIRPARMERGE): `hits_by_strips` on the current pool with the
/// strips cut into `ranges` ranges of equal strip counts (at most one per strip). A range starts
/// each partition's cursor at its first key of the range's first strip (binary search) and
/// writes its hits into its own slice of the spare capacity of `hits`.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn hits_by_strips_on_pool(
    hits: &mut Vec<UngappedHit>,
    states: &[PartState],
    strips: usize,
    ranges: usize,
) {
    let ranges = ranges.clamp(1, strips.max(1));
    let bounds: Vec<usize> = (0..=ranges).map(|r| r * strips / ranges).collect();
    // per range boundary, each partition's first hit at or after it, and the hits before it
    let starts: Vec<Vec<usize>> = bounds
        .iter()
        .map(|&b| {
            let key = (b as u64) << 32;
            states
                .iter()
                .map(|st| {
                    debug_assert_eq!(st.hits.len(), st.out.hit_keys.len());
                    st.out.hit_keys.partition_point(|&k| k < key)
                })
                .collect()
        })
        .collect();
    let before: Vec<usize> = starts.iter().map(|c| c.iter().sum()).collect();
    let total = before[ranges];
    // every key is of one of the stage's strips
    assert_eq!(total, states.iter().map(|st| st.hits.len()).sum::<usize>());
    let len0 = hits.len();
    hits.reserve(total);
    let mut rest = &mut hits.spare_capacity_mut()[..total];
    let mut slices = Vec::with_capacity(ranges);
    for r in 0..ranges {
        let (slice, tail) = rest.split_at_mut(before[r + 1] - before[r]);
        slices.push(slice);
        rest = tail;
    }
    rayon::scope(|scope| {
        for (r, slice) in slices.into_iter().enumerate() {
            let (bounds, starts) = (&bounds, &starts);
            scope.spawn(move |_| {
                let mut cur = starts[r].clone();
                let mut runs = Vec::new();
                let mut k = 0usize;
                hits_of_strips(
                    states,
                    bounds[r]..bounds[r + 1],
                    &mut cur,
                    &mut runs,
                    |hit| {
                        slice[k].write(hit);
                        k += 1;
                    },
                );
                // the range wrote every element of its slice, and only those
                assert!(k == slice.len() && cur == starts[r + 1]);
            });
        }
    });
    // SAFETY: the slices cover `len0 .. len0 + total` of the allocation (`reserve(total)`), and
    // every range wrote each element of its slice (the assert above; a range that panics ends
    // the scope with its panic before this line).
    unsafe { hits.set_len(len0 + total) };
}

/// No NCBI counterpart (LOSAT_X_PAIRPARMERGE): the hits of strips `strips` in key order, passed
/// to `put`: for each strip, every partition's run of keys below the next strip's from `cur[p]`
/// on (which moves past it) is gathered in `runs` and sorted (`hits_by_strips`).
fn hits_of_strips(
    states: &[PartState],
    strips: std::ops::Range<usize>,
    cur: &mut [usize],
    runs: &mut Vec<(u64, u32, u32)>,
    mut put: impl FnMut(UngappedHit),
) {
    for strip in strips {
        let end = ((strip as u64) + 1) << 32;
        runs.clear();
        for (p, st) in states.iter().enumerate() {
            let keys = &st.out.hit_keys;
            let mut j = cur[p];
            while j < keys.len() && keys[j] < end {
                debug_assert!(
                    j == cur[p] || keys[j - 1] < keys[j],
                    "a partition's keys ascend"
                );
                runs.push((keys[j], p as u32, j as u32));
                j += 1;
            }
            cur[p] = j;
        }
        if runs.len() > 1 {
            runs.sort_unstable();
        }
        for &(_, p, j) in runs.iter() {
            put(states[p as usize].hits[j as usize]);
        }
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:428-447
/// The end of a partitioned stage of `strips` strips (module comment): the hits in key order
/// after `hits`; the appended cells numbered `occ0 + rank` of their keys; the chains of the
/// changed buckets of every partition (of all its buckets where its folds did not record them,
/// `PartOut::all_dirty`) written into `hash` with those numbers (each cell at its
/// number, its link and the bucket's head mapped; the partition tables themselves are left as
/// they are); the capacity of the reference's doublings.
fn merge(
    hash: &mut DiagHashTable,
    hits: &mut Vec<UngappedHit>,
    states: &mut [PartState],
    occ0: u32,
    strips: usize,
    order: &mut Vec<(u64, u32, u32)>,
    apps: &mut Vec<(u64, u32, u32)>,
) {
    let parts = states.len();
    // EXPERIMENT (LOSAT_X_PAIRPARMERGE): the same order strip by strip (on the pool in ranges).
    if !hits_by_strips(hits, states, strips, order) {
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
        let mut put = |b: usize| {
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
        };
        if out.all_dirty {
            // EXPERIMENT (LOSAT_X_DCSEEDLEAN, LOSAT_X_BNSEEDLISTS): every bucket of the partition
            // (an unchanged one has the stage's chain, its cells below occ0, so writing it back
            // changes nothing).
            for b in (p..DIAGHASH_NUM_BUCKETS).step_by(parts) {
                put(b);
            }
        } else {
            for &bucket in &out.dirty {
                put(bucket as usize);
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

/// No NCBI counterpart: the fields the LOSAT_TIMING lines of `kind` end with when the switch is
/// on: the partition count of the stages, the number of partitioned stages merged, the time in
/// merges and the number of helpers that left a stalled stage; with LOSAT_X_PAIRPARHOME also the
/// runs folded and the partitions of strips folded away from home (`PartStats::runs`, `away`), and
/// the partitions of strips of the stages (`units`: strips x partitions; with LOSAT_X_SEEDPIECE pieces x
/// partitions). Without the
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
#[path = "x_pair_seed_part_tests.rs"]
mod tests;
