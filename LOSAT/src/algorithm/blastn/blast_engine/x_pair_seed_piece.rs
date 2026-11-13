//! EXPERIMENT (LOSAT_X_SEEDPIECE; a strict-set switch): the partitioned seed stage of the strips
//! paths (`x_pair_seed_part.rs` `seed_stage_part_lists` with LOSAT_X_PAIRPARHOME: dc-megablast's
//! two-hit stage, and the stages of `x_pair_seed.rs` and `na` with LOSAT_X_BNSEEDLISTS) holds a
//! strip's seeds in pieces of at most `PIECE_SEEDS` seeds, each in an arena of a fixed size,
//! instead of in one growable list per partition and strip. No NCBI counterpart (scheduling and
//! memory only).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684
//! ```c
//!     while(s_DetermineScanningOffsets(subject, word_length, lut_word_length, scan_range)) {
//!
//!         hitsfound = scansub(lookup_wrap, subject, offset_pairs, max_hits, &scan_range[1]);
//! ...
//!         hits_extended += extend(offset_pairs, hitsfound, word_params,
//!                                 lookup_wrap, query, subject, matrix,
//!                                 query_info, ewp, init_hitlist, scan_range[2] + lut_word_length);
//!     }
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1447-1454 (`MB_ACCESS_HITS`)
//! ```c
//!    if (PV_TEST(pv, index, pv_array_bts)) {                         \
//!       if (total_hits >= max_hits) {                                 \
//!          if ((Int4)(s_off) > scan_range[0])                         \
//!             scan_range[0] = s_off;                                  \
//!          break;                                                     \
//!       }                                                             \
//! ```
//! NCBI scans the subject into an offset array of `OFFSET_ARRAY_SIZE + longest_chain` pairs per
//! thread, stops when it is full and extends those pairs before the scan goes on from where it
//! stopped; only the diagonal hash carries over from one batch to the next. The strips stage
//! scans strips on any thread and folds every partition's seeds of a strip after those of the
//! strips before, so a strip's seeds wait until every partition has folded them; their lists
//! kept the largest capacity any strip had needed (`Vec` doubling, the buffers recycled), and on
//! dc-megablast's 5 Mb chunks at 16 threads that held about 80 MiB, nearly all of it capacity of
//! partitions that had once held a strip's main-diagonal run (results/s12c8 of line EB).
//!
//! The form here: a strip's scan stops at the first pair a piece has no room for (`Cut`; inside a
//! word's run as well, so no word makes a piece larger) and goes on from that word in a later
//! piece, on any thread (the stop is data, `Resume`, parked with the strip), as NCBI's scan goes
//! on from `scan_range[0]`. A piece is an arena of `PIECE_SEEDS / BLOCK + parts` blocks of
//! `BLOCK` seeds, enough for `PIECE_SEEDS` seeds whatever partitions they go to; each partition's
//! seeds of the piece are a chain of blocks in scan order, taken from a lane of the arena of its
//! own while its share lasts (`Links`), so that they lie in one run. The arenas are a pool of
//! `PIECES_PER_THREAD * threads + RESERVED` slots, allocated at their size once and kept in the
//! subject scratch, so a stage holds at most that many pieces whatever the subject, the query or
//! the density of the seeds (and arenas only for as many pieces as it had at once: a piece takes
//! the slot of one every partition has folded before a slot no piece has used, `Sched::take`, as
//! the strips stage takes a new buffer only when it has none to give back); the strip the
//! slowest partition is in (the lowest strip a partition has not folded all of) alone may take
//! the last `RESERVED` slots, so that it can always go on (no wait for a piece
//! nobody frees). A participant produces one piece per turn and parks the rest of the strip, so
//! that a helper never holds a strip half-scanned (it can leave a stalled stage at any turn).
//!
//! Why the output is the same: a partition folds the pieces of a strip in order and the strips
//! in order, a piece's seeds of the partition in scan order (the blocks of its chain), so each
//! partition folds the seeds of the reference order (`x_pair_seed_part.rs` module comment); the
//! seeds' positions in the strip's scan order count on across the pieces (`PieceSink::k`), so
//! the keys (strip, position), the order of the hits and the numbers of the appended cells are
//! those of the strips stage; the scan gives the same pairs in the same order (`Cut`: a scan
//! started again at a word gives the words from that one on, and the first pairs of that word
//! that an earlier piece took are passed over); the fold's statements do not change.

use super::super::super::{DiagHashTable, UngappedHit, DIAGHASH_NUM_BUCKETS};
use super::super::{Resume, Seed, SeedOut};
use super::{
    bucket_of, lock, merge, pause, read, run_participants, write, AtHome, Failing, InStage,
    OwnLines, PartScratch, PartState, PartStats, Run, HELPER_STALL, IN_STAGE,
};
use std::cell::Cell;
use std::mem::MaybeUninit;
use std::sync::atomic::{AtomicBool, AtomicU32, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Mutex, OnceLock, PoisonError, RwLock};

/// No NCBI counterpart: reads LOSAT_X_SEEDPIECE once (a strict-set switch that carries a number:
/// on unless it is `0`, with pieces of the number of seeds given from 2 to 2^17, else of
/// `PIECE_SEEDS`): the seeds of a piece, `None` when the switch is off. It does not change any
/// value NCBI computes.
pub(super) fn seeds() -> Option<usize> {
    static N: OnceLock<Option<usize>> = OnceLock::new();
    *N.get_or_init(|| crate::utils::xswitch::value("LOSAT_X_SEEDPIECE", 2..=1 << 17, PIECE_SEEDS))
}

/// The seeds of one piece at most (NCBI's `OFFSET_ARRAY_SIZE`, 4,096 pairs per batch, times 4).
pub(super) const PIECE_SEEDS: usize = 16_384;
/// The seeds of a block of a piece (a partition's seeds of a piece are a chain of blocks).
const BLOCK: usize = 64;
/// The pieces of a stage per thread of the stage (the pool of arenas).
pub(super) const PIECES_PER_THREAD: usize = 3;
/// The slots of the pool only the strip the slowest partition is in may take.
pub(super) const RESERVED: usize = 2;
/// No block (the end of a chain, or a partition without a block yet).
const NONE: u16 = u16::MAX;
/// `StripState::published`: the strip's scan is done (the count below it is final).
const DONE: u32 = 1 << 31;

/// One partition's seeds of a piece: the first and the last block of its chain and the seeds in
/// the last one (`BLOCK` before its first seed, so that the first push takes a block).
#[derive(Clone, Copy)]
struct Chain {
    head: u16,
    tail: u16,
    fill: u16,
}

const EMPTY: Chain = Chain {
    head: NONE,
    tail: NONE,
    fill: BLOCK as u16,
};

/// A partition's lane of an arena: its next block and the end of its lane (`Links::take`).
#[derive(Clone, Copy)]
struct Lane {
    next: u16,
    end: u16,
}

/// No NCBI counterpart (scheduling only): the chains of a slot's piece: its blocks of `BLOCK`
/// seeds (`used` of them taken), each partition's blocks linked in its order (`next`, its first
/// and last block in `chains`). A partition takes the blocks of its lane first (`lanes`: the
/// arena's first `parts` runs of equal length, one per partition), then the blocks after the
/// lanes (`over`), then the last free block of another partition's lane (`steal`), so that a
/// partition's seeds of a piece lie in one run of blocks unless it has more than its share, and
/// any `PIECE_SEEDS` seeds still fit. (Blocks taken in the order of the pushes, the partitions'
/// blocks interleaved, made the fold read each partition's seeds in blocks scattered over the
/// arena: the hardware prefetches of the neighbouring blocks of other partitions and the pages
/// of the scattered blocks took the fold +6 % cycles at the same instructions, L2 misses +22 %;
/// dc-megablast cDNA x genome at 4 threads, +2.2 % over the strips stage, line EB s15d / s15e.)
/// Kept in the subject scratch from one stage to the next; the seeds are in the pool's `Store`.
#[derive(Default)]
pub(super) struct Links {
    next: Vec<u16>,
    chains: Vec<Chain>,
    lanes: Vec<Lane>,
    over: u16,
    steal: u16,
    used: u16,
}

impl Links {
    /// A free block for partition p: the next one of its lane, else the next one after the
    /// lanes, else the last free one of another lane.
    #[inline(always)]
    fn take(&mut self, p: usize) -> u16 {
        debug_assert!(p < self.lanes.len());
        // SAFETY: `p` is below `parts`, the number of lanes (`Arena::reset`, `PieceSink::push`).
        let lane = unsafe { self.lanes.get_unchecked_mut(p) };
        let b = if lane.next < lane.end {
            lane.next += 1;
            lane.next - 1
        } else {
            self.take_other()
        };
        self.used += 1;
        b
    }

    /// Empty, for `parts` partitions, pieces of `piece_seeds` seeds and arenas of `blocks` blocks
    /// (`Arena::blocks`; `next` allocated at its size the first time): equal lanes of the full
    /// blocks a piece fills, and the rest after them.
    fn reset(&mut self, parts: usize, piece_seeds: usize, blocks: usize) {
        if self.next.len() != blocks {
            self.next = vec![NONE; blocks];
        }
        self.chains.clear();
        self.chains.resize(parts, EMPTY);
        let lane = piece_seeds / BLOCK / parts;
        self.lanes.clear();
        self.lanes.extend((0..parts).map(|p| Lane {
            next: (p * lane) as u16,
            end: ((p + 1) * lane) as u16,
        }));
        self.over = (parts * lane) as u16;
        self.steal = 0;
        self.used = 0;
    }

    /// `take` once partition p's lane is used up.
    #[cold]
    #[inline(never)]
    fn take_other(&mut self) -> u16 {
        if usize::from(self.over) < self.next.len() {
            self.over += 1;
            return self.over - 1;
        }
        let parts = self.lanes.len();
        for i in 0..parts {
            let q = (usize::from(self.steal) + i) % parts;
            let lane = &mut self.lanes[q];
            if lane.next < lane.end {
                lane.end -= 1;
                self.steal = q as u16;
                return lane.end;
            }
        }
        // (a piece takes at most `PIECE_SEEDS` seeds, `Cut`, and a partition's last block is the
        // only one not full, so the arena's `PIECE_SEEDS / BLOCK + parts` blocks hold them)
        panic!("LOSAT_X_SEEDPIECE: a full arena");
    }

    /// (tests) whether a stage used the slot.
    #[cfg(all(
        test,
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    pub(super) fn allocated(&self) -> bool {
        !self.next.is_empty()
    }
}

/// No NCBI counterpart (scheduling only): the seeds and their positions in the strip's scan
/// order of every slot of the pool, slot `s` at `s * len..(s + 1) * len` (`len`: an arena's
/// blocks of `BLOCK` seeds), kept in the subject scratch from one stage to the next.
///
/// Each buffer is one allocation whose capacity is above glibc's ceiling of the mmap threshold
/// (`utils/x_mapped.rs` `x_unthresholded_capacity`), not initialized: glibc maps it, its free
/// leaves the threshold as it is, and only the pages of the blocks pieces wrote are resident; a
/// dropped store's buffers wait for the next store (`KEPT`). (An arena allocated by itself, 266-328 KiB of seeds, was mapped and its free raised glibc's
/// threshold to its size, so that later blocks of up to that size came from the arenas of the
/// heap and stayed resident there: BLASTN cDNA x genome -subject at 16 threads, peak x1.13-1.15
/// of the strips stage, x1.04 with the threshold held at 128 KiB, line EB s14u / s14w.)
#[derive(Default)]
pub(super) struct Store {
    seeds: Vec<MaybeUninit<Seed>>,
    pos: Vec<MaybeUninit<u32>>,
    /// the seeds of a slot's region in the layout of the buffers (0 before the first stage)
    len: usize,
}

/// No NCBI counterpart (scheduling only): the buffers of a dropped `Store` (`KEPT`).
struct StoreBufs {
    seeds: Vec<MaybeUninit<Seed>>,
    pos: Vec<MaybeUninit<u32>>,
    len: usize,
}

/// No NCBI counterpart (scheduling only): the buffers of the `Store`s dropped so far, for the
/// next stores' first stages (`Store::slots`), as glibc keeps the heap's freed blocks for the
/// next allocations. A scratch lives for one query batch (`run.rs` `SubjectScratch`), and the
/// pages of a mapped buffer go back to the system when it is freed, so each batch's stages wrote
/// the store's pages anew (BLASTN cDNA x genome -subject at 16 threads: 66,000 page faults more
/// than the strips stage, whose buffers come from the heap; line EB s15). A store takes one when
/// it has none, of its layout if one is there, so it holds no more buffers than were alive at
/// once.
static KEPT: Mutex<Vec<StoreBufs>> = Mutex::new(Vec::new());

impl Drop for Store {
    fn drop(&mut self) {
        if self.len == 0 {
            return;
        }
        let bufs = StoreBufs {
            seeds: std::mem::take(&mut self.seeds),
            pos: std::mem::take(&mut self.pos),
            len: self.len,
        };
        lock(&KEPT).push(bufs);
    }
}

impl Store {
    /// The regions of `slots` slots of `len` seeds each (the buffers allocated anew when they
    /// are too small; a store without buffers first takes those of one dropped before, `KEPT`).
    #[allow(clippy::type_complexity)]
    fn slots(
        &mut self,
        slots: usize,
        len: usize,
    ) -> impl Iterator<Item = (&mut [MaybeUninit<Seed>], &mut [MaybeUninit<u32>])> {
        let need = slots * len;
        if self.len == 0 {
            let mut kept = lock(&KEPT);
            let at = kept.iter().rposition(|b| b.len == len);
            let bufs = match at {
                Some(at) => Some(kept.swap_remove(at)),
                None => kept.pop(),
            };
            drop(kept);
            if let Some(b) = bufs {
                self.seeds = b.seeds;
                self.pos = b.pos;
                self.len = b.len;
            }
        }
        // (a stage of another layout allocates anew, so that the pages the regions of the old one
        // touched are not kept; one of the same layout with more slots takes the capacity the
        // buffers have, above glibc's ceiling)
        if self.len != len || self.seeds.capacity() < need || self.pos.capacity() < need {
            // (the old buffers go first, so that the two are not held at once)
            self.seeds = Vec::new();
            self.pos = Vec::new();
            self.seeds = uninit(need);
            self.pos = uninit(need);
            self.len = len;
        } else if self.seeds.len() < need {
            // SAFETY: `MaybeUninit` needs no initialization, and both capacities are at least
            // `need` (the test above).
            unsafe {
                self.seeds.set_len(need);
                self.pos.set_len(need);
            }
        }
        self.seeds[..need]
            .chunks_mut(len)
            .zip(self.pos[..need].chunks_mut(len))
    }
}

/// No NCBI counterpart (scheduling only): one slot of the pool in a stage: one piece's seeds, cut
/// by partition, in its region of the `Store` (blocks of `BLOCK` seeds and their positions) and
/// the chains of its blocks (`Links`). Its size is set for the partition count (`reset`) and
/// does not change with the seeds.
pub(super) struct Arena<'s> {
    // (not initialized, so that only the pages of the blocks a piece takes become resident: a
    // stage of few seeds holds few pages whatever the arena's size; a slot is read only after a
    // push of the same piece wrote it, `for_each_run`)
    seeds: &'s mut [MaybeUninit<Seed>],
    pos: &'s mut [MaybeUninit<u32>],
    links: Links,
}

impl Arena<'_> {
    /// The blocks of an arena for `parts` partitions: `piece_seeds` seeds in any cut, since a
    /// partition's last block is the only one not full.
    fn blocks(parts: usize, piece_seeds: usize) -> usize {
        piece_seeds.div_ceil(BLOCK) + parts
    }

    /// Empty, for `parts` partitions and pieces of `piece_seeds` seeds (its region of the
    /// `Store` is `blocks * BLOCK` seeds; its links allocated at their size the first time).
    fn reset(&mut self, parts: usize, piece_seeds: usize) {
        let blocks = Self::blocks(parts, piece_seeds);
        assert!(
            self.seeds.len() == blocks * BLOCK && self.pos.len() == blocks * BLOCK,
            "LOSAT_X_SEEDPIECE: an arena's region of the store"
        );
        self.links.reset(parts, piece_seeds, blocks);
    }

    /// Partition p's seeds of the piece in scan order, with their positions, through `f`: one
    /// call per run of consecutive blocks of its chain.
    #[inline]
    fn for_each_run(&self, p: usize, mut f: impl FnMut(&[Seed], &[u32])) {
        let chain = self.links.chains[p];
        if chain.head == NONE {
            return;
        }
        let tail = chain.tail as usize;
        let mut b = chain.head as usize;
        let mut first = b;
        loop {
            if b == tail {
                let (a, z) = (first * BLOCK, b * BLOCK + chain.fill as usize);
                // SAFETY: blocks `first..=b` are p's, taken in this piece and filled in order
                // (all of each but the last, `fill` of it), so every slot of `a..z` was written.
                let (seeds, pos) =
                    unsafe { (written(&self.seeds[a..z]), written(&self.pos[a..z])) };
                f(seeds, pos);
                return;
            }
            let next = self.links.next[b] as usize;
            if next != b + 1 {
                let (a, z) = (first * BLOCK, (b + 1) * BLOCK);
                // SAFETY: as above; blocks `first..=b` are full.
                let (seeds, pos) =
                    unsafe { (written(&self.seeds[a..z]), written(&self.pos[a..z])) };
                f(seeds, pos);
                first = next;
            }
            b = next;
        }
    }
}

/// `n` slots, not initialized (allocated, not written: the pages stay untouched until used), in
/// a buffer whose free leaves glibc's mmap threshold as it is (`Store`).
fn uninit<T>(n: usize) -> Vec<MaybeUninit<T>> {
    let mut v = Vec::with_capacity(crate::utils::x_mapped::x_unthresholded_capacity(
        n,
        std::mem::size_of::<T>(),
    ));
    // SAFETY: `MaybeUninit<T>` needs no initialization, and the capacity is at least `n`.
    unsafe { v.set_len(n) };
    v
}

/// The slots `s` as values.
///
/// # Safety
/// Every slot of `s` was written.
unsafe fn written<T>(s: &[MaybeUninit<T>]) -> &[T] {
    // SAFETY: `MaybeUninit<T>` has the layout of `T`, and the caller says every slot holds a value.
    unsafe { std::slice::from_raw_parts(s.as_ptr().cast::<T>(), s.len()) }
}

/// No NCBI counterpart (scheduling only): the receiver of a piece's seeds: each seed to the
/// block chain of its partition (`table`: the partition of each bucket, `bucket % parts`), with
/// its position in the strip's scan order (`k`, counted on from the pieces before), until
/// `limit` (`room`, which `Cut` reads).
pub(in super::super) struct PieceSink<'a> {
    table: &'a [u16; DIAGHASH_NUM_BUCKETS],
    seeds: &'a mut [MaybeUninit<Seed>],
    pos: &'a mut [MaybeUninit<u32>],
    links: &'a mut Links,
    k: u32,
    limit: u32,
}

impl SeedOut for PieceSink<'_> {
    #[inline(always)]
    fn push(&mut self, seed: Seed) {
        let p = self.table[bucket_of(&seed)] as usize;
        let a = &mut *self.links;
        debug_assert!(p < a.chains.len());
        // SAFETY: the table holds `bucket % parts` and there are `parts` chains and lanes
        // (`reset`).
        if unsafe { a.chains.get_unchecked(p) }.fill as usize == BLOCK {
            // (a block of the arena, below `next.len()`: `take` panics when none is free, which
            // keeps the stores below in the arena whatever the caller does)
            let b = a.take(p);
            a.next[b as usize] = NONE;
            // SAFETY: as above.
            let chain = unsafe { a.chains.get_unchecked_mut(p) };
            if chain.head == NONE {
                chain.head = b;
            } else {
                a.next[chain.tail as usize] = b;
            }
            chain.tail = b;
            chain.fill = 0;
        }
        // SAFETY: as above.
        let chain = unsafe { a.chains.get_unchecked_mut(p) };
        let at = chain.tail as usize * BLOCK + chain.fill as usize;
        debug_assert!(at < self.seeds.len());
        // SAFETY: `chain.tail` is a block below `next.len()` (`take`, its index into `next` above
        // when it was taken) and `fill < BLOCK`, so `at` is below `blocks * BLOCK`, the length of
        // both regions (`Arena::reset`).
        unsafe {
            self.seeds.get_unchecked_mut(at).write(seed);
            self.pos.get_unchecked_mut(at).write(self.k);
        }
        chain.fill += 1;
        self.k += 1;
    }

    #[inline(always)]
    fn room(&self) -> usize {
        (self.limit - self.k) as usize
    }

    #[inline(always)]
    fn close(&mut self) {
        self.limit = self.k;
    }
}

/// No NCBI counterpart (scheduling only): a strip's pieces: their slots in the pool, in order,
/// and how many there are, with `DONE` once the strip's scan is done.
#[derive(Default)]
struct StripState {
    pieces: Mutex<Vec<u32>>,
    published: AtomicU32,
}

/// No NCBI counterpart (scheduling only): a strip whose scan stopped at the end of a piece:
/// where it goes on, the next seed's position and the pieces it has.
#[derive(Clone, Copy)]
struct Parked {
    strip: usize,
    from: Resume,
    k: u32,
}

/// No NCBI counterpart (scheduling only): what the participants hand out under one lock: the
/// parked strips, the free slots of the pool (those given back, `free`, and those no piece has
/// used yet, `fresh` up to the pool's size) and the next strip nobody has scanned.
struct Sched {
    parked: Vec<Parked>,
    free: Vec<u32>,
    fresh: u32,
    cursor: usize,
}

impl Sched {
    /// The free slots of a pool of `slots_n`.
    fn avail(&self, slots_n: usize) -> usize {
        self.free.len() + slots_n - self.fresh as usize
    }

    /// A free slot: one given back (the last one first), else the next one no piece has used, so
    /// that a stage holds as many arenas as it ever had pieces at once (as the strips stage
    /// takes a new buffer only when no strip it can give back is left).
    fn take(&mut self) -> u32 {
        self.free.pop().unwrap_or_else(|| {
            self.fresh += 1;
            self.fresh - 1
        })
    }
}

/// The lowest parked strip (`usize::MAX` if none).
fn lowest_parked(parked: &[Parked]) -> usize {
    parked.iter().map(|p| p.strip).min().unwrap_or(usize::MAX)
}

/// A partition's next piece: (strip << 32) | piece.
#[inline(always)]
fn pack(strip: usize, piece: u32) -> u64 {
    ((strip as u64) << 32) | u64::from(piece)
}

#[inline(always)]
fn strip_of(lane: u64) -> usize {
    (lane >> 32) as usize
}

/// No NCBI counterpart (scheduling only): `seed_stage_part_lists`'s stage (after `rebuild`) in
/// pieces (module comment). `scan` scans strip i from the start (`None`) or from where a piece
/// stopped into a piece's sink and returns where the strip goes on (`None`: done); `fold` folds a
/// run of a partition's seeds of a piece (`PartOut::strip` set). The folds are shared as in
/// `seed_stage_home` (home partitions, `Run`), with a piece where that stage has a strip: each
/// partition folds its pieces in order, one folder at a time (its lock), so the output is the
/// same. `ahead` bounds the strips handed out past the slowest partition's (as there), the pool
/// the pieces.
#[allow(clippy::too_many_arguments)]
pub(super) fn seed_stage_pieces<FP, FF>(
    n: usize,
    threads: usize,
    ahead: usize,
    hash: &mut DiagHashTable,
    hits: &mut Vec<UngappedHit>,
    scratch: &mut PartScratch,
    parts: usize,
    occ0: u32,
    piece_seeds: usize,
    scan: &FP,
    fold: &FF,
    timed: bool,
) -> PartStats
where
    FP: Fn(usize, Option<Resume>, &mut PieceSink<'_>) -> Option<Resume> + Sync,
    FF: Fn(&mut PartState, &[Seed], &[u32]) + Sync,
{
    crate::utils::xreach::x_reach!("LOSAT_X_SEEDPIECE", true);
    let threads = threads.max(1);
    let slots_n = PIECES_PER_THREAD * threads + RESERVED;
    let states: Vec<Mutex<OwnLines<PartState>>> = scratch
        .parts
        .drain(..)
        .map(|st| Mutex::new(OwnLines(st)))
        .collect();
    let mut links = std::mem::take(&mut scratch.arenas);
    // (a stage of fewer threads than the one before keeps only its own slots' links)
    links.truncate(slots_n);
    links.resize_with(slots_n, Links::default);
    let mut store = std::mem::take(&mut scratch.store);
    let pool: Vec<RwLock<Arena<'_>>> = store
        .slots(slots_n, Arena::blocks(parts, piece_seeds) * BLOCK)
        .zip(links)
        .map(|((seeds, pos), links)| RwLock::new(Arena { seeds, pos, links }))
        .collect();
    let sched = Mutex::new(Sched {
        parked: Vec::new(),
        free: Vec::with_capacity(slots_n),
        fresh: 0,
        cursor: 0,
    });
    // (hints read without the lock: the parked strips, the free slots, the slots given back, the
    // next strip)
    let parked_n = AtomicUsize::new(0);
    let free_n = AtomicUsize::new(slots_n);
    let given_n = AtomicUsize::new(0);
    let cursor_hint = AtomicUsize::new(0);
    let low_parked = AtomicUsize::new(usize::MAX);
    let lanes: Vec<OwnLines<AtomicU64>> = (0..parts).map(|_| OwnLines(AtomicU64::new(0))).collect();
    let homes: Vec<AtomicBool> = (0..parts).map(|_| AtomicBool::new(false)).collect();
    let strips: Vec<StripState> = (0..n).map(|_| StripState::default()).collect();
    // the pieces below it have gone back to the pool (every partition folded them)
    let freed = AtomicU64::new(0);
    let produced = AtomicU64::new(0);
    let failed = AtomicBool::new(false);
    let left = AtomicUsize::new(0);
    let runs = AtomicU64::new(0);
    let away = AtomicU64::new(0);
    let table: [u16; DIAGHASH_NUM_BUCKETS] = std::array::from_fn(|b| (b % parts) as u16);

    // The slowest partition's next piece (every partition has folded the pieces below it); the
    // pieces below it go back to the pool, by the participant that moves `freed`.
    let slowest = || -> u64 {
        let min_next = lanes
            .iter()
            .map(|lane| lane.0.load(Ordering::Acquire))
            .min()
            .unwrap_or(pack(n, 0));
        // (a load first, as in `seed_stage_home`; lock order: the schedule, then a strip's pieces)
        if freed.load(Ordering::Relaxed) < min_next {
            let mut s = lock(&sched);
            let f = freed.load(Ordering::Relaxed);
            if f < min_next {
                freed.store(min_next, Ordering::Relaxed);
                let (fs, fk) = (strip_of(f), f as u32);
                let (ms, mk) = (strip_of(min_next), min_next as u32);
                for (i, strip) in strips.iter().enumerate().take((ms + 1).min(n)).skip(fs) {
                    let pieces = lock(&strip.pieces);
                    let from = if i == fs { fk as usize } else { 0 };
                    let to = if i == ms { mk as usize } else { pieces.len() };
                    if from < to {
                        s.free.extend_from_slice(&pieces[from..to]);
                    }
                }
                free_n.store(s.avail(slots_n), Ordering::Relaxed);
                given_n.store(s.free.len(), Ordering::Relaxed);
            }
        }
        min_next
    };

    // One piece of the lowest parked strip, or of the next strip within `ahead` strips of the
    // slowest partition's (`low`), if a slot of the pool is free for it; `false` if not.
    let produce = |low: usize| -> bool {
        // (the pieces every partition has folded go back first, so that the piece takes one of
        // their slots before one no piece has used)
        if given_n.load(Ordering::Relaxed) == 0 {
            slowest();
        }
        if free_n.load(Ordering::Relaxed) <= RESERVED {
            slowest();
            // (the strip it would take: the lowest parked one, else the next one)
            let next = if parked_n.load(Ordering::Relaxed) > 0 {
                low_parked.load(Ordering::Relaxed)
            } else {
                cursor_hint.load(Ordering::Relaxed)
            };
            if free_n.load(Ordering::Relaxed) <= RESERVED
                && (next != low || free_n.load(Ordering::Relaxed) == 0)
            {
                return false;
            }
        }
        let (strip, from, k, slot) = {
            let mut s = lock(&sched);
            let lowest = s
                .parked
                .iter()
                .enumerate()
                .min_by_key(|(_, p)| p.strip)
                .map(|(at, p)| (at, *p));
            let (at, strip) = match lowest {
                Some((at, p)) => (Some(at), p.strip),
                None if s.cursor < n && s.cursor < low + ahead => (None, s.cursor),
                None => return false,
            };
            let avail = s.avail(slots_n);
            let may = avail > RESERVED || (strip == low && avail > 0);
            if !may {
                return false;
            }
            let slot = s.take();
            free_n.store(s.avail(slots_n), Ordering::Relaxed);
            given_n.store(s.free.len(), Ordering::Relaxed);
            let (from, k) = match at {
                Some(at) => {
                    let p = s.parked.swap_remove(at);
                    parked_n.store(s.parked.len(), Ordering::Relaxed);
                    low_parked.store(lowest_parked(&s.parked), Ordering::Relaxed);
                    (Some(p.from), p.k)
                }
                None => {
                    s.cursor += 1;
                    cursor_hint.store(s.cursor, Ordering::Relaxed);
                    (None, 0)
                }
            };
            (strip, from, k, slot)
        };
        let (next, k_end, empty) = {
            let mut guard = write(&pool[slot as usize]);
            let arena = &mut *guard;
            arena.reset(parts, piece_seeds);
            let mut sink = PieceSink {
                table: &table,
                seeds: arena.seeds,
                pos: arena.pos,
                links: &mut arena.links,
                k,
                limit: k + piece_seeds as u32,
            };
            let next = scan(strip, from, &mut sink);
            let k_end = sink.k;
            (next, k_end, arena.links.used == 0)
        };
        let st = &strips[strip];
        let count = if empty {
            // (a piece stops only when it is full, so an empty one ends its strip)
            debug_assert!(next.is_none());
            let mut s = lock(&sched);
            s.free.push(slot);
            free_n.store(s.avail(slots_n), Ordering::Relaxed);
            given_n.store(s.free.len(), Ordering::Relaxed);
            st.published.load(Ordering::Relaxed)
        } else {
            let mut pieces = lock(&st.pieces);
            pieces.push(slot);
            pieces.len() as u32
        };
        st.published.store(
            count | if next.is_none() { DONE } else { 0 },
            Ordering::Release,
        );
        produced.fetch_add(1, Ordering::Relaxed);
        if let Some(from) = next {
            let mut s = lock(&sched);
            s.parked.push(Parked {
                strip,
                from,
                k: k_end,
            });
            parked_n.store(s.parked.len(), Ordering::Relaxed);
            low_parked.store(lowest_parked(&s.parked), Ordering::Relaxed);
        }
        true
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
        let (mut scan_ns, mut wait_ns, mut pause_ns) = (0u64, 0u64, 0u64);
        let (mut my_runs, mut my_away) = (0u64, 0u64);
        let mut wait_start: Option<std::time::Instant> = None;
        let mut idle = 0u32;
        let mut stall: Option<(std::time::Instant, u64)> = None;
        // the slowest partition's next piece as this participant last read it (it only grows)
        let mut seen = 0u64;

        // Folds partition p's ready pieces from its next one on, in order, if it can lock p, and
        // returns how many steps (a piece folded, or a strip whose pieces it has all folded left):
        // every ready piece (`Run::Home`, `Run::Tail`), while p has no participant at home
        // (`Run::Orphan`), or one step (`Run::Slowest`).
        let mut run = |p: usize, kind: Run, pause_ns: &mut u64| -> u64 {
            let lane = &lanes[p].0;
            let ready = |v: u64| {
                let i = strip_of(v);
                i < n && {
                    let published = strips[i].published.load(Ordering::Acquire);
                    v as u32 != published & !DONE || published & DONE != 0
                }
            };
            if !ready(lane.load(Ordering::Acquire)) {
                return 0;
            }
            let Ok(mut guard) = states[p].try_lock() else {
                return 0;
            };
            let st = &mut guard.0;
            // (only the participant that holds p's lock moves p's lane)
            let mut v = lane.load(Ordering::Relaxed);
            // (the steps: a piece folded or a strip left; `away` counts the pieces only, so that
            // it stays within `units`, the pieces made times the partitions)
            let (mut folded, mut pieces) = (0u64, 0u64);
            while strip_of(v) < n {
                let (i, k) = (strip_of(v), v as u32);
                let published = strips[i].published.load(Ordering::Acquire);
                if k < published & !DONE {
                    let slot = lock(&strips[i].pieces)[k as usize];
                    let arena = read(&pool[slot as usize]);
                    st.out.strip = i as u64;
                    let t0 = timed.then(std::time::Instant::now);
                    let pause0 = st.scan_pause_ns;
                    arena.for_each_run(p, |seeds, pos| fold(st, seeds, pos));
                    if let Some(t0) = t0 {
                        st.fold_ns += t0.elapsed().as_nanos() as u64;
                    }
                    if owner {
                        *pause_ns += st.scan_pause_ns - pause0;
                    }
                    drop(arena);
                    v = pack(i, k + 1);
                    pieces += 1;
                } else if published & DONE != 0 {
                    v = pack(i + 1, 0);
                } else {
                    break;
                }
                lane.store(v, Ordering::Release);
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
                    my_away += pieces;
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
            let work = |seen: u64| {
                parked_n.load(Ordering::Relaxed) > 0
                    || cursor_hint.load(Ordering::Relaxed) < n.min(strip_of(seen) + ahead)
            };
            if !did || !work(seen) {
                seen = slowest();
                if strip_of(seen) >= n {
                    failing.done = true;
                    break;
                }
            }
            if work(seen) {
                let t0 = (timed && owner).then(std::time::Instant::now);
                if produce(strip_of(seen)) {
                    did = true;
                }
                if let Some(t0) = t0 {
                    scan_ns += t0.elapsed().as_nanos() as u64;
                }
            }
            if !did {
                // nothing of its own: a piece of the slowest partition, or once every strip is
                // handed out and none is parked every ready piece of any partition
                let tail = cursor_hint.load(Ordering::Relaxed) >= n
                    && parked_n.load(Ordering::Relaxed) == 0;
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
                    // the owner can finish the stage alone; the stage's progress: pieces made
                    // plus the partitions' next pieces
                    let now = lanes
                        .iter()
                        .map(|lane| lane.0.load(Ordering::Relaxed))
                        .fold(produced.load(Ordering::Relaxed), u64::wrapping_add);
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
        units: produced.into_inner() * parts as u64,
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
    // (the participant that saw every partition past the last strip gave every piece back)
    debug_assert_eq!(freed.load(Ordering::Relaxed), pack(n, 0));
    debug_assert_eq!(
        sched
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .avail(slots_n),
        slots_n
    );
    scratch.parts = states;
    scratch.arenas = pool
        .into_iter()
        .map(|a| a.into_inner().unwrap_or_else(PoisonError::into_inner).links)
        .collect();
    scratch.store = store;
    stats
}

#[cfg(test)]
mod tests {
    use super::{Arena, Links, BLOCK, NONE};

    /// `Links::take` hands out distinct blocks of the arena for any cut of a piece's seeds by
    /// partition (the pushes of `PieceSink::push`: a block when the partition's last one is full),
    /// up to `piece_seeds` seeds, and a partition within its share keeps its seeds in its lane.
    #[test]
    fn take_gives_each_block_once_and_lanes_to_shares() {
        let mut rng = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = |n: usize| {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            (rng % n as u64) as usize
        };
        for parts in [1usize, 2, 3, 16, 64] {
            for piece_seeds in [1usize, 2, 63, 64, 65, 300, 4000, 16_384] {
                let blocks = Arena::blocks(parts, piece_seeds);
                for cut in 0..4 {
                    let mut links = Links::default();
                    links.reset(parts, piece_seeds, blocks);
                    let mut fill = vec![BLOCK; parts];
                    let mut taken = vec![false; blocks];
                    let mut mine: Vec<Vec<u16>> = vec![Vec::new(); parts];
                    for k in 0..piece_seeds {
                        // all to one partition, round robin, uniform, or mostly one
                        let p = match cut {
                            0 => 0,
                            1 => k % parts,
                            2 => next(parts),
                            _ if next(4) == 0 => next(parts),
                            _ => parts - 1,
                        };
                        if fill[p] == BLOCK {
                            let b = links.take(p);
                            assert!(
                                (b as usize) < blocks && b != NONE,
                                "{parts} {piece_seeds} {cut}"
                            );
                            assert!(
                                !taken[b as usize],
                                "{parts} {piece_seeds} {cut}: block {b} twice"
                            );
                            taken[b as usize] = true;
                            mine[p].push(b);
                            fill[p] = 0;
                        }
                        fill[p] += 1;
                    }
                    assert_eq!(
                        usize::from(links.used),
                        taken.iter().filter(|&&t| t).count()
                    );
                    // round robin: every partition's share in its own lane, one run of blocks
                    let lane = piece_seeds / BLOCK / parts;
                    if cut == 1 && piece_seeds % (BLOCK * parts) == 0 {
                        for (p, bs) in mine.iter().enumerate() {
                            let want: Vec<u16> =
                                (p * lane..(p + 1) * lane).map(|b| b as u16).collect();
                            assert_eq!(bs, &want, "{parts} {piece_seeds}: partition {p}");
                        }
                    }
                }
            }
        }
    }
}
