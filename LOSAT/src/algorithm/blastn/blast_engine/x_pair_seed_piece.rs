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
//! on from `scan_range[0]`. A piece is an arena of blocks of `BLOCK` seeds (`Arena::blocks`),
//! enough for `PIECE_SEEDS` seeds (the switch's number) whatever partitions they go to; each
//! partition's seeds of the piece are a chain of blocks in scan order: its first block, then runs
//! of a page of seeds taken in order from the front of the arena (`Links`), so that a piece
//! touches the pages of its seeds and a partition's seeds after its first block lie in pages of
//! seeds of its own. The arenas are a pool of `PIECES_PER_THREAD * threads + RESERVED` slots,
//! allocated at their size once and kept in the subject scratch, so a stage holds at most that
//! many pieces whatever the subject, the query or the density of the seeds (and arenas only for
//! as many pieces as it had at once: a piece takes
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
    bucket_of, lock, merge, pause, run_participants, AtHome, Failing, InStage, OwnLines,
    PartScratch, PartState, PartStats, Run, HELPER_STALL, IN_STAGE,
};
use std::cell::{Cell, UnsafeCell};
use std::mem::{size_of, MaybeUninit};
use std::sync::atomic::{AtomicBool, AtomicU32, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Mutex, OnceLock, PoisonError};

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
#[cfg_attr(test, derive(Debug, PartialEq))]
#[repr(C)]
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

/// A partition's run of blocks in an arena: its next free block and its end (`Links::take`).
#[derive(Clone, Copy)]
#[cfg_attr(test, derive(Debug, PartialEq))]
#[repr(C)]
struct Span {
    next: u16,
    end: u16,
}

/// The blocks of a run of a partition's blocks after its first: a page of seeds (4 KiB), so that a
/// partition's seeds after its first block lie in pages of seeds of its own, which the hardware
/// prefetches of a fold's reads do not leave (their positions, 4 bytes each, share a page with
/// those of three other runs).
const RUN: usize = {
    let run = PAGE / (BLOCK * size_of::<Seed>());
    assert!(run >= 1);
    run
};

/// No NCBI counterpart (scheduling only): the chains of a slot's piece: its blocks of `BLOCK`
/// seeds (`used` of them taken), each partition's blocks linked in its order (`next`, its first
/// and last block in `chains`). Partition p's first block is block p of the arena; its later
/// blocks come in runs of `RUN` blocks (a page of seeds), each the next run of the blocks after
/// the first ones (`front`), and `spans` holds each partition's run; once the runs are given out,
/// a partition takes the last free block of another partition's run (`take_tail`: a piece of
/// many seeds over many partitions, near the arena's end). So a piece touches the first blocks of
/// its partitions and the pages of its seeds from the front of the arena, a slot keeps the pages
/// of its largest piece, and a partition's seeds lie in its first block and in pages of seeds of
/// its own. (Blocks taken in the order of the pushes, the partitions' blocks
/// interleaved, made the fold read each partition's seeds in blocks scattered over the arena:
/// the hardware prefetches of the neighbouring blocks of other partitions and the pages of the
/// scattered blocks took the fold +6 % cycles at the same instructions, L2 misses +22 %;
/// dc-megablast cDNA x genome at 4 threads, +2.2 % over the strips stage, line EB s15d / s15e.
/// Lanes of the arena of equal length, one per partition, touched a page of every lane a
/// partition had a seed in, whatever its seeds, and over the pieces of a slot the pages of every
/// lane: megablast Sakai x MG1655 -db at 16 threads, the peak 1.3 MiB over the strips stage's
/// with pieces of 16,384 seeds, 0.5 MiB with 8,192 and 1.2 MiB below it with 4,096, line EB
/// s19m2.)
/// Kept in the subject scratch from one stage to the next; the seeds are in the pool's `Store`.
/// Each of its buffers is on cache lines of its own (`Lined`): a participant writes the chains,
/// runs and links of the slot it fills while others read those of the slots they fold, and as
/// small allocations of their own they shared lines with other slots' (BLASTN-short of 30,000
/// short queries x 3 Mb at 16 threads: the scan's loads of a chain's fill and tail were 84 % of
/// its loads served from another core's cache, which were x7 those of the strips stage, and the
/// scan took x1.6 the memory-ordering machine clears, at the same instructions; +1.2 % cycles,
/// line EB s18x).
#[derive(Default)]
pub(super) struct Links {
    next: Lined<u16>,
    chains: Lined<Chain>,
    spans: Lined<Span>,
    used: u16,
    front: u16,
}

/// A cache line of a `Lined` buffer.
#[derive(Clone, Copy)]
#[repr(C, align(64))]
struct Line([u8; 64]);

/// No NCBI counterpart (memory placement only): a buffer of `len` values of `T` on cache lines of
/// its own: its allocation begins on a line and ends on one, so that no other allocation shares
/// a line with it.
struct Lined<T> {
    lines: Vec<Line>,
    len: usize,
    _of: std::marker::PhantomData<T>,
}

/// The values a `Lined` buffer holds: plain data of 2-byte fields, any bytes a value.
///
/// # Safety
/// `T` has no padding, an alignment of at most 64 and every bit pattern of its size is a value.
unsafe trait Plain: Copy {}
// SAFETY: a `u16`, and `repr(C)` structs of `u16` fields (no padding, alignment 2).
unsafe impl Plain for u16 {}
// SAFETY: as above.
unsafe impl Plain for Chain {}
// SAFETY: as above.
unsafe impl Plain for Span {}

impl<T> Default for Lined<T> {
    fn default() -> Self {
        Lined {
            lines: Vec::new(),
            len: 0,
            _of: std::marker::PhantomData,
        }
    }
}

impl<T: Plain> Lined<T> {
    /// `len` values (those it held where its lines are kept, zero bytes in new lines; the caller
    /// sets the ones it reads); its lines allocated anew only when there are too few.
    fn set_len(&mut self, len: usize) {
        let lines = (len * std::mem::size_of::<T>()).div_ceil(std::mem::size_of::<Line>());
        if self.lines.len() < lines {
            self.lines = vec![Line([0; 64]); lines];
        }
        self.len = len;
    }
}

impl<T: Plain> std::ops::Deref for Lined<T> {
    type Target = [T];

    #[inline(always)]
    fn deref(&self) -> &[T] {
        // SAFETY: `len` values of `T` fit in the lines (`set_len`), which are initialized bytes,
        // aligned for `T`, and any bytes are a value of `T` (`Plain`).
        unsafe { std::slice::from_raw_parts(self.lines.as_ptr().cast::<T>(), self.len) }
    }
}

impl<T: Plain> std::ops::DerefMut for Lined<T> {
    #[inline(always)]
    fn deref_mut(&mut self) -> &mut [T] {
        // SAFETY: as above, and `&mut self` is the only reference to the lines.
        unsafe { std::slice::from_raw_parts_mut(self.lines.as_mut_ptr().cast::<T>(), self.len) }
    }
}

impl Links {
    /// A free block for partition p: the next one of its run, else `take_run`.
    #[inline(always)]
    fn take(&mut self, p: usize) -> u16 {
        debug_assert!(p < self.spans.len());
        // SAFETY: `p` is below `parts`, and there are `parts` runs (`Arena::reset`,
        // `PieceSink::push`).
        let span = unsafe { self.spans.get_unchecked_mut(p) };
        let b = if span.next < span.end {
            span.next += 1;
            span.next - 1
        } else {
            self.take_run(p)
        };
        self.used += 1;
        b
    }

    /// Empty, for `parts` partitions and arenas of `blocks` blocks (`Arena::blocks`; `next`
    /// allocated at its size the first time): partition p's first block is block p, and the
    /// runs begin at the first whole run after the first blocks.
    fn reset(&mut self, parts: usize, blocks: usize) {
        if self.next.len() != blocks {
            self.next.set_len(blocks);
            self.next.fill(NONE);
        }
        self.chains.set_len(parts);
        self.chains.fill(EMPTY);
        self.spans.set_len(parts);
        for (p, at) in self.spans.iter_mut().enumerate() {
            *at = Span {
                next: p as u16,
                end: p as u16 + 1,
            };
        }
        self.used = 0;
        self.front = parts.next_multiple_of(RUN) as u16;
    }

    /// `take` once partition p's run is used up: the next run of the arena becomes p's, and p
    /// takes its first block; once the runs are given out, `take_tail`.
    #[cold]
    #[inline(never)]
    fn take_run(&mut self, p: usize) -> u16 {
        let first = self.front;
        if usize::from(first) + RUN > self.next.len() {
            return self.take_tail(p);
        }
        self.front += RUN as u16;
        self.spans[p] = Span {
            next: first + 1,
            end: first + RUN as u16,
        };
        first
    }

    /// `take` once the arena's runs are given out: the last free block of another partition's run
    /// (an unused first block included), which p takes alone.
    #[cold]
    #[inline(never)]
    fn take_tail(&mut self, p: usize) -> u16 {
        // (a piece takes at most `piece_seeds` seeds, `Cut`, and a partition's last block is the
        // only one not full, so the arena's first blocks and runs hold them, `Arena::blocks`)
        let q = self
            .spans
            .iter()
            .position(|span| span.next < span.end)
            .expect("LOSAT_X_SEEDPIECE: a full arena");
        self.spans[q].end -= 1;
        let b = self.spans[q].end;
        self.spans[p] = Span {
            next: b + 1,
            end: b + 1,
        };
        b
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
/// order of every slot of the pool, slot `s` at `s * len..(s + 1) * len` from the first page of
/// each buffer (`len`: an arena's blocks of `BLOCK` seeds), kept in the subject scratch from one
/// stage to the next.
///
/// Each buffer is one allocation whose capacity is above glibc's ceiling of the mmap threshold
/// (`utils/x_mapped.rs` `x_unthresholded_capacity`), not initialized: glibc maps it, its free
/// leaves the threshold as it is, and only the pages of the blocks pieces wrote are resident.
/// (A store kept beyond its scratch, which lives for one query batch, held the union of the
/// pages every batch wrote (a process-wide cache of dropped stores, 7373fb967): BLASTN cDNA x
/// genome -subject at 16 threads, peak x1.13 of the strips stage against x0.90 for one store per
/// batch, line EB s15j.) (An arena allocated by itself,
/// 266-328 KiB of seeds, was mapped and its free raised glibc's threshold to its size, so that
/// later blocks of up to that size came from the arenas of the heap and stayed resident there:
/// BLASTN cDNA x genome -subject at 16 threads, peak x1.13-1.15 of the strips stage, x1.04 with
/// the threshold held at 128 KiB, line EB s14u / s14w.)
#[derive(Default)]
pub(super) struct Store {
    seeds: Vec<MaybeUninit<Seed>>,
    pos: Vec<MaybeUninit<u32>>,
    /// the seeds of a slot's region in the layout of the buffers (0 before the first stage)
    len: usize,
    /// where the regions begin in each buffer: its first page (`page_start`)
    at: (usize, usize),
}

impl Store {
    /// The regions of `slots` slots of `len` seeds each (the buffers allocated anew when they
    /// are too small).
    #[allow(clippy::type_complexity)]
    fn slots(
        &mut self,
        slots: usize,
        len: usize,
    ) -> impl Iterator<Item = (&mut [MaybeUninit<Seed>], &mut [MaybeUninit<u32>])> {
        let need = slots * len;
        // (a page more, so that the seeds' regions can begin on a page: glibc's mapped blocks
        // begin after a header, and a run of a partition's seeds is a page, `RUN`)
        let (seeds_n, pos_n) = (
            need + PAGE / size_of::<Seed>(),
            need + PAGE / size_of::<u32>(),
        );
        // (a stage of another layout allocates anew, so that the pages the regions of the old one
        // touched are not kept; one of the same layout with more slots takes the capacity the
        // buffers have, above glibc's ceiling)
        if self.len != len || self.seeds.capacity() < seeds_n || self.pos.capacity() < pos_n {
            // (the old buffers go first, so that the two are not held at once)
            self.seeds = Vec::new();
            self.pos = Vec::new();
            self.seeds = uninit(seeds_n);
            self.pos = uninit(pos_n);
            self.len = len;
            self.at = (page_start(&self.seeds), page_start(&self.pos));
        } else if self.seeds.len() < seeds_n || self.pos.len() < pos_n {
            // SAFETY: `MaybeUninit` needs no initialization, and the capacities are at least
            // `seeds_n` and `pos_n` (the test above).
            unsafe {
                self.seeds.set_len(seeds_n);
                self.pos.set_len(pos_n);
            }
        }
        let (a, b) = self.at;
        self.seeds[a..a + need]
            .chunks_mut(len)
            .zip(self.pos[b..b + need].chunks_mut(len))
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
    /// The blocks of an arena for `parts` partitions: their first blocks (`parts`, up to a whole
    /// run) and runs of `piece_seeds.div_ceil(BLOCK)` blocks, which hold `piece_seeds` seeds in
    /// any cut, since a partition's last block is the only one not full and every free block can
    /// be taken (`Links::take_tail`). A whole number of runs, so that every slot's seeds begin on
    /// a run.
    fn blocks(parts: usize, piece_seeds: usize) -> usize {
        parts.next_multiple_of(RUN) + piece_seeds.div_ceil(BLOCK).next_multiple_of(RUN)
    }

    /// Empty, for `parts` partitions and pieces of `piece_seeds` seeds (its region of the
    /// `Store` is `blocks * BLOCK` seeds; its links allocated at their size the first time).
    fn reset(&mut self, parts: usize, piece_seeds: usize) {
        let blocks = Self::blocks(parts, piece_seeds);
        assert!(
            self.seeds.len() == blocks * BLOCK && self.pos.len() == blocks * BLOCK,
            "LOSAT_X_SEEDPIECE: an arena's region of the store"
        );
        self.links.reset(parts, blocks);
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

/// A page of memory, the unit the regions of the `Store` and the runs of an arena keep to.
const PAGE: usize = 4096;

/// The first slot of `v` on a page (the first slot if its address cannot be on one).
fn page_start<T>(v: &[T]) -> usize {
    let size = size_of::<T>().max(1);
    let skip = (PAGE - v.as_ptr() as usize % PAGE) % PAGE;
    if skip % size == 0 {
        skip / size
    } else {
        0
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
        // SAFETY: the table holds `bucket % parts`, and there are `parts` chains and `parts` runs
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

/// No NCBI counterpart (scheduling only): how many pieces a strip has, with `DONE` once its scan
/// is done (their slots are in `PieceSlots`).
#[derive(Default)]
struct StripState {
    published: AtomicU32,
}

/// No NCBI counterpart (scheduling only): the slot in the pool of each piece a partition may still
/// fold, which the participants read without a lock: piece k of strip i at `(i % strips) * row
/// + k % pieces`, for `strips` the strips handed out past the slowest partition's (`ahead`)
/// and `pieces` the slots of the pool. Piece k's entry is written by the piece's producer before
/// the strip's count says it is there (`Release`) and read after the count is read (`Acquire`),
/// and it is not written again while a partition has piece k to fold or its slot is to go back:
/// that piece holds its slot until every partition has folded it (`freed`), and so do the later
/// pieces, so a strip's piece k + `pieces` would need one more slot than the pool has; and strip
/// i + `strips` is handed out only once every partition has left strip i and its pieces have gone
/// back (under the schedule's lock, which the producer takes first).
/// (A lock of the strip's list of slots and a read lock of the slot, taken by every partition for
/// every piece, wrote the strip's and the slot's lines that the other partitions read: BLASTN-short
/// of 30,000 short queries x 3 Mb at 16 threads, the participants took x4 the loads served from
/// another core's cache of the strips stage's, 73 % of them at the slot's lock and the strip's
/// count, line EB s19y.)
struct PieceSlots {
    slots: Vec<AtomicU32>,
    strips: usize,
    pieces: usize,
    /// the entries from one strip's row to the next: its pieces and a cache line, so that no two
    /// rows share a line (a producer writes its strip's row while other strips' folds read theirs)
    row: usize,
}

impl PieceSlots {
    fn new(strips: usize, pieces: usize) -> Self {
        let (strips, pieces) = (strips.max(1), pieces.max(1));
        let row = pieces.next_multiple_of(16) + 16;
        PieceSlots {
            slots: (0..strips * row).map(|_| AtomicU32::new(0)).collect(),
            strips,
            pieces,
            row,
        }
    }

    #[inline(always)]
    fn at(&self, strip: usize, k: u32) -> &AtomicU32 {
        &self.slots[(strip % self.strips) * self.row + k as usize % self.pieces]
    }
}

/// No NCBI counterpart (scheduling only): a slot of the pool. Its arena is written by the one
/// participant that took the slot from the schedule, and read by the folds once its piece is
/// published, until every partition has folded it (`Slot::fill`, `Slot::folded`); no lock.
struct Slot<'s>(UnsafeCell<Arena<'s>>);

// SAFETY: the arena is shared between threads only as `Slot::fill` and `Slot::folded` say: one
// writer, or readers after the writer's release; `Arena` itself holds plain data and slices.
unsafe impl Sync for Slot<'_> {}

impl<'s> Slot<'s> {
    /// The arena, to fill a piece in it.
    ///
    /// # Safety
    /// The caller took the slot from the schedule (`Sched::take`), so no other participant
    /// fills it, and every partition folded the piece it held before, its lane stored after the
    /// fold (`Release`) and read before the slot went back to the schedule (`Acquire`).
    #[allow(clippy::mut_from_ref)]
    #[inline(always)]
    unsafe fn fill(&self) -> &mut Arena<'s> {
        // SAFETY: as the caller says, nothing else refers to the arena.
        unsafe { &mut *self.0.get() }
    }

    /// The arena of a published piece, to fold it.
    ///
    /// # Safety
    /// The caller read the strip's count with the piece in it (`Acquire`, after the producer's
    /// `Release`) and has not stored its partition's lane past the piece, so the slot has not
    /// gone back to the schedule.
    #[inline(always)]
    unsafe fn folded(&self) -> &Arena<'s> {
        // SAFETY: as the caller says, no one writes the arena while the reference lives.
        unsafe { &*self.0.get() }
    }
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
    // (each slot's arena on lines of its own: a participant writes the slot it fills while others
    // read their neighbours to fold them)
    let pool: Vec<OwnLines<Slot<'_>>> = store
        .slots(slots_n, Arena::blocks(parts, piece_seeds) * BLOCK)
        .zip(links)
        .map(|((seeds, pos), links)| OwnLines(Slot(UnsafeCell::new(Arena { seeds, pos, links }))))
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
    // (each strip's state on lines of its own: the participant that makes a piece of a strip
    // writes its count while others read those of the strips they fold)
    let strips: Vec<OwnLines<StripState>> =
        (0..n).map(|_| OwnLines(StripState::default())).collect();
    let piece_slots = PieceSlots::new(ahead.min(n), slots_n);
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
        // (a load first, as in `seed_stage_home`)
        if freed.load(Ordering::Relaxed) < min_next {
            let mut s = lock(&sched);
            let f = freed.load(Ordering::Relaxed);
            if f < min_next {
                freed.store(min_next, Ordering::Relaxed);
                let (fs, fk) = (strip_of(f), f as u32);
                let (ms, mk) = (strip_of(min_next), min_next as u32);
                for (i, strip) in strips.iter().enumerate().take((ms + 1).min(n)).skip(fs) {
                    let from = if i == fs { fk } else { 0 };
                    let to = if i == ms {
                        mk
                    } else {
                        strip.0.published.load(Ordering::Acquire) & !DONE
                    };
                    // (the pieces below `min_next`: every partition folded them, after their
                    // producers wrote their slots)
                    s.free
                        .extend((from..to).map(|k| piece_slots.at(i, k).load(Ordering::Relaxed)));
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
            // SAFETY: the slot was taken from the schedule just now, and a slot goes back to it
            // only when every partition's lane is past its piece (`slowest`).
            let arena = unsafe { pool[slot as usize].0.fill() };
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
        let st = &strips[strip].0;
        let count = if empty {
            // (a piece stops only when it is full, so an empty one ends its strip)
            debug_assert!(next.is_none());
            let mut s = lock(&sched);
            s.free.push(slot);
            free_n.store(s.avail(slots_n), Ordering::Relaxed);
            given_n.store(s.free.len(), Ordering::Relaxed);
            st.published.load(Ordering::Relaxed)
        } else {
            // (one participant at a time makes the strip's pieces: the parked strip passes
            // through the schedule's lock, so the count is the one it last stored)
            let k = st.published.load(Ordering::Relaxed) & !DONE;
            piece_slots.at(strip, k).store(slot, Ordering::Relaxed);
            k + 1
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
                    let published = strips[i].0.published.load(Ordering::Acquire);
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
            // (the strip's count as last read: read again only when it says no more)
            let mut known = (usize::MAX, 0u32);
            while strip_of(v) < n {
                let (i, k) = (strip_of(v), v as u32);
                let published = if known.0 == i && k < known.1 & !DONE {
                    known.1
                } else {
                    let published = strips[i].0.published.load(Ordering::Acquire);
                    known = (i, published);
                    published
                };
                if k < published & !DONE {
                    let slot = piece_slots.at(i, k).load(Ordering::Relaxed);
                    // SAFETY: the count read above (`Acquire`) has piece k, and p's lane is
                    // stored past it only after the fold below.
                    let arena = unsafe { pool[slot as usize].0.folded() };
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
        .map(|a| a.0 .0.into_inner().links)
        .collect();
    scratch.store = store;
    stats
}

#[cfg(test)]
mod tests {
    use super::{Arena, Lined, Links, Span, BLOCK, EMPTY, NONE, RUN};

    /// A `Lined` buffer begins on a cache line and holds whole lines, for any length, keeps its
    /// lines when it shrinks and grows within them, and holds the values set.
    #[test]
    fn lined_buffers_are_lines_of_their_own() {
        fn check<T: super::Plain + PartialEq + std::fmt::Debug>(v: T, lens: &[usize]) {
            let mut buf: Lined<T> = Lined::default();
            assert!(buf.is_empty());
            for &len in lens {
                let before = buf.lines.as_ptr();
                let had = buf.lines.len();
                buf.set_len(len);
                buf.fill(v);
                assert_eq!(buf.len(), len);
                assert_eq!(buf.as_ptr() as usize % 64, 0, "{len}");
                assert!(
                    buf.lines.len() * 64 >= len * std::mem::size_of::<T>(),
                    "{len}"
                );
                if had * 64 >= len * std::mem::size_of::<T>() {
                    assert_eq!(buf.lines.as_ptr(), before, "{len}: lines kept");
                }
                assert!(buf.iter().all(|x| *x == v), "{len}");
            }
        }
        check(7u16, &[1, 31, 32, 33, 272, 16, 300]);
        check(EMPTY, &[1, 10, 11, 16, 64, 3]);
        check(Span { next: 3, end: 9 }, &[1, 16, 17, 65, 2]);
    }

    /// `Links::take` hands out distinct blocks of the arena for any cut of a piece's seeds by
    /// partition (the pushes of `PieceSink::push`: a block when the partition's last one is full),
    /// up to `piece_seeds` seeds: partition p's first block is block p, its later ones whole runs
    /// of `RUN` blocks that begin on a run, and the runs taken are the first ones after the first
    /// blocks (a piece touches the front of its arena); once the runs are given out, single blocks
    /// of other partitions' runs.
    #[test]
    fn take_gives_each_block_once_in_runs_from_the_front() {
        let mut rng = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = |n: usize| {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            (rng % n as u64) as usize
        };
        let mut tails = 0;
        for parts in [1usize, 2, 3, 5, 16, 64, 128] {
            for piece_seeds in [1usize, 2, 63, 64, 65, 300, 4000, 4096, 16_384] {
                let blocks = Arena::blocks(parts, piece_seeds);
                assert_eq!(blocks % RUN, 0, "{parts} {piece_seeds}");
                let front = parts.next_multiple_of(RUN);
                for cut in 0..6 {
                    let mut links = Links::default();
                    links.reset(parts, blocks);
                    let mut fill = vec![BLOCK; parts];
                    let mut taken = vec![false; blocks];
                    let mut mine: Vec<Vec<u16>> = vec![Vec::new(); parts];
                    for k in 0..piece_seeds {
                        // all to one partition, round robin, uniform, mostly one, each
                        // partition's seeds one block and one more seed, or that and then a
                        // run and one more seed in turn (the most blocks: every run begun)
                        let p = match cut {
                            0 => 0,
                            1 => k % parts,
                            2 => next(parts),
                            3 if next(4) == 0 => next(parts),
                            3 => parts - 1,
                            4 => (k / (BLOCK + 1)) % parts,
                            _ if k < (BLOCK + 1) * parts => k / (BLOCK + 1),
                            _ => ((k - (BLOCK + 1) * parts) / (BLOCK * RUN)) % parts,
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
                    // (blocks of other runs only once the runs are given out)
                    if usize::from(links.front) == blocks {
                        tails += 1;
                        continue;
                    }
                    let mut runs = 0;
                    for (p, bs) in mine.iter().enumerate() {
                        if let Some(&first) = bs.first() {
                            assert_eq!(usize::from(first), p, "{parts} {piece_seeds} {cut}");
                        }
                        for (i, &b) in bs.iter().enumerate().skip(1) {
                            let b = usize::from(b);
                            if (i - 1) % RUN == 0 {
                                assert!(b >= front && b % RUN == 0, "{parts} {piece_seeds} {cut}");
                                runs += 1;
                            } else {
                                assert_eq!(
                                    b,
                                    usize::from(bs[i - 1]) + 1,
                                    "{parts} {piece_seeds} {cut}"
                                );
                            }
                        }
                    }
                    // the runs taken are the first ones after the first blocks
                    assert!(
                        taken[front + runs * RUN..].iter().all(|&t| !t),
                        "{parts} {piece_seeds} {cut}"
                    );
                    assert_eq!(
                        usize::from(links.front),
                        front + runs * RUN,
                        "{parts} {piece_seeds} {cut}"
                    );
                }
            }
        }
        // (the cuts reach the arena's end)
        assert!(tails > 0);
    }
}
