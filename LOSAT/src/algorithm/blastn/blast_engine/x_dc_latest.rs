//! EXPERIMENT (LOSAT_X_DCLATEST, a strict-set switch since 2026-11-23): the two-hit fold of
//! dc-megablast's strips stage (`x_pair_seed_disc.rs` `fold_two_hit_t`: the owner's fold without a
//! pool or partitions, and each partition's fold of `LOSAT_X_PAIRPARFOLD`) keeps per diagonal only
//! its latest insert, in a map of its own (`DiagLatest`), instead of NCBI's diagonal hash, 512
//! buckets of chains of cells whose stale cells later inserts of the bucket take over. No NCBI
//! counterpart (the layout of the table only).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:361-381
//! (`s_BlastDiagHashRetrieve`)
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:396-449
//! (`s_BlastDiagHashInsert`)
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:825-858,939-941
//! (`s_BlastnDiagHashExtendInitialHit`: what a seed reads from the table and what it inserts)
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:978-982
//! (`s_BlastNaExtendDirect`: a discontiguous table's word length is its template length)
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:174-181
//! (the hash offset after each chunk, and the clear of the table)
//!
//! Names below: S is a seed's `s_off_pos` (its subject offset plus the hash offset), L an entry's
//! level, w the window (`window_size`), wl the word length (the template length, 16, 18 or 21;
//! the lookup width is the same), and W the take-over window of NCBI's insert (w plus Delta
//! plus 1, Delta the smaller of `scan_range` and w - wl, set to 0 where it is negative by the
//! two-hit branch of a seed whose word is a single word).
//!
//! Why the output is the same:
//! - (i) What a seed reads. The seed's own diagonal through the retrieve: the level as
//!   `last_hit` (0 with no cell) and `hit_saved` (false with no cell), read only by the
//!   `last_hit` test (the seed is rejected when S < L) and the window test (the two-hit branch
//!   runs when `hit_saved` is set or S + wl > L + w). The retrieves of other diagonals happen only
//!   in the off-diagonal loop, which runs no iteration with `scan_range == 0` (the scope of the
//!   strips stage), so `hit_len` is read by nothing and the map does not keep it. The three
//!   inserts take their values from the seed and the branch it took, never from what the
//!   retrieve returned.
//! - (ii) Dead entries. An entry decides as no cell at S when L lies below the level bound of
//!   LOSAT_X_DIAGLV (`x_diag_lv.rs` `LvSeeds::bound`: S + 1, with two hits at most S + wl - w; its
//!   proof (ii) covers both tests). The bound only grows with S, so an entry dead at one seed is
//!   dead at every later seed whose S is not smaller. (`LvSeeds` gives no bound, `i32::MIN`, to
//!   positions outside its range, where a dead entry could decide again; the stage takes the map
//!   only where every position of the chunk has a bound, `x_pair_seed_disc.rs` `latest_bound`,
//!   which every chunk has. `put` and `grow` still count no entry as dead without a bound.)
//! - (iii) The map against the chains. NCBI's insert writes into the first cell of the bucket's
//!   chain that has the diagonal or is stale (S' - L > W at the inserting seed S'), else appends a
//!   cell at the head, and the retrieve returns the first cell of the diagonal. So the first cell
//!   of diagonal d holds d's latest insert unless an insert of another diagonal took that cell
//!   over; then the retrieve returns an older insert of d, or nothing. The levels of d's inserts
//!   never decrease: an insert's level is at least its seed's S, and a seed that is not rejected
//!   saw a level at most S (the latest, or a cell taken over, which is older). A take-over at S'
//!   needs L < S' - W, so the cell is dead at S' when W + 1 >= the bound's distance D (S - L >= D
//!   is dead; D = 0 when w < wl, else w - wl + 1). With `scan_range == 0` NCBI's W is w + 1: Delta
//!   is 0 when w >= wl, and when w < wl every seed that is not rejected takes the two-hit branch
//!   (S + wl > L + w when S >= L) with a single word (`s_TypeOfWord` returns at once when the word
//!   length equals the lookup width), which sets Delta to 0 before the insert. So W + 1 >= D, every
//!   taken-over cell is dead, NCBI's chains return d's latest insert or something dead (an older
//!   level, or no cell), and the map (iv) returns the latest insert or something dead: the same
//!   decisions, so the same inserts, seed after seed. (LOSAT's chains take their window once per
//!   chunk without that step, `diag_hash_window`: 2w - wl + 1 when w < wl; see (v).)
//! - (iv) What `find` returns. Every entry holds the values of an insert of its diagonal. When the
//!   latest insert of d is live, `find` returns it: `put` overwrites only the slot `find` found
//!   for the seed's own diagonal, or a slot whose entry is dead (a free slot, (ii)) or empty; so a
//!   live latest is never overwritten, and `grow` moves the live entries and drops only dead ones.
//!   At most one entry of a diagonal is live, and it is the first of that diagonal in `find`'s
//!   order (the home bucket's slots, then the stash): a diagonal enters a home slot only when
//!   `find` found none of it there and none live in the stash, and enters the stash only where
//!   its home bucket has no free slot, over its first entry there if it has one. A stash entry is
//!   passed over only when the seed's S is at or above its home bucket's `stash_until`, which
//!   every stash write raises past the last S at which the entry can be live. So when d's latest
//!   is dead, `find` returns d's latest, an older entry of d (whose level is no higher, (iii)),
//!   or nothing: dead, all of them.
//! - (v) The band 2w + 2 < wl (w at most 6, 7 or 9 for templates 16, 18 or 21) is where LOSAT's
//!   chains, with their window 2w - wl + 1 (iii), can take over a live cell, so there their
//!   decisions depend on the order of their buckets' cells and can differ from NCBI's and the
//!   map's (a seed rejected by one and not the other). No output depends on them: with w < wl a
//!   seed that is not rejected always takes the two-hit branch with a single word, no
//!   off-diagonal cell is read (Delta is 0), `hit_ready` is cleared, and the seed only inserts.
//!   dc-megablast with w < wl extends no seed in NCBI, in the chains or in the map; only the count
//!   of rejected seeds (`[TIMING] x_pairpardisc_seed: skipped=`) can differ there.
//! - (vi) The order. Within a stage the seeds reach each map with S never decreasing: the strips
//!   in order (the owner's `run_strips`; each partition its strips, or its pieces of a strip, in
//!   order with one folder at a time), a strip's words in subject order, a word's two templates
//!   at the same S, a masked subject's ranges in the order they are scanned (`plan_strips`), and
//!   a partition gets a subsequence of the stage's seeds. The ranges' order is the one input
//!   where this fails: a subject of `-db` with `-db_soft_mask` keeps its stored masks' order
//!   (`subject_split.rs`), so its ranges can step back (S goes back to a range's start after a
//!   later range), and NCBI scans them in that order; which cells its chains still hold then
//!   depends on the order of their take-overs, so no per-diagonal map can follow them. The stage
//!   sends such a chunk to the chains (`x_pair_seed_disc.rs` `latest_ranges_ascend`), on every
//!   path that folds on the maps (the owner's fold, the partitions', one thread or many), as it
//!   does a chunk without the bounds of (ii). `put` asserts the order in debug builds.
//! - (vii) Stages, chunks and the offset clear. A map lives from stage to stage in the subject
//!   scratch (`x_pair_seed_part.rs` `PartScratch`, `PartState`). A stage's levels are at most its
//!   offset plus its chunk's length, and NCBI moves the offset past that plus w after every chunk,
//!   so at the next stage every entry lies below the bound of its first position: dead for the
//!   whole stage, like every cell NCBI's chains carry over from earlier chunks (also those that
//!   the reference fold inserted between the stages, which the map never sees; and the map's
//!   inserts never reach the chains, whose cells of the chunk are dead at the next chunk too).
//!   `begin` checks this at each stage's start with `top`, an upper bound of the levels the map
//!   holds, and empties the map where it fails, which happens only after the offset was cleared
//!   (at a quarter of the largest Int4), where NCBI empties its table: there the old levels lie above the new
//!   offset and would reject seeds.
//!
//! The layout: buckets of one cache line, `Bucket`, with 7 slots each (the diagonals, the levels,
//! and in the eighth lanes the bucket's `stash_until` and its slots' `hit_saved` bits); a home
//! bucket chosen by the high bits of the diagonal times a Fibonacci constant (NCBI's bucket is
//! the diagonal modulo 512, the partition of `LOSAT_X_PAIRPARFOLD` that bucket modulo the
//! partitions, so the diagonals of a partition share their low bits), whose slots `find` and `put`
//! compare at once (`lanes_mask`: two 128-bit compares on x86_64 and on wasm32 with simd128, lane
//! by lane elsewhere); a stash of `STASH` entries
//! for a bucket without a free slot, and the buckets doubled when the stash is full and they are
//! fewer than twice the entries, keeping the live entries (`grow`). The fold brings in the home bucket of the seed `AHEAD` places ahead
//! (`prefetch`). A slot is free when its level is below the seed's bound (dead, (ii)) or the slot
//! is empty (level `i32::MIN`). The size is a capacity of the structure (the chains of NCBI's
//! table double too, `DiagHashTable::insert_inl`), the same code for every input; it chooses no
//! path.

use crate::algorithm::blastn::lookup::x_dc_csr::prefetch;
use std::sync::atomic::{AtomicU64, Ordering};

/// The slots of a bucket (lanes 0 to 6; lane 7 holds `stash_until` and the saved bits).
const SLOTS: u32 = 7;
const SLOT_MASK: u32 = (1 << SLOTS) - 1;
/// Lane 7 of `Bucket::diag`: S below which an entry of the stash homed here may be live.
const UNTIL: usize = 7;
/// Lane 7 of `Bucket::lvl`: bit k set when slot k's entry has `hit_saved`.
const SAVED: usize = 7;
/// An empty slot's diagonal and level. No diagonal is `i32::MIN` (`s_off - q_off` with both
/// offsets in `0..i32::MAX`), and no level is below 0.
const NO_DIAG: i32 = i32::MIN;
const NO_LEVEL: i32 = i32::MIN;
/// The entries of the stash at which the buckets double (`put_cold`), while they are fewer than
/// twice the entries the map holds.
const STASH: usize = 32;
/// The buckets of a new map, as a power of two (64 buckets, 4 KiB).
const MIN_BITS: u32 = 6;
/// How many seeds ahead the fold brings in a seed's home bucket.
pub(super) const AHEAD: usize = 8;

/// `find`'s place of an entry: a slot of the home bucket (0 to 6), none (`MISS`), or a stash
/// entry (`STASHED` plus its index).
const MISS: u32 = 7;
const STASHED: u32 = 8;

// No NCBI counterpart: counters of LOSAT_TIMING (the stages folded on the maps, the chunks sent to
// the chains because their scan ranges step back, the map's growths and its largest size).
static STAGES: AtomicU64 = AtomicU64::new(0);
static BACK: AtomicU64 = AtomicU64::new(0);
static GROWS: AtomicU64 = AtomicU64::new(0);
static MOST_BUCKETS: AtomicU64 = AtomicU64::new(0);

/// One cache line: slots 0 to 6 (`diag`, `lvl`), the bucket's `stash_until` (`diag[UNTIL]`) and
/// its slots' saved bits (`lvl[SAVED]`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(C, align(64))]
struct Bucket {
    diag: [i32; 8],
    lvl: [i32; 8],
}

/// A bucket without entries (and with nothing in the stash: `stash_until` is `i32::MIN`).
const EMPTY: Bucket = Bucket {
    diag: [NO_DIAG; 8],
    lvl: [
        NO_LEVEL, NO_LEVEL, NO_LEVEL, NO_LEVEL, NO_LEVEL, NO_LEVEL, NO_LEVEL, 0,
    ],
};

/// An entry of the stash (`NO_DIAG` and `NO_LEVEL` when empty).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Entry {
    diag: i32,
    level: i32,
    saved: bool,
}

/// What `find` returned for a seed (the retrieve's `last_hit` and `hit_saved`, 0 and false for
/// no entry) and where the seed's `put` writes.
#[derive(Clone, Copy, Debug)]
pub(super) struct Found {
    pub(super) level: i32,
    pub(super) saved: bool,
    bucket: u32,
    at: u32,
}

impl Found {
    /// Before a seed's `find` (never written by `put`).
    pub(super) const NONE: Found = Found {
        level: 0,
        saved: false,
        bucket: 0,
        at: MISS,
    };
}

/// No NCBI counterpart: the latest insert of each diagonal (module comment).
pub(super) struct DiagLatest {
    buckets: Vec<Bucket>,
    bits: u32,
    /// `32 - bits`, the shift of `home`.
    shift: u32,
    stash: Vec<Entry>,
    /// The slots and stash entries that hold an entry, live or dead (at least the live entries).
    entries: usize,
    /// At least every level the map holds (`i32::MIN`: none since it was made or emptied).
    top: i32,
    /// Debug builds: the last S `put` saw in this stage (the order of (vi)).
    #[cfg(debug_assertions)]
    last_s: i32,
}

impl Default for DiagLatest {
    fn default() -> Self {
        DiagLatest {
            buckets: Vec::new(),
            bits: 0,
            shift: 32,
            stash: Vec::new(),
            entries: 0,
            top: i32::MIN,
            #[cfg(debug_assertions)]
            last_s: i32::MIN,
        }
    }
}

/// Bit k set where lane k of `lanes` is `v`.
#[inline(always)]
fn eq_mask(lanes: &[i32; 8], v: i32) -> u32 {
    lanes_mask::<false>(lanes, v)
}

/// Bit k set where lane k of `lanes` is below `v`.
#[inline(always)]
fn lt_mask(lanes: &[i32; 8], v: i32) -> u32 {
    lanes_mask::<true>(lanes, v)
}

/// Bit k set where lane k of `lanes` is below `v` (`LT`) or is `v`: two 128-bit compares and one
/// mask where the target has them (x86_64's SSE2, wasm32's simd128), else `scalar_mask`. The
/// compiler does not make the loop over the lanes into the compares (it kept about 30
/// instructions of scalar compares and shuffles per mask, against 10).
#[inline(always)]
fn lanes_mask<const LT: bool>(lanes: &[i32; 8], v: i32) -> u32 {
    #[cfg(target_arch = "x86_64")]
    {
        use core::arch::x86_64::{
            _mm_cmpeq_epi32, _mm_cmplt_epi32, _mm_loadu_si128, _mm_movemask_epi8, _mm_packs_epi16,
            _mm_packs_epi32, _mm_set1_epi32,
        };
        // SAFETY: SSE2 is part of every x86_64 target; the two loads read the 32 bytes of
        // `lanes`.
        unsafe {
            let v = _mm_set1_epi32(v);
            let lo = _mm_loadu_si128(lanes.as_ptr().cast());
            let hi = _mm_loadu_si128(lanes.as_ptr().add(4).cast());
            let (lo, hi) = if LT {
                (_mm_cmplt_epi32(lo, v), _mm_cmplt_epi32(hi, v))
            } else {
                (_mm_cmpeq_epi32(lo, v), _mm_cmpeq_epi32(hi, v))
            };
            // lanes of all ones or zeros, narrowed to one byte each (saturation keeps them)
            let words = _mm_packs_epi32(lo, hi);
            _mm_movemask_epi8(_mm_packs_epi16(words, words)) as u32 & 0xff
        }
    }
    #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
    {
        use core::arch::wasm32::{i32x4_bitmask, i32x4_eq, i32x4_lt, i32x4_splat, v128_load};
        // SAFETY: the two loads read the 32 bytes of `lanes` (no alignment needed).
        let (lo, hi) = unsafe {
            (
                v128_load(lanes.as_ptr().cast()),
                v128_load(lanes.as_ptr().add(4).cast()),
            )
        };
        let v = i32x4_splat(v);
        let (lo, hi) = if LT {
            (i32x4_lt(lo, v), i32x4_lt(hi, v))
        } else {
            (i32x4_eq(lo, v), i32x4_eq(hi, v))
        };
        u32::from(i32x4_bitmask(lo)) | (u32::from(i32x4_bitmask(hi)) << 4)
    }
    #[cfg(not(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    )))]
    {
        scalar_mask::<LT>(lanes, v)
    }
}

/// `lanes_mask` lane by lane (the targets without the compares, and the tests' reference).
#[cfg(any(
    test,
    not(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    ))
))]
#[inline(always)]
fn scalar_mask<const LT: bool>(lanes: &[i32; 8], v: i32) -> u32 {
    let mut mask = 0u32;
    for (k, &lane) in lanes.iter().enumerate() {
        mask |= u32::from(if LT { lane < v } else { lane == v }) << k;
    }
    mask
}

/// The level below which a slot is free at a seed of bound `bound`: the bound, and above
/// `i32::MIN` (an empty slot is always free; with no bound only empty slots are).
#[inline(always)]
fn free_below(bound: i32) -> i32 {
    bound.max(NO_LEVEL + 1)
}

impl DiagLatest {
    /// The start of a stage whose first position has the level bound `bound0`
    /// (`LvSeeds::bound` of the stage's hash offset): a new map gets its buckets; a map whose
    /// levels may not all lie below `bound0` is emptied ((vii): only after the offset's clear).
    pub(super) fn begin(&mut self, bound0: i32) {
        if self.buckets.is_empty() {
            self.alloc(MIN_BITS);
        } else if bound0 == i32::MIN || self.top >= bound0 {
            self.clear();
        }
        #[cfg(debug_assertions)]
        {
            self.last_s = i32::MIN;
        }
    }

    /// `2^bits` empty buckets and an empty stash.
    fn alloc(&mut self, bits: u32) {
        self.buckets.clear();
        self.buckets.resize(1 << bits, EMPTY);
        self.bits = bits;
        self.shift = 32 - bits;
        self.stash.clear();
        self.entries = 0;
        self.top = i32::MIN;
        MOST_BUCKETS.fetch_max(1 << bits, Ordering::Relaxed);
    }

    /// Every entry gone (the buckets kept).
    fn clear(&mut self) {
        self.buckets.fill(EMPTY);
        self.stash.clear();
        self.entries = 0;
        self.top = i32::MIN;
    }

    /// The home bucket of `diag`: the high bits of its product with a Fibonacci constant.
    #[inline(always)]
    fn home(&self, diag: i32) -> usize {
        debug_assert!((1..=32).contains(&self.bits) && self.shift == 32 - self.bits);
        ((diag as u32).wrapping_mul(0x9E37_79B1) >> self.shift) as usize
    }

    /// Brings `diag`'s home bucket toward the core (a hint only).
    #[inline(always)]
    pub(super) fn prefetch(&self, diag: i32) {
        prefetch(self.buckets.as_ptr().wrapping_add(self.home(diag)));
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:828-836
    /// The retrieve of the seed at S (`s`) on `diag` (module comment (iv)): the entry of the
    /// diagonal, which is its latest insert when that is live and something dead otherwise; no
    /// entry gives level 0 and no saved hit, as NCBI's retrieve without a cell. Remembers where
    /// the seed's insert goes.
    #[inline(always)]
    pub(super) fn find(&self, diag: i32, s: i32) -> Found {
        let b = self.home(diag);
        debug_assert!(b < self.buckets.len());
        // SAFETY: `home` is below 2^bits, the number of buckets (`begin` made them).
        let bucket = unsafe { self.buckets.get_unchecked(b) };
        let found = eq_mask(&bucket.diag, diag) & SLOT_MASK;
        if found != 0 {
            let lane = found.trailing_zeros();
            return Found {
                level: bucket.lvl[(lane & 7) as usize],
                saved: (bucket.lvl[SAVED] >> lane) & 1 != 0,
                bucket: b as u32,
                at: lane,
            };
        }
        if s < bucket.diag[UNTIL] {
            if let Some(i) = self.find_stash(diag) {
                let e = &self.stash[i];
                return Found {
                    level: e.level,
                    saved: e.saved,
                    bucket: b as u32,
                    at: STASHED + i as u32,
                };
            }
        }
        Found {
            bucket: b as u32,
            ..Found::NONE
        }
    }

    /// `find` in the stash (the home bucket has a stash entry that may be live): the place of
    /// `diag`'s entry. (A place, not a `Found`: a `Found` returned through memory kept every
    /// seed's `Found` on the stack, written in 4-byte fields and read back 8 bytes at a time.)
    #[cold]
    #[inline(never)]
    fn find_stash(&self, diag: i32) -> Option<usize> {
        self.stash.iter().position(|e| e.diag == diag)
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:939-941
    /// The insert of `(diag, level, saved)` by the seed at S (`s`, of level bound `bound`) whose
    /// `find` gave `found`: over the entry found, else into a free slot of the home bucket, else
    /// into the stash (module comment (iv)).
    #[inline(always)]
    pub(super) fn put(
        &mut self,
        found: Found,
        diag: i32,
        level: i32,
        saved: bool,
        s: i32,
        bound: i32,
    ) {
        #[cfg(debug_assertions)]
        {
            assert!(
                s >= self.last_s,
                "LOSAT_X_DCLATEST: a seed at {s} after one at {}",
                self.last_s
            );
            self.last_s = s;
        }
        debug_assert!(level >= 0);
        if self.top < level {
            self.top = level;
        }
        let b = found.bucket as usize;
        debug_assert!(b < self.buckets.len());
        let lane = if found.at < SLOTS {
            found.at
        } else if found.at == MISS {
            // SAFETY: as in `find`; nothing changed the map since.
            let bucket = unsafe { self.buckets.get_unchecked_mut(b) };
            let free = lt_mask(&bucket.lvl, free_below(bound)) & SLOT_MASK;
            if free == 0 {
                self.put_cold(b, Entry { diag, level, saved }, s, bound);
                return;
            }
            let lane = free.trailing_zeros();
            let slot = &mut bucket.diag[(lane & 7) as usize];
            self.entries += usize::from(*slot == NO_DIAG);
            *slot = diag;
            lane
        } else {
            self.put_stashed(found.at, b, Entry { diag, level, saved }, s, bound);
            return;
        };
        // SAFETY: as in `find`.
        let bucket = unsafe { self.buckets.get_unchecked_mut(b) };
        bucket.lvl[(lane & 7) as usize] = level;
        let bit = 1i32 << lane;
        bucket.lvl[SAVED] = if saved {
            bucket.lvl[SAVED] | bit
        } else {
            bucket.lvl[SAVED] & !bit
        };
    }

    /// `put` over the stash entry `find` found.
    #[cold]
    #[inline(never)]
    fn put_stashed(&mut self, at: u32, b: usize, e: Entry, s: i32, bound: i32) {
        self.stash[(at - STASHED) as usize] = e;
        self.raise_until(b, e.level, s, bound);
    }

    /// `stash_until` of bucket `b` raised past the last S at which an entry of level `level` is
    /// live, for the seeds' bounds (`bound` at `s`: the bound is S plus a constant of the fold's
    /// seeds; with no bound the entry may always be live).
    fn raise_until(&mut self, b: usize, level: i32, s: i32, bound: i32) {
        let until = if bound == i32::MIN {
            i32::MAX
        } else {
            // live while level >= S + (bound - s)
            (i64::from(level) - (i64::from(bound) - i64::from(s)) + 1).min(i64::from(i32::MAX))
                as i32
        };
        let at = &mut self.buckets[b].diag[UNTIL];
        if *at < until {
            *at = until;
        }
    }

    /// `put` where the home bucket `b` has no free slot: the stash entry of the diagonal (dead, or
    /// `find` would have found it), else a free stash entry, else a new one; when the stash is
    /// full and the buckets are fewer than twice the entries, the buckets double (`grow`) and the
    /// entry goes where the new map has room. (So the map holds at most about two buckets per
    /// entry: where many live diagonals share their home at every size, the stash holds them.)
    #[cold]
    #[inline(never)]
    fn put_cold(&mut self, b: usize, e: Entry, s: i32, bound: i32) {
        let free = free_below(bound);
        let at = self
            .stash
            .iter()
            .position(|x| x.diag == e.diag)
            .or_else(|| self.stash.iter().position(|x| x.level < free));
        if let Some(i) = at {
            self.stash[i] = e;
            self.raise_until(b, e.level, s, bound);
            return;
        }
        if self.stash.len() < STASH || self.buckets.len() >= 2 * self.entries {
            self.stash.push(e);
            self.entries += 1;
            self.raise_until(b, e.level, s, bound);
            return;
        }
        self.grow(s, bound);
        self.place(e, s, bound);
    }

    /// Twice the buckets, with the live entries (at the seed at `s` of bound `bound`) moved there
    /// and the dead ones dropped (module comment (iv)).
    #[cold]
    #[inline(never)]
    fn grow(&mut self, s: i32, bound: i32) {
        GROWS.fetch_add(1, Ordering::Relaxed);
        let free = free_below(bound);
        let mut live: Vec<Entry> = Vec::new();
        for bucket in &self.buckets {
            for k in 0..SLOTS as usize {
                if bucket.lvl[k] >= free {
                    live.push(Entry {
                        diag: bucket.diag[k],
                        level: bucket.lvl[k],
                        saved: (bucket.lvl[SAVED] >> k) & 1 != 0,
                    });
                }
            }
        }
        live.extend(self.stash.iter().filter(|e| e.level >= free));
        let top = self.top;
        debug_assert!(self.bits < 32);
        self.alloc(self.bits + 1);
        self.top = top;
        for e in live {
            self.place(e, s, bound);
        }
    }

    /// `e` into a free slot of its home bucket, else into the stash (after `grow`; no entry of its
    /// diagonal is in the map).
    fn place(&mut self, e: Entry, s: i32, bound: i32) {
        let b = self.home(e.diag);
        let bucket = &mut self.buckets[b];
        let free = lt_mask(&bucket.lvl, free_below(bound)) & SLOT_MASK;
        self.entries += 1;
        if free == 0 {
            self.stash.push(e);
            self.raise_until(b, e.level, s, bound);
            return;
        }
        let lane = free.trailing_zeros() as usize;
        bucket.diag[lane] = e.diag;
        bucket.lvl[lane] = e.level;
        let bit = 1i32 << lane;
        bucket.lvl[SAVED] = if e.saved {
            bucket.lvl[SAVED] | bit
        } else {
            bucket.lvl[SAVED] & !bit
        };
    }
}

/// No NCBI counterpart: counts a chunk of the strips stage that has the level bounds
/// (`x_pair_seed_disc.rs` `latest_bound`): folded on the maps when its scan ranges `ascend`, else
/// sent to the chains (`latest_ranges_ascend`).
pub(super) fn count_stage(ascend: bool) {
    let counter = if ascend { &STAGES } else { &BACK };
    counter.fetch_add(1, Ordering::Relaxed);
}

/// No NCBI counterpart: the fields of the `[TIMING] x_dc_latest:` line (after the stage's
/// `x_pairpardisc_seed:` line while the switch is on): the stages folded on the maps, the chunks
/// whose scan ranges step back (sent to the chains), the maps' growths and the most buckets a map
/// had.
pub(super) fn timing_fields() -> String {
    format!(
        " stages={} back={} grows={} buckets={}",
        STAGES.load(Ordering::Relaxed),
        BACK.load(Ordering::Relaxed),
        GROWS.load(Ordering::Relaxed),
        MOST_BUCKETS.load(Ordering::Relaxed)
    )
}

#[cfg(test)]
#[path = "x_dc_latest_tests.rs"]
mod tests;
