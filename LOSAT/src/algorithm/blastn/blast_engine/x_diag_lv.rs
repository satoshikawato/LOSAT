//! EXPERIMENT (LOSAT_X_DIAGLV, a strict-set switch since 2026-10-17; shadow
//! LOSAT_X_DIAGLVSHADOW, off unless set to any value): a seed's retrieve of its own diagonal from the diagonal hash skips its
//! read-only chain walk when every cell of the diagonal's bucket has a level below the smallest
//! level that could change a decision of the seed. The table keeps, per bucket, an upper bound
//! of the levels of its cells (`DiagHashTable::max_level`, 512 values).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:361-382
//! (`s_BlastDiagHashRetrieve`)
//! ```c
//!     Uint4 bucket = ((Uint4) diag * 0x9E370001) % DIAGHASH_NUM_BUCKETS;
//!     Uint4 index = table->backbone[bucket];
//!
//!     while (index) {
//!         if (table->chain[index].diag == diag) {
//!             *level = table->chain[index].level;
//!             *hit_len = table->chain[index].hit_len;
//!             *hit_saved = table->chain[index].hit_saved;
//!             return 1;
//!         }
//!
//!         index = table->chain[index].next;
//!     }
//!     return 0;
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:828-841,887
//! (`s_BlastnDiagHashExtendInitialHit`: the two-hit and the one-hit branch)
//! ```c
//!     diag = s_off - q_off;
//!     s_end = s_off + word_length;
//!     s_off_pos = s_off + hash_table->offset;
//!     s_end_pos = s_end + hash_table->offset;
//!
//!     rc = s_BlastDiagHashRetrieve(hash_table, diag, &last_hit, &s_l, &hit_saved);
//!
//!     /* if there is no record in hashtable, we set last_hit to be a very negative number */
//!     if(!rc)  last_hit = 0;
//!
//!     /* hit within the explored area should be rejected*/
//!     if (s_off_pos < last_hit) return 0;
//!
//!     if (two_hits && (hit_saved || s_end_pos > last_hit + window_size )) {
//! ...
//!     }  else if (check_masks) {
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:396-450
//! (`s_BlastDiagHashInsert`)
//! ```c
//!     Uint4 bucket = ((Uint4) diag * 0x9E370001) % DIAGHASH_NUM_BUCKETS;
//!     Uint4 index = table->backbone[bucket];
//! ...
//!         if (table->chain[index].diag == diag) {
//!             table->chain[index].level = level;
//! ...
//!             if ( s_off - table->chain[index].level > window_size) {
//!                 table->chain[index].diag = diag;
//!                 table->chain[index].level = level;
//! ...
//!     cell = table->chain + table->occupancy;
//!     cell->diag = diag;
//!     cell->level = level;
//! ...
//!     cell->next = table->backbone[bucket];
//!     table->backbone[bucket] = table->occupancy;
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:174-181
//! ```c
//!     } else if (ewp->hash_table) {
//!         if (ewp->hash_table->offset >= INT4_MAX / 4) {
//!             ewp->hash_table->occupancy = 1;
//!             ewp->hash_table->offset = ewp->hash_table->window;
//!             memset(ewp->hash_table->backbone, 0,
//!                    ewp->hash_table->num_buckets * sizeof(Int4));
//!         } else {
//!             ewp->hash_table->offset += subject_length + ewp->hash_table->window;
//! ```
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:939-941
//! (`s_BlastnDiagHashExtendInitialHit`: the insert at the end of every seed)
//! ```c
//!     s_BlastDiagHashInsert(hash_table, diag, s_end_pos,
//!                           (hit_ready) ? 0 : s_end_pos - s_off_pos,
//!                           hit_ready, s_off_pos, window_size + Delta + 1);
//! ```
//!
//! The bound of a seed at `s_off_pos` (`LvSeeds::bound`; `level_bound` for one seed) is
//! `s_off_pos + 1` (the `last_hit` test) and, with two hits, at most `s_end_pos - window_size`
//! = `s_off_pos + word_length - window_size` (the window test). Both are `s_off_pos` plus a
//! constant of the fold's seeds, so `LvSeeds::new` takes that constant (and the positions that
//! have a bound) once per chunk or fold call, and a seed's bound is one add and two compares.
//! For dc-megablast (window 40, template 16, 18 or 21) the window term binds, `s_off_pos - 24`
//! to `s_off_pos - 19`; with one hit (megablast, blastn, blastn-short) the bound is
//! `s_off_pos + 1`. `retrieve_below` returns `None` without a walk when the bucket's
//! `max_level` is below the bound, and is `retrieve` otherwise. The callers: the two retrieves
//! of `process_offset_pairs` (`run.rs`; not while the debug, trace or window diagnostics print
//! what they read), the one-hit fold of PAIRPAR (`x_pair_seed.rs` `fold`) and the
//! two-hit fold of `LOSAT_X_PAIRPARDISC` (`x_pair_seed_disc.rs` `fold_two_hit`). The
//! off-diagonal retrieves (`scan_range > 0`, an unported option) stay plain walks.
//!
//! Why the output is the same:
//! - (i) `max_level[b]` is at least the level of every cell reachable from the head of bucket
//!   `b`, at all times. Every level is written by `insert`, which raises `max_level` of the
//!   diagonal's bucket to `level` before its walk and writes `level` into one cell of that
//!   bucket's chain only (the cell of the same diagonal, a stale cell it takes over, both found
//!   by walking this chain, or the cell it appends and links at this bucket's head; a take-over
//!   keeps the cell in its bucket, since the walk visits only this bucket). The clear when the
//!   offset reaches `INT4_MAX / 4` empties every chain and resets every bound to `i32::MIN`;
//!   moving the offset between subjects writes no level; the cells the chain grows by are not
//!   linked; `copy_hash` copies the bounds with the chains. The bound never decreases between
//!   clears, so it stays an upper bound when a cell's level is lowered.
//! - (ii) When `max_level[b] < bound`, a cell that the walk would return has a level
//!   `L <= max_level[b] < bound`, and it makes the decisions of no cell (`last_hit = 0`,
//!   `hit_saved = 0`):
//!   - the `last_hit` test: `L < s_off_pos + 1`, so `s_off_pos < L` is false; with no cell
//!     `s_off_pos < 0` is false (`s_off_pos >= 0`; positions beyond `i32::MAX`, where the
//!     decisions' `as i32` would wrap, get no bound);
//!   - with two hits, the window test `hit_saved || s_end_pos > last_hit + window_size`, where
//!     `s_end_pos = s_off_pos + word_length` at this point: `L < s_end_pos - window_size`
//!     makes it true whatever `hit_saved` says, and with no cell it is true because
//!     `s_end_pos > window_size` (the hash offset starts at `window_size` and only grows or
//!     returns to it, so `s_off_pos >= window_size`; positions where that fails get no
//!     bound); with one hit (`window_size == 0`) the test is false either way and the
//!     `check_masks` branch runs;
//!   - nothing else reads what the retrieve returns: `s_l` (`hit_len`) is never read, and
//!     `last_hit` and `hit_saved` only by these two tests (na_ungapped.c:839, 841; the
//!     off-diagonal block at 853-886 makes retrieves of its own).
//!
//!   So the bound needs no term for the insert. The prototype's third term, the insert's stale
//!   limit `s_off_pos - (window_size + Delta + 1)`, only lowered the bound (it came from a
//!   variant whose skip also covered the insert's walk): the insert at the end of the seed
//!   (na_ungapped.c:939-941) takes `s_end_pos`, `hit_ready`, `s_off_pos` and the window, all of
//!   this seed and its extension, and `hit_ready` depends on the retrieve only through the two
//!   tests above; `s_BlastDiagHashInsert` walks the chain itself to find the cell it updates,
//!   takes over or appends, whether the retrieve walked or not; and the next seed reads only
//!   the table, which no retrieve writes.
//!
//!   The bound is tight: a cell at `s_off_pos + 1` is rejected where no cell is not, and with
//!   two hits a cell at `s_end_pos - window_size` with `hit_saved = 0` fails the window test
//!   where no cell passes it (test `the_bound_is_tight_at_every_reachable_setting`).
//! - (iii) `insert` is not changed and still walks the chain, so the table is the reference's
//!   cell for cell; only the retrieve's read-only walk is skipped. The decisions, the ungapped
//!   extensions, the hit list and the table are therefore the reference's element for element.
//!
//! Which programs reach it: BLASTN only (every task: megablast, dc-megablast, blastn,
//! blastn-short), `-subject` and `-db`, when the query block is longer than 8000 letters (the
//! diagonal hash; a shorter block uses the diagonal array, which this does not touch). No other
//! program has a `DiagHashTable` (BLASTP, BLASTX, TBLASTN and TBLASTX keep diagonal arrays of
//! their own). The Web reactors receive no environment, so they run the default set: this switch
//! and LOSAT_X_DIAGLV2 (below) are on there.
//!
//! With the switch off the table has no bounds (`max_level` empty), `insert` does one length
//! test of the empty vector, the clear fills an empty vector, and the callers call `retrieve`.
//!
//! The shadow (`LOSAT_X_DIAGLVSHADOW`) checks both parts on every skipped retrieve: it walks the
//! bucket and checks (i) for each cell against the bound, then makes the walk of the retrieve and
//! checks (ii), that the decisions the fold takes from the cell it returns (the `last_hit` test,
//! then with two hits the window test) are those from no cell.
//!
//! EXPERIMENT (LOSAT_X_DIAGLV2, a strict-set switch since 2026-10-22; its shadow is
//! LOSAT_X_DIAGLVSHADOW): a retrieve that the bound does not skip is answered without the walk
//! when every cell of the bucket but the one the latest insert wrote lies below the bound. A
//! two-hit table (window above 0, dc-megablast) with level bounds keeps, per bucket, the latest
//! insert's `(diag, level, hit_len, hit_saved)` and `other`, an upper bound of the levels of the
//! bucket's other cells (`DiagHashTable::x_last`, `LvLast`). Where `other` is below the seed's
//! bound, the retrieve returns the latest insert's values for its diagonal and no cell for any
//! other diagonal. Only the strips folds of LOSAT_X_PAIRPARDISC keep and read the memo: each of
//! their stages makes it again from the bounds at its start (`x_lv2_reset`: no latest insert,
//! `other` the bucket's bound), so the reference fold (`process_offset_pairs`: one thread, and
//! every worker with many subjects) runs the code before the switch and need not keep it.
//!
//! Why the output is the same:
//! - (iv) The latest insert's cell is the first cell of its diagonal in the bucket's chain, with
//!   the insert's values. `insert` walks the chain from its head and stops at the first cell of
//!   its diagonal or at the first stale cell, whichever comes first, and writes its values there
//!   (a take-over gives the stale cell the diagonal, and no cell before it has that diagonal), or
//!   appends a cell at the head. Only an insert of the same bucket changes a cell of the bucket
//!   or links a cell into its chain, and there is none after the latest; the growth of the chain
//!   and the merge of `LOSAT_X_PAIRPARFOLD` (which renumbers cells but keeps each chain's order)
//!   change neither. So the walk of the retrieve for that diagonal returns the latest insert's
//!   values, which the memo returns.
//! - (v) `other` is at least the level of every other cell of the bucket. Every insert first
//!   raises `other` to the level of the previous latest insert and then becomes the latest; a
//!   cell that the latest insert did not write was written by an earlier insert (whose level
//!   `other` took) or holds what it held before the bucket's first insert of the stage, which is
//!   at or below the bucket's bound when the stage made the memo (`x_lv2_reset`; nothing after
//!   the table was made or cleared: `LvLast::EMPTY`, `other = i32::MIN`); the tests'
//!   tables made from cells start at `LvLast::UNKNOWN` (`other = i32::MAX`), which never
//!   answers. The stage copies of `LOSAT_X_PAIRPARFOLD` (`rebuild`) copy the memo with the
//!   chains. `copy_hash`, the merge (per changed bucket, from the partition that owns it) and the
//!   clear keep it too, though no fold reads it before the next stage's reset: so that the tests
//!   find the stage's table equal to the serial fold's, memo included.
//! - (vi) For any other diagonal, the cell the walk would return is not the latest insert's, so
//!   by (v) its level is below the bound, and by (ii) it makes the decisions of no cell. The
//!   retrieve's other values are read by no decision ((ii), last item). `insert_memo` is
//!   `insert_plain` and the memo's upkeep, so the table is the reference's cell for cell.
//!
//! One-hit tables (megablast, blastn, blastn-short: window 0) keep no memo: their inserts would
//! pay for it where no measured fold gains (the memo is a dc-megablast measurement, line R
//! REPORT_R12 §2.3). The reference fold does not keep it either: on many short subjects (100
//! 16S rRNA queries against a 16S database, every worker in the reference fold) keeping and
//! testing it cost x1.008 of the wall time where the walks it saves are one cell long
//! (v4-A15 ab5). The reference fold, the one-hit strips folds, and with the switch off (or
//! LOSAT_X_DIAGLV off: `x_last` empty) the dc-megablast strips folds call
//! `retrieve_below::<false>` and `DiagHashTable::insert_plain`, the code before the switch.
//! (Debug builds assert the memo's validity, `DiagHashTable::x_last_fresh`.) The shadow
//! (LOSAT_X_DIAGLVSHADOW) checks every answer: the bucket's walk must give (iv) (the first cell of the latest insert's
//! diagonal holds its values) and (v) (every other cell at or below `other`), and the walk of the
//! retrieve must return what the memo returned for that diagonal, or for any other diagonal a
//! cell with the decisions of no cell; the summary line counts the answers (`memo=`) and those
//! for the latest insert's diagonal (`memofound=`).

use super::{DiagHashTable, DIAGHASH_NUM_BUCKETS};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

/// The switch's three states.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Mode {
    Off,
    On,
    Shadow,
}

/// No NCBI counterpart: reads LOSAT_X_DIAGLVSHADOW (off unless set to any value) and LOSAT_X_DIAGLV
/// (a strict-set switch: on unless it is `0`, or `LOSAT_X_DEFAULT=0` and it is not set) once; the
/// shadow implies the switch. They do not change any
/// value NCBI computes.
fn mode() -> Mode {
    static MODE: OnceLock<Mode> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_DIAGLVSHADOW").is_some() {
            Mode::Shadow
        } else if crate::utils::xswitch::on("LOSAT_X_DIAGLV") {
            Mode::On
        } else {
            Mode::Off
        }
    })
}

// No NCBI counterpart: the counters of the shadow's summary line: the retrieves, the skipped
// ones, those whose bucket passed the check of (i), those whose decisions passed the check of
// (ii), and of those, the ones whose walk found a cell of the diagonal (below the bound).
static RETRIEVES: AtomicU64 = AtomicU64::new(0);
static SKIPPED: AtomicU64 = AtomicU64::new(0);
static VERIFIED: AtomicU64 = AtomicU64::new(0);
static DECISIONS: AtomicU64 = AtomicU64::new(0);
static CELLS: AtomicU64 = AtomicU64::new(0);
// No NCBI counterpart (LOSAT_X_DIAGLV2): the memo's answers that the shadow checked, and of those
// the ones for the latest insert's diagonal.
static MEMO_ANSWERS: AtomicU64 = AtomicU64::new(0);
static MEMO_FOUND: AtomicU64 = AtomicU64::new(0);

#[cfg(test)]
thread_local! {
    // The walks this thread skipped (tests only: a fold's skips, counted on its thread).
    static TEST_SKIPS: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
    // The memo's answers on this thread (tests only, LOSAT_X_DIAGLV2).
    static TEST_MEMOS: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
}

/// No NCBI counterpart: the bounds and the shadow flag of a new table (`DiagHashTable::new`):
/// one `i32::MIN` per bucket while the switch is on, nothing while it is off.
pub(super) fn new_levels() -> (Vec<i32>, bool) {
    match mode() {
        Mode::Off => (Vec::new(), false),
        Mode::On => (vec![i32::MIN; DIAGHASH_NUM_BUCKETS], false),
        Mode::Shadow => (vec![i32::MIN; DIAGHASH_NUM_BUCKETS], true),
    }
}

/// No NCBI counterpart (LOSAT_X_DIAGLV2, module comment): the latest insert of one bucket and an
/// upper bound of the levels of the bucket's other cells.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct LvLast {
    pub(super) diag: i32,
    pub(super) level: i32,
    pub(super) hit_len: i32,
    pub(super) other: i32,
    pub(super) hit_saved: bool,
}

impl LvLast {
    /// No diagonal (`diag = s_off - q_off` with both offsets in `0..i32::MAX` is above
    /// `i32::MIN`): the latest of a bucket that had no insert since the table was made, cleared or
    /// its memo made again (`x_lv2_reset`). The output does not depend on it (such a bucket's
    /// `other` is its bound, so its memo answers only where the bound skips); the shadow's check
    /// of (iv) does not expect a cell for it.
    pub(super) const NO_DIAG: i32 = i32::MIN;
    /// A bucket with no cell: no insert yet (its `level` raises `other` to `i32::MIN` only).
    pub(super) const EMPTY: LvLast = LvLast {
        diag: LvLast::NO_DIAG,
        level: i32::MIN,
        hit_len: 0,
        other: i32::MIN,
        hit_saved: false,
    };
    /// Tests only: a bucket whose cells were not made by inserts that the memo saw; never
    /// answers (`other` is never below a bound).
    #[cfg(test)]
    pub(super) const UNKNOWN: LvLast = LvLast {
        other: i32::MAX,
        ..LvLast::EMPTY
    };

    /// The insert of `(diag, level, hit_len, hit_saved)` into this bucket becomes the latest;
    /// the previous latest is one of the other cells now (module comment, (v)).
    #[inline(always)]
    pub(super) fn note(&mut self, diag: i32, level: i32, hit_len: i32, hit_saved: bool) {
        self.other = self.other.max(self.level);
        self.diag = diag;
        self.level = level;
        self.hit_len = hit_len;
        self.hit_saved = hit_saved;
    }
}

/// No NCBI counterpart: the memo of a new table (`DiagHashTable::new`, window `window`): one
/// `LvLast::EMPTY` per bucket while LOSAT_X_DIAGLV and LOSAT_X_DIAGLV2 (a strict-set switch: on
/// unless it is `0`, or `LOSAT_X_DEFAULT=0` and it is not set; read once) are on and the table
/// takes two hits, nothing otherwise. The switch does not change any value NCBI computes.
pub(super) fn new_last(window: i32) -> Vec<LvLast> {
    if lv2_on() && window > 0 {
        vec![LvLast::EMPTY; DIAGHASH_NUM_BUCKETS]
    } else {
        Vec::new()
    }
}

/// No NCBI counterpart: LOSAT_X_DIAGLV2 and LOSAT_X_DIAGLV, which it needs, are on (read once).
pub(super) fn lv2_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| mode() != Mode::Off && crate::utils::xswitch::on("LOSAT_X_DIAGLV2"))
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:828-841
/// ```c
///     s_end = s_off + word_length;
///     s_off_pos = s_off + hash_table->offset;
///     s_end_pos = s_end + hash_table->offset;
/// ...
///     if (s_off_pos < last_hit) return 0;
///
///     if (two_hits && (hit_saved || s_end_pos > last_hit + window_size )) {
/// ```
/// What the level bound and the shadow's decisions read of a fold's seeds, fixed for a chunk
/// (the reference fold) or a fold call (the strips folds) and made there once: whether the
/// seeds take two hits, the word length that `s_end` adds to `s_off`, and the window. A seed's
/// bound (module comment, (ii)) is `s_off_pos + add`, where `add` is 1 (the `last_hit` test),
/// with two hits at most `word_length - window_size` (the window test), for `s_off_pos` in
/// `lo..=hi`, and `i32::MIN` (nothing is below it) elsewhere: where `s_end_pos` passes
/// `i32::MAX` (the decisions' `as i32` would wrap) and, with two hits, where `s_end_pos <=
/// window_size` (with no cell the window test would be false).
#[derive(Clone, Copy, Debug)]
pub(super) struct LvSeeds {
    two_hits: bool,
    word_length: usize,
    window_size: usize,
    lo: usize,
    hi: usize,
    add: i64,
}

impl LvSeeds {
    /// The seeds' constants: `two_hits` (`window_size > 0`), `word_length` (the length that
    /// `s_end` adds to `s_off`: the two-stage word length, or the lookup word without the
    /// two-stage table) and `window_size`.
    pub(super) fn new(two_hits: bool, word_length: usize, window_size: usize) -> Self {
        // The positions with a bound: `s_end_pos = s_off_pos + word_length` at most `i32::MAX`
        // (`word_length.max(1)`, so that the bound `s_off_pos + 1` fits too) and, with two hits,
        // above `window_size`. No position has one when `lo > hi`.
        let (lo, hi) = match (i32::MAX as usize).checked_sub(word_length.max(1)) {
            Some(hi) if two_hits => ((window_size + 1).saturating_sub(word_length), hi),
            Some(hi) => (0, hi),
            None => (1, 0),
        };
        let add = if two_hits {
            1.min(word_length as i64 - window_size as i64)
        } else {
            1
        };
        Self {
            two_hits,
            word_length,
            window_size,
            lo,
            hi,
            add,
        }
    }

    /// The level bound of a seed at `s_off_pos` (`s_off + hash_table->offset`): a cell whose
    /// level lies below it makes the decisions of no cell (module comment, (ii)). In `lo..=hi`
    /// it lies in `1..=i32::MAX` (`s_off_pos + 1 <= i32::MAX`; with two hits and `add < 1`,
    /// `s_off_pos + word_length - window_size >= 1` because `s_end_pos > window_size`).
    #[inline(always)]
    pub(super) fn bound(&self, s_off_pos: usize) -> i32 {
        if s_off_pos < self.lo || s_off_pos > self.hi {
            return i32::MIN;
        }
        let bound = s_off_pos as i64 + self.add;
        debug_assert!((1..=i32::MAX as i64).contains(&bound));
        bound as i32
    }

    /// The decisions of `s_BlastnDiagHashExtendInitialHit` that read the retrieve (module
    /// comment, (ii)), as the four callers take them from `last_hit` and `hit_saved` (`0` and
    /// `false` with no cell) before the seed is extended: rejected (`s_off_pos < last_hit`),
    /// and if not, whether the two-hit branch runs (`two_hits && (hit_saved || s_end_pos >
    /// last_hit + window_size)`, `s_end_pos = s_off_pos + word_length`).
    pub(super) fn decisions(
        &self,
        found: Option<(i32, i32, bool)>,
        s_off_pos: usize,
    ) -> (bool, bool) {
        let (last_hit, _hit_len, hit_saved) = found.unwrap_or((0, 0, false));
        let rejected = (s_off_pos as i32) < last_hit;
        let s_end_pos_i32 = (s_off_pos + self.word_length) as i32;
        let window_end = last_hit.wrapping_add(self.window_size as i32);
        let branch = self.two_hits && (hit_saved || s_end_pos_i32 > window_end);
        (rejected, !rejected && branch)
    }
}

/// The level bound of one seed (`LvSeeds::bound`; the folds make their `LvSeeds` once).
#[cfg(test)]
pub(super) fn level_bound(
    s_off_pos: usize,
    two_hits: bool,
    word_length: usize,
    window_size: usize,
) -> i32 {
    LvSeeds::new(two_hits, word_length, window_size).bound(s_off_pos)
}

impl DiagHashTable {
    /// No NCBI counterpart: whether this table keeps level bounds (LOSAT_X_DIAGLV was on when
    /// it was made).
    #[inline(always)]
    pub(super) fn x_lv_on(&self) -> bool {
        !self.max_level.is_empty()
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:368-381
    /// ```c
    ///     Uint4 bucket = ((Uint4) diag * 0x9E370001) % DIAGHASH_NUM_BUCKETS;
    ///     Uint4 index = table->backbone[bucket];
    ///
    ///     while (index) {
    ///         if (table->chain[index].diag == diag) {
    /// ```
    /// `retrieve` for a seed of `seeds` at `s_off_pos`, or `None` without a walk when every cell
    /// of `diag`'s bucket has a level below the seed's bound (`max_level`). A table without
    /// bounds always walks. With LOSAT_X_DIAGLVSHADOW every skipped walk is done and checked
    /// (`x_lv_verify`). With `MEMO` (the folds of a dc-megablast strips stage only, after the
    /// stage's `x_lv2_reset`: the field comment of `x_last`) the memo of LOSAT_X_DIAGLV2 may
    /// answer without the walk.
    #[inline(always)]
    pub(super) fn retrieve_below<const MEMO: bool>(
        &self,
        diag: i32,
        seeds: &LvSeeds,
        s_off_pos: usize,
    ) -> Option<(i32, i32, bool)> {
        let bucket =
            ((diag as u32).wrapping_mul(0x9E370001) % DIAGHASH_NUM_BUCKETS as u32) as usize;
        let bound = seeds.bound(s_off_pos);
        let below = matches!(self.max_level.get(bucket), Some(&max_level) if max_level < bound);
        if self.x_lv_shadow {
            self.x_lv_verify(diag, bucket, seeds, s_off_pos, bound, below);
        }
        if below {
            #[cfg(test)]
            TEST_SKIPS.with(|skips| skips.set(skips.get() + 1));
            return None;
        }
        // EXPERIMENT (LOSAT_X_DIAGLV2, module comment (iv)-(vi)): every cell but the latest
        // insert's lies below the bound; answer from the memo without the walk. (`insert_memo`
        // keeps the memo only on a table with bounds.) Without `MEMO` (the reference fold, the one-hit
        // strips folds, and dc-megablast's strips fold on a table without the memo) there is no
        // test, and the memo is not read.
        debug_assert!(self.x_last.is_empty() || self.x_lv_on());
        if MEMO {
            #[cfg(debug_assertions)]
            assert!(
                self.x_last.is_empty() || self.x_last_fresh,
                "LOSAT_X_DIAGLV2: a retrieve reads a stale memo (no x_lv2_reset since an insert without it)"
            );
            if let Some(last) = self.x_last.get(bucket) {
                if last.other < bound {
                    let found =
                        (last.diag == diag).then_some((last.level, last.hit_len, last.hit_saved));
                    if self.x_lv_shadow {
                        self.x_lv2_verify(diag, bucket, seeds, s_off_pos, found);
                    }
                    #[cfg(test)]
                    TEST_MEMOS.with(|memos| memos.set(memos.get() + 1));
                    return found;
                }
            }
        }
        self.retrieve(diag)
    }

    /// No NCBI counterpart: LOSAT_X_DIAGLVSHADOW for one answer of the memo (LOSAT_X_DIAGLV2).
    /// The bucket's walk must find the latest insert's values in the first cell of its diagonal
    /// ((iv)) and every other cell at or below `other` ((v)); the retrieve's walk must return
    /// the memo's answer for the latest insert's diagonal, and for any other diagonal a cell with
    /// the decisions of no cell ((vi)). Anything else prints the summary line with `DIFFERENT`
    /// and panics.
    #[cold]
    #[inline(never)]
    fn x_lv2_verify(
        &self,
        diag: i32,
        bucket: usize,
        seeds: &LvSeeds,
        s_off_pos: usize,
        found: Option<(i32, i32, bool)>,
    ) {
        let last = self.x_last[bucket];
        let mut first = true;
        let mut index = self.backbone[bucket];
        while index != 0 {
            let cell = &self.chain[index as usize];
            let ok = if first && cell.diag == last.diag {
                first = false;
                (cell.level, cell.hit_len, cell.hit_saved)
                    == (last.level, last.hit_len, last.hit_saved)
            } else {
                cell.level <= last.other
            };
            if !ok {
                print_line("DIFFERENT");
                panic!(
                    "LOSAT_X_DIAGLV2 (shadow): bucket {bucket} cell {index} (diagonal {}, level {}) breaks the memo {last:?}",
                    cell.diag, cell.level
                );
            }
            index = cell.next;
        }
        // (iv): a bucket with cells has had an insert since the clear, and its cell (unless the
        // memo was rebuilt from the bounds since, `x_lv2_reset`, and no insert came after)
        if first && self.backbone[bucket] != 0 && last.diag != LvLast::NO_DIAG {
            print_line("DIFFERENT");
            panic!(
                "LOSAT_X_DIAGLV2 (shadow): bucket {bucket} has no cell of the memo's diagonal, which breaks the memo {last:?}"
            );
        }
        let walked = self.retrieve(diag);
        let same = if found.is_some() {
            walked == found
        } else {
            seeds.decisions(walked, s_off_pos) == seeds.decisions(None, s_off_pos)
        };
        if !same {
            print_line("DIFFERENT");
            panic!(
                "LOSAT_X_DIAGLV2 (shadow): diagonal {diag} at s_off_pos {s_off_pos} ({seeds:?}): the memo {last:?} answers {found:?}, the walk {walked:?}"
            );
        }
        MEMO_ANSWERS.fetch_add(1, Ordering::Relaxed);
        if found.is_some() {
            MEMO_FOUND.fetch_add(1, Ordering::Relaxed);
        }
    }

    /// No NCBI counterpart: LOSAT_X_DIAGLVSHADOW for one retrieve. A skipped walk is walked
    /// here: every cell of the bucket must lie below the bound and at or below the bucket's
    /// `max_level` ((i): the walk would have returned no cell or a cell below the bound), and
    /// the decisions the seed takes from the cell the retrieve's walk returns must be those it
    /// takes from no cell ((ii): the bound is high enough to skip only what changes nothing).
    /// Anything else prints the summary line with `DIFFERENT` and panics.
    #[cold]
    #[inline(never)]
    fn x_lv_verify(
        &self,
        diag: i32,
        bucket: usize,
        seeds: &LvSeeds,
        s_off_pos: usize,
        bound: i32,
        below: bool,
    ) {
        RETRIEVES.fetch_add(1, Ordering::Relaxed);
        if !below {
            return;
        }
        SKIPPED.fetch_add(1, Ordering::Relaxed);
        let max_level = self.max_level[bucket];
        let mut index = self.backbone[bucket];
        while index != 0 {
            let cell = &self.chain[index as usize];
            if !(cell.level < bound && cell.level <= max_level) {
                print_line("DIFFERENT");
                panic!(
                    "LOSAT_X_DIAGLVSHADOW: bucket {bucket} cell {index} (diagonal {}) has level {}, not below the bound {bound} (bucket bound {max_level})",
                    cell.diag, cell.level
                );
            }
            index = cell.next;
        }
        VERIFIED.fetch_add(1, Ordering::Relaxed);
        let found = self.retrieve(diag);
        let walked = seeds.decisions(found, s_off_pos);
        let none = seeds.decisions(None, s_off_pos);
        if walked != none {
            print_line("DIFFERENT");
            panic!(
                "LOSAT_X_DIAGLVSHADOW: diagonal {diag} at s_off_pos {s_off_pos} (bound {bound}, {seeds:?}): the walk's cell {found:?} gives the decisions (rejected, two-hit branch) {walked:?}, no cell {none:?}"
            );
        }
        DECISIONS.fetch_add(1, Ordering::Relaxed);
        if found.is_some() {
            CELLS.fetch_add(1, Ordering::Relaxed);
        }
    }

    /// Tests only: bounds (and the shadow flag) for a table made without them, from its cells.
    #[cfg(test)]
    pub(super) fn x_lv_start(&mut self, shadow: bool) {
        self.max_level = vec![i32::MIN; self.num_buckets];
        for bucket in 0..self.num_buckets {
            let mut index = self.backbone[bucket];
            while index != 0 {
                let cell = &self.chain[index as usize];
                self.max_level[bucket] = self.max_level[bucket].max(cell.level);
                index = cell.next;
            }
        }
        self.x_lv_shadow = shadow;
    }

    /// No NCBI counterpart (LOSAT_X_DIAGLV2): the memo of every bucket made again from the level
    /// bounds of LOSAT_X_DIAGLV, at the start of a strips stage of dc-megablast (the only folds
    /// that keep and read it; the reference fold may have changed the table since the last
    /// stage): no latest insert (`LvLast::NO_DIAG`) and `other` the bucket's bound, which every
    /// cell's level is at or below. Such a memo answers no cell where `other` lies below the
    /// seed's bound, which is LOSAT_X_DIAGLV's skip ((i), (ii)); the stage's inserts make it
    /// the latest-insert memo of (iv) to (vi) again.
    pub(super) fn x_lv2_reset(&mut self) {
        debug_assert!(self.x_lv_on() && self.x_last.len() == self.max_level.len());
        for (last, &max_level) in self.x_last.iter_mut().zip(&self.max_level) {
            *last = LvLast {
                other: max_level,
                ..LvLast::EMPTY
            };
        }
        #[cfg(debug_assertions)]
        {
            self.x_last_fresh = true;
        }
    }

    /// Tests only: a memo (LOSAT_X_DIAGLV2) for a table with bounds, `LvLast::EMPTY` for a
    /// bucket without cells and `LvLast::UNKNOWN` (never answers) for one with cells.
    #[cfg(test)]
    pub(super) fn x_lv2_start(&mut self) {
        assert!(self.x_lv_on());
        #[cfg(debug_assertions)]
        {
            self.x_last_fresh = true;
        }
        self.x_last = (0..self.num_buckets)
            .map(|bucket| {
                if self.backbone[bucket] == 0 {
                    LvLast::EMPTY
                } else {
                    LvLast::UNKNOWN
                }
            })
            .collect();
    }
}

/// No NCBI counterpart: the shadow's line, `identical` at the end of a search.
fn print_line(verdict: &str) {
    eprintln!(
        "[X_SHADOW] DIAGLV retrieves={} skipped={} verified={} decisions={} cells={} memo={} memofound={} {verdict}",
        RETRIEVES.load(Ordering::Relaxed),
        SKIPPED.load(Ordering::Relaxed),
        VERIFIED.load(Ordering::Relaxed),
        DECISIONS.load(Ordering::Relaxed),
        CELLS.load(Ordering::Relaxed),
        MEMO_ANSWERS.load(Ordering::Relaxed),
        MEMO_FOUND.load(Ordering::Relaxed)
    );
}

/// No NCBI counterpart: the one summary line of LOSAT_X_DIAGLVSHADOW on stderr at the end of
/// the search (`run.rs` `search`); nothing without the shadow.
pub(super) fn summary() {
    if mode() == Mode::Shadow {
        print_line("identical");
    }
}

/// Tests only: the walks this thread has skipped so far.
#[cfg(test)]
pub(super) fn test_skips() -> u64 {
    TEST_SKIPS.with(|skips| skips.get())
}

/// Tests only: the retrieves the memo (LOSAT_X_DIAGLV2) answered on this thread so far.
#[cfg(test)]
pub(super) fn test_memos() -> u64 {
    TEST_MEMOS.with(|memos| memos.get())
}

#[cfg(test)]
mod tests {
    use super::super::{advance_diag_table_offset, diag_hash_insert_window, DiagStruct};
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

    fn bucket_of(diag: i32) -> usize {
        ((diag as u32).wrapping_mul(0x9E370001) % DIAGHASH_NUM_BUCKETS as u32) as usize
    }

    /// A table without bounds, or with bounds (and the shadow) from the start.
    fn table(window: i32, levels: bool, shadow: bool) -> DiagHashTable {
        let mut h = DiagHashTable::new(window);
        h.max_level = Vec::new();
        h.x_lv_shadow = false;
        h.x_last = Vec::new();
        if levels {
            h.x_lv_start(shadow);
        }
        h
    }

    type HashState = (Vec<u32>, Vec<(i32, i32, bool, i32, u32)>, u32, u32, i32);

    fn hash_state(h: &DiagHashTable) -> HashState {
        (
            h.backbone.clone(),
            h.chain[..h.occupancy as usize]
                .iter()
                .map(|c| (c.diag, c.level, c.hit_saved, c.hit_len, c.next))
                .collect(),
            h.occupancy,
            h.capacity,
            h.offset,
        )
    }

    /// (i): `max_level` bounds every reachable cell's level; a bucket's walk from its head.
    fn check_bounds(h: &DiagHashTable, step: usize) {
        for bucket in 0..h.num_buckets {
            let mut index = h.backbone[bucket];
            while index != 0 {
                let cell = &h.chain[index as usize];
                assert_eq!(bucket_of(cell.diag), bucket, "step {step}");
                assert!(
                    cell.level <= h.max_level[bucket],
                    "step {step}: bucket {bucket} cell {index} level {} above its bound {}",
                    cell.level,
                    h.max_level[bucket]
                );
                index = cell.next;
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:396-450
    // ```c
    //         if (table->chain[index].diag == diag) {
    //             table->chain[index].level = level;
    // ...
    //             if ( s_off - table->chain[index].level > window_size) {
    //                 table->chain[index].diag = diag;
    // ```
    // Random inserts (diagonals d + 512 k share a bucket, so cells are matched, taken over
    // and appended, the chain grows past 256 cells, levels go up and down, some negative),
    // offset advances between "subjects" and clears at INT4_MAX / 4: after every step the
    // bounds hold for every bucket, the table is the table without bounds, and after a clear
    // every bound is `i32::MIN`.
    #[test]
    fn the_bounds_hold_under_inserts_advances_and_clears() {
        let mut rng = Lcg(0xD1A6_17C0);
        for case in 0..30 {
            let window = [0, 16, 40][case % 3];
            let mut on = table(window, true, false);
            let mut off = table(window, false, false);
            let base = rng.below(2000) as i32 - 1000;
            // 24 diagonals in three buckets, and (in every other case) 600 more, so that the
            // chain grows past its 256 cells between clears
            let spread = if case % 2 == 0 { 600 } else { 16 };
            let diags: Vec<i32> = (0..24 + spread)
                .map(|k| {
                    if k < 24 {
                        base + 512 * (k % 8) + (k / 8) * 7
                    } else {
                        rng.below(1_000_000) as i32 - 500_000
                    }
                })
                .collect();
            let mut grown = false;
            let (mut hl_on, mut hl_off) = (Vec::<DiagStruct>::new(), Vec::<DiagStruct>::new());
            let (mut ln_on, mut ln_off) = (Vec::<u8>::new(), Vec::<u8>::new());
            let (mut dt_on, mut dt_off) = (window, window);
            let mut clears = 0;
            for step in 0..4000 {
                match rng.below(100) {
                    0..=1 => {
                        // the next subject; sometimes past INT4_MAX / 4 (the clear)
                        let len = if rng.below(4) == 0 {
                            (i32::MAX / 4) as usize
                        } else {
                            rng.below(5000) as usize
                        };
                        let cleared = on.offset >= i32::MAX / 4;
                        advance_diag_table_offset(
                            &mut dt_on, window, len, false, &mut hl_on, &mut ln_on, &mut on,
                        );
                        advance_diag_table_offset(
                            &mut dt_off,
                            window,
                            len,
                            false,
                            &mut hl_off,
                            &mut ln_off,
                            &mut off,
                        );
                        if cleared {
                            clears += 1;
                            assert_eq!(on.occupancy, 1);
                            assert!(on.max_level.iter().all(|&m| m == i32::MIN));
                        }
                    }
                    _ => {
                        let diag = diags[rng.below(diags.len() as u64) as usize];
                        let s_off = on.offset.saturating_add(rng.below(3000) as i32);
                        let level = s_off + rng.below(400) as i32 - 100;
                        let hit_len = rng.below(30) as i32;
                        let saved = rng.below(2) == 0;
                        let stale = rng.below(60) as i32 - 10;
                        on.insert_memo(diag, level, hit_len, saved, s_off, stale);
                        off.insert_memo(diag, level, hit_len, saved, s_off, stale);
                    }
                }
                grown |= on.capacity > 256;
                check_bounds(&on, step);
                assert_eq!(hash_state(&on), hash_state(&off), "case {case} step {step}");
                assert!(off.max_level.is_empty());
            }
            assert!(clears > 0, "case {case}: no clear");
            assert!(grown || spread < 600, "case {case}: the chain never grew");
        }
    }

    /// The decisions of `s_BlastnDiagHashExtendInitialHit` that read the retrieve (module
    /// comment, (ii)): rejected, and the two-hit branch taken.
    fn decisions(
        found: Option<(i32, i32, bool)>,
        s_off_pos: usize,
        two_hits: bool,
        word_length: usize,
        window_size: usize,
    ) -> (bool, bool) {
        let (last_hit, _hit_len, hit_saved) = found.unwrap_or((0, 0, false));
        let s_off_pos_i32 = s_off_pos as i32;
        let rejected = s_off_pos_i32 < last_hit;
        let s_end_pos_i32 = (s_off_pos + word_length) as i32;
        let branch = two_hits && (hit_saved || s_end_pos_i32 > last_hit + window_size as i32);
        (rejected, branch)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:833-841
    // ```c
    //     rc = s_BlastDiagHashRetrieve(hash_table, diag, &last_hit, &s_l, &hit_saved);
    //     if(!rc)  last_hit = 0;
    //     if (s_off_pos < last_hit) return 0;
    //     if (two_hits && (hit_saved || s_end_pos > last_hit + window_size )) {
    // ```
    // Random tables filled as a fold fills them (seed positions rising from the hash offset,
    // levels at or past the seed), probed at random diagonals and positions with the bounds
    // of one hit (megablast and blastn word lengths, window 0) and of two hits
    // (dc-megablast's window 40 and template lengths, `scan_range` 0 and 2): the decisions
    // with `retrieve_below` are those with `retrieve`, a walk that is not skipped returns the
    // same cell, and walks are skipped.
    #[test]
    fn retrieve_below_makes_the_decisions_of_retrieve() {
        let mut rng = Lcg(0x5EED_D1A6);
        let configs: [(bool, usize, usize, usize); 6] = [
            (false, 28, 0, 0),
            (false, 11, 0, 0),
            (true, 18, 40, 0),
            (true, 21, 40, 0),
            (true, 16, 40, 2),
            (true, 11, 40, 2),
        ];
        for (case, &(two_hits, word_length, window_size, scan_range)) in configs.iter().enumerate()
        {
            for shadow in [false, true] {
                let window = window_size as i32;
                let diag_hash_window =
                    diag_hash_insert_window(window_size, scan_range, word_length);
                let seeds = LvSeeds::new(two_hits, word_length, window_size);
                let mut h = table(window, true, shadow);
                h.offset = window + rng.below(100_000) as i32;
                let base = h.offset as usize;
                let diags: Vec<i32> = (0..3000)
                    .map(|_| rng.below(20_000) as i32 - 10_000)
                    .collect();
                let (mut skipped, mut probes, mut found_below) = (0u64, 0u64, 0u64);
                let mut s = base;
                for _ in 0..20_000 {
                    s += rng.below(3) as usize;
                    // a probe at the current position (as the fold's seeds come in order)
                    let diag = diags[rng.below(diags.len() as u64) as usize];
                    let bound = seeds.bound(s);
                    let full = h.retrieve(diag);
                    let before = test_skips();
                    let fast = h.retrieve_below::<true>(diag, &seeds, s);
                    let skip = test_skips() > before;
                    probes += 1;
                    if skip {
                        skipped += 1;
                        assert_eq!(fast, None);
                        if let Some((level, _, _)) = full {
                            assert!(level < bound);
                            found_below += 1;
                        }
                    } else {
                        assert_eq!(fast, full, "case {case}: a walk returned another cell");
                    }
                    assert_eq!(
                        decisions(fast, s, two_hits, word_length, window_size),
                        decisions(full, s, two_hits, word_length, window_size),
                        "case {case} diag {diag} s_off_pos {s} bound {bound}: {full:?}"
                    );
                    // then an insert as the fold makes one (first hit or extension end)
                    let reach = if rng.below(4) == 0 { 300 } else { 3 };
                    let level = (s + word_length + rng.below(reach) as usize) as i32;
                    let saved = rng.below(3) == 0;
                    let hit_len = if saved { 0 } else { word_length as i32 };
                    h.insert_memo(diag, level, hit_len, saved, s as i32, diag_hash_window);
                }
                assert!(
                    skipped * 5 > probes && found_below > 0,
                    "case {case}: too few skips ({skipped} of {probes}, {found_below} with a cell)"
                );
            }
        }
    }

    /// LOSAT_X_DIAGLV2 (iv) and (v): in every bucket the first cell of the latest insert's
    /// diagonal holds the latest insert's values and every other cell lies at or below `other`.
    fn check_memo(h: &DiagHashTable, at: &str) {
        for bucket in 0..h.num_buckets {
            let last = h.x_last[bucket];
            let mut first = true;
            let mut index = h.backbone[bucket];
            while index != 0 {
                let cell = &h.chain[index as usize];
                if first && cell.diag == last.diag {
                    first = false;
                    assert_eq!(
                        (cell.level, cell.hit_len, cell.hit_saved),
                        (last.level, last.hit_len, last.hit_saved),
                        "{at}: bucket {bucket} cell {index}"
                    );
                } else {
                    assert!(
                        cell.level <= last.other,
                        "{at}: bucket {bucket} cell {index} level {} above {last:?}",
                        cell.level
                    );
                }
                index = cell.next;
            }
            assert!(
                !first || h.backbone[bucket] == 0 || last.diag == LvLast::NO_DIAG,
                "{at}: bucket {bucket} has no cell of the memo's diagonal ({last:?})"
            );
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:368-381,396-450
    // ```c
    //     while (index) {
    //         if (table->chain[index].diag == diag) {
    // ...
    //             if ( s_off - table->chain[index].level > window_size) {
    // ```
    // LOSAT_X_DIAGLV2: tables filled as the dc-megablast fold fills them (seeds in rising
    // order, an insert after every probe; diagonals d + 512 k share a bucket, so cells are
    // matched, taken over and appended), with offset advances between "subjects" and clears at
    // INT4_MAX / 4, at the template lengths of dc-megablast and `scan_range` 0 and 2, with runs
    // of the reference fold (no memo read or kept: `retrieve_below::<false>`, `insert_plain`)
    // between strips stages (each makes the memo again, `x_lv2_reset`): the memo holds (iv) and
    // (v) in every stage, the table is the table without a memo, a memo answer for the latest
    // insert's diagonal is the walk's cell, every answer makes the decisions of the walk, and the
    // memo answers (also for the latest insert's diagonal) where the bound does not skip.
    #[test]
    fn the_memo_makes_the_decisions_of_retrieve() {
        let mut rng = Lcg(0x1A57_D1A6);
        let configs: [(usize, usize); 4] = [(16, 0), (18, 0), (21, 0), (11, 2)];
        for (case, &(word_length, scan_range)) in configs.iter().enumerate() {
            for shadow in [false, true] {
                let window_size = 40usize;
                let window = window_size as i32;
                let diag_hash_window =
                    diag_hash_insert_window(window_size, scan_range, word_length);
                let seeds = LvSeeds::new(true, word_length, window_size);
                let mut h = table(window, true, shadow);
                h.x_lv2_start();
                let mut plain = table(window, false, false);
                let base = rng.below(2000) as i32 - 1000;
                let diags: Vec<i32> = (0..2400)
                    .map(|k| {
                        if k < 48 {
                            base + 512 * (k % 16) + (k / 16) * 7
                        } else {
                            rng.below(20_000) as i32 - 10_000
                        }
                    })
                    .collect();
                let (mut hl_on, mut hl_off) = (Vec::<DiagStruct>::new(), Vec::<DiagStruct>::new());
                let (mut ln_on, mut ln_off) = (Vec::<u8>::new(), Vec::<u8>::new());
                let (mut dt_on, mut dt_off) = (window, window);
                let (mut memos, mut memo_found, mut clears) = (0u64, 0u64, 0u64);
                // steps left in a run of the reference fold, and the runs and stage starts
                let (mut reference_left, mut runs, mut resets) = (0usize, 0u64, 0u64);
                let mut s = h.offset as usize;
                for step in 0..30_000 {
                    if reference_left == 0 && rng.below(700) == 0 {
                        reference_left = 1 + rng.below(400) as usize;
                        runs += 1;
                    } else if reference_left == 0 && rng.below(1000) == 0 {
                        // a stage after a stage
                        h.x_lv2_reset();
                        resets += 1;
                        check_memo(&h, &format!("case {case} step {step} reset"));
                    }
                    if rng.below(1000) == 0 {
                        // the next subject; sometimes past INT4_MAX / 4 (the clear)
                        let len = if rng.below(3) == 0 {
                            (i32::MAX / 4) as usize
                        } else {
                            rng.below(5000) as usize
                        };
                        let cleared = h.offset >= i32::MAX / 4;
                        advance_diag_table_offset(
                            &mut dt_on, window, len, false, &mut hl_on, &mut ln_on, &mut h,
                        );
                        advance_diag_table_offset(
                            &mut dt_off,
                            window,
                            len,
                            false,
                            &mut hl_off,
                            &mut ln_off,
                            &mut plain,
                        );
                        if cleared {
                            clears += 1;
                            assert!(h.x_last.iter().all(|&l| l == LvLast::EMPTY));
                        }
                        s = h.offset as usize;
                    }
                    s += rng.below(3) as usize;
                    let diag = diags[rng.below(diags.len() as u64) as usize];
                    let full = plain.retrieve(diag);
                    assert_eq!(h.retrieve(diag), full, "case {case} step {step}");
                    let before = test_memos();
                    let fast = if reference_left > 0 {
                        h.retrieve_below::<false>(diag, &seeds, s)
                    } else {
                        h.retrieve_below::<true>(diag, &seeds, s)
                    };
                    if reference_left > 0 {
                        assert_eq!(test_memos(), before, "case {case} step {step}: memo read");
                    } else if test_memos() > before {
                        memos += 1;
                        if fast.is_some() {
                            memo_found += 1;
                            assert_eq!(fast, full, "case {case} step {step}: the memo's cell");
                        }
                    }
                    assert_eq!(
                        seeds.decisions(fast, s),
                        seeds.decisions(full, s),
                        "case {case} step {step} diag {diag} s_off_pos {s}: {fast:?} {full:?}"
                    );
                    let reach = if rng.below(4) == 0 { 300 } else { 3 };
                    let level = (s + word_length + rng.below(reach) as usize) as i32;
                    let saved = rng.below(3) == 0;
                    let hit_len = if saved { 0 } else { word_length as i32 };
                    if reference_left > 0 {
                        h.insert_plain(diag, level, hit_len, saved, s as i32, diag_hash_window);
                        reference_left -= 1;
                        if reference_left == 0 {
                            // the next strips stage
                            h.x_lv2_reset();
                            resets += 1;
                            check_memo(&h, &format!("case {case} step {step} reset"));
                        }
                    } else {
                        h.insert_memo(diag, level, hit_len, saved, s as i32, diag_hash_window);
                    }
                    plain.insert_plain(diag, level, hit_len, saved, s as i32, diag_hash_window);
                    if step % 97 == 0 {
                        if reference_left == 0 {
                            check_memo(&h, &format!("case {case} step {step}"));
                        }
                        assert_eq!(
                            hash_state(&h),
                            hash_state(&plain),
                            "case {case} step {step}"
                        );
                    }
                }
                assert!(clears > 0, "case {case}: no clear");
                assert!(
                    runs > 10 && resets > runs,
                    "case {case}: {runs} runs of the reference fold, {resets} resets"
                );
                assert!(
                    memos > 1000 && memo_found > 100,
                    "case {case}: too few memo answers ({memos}, {memo_found} for the latest insert)"
                );
            }
        }
    }

    // LOSAT_X_DIAGLV2, debug builds: a read of the memo after an insert that did not keep it
    // (the reference fold's `insert_plain`) and before a stage made it again (`x_lv2_reset`)
    // panics; after the reset it answers.
    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "stale memo")]
    fn a_memo_read_after_a_memo_free_insert_panics_in_debug_builds() {
        let seeds = LvSeeds::new(true, 18, 40);
        let mut h = table(40, true, false);
        h.x_lv2_start();
        h.insert_memo(7, 500, 0, true, 480, 41);
        h.insert_plain(7 + 512, 490, 0, true, 470, 41);
        h.x_lv2_reset();
        let _ = h.retrieve_below::<true>(7, &seeds, 480);
        h.insert_plain(7, 520, 0, true, 500, 41);
        let _ = h.retrieve_below::<true>(7, &seeds, 480);
    }

    #[test]
    fn the_bound_is_the_smallest_limit_and_none_where_the_argument_fails() {
        // one hit: s_off_pos + 1
        assert_eq!(level_bound(1000, false, 28, 0), 1001);
        assert_eq!(level_bound(1000, false, 11, 0), 1001);
        // two hits (window 40): the window limit s_end_pos - window_size lies below
        // s_off_pos + 1 (template 18: 1000 + 18 - 40; word 11: 1000 + 11 - 40)
        assert_eq!(level_bound(1000, true, 18, 40), 978);
        assert_eq!(level_bound(1000, true, 11, 40), 971);
        // a word longer than the window (not reachable): s_off_pos + 1
        assert_eq!(level_bound(1000, true, 50, 40), 1001);
        // s_end_pos <= window_size: no bound
        assert_eq!(level_bound(10, true, 18, 40), i32::MIN);
        assert_eq!(level_bound(22, true, 18, 40), i32::MIN);
        assert_eq!(level_bound(23, true, 18, 40), 1);
        // past i32::MAX: no bound
        assert_eq!(level_bound(i32::MAX as usize, false, 28, 0), i32::MIN);
        assert_eq!(level_bound(usize::MAX, false, 28, 0), i32::MIN);
        assert_eq!(level_bound(i32::MAX as usize - 27, false, 28, 0), i32::MIN);
        assert_eq!(
            level_bound(i32::MAX as usize - 28, false, 28, 0),
            i32::MAX - 27
        );
        assert_eq!(
            level_bound(i32::MAX as usize - 28, true, 28, 40),
            i32::MAX - 40
        );
        // an empty bucket (`i32::MIN`) is never below `i32::MIN`
        let h = table(40, true, false);
        assert_eq!(
            h.retrieve_below::<true>(5, &LvSeeds::new(true, 18, 40), 10),
            None
        );
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:839-841,939-941
    // ```c
    //     if (s_off_pos < last_hit) return 0;
    //
    //     if (two_hits && (hit_saved || s_end_pos > last_hit + window_size )) {
    // ...
    //     s_BlastDiagHashInsert(hash_table, diag, s_end_pos,
    //                           (hit_ready) ? 0 : s_end_pos - s_off_pos,
    //                           hit_ready, s_off_pos, window_size + Delta + 1);
    // ```
    // The bound at its edge, for every setting the command line reaches: two hits with window
    // 40 and templates 16, 18 and 21 (`scan_range` 0), and one hit with words 7, 11 and 28. One
    // cell on the probed diagonal at `bound - 1` (with `hit_saved` 0 and 1) is skipped, also
    // by the shadow, and decides as no cell; at `bound` it is walked, and with `hit_saved` 0
    // (with one hit, either) it decides otherwise than no cell, so a bound one higher would
    // change a decision. The insert the fold makes at that seed then leaves the same table with
    // the bounds as without, and the bounds hold.
    #[test]
    fn the_bound_is_tight_at_every_reachable_setting() {
        let configs: [(bool, usize, usize); 6] = [
            (true, 16, 40),
            (true, 18, 40),
            (true, 21, 40),
            (false, 7, 0),
            (false, 11, 0),
            (false, 28, 0),
        ];
        let effective = |(rejected, branch): (bool, bool)| (rejected, !rejected && branch);
        for &(two_hits, word_length, window_size) in &configs {
            let seeds = LvSeeds::new(two_hits, word_length, window_size);
            let diag_hash_window = diag_hash_insert_window(window_size, 0, word_length);
            let window = window_size as i32;
            for s_off_pos in [window_size, window_size + 1, 5_000, 1 << 20] {
                let bound = seeds.bound(s_off_pos);
                let expected = if two_hits {
                    s_off_pos + word_length - window_size
                } else {
                    s_off_pos + 1
                };
                assert_eq!(bound as usize, expected);
                assert_eq!(
                    bound,
                    level_bound(s_off_pos, two_hits, word_length, window_size)
                );
                for (level, hit_saved) in [
                    (bound - 1, false),
                    (bound - 1, true),
                    (bound, false),
                    (bound, true),
                ] {
                    let at = format!(
                        "two_hits {two_hits} word {word_length} s_off_pos {s_off_pos} bound {bound} level {level} hit_saved {hit_saved}"
                    );
                    let diag = 4321;
                    let hit_len = if hit_saved { 0 } else { word_length as i32 };
                    // one cell, written by an earlier seed of this diagonal
                    let earlier = (level as usize).saturating_sub(word_length) as i32;
                    let mut on = table(window, true, true);
                    let mut off = table(window, false, false);
                    on.insert_memo(diag, level, hit_len, hit_saved, earlier, diag_hash_window);
                    off.insert_memo(diag, level, hit_len, hit_saved, earlier, diag_hash_window);
                    let full = off.retrieve(diag);
                    assert_eq!(full, Some((level, hit_len, hit_saved)), "{at}");
                    let before = test_skips();
                    let fast = on.retrieve_below::<true>(diag, &seeds, s_off_pos);
                    let skipped = test_skips() > before;
                    let walked = decisions(full, s_off_pos, two_hits, word_length, window_size);
                    let none = decisions(None, s_off_pos, two_hits, word_length, window_size);
                    // the shadow's decisions are those of the test's
                    assert_eq!(seeds.decisions(full, s_off_pos), effective(walked), "{at}");
                    assert_eq!(seeds.decisions(None, s_off_pos), effective(none), "{at}");
                    if level < bound {
                        assert!(skipped, "{at}: not skipped");
                        assert_eq!(fast, None, "{at}");
                        assert_eq!(walked, none, "{at}");
                    } else {
                        assert!(!skipped, "{at}: skipped");
                        assert_eq!(fast, full, "{at}");
                        if !hit_saved || !two_hits {
                            assert_ne!(effective(walked), effective(none), "{at}: not tight");
                        }
                    }
                    // the insert at the end of this seed (unless it was rejected)
                    if !walked.0 {
                        let s_end_pos = (s_off_pos + word_length) as i32;
                        let ready = walked.1;
                        let len = if ready { 0 } else { word_length as i32 };
                        for h in [&mut on, &mut off] {
                            h.insert_memo(
                                diag,
                                s_end_pos,
                                len,
                                ready,
                                s_off_pos as i32,
                                diag_hash_window,
                            );
                        }
                    }
                    check_bounds(&on, 0);
                    assert_eq!(hash_state(&on), hash_state(&off), "{at}");
                }
            }
        }
    }

    #[test]
    #[should_panic(expected = "LOSAT_X_DIAGLVSHADOW")]
    fn the_shadow_panics_on_a_cell_above_the_bound() {
        let mut h = table(40, true, true);
        h.insert_memo(7, 500, 0, true, 480, 41);
        // a bound that lies wrong about the bucket
        h.max_level[bucket_of(7)] = 100;
        // one hit at 199: the seed's bound is 200
        let _ = h.retrieve_below::<true>(7 + 512, &LvSeeds::new(false, 28, 0), 199);
    }

    // LOSAT_X_DIAGLV2's shadow: a memo whose `other` lies wrong about an older cell of the
    // bucket (a cell at 980 above `other` 900) is caught on the first answer.
    #[test]
    #[should_panic(expected = "breaks the memo")]
    fn the_shadow_panics_on_a_memo_that_misses_a_cell() {
        let mut h = table(40, true, true);
        h.x_lv2_start();
        h.insert_memo(7, 980, 0, true, 960, 41);
        // same bucket; at 1000 the cell at 980 is not stale (1000 - 980 <= 41): appended
        h.insert_memo(7 + 512, 1500, 0, true, 1000, 41);
        h.x_last[bucket_of(7)].other = 900;
        // two hits, template 18, at 1000: the bound is 978, above `other`, below max_level
        let _ = h.retrieve_below::<true>(7, &LvSeeds::new(true, 18, 40), 1000);
    }

    // The shadow's check of (ii): with a one-hit bound one too high (`s_off_pos + 2`), a cell
    // at `s_off_pos + 1` passes the check of (i) but is rejected where no cell is not.
    #[test]
    #[should_panic(expected = "gives the decisions")]
    fn the_shadow_panics_on_a_bound_that_changes_a_decision() {
        let seeds = LvSeeds {
            add: 2,
            ..LvSeeds::new(false, 28, 0)
        };
        let mut h = table(0, true, true);
        h.insert_memo(7, 1001, 0, true, 980, 1);
        let _ = h.retrieve_below::<true>(7, &seeds, 1000);
    }

    // The same with a two-hit bound one too high (`s_end_pos - window_size + 1`, template 18):
    // a cell at `s_end_pos - window_size` with `hit_saved` 0 fails the window test where no
    // cell passes it.
    #[test]
    #[should_panic(expected = "gives the decisions")]
    fn the_shadow_panics_on_a_two_hit_bound_that_changes_a_decision() {
        let right = LvSeeds::new(true, 18, 40);
        let seeds = LvSeeds {
            add: right.add + 1,
            ..right
        };
        let mut h = table(40, true, true);
        h.insert_memo(7, 1000 + 18 - 40, 18, false, 960, 41);
        let _ = h.retrieve_below::<true>(7, &seeds, 1000);
    }
}
