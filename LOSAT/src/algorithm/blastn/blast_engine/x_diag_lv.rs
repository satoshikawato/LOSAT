//! EXPERIMENT (LOSAT_X_DIAGLV, a strict-set switch since 2026-10-17): a seed's retrieve of its own
//! diagonal from the diagonal hash skips its read-only chain walk when every cell of the diagonal's
//! bucket has a level below the smallest level that could change a decision of the seed. The table
//! keeps, per bucket, an upper bound of the levels of its cells (`DiagHashTable::max_level`, 512
//! values).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:361-382
//! (`s_BlastDiagHashRetrieve`)
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:828-841,887
//! (`s_BlastnDiagHashExtendInitialHit`: the two-hit and the one-hit branch)
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:396-450
//! (`s_BlastDiagHashInsert`)
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:174-181
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:939-941
//! (`s_BlastnDiagHashExtendInitialHit`: the insert at the end of every seed)
//!
//! The bound of a seed at `s_off_pos` (`LvSeeds::bound`; `level_bound` for one seed) is
//! `s_off_pos + 1` (the `last_hit` test) and, with two hits, at most `s_end_pos - window_size`
//! = `s_off_pos + word_length - window_size` (the window test). Both are `s_off_pos` plus a
//! constant of the fold's seeds, so `LvSeeds::new` takes that constant (and the positions that
//! have a bound) once per chunk or fold call, and a seed's bound is one add and two compares.
//! For dc-megablast (window 40, template 16, 18 or 21) the window term binds, `s_off_pos - 24`
//! to `s_off_pos - 19`; with one hit (megablast, blastn, blastn-short) the bound is
//! `s_off_pos + 1`. `retrieve_below` returns `None` without a walk when the bucket's
//! `max_level` is below the bound, and is `retrieve` otherwise. The callers: the retrieves of the
//! two reference scans, `extend_offset_pairs` (`pair_extend.rs`, two-stage tables) and
//! `scan_na_lookup` (`na_scan.rs`, the other tables), not while the debug, trace or window
//! diagnostics print what they read; the one-hit fold of PAIRPAR (`x_pair_seed.rs` `fold`), the
//! two-hit fold of `LOSAT_X_PAIRPARDISC` (`x_pair_seed_disc.rs` `fold_two_hit`), the na fold
//! of `LOSAT_X_PAIRPARNA` (`x_pair_seed_na.rs` `fold_na_t`) and the bucket-split folds of
//! `LOSAT_X_PAIRPARFOLD` (`x_pair_seed_part.rs` `last_hit`). The
//! off-diagonal retrieves (`scan_range > 0`, an unported option) stay plain walks.
//!
//! Why the output is the same:
//! - (i) `max_level[b]` is at least the level of every cell reachable from the head of bucket
//!   `b`, at all times. Every level is written by `insert`, which raises `max_level` of the
//!   diagonal's bucket to `level` before its walk and writes `level` into one cell of that
//!   bucket's chain only (the cell of the same diagonal, a stale cell it takes over, both found
//!   by walking this chain, or the cell it appends and links at this bucket's head; a take-over
//!   keeps the cell in its bucket, since the walk visits only this bucket). The clear when the
//!   offset reaches `INT4_MAX` divided by 4 empties every chain and resets every bound to `i32::MIN`;
//!   moving the offset between subjects writes no level; the cells the chain grows by are not
//!   linked; `copy_hash` copies the bounds with the chains. The bound never decreases between
//!   clears, so it stays an upper bound when a cell's level is lowered.
//! - (ii) When `max_level[b] < bound`, a cell that the walk would return has a level
//!   `L <= max_level[b] < bound`, and it makes the decisions of no cell (`last_hit` and
//!   `hit_saved` both 0):
//!   - the `last_hit` test: `L < s_off_pos + 1`, so `s_off_pos < L` is false; with no cell
//!     `s_off_pos < 0` is false (`s_off_pos >= 0`; positions beyond `i32::MAX`, where the
//!     decisions' `as i32` would wrap, get no bound);
//!   - with two hits, the window test (`hit_saved` set, or `s_end_pos` past `last_hit` plus
//!     `window_size`), where `s_end_pos` is `s_off_pos` plus `word_length` at this point: a level
//!     `L` below `s_end_pos` minus `window_size` makes it true whatever `hit_saved` says, and with
//!     no cell it is true because `s_end_pos > window_size` (the hash offset starts at
//!     `window_size` and only grows or returns to it, so `s_off_pos >= window_size`; positions
//!     where that fails get no bound); with one hit (`window_size == 0`) the test is false either
//!     way and the `check_masks` branch runs;
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
//! EXPERIMENT (LOSAT_X_DIAGLV2, a strict-set switch since 2026-10-22): a retrieve that the bound
//! does not skip is answered without the walk when every cell of the bucket but the one the latest
//! insert wrote lies below the bound. A two-hit table (window above 0, dc-megablast) with level
//! bounds keeps, per bucket, the latest insert's `(diag, level, hit_len, hit_saved)` and `other`,
//! an upper bound of the levels of the bucket's other cells (`DiagHashTable::x_last`, `LvLast`).
//! Where `other` is below the seed's bound, the retrieve returns the latest insert's values for its
//! diagonal and no cell for any other diagonal. Only the strips folds of LOSAT_X_PAIRPARDISC keep
//! and read the memo: each of their stages makes it again from the bounds at its start
//! (`x_lv2_reset`: no latest insert, `other` the bucket's bound), so the reference fold
//! (`extend_offset_pairs`: a chunk without a pool with LOSAT_X_BNDCONE `0`, a masked chunk, and
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
//! (Debug builds assert the memo's validity, `DiagHashTable::x_last_fresh`.)

use super::{DiagHashTable, DIAGHASH_NUM_BUCKETS};
use std::sync::OnceLock;

/// No NCBI counterpart: reads LOSAT_X_DIAGLV (a strict-set switch: on unless it is `0`) once. It
/// does not change any value NCBI computes.
fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_DIAGLV"))
}

#[cfg(test)]
thread_local! {
    // The walks this thread skipped (tests only: a fold's skips, counted on its thread).
    static TEST_SKIPS: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
    // The memo's answers on this thread (tests only, LOSAT_X_DIAGLV2).
    static TEST_MEMOS: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
}

/// No NCBI counterpart: the bounds of a new table (`DiagHashTable::new`): one `i32::MIN` per
/// bucket while the switch is on, nothing while it is off.
pub(super) fn new_levels() -> Vec<i32> {
    if on() {
        vec![i32::MIN; DIAGHASH_NUM_BUCKETS]
    } else {
        Vec::new()
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
    /// `other` is its bound, so its memo answers only where the bound skips); the tests' check of
    /// (iv) does not expect a cell for it.
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
/// unless it is `0`; read once) are on and the table
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
    *ON.get_or_init(|| on() && crate::utils::xswitch::on("LOSAT_X_DIAGLV2"))
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:828-841
/// What the level bound and the decisions read of a fold's seeds, fixed for a chunk
/// (the reference fold) or a fold call (the strips folds) and made there once: whether the
/// seeds take two hits, the word length that `s_end` adds to `s_off`, and the window. A seed's
/// bound (module comment, (ii)) is `s_off_pos + add`, where `add` is 1 (the `last_hit` test),
/// with two hits at most `word_length - window_size` (the window test), for `s_off_pos` in
/// `lo..=hi`, and `i32::MIN` (nothing is below it) elsewhere: where `s_end_pos` passes
/// `i32::MAX` (the decisions' `as i32` would wrap) and, with two hits, where `s_end_pos <=
/// window_size` (with no cell the window test would be false).
#[derive(Clone, Copy, Debug)]
pub(super) struct LvSeeds {
    #[cfg(test)]
    two_hits: bool,
    #[cfg(test)]
    word_length: usize,
    #[cfg(test)]
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
            #[cfg(test)]
            two_hits,
            #[cfg(test)]
            word_length,
            #[cfg(test)]
            window_size,
            lo,
            hi,
            add,
        }
    }

    /// The level bound of a seed at `s_off_pos` (`s_off` plus `hash_table->offset`): a cell whose
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

    /// Tests only: the decisions of `s_BlastnDiagHashExtendInitialHit` that read the retrieve
    /// (module comment, (ii)), as the four callers take them from `last_hit` and `hit_saved` (`0`
    /// and `false` with no cell) before the seed is extended: rejected (`s_off_pos < last_hit`),
    /// and if not, whether the two-hit branch runs (`two_hits && (hit_saved || s_end_pos >
    /// last_hit + window_size)`, `s_end_pos = s_off_pos + word_length`).
    #[cfg(test)]
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
    /// `retrieve` for a seed of `seeds` at `s_off_pos`, or `None` without a walk when every cell
    /// of `diag`'s bucket has a level below the seed's bound (`max_level`). A table without
    /// bounds always walks. With `MEMO` (the folds of a dc-megablast strips stage only, after the
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
                    #[cfg(test)]
                    TEST_MEMOS.with(|memos| memos.set(memos.get() + 1));
                    return found;
                }
            }
        }
        self.retrieve(diag)
    }

    /// Tests only: bounds for a table made without them, from its cells.
    #[cfg(test)]
    pub(super) fn x_lv_start(&mut self) {
        self.max_level = vec![i32::MIN; self.num_buckets];
        for bucket in 0..self.num_buckets {
            let mut index = self.backbone[bucket];
            while index != 0 {
                let cell = &self.chain[index as usize];
                self.max_level[bucket] = self.max_level[bucket].max(cell.level);
                index = cell.next;
            }
        }
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
#[path = "x_diag_lv_tests.rs"]
mod tests;
