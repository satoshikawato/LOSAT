//! EXPERIMENT (LOSAT_X_DIAGLV, shadow LOSAT_X_DIAGLVSHADOW; off unless set to any value, not in
//! the strict set): a seed's retrieve of its own diagonal from the diagonal hash skips its
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
//! The bound of a seed at `s_off_pos` (`level_bound`) is the smallest of `s_off_pos + 1` (the
//! `last_hit` test), `s_off_pos - (window_size + Delta + 1)` (the insert's stale limit; it is
//! not needed for the argument below and only lowers the bound, so fewer walks are skipped)
//! and, with two hits, `s_off_pos + word_length - window_size` (the window test), in i64 and
//! clamped to i32. `retrieve_below` returns `None` without a walk when the bucket's
//! `max_level` is below the bound, and is `retrieve` otherwise. The callers: the two retrieves
//! of `process_offset_pairs` (`run.rs`; not while the debug, trace or window diagnostics print
//! what they read), the one-hit fold of `LOSAT_X_PAIRPAR` (`x_pair_seed.rs` `fold`) and the
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
//!     `s_end_pos = s_off_pos + word_length` at this point: `L < s_off_pos + word_length -
//!     window_size` makes it true whatever `hit_saved` says, and with no cell it is true
//!     because `s_end_pos > window_size` (the hash offset starts at `window_size` and only
//!     grows or returns to it, so `s_off_pos >= window_size`; `level_bound` gives no bound
//!     where that fails); with one hit (`window_size == 0`) the test is false either way and
//!     the `check_masks` branch runs;
//!   - nothing else reads this retrieve's `last_hit`, `hit_saved` or `hit_len` (the
//!     off-diagonal retrieves are their own walks, and the insert walks the chain itself).
//! - (iii) `insert` is not changed and still walks the chain, so the table is the reference's
//!   cell for cell; only the retrieve's read-only walk is skipped. The decisions, the ungapped
//!   extensions, the hit list and the table are therefore the reference's element for element.
//!
//! Which programs reach it: BLASTN only (every task: megablast, dc-megablast, blastn,
//! blastn-short), `-subject` and `-db`, when the query block is longer than 8000 letters (the
//! diagonal hash; a shorter block uses the diagonal array, which this does not touch). No other
//! program has a `DiagHashTable` (BLASTP, BLASTX, TBLASTN and TBLASTX keep diagonal arrays of
//! their own). The Web reactors receive no environment, so the switch is never on there.
//!
//! With the switch off the table has no bounds (`max_level` empty), `insert` does one length
//! test of the empty vector, the clear fills an empty vector, and the callers call `retrieve`.

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

/// No NCBI counterpart: reads LOSAT_X_DIAGLVSHADOW and LOSAT_X_DIAGLV once (off unless set to
/// any value; not in the strict set); the shadow implies the switch. They do not change any
/// value NCBI computes.
fn mode() -> Mode {
    static MODE: OnceLock<Mode> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_DIAGLVSHADOW").is_some() {
            Mode::Shadow
        } else if std::env::var_os("LOSAT_X_DIAGLV").is_some() {
            Mode::On
        } else {
            Mode::Off
        }
    })
}

// No NCBI counterpart: the counters of the shadow's summary line.
static RETRIEVES: AtomicU64 = AtomicU64::new(0);
static SKIPPED: AtomicU64 = AtomicU64::new(0);
static VERIFIED: AtomicU64 = AtomicU64::new(0);

#[cfg(test)]
thread_local! {
    // The walks this thread skipped (tests only: a fold's skips, counted on its thread).
    static TEST_SKIPS: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
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

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:825,828-841,939-941
/// ```c
///     Int4 Delta = MIN(word_params->options->scan_range, window_size - word_length);
/// ...
///     s_off_pos = s_off + hash_table->offset;
///     s_end_pos = s_end + hash_table->offset;
/// ...
///     if (s_off_pos < last_hit) return 0;
///
///     if (two_hits && (hit_saved || s_end_pos > last_hit + window_size )) {
/// ...
///     s_BlastDiagHashInsert(hash_table, diag, s_end_pos,
///                           (hit_ready) ? 0 : s_end_pos - s_off_pos,
///                           hit_ready, s_off_pos, window_size + Delta + 1);
/// ```
/// The bound of a seed at `s_off_pos` (`s_off + hash_table->offset`): a cell whose level lies
/// below it makes the decisions of no cell (module comment, (ii)). `diag_hash_window` is the
/// insert's `window_size + Delta + 1` (`diag_hash_insert_window`), `word_length` the length
/// `s_end` adds to `s_off`. `i32::MIN` (nothing is below it) where the argument does not hold:
/// a position past `i32::MAX`, or two hits with `s_end_pos <= window_size`.
#[inline(always)]
pub(super) fn level_bound(
    s_off_pos: usize,
    diag_hash_window: i32,
    two_hits: bool,
    word_length: usize,
    window_size: usize,
) -> i32 {
    let s_end_pos = s_off_pos.saturating_add(word_length);
    if s_end_pos > i32::MAX as usize {
        return i32::MIN;
    }
    let s_off_pos = s_off_pos as i64;
    let mut bound = (s_off_pos + 1).min(s_off_pos - diag_hash_window as i64);
    if two_hits {
        let window_bound = s_end_pos as i64 - window_size as i64;
        if window_bound <= 0 {
            return i32::MIN;
        }
        bound = bound.min(window_bound);
    }
    bound.clamp(i32::MIN as i64, i32::MAX as i64) as i32
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
    /// `retrieve`, or `None` without a walk when every cell of `diag`'s bucket has a level
    /// below `bound` (`max_level`). A table without bounds always walks. With
    /// LOSAT_X_DIAGLVSHADOW every skipped walk is done and checked (`x_lv_verify`).
    #[inline(always)]
    pub(super) fn retrieve_below(&self, diag: i32, bound: i32) -> Option<(i32, i32, bool)> {
        let bucket =
            ((diag as u32).wrapping_mul(0x9E370001) % DIAGHASH_NUM_BUCKETS as u32) as usize;
        let below = matches!(self.max_level.get(bucket), Some(&max_level) if max_level < bound);
        if self.x_lv_shadow {
            self.x_lv_verify(bucket, bound, below);
        }
        if below {
            #[cfg(test)]
            TEST_SKIPS.with(|skips| skips.set(skips.get() + 1));
            return None;
        }
        self.retrieve(diag)
    }

    /// No NCBI counterpart: LOSAT_X_DIAGLVSHADOW for one retrieve. A skipped walk is walked
    /// here: every cell of the bucket must lie below the bound and at or below the bucket's
    /// `max_level` (so the walk would have returned no cell or a cell below the bound); any
    /// other cell prints the summary line with `DIFFERENT` and panics.
    #[cold]
    #[inline(never)]
    fn x_lv_verify(&self, bucket: usize, bound: i32, below: bool) {
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
}

/// No NCBI counterpart: the shadow's line, `identical` at the end of a search.
fn print_line(verdict: &str) {
    eprintln!(
        "[X_SHADOW] DIAGLV retrieves={} skipped={} verified={} {verdict}",
        RETRIEVES.load(Ordering::Relaxed),
        SKIPPED.load(Ordering::Relaxed),
        VERIFIED.load(Ordering::Relaxed)
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
                        on.insert(diag, level, hit_len, saved, s_off, stale);
                        off.insert(diag, level, hit_len, saved, s_off, stale);
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
                    let bound =
                        level_bound(s, diag_hash_window, two_hits, word_length, window_size);
                    let full = h.retrieve(diag);
                    let before = test_skips();
                    let fast = h.retrieve_below(diag, bound);
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
                    h.insert(diag, level, hit_len, saved, s as i32, diag_hash_window);
                }
                assert!(
                    skipped * 5 > probes && found_below > 0,
                    "case {case}: too few skips ({skipped} of {probes}, {found_below} with a cell)"
                );
            }
        }
    }

    #[test]
    fn the_bound_is_the_smallest_limit_and_none_where_the_argument_fails() {
        // one hit: s_off_pos + 1 (the insert's limit 1 - word_length lies above it)
        let w28 = diag_hash_insert_window(0, 0, 28);
        assert_eq!(level_bound(1000, w28, false, 28, 0), 1001);
        // two hits (window 40, template 18): the insert's limit, 41 below
        let w18 = diag_hash_insert_window(40, 0, 18);
        assert_eq!(w18, 41);
        assert_eq!(level_bound(1000, w18, true, 18, 40), 959);
        // two hits with a window limit below the insert's (word 11, window 40, Delta 2 -> 43;
        // window limit 1000 + 11 - 40 = 971 > 957)
        assert_eq!(
            level_bound(1000, diag_hash_insert_window(40, 2, 11), true, 11, 40),
            957
        );
        // and where only the window limit is smallest (a window above the insert's window)
        assert_eq!(level_bound(1000, -100, true, 18, 40), 978);
        // s_end_pos <= window_size: no bound
        assert_eq!(level_bound(10, w18, true, 18, 40), i32::MIN);
        assert_eq!(level_bound(22, w18, true, 18, 40), i32::MIN);
        assert_ne!(level_bound(23, w18, true, 18, 40), i32::MIN);
        // past i32::MAX: no bound
        assert_eq!(level_bound(i32::MAX as usize, w28, false, 28, 0), i32::MIN);
        assert_eq!(level_bound(usize::MAX, w28, false, 28, 0), i32::MIN);
        assert_eq!(
            level_bound(i32::MAX as usize - 28, w28, false, 28, 0),
            i32::MAX - 27
        );
        // an empty bucket (`i32::MIN`) is never below `i32::MIN`
        let h = table(40, true, false);
        assert_eq!(h.retrieve_below(5, i32::MIN), None);
    }

    #[test]
    #[should_panic(expected = "LOSAT_X_DIAGLVSHADOW")]
    fn the_shadow_panics_on_a_cell_above_the_bound() {
        let mut h = table(40, true, true);
        h.insert(7, 500, 0, true, 480, 41);
        // a bound that lies wrong about the bucket
        h.max_level[bucket_of(7)] = 100;
        let _ = h.retrieve_below(7 + 512, 200);
    }
}
