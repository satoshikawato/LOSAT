//! The seed-stage state of a TBLASTX query batch is freed as
//! soon as nothing can read it any more, before the HSP list of the batch's last subject is
//! merged and linked: the group tables of the unit path when the units of that subject have run,
//! the batch's lookup table and the worker's diagonal table when its seed stage has ended.
//!
//! This file is a child module of `blast_engine/run_impl.rs` (declared there with `#[path]`, as
//! `x_pair_par.rs` is), so it names the unit path's type without changing it.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1002-1004
//! ```c
//! if ((status = BlastExtendWordNew(query->length, word_params,
//!                                 &aux_struct->ewp)) != 0)
//!   return status;
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1604-1605
//! ```c
//! word_params = BlastInitialWordParametersFree(word_params);
//! s_BlastCoreAuxStructFree(aux_struct);
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:110-113
//! ```c
//! static BlastCoreAuxStruct*
//! s_BlastCoreAuxStructFree(BlastCoreAuxStruct* aux_struct)
//! {
//!     BlastExtendWordFree(aux_struct->ewp);
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:214-215
//! ```c
//! retval->m_InternalData->m_LookupTable.Reset
//!     (new TLookupTableWrap(lut, LookupTableWrapFree));
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:239-242
//! ```c
//! case eAaLookupTable:
//!    lookup->lut = (void*)
//!       BlastAaLookupTableDestruct((BlastAaLookupTable*)lookup->lut);
//!    break;
//! ```
//! NCBI creates the diagonal table (`BlastCoreAuxStruct.ewp`) before the subject loop of the
//! preliminary search and frees it after the loop (blast_engine.c:1605, after every subject has
//! been linked inside the loop); its lookup table belongs to the search's internal data and is
//! freed with it by `LookupTableWrapFree`. Both live through every `BLAST_LinkHsps` call. When
//! LOSAT frees its copies is not a value NCBI computes, and nothing NCBI reports depends on it.
//!
//! # What it does
//!
//! In LOSAT the HSP list of a subject is made and linked inside the per-subject closure of the
//! batch function (`process_subject`): the units of the unit path (or the reference frame loop)
//! find the HSPs, the unit lists are merged into the subject's list, and the list is linked. The
//! heap of a large search peaks in that merge and in the first linking call, and neither reads a
//! table the scan used. AP027131 x AP027133 (genetic code 4, one thread), bytes requested under
//! massif, before -> after this module: end of the units 110.4 -> 102.4 MiB, the merge 155.8 ->
//! 81.0 MiB, the first linking call 167.6 -> 107.3 MiB, which is the peak of the run in both. Of
//! the 167.6 MiB, 51.8 MiB are group tables and 8.0 MiB the worker's diagonal table.
//!
//! The batch function moves its tables into two holders (`XSeedBatch`), each an `XSeedState`
//! with the same rule: `lookup` holds the batch's lookup table (the full table of the reference
//! frame loop, or the stub of PAIRLUT), `units` the unit path with its group tables.
//! Every call of `process_subject` opens one stage of each (`enter`) as its first statements. It
//! ends the units stage (`XSeedStage::end`) when the subject's units have run, before their lists
//! are merged, and the lookup stage after the seed stage, before the first linking call; a stage
//! dropped without `end` (the return for a subject too short to search, an unwinding panic) ends
//! in `Drop`. Each state counts its stages that have not ended, starting from the number of
//! subjects of the batch. The stage end that brings a count to zero takes that state's tables
//! out of it, which frees them unless another stage still holds them.
//!
//! The worker's diagonal table (`WorkerState.diag_array`) is read and written only by the
//! reference frame loop. (a) A batch that the unit path takes
//! never runs that loop, so its workers are created with a table of zero cells
//! (`XSeedBatch::worker_table_cells`); (b) otherwise the worker whose stage end released the
//! batch's lookup table drops its own table there (`drop_worker_table`). There is no per-subject
//! release in the reference loop: a batch of many subjects would pay an allocation of the whole
//! table per subject.
//!
//! Known limits, accepted:
//! - The gain is for the last subject of a batch (the only one of a pairwise search). While
//!   earlier subjects are merged and linked, later ones still need the tables.
//! - The lookup stage ends after the frame loop, which names the lookup table. Where PAIRLUT's
//!   fallback builds the groups from the full table (`lut_build` gave none for a batch that takes
//!   the units), that table is the full one, not a stub, and it stays alive through the merge of
//!   the unit lists: that path keeps the higher seed-stage peak.
//!
//! # Why the result is the same, and why there is no comparison mode
//!
//! 1. The release changes no value and no order of computation: it only decides when three
//!    allocations are freed. The searches run the same statements on the same tables.
//! 2. The tables are reachable only through a stage. `XSeedBatch::new` takes them by value (the
//!    batch function keeps no binding of them), each slot is private to this module, and the
//!    only access a state hands out is `XSeedStage::tables`.
//! 3. No read after an end inside `process_subject`: `tables` borrows the stage and `end` takes
//!    the stage by value, so the borrow checker rejects any use of `x_pair_par_ref` (or of
//!    anything borrowed from it) after `x_units_stage.end()`, and any use of `lookup_ref` after
//!    `x_seed_stage.end()`. The merge that follows the units stage is a function without `self`
//!    (`XPairPar::merge_subject`): it gets the units' own lists by value and cannot name a group
//!    table.
//! 4. A stage that is still open keeps its tables alive through its own `Arc`, whatever the
//!    count says: taking the state's `Arc` out of the slot cannot free tables another call of
//!    `process_subject` is reading.
//! 5. A call that arrives after a release finds the slot empty, and `enter` stops the run with
//!    a message naming the subject. It cannot search without tables.
//! 6. So a wrong count can only keep tables longer (too high: they are freed with the state;
//!    too low with stages open at once: freed when the last of those ends) or stop the run (too
//!    low with a later call). It cannot change output. The two counts are independent.
//! 7. The diagonal table is indexed, through a raw pointer, only by the reference frame loop.
//!    That loop runs zero times exactly when the subject's frame slice is empty, and
//!    `process_subject` asserts before the loop, unconditionally, that the slice is empty or the
//!    table has its full size. A worker whose table was dropped by (b) gets no further subject
//!    unless the count was too low, and then `enter` (item 5) stops the run before anything else.
//!
//! Items 2-5 and 7 hold by ownership and by assertions that always run, not by comparing two
//! computations, and there is no second result to compare with: the same code produced it.
//! Hence no shadow mode. `LOSAT_X_STATS=1` prints one `[X_STATS] SEEDRELEASE` line at exit
//! (`print_summary`, called by `main`).

use super::x_pair_par::XPairPar;
use crate::algorithm::tblastx::lookup::BlastAaLookupTable;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering as AtomicOrdering};
use std::sync::{Arc, Mutex, MutexGuard, PoisonError};

// No NCBI counterpart: counters for LOSAT_X_STATS (batches searched; batches
// whose lookup table a stage end released; batches whose unit path, with its group tables, a
// stage end released; batches whose workers hold no diagonal table; worker tables dropped at the
// release). They do not change any value NCBI computes.
static BATCHES: AtomicU64 = AtomicU64::new(0);
static TABLES_RELEASED: AtomicU64 = AtomicU64::new(0);
static UNIT_TABLES_RELEASED: AtomicU64 = AtomicU64::new(0);
static BATCHES_WITHOUT_WORKER_TABLE: AtomicU64 = AtomicU64::new(0);
static WORKER_TABLES_RELEASED: AtomicU64 = AtomicU64::new(0);

/// The `LOSAT_X_STATS` line (printed by `main` at exit).
// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
pub fn print_summary() {
    if std::env::var_os("LOSAT_X_STATS").is_some() {
        let load = |c: &AtomicU64| c.load(AtomicOrdering::Relaxed);
        eprintln!(
            "[X_STATS] SEEDRELEASE batches={} tables_released={} unit_tables_released={} batches_without_worker_table={} worker_tables_released={}",
            load(&BATCHES),
            load(&TABLES_RELEASED),
            load(&UNIT_TABLES_RELEASED),
            load(&BATCHES_WITHOUT_WORKER_TABLE),
            load(&WORKER_TABLES_RELEASED)
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:503-505
// ```c
//     /* scan the subject sequence for hits */
//     hits = scansub(lookup_wrap, subject,
//                               offset_pairs, array_size, scan_range);
// ```
// The tables the seed stage of a batch scans with, in two holders: `lookup` is NCBI's
// `lookup_wrap->lut` (the full table of the reference frame loop, or the stub of
// PAIRLUT), read by the frame loop; `units` is the unit path with its per-group tables
// (x_pair_par.rs), read while the units of a subject run. Nothing after the seed stage reads
// either, and nothing after a subject's units reads the second.
/// The seed-stage tables of one query batch.
pub(super) struct XSeedBatch {
    /// The batch's lookup table.
    pub(super) lookup: XSeedState<BlastAaLookupTable>,
    /// The unit path of the batch (`None`: the reference frame loop takes every subject).
    pub(super) units: XSeedState<Option<XPairPar>>,
    /// The unit path takes every subject of the batch: the reference frame
    /// loop, the only reader of a worker's diagonal table, runs for no subject.
    units_only: bool,
}

impl XSeedBatch {
    // No NCBI counterpart: whether the reference frame loop can run for this batch is read
    // from the unit path before it moves into its holder (dispatch point 2 of PAIRPAR in
    // run_impl.rs: the loop runs zero times exactly when the unit path exists). One batch is
    // counted once, whatever the number of holders; the release of the unit path is counted only
    // for a batch that has one.
    /// The holders of a batch of `subjects` subjects.
    pub(super) fn new(
        lookup: BlastAaLookupTable,
        pair_par: Option<XPairPar>,
        subjects: usize,
    ) -> Self {
        let units_only = pair_par.is_some();
        BATCHES.fetch_add(1, AtomicOrdering::Relaxed);
        if units_only {
            BATCHES_WITHOUT_WORKER_TABLE.fetch_add(1, AtomicOrdering::Relaxed);
        }
        let units_released = pair_par.is_some().then_some(&UNIT_TABLES_RELEASED);
        Self {
            lookup: XSeedState::new(lookup, subjects, Some(&TABLES_RELEASED)),
            units: XSeedState::new(pair_par, subjects, units_released),
            units_only,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-60
    // ```c
    // diag_array_length = 1;
    // /* What power of 2 is just longer than the query? */
    // while (diag_array_length < (qlen+window_size))
    // {
    //         diag_array_length = diag_array_length << 1;
    // }
    // ...
    // diag_table->diag_array_length = diag_array_length;
    // ```
    // `cells` is that length. NCBI always allocates it (the word finder always runs on it). In
    // LOSAT the unit path scans on tables of its own, so a worker of a batch that never runs the
    // reference frame loop holds a table nobody reads; it is not allocated.
    /// The number of cells a worker's diagonal table is created with: `cells`, or 0 when no
    /// worker of the batch runs the reference frame loop.
    pub(super) fn worker_table_cells(&self, cells: usize) -> usize {
        worker_table_cells(self.units_only, cells)
    }
}

// No NCBI counterpart: the rule of `XSeedBatch::worker_table_cells` as a function of its two
// inputs, so that it is tested without a unit path.
fn worker_table_cells(units_only: bool, cells: usize) -> usize {
    if units_only {
        0
    } else {
        cells
    }
}

// No NCBI counterpart: NCBI's lookup table belongs to the search's internal data for the whole
// search (blast_aux_priv.cpp:214-215); this holder exists so that LOSAT can free its tables
// earlier (module doc). Generic over the tables: the batch uses it twice (`XSeedBatch`), and
// its rules are tested without a lookup table.
/// Tables of one batch and the count of their stages that have not ended.
pub(super) struct XSeedState<T> {
    /// The batch's hold on the tables; `None` once the last stage end has taken it.
    slot: Mutex<Option<Arc<T>>>,
    /// Calls of `process_subject` that have not ended their stage of this state.
    pending: AtomicUsize,
    /// The LOSAT_X_STATS counter of this state's release, if it has one.
    released: Option<&'static AtomicU64>,
}

impl<T> XSeedState<T> {
    // No NCBI counterpart: see the type.
    /// A state holding `tables` for `stages` stages; `released` is the counter that the release
    /// adds one to.
    pub(super) fn new(tables: T, stages: usize, released: Option<&'static AtomicU64>) -> Self {
        Self {
            slot: Mutex::new(Some(Arc::new(tables))),
            pending: AtomicUsize::new(stages),
            released,
        }
    }

    // No NCBI counterpart: the slot is only cloned or taken under the lock, and neither can
    // panic, so a poisoned lock still holds a valid value; taking it avoids a second panic
    // when a stage ends while its thread unwinds.
    fn lock(&self) -> MutexGuard<'_, Option<Arc<T>>> {
        self.slot.lock().unwrap_or_else(PoisonError::into_inner)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1429-1432
    // ```c
    // if (seq_arg.seq->length < min_subj_seq_length) {
    //     BlastSeqSrcReleaseSequence(seq_src, &seq_arg);
    //     continue;
    // }
    // ```
    // One stage per iteration of NCBI's subject loop, opened before this length test (the stage
    // of a subject that is not searched ends when it is dropped). NCBI's iterations all read the
    // one lookup table of the search; here each takes its own hold on it.
    /// Opens the stage of subject `s_idx`: a hold on the tables until the stage ends.
    /// Stops the run when the tables were already released (module doc, item 5).
    pub(super) fn enter(&self, s_idx: usize) -> XSeedStage<'_, T> {
        // The guard is a temporary of this statement: the lock is free again before the check.
        let tables = self.lock().clone();
        let Some(tables) = tables else {
            panic!(
                "the seed-stage tables of the batch were released before the seed stage of subject {s_idx} started"
            );
        };
        XSeedStage {
            state: self,
            tables: Some(tables),
        }
    }

    // No NCBI counterpart: NCBI frees its tables after the subject loop (blast_engine.c:1605;
    // lookup_wrap.c:239-242). This is the earlier release of LOSAT's copies: the count goes down
    // by one, never below zero, and the stage end that leaves no announced stage pending takes
    // the batch's hold out of the slot. The tables are freed when the last hold is dropped: here,
    // outside the lock, unless a stage beyond the announced number is still open.
    /// One stage has ended. True when this call took the batch's hold on the tables.
    fn stage_ended(&self) -> bool {
        let before = match self.pending.fetch_update(
            AtomicOrdering::SeqCst,
            AtomicOrdering::SeqCst,
            |pending| Some(pending.saturating_sub(1)),
        ) {
            Ok(before) | Err(before) => before,
        };
        if before > 1 {
            return false;
        }
        let taken = self.lock().take();
        let released = taken.is_some();
        drop(taken);
        if released {
            if let Some(counter) = self.released {
                counter.fetch_add(1, AtomicOrdering::Relaxed);
            }
        }
        released
    }
}

// No NCBI counterpart: see `XSeedState::enter`.
/// One subject's stage of a state: the only way to reach that state's tables.
pub(super) struct XSeedStage<'a, T> {
    state: &'a XSeedState<T>,
    /// This stage's hold; `None` once the stage has ended.
    tables: Option<Arc<T>>,
}

impl<T> XSeedStage<'_, T> {
    // No NCBI counterpart: the borrow of the stage is what ends every use of the tables at
    // `end` (module doc, item 3).
    /// The tables, borrowed from the stage.
    pub(super) fn tables(&self) -> &T {
        self.tables
            .as_deref()
            .expect("a stage holds the tables until it ends")
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-873
    // ```c
    // if (hit_params->link_hsp_params) {
    //     status = BLAST_LinkHsps(program_number, hsp_list_out, query_info,
    //               subject->length, gap_align->sbp, hit_params->link_hsp_params,
    //               score_options->gapped_calculation);
    // ```
    // Called before this linking call, which reads neither the lookup table nor the diagonal
    // table (NCBI keeps both allocated through it): the units stage when the subject's units
    // have run, the lookup stage after the seed stage.
    /// Ends the stage. True when this was the stage end that released the state's tables (the
    /// last one of the batch); for the lookup stage the caller then drops
    /// its diagonal table.
    pub(super) fn end(mut self) -> bool {
        self.finish()
    }

    // No NCBI counterpart: the stage's own hold goes first, so that taking the batch's hold in
    // `stage_ended` frees the tables when no other stage is open.
    fn finish(&mut self) -> bool {
        match self.tables.take() {
            Some(own) => {
                drop(own);
                self.state.stage_ended()
            }
            None => false,
        }
    }
}

// No NCBI counterpart: a stage that is dropped without `end` (the return for a subject shorter
// than a word of codons plus two, blast_engine.c:1429-1432; an unwinding panic) ends here.
impl<T> Drop for XSeedStage<'_, T> {
    fn drop(&mut self) {
        self.finish();
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:203-212
// ```c
// Blast_ExtendWord *BlastExtendWordFree(Blast_ExtendWord * ewp)
// {
//
//     if (ewp == NULL)
//         return NULL;
//
//     s_BlastDiagTableFree(ewp->diag_table);
//     s_BlastDiagHashFree(ewp->hash_table);
//     sfree(ewp);
//     return NULL;
// ```
// NCBI frees the diagonal table after the subject loop (blast_engine.c:1605). The worker that
// ended the last lookup stage of the batch frees its table here, before it links that subject;
// it gets no further subject (module doc, item 7).
/// Drops a worker's diagonal table (called when the lookup stage's `end` returned true).
pub(super) fn drop_worker_table<C>(table: &mut Vec<C>) {
    if table.capacity() != 0 {
        WORKER_TABLES_RELEASED.fetch_add(1, AtomicOrdering::Relaxed);
    }
    *table = Vec::new();
}

#[cfg(test)]
mod tests {
    // No NCBI counterpart: NCBI frees its tables after the subject loop and has no unit test of
    // an earlier release. These tests check the ownership rules of the module doc on a payload
    // that counts its drops.
    use super::*;

    /// A stand-in for the tables: counts how often it is dropped.
    struct Payload {
        drops: Arc<AtomicUsize>,
        value: u32,
    }

    impl Drop for Payload {
        fn drop(&mut self) {
            self.drops.fetch_add(1, AtomicOrdering::SeqCst);
        }
    }

    fn new_state(stages: usize) -> (XSeedState<Payload>, Arc<AtomicUsize>) {
        let drops = Arc::new(AtomicUsize::new(0));
        let payload = Payload {
            drops: Arc::clone(&drops),
            value: 7,
        };
        (XSeedState::new(payload, stages, None), drops)
    }

    fn freed(drops: &AtomicUsize) -> usize {
        drops.load(AtomicOrdering::SeqCst)
    }

    #[test]
    fn one_subject_holds_the_tables_until_its_stage_ends() {
        let (state, drops) = new_state(1);
        let stage = state.enter(0);
        assert_eq!(stage.tables().value, 7);
        assert_eq!(freed(&drops), 0);
        assert!(stage.end());
        assert_eq!(freed(&drops), 1);
        drop(state);
        assert_eq!(freed(&drops), 1);
    }

    #[test]
    fn subjects_in_order_free_the_tables_at_the_last_end() {
        let (state, drops) = new_state(4);
        for s_idx in 0..4 {
            let stage = state.enter(s_idx);
            assert_eq!(stage.tables().value, 7);
            assert_eq!(freed(&drops), 0);
            assert_eq!(stage.end(), s_idx == 3, "subject {s_idx}");
            assert_eq!(freed(&drops), usize::from(s_idx == 3), "subject {s_idx}");
        }
    }

    #[test]
    fn a_stage_dropped_without_end_counts() {
        // Not the last one: the next stage's `end` is the release.
        let (state, drops) = new_state(2);
        drop(state.enter(0));
        assert_eq!(freed(&drops), 0);
        assert!(state.enter(1).end());
        assert_eq!(freed(&drops), 1);

        // The last one: its drop is the release.
        let (state, drops) = new_state(2);
        assert!(!state.enter(0).end());
        assert_eq!(freed(&drops), 0);
        {
            let stage = state.enter(1);
            assert_eq!(stage.tables().value, 7);
            assert_eq!(freed(&drops), 0);
        }
        assert_eq!(freed(&drops), 1);
    }

    #[test]
    fn a_stage_dropped_by_a_panic_counts() {
        let (state, drops) = new_state(1);
        let unwound = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let stage = state.enter(0);
            assert_eq!(stage.tables().value, 7);
            panic!("a panic inside the seed stage");
        }));
        assert!(unwound.is_err());
        assert_eq!(freed(&drops), 1);
    }

    #[test]
    fn stages_open_at_once_free_the_tables_once_after_the_last() {
        let (state, drops) = new_state(3);
        let a = state.enter(0);
        let b = state.enter(1);
        let c = state.enter(2);
        assert!(!c.end());
        assert_eq!(freed(&drops), 0);
        assert!(!a.end());
        assert_eq!(freed(&drops), 0);
        assert_eq!(b.tables().value, 7);
        assert!(b.end());
        assert_eq!(freed(&drops), 1);
        drop(state);
        assert_eq!(freed(&drops), 1);
    }

    #[cfg(any(not(target_arch = "wasm32"), feature = "wasm-threads"))]
    #[test]
    fn stages_on_several_threads_free_the_tables_once_after_the_last() {
        const THREADS: usize = 8;
        for round in 0..200 {
            let (state, drops) = new_state(THREADS);
            // Even rounds: every stage is open before the first one ends. Odd rounds: the
            // threads enter and end as they come.
            let all_open = std::sync::Barrier::new(if round % 2 == 0 { THREADS } else { 1 });
            let released = AtomicUsize::new(0);
            std::thread::scope(|scope| {
                for s_idx in 0..THREADS {
                    let (state, drops, all_open, released) = (&state, &drops, &all_open, &released);
                    scope.spawn(move || {
                        let stage = state.enter(s_idx);
                        all_open.wait();
                        assert_eq!(stage.tables().value, 7);
                        // This stage is open, so the tables are alive.
                        assert_eq!(freed(drops), 0);
                        if stage.end() {
                            // The release is the last end: every other stage has ended.
                            assert_eq!(freed(drops), 1);
                            released.fetch_add(1, AtomicOrdering::SeqCst);
                        }
                    });
                }
            });
            assert_eq!(released.load(AtomicOrdering::SeqCst), 1, "round {round}");
            assert_eq!(freed(&drops), 1, "round {round}");
        }
    }

    #[test]
    #[should_panic(
        expected = "the seed-stage tables of the batch were released before the seed stage of subject 5 started"
    )]
    fn enter_after_the_release_stops_the_run() {
        let (state, _drops) = new_state(1);
        assert!(state.enter(0).end());
        let _stage = state.enter(5);
    }

    #[test]
    fn a_count_that_is_too_low_keeps_the_tables_until_the_last_open_stage_ends() {
        // One stage announced, two open at once.
        let (state, drops) = new_state(1);
        let a = state.enter(0);
        let b = state.enter(1);
        // The announced count is used up: the batch's hold goes, the tables stay.
        assert!(a.end());
        assert_eq!(freed(&drops), 0);
        assert_eq!(b.tables().value, 7);
        // The count stays at zero; the last open stage frees the tables and reports no release.
        assert!(!b.end());
        assert_eq!(freed(&drops), 1);
    }

    #[test]
    fn a_count_that_is_too_high_keeps_the_tables_until_the_state_is_dropped() {
        let (state, drops) = new_state(3);
        assert!(!state.enter(0).end());
        assert!(!state.enter(1).end());
        assert_eq!(freed(&drops), 0);
        drop(state);
        assert_eq!(freed(&drops), 1);
    }

    #[test]
    fn a_batch_without_subjects_keeps_the_tables_until_the_state_is_dropped() {
        let (state, drops) = new_state(0);
        assert_eq!(freed(&drops), 0);
        drop(state);
        assert_eq!(freed(&drops), 1);
    }

    #[test]
    fn workers_get_no_diagonal_table_only_when_the_batch_is_units_only() {
        let cells = 1usize << 21;
        for (units_only, expected) in [(true, 0), (false, cells)] {
            assert_eq!(
                worker_table_cells(units_only, cells),
                expected,
                "units_only {units_only}"
            );
        }
    }

    #[test]
    fn drop_worker_table_leaves_no_allocation() {
        let mut table = vec![0u32; 1024];
        drop_worker_table(&mut table);
        assert!(table.is_empty());
        assert_eq!(table.capacity(), 0);
        // A table of zero cells (no reference frame loop in the batch) stays as it is.
        drop_worker_table(&mut table);
        assert_eq!(table.capacity(), 0);
    }

    #[test]
    fn the_release_adds_one_to_the_counter_of_its_state() {
        fn counter() -> &'static AtomicU64 {
            Box::leak(Box::new(AtomicU64::new(0)))
        }
        let count = |c: &AtomicU64| c.load(AtomicOrdering::SeqCst);

        // Once, at the release, also when further stages were open.
        let released = counter();
        let state = XSeedState::new(7u32, 2, Some(released));
        let a = state.enter(0);
        let b = state.enter(1);
        let c = state.enter(2);
        assert!(!a.end());
        assert_eq!(count(released), 0);
        assert!(b.end());
        assert_eq!(count(released), 1);
        assert!(!c.end());
        drop(state);
        assert_eq!(count(released), 1);
    }

    /// A lookup table without entries (every field of the type is public).
    fn empty_lookup() -> BlastAaLookupTable {
        BlastAaLookupTable {
            backbone: Vec::new(),
            overflow: Vec::new(),
            pv: Vec::new(),
            frame_bases: Vec::new(),
            num_contexts: 0,
            query_length: 0,
            word_length: 3,
            alphabet_size: 28,
            charsize: 5,
            mask: 0,
            longest_chain: 11,
            threshold: 13,
            row_max: Vec::new(),
        }
    }

    #[test]
    fn the_two_holders_of_a_batch_count_their_stages_apart() {
        // A batch without a unit path (the reference frame loop takes every subject): its
        // workers keep their diagonal table.
        let batch = XSeedBatch::new(empty_lookup(), None, 2);
        assert_eq!(batch.worker_table_cells(64), 64);

        // First subject: the units stage ends first, then the lookup stage.
        let lookup_stage = batch.lookup.enter(0);
        let units_stage = batch.units.enter(0);
        assert!(units_stage.tables().is_none());
        assert!(!units_stage.end());
        assert_eq!(lookup_stage.tables().longest_chain, 11);
        assert!(!lookup_stage.end());

        // Last subject: the end of its units stage releases the unit holder only; the lookup
        // table is still read through the stage that is open.
        let lookup_stage = batch.lookup.enter(1);
        let units_stage = batch.units.enter(1);
        assert!(units_stage.end());
        assert_eq!(lookup_stage.tables().threshold, 13);
        assert!(lookup_stage.end());
    }

    #[test]
    #[should_panic(
        expected = "the seed-stage tables of the batch were released before the seed stage of subject 1 started"
    )]
    fn a_unit_stage_after_the_release_of_the_unit_holder_stops_the_run() {
        let batch = XSeedBatch::new(empty_lookup(), None, 1);
        let lookup_stage = batch.lookup.enter(0);
        assert!(batch.units.enter(0).end());
        // The lookup holder still serves the open stage, and would serve a new one.
        assert_eq!(lookup_stage.tables().longest_chain, 11);
        let _units_stage = batch.units.enter(1);
    }

    #[test]
    fn the_holders_are_shared_between_the_threads_of_a_search() {
        fn shared<S: Send + Sync>() {}
        shared::<XSeedBatch>();
    }
}
