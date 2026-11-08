//! EXPERIMENT (LOSAT_X_BPLKPAHEAD): in the blastp block search
//! (`LOSAT_X_BPQBLOCK`, `x_query_blocks.rs`), the next block's query frames, lookup table and query
//! contexts are built during the block before it, by a task on the search's pool, instead of at the
//! start of the next block's own search.
//!
//! Without the switch each block's search (`run_resolved_in_pool`) encodes its queries
//! (`blastp_query_frames`) and builds its table and contexts (`blastp_compressed_lookup` for
//! blastp-fast, `blastp_aa_lookup` for exact blastp) before it scans the subjects, and its serial
//! stages (set-up, final traceback, output) and the tail of its redo windows leave the pool's other
//! threads idle. Exact blastp's table is built on the block's thread alone (11-15 ms per block of
//! 50,000 letters at 16 threads on the workstation, line R `REPORT_R3.md` §5.2), while the other 15
//! threads wait. In the human x yeast blastp-fast job at 16
//! threads (workstation, `LOSAT_X_BPQBLOCK=50000`, with `LOSAT_X_BPLKPFLAT` and
//! `LOSAT_X_BPSUBJRUN`) the lookup stage took 1.91 s of the 6.88 s, with 1.87 s of idle pool time in
//! the job (perf-mb-20261014 `logs/split2_report.md`). With the switch on, block i's search starts
//! block i + 1's build as soon as block i's own table is ready (at the block's start when it was
//! itself built ahead), as one task spawned in a `rayon::in_place_scope` around block i's search
//! (`XLookupAhead::run_block`); block i + 1's search then takes the values instead of computing them.
//!
//! Item SC2d: blastp-fast's flat table is built ahead in two parts, so that one table is alive at a
//! time: part 1 (`x_flat_ahead_appends`: the frames, the lookup query, NCBI's tables, the build's
//! appends pass with the contexts' Karlin blocks as its extra tasks, `tblastx/lookup/x_clut_par.rs`)
//! starts when block i's table is ready, beside block i's scan; part 2 (`x_flat_ahead_table`: the
//! stripes pass, which writes the table) runs on the vectors of block i's table, which block i's
//! thread gives to the second spare after block i's preliminary stage (`XAheadBlock::after_prelim`,
//! `x_make_room_flat`), when nothing of block i reads its table any more. The two meet without a
//! wait (`x_run_with_split_build_ahead`): the side that comes second starts part 2 (part 1's task
//! goes on into it, or block i's thread spawns it as one more task of the scope).
//!
//! Why the reports do not change: block i + 1's frames, table and contexts are functions of block i +
//! 1's query records and the search's arguments only. The task computes them by the calls the
//! block's own search makes, with the same inputs: `blastp_query_frames` of the block's records,
//! and `blastp_compressed_lookup` with the block's batches (`QueryInput::whole` of its records cut
//! by the search's batch size, as the block loop gives them to `run_resolved_in_pool`) and lengths,
//! and the same pool (the flat table and the Karlin blocks do not depend on the threads that
//! compute them, `tblastx/lookup/x_clut_par.rs` and `x_query_blocks.rs`); for exact blastp the call
//! of `blastp_aa_lookup` (`blastp_aa_lookup_ahead`: the same threshold, matrix, word length and
//! `check_ideal`), the in-place builder on the block search's spare vectors
//! (`tblastx/lookup/x_lut_direct.rs` `XAaSpares`: cleared and zero-filled to the table's sizes, as
//! the C `calloc` gives them) with neighbour lists memoised under the same key (each list is a
//! function of the word, the matrix, the threshold, the word length and the alphabet; a list
//! memoised by another thread is the list the recursion gives again). It is a leaf task: no rayon
//! inside, no lock but its spare's, held around a move. Nothing of block i reads them, and they
//! read nothing that block i writes: the only shared state is the block search's spare vectors,
//! whose values no build reads (`x_clut_par.rs`; the exact build writes every entry it keeps), and
//! the task uses its own (the second spare, `x_query_blocks::x_bp_lookup_spares`, or the exact
//! table's `XAaSpares::x_make_room`), and the exact build's neighbour lists, which only one build
//! uses at a time. The flat build's part 2 writes the vectors of block i's table only after block
//! i's thread has given them to the second spare, after block i's last read of the table (its
//! preliminary stage); no entry of them is read before the stripes pass writes it (`x_clut_par.rs`).
//! Its part 1's values are owned (`x_clut_par::XFlatAppends`): what part 1 computes does not depend
//! on the thread or the time at which part 2 runs. No floating-point operation changes order: each context's Karlin block is the
//! same call on the same letters. The search of block i + 1 checks that the values are those of its
//! own records (address and count of the slice).
//!
//! Where it acts: the block search, with `-subject` or `-db` (the block loop and the lookup region
//! are the same; `blast_engine/database.rs` only reads the subjects), on a pool of more than one
//! thread: for blastp-fast (`CompressedAaLookupTable`) where `LOSAT_X_BPLKPFLAT` acts on the flat
//! table of `LOSAT_X_BPCLUT` (on, `x_query_blocks::x_bp_lookup_spare`); for exact blastp
//! (`AaLookupTable`) where its table is built in place (`x_query_blocks::x_bp_aa_spares`), not
//! with `LOSAT_DIAGNOSTICS` (whose reference builder prints its phases, which would come earlier).
//! Not at one thread, for the search of all queries, or for the query chunks of a split batch (each builds its own table inside its block's search). The last block
//! starts no build.
//!
//! Errors and refusals: the task keeps nothing but a table. A build that fails or refuses (its
//! `Err` is dropped), that the reference builder would make (`TooLarge`), or whose vectors lack the
//! room (below) gives no values, and block i + 1 builds inline as without the switch, so its error
//! is raised by the same call at the same point of its search: after block i's rows, before any of
//! block i + 1's. In the block search such errors do not arise (`x_word_blocks` sends a query set
//! whose tables could refuse to the search of all queries before any block); the rule stands
//! regardless.
//!
//! Panics: the task catches a panic of its build (`catch_unwind`) and keeps it, and the block loop
//! resumes it (`resume_unwind`) when block i + 1's search would start, after block i's rows, at the
//! point where block i + 1's inline build would have raised it (nothing is written in between), so
//! the run ends with the same payload and exit status (101). The panic hook prints its message when
//! the panic happens, on the pool thread that ran the task, so that line can come earlier in standard
//! error (during block i) and name another thread. If block i itself fails, its error ends the run and
//! the kept panic is dropped (its message was printed). With `panic = "abort"` (the release profile)
//! a panic cannot be kept: it ends the process when it happens, as a panic of any pool task does
//! today, possibly before block i's rows are written.
//!
//! No deadlock: the task waits for nothing but its own nested parallel passes (the one-latch passes
//! of the flat table and its Karlin tasks; the exact build has none), which wait for nothing else,
//! and block i's search waits for nothing of the task until the scope ends, after block i's search
//! has returned. The flat build's two parts meet in one slot (`XHandoff`) whose lock is held only
//! around a move: neither block i's thread nor part 1's task waits for the other side; part 2 is
//! a task of the same scope, like part 1. A thread that waits inside block i's nested parallel work may steal and run the
//! task (or one of its passes) to the end, and a thread that waits inside the task's passes may run
//! one of block i's jobs to the end; neither waits for the frame it interrupted, so both waits end
//! (memory `rayon-waiting-helper-deadlock`: the hazard there was a spawned loop that waited for
//! progress made elsewhere). The scope's owner, the block's thread, runs the task itself if no
//! other thread took it. The locks the task takes (the spare's, the result's slot) are held only
//! around a move.
//!
//! Memory: for exact blastp two tables are alive at once (block i's, which its search scans, and
//! block i + 1's), so the peak holds one more table's backbone (0.44 MiB), overflow (about 7-8 MiB
//! at `LOSAT_X_BPQBLOCK=50000`) and PV array (3.5 KiB), plus block i + 1's query frames and
//! contexts. For blastp-fast (item SC2d) block i + 1's table is written on block i's vectors after
//! block i's preliminary stage: beside block i's table and scan there are only block i + 1's appends
//! (as many as its table's hits) and the appends pass's scratch, and after the preliminary stage one
//! table and the appends. Building the whole table during block i's scan (before SC2d) held two
//! tables and the appends at the end of the build. All of the table's vectors are allocated by the block
//! search's thread, never by the pool's task (a pool thread's allocations stay in its own glibc
//! arena, which keeps them resident after they are freed): that thread makes the second spare hold
//! vectors with the room for the largest table so far plus a half (`x_room_of`): for blastp-fast
//! after block i's preliminary stage, block i's table's vectors (`x_make_room_flat`,
//! `XFlatSpare::x_make_room`: `hits` grown in place where short, to twice the largest table's hits;
//! a PV array where the spare has none); for exact blastp when block i's table is ready
//! (`XAaSpares::x_make_room`, the same rule for the overflow entries; each block's table is given
//! back at the end of its search). The task's build takes them all (`XFlatArgs::x_ahead`,
//! `XAaSpares::x_build_ahead`). A build whose `hits` (overflow) would not fit gives no values and the
//! block builds its table itself, on its own thread (for blastp-fast on the vectors the second
//! spare got, moved to the block spare: `XFlatSpare::x_give_to`). Untouched capacity is not resident. The exact build's neighbour lists (about
//! 1.5-1.9 MiB) are the block search's thread's, moved to `XAaSpares` once and used by every build
//! ahead, so no pool thread keeps a set of its own.
//!
//! `LOSAT_X_BPLKPAHEAD` is a switch of the strict set (on by default since 2026-10-17): on unless it is `0`.
//!
//! With `LOSAT_X_STATS`, `main` prints `[X_STATS] BPLKPAHEAD blocks= prebuilt= dropped=` at exit
//! when a counter is not zero (`blocks`: blocks searched where the switch acts; `prebuilt`: blocks
//! that took values built ahead; `dropped`: builds ahead that gave no values). The unit test
//! `values_built_ahead_are_the_values_the_block_builds_itself` builds a block's values ahead and
//! inline, for blastp-fast (its two parts in either order) and exact blastp, and compares the query
//! frames, the table and the contexts (the Karlin blocks bit for bit).
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:257-276
//! ```c++
//!         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//!             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//!             CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
//!             ...
//!                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
//!                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
//!                 results = lcl_blast.Run();
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:206-215
//! ```c++
//!     // 5. Create the lookup table
//!     if ( !retval->m_QuerySplitter->IsQuerySplit() ) {
//!         LookupTableWrap* lut =
//!             CSetupFactory::CreateLookupTable(query_data, opts_memento.get(),
//!                                              sbp, lookup_segments_wrap,
//!                                              retval->m_InternalData->m_RpsData,
//!                                              seqsrc,
//!                                              num_threads);
//! ```
//! NCBI makes each batch's lookup table when that batch's search is set up, after the previous
//! batch's report. The table is a function of the batch's queries and the options only, so building
//! it earlier, beside the previous block's search, gives the same table. Scheduling only.

use std::cell::Cell;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use std::sync::{Mutex, PoisonError};

use super::{blastp_aa_lookup_ahead, blastp_query_frames, BlastpRuntimeLookup};
use crate::algorithm::blastp::args::{BlastpLookupTableType, ResolvedBlastpArgs};
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use crate::algorithm::blastp::x_query_blocks::XLookupFor;
use crate::algorithm::blastp::x_query_blocks::{
    x_bp_aa_block_spare, x_bp_aa_spares, x_bp_lkreuse, x_bp_lookup_spares,
};
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::{
    x_bp_flat_ahead_appends, x_bp_flat_ahead_table, XFlatAppends,
};
use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::{
    x_bpclut_on, XFlatCompressedLookup,
};
use crate::algorithm::tblastx::lookup::{QueryContext, XAaSpares};
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::blastinput::fasta_reader::FastaRecord;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use crate::blastinput::seq_range;
use crate::utils::threading::SearchPool;

// No NCBI counterpart: reads the LOSAT_X_BPLKPAHEAD switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPLKPAHEAD): a strict-set switch, on unless it is `0`.
pub(crate) fn x_bp_lkp_ahead() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPLKPAHEAD"))
}

// No NCBI counterpart: counters of the `[X_STATS] BPLKPAHEAD` line; they change no value NCBI
// computes.
static X_BLOCKS: AtomicU64 = AtomicU64::new(0);
static X_PREBUILT: AtomicU64 = AtomicU64::new(0);
static X_DROPPED: AtomicU64 = AtomicU64::new(0);

// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPLKPAHEAD): with `LOSAT_X_STATS`, the counters' line at exit (printed by
/// `main`) when a counter is not zero: `[X_STATS] BPLKPAHEAD blocks= prebuilt= dropped=`.
pub(crate) fn x_bp_lkp_ahead_print_summary() {
    if std::env::var_os("LOSAT_X_STATS").is_some() {
        let load = |counter: &AtomicU64| counter.load(Ordering::Relaxed);
        let (blocks, prebuilt, dropped) = (load(&X_BLOCKS), load(&X_PREBUILT), load(&X_DROPPED));
        if blocks != 0 || prebuilt != 0 || dropped != 0 {
            eprintln!("[X_STATS] BPLKPAHEAD blocks={blocks} prebuilt={prebuilt} dropped={dropped}");
        }
    }
}

// No NCBI counterpart: identifies a block's records; it changes no value NCBI computes.
/// The address and count of a block's query records (the slice the block loop gives its search).
fn x_key(query_records: &[FastaRecord]) -> (usize, usize) {
    (query_records.as_ptr() as usize, query_records.len())
}

/// EXPERIMENT (LOSAT_X_BPLKPAHEAD): a block's query frames, lookup table (blastp-fast's flat table
/// or exact blastp's table) and query contexts, built ahead of its search (`x_build_ahead`).
pub(super) struct XAheadLookup {
    /// The block's query records (`x_key`).
    key: (usize, usize),
    query_frames: Vec<Vec<QueryFrame>>,
    table: BlastpRuntimeLookup,
    contexts: Vec<QueryContext>,
}

impl XAheadLookup {
    // No NCBI counterpart: hands over values built ahead; it changes no value NCBI computes.
    /// The query frames (taken; the search keeps them as its own).
    pub(super) fn take_query_frames(&mut self) -> Vec<Vec<QueryFrame>> {
        std::mem::take(&mut self.query_frames)
    }
}

/// EXPERIMENT (LOSAT_X_BPLKPAHEAD): what a block's search (`run_resolved_in_pool`) gets from the block
/// loop: the values built ahead for it, if any, and the start of the next block's build, called
/// once on the block's thread when the block's lookup table is ready (`lookup_ready`). Every other
/// search gets `none()`.
pub(super) struct XAheadBlock<'a> {
    prebuilt: Option<XAheadLookup>,
    start_next: Option<&'a dyn Fn(&BlastpRuntimeLookup)>,
    /// Item SC2d: called once on the block's thread after the block's preliminary stage with the
    /// block's flat table, whose vectors the second part of the next block's build takes
    /// (`after_prelim`).
    after_prelim_flat: Option<&'a dyn Fn(XFlatCompressedLookup)>,
    /// EXPERIMENT (LOSAT_X_BPLKDROP): whether a later block's build may take the vectors of this
    /// block's table from the block search's spares (`x_table_taken_later`).
    table_taken_later: bool,
}

impl XAheadBlock<'_> {
    // No NCBI counterpart: chooses LOSAT's path; it changes no value NCBI computes.
    /// Nothing built ahead and nothing to start: the search as without the switch.
    pub(super) fn none() -> Self {
        Self {
            prebuilt: None,
            start_next: None,
            after_prelim_flat: None,
            table_taken_later: true,
        }
    }

    // No NCBI counterpart: hands over values built ahead (module header); it changes no value NCBI
    // computes.
    /// The values built ahead for `query_records` (the same slice), counted as taken; values built
    /// for another slice (which the block loop does not make) are dropped, counted as such, and the
    /// search builds its own.
    pub(super) fn take_prebuilt(&mut self, query_records: &[FastaRecord]) -> Option<XAheadLookup> {
        let prebuilt = self.prebuilt.take()?;
        if prebuilt.key != x_key(query_records) {
            X_DROPPED.fetch_add(1, Ordering::Relaxed);
            return None;
        }
        X_PREBUILT.fetch_add(1, Ordering::Relaxed);
        Some(prebuilt)
    }

    // No NCBI counterpart: starts the next block's build (module header); it changes no value NCBI
    // computes.
    /// The block's lookup table is `lookup`: starts the next block's build, if the block loop gave
    /// one to start.
    pub(super) fn lookup_ready(&self, lookup: &BlastpRuntimeLookup) {
        if let Some(start_next) = self.start_next {
            start_next(lookup);
        }
    }

    // No NCBI counterpart: hands the vectors of a table that nothing reads any more to the next
    // block's build (module header, item SC2d); it changes no value NCBI computes.
    /// The block's preliminary stage has ended and nothing reads its lookup table `lookup` any more:
    /// where the block loop gave the second part of the next block's flat table to start, a flat
    /// table's vectors go to it (`None`: the table was taken); otherwise the table is given back.
    pub(super) fn after_prelim(&self, lookup: BlastpRuntimeLookup) -> Option<BlastpRuntimeLookup> {
        match (lookup, self.after_prelim_flat) {
            (BlastpRuntimeLookup::CompressedFlat(table), Some(after_prelim_flat)) => {
                after_prelim_flat(table);
                None
            }
            (lookup, _) => Some(lookup),
        }
    }

    // No NCBI counterpart: memory only (LOSAT_X_BPLKDROP); it changes no value NCBI computes.
    /// Whether a later block's build may take the vectors of this block's lookup table, which the
    /// search gives to the block search's spares at its end (`x_query_blocks::x_bp_lookup_spare`,
    /// `x_bp_aa_spares`): the next block's own build where the switch does not act, the build
    /// ahead of the block after the next where it does (`x_make_room`, called when the next block's
    /// table is ready). True for `none()` (a search that is not a block of the block loop; outside
    /// the block search there are no spares).
    pub(super) fn x_table_taken_later(&self) -> bool {
        self.table_taken_later
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:206-215 (module header)
// ```c++
//     // 5. Create the lookup table
//     if ( !retval->m_QuerySplitter->IsQuerySplit() ) {
//         LookupTableWrap* lut =
//             CSetupFactory::CreateLookupTable(query_data, opts_memento.get(),
// ```
/// EXPERIMENT (LOSAT_X_BPLKPAHEAD): the lookup table and contexts of the lookup region from values
/// built ahead (`prebuilt`, whose query frames the search took with `take_query_frames`). They are
/// the values the block's own build gives (module header; the unit test
/// `values_built_ahead_are_the_values_the_block_builds_itself` compares them).
pub(super) fn x_use_prebuilt(prebuilt: XAheadLookup) -> (BlastpRuntimeLookup, Vec<QueryContext>) {
    let XAheadLookup {
        table, contexts, ..
    } = prebuilt;
    (table, contexts)
}

/// What the build ahead of a block gave (`V`: its values).
#[cfg_attr(
    not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )),
    expect(
        dead_code,
        reason = "Values, Dropped and Panicked are constructed only by x_run_with_build_ahead (used only with `parallel` threads)"
    )
)]
enum XBuilt<V> {
    Values(V),
    /// No values: the block builds inline (module header).
    Dropped,
    /// The build panicked: resumed where the block's search would start.
    Panicked(Box<dyn std::any::Any + Send>),
}

// No NCBI counterpart: hands over the outcome of a build ahead (module header); it changes no value
// NCBI computes.
/// The values that the outcome of the build ahead of a block gives that block: its values; none for
/// a build that gave none (counted), so that the block builds inline; a kept panic is resumed here,
/// before the block's search (module header).
fn x_take_built<V>(built: Option<XBuilt<V>>) -> Option<V> {
    match built? {
        XBuilt::Values(values) => Some(values),
        XBuilt::Dropped => {
            X_DROPPED.fetch_add(1, Ordering::Relaxed);
            None
        }
        XBuilt::Panicked(payload) => std::panic::resume_unwind(payload),
    }
}

// No NCBI counterpart: hands over the outcome of a build ahead (module header); it changes no value
// NCBI computes.
/// The outcome of a build ahead from its result under `catch_unwind`.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn x_built_of<V>(outcome: std::thread::Result<Option<V>>) -> XBuilt<V> {
    match outcome {
        Ok(Some(values)) => XBuilt::Values(values),
        Ok(None) => XBuilt::Dropped,
        Err(payload) => XBuilt::Panicked(payload),
    }
}

// No NCBI counterpart: scheduling only (module header: why the scope cannot deadlock); it changes no
// value NCBI computes.
/// The scheduling of `XLookupAhead::run_block`: runs `run` on this thread inside a
/// `rayon::in_place_scope` of the current pool, giving it `spawn`, which `run` calls on this thread
/// when the next block's build may start: the first call spawns `build` as one task of the scope
/// (later calls do nothing). Returns `run`'s result and, once every task of the scope has ended, the
/// build's outcome: its values, none, or its panic (caught; `None` if it was never spawned).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn x_run_with_build_ahead<V: Send, R>(
    build: impl FnOnce() -> Option<V> + Send,
    run: impl FnOnce(&dyn Fn()) -> R,
) -> (R, Option<XBuilt<V>>) {
    let built: Mutex<Option<XBuilt<V>>> = Mutex::new(None);
    let result = rayon::in_place_scope(|scope| {
        let build = Cell::new(Some(build));
        let built = &built;
        let spawn = || {
            if let Some(build) = build.take() {
                scope.spawn(move |_| {
                    let outcome = x_built_of(std::panic::catch_unwind(
                        std::panic::AssertUnwindSafe(build),
                    ));
                    *built.lock().unwrap_or_else(PoisonError::into_inner) = Some(outcome);
                });
            }
        };
        run(&spawn)
    });
    (
        result,
        built.into_inner().unwrap_or_else(PoisonError::into_inner),
    )
}

/// Item SC2d: where the two parts of a split build ahead meet (`x_run_with_split_build_ahead`).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
enum XHandoff<P> {
    /// Neither part 1's values nor the vectors yet.
    Empty,
    /// Part 1's values, before the vectors.
    Part(P),
    /// The vectors, before part 1's values.
    Vectors,
    /// Part 2 has started, or part 1 gave no values: nothing more is handed over.
    Closed,
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
impl<P> XHandoff<P> {
    // No NCBI counterpart: scheduling only; it changes no value NCBI computes.
    /// Part 1 gave `part`: `Some(part)` when the vectors came first (part 2 runs now, on the caller's
    /// task), else `None` (kept until the vectors come). The lock is held only around the move.
    fn part(slot: &Mutex<Self>, part: P) -> Option<P> {
        let mut slot = slot.lock().unwrap_or_else(PoisonError::into_inner);
        match std::mem::replace(&mut *slot, Self::Closed) {
            Self::Vectors => Some(part),
            Self::Empty => {
                *slot = Self::Part(part);
                None
            }
            Self::Part(_) | Self::Closed => unreachable!("part 1 runs once"),
        }
    }

    // No NCBI counterpart: scheduling only; it changes no value NCBI computes.
    /// The vectors are ready: part 1's values when they came first (part 2 starts now), else `None`
    /// (part 2 starts when part 1 ends, or never: part 1 gave no values, or this is a second call).
    fn vectors(slot: &Mutex<Self>) -> Option<P> {
        let mut slot = slot.lock().unwrap_or_else(PoisonError::into_inner);
        match std::mem::replace(&mut *slot, Self::Closed) {
            Self::Part(part) => Some(part),
            Self::Empty => {
                *slot = Self::Vectors;
                None
            }
            state @ (Self::Vectors | Self::Closed) => {
                *slot = state;
                None
            }
        }
    }

    // No NCBI counterpart: scheduling only; it changes no value NCBI computes.
    /// Part 1 gave no values (or panicked): the vectors start nothing.
    fn close(slot: &Mutex<Self>) {
        *slot.lock().unwrap_or_else(PoisonError::into_inner) = Self::Closed;
    }
}

// No NCBI counterpart: scheduling only (module header: why nothing waits); it changes no value NCBI
// computes.
/// Item SC2d: the scheduling of a build ahead in two parts (`XLookupAhead::run_block` for the flat
/// table): runs `run` on this thread inside a `rayon::in_place_scope` of the current pool, giving it
/// `start` and `vectors`, which `run` calls on this thread. The first call of `start` spawns `first`
/// (part 1) as one task of the scope; `vectors` says that what part 2 takes is ready. Whichever of
/// the two comes second starts `second` (part 2) with part 1's values: the task of part 1 goes on
/// into it, or `vectors` spawns it as one more task of the scope. Neither side waits for the other.
/// Returns `run`'s result and, once every task of the scope has ended, the build's outcome: part 2's
/// values, none (part 1 or part 2 gave none, or part 1's values never got the vectors), or the
/// panic of either part (caught); `None` if `start` was never called.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn x_run_with_split_build_ahead<P: Send, V: Send, R>(
    first: impl FnOnce() -> Option<P> + Send,
    second: impl Fn(P) -> Option<V> + Sync,
    run: impl FnOnce(&dyn Fn(), &dyn Fn()) -> R,
) -> (R, Option<XBuilt<V>>) {
    let built: Mutex<Option<XBuilt<V>>> = Mutex::new(None);
    let handoff: Mutex<XHandoff<P>> = Mutex::new(XHandoff::Empty);
    let result = rayon::in_place_scope(|scope| {
        let first = Cell::new(Some(first));
        let (built, handoff, second) = (&built, &handoff, &second);
        let keep = move |outcome: XBuilt<V>| {
            *built.lock().unwrap_or_else(PoisonError::into_inner) = Some(outcome);
        };
        // Part 2, once.
        let finish = move |part: P| {
            keep(x_built_of(std::panic::catch_unwind(
                std::panic::AssertUnwindSafe(|| second(part)),
            )));
        };
        let start = || {
            if let Some(first) = first.take() {
                scope.spawn(move |_| {
                    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(first)) {
                        Ok(Some(part)) => {
                            if let Some(part) = XHandoff::part(handoff, part) {
                                finish(part);
                            }
                        }
                        outcome => {
                            XHandoff::close(handoff);
                            keep(x_built_of(outcome.map(|_| None)));
                        }
                    }
                });
            }
        };
        let vectors = || {
            if let Some(part) = XHandoff::vectors(handoff) {
                scope.spawn(move |_| finish(part));
            }
        };
        run(&start, &vectors)
    });
    let built = match handoff.into_inner().unwrap_or_else(PoisonError::into_inner) {
        // Part 1's values that never got the vectors.
        XHandoff::Part(_) => Some(XBuilt::Dropped),
        _ => built.into_inner().unwrap_or_else(PoisonError::into_inner),
    };
    (result, built)
}

// No NCBI counterpart: the floor of the room for the next table's hits (module header); untouched
// capacity is not resident. It changes no value NCBI computes.
/// The fewest hits (exact blastp: overflow entries) the second spare has room for: 2^18 (1 MiB).
const X_MIN_HITS_ROOM: usize = 1 << 18;

// No NCBI counterpart: the room for the next table's hits (module header); it changes no value NCBI
// computes.
/// The room for the next block's table, from a table of `hits` hits (exact blastp: overflow entries)
/// and the most of a table so far (`largest_hits`, raised to `hits`): the room asked (the largest
/// plus a half, at least `X_MIN_HITS_ROOM`) and the capacity a vector short of it is given (twice
/// the largest, at least the room).
fn x_room_of(largest_hits: &Cell<usize>, hits: usize) -> (usize, usize) {
    let largest = largest_hits.get().max(hits);
    largest_hits.set(largest);
    let room = (largest + largest / 2).max(X_MIN_HITS_ROOM);
    (room, (2 * largest).max(room))
}

// No NCBI counterpart: which thread allocates the next table's vectors (module header); it changes
// no value NCBI computes.
/// On the block's thread, when its exact table `lookup` is ready: makes `XAaSpares` hold vectors with
/// the room for the next block's table (`x_room_of`), allocating here what is short. `false` where
/// the block's table is not an exact one with its spares (no build ahead then; a flat table's room is
/// made after the block's preliminary stage, `x_make_room_flat`).
#[cfg_attr(
    not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )),
    expect(
        dead_code,
        reason = "called only by run_block's build ahead (used only with `parallel` threads)"
    )
)]
fn x_make_room(
    pool: &SearchPool<'_>,
    lookup: &BlastpRuntimeLookup,
    largest_hits: &Cell<usize>,
) -> bool {
    match lookup {
        // Exact blastp: the backbone and the pv array have the sizes of the block's own table (fixed
        // by the word length and the alphabet).
        BlastpRuntimeLookup::Aa(table) => {
            let Some(spares) = x_bp_aa_spares(pool) else {
                return false;
            };
            let (room, grown) = x_room_of(largest_hits, table.overflow.len());
            spares.x_make_room(table.backbone.len(), table.pv.len(), room, grown);
            true
        }
        // EXPERIMENT (LOSAT_X_BPLK16): the same for a table of 16-bit hits; the lengths of its
        // vectors count their padding, which the next table of that width has as well.
        BlastpRuntimeLookup::Aa16(table) => {
            let Some(spares) = x_bp_aa_spares(pool) else {
                return false;
            };
            let (room, grown) = x_room_of(largest_hits, table.overflow.len());
            spares.x_make_room16(table.backbone.len(), table.pv.len(), room, grown);
            true
        }
        BlastpRuntimeLookup::CompressedFlat(_) | BlastpRuntimeLookup::Compressed(_) => false,
    }
}

// No NCBI counterpart: which vectors the next table takes (module header, item SC2d); it changes no
// value NCBI computes.
/// On the block's thread, after its preliminary stage: the second spare takes the `starts` and
/// `hits` of the block's flat table `table` (which nothing reads any more) and gets the room for the
/// next block's table as `x_make_room` makes it (the hits grown in place where short). `false` where
/// the block search has no spares (the table is dropped).
#[cfg_attr(
    not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )),
    expect(
        dead_code,
        reason = "called only by run_block's build ahead (used only with `parallel` threads)"
    )
)]
fn x_make_room_flat(
    pool: &SearchPool<'_>,
    table: XFlatCompressedLookup,
    largest_hits: &Cell<usize>,
) -> bool {
    let Some((block, ahead)) = x_bp_lookup_spares(pool) else {
        return false;
    };
    let backbone_size = table.starts.len().saturating_sub(1);
    let (room, grown) = x_room_of(largest_hits, table.hits.len());
    ahead.give(table);
    ahead.x_make_room(block, backbone_size, room, grown);
    true
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:206-215 (module header)
// ```c++
//         LookupTableWrap* lut =
//             CSetupFactory::CreateLookupTable(query_data, opts_memento.get(),
//                                              sbp, lookup_segments_wrap,
// ```
/// EXPERIMENT (LOSAT_X_BPLKPAHEAD): the query frames, table and contexts of a block of
/// `query_records` for exact blastp, as its search builds them (module header: the same calls with
/// the same inputs), its table on the spare vectors of `XAaSpares`. `None` where the build fails,
/// refuses or makes no table on the spare's vectors, and for blastp-fast, whose flat table is built
/// ahead in two parts (`x_flat_ahead_appends`, `x_flat_ahead_table`; item SC2d).
#[cfg_attr(
    not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )),
    expect(
        dead_code,
        reason = "called only by run_block's build ahead (used only with `parallel` threads)"
    )
)]
fn x_build_ahead(
    args: &ResolvedBlastpArgs,
    query_records: &[FastaRecord],
    pool: &SearchPool<'_>,
) -> Option<XAheadLookup> {
    if args.lookup_table_type != BlastpLookupTableType::AaLookupTable {
        return None;
    }
    let query_frames = blastp_query_frames(args, query_records);
    // `blastp_aa_lookup`'s call on the block search's spare vectors.
    let (table, contexts) = blastp_aa_lookup_ahead(args, &query_frames, x_bp_aa_spares(pool)?)?;
    Some(XAheadLookup {
        key: x_key(query_records),
        query_frames,
        table,
        contexts,
    })
}

/// Item SC2d: what part 1 of the build ahead of a block's flat table gives part 2
/// (`x_flat_ahead_appends`).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
struct XFlatAheadPart {
    key: (usize, usize),
    query_frames: Vec<Vec<QueryFrame>>,
    contexts: Vec<QueryContext>,
    /// The lookup query and its segments, which part 2 reads only for a table's capacity test.
    query: Vec<u8>,
    locations: Vec<(i32, i32)>,
    appends: XFlatAppends,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:206-215 (module header)
// ```c++
//         LookupTableWrap* lut =
//             CSetupFactory::CreateLookupTable(query_data, opts_memento.get(),
//                                              sbp, lookup_segments_wrap,
// ```
/// Item SC2d (LOSAT_X_BPLKPAHEAD): part 1 of the build ahead of a block of `query_records` for
/// blastp-fast (BLOSUM62's flat table), on a pool task once the block before it has its table: the
/// calls of `blastp_compressed_lookup` with the same inputs (`batch_size` is the search's query
/// batch size; `XLookupFor::Ahead`) up to the flat build's appends pass (`x_bp_flat_ahead_appends`,
/// with the contexts' Karlin blocks as its extra tasks), then the contexts' Karlin blocks set
/// (`XLookupFlat::finish`). `None` where the search's build would fail or the build ends at its
/// appends pass (module header).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn x_flat_ahead_appends(
    args: &ResolvedBlastpArgs,
    query_records: &[FastaRecord],
    batch_size: u32,
    pool: &SearchPool<'_>,
) -> Option<XFlatAheadPart> {
    if args.lookup_table_type != BlastpLookupTableType::CompressedAaLookupTable
        || args.scoring.matrix != crate::config::ScoringMatrix::Blosum62
    {
        return None;
    }
    let query_frames = blastp_query_frames(args, query_records);
    let query_lengths: Vec<usize> = query_records
        .iter()
        .map(|record| record.seq().len())
        .collect();
    let query_batches = seq_range::QueryInput::whole(query_records).batches(batch_size);
    let (x_flat, query, locations, mut contexts) =
        crate::algorithm::blastp::x_query_blocks::XLookupFlat::lookup_query(
            &query_frames,
            args.word_size,
            pool,
            XLookupFor::Ahead,
        );
    let ncbi_tables = super::ncbi_compressed_lookup_tables(
        &query_batches,
        &query_lengths,
        &contexts,
        &locations,
        false,
    )
    .ok()?;
    let appends = x_bp_flat_ahead_appends(
        args.word_size,
        crate::core::blast_util::ncbi_int4_from_double(100.0 * args.threshold),
        &query,
        &locations,
        &ncbi_tables,
        pool,
        x_flat.build_args(crate::algorithm::blastp::x_query_blocks::x_block_appends(
            pool,
            query_records,
        )),
    )?;
    x_flat.finish(&mut contexts);
    Some(XFlatAheadPart {
        key: x_key(query_records),
        query_frames,
        contexts,
        query,
        locations,
        appends,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:206-215 (module header)
// ```c++
//         LookupTableWrap* lut =
//             CSetupFactory::CreateLookupTable(query_data, opts_memento.get(),
// ```
/// Item SC2d (LOSAT_X_BPLKPAHEAD): part 2 of the build ahead of a block's flat table, once part 1
/// gave `part` and the block before it has given its table's vectors to the second spare after its
/// preliminary stage (`x_make_room_flat`): the stripes pass on those vectors
/// (`x_bp_flat_ahead_table`). `None` where they lack the room (the second spare keeps them) or a
/// table refuses.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn x_flat_ahead_table(part: XFlatAheadPart, pool: &SearchPool<'_>) -> Option<XAheadLookup> {
    let (_, spare) = x_bp_lookup_spares(pool)?;
    let XFlatAheadPart {
        key,
        query_frames,
        contexts,
        query,
        locations,
        appends,
    } = part;
    let table = x_bp_flat_ahead_table(appends, &query, &locations, pool, spare)?;
    crate::utils::xreach::x_reach!("LOSAT_X_BPCLUT", true);
    Some(XAheadLookup {
        key,
        query_frames,
        table: BlastpRuntimeLookup::CompressedFlat(table),
        contexts,
    })
}

/// EXPERIMENT (LOSAT_X_BPLKPAHEAD): the block loop's side of the switch (module header): whether it
/// acts for this block search, and what the build ahead of the next block gave.
pub(super) struct XLookupAhead<'a, 'p> {
    args: &'a ResolvedBlastpArgs,
    pool: &'a SearchPool<'p>,
    #[cfg_attr(
        not(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        )),
        expect(
            dead_code,
            reason = "read only by run_block's build ahead (used only with `parallel` threads)"
        )
    )]
    batch_size: u32,
    acts: bool,
    next: Option<XBuilt<XAheadLookup>>,
    /// The most hits of a table so far (`x_make_room`).
    #[cfg_attr(
        not(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        )),
        expect(
            dead_code,
            reason = "read only by run_block's build ahead (used only with `parallel` threads)"
        )
    )]
    largest_hits: Cell<usize>,
}

impl<'a, 'p> XLookupAhead<'a, 'p> {
    // No NCBI counterpart: chooses LOSAT's path; it does not change any value NCBI computes.
    /// The block loop's state for the block search of `args` on `pool` (the pool with the block
    /// search's subjects, `x_query_blocks::XBlockSubjects`); `batch_size` is the search's query
    /// batch size. It acts where the module header says.
    pub(super) fn new(
        args: &'a ResolvedBlastpArgs,
        pool: &'a SearchPool<'p>,
        batch_size: u32,
    ) -> Self {
        let acts = cfg!(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        )) && x_bp_lkp_ahead()
            && match args.lookup_table_type {
                // The flat table's (BLOSUM62's; another matrix's flat builders take no spare
                // vectors: item SC2b).
                BlastpLookupTableType::CompressedAaLookupTable => {
                    x_bpclut_on()
                        && args.scoring.matrix == crate::config::ScoringMatrix::Blosum62
                        && x_bp_lookup_spares(pool).is_some()
                }
                BlastpLookupTableType::AaLookupTable => x_bp_aa_spares(pool).is_some(),
            };
        if x_bp_lkp_ahead() {
            crate::utils::xreach::x_reach!("LOSAT_X_BPLKPAHEAD", acts);
        }
        // EXPERIMENT (LOSAT_X_BPLKREUSE): where exact blastp's tables are not built ahead, each
        // block's build takes the vectors of the block before (`x_bp_aa_block_spare`).
        if x_bp_lkreuse()
            && matches!(args.lookup_table_type, BlastpLookupTableType::AaLookupTable)
            && XAaSpares::x_applies()
        {
            crate::utils::xreach::x_reach!(
                "LOSAT_X_BPLKREUSE",
                x_bp_aa_block_spare(pool).is_some()
            );
        }
        Self {
            args,
            pool,
            batch_size,
            acts,
            next: None,
            largest_hits: Cell::new(0),
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:257-276 (module header)
    // ```c++
    //         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
    //             ...
    //                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
    //                 results = lcl_blast.Run();
    // ```
    /// Searches a block with `run` (its `run_resolved_in_pool`), giving it the values built ahead for
    /// it (its search checks that they are its records', `XAheadBlock::take_prebuilt`); with `next`
    /// (the next block's records) and where the switch acts,
    /// the next block's build is spawned in a `rayon::in_place_scope` around the search when the
    /// block's table is ready, and its outcome kept for the next call. A kept panic of the build for
    /// this block is resumed here, before its search (module header). `later_blocks` is the number
    /// of blocks after this one (`XAheadBlock::x_table_taken_later`).
    pub(super) fn run_block<R>(
        &mut self,
        next: Option<&[FastaRecord]>,
        later_blocks: usize,
        run: impl FnOnce(XAheadBlock<'_>) -> anyhow::Result<R>,
    ) -> anyhow::Result<R> {
        let prebuilt = x_take_built(self.next.take());
        if !self.acts {
            return run(XAheadBlock {
                prebuilt,
                start_next: None,
                after_prelim_flat: None,
                table_taken_later: later_blocks >= 1,
            });
        }
        X_BLOCKS.fetch_add(1, Ordering::Relaxed);
        // Exact blastp: where this block's thread may build its table itself (nothing built ahead
        // for it) and at the last block (after which that thread may search again), the neighbour
        // lists kept for the builds ahead go back to it; at the last block (no build ahead follows)
        // the spare vectors are dropped.
        if prebuilt.is_none() || next.is_none() {
            if let (BlastpLookupTableType::AaLookupTable, Some(spares)) =
                (self.args.lookup_table_type, x_bp_aa_spares(self.pool))
            {
                spares.x_home(next.is_none());
            }
        }
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        if let Some(next) = next {
            let (args, pool, batch_size, largest_hits) =
                (self.args, self.pool, self.batch_size, &self.largest_hits);
            // Item SC2d: blastp-fast's flat table in two parts (module header).
            if args.lookup_table_type == BlastpLookupTableType::CompressedAaLookupTable {
                let (result, built) = x_run_with_split_build_ahead(
                    move || x_flat_ahead_appends(args, next, batch_size, pool),
                    |part| x_flat_ahead_table(part, pool),
                    |start, vectors| {
                        // On this thread, when the block's flat table is ready: part 1 (the appends
                        // pass) is spawned.
                        let start_next = |lookup: &BlastpRuntimeLookup| {
                            if matches!(lookup, BlastpRuntimeLookup::CompressedFlat(_)) {
                                start();
                            }
                        };
                        // After the block's preliminary stage: its table's vectors go to part 2
                        // (the stripes pass).
                        let after_prelim_flat = |table: XFlatCompressedLookup| {
                            if x_make_room_flat(pool, table, largest_hits) {
                                vectors();
                            }
                        };
                        run(XAheadBlock {
                            prebuilt,
                            start_next: Some(&start_next),
                            after_prelim_flat: Some(&after_prelim_flat),
                            table_taken_later: later_blocks >= 2,
                        })
                    },
                );
                // No table ahead: the next block builds its own, on the vectors the second spare
                // got (the block spare's).
                if !matches!(built, Some(XBuilt::Values(_))) {
                    if let Some((block, ahead)) = x_bp_lookup_spares(pool) {
                        ahead.x_give_to(block);
                    }
                }
                self.next = built;
                return result;
            }
            let (result, built) = x_run_with_build_ahead(
                move || x_build_ahead(args, next, pool),
                |spawn| {
                    // On this thread, when the block's table is ready: the next table's vectors are
                    // made ready here, then its build is spawned.
                    let start_next = |lookup: &BlastpRuntimeLookup| {
                        if x_make_room(pool, lookup, largest_hits) {
                            spawn();
                        }
                    };
                    run(XAheadBlock {
                        prebuilt,
                        start_next: Some(&start_next),
                        after_prelim_flat: None,
                        table_taken_later: later_blocks >= 2,
                    })
                },
            );
            self.next = built;
            return result;
        }
        let _ = next;
        run(XAheadBlock {
            prebuilt,
            start_next: None,
            after_prelim_flat: None,
            table_taken_later: later_blocks >= 2,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    use crate::algorithm::blastp::args::{BlastpApp, BlastpArgs};
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    use crate::algorithm::blastp::x_query_blocks::XBlockSubjects;
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::x_assert_same_flat_tables_as;
    use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::XFlatCompressedLookup;
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    use crate::algorithm::tblastx::lookup::x_pair_lut::x_assert_same_contexts_as;
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    use crate::algorithm::tblastx::lookup::{x_assert_same_table_as, x_same_table};
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    use std::path::PathBuf;

    // No NCBI counterpart: a check of the unit tests; it changes no value NCBI computes.
    /// The first difference between two query frame sets (every field), or `None`.
    fn x_frames_difference(got: &[Vec<QueryFrame>], want: &[Vec<QueryFrame>]) -> Option<String> {
        if got.len() != want.len() {
            return Some(format!(
                "{} queries' frames built ahead, {} inline",
                got.len(),
                want.len()
            ));
        }
        for (query, (got, want)) in got.iter().zip(want).enumerate() {
            let same = got.len() == want.len()
                && got.iter().zip(want).all(|(a, b)| {
                    a.frame == b.frame
                        && a.aa_seq == b.aa_seq
                        && a.aa_seq_nomask == b.aa_seq_nomask
                        && a.aa_len == b.aa_len
                        && a.orig_len == b.orig_len
                        && a.seg_masks == b.seg_masks
                });
            if !same {
                return Some(format!("the frames of query {query} differ"));
            }
        }
        None
    }

    // No NCBI counterpart: a check of the unit tests; it changes no value NCBI computes.
    /// Whether two flat tables are equal, field for field and slot for slot
    /// (`x_assert_same_flat_tables_as` names the first difference).
    fn x_same_flat_tables(a: &XFlatCompressedLookup, b: &XFlatCompressedLookup) -> bool {
        a.word_length == b.word_length
            && a.compressed_alphabet_size == b.compressed_alphabet_size
            && a.reciprocal_alphabet_size == b.reciprocal_alphabet_size
            && a.longest_chain == b.longest_chain
            && a.pv_array_bts == b.pv_array_bts
            && a.scaled_compress_table_u8 == b.scaled_compress_table_u8
            && a.pv == b.pv
            && a.starts == b.starts
            && a.hits == b.hits
    }

    // No NCBI counterpart: a check of the unit tests; it changes no value NCBI computes.
    /// Whether two context sets are equal, every field, the Karlin blocks bit for bit
    /// (`x_assert_same_contexts_as` names the first difference).
    fn x_same_contexts(a: &[QueryContext], b: &[QueryContext]) -> bool {
        let bits = |context: &QueryContext| {
            let p = &context.karlin_params;
            [
                p.lambda.to_bits(),
                p.k.to_bits(),
                p.h.to_bits(),
                p.alpha.to_bits(),
                p.beta.to_bits(),
            ]
        };
        a.len() == b.len()
            && a.iter().zip(b).all(|(a, b)| {
                a.q_idx == b.q_idx
                    && a.f_idx == b.f_idx
                    && a.frame == b.frame
                    && a.aa_seq == b.aa_seq
                    && a.aa_seq_nomask == b.aa_seq_nomask
                    && a.aa_len == b.aa_len
                    && a.orig_len == b.orig_len
                    && a.frame_base == b.frame_base
                    && a.is_valid == b.is_valid
                    && bits(a) == bits(b)
            })
    }

    // Nested parallel work (two levels of `par_iter`) of block `block`, as a block's search and a
    // build run it on the pool.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    fn nested_work(block: usize) -> u64 {
        use rayon::prelude::*;
        (0..48u64)
            .into_par_iter()
            .map(|x| {
                (0..x + block as u64)
                    .into_par_iter()
                    .map(|y| std::hint::black_box(y * y))
                    .sum::<u64>()
            })
            .sum()
    }

    /// How `fake_loop` builds a block's "table": inline, ahead in one task
    /// (`x_run_with_build_ahead`), or ahead in two parts (`x_run_with_split_build_ahead`, item SC2d)
    /// with the vectors given after the block's nested work (part 1 may end before or after), before
    /// part 1 starts (the vectors come first), or once part 1 has ended (part 1 comes first, mostly).
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[derive(Clone, Copy, Debug, PartialEq)]
    enum XFakeAhead {
        Inline,
        Whole,
        Split,
        SplitVectorsFirst,
        SplitPartFirst,
    }

    /// Sets the flag when dropped (also on a panic).
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    struct XSetOnDrop<'a>(&'a std::sync::atomic::AtomicBool);

    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    impl Drop for XSetOnDrop<'_> {
        fn drop(&mut self) {
            self.0.store(true, Ordering::Release);
        }
    }

    /// A block loop with fake blocks, as `search` drives the blocks: block b's "table" is `build(b)`
    /// (an `Err` for a failing or refused build; it may panic), made by its search inline, or ahead,
    /// during block b - 1, as `ahead` says (in two parts: `build(b)` in part 1, nested work in part
    /// 2); the search then starts the next block's build, runs nested parallel work (its "scan"),
    /// gives the vectors to part 2 of a split build, fails where `fails` says, and writes its row to
    /// `rows`. Returns the loop's outcome (the first error ends it); a panic leaves the loop.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    fn fake_loop(
        threads: usize,
        blocks: usize,
        ahead: XFakeAhead,
        build: &(dyn Fn(usize) -> Result<u64, String> + Sync),
        fails: Option<usize>,
        rows: &Mutex<Vec<String>>,
    ) -> Result<(), String> {
        use std::sync::atomic::AtomicBool;
        crate::utils::threading::with_search_pool(threads, "blastp", |_pool| {
            let mut next: Option<XBuilt<u64>> = None;
            for block in 0..blocks {
                let prebuilt = x_take_built(next.take());
                let took = prebuilt.is_some();
                let part_ended = AtomicBool::new(false);
                // `spawn`: the build's start, and for a split build the vectors' hand-over.
                let search = |spawn: Option<(&dyn Fn(), Option<&dyn Fn()>)>| -> Result<(), String> {
                    let table = match prebuilt {
                        Some(table) => table,
                        None => build(block)?,
                    };
                    let (start, vectors) = spawn.unzip();
                    let vectors = vectors.flatten();
                    if let (XFakeAhead::SplitVectorsFirst, Some(vectors)) = (ahead, vectors) {
                        vectors();
                    }
                    if let Some(start) = start {
                        start();
                    }
                    if ahead == XFakeAhead::SplitPartFirst && vectors.is_some() {
                        while !part_ended.load(Ordering::Acquire) {
                            if !matches!(rayon::yield_now(), Some(rayon::Yield::Executed)) {
                                std::thread::yield_now();
                            }
                        }
                    }
                    let work = nested_work(block);
                    if let (true, Some(vectors)) = (ahead != XFakeAhead::SplitVectorsFirst, vectors)
                    {
                        vectors();
                    }
                    if fails == Some(block) {
                        return Err(format!("block {block} failed"));
                    }
                    rows.lock()
                        .unwrap()
                        .push(format!("block {block} table {table} work {work}"));
                    Ok(())
                };
                let outcome = match ahead {
                    _ if block + 1 == blocks => search(None),
                    XFakeAhead::Inline => search(None),
                    XFakeAhead::Whole => {
                        let (outcome, built) = x_run_with_build_ahead(
                            move || {
                                nested_work(block + 1);
                                build(block + 1).ok()
                            },
                            |spawn| search(Some((spawn, None))),
                        );
                        next = built;
                        outcome
                    }
                    XFakeAhead::Split
                    | XFakeAhead::SplitVectorsFirst
                    | XFakeAhead::SplitPartFirst => {
                        let part_ended = &part_ended;
                        let (outcome, built) = x_run_with_split_build_ahead(
                            move || {
                                let _ended = XSetOnDrop(part_ended);
                                nested_work(block + 1);
                                build(block + 1).ok()
                            },
                            |table| {
                                nested_work(block + 2);
                                Some(table)
                            },
                            |start, vectors| search(Some((start, Some(vectors)))),
                        );
                        next = built;
                        outcome
                    }
                };
                assert!(!took || block > 0, "the first block is never built ahead");
                if let Err(error) = outcome {
                    return Ok(Err(error));
                }
            }
            Ok(Ok(()))
        })
        .expect("pool")
    }

    /// The loop of `fake_loop` without and with the builds ahead (whole and split, item SC2d): the
    /// rows written, the outcome, and a panic's message (`None` without one); all runs must agree.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    fn same_loops(
        threads: usize,
        blocks: usize,
        build: &(dyn Fn(usize) -> Result<u64, String> + Sync),
        fails: Option<usize>,
    ) -> (Vec<String>, Result<(), String>, Option<String>) {
        let run = |ahead: XFakeAhead| {
            let rows = Mutex::new(Vec::new());
            let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                fake_loop(threads, blocks, ahead, build, fails, &rows)
            }));
            let rows = rows.into_inner().unwrap();
            match outcome {
                Ok(outcome) => (rows, outcome, None),
                Err(payload) => {
                    let message = payload
                        .downcast_ref::<String>()
                        .cloned()
                        .or_else(|| payload.downcast_ref::<&str>().map(|text| text.to_string()))
                        .unwrap_or_default();
                    (rows, Ok(()), Some(message))
                }
            }
        };
        let off = run(XFakeAhead::Inline);
        for ahead in [
            XFakeAhead::Whole,
            XFakeAhead::Split,
            XFakeAhead::SplitVectorsFirst,
            XFakeAhead::SplitPartFirst,
        ] {
            let on = run(ahead);
            assert_eq!(
                on, off,
                "threads {threads}, {blocks} blocks, fails {fails:?}, {ahead:?}"
            );
        }
        off
    }

    // Every block after the first takes the values built ahead during the block before it; the rows
    // are those of the loop without the build ahead, at 2, 3 and 7 threads, with one block (nothing
    // built ahead), two and six; nested parallel work in the searches and the builds ends.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn blocks_take_the_values_built_ahead_and_write_the_same_rows() {
        let build = |block: usize| -> Result<u64, String> { Ok(10 * block as u64 + 7) };
        for threads in [2, 3, 7] {
            for blocks in [1, 2, 6] {
                let (rows, outcome, panic) = same_loops(threads, blocks, &build, None);
                assert_eq!(rows.len(), blocks);
                assert_eq!(outcome, Ok(()));
                assert_eq!(panic, None);
            }
        }
    }

    // A build that fails (an error, a refusal) gives no values ahead; its block builds inline and its
    // error ends the loop where it ends without the build ahead: after the rows of the blocks before
    // it, before its own. The first block's error, and the last block's, too.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn a_failing_build_ahead_fails_in_its_block_after_the_rows_before_it() {
        for failing in [0, 1, 3, 5] {
            let build = move |block: usize| -> Result<u64, String> {
                if block == failing {
                    Err(format!("build {block} refused"))
                } else {
                    Ok(block as u64)
                }
            };
            for threads in [2, 7] {
                let (rows, outcome, panic) = same_loops(threads, 6, &build, None);
                assert_eq!(rows.len(), failing);
                assert_eq!(outcome, Err(format!("build {failing} refused")));
                assert_eq!(panic, None);
            }
        }
    }

    // A build that panics: its panic is kept and resumed where its block's search would start, after
    // the rows of the blocks before it, with the same payload as the inline build's panic.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn a_panic_of_a_build_ahead_surfaces_where_its_block_would_raise_it() {
        for panicking in [1, 4] {
            let build = move |block: usize| -> Result<u64, String> {
                if block == panicking {
                    panic!("build {block} panicked");
                }
                Ok(block as u64)
            };
            for threads in [2, 3] {
                let (rows, outcome, panic) = same_loops(threads, 5, &build, None);
                assert_eq!(rows.len(), panicking);
                assert_eq!(outcome, Ok(()));
                assert_eq!(panic, Some(format!("build {panicking} panicked")));
            }
        }
    }

    // A block whose search fails ends the loop with its error; the values built ahead for the next
    // block, an error of that build or its panic are dropped (no panic is resumed).
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn a_failing_block_drops_what_was_built_ahead_for_the_next() {
        let fine = |block: usize| -> Result<u64, String> { Ok(block as u64) };
        let failing = |block: usize| -> Result<u64, String> {
            if block == 3 {
                Err("build 3 refused".to_string())
            } else {
                Ok(block as u64)
            }
        };
        let panicking = |block: usize| -> Result<u64, String> {
            if block == 3 {
                panic!("build 3 panicked");
            }
            Ok(block as u64)
        };
        let builds: [&(dyn Fn(usize) -> Result<u64, String> + Sync); 3] =
            [&fine, &failing, &panicking];
        for build in builds {
            for threads in [2, 7] {
                let (rows, outcome, panic) = same_loops(threads, 5, build, Some(2));
                assert_eq!(rows.len(), 2);
                assert_eq!(outcome, Err("block 2 failed".to_string()));
                assert_eq!(panic, None);
            }
        }
    }

    // `x_run_with_build_ahead` spawns the build once, on the first call of `spawn` only, and gives no
    // outcome when `spawn` is not called.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn the_build_ahead_is_spawned_once_when_asked() {
        use std::sync::atomic::AtomicUsize;
        crate::utils::threading::with_search_pool(3, "blastp", |_pool| {
            let builds = AtomicUsize::new(0);
            let (result, built) = x_run_with_build_ahead(
                || Some(builds.fetch_add(1, Ordering::Relaxed)),
                |spawn| {
                    spawn();
                    spawn();
                    5
                },
            );
            assert_eq!(result, 5);
            assert!(matches!(built, Some(XBuilt::Values(0))));
            assert_eq!(builds.load(Ordering::Relaxed), 1);
            let (_, built) = x_run_with_build_ahead(|| Some(1u8), |_spawn| ());
            assert!(built.is_none());
            let (_, built) = x_run_with_build_ahead(|| None::<u8>, |spawn| spawn());
            assert!(matches!(built, Some(XBuilt::Dropped)));
            Ok(())
        })
        .expect("pool");
    }

    // Item SC2d: the hand-off of a split build ahead starts part 2 exactly once, on the side that
    // comes second, in either order: part 1's values first (the vectors take them), or the vectors
    // first (part 1's task goes on); a second call of the vectors, and the vectors after part 1 gave
    // no values, start nothing.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn the_hand_off_starts_part_two_once_on_the_side_that_comes_second() {
        let slot = Mutex::new(XHandoff::Empty);
        assert_eq!(XHandoff::part(&slot, 7u8), None, "part 1 first: kept");
        assert_eq!(
            XHandoff::vectors(&slot),
            Some(7),
            "then the vectors start part 2"
        );
        assert_eq!(
            XHandoff::vectors(&slot),
            None,
            "a second call starts nothing"
        );
        let slot = Mutex::new(XHandoff::Empty);
        assert_eq!(XHandoff::vectors(&slot), None, "the vectors first: kept");
        assert_eq!(XHandoff::vectors(&slot), None);
        assert_eq!(
            XHandoff::part(&slot, 9u8),
            Some(9),
            "then part 1's task starts part 2"
        );
        assert_eq!(XHandoff::vectors(&slot), None);
        let slot = Mutex::new(XHandoff::<u8>::Empty);
        XHandoff::close(&slot);
        assert_eq!(XHandoff::vectors(&slot), None, "part 1 gave no values");
        assert!(matches!(slot.into_inner().unwrap(), XHandoff::Closed));
    }

    // Item SC2d: `x_run_with_split_build_ahead` runs each part once, with part 1's values, in either
    // order of the hand-over; gives no values where part 1 or part 2 gives none or the vectors never
    // come (part 1's values dropped), keeps a panic of either part, and gives no outcome when part 1
    // is not started.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn the_split_build_ahead_runs_each_part_once_and_keeps_its_outcome() {
        use std::sync::atomic::AtomicUsize;
        crate::utils::threading::with_search_pool(3, "blastp", |_pool| {
            for vectors_first in [false, true] {
                let (firsts, seconds) = (AtomicUsize::new(0), AtomicUsize::new(0));
                let (result, built) = x_run_with_split_build_ahead(
                    || Some(10 + firsts.fetch_add(1, Ordering::Relaxed)),
                    |part| Some(part * 100 + seconds.fetch_add(1, Ordering::Relaxed)),
                    |start, vectors| {
                        if vectors_first {
                            vectors();
                        }
                        start();
                        start();
                        nested_work(3);
                        vectors();
                        vectors();
                        5
                    },
                );
                assert_eq!(result, 5);
                assert!(
                    matches!(built, Some(XBuilt::Values(1000))),
                    "{vectors_first}"
                );
                assert_eq!(firsts.load(Ordering::Relaxed), 1);
                assert_eq!(seconds.load(Ordering::Relaxed), 1);
            }
            let both = |start: &dyn Fn(), vectors: &dyn Fn()| {
                start();
                vectors();
            };
            let (_, built) = x_run_with_split_build_ahead(|| None::<u8>, Some, both);
            assert!(matches!(built, Some(XBuilt::Dropped)), "part 1 gives none");
            let (_, built) = x_run_with_split_build_ahead(|| Some(1u8), |_| None::<u8>, both);
            assert!(matches!(built, Some(XBuilt::Dropped)), "part 2 gives none");
            let (_, built) = x_run_with_split_build_ahead(|| Some(1u8), Some, |start, _| start());
            assert!(matches!(built, Some(XBuilt::Dropped)), "no vectors");
            let (_, built) =
                x_run_with_split_build_ahead(|| Some(1u8), Some, |_, vectors| vectors());
            assert!(built.is_none(), "part 1 not started");
            let panics = |built: Option<XBuilt<u8>>, what: &str| match built {
                Some(XBuilt::Panicked(payload)) => {
                    assert_eq!(payload.downcast_ref::<&str>(), Some(&what));
                }
                _ => panic!("{what}: no panic kept"),
            };
            let (_, built) = x_run_with_split_build_ahead(
                || -> Option<u8> { panic!("part 1 panicked") },
                Some,
                both,
            );
            panics(built, "part 1 panicked");
            for vectors_first in [false, true] {
                let (_, built) = x_run_with_split_build_ahead(
                    || Some(1u8),
                    |_| -> Option<u8> { panic!("part 2 panicked") },
                    |start, vectors| {
                        if vectors_first {
                            vectors();
                        }
                        start();
                        nested_work(2);
                        vectors();
                    },
                );
                panics(built, "part 2 panicked");
            }
            Ok(())
        })
        .expect("pool");
    }

    fn context(lambda: f64) -> QueryContext {
        QueryContext {
            q_idx: 1,
            f_idx: 0,
            frame: 0,
            aa_seq: vec![0, 4, 7, 0],
            aa_seq_nomask: None,
            aa_len: 2,
            orig_len: 2,
            frame_base: 3,
            is_valid: true,
            karlin_params: crate::stats::KarlinParams {
                lambda,
                ..Default::default()
            },
        }
    }

    fn table(hits: Vec<i32>) -> XFlatCompressedLookup {
        XFlatCompressedLookup {
            word_length: 5,
            compressed_alphabet_size: 15,
            reciprocal_alphabet_size: 7,
            longest_chain: 2,
            pv: vec![1, 0],
            pv_array_bts: 5,
            scaled_compress_table_u8: [-1; 256],
            starts: vec![0, 2, 2],
            hits,
        }
    }

    // The comparisons see a difference in one slot of the table's hits, in one bit of a Karlin
    // block, and in one letter of a frame; equal values pass.
    #[test]
    fn the_comparisons_see_one_slot_one_bit_and_one_letter() {
        assert!(x_same_flat_tables(&table(vec![3, 9]), &table(vec![3, 9])));
        assert!(!x_same_flat_tables(&table(vec![3, 9]), &table(vec![3, 8])));
        let lambda = 0.3176;
        assert!(x_same_contexts(&[context(lambda)], &[context(lambda)]));
        let next = f64::from_bits(lambda.to_bits() + 1);
        assert!(!x_same_contexts(&[context(lambda)], &[context(next)]));
        assert!(!x_same_contexts(&[context(lambda)], &[]));
        let frame = crate::algorithm::blastp::encoding::encode_protein_query_frame_with_seg(
            b"MKVLAAGIVG",
            None,
        );
        let mut other = frame.clone();
        assert!(x_frames_difference(&[vec![frame.clone()]], &[vec![other.clone()]]).is_none());
        other.aa_seq[3] ^= 1;
        assert!(x_frames_difference(&[vec![frame.clone()]], &[vec![other]]).is_some());
        assert!(x_frames_difference(&[vec![frame]], &[]).is_some());
    }

    /// The arguments of `losat blastp -task <task>` (SEG on) as the search resolves them.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    fn test_args(task: &str) -> ResolvedBlastpArgs {
        BlastpArgs {
            query: PathBuf::from("query.faa"),
            subject: Some(PathBuf::from("subject.faa")),
            db: None,
            db_filter: Default::default(),
            dbsize: None,
            show_gis: false,
            parse_deflines: false,
            id_scope: None,
            searchsp: None,
            task: task.to_string(),
            evalue: None,
            threshold: None,
            word_size: None,
            num_threads: 1,
            mt_mode: None,
            out: None,
            max_target_seqs: None,
            num_descriptions: None,
            num_alignments: None,
            max_hsps_per_subject: None,
            qcov_hsp_perc: None,
            culling_limit: None,
            ungapped: false,
            window_size: None,
            matrix: None,
            gap_open: None,
            gap_extend: None,
            comp_based_stats: "2".to_string(),
            seg: Some("yes".to_string()),
            soft_masking: None,
            lcase_masking: false,
            use_sw_tback: false,
            query_loc: None,
            subject_loc: None,
            outfmt: "6".to_string(),
            archive: None,
            app: BlastpApp::Blastp,
            inclusion_ethresh: None,
            first_query_ordinal: 0,
            psiblast_pssm: None,
            rps: None,
        }
        .resolve()
        .expect("resolved arguments")
    }

    /// `count` protein records of 40 to 239 residues from a xorshift generator; every third has a
    /// run of one residue (a SEG mask).
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    fn protein_records(seed: u64, count: usize) -> Vec<FastaRecord> {
        let letters = b"ACDEFGHIKLMNPQRSTVWY";
        let mut state = seed;
        let mut below = |n: usize| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % n as u64) as usize
        };
        (0..count)
            .map(|i| {
                let length = 40 + below(200);
                let mut sequence: Vec<u8> = (0..length).map(|_| letters[below(20)]).collect();
                if i % 3 == 0 {
                    let at = below(length - 30);
                    let run = letters[below(20)];
                    sequence[at..at + 24].fill(run);
                }
                FastaRecord::new(format!("Query_{}", i + 1), b"q", &sequence)
            })
            .collect()
    }

    /// The table and contexts of a block as its own search builds them (`run_resolved_in_pool`):
    /// `blastp_aa_lookup` for exact blastp, `blastp_compressed_lookup` with the block's batches and
    /// lengths for blastp-fast.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    fn own_build(
        args: &ResolvedBlastpArgs,
        records: &[FastaRecord],
        frames: &[Vec<QueryFrame>],
        batch_size: u32,
        pool: &SearchPool<'_>,
    ) -> (BlastpRuntimeLookup, Vec<QueryContext>) {
        match args.lookup_table_type {
            BlastpLookupTableType::AaLookupTable => {
                super::super::blastp_aa_lookup(args, frames, false)
            }
            BlastpLookupTableType::CompressedAaLookupTable => {
                let lengths: Vec<usize> = records.iter().map(|record| record.seq().len()).collect();
                let batches = seq_range::QueryInput::whole(records).batches(batch_size);
                super::super::blastp_compressed_lookup(
                    args,
                    frames,
                    &batches,
                    &lengths,
                    pool,
                    XLookupFor::Search,
                    0,
                )
                .expect("the block's own build")
            }
        }
    }

    // A block's values built ahead (as tasks of the pool, on the second spare's vectors made ready
    // as the block loop makes them) are the values the block's own search builds: the query frames,
    // the table and the contexts, with the Karlin blocks bit for bit, on a pool of three threads, SEG
    // masks included. Exact blastp: `x_build_ahead`, once the block before it has its table.
    // blastp-fast (item SC2d): part 1 (`x_flat_ahead_appends`) and part 2 (`x_flat_ahead_table`) on
    // the vectors of the first block's table, given to the second spare as after its preliminary
    // stage (`x_make_room_flat`; the table holds its `starts` buffer): part 1 before the vectors,
    // the vectors before part 1, and both through `x_run_with_split_build_ahead`.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn values_built_ahead_are_the_values_the_block_builds_itself() {
        let batch_size = 3000;
        for (task, order) in [
            ("blastp-fast", 0),
            ("blastp-fast", 1),
            ("blastp-fast", 2),
            ("blastp", 0),
        ] {
            let args = test_args(task);
            let first = protein_records(0x9e37_79b9_7f4a_7c15, 24);
            let next = protein_records(0xd1b5_4a32_d192_ed03, 30);
            crate::utils::threading::with_search_pool(3, "blastp", |pool| {
                let subjects = XBlockSubjects::new(&[], Vec::new());
                let pool = pool.x_with_ext(&subjects);
                let no_spares = "no spares in the pool (LOSAT_X_BPCLUT or LOSAT_X_BPLKPFLAT is \
                                 set to 0)";
                // The block loop: the first block's own table is ready.
                let first_frames = blastp_query_frames(&args, &first);
                let (first_lookup, _) = own_build(&args, &first, &first_frames, batch_size, &pool);
                let mut first_starts = None;
                let built = match first_lookup {
                    BlastpRuntimeLookup::CompressedFlat(first_table) => {
                        first_starts = Some(first_table.starts.as_ptr() as usize);
                        let largest = Cell::new(0);
                        let values = |built: Option<XAheadLookup>| {
                            Some(built.map_or(XBuilt::Dropped, XBuilt::Values))
                        };
                        match order {
                            // Part 1, then the first table's vectors, then part 2.
                            0 => {
                                let part = x_flat_ahead_appends(&args, &next, batch_size, &pool);
                                assert!(
                                    x_make_room_flat(&pool, first_table, &largest),
                                    "{no_spares}"
                                );
                                values(part.and_then(|part| x_flat_ahead_table(part, &pool)))
                            }
                            // The vectors, then part 1 and part 2.
                            1 => {
                                assert!(
                                    x_make_room_flat(&pool, first_table, &largest),
                                    "{no_spares}"
                                );
                                let part = x_flat_ahead_appends(&args, &next, batch_size, &pool);
                                values(part.and_then(|part| x_flat_ahead_table(part, &pool)))
                            }
                            _ => {
                                x_run_with_split_build_ahead(
                                    || x_flat_ahead_appends(&args, &next, batch_size, &pool),
                                    |part| x_flat_ahead_table(part, &pool),
                                    |start, vectors| {
                                        start();
                                        assert!(
                                            x_make_room_flat(&pool, first_table, &largest),
                                            "{no_spares}"
                                        );
                                        vectors();
                                    },
                                )
                                .1
                            }
                        }
                    }
                    // Exact blastp: the vectors for the next block's table are made ready, and its
                    // build is spawned.
                    first_lookup => {
                        assert!(
                            x_make_room(&pool, &first_lookup, &Cell::new(0)),
                            "{no_spares}"
                        );
                        x_run_with_build_ahead(
                            || x_build_ahead(&args, &next, &pool),
                            |spawn| spawn(),
                        )
                        .1
                    }
                };
                let Some(XBuilt::Values(ahead)) = built else {
                    panic!(
                        "{task}: no values built ahead (LOSAT_X_BPLKPAHEAD is 0, LOSAT_X_BPQBLOCK \
                         is changed or LOSAT_DIAGNOSTICS is set)"
                    );
                };
                // The next block's own search.
                let frames = blastp_query_frames(&args, &next);
                assert!(
                    frames
                        .iter()
                        .flatten()
                        .any(|frame| !frame.seg_masks.is_empty()),
                    "{task}: no SEG mask in the fixture"
                );
                if let Some(difference) = x_frames_difference(&ahead.query_frames, &frames) {
                    panic!("{task}: {difference}");
                }
                let (table, contexts) = own_build(&args, &next, &frames, batch_size, &pool);
                assert!(!contexts.is_empty(), "{task}: no contexts");
                match (&ahead.table, &table) {
                    (
                        BlastpRuntimeLookup::CompressedFlat(got),
                        BlastpRuntimeLookup::CompressedFlat(want),
                    ) => {
                        assert!(!want.hits.is_empty(), "{task}: the flat table has no hits");
                        assert_eq!(
                            Some(got.starts.as_ptr() as usize),
                            first_starts,
                            "{task}, order {order}: the table is not on the first table's vectors"
                        );
                        if !x_same_flat_tables(got, want) {
                            x_assert_same_flat_tables_as("ahead:", got, want);
                            panic!("{task}: the tables differ");
                        }
                    }
                    (BlastpRuntimeLookup::Aa(got), BlastpRuntimeLookup::Aa(want)) => {
                        assert!(
                            want.backbone.iter().any(|cell| cell.num_used > 0)
                                || !want.overflow.is_empty(),
                            "{task}: the table has no hits"
                        );
                        if !x_same_table(got, want) {
                            x_assert_same_table_as("ahead:", got, want);
                            panic!("{task}: the tables differ");
                        }
                    }
                    // LOSAT_X_BPLK16: the tables of 16-bit hits, field by field as
                    // `x_same_table`.
                    (BlastpRuntimeLookup::Aa16(got), BlastpRuntimeLookup::Aa16(want)) => {
                        assert!(
                            want.backbone.iter().any(|cell| cell.num_used > 0)
                                || !want.overflow.is_empty(),
                            "{task}: the table has no hits"
                        );
                        assert!(
                            got.backbone.len() == want.backbone.len()
                                && got.backbone.iter().zip(&want.backbone).all(|(a, b)| {
                                    a.num_used == b.num_used && a.entries == b.entries
                                })
                                && got.overflow == want.overflow
                                && got.pv == want.pv
                                && got.frame_bases == want.frame_bases
                                && got.num_contexts == want.num_contexts
                                && got.query_length == want.query_length
                                && got.word_length == want.word_length
                                && got.alphabet_size == want.alphabet_size
                                && got.charsize == want.charsize
                                && got.mask == want.mask
                                && got.longest_chain == want.longest_chain
                                && got.threshold == want.threshold
                                && got.row_max == want.row_max,
                            "{task}: the tables of 16-bit hits differ"
                        );
                    }
                    _ => panic!("{task}: the tables are of different kinds"),
                }
                if !x_same_contexts(&ahead.contexts, &contexts) {
                    x_assert_same_contexts_as("ahead:", &ahead.contexts, &contexts);
                    panic!("{task}: the contexts differ");
                }
                Ok(())
            })
            .expect("pool");
        }
    }

    // A block's search takes values built ahead only for its own records (the same slice); others are
    // dropped and the search builds its own.
    #[test]
    fn values_built_ahead_are_taken_only_by_their_block() {
        let records = FastaRecord::new("Query_1", b"q1", b"MKVLAAGIVG");
        let records = vec![records.clone(), records];
        let values = |key| XAheadLookup {
            key,
            query_frames: Vec::new(),
            table: BlastpRuntimeLookup::CompressedFlat(table(Vec::new())),
            contexts: Vec::new(),
        };
        let mut block = XAheadBlock {
            prebuilt: Some(values(x_key(&records))),
            start_next: None,
            after_prelim_flat: None,
            table_taken_later: true,
        };
        assert!(block.take_prebuilt(&records[..1]).is_none());
        let mut block = XAheadBlock {
            prebuilt: Some(values(x_key(&records))),
            start_next: None,
            after_prelim_flat: None,
            table_taken_later: true,
        };
        assert!(block.take_prebuilt(&records).is_some());
        assert!(block.take_prebuilt(&records).is_none());
        assert!(XAheadBlock::none().take_prebuilt(&records).is_none());
    }
}
