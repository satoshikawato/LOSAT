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
//! Why the reports do not change: block i + 1's frames, table and contexts are functions of block i
//! + 1's query records and the search's arguments only. The task computes them by the calls the
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
//! uses at a time. No floating-point operation changes order: each context's Karlin block is the
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
//! block i + 1's. In the block search such errors do not arise (`x_block_words_fit` sends a query set
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
//! has returned. A thread that waits inside block i's nested parallel work may steal and run the
//! task (or one of its passes) to the end, and a thread that waits inside the task's passes may run
//! one of block i's jobs to the end; neither waits for the frame it interrupted, so both waits end
//! (memory `rayon-waiting-helper-deadlock`: the hazard there was a spawned loop that waited for
//! progress made elsewhere). The scope's owner, the block's thread, runs the task itself if no
//! other thread took it. The locks the task takes (the spare's, the result's slot) are held only
//! around a move.
//!
//! Memory: two tables are alive at once (block i's, which its search scans, and block i + 1's), so
//! the peak holds one more table's `starts`, `hits` and PV array (at `LOSAT_X_BPQBLOCK=50000` about
//! 3 MiB + up to 19 MiB; exact blastp: backbone 0.44 MiB, overflow about 7-8 MiB, PV 3.5 KiB), plus
//! block i + 1's query frames and contexts. All of the table's vectors are allocated by the block
//! search's thread, never by the pool's task (a pool thread's allocations stay in its own glibc
//! arena, which keeps them resident after they are freed): when block i's table is ready, that
//! thread makes the second spare hold vectors with the room for the largest table so far plus a
//! half (`XFlatSpare::x_make_room`: the second spare's own, else those block i - 1 gave back, else
//! new ones of twice the largest table's hits; for exact blastp `XAaSpares::x_make_room`, the same
//! rule for the overflow entries, and each block's table is given back at the end of its search),
//! and the task's build takes them all (`XFlatArgs::x_ahead`, `XAaSpares::x_build_ahead`). A build
//! whose `hits` (overflow) would not fit gives no values and the block builds its table itself, on
//! its own thread. Untouched capacity is not resident. The exact build's neighbour lists (about
//! 1.5-1.9 MiB) are the block search's thread's, moved to `XAaSpares` once and used by every build
//! ahead, so no pool thread keeps a set of its own.
//!
//! `LOSAT_X_BPLKPAHEAD` is a switch of the strict set (on by default since 2026-10-17): on unless it is `0`.
//!
//! With `LOSAT_X_STATS`, `main` prints `[X_STATS] BPLKPAHEAD blocks= prebuilt= dropped=` at exit
//! when a counter is not zero (`blocks`: blocks searched where the switch acts; `prebuilt`: blocks
//! that took values built ahead; `dropped`: builds ahead that gave no values). The unit test
//! `values_built_ahead_are_the_values_the_block_builds_itself` builds a block's values ahead and
//! inline, for blastp-fast and exact blastp, and compares the query frames, the table and the
//! contexts (the Karlin blocks bit for bit).
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
use std::sync::{Mutex, OnceLock, PoisonError};

use super::{
    blastp_aa_lookup_ahead, blastp_compressed_lookup, blastp_query_frames, BlastpRuntimeLookup,
};
use crate::algorithm::blastp::args::{BlastpLookupTableType, ResolvedBlastpArgs};
use crate::algorithm::blastp::x_query_blocks::{
    x_bp_aa_block_spare, x_bp_aa_spares, x_bp_lkreuse, x_bp_lookup_spares, XLookupFor,
};
use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::x_bpclut_on;
use crate::algorithm::tblastx::lookup::{QueryContext, XAaSpares};
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::blastinput::fasta_reader::FastaRecord;
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
                    let outcome =
                        match std::panic::catch_unwind(std::panic::AssertUnwindSafe(build)) {
                            Ok(Some(values)) => XBuilt::Values(values),
                            Ok(None) => XBuilt::Dropped,
                            Err(payload) => XBuilt::Panicked(payload),
                        };
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

// No NCBI counterpart: the floor of the room for the next table's hits (module header); untouched
// capacity is not resident. It changes no value NCBI computes.
/// The fewest hits (exact blastp: overflow entries) the second spare has room for: 2^18 (1 MiB).
const X_MIN_HITS_ROOM: usize = 1 << 18;

// No NCBI counterpart: which thread allocates the next table's vectors (module header); it changes
// no value NCBI computes.
/// On the block's thread, when its table `lookup` is ready: makes the second spare (exact blastp:
/// `XAaSpares`) hold vectors with the room for the next block's table (the largest table so far,
/// `largest_hits` (exact blastp: overflow entries), plus a half; at least `X_MIN_HITS_ROOM`),
/// allocating here what is short (twice the largest). `false` where the block's table is neither a
/// flat one nor an exact one with its spares (no build ahead then).
fn x_make_room(
    pool: &SearchPool<'_>,
    lookup: &BlastpRuntimeLookup,
    largest_hits: &Cell<usize>,
) -> bool {
    match lookup {
        BlastpRuntimeLookup::CompressedFlat(table) => {
            let Some((block, ahead)) = x_bp_lookup_spares(pool) else {
                return false;
            };
            let largest = largest_hits.get().max(table.hits.len());
            largest_hits.set(largest);
            let room = (largest + largest / 2).max(X_MIN_HITS_ROOM);
            ahead.x_make_room(
                block,
                table.starts.len().saturating_sub(1),
                room,
                (2 * largest).max(room),
            );
            true
        }
        // Exact blastp: the same rule for the overflow entries; the backbone and the pv array have
        // the sizes of the block's own table (fixed by the word length and the alphabet).
        BlastpRuntimeLookup::Aa(table) => {
            let Some(spares) = x_bp_aa_spares(pool) else {
                return false;
            };
            let largest = largest_hits.get().max(table.overflow.len());
            largest_hits.set(largest);
            let room = (largest + largest / 2).max(X_MIN_HITS_ROOM);
            spares.x_make_room(
                table.backbone.len(),
                table.pv.len(),
                room,
                (2 * largest).max(room),
            );
            true
        }
        // EXPERIMENT (LOSAT_X_BPLK16): the same for a table of 16-bit hits; the lengths of its
        // vectors count their padding, which the next table of that width has as well.
        BlastpRuntimeLookup::Aa16(table) => {
            let Some(spares) = x_bp_aa_spares(pool) else {
                return false;
            };
            let largest = largest_hits.get().max(table.overflow.len());
            largest_hits.set(largest);
            let room = (largest + largest / 2).max(X_MIN_HITS_ROOM);
            spares.x_make_room16(
                table.backbone.len(),
                table.pv.len(),
                room,
                (2 * largest).max(room),
            );
            true
        }
        BlastpRuntimeLookup::Compressed(_) => false,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:206-215 (module header)
// ```c++
//         LookupTableWrap* lut =
//             CSetupFactory::CreateLookupTable(query_data, opts_memento.get(),
//                                              sbp, lookup_segments_wrap,
// ```
/// EXPERIMENT (LOSAT_X_BPLKPAHEAD): the query frames, table and contexts of a block of
/// `query_records`, as its search builds them (module header: the same calls with the same inputs;
/// `batch_size` is the search's query batch size), its table on the second spare's vectors
/// (`XLookupFor::Ahead`; exact blastp: `XAaSpares`). `None` where the build fails, refuses or makes
/// no flat table (exact blastp: makes no table on the spare's vectors).
fn x_build_ahead(
    args: &ResolvedBlastpArgs,
    query_records: &[FastaRecord],
    batch_size: u32,
    pool: &SearchPool<'_>,
) -> Option<XAheadLookup> {
    let query_frames = blastp_query_frames(args, query_records);
    let (table, contexts) = match args.lookup_table_type {
        // Exact blastp: `blastp_aa_lookup`'s call on the block search's spare vectors.
        BlastpLookupTableType::AaLookupTable => {
            blastp_aa_lookup_ahead(args, &query_frames, x_bp_aa_spares(pool)?)?
        }
        BlastpLookupTableType::CompressedAaLookupTable => {
            let query_lengths: Vec<usize> = query_records
                .iter()
                .map(|record| record.seq().len())
                .collect();
            let query_batches = seq_range::QueryInput::whole(query_records).batches(batch_size);
            let (lookup, contexts) = blastp_compressed_lookup(
                args,
                &query_frames,
                &query_batches,
                &query_lengths,
                pool,
                XLookupFor::Ahead,
            )
            .ok()?;
            if !matches!(lookup, BlastpRuntimeLookup::CompressedFlat(_)) {
                return None;
            }
            (lookup, contexts)
        }
    };
    Some(XAheadLookup {
        key: x_key(query_records),
        query_frames,
        table,
        contexts,
    })
}

/// EXPERIMENT (LOSAT_X_BPLKPAHEAD): the block loop's side of the switch (module header): whether it
/// acts for this block search, and what the build ahead of the next block gave.
pub(super) struct XLookupAhead<'a, 'p> {
    args: &'a ResolvedBlastpArgs,
    pool: &'a SearchPool<'p>,
    batch_size: u32,
    acts: bool,
    next: Option<XBuilt<XAheadLookup>>,
    /// The most hits of a table so far (`x_make_room`).
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
                BlastpLookupTableType::CompressedAaLookupTable => {
                    x_bpclut_on() && x_bp_lookup_spares(pool).is_some()
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
            let (result, built) = x_run_with_build_ahead(
                move || x_build_ahead(args, next, batch_size, pool),
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
            table_taken_later: later_blocks >= 2,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::blastp::args::{BlastpApp, BlastpArgs};
    use crate::algorithm::blastp::x_query_blocks::XBlockSubjects;
    use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::{
        x_assert_same_flat_tables_as, XFlatCompressedLookup,
    };
    use crate::algorithm::tblastx::lookup::x_pair_lut::x_assert_same_contexts_as;
    use crate::algorithm::tblastx::lookup::{x_assert_same_table_as, x_same_table};
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
            && a.used_cells == b.used_cells
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

    /// A block loop with fake blocks, as `search` drives the blocks: block b's "table" is `build(b)`
    /// (an `Err` for a failing or refused build; it may panic), made by its search inline, or ahead,
    /// during block b - 1, with `ahead`; the search then starts the next block's build, runs nested
    /// parallel work, fails where `fails` says, and writes its row to `rows`. Returns the loop's
    /// outcome (the first error ends it); a panic leaves the loop.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    fn fake_loop(
        threads: usize,
        blocks: usize,
        ahead: bool,
        build: &(dyn Fn(usize) -> Result<u64, String> + Sync),
        fails: Option<usize>,
        rows: &Mutex<Vec<String>>,
    ) -> Result<(), String> {
        crate::utils::threading::with_search_pool(threads, "blastp", |_pool| {
            let mut next: Option<XBuilt<u64>> = None;
            for block in 0..blocks {
                let prebuilt = x_take_built(next.take());
                let took = prebuilt.is_some();
                let search = |spawn: Option<&dyn Fn()>| -> Result<(), String> {
                    let table = match prebuilt {
                        Some(table) => table,
                        None => build(block)?,
                    };
                    if let Some(spawn) = spawn {
                        spawn();
                    }
                    let work = nested_work(block);
                    if fails == Some(block) {
                        return Err(format!("block {block} failed"));
                    }
                    rows.lock()
                        .unwrap()
                        .push(format!("block {block} table {table} work {work}"));
                    Ok(())
                };
                let outcome = if ahead && block + 1 < blocks {
                    let (outcome, built) = x_run_with_build_ahead(
                        move || {
                            nested_work(block + 1);
                            build(block + 1).ok()
                        },
                        |spawn| search(Some(spawn)),
                    );
                    next = built;
                    outcome
                } else {
                    search(None)
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

    /// The loop of `fake_loop` without and with the build ahead: the rows written, the outcome, and
    /// a panic's message (`None` without one); both runs must agree.
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
        let run = |ahead: bool| {
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
        let off = run(false);
        let on = run(true);
        assert_eq!(
            on, off,
            "threads {threads}, {blocks} blocks, fails {fails:?}"
        );
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
            used_cells: 1,
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
                blastp_compressed_lookup(args, frames, &batches, &lengths, pool, XLookupFor::Search)
                    .expect("the block's own build")
            }
        }
    }

    // A block's values built ahead (`x_build_ahead` as a task of the pool, on the second spare's
    // vectors made ready as the block loop makes them, once the block before it has its table) are
    // the values the block's own search builds: the query frames, the table and the contexts, with
    // the Karlin blocks bit for bit, for blastp-fast (the flat table) and for exact blastp, on a
    // pool of three threads, SEG masks included.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn values_built_ahead_are_the_values_the_block_builds_itself() {
        let batch_size = 3000;
        for task in ["blastp-fast", "blastp"] {
            let args = test_args(task);
            let first = protein_records(0x9e37_79b9_7f4a_7c15, 24);
            let next = protein_records(0xd1b5_4a32_d192_ed03, 30);
            crate::utils::threading::with_search_pool(3, "blastp", |pool| {
                let subjects = XBlockSubjects::new(&[], Vec::new());
                let pool = pool.x_with_ext(&subjects);
                // The block loop: the first block's own table is ready, the vectors for the next
                // block's table are made ready, and its build is spawned.
                let first_frames = blastp_query_frames(&args, &first);
                let (first_lookup, _) = own_build(&args, &first, &first_frames, batch_size, &pool);
                assert!(
                    x_make_room(&pool, &first_lookup, &Cell::new(0)),
                    "{task}: no spares in the pool (LOSAT_X_BPCLUT or LOSAT_X_BPLKPFLAT is set \
                     to 0)"
                );
                let (_, built) = x_run_with_build_ahead(
                    || x_build_ahead(&args, &next, batch_size, &pool),
                    |spawn| spawn(),
                );
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
            table_taken_later: true,
        };
        assert!(block.take_prebuilt(&records[..1]).is_none());
        let mut block = XAheadBlock {
            prebuilt: Some(values(x_key(&records))),
            start_next: None,
            table_taken_later: true,
        };
        assert!(block.take_prebuilt(&records).is_some());
        assert!(block.take_prebuilt(&records).is_none());
        assert!(XAheadBlock::none().take_prebuilt(&records).is_none());
    }
}
