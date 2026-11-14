//! RPSTBLASTN's search of BLASTX's query batches (series PR, item PS12a, DECISIONS 38 D1): the six
//! frames of each nucleotide query, set up as BLASTX sets them up (`prepare_unit_queries`: frames,
//! strands, genetic codes, masks, `-query_loc`; `score_block` and `effective_lengths`: the frames'
//! Karlin-Altschul blocks and search spaces), are the "subjects" of the reversed RPS search
//! (`rps_search::translated`); its final HSPs become BLASTX's results (`blastx_lists`), which
//! BLASTX's reports write. The search takes the place of BLASTX's own (`RpsTblastnSearch`).
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/rpstblastn_app.cpp:191-209
//! ```c++
//!         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//!
//!             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//!             CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
//!
//!             SaveSearchStrategy(args, m_CmdLineArgs, queries, opts_hndl);
//!
//!             CRef<CSearchResultSet> results;
//!
//!             if (m_CmdLineArgs->ExecuteRemotely()) {
//!                 ...
//!             } else {
//!             	CLocalRPSBlast  local_search (query_batch, db_args->GetDatabaseName(), opts_hndl, args[kArgNumThreads].AsInteger() );
//!             	results = local_search.Run();
//!             }
//! ```

use std::io::Write;
use std::sync::Arc;

use anyhow::{bail, Result};

use crate::algorithm::blastp::hsp::{BlastpHitList, BlastpHsp};
use crate::algorithm::blastx::args::{ResolvedOptions, RpsTblastnSearch};
use crate::algorithm::blastx::database::BlastxDatabase;
use crate::algorithm::blastx::input::FastaRecord;
use crate::algorithm::blastx::native::{engine, BatchRender, QueryBatchReader};
use crate::algorithm::blastx::parameters::{
    effective_lengths, gapped_karlin, gumbel_block, score_block, ContextParameters,
};
use crate::algorithm::blastx::preliminary::PreliminaryHsp;
use crate::algorithm::blastx::query_setup::{
    prepare_unit_queries, unknown_matrix_message, PreparedQueryBatch,
};
use crate::algorithm::blastx::results::{Hsp, HspList};
use crate::algorithm::blastx::runtime::BatchResults;
use crate::algorithm::rps_search::engine::{
    RpsFrame, RpsQuery, RpsRun, RpsSearch, RpsWorkers, NUM_FRAMES,
};
use crate::algorithm::rps_search::volumes::{
    mask_unsupported_letters, merge_volume, per_volume_search_spaces, shift_ordinal_ids, Consensus,
    RpsVolumeError, RpsVolumeSet, RpsVolumeSetup,
};
use crate::core::blast_util::frame_to_context;
use crate::stats::KarlinParams;

/// The program of LOSAT's refusals.
const PROGRAM: &str = "RPSTBLASTN";

/// The RPS search of an RPSTBLASTN run: `RpsSearch` (rpsblast's options of the RPS engine),
/// `-ungapped`, and what the run's database tells (`for_database`). The formatting options are
/// BLASTX's (`BlastxArgs::formatting`; `ResolvedOptions::report_sort` keeps NCBI's no-sort).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct TranslatedRpsSearch {
    pub(crate) search: RpsSearch,
    /// `-ungapped` (with composition-based statistics 0: NCBI's options error at the others is
    /// BLASTX's `resolve`).
    pub(crate) ungapped: bool,
    /// `-use_sw_tback` and `-comp_based_stats` as given: the composition-based statistics of
    /// later items are refused per run (`check_options`).
    pub(crate) use_sw_tback: bool,
    pub(crate) comp_based_stats: String,
    /// Whether the database has several volumes, whose results NCBI merges
    /// (`RpsTblastnSearch::merges_volumes`); set by `for_database`.
    pub(crate) several_volumes: bool,
    /// The volumes' matrix and gap costs as their set-up read them (`RpsVolumeSetup::costs`), for
    /// the search strategy's requests (`handle_sets`); set by `for_database` of a run that writes
    /// one.
    pub(crate) costs: Option<(String, i32, i32)>,
}

impl RpsTblastnSearch for TranslatedRpsSearch {
    fn for_database(
        &self,
        options: &ResolvedOptions,
        database: Option<&BlastxDatabase>,
        strategy: bool,
    ) -> Result<Arc<dyn RpsTblastnSearch>> {
        let Some(database) = database else {
            bail!("a search without -db is not supported by LOSAT's {PROGRAM}");
        };
        // LOSAT's refusals of the volumes' files come before the first byte, as rpsblast's do;
        // NCBI's own errors of the files stay where NCBI's search of a batch raises them
        // (`search_window`), and the options' refusals below wait for the files to open.
        let volumes = RpsVolumeSet::of(database.x_seqdb());
        let mut setup = RpsVolumeSetup::new(PROGRAM, options.window_size);
        let mut faulty_files = false;
        for &index in &volumes.order {
            if setup.check(&volumes.volumes[index], false)?.is_err() {
                faulty_files = true;
            }
        }
        if !faulty_files {
            self.check_options()?;
        }
        // The volumes' costs are read for the search strategy's requests only (`handle_sets`;
        // after volume files that failed NCBI's checks, `costs` opens the first volume's again).
        Ok(Arc::new(Self {
            several_volumes: volumes.order.len() > 1,
            costs: strategy.then(|| setup.costs(volumes.first())).flatten(),
            ..self.clone()
        }))
    }

    fn run_query_batches(
        &self,
        options: &ResolvedOptions,
        database: Option<&BlastxDatabase>,
        read: &mut QueryBatchReader<'_>,
        diagnostics: &mut dyn Write,
        search_error: &dyn Fn(anyhow::Error) -> anyhow::Error,
        render: &mut BatchRender<'_>,
        take: &mut dyn FnMut(&[FastaRecord], usize) -> Result<()>,
    ) -> Result<usize> {
        let Some(database) = database else {
            bail!("a search without -db is not supported by LOSAT's {PROGRAM}");
        };
        // The queries of the batches replayed so far: each batch's search strategy request and
        // render take the number before it.
        let processed = std::cell::Cell::new(0usize);
        let diagnostics = std::cell::RefCell::new(diagnostics);
        let volumes = RpsVolumeSet::of(database.x_seqdb());
        let setup = RpsVolumeSetup::new(PROGRAM, options.window_size);
        let lengths = database_lengths(database);
        // `-num_threads 0`: one thread for each volume (NCBI's "auto", as rpsblast's).
        let threads = if self.search.auto_threads {
            volumes.order.len()
        } else {
            usize::try_from(options.num_threads.max(1)).unwrap_or(1)
        };
        // The letters the reader reads ahead of the searches (D1 of PS12s): four reader batches
        // (`BATCH_SIZE` letters each, NCBI's batch size of rpstblastn by default) for each
        // thread. The bound only sets how far the reader reads ahead: every input and thread
        // count takes the same path, and a long record, a batch of its own, makes a window of at
        // most the bound plus that batch. (`native::run` refused a `BATCH_SIZE` it cannot read.)
        let batch_size = crate::blastinput::app::query_batch_size(
            crate::algorithm::blastx::query_setup::RPSTBLASTN_BATCH_SIZE as u32,
        )? as usize;
        let bound = batch_size.saturating_mul(threads).saturating_mul(4);
        crate::utils::threading::x_with_search_pool_local(threads, PROGRAM, |pool| {
            // The workers' scratch, kept from window to window and volume to volume
            // (`RpsWorkers`).
            let workers = RpsWorkers::new(pool.threads());
            let window = std::cell::RefCell::new(RpsWindow::default());
            let failed = std::cell::Cell::new(false);
            // The window's batches searched together (`search_window`), then its events
            // replayed in the reader's order (`replay`): what NCBI's one-batch loop writes, in
            // its order, up to the first batch that fails, which ends the run with its error.
            // Each batch's search strategy request (`take`, `-export_search_strategy`) is
            // written in the replay, before the batch's report or error, as NCBI writes it
            // before the batch's search (`SaveSearchStrategy`, quoted in the module's header):
            // after the batches before it are reported (their search spaces, which the request
            // holds, `StrategyExport::note_searched`), and none for a batch after one that fails.
            let mut flush = || -> Result<()> {
                let events = std::mem::take(&mut *window.borrow_mut()).events;
                let batches: Vec<&[FastaRecord]> = events
                    .iter()
                    .filter_map(|event| match event {
                        RpsEvent::Batch(queries) => Some(queries.as_slice()),
                        RpsEvent::Warn(_) => None,
                    })
                    .collect();
                let searched = self.search_window(
                    &batches, options, database, &volumes, &setup, lengths, pool, &workers,
                );
                replay(
                    &events,
                    searched,
                    &mut |bytes| Ok(diagnostics.borrow_mut().write_all(bytes)?),
                    &mut |queries| take(queries, processed.get()),
                    &mut |queries, results| {
                        render(
                            &mut **diagnostics.borrow_mut(),
                            queries,
                            &mut vec![results],
                            processed.get(),
                        )?;
                        processed.set(processed.get() + queries.len());
                        Ok(())
                    },
                )
                .map_err(|error| match error.downcast::<RpsBatchError>() {
                    // NCBI's errors of the set-up of the batch (`BLAST engine error`).
                    Ok(RpsBatchError::Setup(error)) => search_error(error),
                    Ok(RpsBatchError::Engine(text)) => engine(text),
                    Err(error) => error,
                })
                .inspect_err(|_| failed.set(true))
            };
            let read_result = read(
                &mut |queries| {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
                    // ```c++
                    // CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
                    //     : m_QueryVector(& queries)
                    // {
                    //     if (queries.Empty()) {
                    //         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
                    //     }
                    // ```
                    // (After the batches before it, which NCBI's loop has written.)
                    if queries.is_empty() {
                        flush()?;
                        return Err(engine("Empty CBlastQueryVector"));
                    }
                    let full = {
                        let mut window = window.borrow_mut();
                        window.letters = queries.iter().fold(window.letters, |sum, q| {
                            sum.saturating_add(q.sequence.len())
                        });
                        window.events.push(RpsEvent::Batch(queries.to_vec()));
                        window.letters >= bound
                    };
                    if full {
                        flush()?;
                    }
                    Ok(())
                },
                &mut |bytes: &[u8]| {
                    window
                        .borrow_mut()
                        .events
                        .push(RpsEvent::Warn(bytes.to_vec()));
                    Ok(())
                },
            );
            // What the reader delivered before its end or its error, searched and replayed
            // before that error; nothing after a batch that failed.
            if !failed.get() {
                flush()?;
            }
            read_result
        })?;
        Ok(processed.get())
    }

    fn merges_volumes(&self) -> bool {
        self.several_volumes
    }

    fn handle_sets(
        &self,
        database: Option<&BlastxDatabase>,
        searchsp: Option<i64>,
        dbsize: Option<i64>,
    ) -> Vec<crate::report::archive::options::Param> {
        database.map_or_else(Vec::new, |database| {
            crate::algorithm::rps_search::volumes::search_sets(
                database.x_seqdb(),
                self.costs.clone(),
                searchsp,
                dbsize,
            )
        })
    }
}

/// What ends the search of a query batch besides LOSAT's refusals: NCBI's error of the batch's
/// set-up (BLASTX's, `prepare_unit_queries`), or NCBI's engine error with its text.
#[derive(Debug)]
enum RpsBatchError {
    Setup(anyhow::Error),
    Engine(String),
}

impl std::fmt::Display for RpsBatchError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Setup(error) => write!(f, "{error}"),
            Self::Engine(text) => f.write_str(text),
        }
    }
}

impl std::error::Error for RpsBatchError {}

/// The database's length and number of sequences for the effective lengths of the frames: the
/// whole profile database's (every volume's search reads the totals, `x_AdjustDbSize` without
/// `-dbsize`); the database is protein, so the length is not divided by 3.
fn database_lengths(database: &BlastxDatabase) -> (i64, i64) {
    let pairwise = database.pairwise_database();
    (pairwise.total_letters as i64, pairwise.num_sequences as i64)
}

impl TranslatedRpsSearch {
    /// The searches of a window of query batches (D1 of PS12s), each batch in NCBI's order of
    /// its set-up (`blast_aux_priv.cpp` `BlastSetupPreliminarySearchEx`): the queries
    /// (`prepare_unit_queries`), then for each volume its RPS files (NCBI's errors) and the score
    /// blocks of the frames (`score_block`), the batch searched in the volume only where a frame
    /// is valid (NCBI's `CheckInternalData`). A batch without a valid frame keeps BLASTX's
    /// unsearched results: with one volume its queries get the Karlin-Altschul warning, with
    /// several the merge drops the warnings and gives each query an empty alignment set
    /// (`ResolvedOptions::merged_volumes`).
    ///
    /// The window is searched volume by volume, each volume mapped once for every batch of the
    /// window, whose queries search it together (`RpsRun::search_volume`: every frame of every
    /// query a unit of the run's pool), as NCBI's `CLocalRPSBlast` searches each batch alone:
    /// every value a batch holds (its queries' frames, statistics, masks, `-qcov_hsp_perc`'s
    /// first query, merged lists) is the batch's own. Returns each batch's results or the error
    /// NCBI's loop meets first in it (its set-up, then volume by volume: the files, the search
    /// spaces, the volume's refusals, then the first query in order that fails); the batches
    /// after the first that fails are not searched (`None`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:158-164
    /// ```c++
    ///     // 2. Take care of any rps information
    ///     if (Blast_ProgramIsRpsBlast(options->GetProgramType())) {
    ///         const char* name = BlastSeqSrcGetName(seqsrc);
    ///         const string rps_dbname(name ? name : "");
    ///         retval->m_InternalData->m_RpsData =
    ///             CSetupFactory::CreateRpsStructures(rps_dbname, options);
    ///     }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:176-181
    /// ```c++
    ///
    ///     int status = m_PrelimSearch->CheckInternalData();
    ///     if (status != 0)
    ///     {
    ///          // Search was not run, but we send back an empty CSearchResultSet.
    ///          CRef<ILocalQueryData> local_query_data = m_QueryFactory->MakeLocalQueryData(m_Opts);
    /// ```
    #[allow(clippy::too_many_arguments)]
    fn search_window(
        &self,
        batches: &[&[FastaRecord]],
        options: &ResolvedOptions,
        database: &BlastxDatabase,
        volumes: &RpsVolumeSet<'_>,
        setup: &RpsVolumeSetup,
        (db_length, db_count): (i64, i64),
        pool: &crate::utils::threading::SearchPool<'_>,
        workers: &RpsWorkers,
    ) -> Vec<Option<Result<BatchResults>>> {
        let mut outcome: Vec<Option<Result<BatchResults>>> = batches.iter().map(|_| None).collect();
        if batches.is_empty() {
            return outcome;
        }
        // The batches [0, states.len()) have not failed; a failing batch ends the window's search
        // there (`fail`), as its error ends NCBI's loop.
        let mut states: Vec<WindowBatch> = Vec::with_capacity(batches.len());
        let fail = |states: &mut Vec<WindowBatch>,
                    outcome: &mut Vec<Option<Result<BatchResults>>>,
                    batch: usize,
                    error: anyhow::Error| {
            states.truncate(batch);
            outcome[batch] = Some(Err(error));
            // A later batch that failed before is not reached (NCBI's loop ends here).
            for later in &mut outcome[batch + 1..] {
                *later = None;
            }
        };
        for (batch, prepared) in prepare_batches(batches, options, pool)
            .into_iter()
            .enumerate()
        {
            match prepared {
                Ok(prepared) => states.push(WindowBatch::new(prepared, batches[batch].len())),
                Err(error) => {
                    fail(
                        &mut states,
                        &mut outcome,
                        batch,
                        RpsBatchError::Setup(error).into(),
                    );
                    break;
                }
            }
        }
        // BLOSUM62 (the `.aux`'s matrix, `RpsVolumeSetup::check`): BLASTX's ideal block, which
        // `check_ideal` takes for RPSTBLASTN as for BLASTX.
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2746-2748
        // ```c
        //    Boolean check_ideal =
        //       (program == eBlastTypeBlastx || program == eBlastTypeTblastx ||
        //        program == eBlastTypeRpsTblastn);
        // ```
        let gapped = gapped_karlin(options);
        let gumbel = gumbel_block(options, db_length);
        // `avg_subject_length` of the hit cutoff: the shortest profile (`SeqDb::min_length`, as
        // rpsblast's `cutoff_subject_length`).
        let cutoff_subject_length = database.x_seqdb().min_length();
        let x_drops = self.search.x_drop_bits();
        let per_volume = per_volume_search_spaces(
            options.db_length,
            options.effective_search_space,
            volumes.order.len(),
        );
        // Volume after volume, each mapped for the window's search alone, as NCBI's
        // `CLocalRPSBlast` searches each query batch.
        for (position, &index) in volumes.order.iter().enumerate() {
            if states.is_empty() {
                break;
            }
            let volume = &volumes.volumes[index];
            // The volume's files: the same for every batch, so their error is the first batch's.
            let files = match setup.reopen(volume, false) {
                Ok(Ok(files)) => files,
                Ok(Err(RpsVolumeError::Files(text))) => {
                    fail(
                        &mut states,
                        &mut outcome,
                        0,
                        RpsBatchError::Engine(text).into(),
                    );
                    break;
                }
                // The empty matrix name fails every query's score block (blastx's unknown
                // matrix: one "Unknown error code -1" per query).
                Ok(Err(RpsVolumeError::NoMatrix)) => {
                    let text = unknown_matrix_message(batches[0]);
                    fail(
                        &mut states,
                        &mut outcome,
                        0,
                        RpsBatchError::Engine(text).into(),
                    );
                    break;
                }
                Err(refusal) => {
                    fail(&mut states, &mut outcome, 0, refusal);
                    break;
                }
            };
            // Each batch's set-up for the volume, in the reader's order: the batches with a
            // valid frame search it, with their frames' letters masked for its alphabet.
            let mut searching: Vec<usize> = Vec::new();
            let mut masks: Vec<Vec<Option<Vec<u8>>>> = Vec::new();
            for batch in 0..states.len() {
                let WindowBatch {
                    prepared,
                    parameters,
                    reported,
                    searched,
                    ..
                } = &mut states[batch];
                let parameters = parameters.get_or_insert_with(|| score_block(prepared));
                if !parameters.iter().any(|p| p.valid) {
                    continue;
                }
                // `-dbsize` on several volumes: each volume's statistics read the volume's own
                // number of sequences (`per_volume_search_spaces`, as rpsblast's).
                let lengths = if per_volume {
                    let num_seqs = i64::from(volume.num_oids());
                    effective_lengths(prepared, parameters, options, db_length, num_seqs, None)
                        .inspect(|()| {
                            // The merged result keeps the first volume's statistics (in path
                            // order), which the reports print (`RpsVolumeSet::first`, as
                            // rpsblast's footers).
                            if position == 0 {
                                *reported = Some(parameters.clone());
                            }
                        })
                } else if !*searched {
                    effective_lengths(prepared, parameters, options, db_length, db_count, None)
                } else {
                    Ok(())
                };
                if let Err(error) = lengths {
                    fail(
                        &mut states,
                        &mut outcome,
                        batch,
                        RpsBatchError::Setup(error).into(),
                    );
                    break;
                }
                // `Blast_MaskUnsupportedAA` of a 26-column database, per frame.
                masks.push(
                    prepared
                        .contexts
                        .iter()
                        .map(|context| {
                            (context.length > 0).then(|| {
                                let letters = &prepared.sequence_start
                                    [context.offset..context.offset + context.length + 2];
                                mask_unsupported_letters(letters, files.alphabet)
                            })?
                        })
                        .collect(),
                );
                searching.push(batch);
            }
            if searching.is_empty() {
                continue;
            }
            // The queries of the searching batches, batch after batch, in one run.
            let mut queries = Vec::new();
            let mut spans = Vec::with_capacity(searching.len());
            for (&batch, masked) in searching.iter().zip(&masks) {
                let state = &states[batch];
                let start = queries.len();
                queries.extend(frame_queries(
                    &state.prepared,
                    state
                        .parameters
                        .as_deref()
                        .expect("a searching batch has its score blocks"),
                    masked,
                ));
                spans.push(queries.len() - start);
            }
            let run = RpsRun {
                queries,
                gapped,
                gumbel: gumbel.as_ref(),
                cutoff_subject_length,
                evalue: options.evalue,
                prelim_evalue: options.evalue,
                cbs: 0,
                sbp: None,
                hitlist_size: options.hitlist_size.max(1) as usize,
                prelim_hitlist_size: crate::algorithm::blastx::results::prelim_hitlist_size(
                    options,
                ),
                window: options.window_size,
                x_drop_ungapped_bits: x_drops.ungapped,
                x_drop_gapped_bits: x_drops.gapped,
                x_drop_final_bits: x_drops.final_,
                gap_trigger_bits: x_drops.gap_trigger,
                options: &self.search,
                program: PROGRAM,
            };
            let traced = run.search_volume(
                &files,
                &Consensus::new(setup, volume, &files),
                pool,
                workers,
            );
            drop(run);
            drop(masks);
            drop(files);
            let traced = match traced {
                Ok(traced) => traced,
                // The volume's refusal: the first batch that searches it meets it. The batches
                // before it (without a valid frame) still set up the later volumes, whose NCBI
                // errors they meet first (review of PS12s, B1).
                Err(refusal) => {
                    fail(&mut states, &mut outcome, searching[0], refusal);
                    continue;
                }
            };
            // Each batch's lists, or its error in NCBI's order: the first refusal of its queries'
            // cutoffs, the volume's refusal after them, then the first error of its queries'
            // searches in query order (`RpsVolumeResults::batches`).
            for (&batch, lists) in searching.iter().zip(traced.batches(&spans)) {
                let lists = match lists {
                    Ok(lists) => lists,
                    Err(error) => {
                        fail(&mut states, &mut outcome, batch, error);
                        break;
                    }
                };
                let mut lists: Vec<Option<BlastpHitList>> = lists
                    .into_iter()
                    .map(|traced| {
                        let mut list = traced.hit_list;
                        if let Some(list) = list.as_mut() {
                            shift_ordinal_ids(list, volumes.starts[index]);
                            list.hsplist_array.truncate(list.hsplist_count);
                        }
                        list
                    })
                    .collect();
                let state = &mut states[batch];
                if position == 0 {
                    std::mem::swap(&mut state.hit_lists, &mut lists);
                } else {
                    merge_volume(&mut state.hit_lists, lists);
                }
                state.searched = true;
            }
        }
        for (batch, state) in states.into_iter().enumerate() {
            outcome[batch] = Some(Ok(state.results()));
        }
        outcome
    }

    /// LOSAT's refusals of the options before the first byte of the reports: what NCBI BLAST+
    /// 2.17.0's rpstblastn does in its search with them has no result that can be checked.
    pub(crate) fn check_options(&self) -> Result<()> {
        let requested_cbs = self.search.requested_cbs;
        if self.use_sw_tback && requested_cbs > 0 {
            bail!(
                "-use_sw_tback with composition-based statistics is not supported by LOSAT's {PROGRAM} yet"
            );
        }
        // NCBI's options error of an ungapped search with composition-based statistics comes
        // from BLASTX's `resolve`; the refusal of the modes LOSAT does not run yet follows it.
        if requested_cbs > 0 && !self.ungapped {
            bail!(
                "-comp_based_stats {} is not supported by LOSAT's {PROGRAM} yet",
                self.comp_based_stats
            );
        }
        if self.search.sum_stats {
            bail!(
                "-sum_stats true (NCBI BLAST+ 2.17.0's rpstblastn links HSPs with the profiles' indexes as query indexes and reads past its arrays: its output depends on memory contents or it crashes) is not supported by LOSAT's {PROGRAM}"
            );
        }
        // DECISIONS 38 D3: every probed query (`out/u.*` of the design) segfaults after the
        // search, valgrind: a read at address 0 in `BLAST_ComputeTraceback_MT`.
        if self.ungapped {
            bail!("{}", crate::cli::run_refusal(PROGRAM, "ungapped"));
        }
        Ok(())
    }
}

/// One of the reader's events held in a window (`run_query_batches`, D1 of PS12s): a warning as
/// its bytes, a query batch as a copy of its records.
enum RpsEvent {
    Warn(Vec<u8>),
    Batch(Vec<FastaRecord>),
}

/// The reader's events not yet searched and written, and the letters of their batches.
#[derive(Default)]
struct RpsWindow {
    events: Vec<RpsEvent>,
    letters: usize,
}

/// The events of a window in the reader's order (D1 of PS12s): each warning written (`warn`),
/// each batch's search strategy request taken (`take`), then its results reported (`report`)
/// with its records, from `searched` (one entry per batch, `search_window`); a batch that failed
/// ends the replay with its error, after its request and the events before it, as it ended
/// NCBI's loop of one batch at a time.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/rpstblastn_app.cpp:193-208
/// ```c++
///             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
///             CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
///
///             SaveSearchStrategy(args, m_CmdLineArgs, queries, opts_hndl);
///             ...
///             	CLocalRPSBlast  local_search (query_batch, db_args->GetDatabaseName(), opts_hndl, args[kArgNumThreads].AsInteger() );
///             	results = local_search.Run();
/// ```
fn replay<T>(
    events: &[RpsEvent],
    searched: Vec<Option<Result<T>>>,
    warn: &mut dyn FnMut(&[u8]) -> Result<()>,
    take: &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
    report: &mut dyn FnMut(&[FastaRecord], T) -> Result<()>,
) -> Result<()> {
    let mut searched = searched.into_iter();
    for event in events {
        match event {
            RpsEvent::Warn(bytes) => warn(bytes)?,
            RpsEvent::Batch(queries) => {
                // The batch's request, before its search's result or error is looked at.
                take(queries)?;
                let results = searched
                    .next()
                    .flatten()
                    .expect("each batch up to the first that failed is searched")?;
                report(queries, results)?;
            }
        }
    }
    Ok(())
}

/// What the search of a window keeps for each of its query batches (D1 of PS12s): the values
/// NCBI's search of the batch alone holds (`search_window`).
struct WindowBatch {
    prepared: PreparedQueryBatch,
    /// The frames' score blocks (`score_block`, after the first volume's files open) and search
    /// spaces (`effective_lengths`).
    parameters: Option<Vec<ContextParameters>>,
    /// With `-dbsize` on several volumes: the first volume's (the merged result's statistics).
    reported: Option<Vec<ContextParameters>>,
    /// Whether a volume searched the batch.
    searched: bool,
    hit_lists: Vec<Option<BlastpHitList>>,
}

impl WindowBatch {
    fn new(prepared: PreparedQueryBatch, queries: usize) -> Self {
        Self {
            prepared,
            parameters: None,
            reported: None,
            searched: false,
            hit_lists: vec![None; queries],
        }
    }

    /// The batch's results (`unsearched_batch` where no volume searched it).
    fn results(self) -> BatchResults {
        let parameters = self
            .reported
            .or(self.parameters)
            .unwrap_or_else(|| score_block(&self.prepared));
        if !self.searched {
            return unsearched_batch(&self.prepared, parameters);
        }
        BatchResults {
            query_ordinal: 0,
            queries: self
                .hit_lists
                .into_iter()
                .enumerate()
                .map(|(query, list)| blastx_lists(list, query))
                .collect(),
            prepared: self.prepared.x_for_results(),
            parameters,
            search_skipped: false,
        }
    }
}

/// Each batch of a window prepared (`prepare_batch`), the batches at the same time on the run's
/// pool; each batch's result is its own.
fn prepare_batches(
    batches: &[&[FastaRecord]],
    options: &ResolvedOptions,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> Vec<Result<PreparedQueryBatch>> {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if pool.enabled() && pool.threads() > 1 {
        use rayon::prelude::*;
        return pool.install(|| {
            batches
                .par_iter()
                .with_max_len(1)
                .map(|queries| prepare_batch(queries, options))
                .collect()
        });
    }
    let _ = pool;
    batches
        .iter()
        .map(|queries| prepare_batch(queries, options))
        .collect()
}

/// A batch's set-up before its RPS structures: the query splitter's `CHUNK_SIZE`, then the
/// queries (`prepare_unit_queries`: frames, masks, `-query_loc`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:59-61
/// ```c++
///     char* chunk_sz_str = getenv("CHUNK_SIZE");
///     if (chunk_sz_str && !NStr::IsBlank(chunk_sz_str)) {
///         retval = NStr::StringToInt(chunk_sz_str);
/// ```
/// NCBI reads it in the set-up of every query batch (`CQuerySplitter`), before its RPS
/// structures; RPSTBLASTN never splits a query (`split_query_aux_priv.cpp:88`), so a value it
/// reads changes nothing, and one it cannot convert is its `CStringException`
/// (`app::read_query_split_variable`), which the batch's error frame passes on as it is
/// (`native::database_search_error`): the first batch stops here, after the outfmt 0 prolog,
/// the messages of its records and its search strategy request (`replay`).
fn prepare_batch(queries: &[FastaRecord], options: &ResolvedOptions) -> Result<PreparedQueryBatch> {
    crate::blastinput::app::read_query_split_variable("CHUNK_SIZE")?;
    prepare_unit_queries(queries, None, options)
}

/// The results of a batch that is not searched: no list for any query, BLASTX's unsearched
/// batch (its reports warn for each query, unless the run merges volumes).
fn unsearched_batch(
    prepared: &PreparedQueryBatch,
    parameters: Vec<ContextParameters>,
) -> BatchResults {
    BatchResults {
        query_ordinal: 0,
        queries: vec![Vec::new(); prepared.original_lengths.len()],
        prepared: prepared.x_for_results(),
        parameters,
        search_skipped: true,
    }
}

/// The frames of each query of a prepared batch as the RPS engine reads them, with the
/// statistics of their contexts (`score_block`, `effective_lengths`): BLASTX's layout of
/// `sequence_start` holds each frame with letters between two sentinels at its context's offset;
/// `masked` holds a frame's letters where a 26-column database masks some
/// (`mask_unsupported_letters`). `-qcov_hsp_perc` reads the batch's first query's frames.
pub(crate) fn frame_queries<'a>(
    prepared: &'a PreparedQueryBatch,
    parameters: &[ContextParameters],
    masked: &'a [Option<Vec<u8>>],
) -> Vec<RpsQuery<'a>> {
    let frame = |context: usize| -> RpsFrame<'a> {
        let c = &prepared.contexts[context];
        let p = &parameters[context];
        let letters = |buffer: &'a [u8]| -> &'a [u8] {
            if c.length == 0 {
                &[]
            } else {
                &buffer[c.offset..c.offset + c.length + 2]
            }
        };
        let searched = letters(&prepared.sequence_start);
        let nomask = letters(&prepared.sequence_start_nomask);
        let aa_seq = masked[context].as_deref().unwrap_or(searched);
        RpsFrame {
            frame: c.frame,
            aa_seq,
            aa_seq_nomask: (nomask != aa_seq).then_some(nomask),
            valid: p.valid && c.length > 0,
            kbp_std: p.ungapped,
            search_space: p.search_space,
            first_query_length: prepared.contexts[context % NUM_FRAMES].length,
        }
    };
    let batch_kbp = parameters
        .iter()
        .find(|p| p.valid)
        .map_or_else(KarlinParams::default, |p| p.ungapped);
    (0..prepared.original_lengths.len())
        .map(|query| RpsQuery {
            frames: (0..NUM_FRAMES)
                .map(|i| frame(query * NUM_FRAMES + i))
                .collect(),
            batch_kbp,
        })
        .collect()
}

/// An RPS final HSP as BLASTX's results hold it: the roles are back (the query side is the
/// user query's frame, in protein letters; the subject side the profile's consensus), so the
/// coordinates become BLASTX's 0-based half-open ones on the frame and the profile, the context
/// is the query's frame context, and the scores, the e-value, the counts, the edit script and
/// the composition-based method are the traceback's.
fn blastx_hsp(hsp: &BlastpHsp, query_index: usize) -> Hsp {
    let to_i32 = |value: usize| i32::try_from(value).expect("a position below 2^31");
    Hsp {
        hsp: PreliminaryHsp {
            context: query_index * NUM_FRAMES + frame_to_context(hsp.query_frame),
            frame: hsp.query_frame as i8,
            score: hsp.raw_score,
            q_start: to_i32(hsp.q_start) - 1,
            q_end: to_i32(hsp.q_end),
            q_gapped_start: hsp.gapped_q_start,
            s_start: to_i32(hsp.s_start) - 1,
            s_end: to_i32(hsp.s_end),
            s_gapped_start: hsp.gapped_s_start,
        },
        // No sum statistics (`do_sum_stats` is FALSE for RPSTBLASTN): no linked set.
        num: 0,
        evalue: hsp.e_value,
        bit_score: hsp.bit_score,
        identity: hsp.num_ident,
        positive: hsp.num_positives,
        edit_script: hsp.gap_info.clone().unwrap_or_default(),
        composition_method: i32::from(hsp.comp_adjustment_method),
    }
}

/// A query's final hit list as BLASTX's results: one list per profile in the hit list's order,
/// each with its HSPs in their order, named by the profile's ordinal id in the database.
fn blastx_lists(hit_list: Option<BlastpHitList>, query_index: usize) -> Vec<HspList> {
    let Some(hit_list) = hit_list else {
        return Vec::new();
    };
    hit_list
        .hsplist_array
        .into_iter()
        .take(hit_list.hsplist_count)
        .filter(|list| !list.hsps.is_empty())
        .map(|list| {
            let mut converted = HspList {
                query_index,
                oid: list.oid as usize,
                hsps: list
                    .hsps
                    .iter()
                    .map(|hsp| blastx_hsp(hsp, query_index))
                    .collect(),
                best_evalue: list.best_evalue,
            };
            converted.refresh();
            converted
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::blastp::hsp::BlastpHspList;
    use crate::common::GapEditOp;

    fn options(extra: &[&str]) -> ResolvedOptions {
        let mut argv = vec!["losat", "rpstblastn", "-query", "q.fna", "-db", "cdd"];
        argv.extend_from_slice(&["-comp_based_stats", "0"]);
        argv.extend_from_slice(extra);
        let cli: crate::cli::Cli = crate::cli::try_parse_from(argv).unwrap();
        let crate::cli::Commands::Rpstblastn(args) = cli.command else {
            unreachable!()
        };
        args.into_blastx_args().unwrap().resolve().unwrap()
    }

    fn record(id: &str, letters: &[u8]) -> FastaRecord {
        FastaRecord {
            internal_id: id.to_string(),
            title: String::new(),
            sequence: letters.to_vec(),
            lowercase_masks: Vec::new(),
            warnings: Vec::new(),
            source: None,
        }
    }

    /// 60 letters of a coding sequence (20 codons) whose six frames all have letters.
    const DNA: &[u8] = b"ATGGCTAAAGGTGAACGTCTGATCGCTGTTGCTAAAGAACTGGGTCTGACCGCTAAAGAA";

    // blast_engine.c:805-813 (quoted at `RpsFrame`): the frames of each query in NCBI's context
    // order, their letters at BLASTX's offsets between two sentinels, and the frames of the strand
    // that `-strand` leaves out without letters.
    #[test]
    fn the_frames_of_a_query_are_blastxs_contexts() {
        let options = options(&[]);
        let records = [record("a", DNA), record("b", &DNA[..31])];
        let prepared = prepare_unit_queries(&records, None, &options).unwrap();
        let parameters = score_block(&prepared);
        let masked = vec![None; prepared.contexts.len()];
        let queries = frame_queries(&prepared, &parameters, &masked);
        assert_eq!(queries.len(), 2);
        let frames: Vec<i8> = queries[0].frames.iter().map(|f| f.frame).collect();
        assert_eq!(frames, [1, 2, 3, -1, -2, -3]);
        let lengths: Vec<usize> = queries[0].frames.iter().map(RpsFrame::length).collect();
        assert_eq!(lengths, [20, 19, 19, 20, 19, 19]);
        let lengths: Vec<usize> = queries[1].frames.iter().map(RpsFrame::length).collect();
        assert_eq!(lengths, [10, 10, 9, 10, 10, 9]);
        for (query, frame_query) in queries.iter().enumerate() {
            for (i, frame) in frame_query.frames.iter().enumerate() {
                let context = &prepared.contexts[query * NUM_FRAMES + i];
                assert_eq!(frame.aa_seq.first(), Some(&0));
                assert_eq!(frame.aa_seq.last(), Some(&0));
                assert_eq!(
                    &frame.aa_seq[1..frame.aa_seq.len() - 1],
                    &prepared.sequence_start
                        [context.offset + 1..context.offset + 1 + context.length]
                );
                assert!(frame.valid);
                assert_eq!(frame.aa_seq_nomask, None);
                assert_eq!(frame_to_context(i32::from(frame.frame)), i);
            }
        }
        // The first frame translates the coding sequence (BLASTX's genetic code 1).
        let translated: Vec<u8> = queries[0].frames[0].aa_seq[1..21]
            .iter()
            .map(|&letter| b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ"[usize::from(letter)])
            .collect();
        assert_eq!(translated, b"MAKGERLIAVAKELGLTAKE");
        let plus = options_and_frames(&["-strand", "plus"]);
        assert_eq!(plus, [20, 19, 19, 0, 0, 0]);
        let minus = options_and_frames(&["-strand", "minus"]);
        assert_eq!(minus, [0, 0, 0, 20, 19, 19]);
    }

    fn options_and_frames(extra: &[&str]) -> Vec<usize> {
        let options = options(extra);
        let records = [record("a", DNA)];
        let prepared = prepare_unit_queries(&records, None, &options).unwrap();
        let parameters = score_block(&prepared);
        let masked = vec![None; prepared.contexts.len()];
        let queries = frame_queries(&prepared, &parameters, &masked);
        for frame in &queries[0].frames {
            assert_eq!(frame.valid, frame.length() > 0);
        }
        queries[0].frames.iter().map(RpsFrame::length).collect()
    }

    // blast_engine.c:813 (quoted at `RpsFrame`): the last frame with letters.
    #[test]
    fn the_preliminary_stat_length_is_the_last_frame_with_letters() {
        let options_both = options(&[]);
        let records = [record("a", &DNA[..59])];
        let prepared = prepare_unit_queries(&records, None, &options_both).unwrap();
        let parameters = score_block(&prepared);
        let masked = vec![None; prepared.contexts.len()];
        let queries = frame_queries(&prepared, &parameters, &masked);
        // Frames of 59 letters: +1 19, +2 19, +3 19, -1 19, -2 19, -3 19 ((59 - i) / 3).
        assert_eq!(queries[0].stat_length(), queries[0].frames[5].length());
        let options_plus = options(&["-strand", "plus"]);
        let prepared = prepare_unit_queries(&[record("b", DNA)], None, &options_plus).unwrap();
        let parameters = score_block(&prepared);
        let masked = vec![None; prepared.contexts.len()];
        let queries = frame_queries(&prepared, &parameters, &masked);
        assert_eq!(queries[0].stat_length(), 19);
        // A query of 4 letters: one letter in frames +1, -1 only (frames +2, +3: (4 - 1) / 3 =
        // 1, (4 - 2) / 3 = 0 letters).
        let prepared = prepare_unit_queries(&[record("c", b"ACGT")], None, &options_both).unwrap();
        let parameters = score_block(&prepared);
        let masked = vec![None; prepared.contexts.len()];
        let queries = frame_queries(&prepared, &parameters, &masked);
        let lengths: Vec<usize> = queries[0].frames.iter().map(RpsFrame::length).collect();
        assert_eq!(lengths, [1, 1, 0, 1, 1, 0]);
        assert_eq!(queries[0].stat_length(), 1);
    }

    // blast_setup_cxx.cpp:72-77 and blast_query_info.c:379-381 (quoted at
    // `RpsQuery::translation_length`): each context starts one sentinel after the previous one,
    // an empty context takes no byte, and the buffer's length is the last context's offset plus
    // its length; BLASTX's layout (`prepare_unit_queries`) is NCBI's.
    #[test]
    fn the_translation_buffer_has_ncbis_layout() {
        let layout = |extra: &[&str], letters: &[u8]| -> (Vec<(usize, usize)>, usize) {
            let options = options(extra);
            let prepared = prepare_unit_queries(&[record("a", letters)], None, &options).unwrap();
            let parameters = score_block(&prepared);
            let masked = vec![None; prepared.contexts.len()];
            let queries = frame_queries(&prepared, &parameters, &masked);
            let first = prepared.contexts[0].offset;
            let contexts: Vec<(usize, usize)> = prepared.contexts[..NUM_FRAMES]
                .iter()
                .map(|c| (c.offset - first, c.length))
                .collect();
            let (offset, length) = contexts[NUM_FRAMES - 1];
            assert_eq!(
                queries[0].translation_length(),
                offset + length,
                "{extra:?}"
            );
            (contexts, queries[0].translation_length())
        };
        // Both strands: six frames of 20, 19, 19 letters.
        let (both, length) = layout(&[], DNA);
        assert_eq!(
            both,
            [(0, 20), (21, 19), (41, 19), (61, 20), (82, 19), (102, 19)]
        );
        assert_eq!(length, 121);
        // -strand plus: the minus frames are empty and start where the last plus frame ends.
        let (plus, length) = layout(&["-strand", "plus"], DNA);
        assert_eq!(
            plus,
            [(0, 20), (21, 19), (41, 19), (61, 0), (61, 0), (61, 0)]
        );
        assert_eq!(length, 61);
        // -strand minus: the plus frames take no byte.
        let (minus, length) = layout(&["-strand", "minus"], DNA);
        assert_eq!(minus, [(0, 0), (0, 0), (0, 0), (0, 20), (21, 19), (41, 19)]);
        assert_eq!(length, 60);
        // A query of 4 letters: frames +3 and -3 are empty.
        let (short, length) = layout(&[], b"ACGT");
        assert_eq!(short, [(0, 1), (2, 1), (4, 0), (4, 1), (6, 1), (8, 0)]);
        assert_eq!(length, 8);
    }

    fn hsp(frame: i32, score: i32, q: (usize, usize), s: (usize, usize)) -> BlastpHsp {
        BlastpHsp {
            num: 1,
            identity: 50.0,
            length: q.1 - q.0 + 2,
            mismatch: 1,
            gapopen: 1,
            q_start: q.0,
            q_end: q.1,
            s_start: s.0,
            s_end: s.1,
            gapped_q_start: (q.0 + 2) as i32,
            gapped_s_start: (s.0 + 2) as i32,
            e_value: 1e-5,
            bit_score: 30.5,
            num_ident: 7,
            query_context: 0,
            query_frame: frame,
            subject_frame: 0,
            query_length: 100,
            q_idx: 0,
            s_idx: 0,
            raw_score: score,
            // The query's columns (1-based closed), with one deletion after the fifth.
            gap_info: Some(match (q.1 - q.0 + 1) as u32 {
                columns if columns > 5 => vec![
                    GapEditOp::Sub(5),
                    GapEditOp::Del(1),
                    GapEditOp::Sub(columns - 5),
                ],
                columns => vec![GapEditOp::Sub(columns)],
            }),
            num_positives: 9,
            comp_adjustment_method: 0,
        }
    }

    // The conversion of the RPS traceback's final HSPs (1-based closed coordinates after the role
    // swap, `traceback::swap_roles`) into BLASTX's (0-based half-open, context = query * 6 +
    // the frame's context).
    #[test]
    fn rps_hsps_become_blastx_hsps() {
        let converted = blastx_hsp(&hsp(-2, 41, (11, 30), (3, 23)), 2);
        assert_eq!(
            converted.hsp,
            PreliminaryHsp {
                context: 2 * 6 + 4,
                frame: -2,
                score: 41,
                q_start: 10,
                q_end: 30,
                q_gapped_start: 13,
                s_start: 2,
                s_end: 23,
                s_gapped_start: 5,
            }
        );
        assert_eq!(
            (
                converted.evalue,
                converted.bit_score,
                converted.identity,
                converted.positive,
                converted.num,
                converted.composition_method
            ),
            (1e-5, 30.5, 7, 9, 0, 0)
        );
        assert_eq!(
            converted.edit_script,
            [GapEditOp::Sub(5), GapEditOp::Del(1), GapEditOp::Sub(15)]
        );
        // Every frame's context.
        for (frame, context) in [(1, 0), (2, 1), (3, 2), (-1, 3), (-2, 4), (-3, 5)] {
            assert_eq!(
                blastx_hsp(&hsp(frame, 1, (1, 2), (1, 2)), 0).hsp.context,
                context
            );
        }
        let mut method = hsp(1, 1, (1, 2), (1, 2));
        method.comp_adjustment_method = 1;
        assert_eq!(blastx_hsp(&method, 0).composition_method, 1);
    }

    // A query's hit list: one BLASTX list per profile with HSPs, in the hit list's order, named
    // by the profile's ordinal id; the best e-value refreshed (`HspList::refresh`).
    #[test]
    fn a_hit_list_becomes_blastx_lists_in_its_order() {
        let list = |oid: u32, evalues: &[f64]| BlastpHspList {
            oid,
            query_index: 0,
            hsps: evalues
                .iter()
                .map(|&e| BlastpHsp {
                    num: 1,
                    e_value: e,
                    ..hsp(1, 10, (1, 9), (1, 9))
                })
                .collect(),
            best_evalue: 0.0,
        };
        let hit_list = BlastpHitList {
            hsplist_count: 3,
            hsplist_max: 3,
            worst_evalue: 0.0,
            low_score: 0,
            heapified: false,
            hsplist_array: vec![
                list(5, &[1e-3, 1e-9]),
                list(1, &[]),
                list(2, &[0.5]),
                list(9, &[1e-30]),
            ],
            hsplist_current: 0,
            num_hits: 0,
        };
        let lists = blastx_lists(Some(hit_list), 3);
        let oids: Vec<usize> = lists.iter().map(|l| l.oid).collect();
        assert_eq!(oids, [5, 2]);
        assert!(lists.iter().all(|l| l.query_index == 3));
        assert_eq!(lists[0].best_evalue, 1e-9);
        assert_eq!(lists[0].hsps.len(), 2);
        assert_eq!(lists[0].hsps[0].hsp.context, 3 * 6);
        assert!(blastx_lists(None, 0).is_empty());
    }

    // local_blast.cpp:176-181 (quoted at `search_window`): a batch without a valid frame is not
    // searched; its results are BLASTX's skipped batch.
    #[test]
    fn a_batch_without_a_valid_frame_is_not_searched() {
        let options = options(&[]);
        let records = [record("n", b"NNNNNN"), record("a", b"AC")];
        let prepared = prepare_unit_queries(&records, None, &options).unwrap();
        let parameters = score_block(&prepared);
        assert!(!parameters.iter().any(|p| p.valid));
        let results = unsearched_batch(&prepared, parameters);
        assert!(results.search_skipped);
        assert_eq!(results.queries.len(), 2);
        assert!(results.queries.iter().all(Vec::is_empty));
        assert_eq!(results.parameters.len(), 12);
    }

    // The refusals before the reports (`for_database`).
    #[test]
    fn sum_statistics_and_an_ungapped_search_are_refused_before_the_reports() {
        let mut search = TranslatedRpsSearch {
            search: RpsSearch::default(),
            ungapped: false,
            use_sw_tback: false,
            comp_based_stats: "0".to_string(),
            several_volumes: false,
            costs: None,
        };
        assert!(search.check_options().is_ok());
        search.search.sum_stats = true;
        assert!(search
            .check_options()
            .unwrap_err()
            .to_string()
            .contains("-sum_stats true"));
        search.search.sum_stats = false;
        search.ungapped = true;
        assert!(search
            .check_options()
            .unwrap_err()
            .to_string()
            .contains("crashes after an ungapped search"));
        assert!(!search.merges_volumes());
    }

    // D1 of PS12s: a window's events are replayed in the reader's order, each warning written and
    // each batch reported, up to the first batch that failed, whose error ends the replay after the
    // events before it (the warnings and batches after it are not written, as NCBI's loop never
    // reads them); a window without a failure writes everything.
    #[test]
    fn a_window_replays_the_readers_events_up_to_the_first_failed_batch() {
        let events = vec![
            RpsEvent::Warn(b"w0\n".to_vec()),
            RpsEvent::Batch(vec![record("a", DNA)]),
            RpsEvent::Batch(vec![record("b", DNA), record("c", DNA)]),
            RpsEvent::Warn(b"w1\n".to_vec()),
            RpsEvent::Batch(vec![record("d", DNA)]),
            RpsEvent::Warn(b"w2\n".to_vec()),
            RpsEvent::Batch(vec![record("e", DNA)]),
        ];
        // The log of a replay: warnings, each batch's request (`take:` and its first query) and
        // report; `take_fails` names the first query of a batch whose request fails.
        let run = |searched: Vec<Option<Result<u32>>>,
                   take_fails: Option<&str>|
         -> (Vec<String>, Option<String>) {
            let log = std::cell::RefCell::new(Vec::new());
            let result = replay(
                &events,
                searched,
                &mut |bytes| {
                    log.borrow_mut()
                        .push(String::from_utf8_lossy(bytes).trim().to_string());
                    Ok(())
                },
                &mut |queries| {
                    let first = queries[0].internal_id.as_str();
                    if take_fails == Some(first) {
                        anyhow::bail!("request {first} failed");
                    }
                    log.borrow_mut().push(format!("take:{first}"));
                    Ok(())
                },
                &mut |queries, results| {
                    let ids: Vec<&str> = queries.iter().map(|q| q.internal_id.as_str()).collect();
                    log.borrow_mut()
                        .push(format!("{}={results}", ids.join("+")));
                    Ok(())
                },
            );
            (
                log.into_inner(),
                result.err().map(|error| error.to_string()),
            )
        };
        let lines =
            |lines: &[&str]| -> Vec<String> { lines.iter().map(|line| line.to_string()).collect() };
        let all = run(
            vec![Some(Ok(1)), Some(Ok(2)), Some(Ok(3)), Some(Ok(4))],
            None,
        );
        let expected = lines(&[
            "w0", "take:a", "a=1", "take:b", "b+c=2", "w1", "take:d", "d=3", "w2", "take:e", "e=4",
        ]);
        assert_eq!(all, (expected, None));
        // A failing batch has its request; the batches after it have none.
        let failed = run(
            vec![
                Some(Ok(1)),
                Some(Ok(2)),
                Some(Err(anyhow::anyhow!("batch 3 failed"))),
                None,
            ],
            None,
        );
        assert_eq!(
            failed,
            (
                lines(&["w0", "take:a", "a=1", "take:b", "b+c=2", "w1", "take:d"]),
                Some("batch 3 failed".to_string())
            )
        );
        let first = run(
            vec![
                Some(Err(anyhow::anyhow!("batch 1 failed"))),
                None,
                None,
                None,
            ],
            None,
        );
        assert_eq!(
            first,
            (lines(&["w0", "take:a"]), Some("batch 1 failed".to_string()))
        );
        // A batch's failing request comes after the batches before it and before the batch's
        // own error.
        let request = run(
            vec![
                Some(Ok(1)),
                Some(Ok(2)),
                Some(Err(anyhow::anyhow!("batch 3 failed"))),
                None,
            ],
            Some("d"),
        );
        assert_eq!(
            request,
            (
                lines(&["w0", "take:a", "a=1", "take:b", "b+c=2", "w1"]),
                Some("request d failed".to_string())
            )
        );
    }

    /// The files of the fixture database `name` (tests/db_search/rps/NAME.*.gz) inflated into a
    /// directory of this test but those whose names end with one of `left_out`, and the database
    /// opened.
    fn rps_database(
        name: &str,
        left_out: &[&str],
        test: &str,
    ) -> (std::path::PathBuf, BlastxDatabase) {
        use std::io::Read;
        let fixtures = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/rps");
        let dir = std::env::temp_dir().join(format!("losat-rpstbn-{test}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        for entry in std::fs::read_dir(&fixtures).unwrap() {
            let file = entry.unwrap().file_name().into_string().unwrap();
            let Some(inflated) = file.strip_suffix(".gz") else {
                continue;
            };
            let Some(extension) = inflated.strip_prefix(name) else {
                continue;
            };
            if !extension.starts_with('.') || left_out.contains(&extension) {
                continue;
            }
            let mut bytes = Vec::new();
            flate2::read::GzDecoder::new(std::fs::read(fixtures.join(&file)).unwrap().as_slice())
                .read_to_end(&mut bytes)
                .unwrap();
            std::fs::write(dir.join(inflated), bytes).unwrap();
        }
        let database = BlastxDatabase::open(
            dir.join(name).to_str().unwrap(),
            &crate::blastinput::ncbi_environment::ApplicationSettings::default(),
        )
        .unwrap();
        (dir, database)
    }

    /// `search_window` of `batches` with the RPS options `search` on `database`, serially.
    fn window_errors(
        search: &TranslatedRpsSearch,
        database: &BlastxDatabase,
        batches: &[&[FastaRecord]],
    ) -> Vec<Option<std::result::Result<(), String>>> {
        let options = options(&[]);
        let volumes = RpsVolumeSet::of(database.x_seqdb());
        let setup = RpsVolumeSetup::new(PROGRAM, options.window_size);
        search
            .search_window(
                batches,
                &options,
                database,
                &volumes,
                &setup,
                database_lengths(database),
                &crate::utils::threading::SearchPool::x_serial(),
                &RpsWorkers::new(1),
            )
            .into_iter()
            .map(|outcome| outcome.map(|result| result.map(|_| ()).map_err(|e| e.to_string())))
            .collect()
    }

    /// The record `w45` of nwin.fna (build_nq.py: a back-translated window of a profile's
    /// fragment, which hits `rps_mini`).
    fn hit_record() -> FastaRecord {
        let text = std::fs::read_to_string(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/rps/nwin.fna"),
        )
        .unwrap();
        let block = text
            .split('>')
            .find(|block| block.starts_with("w45 "))
            .unwrap();
        let letters: String = block.lines().skip(1).collect();
        record("w45", letters.as_bytes())
    }

    // B1 of the review of PS12s: a volume's refusal met by a later batch of the window does not
    // end the volume loop for the batches before it (without a valid frame), which set up every
    // volume as NCBI's loop of one batch at a time does and meet NCBI's error of a later volume's
    // missing file first (rpsblast_local.cpp:384-399: one `CLocalBlast` per volume, whose set-up
    // maps the RPS files; probe e1 of the review: NCBI rc 3).
    #[test]
    fn a_volumes_refusal_leaves_the_earlier_batches_to_meet_later_volumes() {
        let (dir, database) = rps_database("rps_mini3v", &[".01.loo"], "b1");
        let mut search = TranslatedRpsSearch {
            search: RpsSearch::default(),
            ungapped: false,
            use_sw_tback: false,
            comp_based_stats: "0".to_string(),
            several_volumes: true,
            costs: None,
        };
        search.search.x_drop_gapped_bits = Some(1e9);
        let n_only = [record("n", &[b'N'; 210])];
        let valid = [hit_record()];
        let outcome = window_errors(&search, &database, &[&n_only, &valid]);
        let Some(Err(error)) = &outcome[0] else {
            panic!("{outcome:?}");
        };
        assert!(
            error.contains("Cannot memory map RPS-BLAST database file")
                && error.contains("rps_mini3v.01.loo"),
            "{error}"
        );
        assert!(outcome[1].is_none(), "{outcome:?}");
        // With every file the X-drop's refusal is the valid batch's, after the first batch's
        // results.
        let (whole_dir, whole) = rps_database("rps_mini3v", &[], "b1w");
        let outcome = window_errors(&search, &whole, &[&n_only, &valid]);
        assert!(matches!(outcome[0], Some(Ok(()))), "{outcome:?}");
        assert!(
            matches!(&outcome[1], Some(Err(error)) if error.contains("an -xdrop_gap of 1000000000.0 bits")),
            "{outcome:?}"
        );
        let _ = std::fs::remove_dir_all(dir);
        let _ = std::fs::remove_dir_all(whole_dir);
    }

    // S1 of the review of PS12s: in a batch, a query's refusal of its cutoffs comes before an
    // earlier query's search error (base's and NCBI's order, blast_engine.c:1372-1374). `w45`'s
    // traceback reads the consensus letters, whose file is missing (NCBI's error); `p_tgg`'s frames
    // have a Lambda whose ungapped X-drop of 9e6 bits leaves 32 bits (LOSAT's refusal).
    #[test]
    fn a_batchs_cutoff_refusal_comes_before_an_earlier_querys_search_error() {
        let (dir, database) = rps_database("rps_mini", &[".psq"], "s1");
        let mut search = TranslatedRpsSearch {
            search: RpsSearch::default(),
            ungapped: false,
            use_sw_tback: false,
            comp_based_stats: "0".to_string(),
            several_volumes: false,
            costs: None,
        };
        search.search.x_drop_ungapped_bits = Some(9e6);
        let repeat = b"TGGTGTTGC".repeat(20);
        let both = [hit_record(), record("p_tgg", &repeat)];
        let outcome = window_errors(&search, &database, &[&both]);
        assert!(
            matches!(&outcome[0], Some(Err(error)) if error.contains("an -xdrop_ungap of 9000000.0 bits")),
            "{outcome:?}"
        );
        // In batches of their own the hit record's batch fails first, with NCBI's error.
        let hit = [hit_record()];
        let tgg = [record("p_tgg", &repeat)];
        let outcome = window_errors(&search, &database, &[&hit, &tgg]);
        assert!(
            matches!(&outcome[0], Some(Err(error)) if error.contains("rps_mini.psq")),
            "{outcome:?}"
        );
        assert!(outcome[1].is_none(), "{outcome:?}");
        let _ = std::fs::remove_dir_all(dir);
    }
}
