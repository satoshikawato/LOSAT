//! RPSTBLASTN's search of BLASTX's query batches (series PR, item PS12a, DECISIONS 38 D1): the six
//! frames of each nucleotide query, set up as BLASTX sets them up (`prepare_unit_queries`: frames,
//! strands, genetic codes, masks, `-query_loc`; `score_block` and `effective_lengths`: the frames'
//! Karlin-Altschul blocks and search spaces), are the "subjects" of the reversed RPS search
//! (`rps_search::translated`); its final HSPs become BLASTX's results (`blastx_lists`), which
//! BLASTX's reports write. The search takes the place of BLASTX's own (`RpsTblastnSearch`).
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/rpstblastn_app.cpp:191-209

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
    batch_gapped_k, RpsFrame, RpsQuery, RpsRun, RpsSearch, RpsWorkers, NUM_FRAMES,
};
use crate::algorithm::rps_search::volumes::{
    mask_unsupported_letters, merge_volume, report_ordinal_ids, shift_ordinal_ids, Consensus,
    RpsUnits, RpsVolumeError, RpsVolumeSet, RpsVolumeSetup,
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
    /// Whether the database has several volumes, whose results NCBI merges
    /// (`RpsTblastnSearch::merges_volumes`); set by `for_database`.
    pub(crate) several_volumes: bool,
    /// The volumes' matrix and gap costs as their set-up read them (`RpsVolumeSetup::costs`), for
    /// the requests of the archive and of the search strategy (`handle_sets`); set by `for_database`
    /// of a run that writes one.
    pub(crate) costs: Option<(String, i32, i32)>,
    /// The database length that `x_AdjustDbSize` writes into the options handle
    /// (`RpsUnits::strategy_db_length`), for the same requests; set with `costs`.
    pub(crate) strategy_db_length: Option<i64>,
    /// The search units of the database and their statistics (`RpsUnits`), set by
    /// `for_database`.
    pub(crate) units: Option<RpsUnits>,
    /// Whether the run writes a BLAST archive or a search strategy, the readers of each batch's
    /// `EffectiveSearchSpace` option (`WindowBatch::request_space`); set by `for_database`.
    pub(crate) strategy: bool,
    /// The set-up that checked the volumes in `for_database` (series PR, item PS12m): the run
    /// maps them again for each batch with it (`RpsVolumeSetup::reopen`), which knows the volumes
    /// whose cells were regrouped; set by `for_database`.
    pub(crate) setup: Option<Arc<RpsVolumeSetup>>,
}

impl RpsTblastnSearch for TranslatedRpsSearch {
    fn for_database(
        &self,
        options: &ResolvedOptions,
        database: Option<&BlastxDatabase>,
        requests: bool,
    ) -> Result<Arc<dyn RpsTblastnSearch>> {
        let Some(database) = database else {
            bail!("a search without -db is not supported by LOSAT's {PROGRAM}");
        };
        // The search units and their statistics (one owner with rpsblast's, line RG D1a-4),
        // whose refusals come at each unit's search of a batch (`search_window`).
        let units = RpsUnits::of(
            database.x_seqdb(),
            options.effective_search_space,
            options.db_length,
            PROGRAM,
        )?;
        let strategy_db_length = requests.then(|| units.strategy_db_length()).flatten();
        // LOSAT's refusals of the volumes' files come before the first byte, as rpsblast's do;
        // NCBI's own errors of the files stay where NCBI's search of a batch raises them
        // (`search_window`), and the options' refusals below wait for the files to open.
        let volumes = RpsVolumeSet::with_order(database.x_seqdb(), units.order().to_vec());
        let mut setup =
            RpsVolumeSetup::new(PROGRAM).with_database(database.name(), volumes.volumes.len());
        let mut faulty_files = false;
        for &index in &volumes.order {
            if setup.check(&volumes.volumes[index], false)?.is_err() {
                faulty_files = true;
            }
        }
        if !faulty_files {
            self.check_options()?;
        }
        // The volumes' costs are read for the requests of the archive and of the search strategy only (`handle_sets`;
        // after volume files that failed NCBI's checks, `costs` opens the first volume's again).
        let costs = requests.then(|| setup.costs(volumes.first())).flatten();
        Ok(Arc::new(Self {
            several_volumes: volumes.order.len() > 1,
            costs,
            strategy_db_length,
            units: Some(units),
            strategy: requests,
            setup: Some(Arc::new(setup)),
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
        take: &mut dyn FnMut(&mut dyn Write, &[FastaRecord], usize) -> Result<()>,
    ) -> Result<usize> {
        let Some(database) = database else {
            bail!("a search without -db is not supported by LOSAT's {PROGRAM}");
        };
        // The queries of the batches replayed so far: each batch's search strategy request and
        // render take the number before it.
        let processed = std::cell::Cell::new(0usize);
        let diagnostics = std::cell::RefCell::new(diagnostics);
        let (Some(units), Some(setup)) = (self.units.as_ref(), self.setup.as_deref()) else {
            bail!("LOSAT internal error: {PROGRAM}'s query batches are searched before its database is set up");
        };
        let volumes = RpsVolumeSet::with_order(database.x_seqdb(), units.order().to_vec());
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
            // before the batch's search (`SaveSearchStrategy`, cited in the module's header):
            // after the batches before it are reported (their search spaces, which the request
            // holds, `StrategyExport::note_searched`), and none for a batch after one that fails.
            // NCBI's posts of the bytes replaced in it go to the diagnostics there, after the
            // batch's reader messages (its `Warn` events, replayed before it; K4b-2 of line RG).
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
                    &batches, options, database, &volumes, units, setup, pool, &workers,
                );
                replay(
                    &events,
                    searched,
                    &mut |bytes| Ok(diagnostics.borrow_mut().write_all(bytes)?),
                    &mut |queries| take(&mut **diagnostics.borrow_mut(), queries, processed.get()),
                    &mut |queries, results| {
                        // The window's searches are done: the pool is free to format the batch's
                        // reports.
                        render(
                            &mut **diagnostics.borrow_mut(),
                            queries,
                            &mut vec![results],
                            processed.get(),
                            Some(pool),
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

    fn handle_sets(&self) -> Vec<crate::report::archive::options::Param> {
        crate::algorithm::rps_search::volumes::search_sets(
            self.strategy_db_length,
            self.costs.clone(),
        )
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
    /// spaces, the volume's refusals, then the first query in order that fails; after the
    /// volumes, LOSAT's refusal of a profile that the reports would not find,
    /// `report_ordinal_ids`); the batches after the first that fails are not searched (`None`).
    /// Each volume is a search unit of `units` (`RpsUnits`), whose statistics its search reads:
    /// the frames' search spaces, the Gumbel block and the shortest length of the hit cutoff.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:158-164
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:176-181
    #[allow(clippy::too_many_arguments)]
    fn search_window(
        &self,
        batches: &[&[FastaRecord]],
        options: &ResolvedOptions,
        database: &BlastxDatabase,
        volumes: &RpsVolumeSet<'_>,
        units: &RpsUnits,
        setup: &RpsVolumeSetup,
        pool: &crate::utils::threading::SearchPool<'_>,
        workers: &RpsWorkers,
    ) -> Vec<Option<Result<BatchResults>>> {
        let mut outcome: Vec<Option<Result<BatchResults>>> = batches.iter().map(|_| None).collect();
        if batches.is_empty() {
            return outcome;
        }
        // The composition-based statistics mode (`-comp_based_stats`; 1 by default), and with
        // it BLOSUM62's score block for the redo's matrix information (rpsblast's,
        // `rps_hooks::rps_search`), before any batch.
        let cbs = self.search.requested_cbs;
        let gapped = gapped_karlin(options);
        let sbp = if cbs > 0 {
            match crate::algorithm::rps_search::redo::redo_score_block(
                options.scoring,
                crate::core::psi::KarlinBlk {
                    lambda: gapped.lambda,
                    k: gapped.k,
                    h: gapped.h,
                },
            ) {
                Ok(sbp) => Some(sbp),
                Err(error) => {
                    outcome[0] = Some(Err(error));
                    return outcome;
                }
            }
        } else {
            None
        };
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
        let x_drops = self.search.x_drop_bits();
        // Volume after volume, each mapped for the window's search alone, as NCBI's
        // `CLocalRPSBlast` searches each query batch.
        for ((position, &index), unit) in volumes.order.iter().enumerate().zip(units.units()) {
            if states.is_empty() {
                break;
            }
            let volume = &volumes.volumes[index];
            // The unit's Gumbel block and `avg_subject_length` of the hit cutoff (the shortest
            // length, as rpsblast's `cutoff_subject_length`).
            let gumbel = gumbel_block(options, units.gumbel_length(unit));
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:932-940
            // (The stretch is NCBI's only inside the Gumbel block's branch: blast_parameters.c:922-945.)
            let prelim_evalue = if cbs > 1 && gumbel.is_some() {
                5.0 * options.evalue
            } else {
                options.evalue
            };
            let cutoff_subject_length = unit.min_length;
            let (db_length, db_count) = units.eff_lengths(unit);
            // The volume's files: the same for every batch, so their error is the first batch's.
            // With composition-based statistics the volume's frequency ratios (`.freq`) open too.
            let files = match setup.reopen(volume, cbs > 0) {
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
            // The preliminary stage's lengths where they are not the setup's
            // (`RpsUnits::prelim_lengths`: one volume whose total is negative), and its spaces of
            // every frame of every searching batch's queries, in their order
            // (`RpsRun::prelim_spaces`).
            let prelim_lengths = units.prelim_lengths(unit, options.effective_search_space);
            let mut prelim_spaces: Option<Vec<i64>> = prelim_lengths.map(|_| Vec::new());
            for batch in 0..states.len() {
                let WindowBatch {
                    prepared,
                    parameters,
                    reported,
                    lengths: held,
                    ..
                } = &mut states[batch];
                let parameters = parameters.get_or_insert_with(|| score_block(prepared));
                if !parameters.iter().any(|p| p.valid) {
                    continue;
                }
                // The unit's own search spaces (`RpsUnits::eff_lengths`), computed where the
                // parameters hold another unit's lengths (in the common case every unit's are
                // the first's).
                let lengths = if *held == Some((db_length, db_count)) {
                    Ok(())
                } else {
                    // The merged result keeps the first unit's statistics (in path order), which
                    // the reports print (`RpsVolumeSet::first`, as rpsblast's footers): kept
                    // before a later unit's replace them.
                    if position > 0 && reported.is_none() {
                        *reported = Some(parameters.clone());
                    }
                    effective_lengths(prepared, parameters, options, db_length, db_count, None)
                        .inspect(|()| *held = Some((db_length, db_count)))
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
                // The frames' search spaces of the preliminary stage (`prelim_spaces`).
                if let (Some((length, count)), Some(spaces)) =
                    (prelim_lengths, prelim_spaces.as_mut())
                {
                    let mut prelim = parameters.clone();
                    if let Err(error) =
                        effective_lengths(prepared, &mut prelim, options, length, count, None)
                    {
                        fail(
                            &mut states,
                            &mut outcome,
                            batch,
                            RpsBatchError::Setup(error).into(),
                        );
                        break;
                    }
                    spaces.extend(prelim.iter().map(|context| context.search_space));
                }
                // The rest of the unit's setup (`BLAST_GapAlignSetUp` after its effective lengths)
                // and its engine's start, the same for every batch: the first batch that searches
                // the unit meets them. A count of 0 without the Gumbel block and a concatenated
                // length NCBI cannot search are LOSAT's refusals (`RpsUnit::check_setup_count`,
                // `check_engine_length`), an average below one letter NCBI's error (review A of
                // D1a-4, MF5; review B, MF1: after the batch's query set-up).
                let gumbel_block = options.has_gumbel_block();
                let unit_checked = unit
                    .check_setup_count(gumbel_block, PROGRAM)
                    .and_then(|()| match unit.setup_subject_length(gumbel_block) {
                        Some(_) => Ok(()),
                        None => Err(RpsBatchError::Setup(anyhow::anyhow!(
                            crate::stats::setup::SUBJECT_LENGTH_TOO_SHORT
                        ))
                        .into()),
                    })
                    .and_then(|()| unit.check_engine_length(options.window_size, PROGRAM));
                if let Err(error) = unit_checked {
                    fail(&mut states, &mut outcome, batch, error);
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
                    !crate::blastinput::app::ReportFormat::from_number(i32::from(options.outfmt))
                        .is_some_and(crate::blastinput::app::ReportFormat::writes_seqaligns),
                ));
                spans.push(queries.len() - start);
            }
            let run = RpsRun {
                queries,
                gapped,
                gumbel: gumbel.as_ref(),
                cutoff_subject_length,
                evalue: options.evalue,
                prelim_evalue,
                cbs,
                sbp: sbp.as_ref(),
                hitlist_size: options.hitlist_size.max(1) as usize,
                prelim_hitlist_size: crate::algorithm::blastx::results::prelim_hitlist_size(
                    options,
                ),
                window: options.window_size,
                // The unit's concatenated length `Lc` (an alias file's `LENGTH` and `NSEQ` over
                // one volume).
                concat_length: unit.concat_length(),
                sorts_results: unit.sorts_results,
                prelim_spaces: prelim_spaces.as_deref(),
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
                let k_marks: Vec<Option<u32>> = lists.iter().map(|traced| traced.k_mark).collect();
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
                    // The K of the footers: volume 0's marks and `.aux` Ks, as rpsblast's
                    // (`thread_results[0]`).
                    if cbs > 0 {
                        state.gapped_k = Some(footer_gapped_k(
                            &k_marks,
                            &files.aux.karlin_k,
                            gapped.k,
                            state
                                .parameters
                                .as_deref()
                                .is_some_and(|parameters| parameters[0].valid),
                        ));
                    }
                } else {
                    merge_volume(&mut state.hit_lists, lists);
                }
                state.searched = true;
            }
            // The volume's files close after the batches' footer K, read from `files.aux` in place.
            drop(files);
        }
        // The reports look the profiles up by ordinal id (`report_ordinal_ids`), batch after batch
        // as NCBI's loop reports them: a profile that they would not find among those they report
        // (the first `-max_target_seqs` of each query's merged list, `ResolvedOptions::
        // shown_lists`) ends the run at its batch; one past that cut leaves the list, its HSPs
        // counted (`BatchResults::unlisted`). Then each batch's `EffectiveSearchSpace` option that
        // its search leaves in the options handle (`WindowBatch::request_space`).
        let db = database.x_seqdb();
        for batch in 0..states.len() {
            let state = &mut states[batch];
            let checked = report_ordinal_ids(
                db,
                &mut state.hit_lists,
                options.shown_lists(),
                PROGRAM,
                units.searches_volumes_bare(),
            )
            .and_then(|unlisted| {
                state.unlisted = unlisted.unwrap_or_default();
                if self.strategy {
                    state.request_space(units, options, &gapped)
                } else {
                    Ok(None)
                }
            });
            match checked {
                Ok(request_space) => state.request_space = request_space,
                Err(error) => {
                    fail(&mut states, &mut outcome, batch, error);
                    break;
                }
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
    /// The first unit's (the merged result's statistics), kept where a later unit's effective
    /// lengths replace them in `parameters`.
    reported: Option<Vec<ContextParameters>>,
    /// Whether a volume searched the batch.
    searched: bool,
    hit_lists: Vec<Option<BlastpHitList>>,
    /// The `EffectiveSearchSpace` option that the batch's search leaves (`request_space`).
    request_space: Option<i64>,
    /// The effective lengths (`RpsUnits::eff_lengths`) that `parameters`' search spaces were
    /// computed with.
    lengths: Option<(i64, i32)>,
    /// The HSPs of each query's profiles that left its merged list past the formatter's cut
    /// (`report_ordinal_ids`; `BatchResults::unlisted`).
    unlisted: Vec<usize>,
    /// The queries' footer K with composition-based statistics (`footer_gapped_k`).
    gapped_k: Option<Vec<f64>>,
}

impl WindowBatch {
    fn new(prepared: PreparedQueryBatch, queries: usize) -> Self {
        Self {
            prepared,
            parameters: None,
            reported: None,
            searched: false,
            hit_lists: vec![None; queries],
            request_space: None,
            lengths: None,
            unlisted: Vec::new(),
            gapped_k: None,
        }
    }

    /// The `EffectiveSearchSpace` option that the search of the batch leaves in NCBI's options
    /// handle, which the search strategy's next request and the BLAST archive hold: `-searchsp`
    /// as given, else the last unit's value for its context 0 (`RpsUnits::request_space`, with the
    /// context's Karlin block as its search space reads it); `None` for a batch that no volume
    /// searched.
    fn request_space(
        &self,
        units: &RpsUnits,
        options: &ResolvedOptions,
        gapped: &KarlinParams,
    ) -> Result<Option<i64>> {
        if !self.searched {
            return Ok(None);
        }
        if options.effective_search_space != 0 {
            return Ok(Some(options.effective_search_space));
        }
        let Some(context0) = self.parameters.as_deref().and_then(<[_]>::first) else {
            return Ok(None);
        };
        let length = if context0.valid {
            self.prepared.contexts[0].length as i64
        } else {
            0
        };
        let ka = if options.gapped {
            *gapped
        } else {
            context0.ungapped
        };
        units.request_space(length, &ka, PROGRAM).map(Some)
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
            request_space: self.request_space,
            unlisted: self.unlisted,
            shown: Vec::new(),
            gapped_k: self.gapped_k,
        }
    }
}

/// The gapped K of the footers of a batch's queries with composition-based statistics: the
/// redo of each list writes `kbp_gap[0]`, and the stream of a translated query's lists reads
/// them by ordinal id (`stream_lists`), so the batch's highest mark leaves it
/// (`RpsTraced::k_mark`, `batch_gapped_k` of the batch). Only a footer that reads `kbp_gap[0]`
/// shows it: the first query's, where its first valid context is 0 (frame +1,
/// `first_frame_valid`); every other footer shows the table's K `table_k`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:82-104
fn footer_gapped_k(
    k_marks: &[Option<u32>],
    karlin_k: &[f64],
    table_k: f64,
    first_frame_valid: bool,
) -> Vec<f64> {
    if !first_frame_valid {
        return vec![table_k; k_marks.len()];
    }
    batch_gapped_k(&[0..k_marks.len()], k_marks, karlin_k, table_k, false)
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
        request_space: None,
        unlisted: Vec::new(),
        shown: Vec::new(),
        gapped_k: None,
    }
}

/// The frames of each query of a prepared batch as the RPS engine reads them, with the
/// statistics of their contexts (`score_block`, `effective_lengths`): BLASTX's layout of
/// `sequence_start` holds each frame with letters between two sentinels at its context's offset;
/// `masked` holds a frame's letters where a 26-column database masks some
/// (`mask_unsupported_letters`). `-qcov_hsp_perc` reads the batch's first query's frames.
/// `recount_identities` is false for the Seq-aligns as they are (outfmt 8, 9, 11, 12:
/// `ReportFormat::writes_seqaligns`), which show the traceback's counts on the searched letters,
/// not the formatter's counts on the letters as read (`traceback::formatter_identities`).
pub(crate) fn frame_queries<'a>(
    prepared: &'a PreparedQueryBatch,
    parameters: &[ContextParameters],
    masked: &'a [Option<Vec<u8>>],
    recount_identities: bool,
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
            aa_seq_nomask: (recount_identities && nomask != aa_seq).then_some(nomask),
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
/// the composition-based method are the traceback's. The ends stay the HSP's, which NCBI's sort
/// before the Seq-aligns reads; the reports print the Seq-align's (`Hsp::seqalign_ends`).
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
#[path = "search_tests.rs"]
mod tests;
