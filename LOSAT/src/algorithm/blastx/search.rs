//! Internal native preliminary boundary. Linking/Kappa/traceback/report remain separate.
use super::{
    args::ResolvedOptions,
    input::FastaRecord,
    parameters::{
        effective_lengths, ka_once_on, score_block_for, subject_parameters, ContextParameters,
    },
    preliminary::{purge_endpoints, sort_by_score, ungapped, PreliminaryHsp},
    query_setup::{unit_ranges, PreparedQueryBatch},
    reader_batch::{prepare_reader_batch, UnitPrep},
    seed::{BxLookup, Diagonals},
    split::{split_queries, split_sizes, QueryChunk},
    subjects::Subjects,
};
use anyhow::{ensure, Result};

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:235-271
pub(crate) trait PreliminarySink {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:176-181,202-210
    fn unsearched_batch(
        &mut self,
        _ordinal: usize,
        _batch: &PreparedQueryBatch,
        _params: &[ContextParameters],
    ) -> Result<()> {
        Ok(())
    }
    fn batch_start(
        &mut self,
        _ordinal: usize,
        _batch: &PreparedQueryBatch,
        _params: &[ContextParameters],
    ) -> Result<()> {
        Ok(())
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3916-3921
    fn containment(&mut self, _h: &PreliminaryHsp, _contained: bool) {}
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3917-3921
    // NCBI tests each candidate HSP for containment in the interval tree. The containment callback
    // only reports that test to an observer. This method says whether the callback reads its
    // arguments, so that the search can skip the copy of every candidate HSP kept for it. No value
    // NCBI computes changes.
    /// False when `containment` ignores its arguments, so that the search need not keep a copy of
    /// every candidate HSP for it.
    fn x_wants_containment(&self) -> bool {
        true
    }
    fn chunk_start(&mut self, _chunk: &QueryChunk) -> Result<()> {
        Ok(())
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1450-1458
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:550-552
    fn preliminary_purge(
        &mut self,
        _oid: usize,
        _before: &[PreliminaryHsp],
        _after: &[PreliminaryHsp],
    ) -> Result<()> {
        Ok(())
    }
    fn subject_start(
        &mut self,
        _oid: usize,
        _batch: &PreparedQueryBatch,
        _params: &[ContextParameters],
        _prepared_link: Option<crate::core::link_hsps::LinkCutoffs>,
    ) -> Result<()> {
        Ok(())
    }
    fn subject(
        &mut self,
        _oid: usize,
        _hsps: &[PreliminaryHsp],
        _batch: &PreparedQueryBatch,
        _params: &[ContextParameters],
    ) -> Result<()> {
        Ok(())
    }
    fn chunk_end(&mut self, _chunk: &QueryChunk, _is_split: bool) -> Result<()> {
        Ok(())
    }
    fn batch_end(
        &mut self,
        _batch: &PreparedQueryBatch,
        _params: &[ContextParameters],
    ) -> Result<()> {
        Ok(())
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:491-505
// The sink of `search_preliminary` (the tests'): it keeps nothing.
#[cfg(test)]
struct Discard;
#[cfg(test)]
impl PreliminarySink for Discard {}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:219-246
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1428-1443
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:491-555
// Only the tests call it (not `pub`: a subject's `lowercase_masks` past its letters would give the
// scan ranges outside the subject).
#[cfg(test)]
pub(crate) fn search_preliminary(
    records: &[FastaRecord],
    subjects: &[FastaRecord],
    options: &ResolvedOptions,
) -> Result<()> {
    crate::utils::threading::with_search_pool(options.num_threads as usize, "BLASTX", |pool| {
        search_core(records, subjects.into(), options, &mut Discard, pool, None)
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:491-526
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1413
// NCBI scans the subjects one after another in this loop. A pool of 2+ threads scans this many
// consecutive subjects per worker in a wave (on its own diagonal table). Only the unit of parallel
// work changes.
/// Subjects per worker in one wave (also the default size of a task of the subject flow,
/// `x_subject_flow.rs`).
const X_BX_CHUNK: usize = 64;

// EXPERIMENT (LOSAT_X_BXSUBJFLOW): the subject loop as a flow of tasks (a child module, so that
// it uses compute_subject and SubjectComputation).
#[cfg(feature = "parallel")]
#[path = "x_subject_flow.rs"]
pub(crate) mod x_subject_flow;

// LOSAT_X_BXDIAGKEEP, LOSAT_X_BXUNITCAP: the workers' diagonal tables kept one per thread, and at
// most K units searched at a time (a child module, so that it uses compute_subject and
// SubjectComputation).
#[cfg(feature = "parallel")]
#[path = "x_bx_mem.rs"]
pub(crate) mod x_bx_mem;

// No NCBI counterpart (LOSAT's check of a `-db` stream, `XStreamDb::x_check`): the point where a
// database refused before the search fails the first batch with a valid context
// (`Subjects::Unreadable`, below) is marked, and an error that the call raises after it, while the
// stream's reads have not checked every sequence, is an `XPastHead`: the run reports the
// database's first error before it (`native::search_error`), as when the whole database was read
// before the first batch.
pub(crate) fn search_core(
    records: &[FastaRecord],
    subjects: Subjects<'_>,
    options: &ResolvedOptions,
    sink: &mut impl PreliminarySink,
    pool: &crate::utils::threading::SearchPool<'_>,
    prep: Option<UnitPrep>,
) -> Result<()> {
    let past_head = std::cell::Cell::new(false);
    x_search_core(records, subjects, options, sink, pool, prep, &past_head).map_err(|error| {
        if past_head.get() {
            subjects.x_past_head(error)
        } else {
            error
        }
    })
}

/// `search_core`, which sets `past_head` where the database refused before the search fails.
#[allow(clippy::too_many_arguments)]
fn x_search_core(
    records: &[FastaRecord],
    subjects: Subjects<'_>,
    options: &ResolvedOptions,
    sink: &mut impl PreliminarySink,
    pool: &crate::utils::threading::SearchPool<'_>,
    prep: Option<UnitPrep>,
    past_head: &std::cell::Cell<bool>,
) -> Result<()> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_multiseq.cpp:175-181
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:53-54
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:246-250
    // LOSAT does not split a subject; a batch that is searched refuses longer ones (below).
    let longest_subject = subjects.longest_length();
    let numbers = subjects.numbers()?;
    let (db_length, db_count) = (numbers.total, numbers.count);
    // No NCBI counterpart: decides whether candidate HSPs are copied for the containment callback;
    // it does not change any value NCBI computes.
    let x_containment = sink.x_wants_containment();
    use super::x_timing::Batch as XBatch;
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:138-150
    // NCBI's query splitter reads its sizes (`SBlastSetupData`, `CQuerySplitter`) before the
    // batch's queries are set up, so their errors come before the setup's (an empty query's
    // `Sequence contains no data`), also for a batch that is not searched. `records` are one
    // query batch, or one prepared unit of it (`prep`, whose shape counts the batch's queries).
    let sizes = split_sizes(
        prep.as_ref()
            .and_then(|prep| prep.prepared.reader_shape.as_deref())
            .map_or(records.len(), |shape| shape.query_count),
    )?;
    let ranges = unit_ranges(records, options, subjects.seqdb());
    // Dispatch point (LOSAT_X_BXCAPFILL, LOSAT_X_BXBUILDPAR, in the strict set): a reader batch
    // searched whole is left the threads of its places (x_bx_mem::whole_batch_threads) for its
    // contexts' blocks, its lookup table's build and its subjects. `None` for every other unit.
    #[cfg(feature = "parallel")]
    let x_whole = x_bx_mem::whole_batch_threads(
        pool,
        options,
        subjects.len(),
        super::query_setup::unit_places(records, options, subjects.seqdb()),
    );
    #[cfg(not(feature = "parallel"))]
    let x_whole: Option<usize> = None;
    // The units of a reader batch (`records`) are prepared first, together, and the values NCBI
    // takes of the whole query batch are reduced (`reader_batch.rs`); each unit is then searched
    // from its preparation. A caller that searches the units one by one has prepared them
    // already (`prep`, the only unit of `records`). A batch of one unit is its own batch.
    let mut preps: std::collections::VecDeque<UnitPrep> = match prep {
        Some(prep) => {
            ensure!(ranges.len() == 1, "a prepared unit is searched alone");
            std::collections::VecDeque::from([prep])
        }
        None if ranges.len() > 1 => {
            prepare_reader_batch(records, &ranges, options, subjects, pool)?.into()
        }
        None => std::collections::VecDeque::new(),
    };
    for range in ranges {
        // Diagnostic (LOSAT_TIMING): the clocks of the batch's stages (x_timing.rs). No NCBI
        // counterpart.
        let x_t_batch = super::x_timing::lap();
        let x_t = x_t_batch;
        let queries = &records[range.clone()];
        let (mut full, mut full_params) = match preps.pop_front() {
            Some(UnitPrep {
                prepared, params, ..
            }) => (prepared, params),
            None => {
                let UnitPrep {
                    prepared, params, ..
                } = super::reader_batch::prepare_lone_unit(queries, options, x_whole.is_some())?;
                (prepared, params)
            }
        };
        super::x_timing::batch(XBatch::QuerySetup, x_t);
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:176-181,202-210
        if !full_params.iter().any(|p| p.valid) {
            for context in &mut full.contexts {
                context.is_valid = false;
            }
            sink.unsearched_batch(range.start, &full, &full_params)?;
            continue;
        }
        // NCBI reads the database's sequences in the preliminary search of the first batch
        // with a valid context (batches without one are not searched, above). The subjects are
        // checked there too, after the split sizes and the setup of the queries, whose errors
        // come first at any thread count (`native.rs`, `XUnit::take_prep`).
        // A database may include no sequence (a list or an alias file's range that keeps none):
        // its batches are searched without subjects.
        ensure!(
            !subjects.is_empty() || !matches!(subjects, Subjects::Fasta(_)),
            "BLASTX preliminary requires a protein subject"
        );
        // The refusal of a database sequence (`database::too_long`): the `-subject` record is
        // named by its number in the file (`Subject_N`).
        if longest_subject > super::database::MAX_DBSEQ_LEN {
            return Err(match &subjects {
                Subjects::Fasta(records) => {
                    let (index, length) = records
                        .iter()
                        .map(|record| record.sequence.len())
                        .enumerate()
                        .find(|(_, length)| *length > super::database::MAX_DBSEQ_LEN)
                        .expect("the longest subject is above the limit");
                    super::database::too_long(
                        format_args!("subject sequence {}", index + 1),
                        length,
                    )
                }
                _ => super::database::too_long("a database sequence", longest_subject),
            });
        }
        if subjects.is_unreadable() {
            return Err(super::subjects::UnreadableSubjects.into());
        }
        past_head.set(true);
        let x_t = super::x_timing::lap();
        // Dispatch point (LOSAT_X_BXKAONCE, in the strict set), taking the blocks: `score_block`'s
        // raw blocks, before the lengths, kept for the one chunk of a batch that is not split (the
        // chunk loop below). No NCBI counterpart.
        let mut x_raw_params = ka_once_on().then(|| full_params.clone());
        effective_lengths(&full, &mut full_params, options, db_length, db_count, None)?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2781-2787
        for (context, p) in full.contexts.iter_mut().zip(&full_params) {
            context.is_valid = p.valid;
        }
        // The options' spaces that each chunk's setup reads (`SplitQuery_SetEffectiveSearchSpace`):
        // the batch's own, unless the statistics' total before the fallback to OID 0 is not
        // positive (`parameters::preliminary_spaces`).
        let spaces: Vec<_> = match super::parameters::preliminary_spaces(
            &full,
            &full_params,
            options,
            numbers,
            matches!(subjects, Subjects::Stream(_)),
        )? {
            Some(spaces) => spaces,
            None => full_params.iter().map(|p| p.search_space).collect(),
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:219-222
        sink.batch_start(range.start, &full, &full_params)?;
        let mut chunks = split_queries(queries, &full, options, sizes)?;
        // The statistics of each subject read the options' spaces of the whole query by each
        // chunk's context index (`s_GetEffectiveSearchSpaceForContext`, blast_setup.c:676-697).
        if subjects.per_subject() {
            let option_spaces: std::sync::Arc<[i64]> = spaces.clone().into();
            for chunk in &mut chunks {
                chunk.prepared.subject_option_spaces = Some(option_spaces.clone());
            }
        }
        super::x_timing::batch(XBatch::Setup, x_t);
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:233-236
        let is_split = chunks.len() > 1;
        // With a pool of 2+ threads: one job per query chunk instead of one per
        // subject.  NCBI c++/src/algo/blast/api/prelim_stage.cpp:235-271 searches
        // the chunks one after another and only `BlastHSPStreamMerge` relates
        // them; here a wave of chunks is searched on the pool and then written
        // to the sink in chunk order, subject by subject, as before.
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:237-271
        // NCBI searches the query chunks one after another and relates them only in
        // BlastHSPStreamMerge. With a pool of 2+ threads a wave of chunks is searched on the pool
        // (x_compute_chunk) and written to the sink in chunk order below. Each chunk is searched
        // exactly as before.
        // Dispatch point (LOSAT_X_BXCHUNKFLOW, in the strict set): where the subject flow would
        // run and there are enough subjects for the pool (x_subject_flow::chunks_in_order), the
        // chunks stay in the serial chunk loop below, one after another, and the flow searches
        // each chunk's subjects on all threads.
        let x_t = super::x_timing::lap();
        #[cfg(feature = "parallel")]
        let x_chunks_par = super::runtime::x_bx_par_rule(|| {
            pool.enabled()
                && chunks.len() > 1
                && !x_subject_flow::chunks_in_order(pool, options, subjects)
        });
        #[cfg(feature = "parallel")]
        let chunks = if x_chunks_par {
            use rayon::prelude::*;
            crate::utils::threading::report_stage("blastx", "query_chunks", chunks.len(), true);
            // Dispatch point (LOSAT_X_BXUNITCAP, in the strict set): waves of at most K chunks
            // (fewer inside the capped batch flow, x_bx_mem::chunk_wave), each chunk's subjects as
            // tasks on the pool (x_bx_mem.rs).
            let x_cap = x_bx_mem::chunk_unit_cap(pool, options, subjects);
            // Dispatch point (LOSAT_X_BXBUILDPAR, in the strict set): a capped wave's chunks build
            // their tables on the threads the cap leaves each (x_bx_mem.rs).
            let x_parts =
                x_bx_mem::build_parts(pool, options, x_cap.is_some(), subjects.len(), false, None);
            let mut x_next = 0;
            while x_next < chunks.len() {
                let x_wave = x_cap.map(|(k, _)| x_bx_mem::chunk_wave(pool, k));
                let wave = x_wave
                    .as_ref()
                    .map_or(pool.threads().max(2), x_bx_mem::XChunkWave::size);
                let wave_chunks = &chunks[x_next..chunks.len().min(x_next + wave)];
                x_next += wave_chunks.len();
                let computed: Vec<_> = pool.install(|| {
                    wave_chunks
                        .par_iter()
                        .with_max_len(1)
                        .map(|chunk| {
                            x_compute_chunk(
                                chunk,
                                subjects,
                                options,
                                &spaces,
                                numbers,
                                x_containment,
                                x_cap.map(|(_, per_task)| per_task),
                                x_parts,
                            )
                        })
                        .collect()
                });
                // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:269-271
                // Replay in chunk order: the sink sees chunk_start, each subject and chunk_end in
                // the same order as the serial chunk loop below.
                for (chunk, computed) in wave_chunks.iter().zip(computed) {
                    sink.chunk_start(chunk)?;
                    let Some((params, results)) = computed? else {
                        continue;
                    };
                    let batch = &chunk.prepared;
                    let mut kept = results.into_iter().peekable();
                    for oid in 0..subjects.len() {
                        let computed = match kept.next_if(|(kept_oid, _)| *kept_oid == oid) {
                            Some((_, computed)) => computed,
                            None => SubjectComputation::empty(),
                        };
                        sink.subject_start(oid, batch, &params, computed.prepared_link)?;
                        for (h, contained) in &computed.containment {
                            sink.containment(h, *contained);
                        }
                        let SubjectComputation { raw, purged, .. } = computed;
                        if options.gapped && !raw.is_empty() {
                            sink.preliminary_purge(oid, &raw, &purged)?;
                        }
                        sink.subject(oid, &purged, batch, &params)?;
                    }
                    sink.chunk_end(chunk, is_split)?;
                }
            }
            Vec::new()
        } else {
            chunks
        };
        super::x_timing::batch(XBatch::ChunksPar, x_t);
        // Dispatch point (LOSAT_X_BXCHUNKCROSS, in the strict set): where the loop below would
        // search two or more chunks one after another, each with the subject flow, they run as one
        // flow of tasks across the chunks; each chunk's setup is made by the loop's statements
        // below, the sink sees the chunks and subjects in the loop's order (x_subject_flow.rs).
        #[cfg(feature = "parallel")]
        let chunks = match x_subject_flow::crosses_chunks(pool, options, subjects, chunks.len()) {
            Some(x_k) => {
                if x_raw_params.take().is_some() {
                    crate::utils::xreach::x_reach!("LOSAT_X_BXKAONCE", false);
                }
                let x_parts = x_bx_mem::cross_build_parts(pool);
                let fixed_spaces = &spaces;
                let x_chunks = &chunks;
                let x_make =
                    |c: usize, x_parts: usize| -> Result<Option<x_subject_flow::XChunkSetup>> {
                        let batch = &x_chunks[c].prepared;
                        let mut params = score_block_for(batch, options.scoring);
                        if !params.iter().any(|p| p.valid) {
                            return Ok(None);
                        }
                        let lookup = BxLookup::x_new_in_parts(batch, &params, options, x_parts)?;
                        let last = batch.contexts.last().expect("BLASTX contexts");
                        let length = last.offset + last.length;
                        let table = Diagonals::new(length, options.window_size)?;
                        let xdrop = subject_parameters(
                            batch,
                            &mut params,
                            options,
                            numbers,
                            Some(fixed_spaces),
                        )?;
                        Ok(Some(x_subject_flow::XChunkSetup {
                            params,
                            lookup,
                            xdrop,
                            table,
                            length,
                        }))
                    };
                let x_t = super::x_timing::lap();
                x_subject_flow::run_chunks(
                    x_k,
                    x_chunks,
                    &x_make,
                    x_parts,
                    subjects,
                    options,
                    pool.threads() - 1,
                    x_bx_mem::diag_keep(options, subjects),
                    x_containment,
                    x_bx_mem::in_capped_batch(pool),
                    is_split,
                    sink,
                )?;
                super::x_timing::batch(XBatch::ScanPar, x_t);
                Vec::new()
            }
            None => chunks,
        };
        for chunk in chunks {
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:237-244
            sink.chunk_start(&chunk)?;
            let batch = &chunk.prepared;
            // Dispatch point (LOSAT_X_BXKAONCE, in the strict set): the one chunk of a batch that
            // is not split is a clone of `full` (`split_queries`) whose contexts' `is_valid` are
            // the raw blocks' `valid` (set above); `score_block` on it would compute the same
            // compositions and so the same Karlin blocks, and `effective_lengths` writes only
            // `length_adjustment` and `search_space`, never `valid` or `ungapped`, so the raw
            // blocks kept above are its result. The chunks of a split batch have sequences of
            // their own and keep `score_block`.
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:243-246
            // NCBI computes the query's score block once in setup, and a chunk's own only for a
            // split query (the loop this is in runs only when `query_splitter->IsQuerySplit()` is true).
            let mut params = match x_raw_params.take() {
                Some(raw) if !is_split => {
                    crate::utils::xreach::x_reach!("LOSAT_X_BXKAONCE", true);
                    raw
                }
                other => {
                    if other.is_some() {
                        crate::utils::xreach::x_reach!("LOSAT_X_BXKAONCE", false);
                    }
                    score_block_for(batch, options.scoring)
                }
            };
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:277-288
            if !params.iter().any(|p| p.valid) {
                continue;
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:685-696
            let x_t = super::x_timing::lap();
            let fixed_spaces = &spaces;
            // Dispatch point (LOSAT_X_BXBUILDPAR, in the strict set): a batch of the capped batch
            // flow builds its table on the threads the cap leaves it (x_bx_mem.rs).
            #[cfg(feature = "parallel")]
            let x_parts =
                x_bx_mem::build_parts(pool, options, false, subjects.len(), true, x_whole);
            #[cfg(not(feature = "parallel"))]
            let x_parts = 1;
            let lookup = BxLookup::x_new_in_parts(batch, &params, options, x_parts)?;
            let last = batch.contexts.last().expect("BLASTX contexts");
            let xdrop =
                subject_parameters(batch, &mut params, options, numbers, Some(fixed_spaces))?;
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:172-180
            // utils::threading owns the search pool; each slot retains its own
            // diagonal table across bounded O(worker-count) reduction windows.
            // Caller/collector/Kappa state remains serial in the original OID order.
            // With a pool of 2+ threads a wave is `X_BX_CHUNK` consecutive subjects per worker,
            // and a worker scans its `X_BX_CHUNK` subjects one after another on its own diagonal
            // table, as the serial loop does on the single table. One subject per worker per wave
            // would mean a fork-join for every `threads` subjects, each shorter than the
            // fork-join itself.
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1413
            // A worker scans `X_BX_CHUNK` consecutive subjects one after another, as the
            // serial loop does on the single table. Scheduling only.
            let wave_size = if pool.enabled() {
                pool.threads() * X_BX_CHUNK
            } else {
                1
            };
            #[cfg(feature = "parallel")]
            let x_flow_k = x_subject_flow::chosen(pool, options, subjects, x_whole);
            // Dispatch point (LOSAT_X_BXDIAGKEEP, in the strict set): the flow's helpers run on the
            // tables their threads keep (x_bx_mem.rs); only where the flow searches the subjects
            // (the waves below then have nothing to do).
            #[cfg(feature = "parallel")]
            let x_keep =
                pool.enabled() && x_flow_k.is_some() && x_bx_mem::diag_keep(options, subjects);
            // The loop's own table: where the helpers run on the tables their threads keep, the
            // one this thread keeps too (x_bx_mem.rs `XDiagSlot`), not a new table for every
            // batch or chunk the thread searches.
            let mut x_new_table = None;
            #[cfg(feature = "parallel")]
            let mut x_own_slot = x_keep
                .then(|| x_bx_mem::XDiagSlot::new(last.offset + last.length, options.window_size));
            #[cfg(feature = "parallel")]
            let diagonals: &mut Diagonals = match x_own_slot.as_mut() {
                Some(slot) => slot.get()?,
                None => x_new_table.insert(Diagonals::new(
                    last.offset + last.length,
                    options.window_size,
                )?),
            };
            #[cfg(not(feature = "parallel"))]
            let diagonals: &mut Diagonals = x_new_table.insert(Diagonals::new(
                last.offset + last.length,
                options.window_size,
            )?);
            #[cfg(feature = "parallel")]
            let mut worker_diagonals = if pool.enabled() && !x_keep {
                (0..subjects.len().min(pool.threads()))
                    .map(|_| Diagonals::new(last.offset + last.length, options.window_size))
                    .collect::<Result<Vec<_>>>()?
            } else {
                Vec::new()
            };
            super::x_timing::batch(XBatch::Lookup, x_t);
            // Dispatch point (LOSAT_X_BXSUBJFLOW[=k], in the strict set: on unless `=0`): the
            // subjects as a flow of tasks of k consecutive subjects, handed to the sink in subject
            // order by this thread
            // (x_subject_flow.rs); the waves below then have nothing left to do.
            #[cfg(feature = "parallel")]
            let x_flowed = if let Some(x_k) = x_flow_k {
                let x_t = super::x_timing::lap();
                // Dispatch point (LOSAT_X_BXCAPFILL, in the strict set): a batch of one task of the
                // capped batch flow has helpers for the threads the cap leaves it, a reader batch
                // searched whole for the threads of its places (x_bx_mem::flow_threads).
                let x_threads = x_bx_mem::flow_threads(pool, options, subjects.len(), x_whole);
                let helper_tables = if x_keep {
                    (1..x_threads)
                        .map(|_| {
                            x_bx_mem::XScratch::Slot(x_bx_mem::XDiagSlot::new(
                                last.offset + last.length,
                                options.window_size,
                            ))
                        })
                        .collect()
                } else {
                    std::mem::take(&mut worker_diagonals)
                        .into_iter()
                        .take(x_threads - 1)
                        .map(x_bx_mem::XScratch::Owned)
                        .collect()
                };
                x_subject_flow::run(
                    x_k,
                    subjects,
                    batch,
                    &params,
                    options,
                    &lookup,
                    &mut *diagonals,
                    helper_tables,
                    xdrop,
                    x_containment,
                    x_bx_mem::in_capped_batch(pool) || x_whole.is_some(),
                    sink,
                )?;
                super::x_timing::batch(XBatch::ScanPar, x_t);
                true
            } else {
                false
            };
            #[cfg(not(feature = "parallel"))]
            let x_flowed = false;
            // The flow (LOSAT_X_BXSUBJFLOW) leaves the waves, and so the chunks, nothing to do.
            if !x_flowed {
                crate::utils::xreach::x_reach!("LOSAT_X_BXCHUNK", pool.enabled());
            }
            let x_waves = if x_flowed {
                0
            } else {
                subjects.len().div_ceil(wave_size)
            };
            // The blocks of the serial pass.
            let mut x_pass = subjects.x_pass();
            for wave_index in 0..x_waves {
                let wave_start = wave_index * wave_size;
                let wave: Vec<usize> =
                    (wave_start..subjects.len().min(wave_start + wave_size)).collect();
                // The subjects of a wave that the pool searches,
                // read before its fork.
                let x_wave_block = if pool.enabled() {
                    subjects.x_read_range(wave_start, wave_start + wave.len())?
                } else {
                    None
                };
                let x_wave_subjects = subjects.x_with_block(x_wave_block.as_ref());
                let x_t = super::x_timing::lap();
                // Diagnostic (LOSAT_TIMING): the wave's tasks (x_timing::SubjectWave). No NCBI
                // counterpart.
                let x_wave = super::x_timing::SubjectWave::begin();
                let mut scheduled = {
                    #[cfg(feature = "parallel")]
                    {
                        use rayon::prelude::*;
                        if pool.enabled() {
                            pool.install(|| {
                                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1413,1469-1470
                                // Each worker runs compute_subject for its `X_BX_CHUNK`
                                // subjects in order; the results are flattened in subject order.
                                let x_k = X_BX_CHUNK;
                                wave.par_chunks(x_k)
                                    .zip(worker_diagonals.par_iter_mut())
                                    .enumerate()
                                    .flat_map_iter(|(x_slot, (x_part, scratch))| {
                                        x_wave.task(|| {
                                            x_part
                                                .iter()
                                                .enumerate()
                                                .map(|(x_within, &oid)| {
                                                    debug_assert_eq!(
                                                        oid,
                                                        wave_index * wave_size
                                                            + x_slot * x_k
                                                            + x_within
                                                    );
                                                    compute_subject_any(
                                                        x_wave_subjects,
                                                        batch,
                                                        &params,
                                                        options,
                                                        &lookup,
                                                        scratch,
                                                        xdrop,
                                                        oid,
                                                        x_containment,
                                                    )
                                                })
                                                .collect::<Vec<_>>()
                                        })
                                    })
                                    .collect::<Vec<_>>()
                            })
                        } else {
                            Vec::new()
                        }
                    }
                    #[cfg(not(feature = "parallel"))]
                    {
                        Vec::<Result<SubjectComputation>>::new()
                    }
                }
                .into_iter();
                x_wave.end();
                super::x_timing::batch(XBatch::ScanPar, x_t);
                let x_t_replay = super::x_timing::lap();
                for &oid in &wave {
                    let computed = if let Some(computed) = scheduled.next() {
                        computed?
                    } else {
                        let x_t = super::x_timing::lap();
                        let x_subjects = if x_wave_block.is_some() {
                            x_wave_subjects
                        } else {
                            x_pass.view(oid)?
                        };
                        let computed = compute_subject_any(
                            x_subjects,
                            batch,
                            &params,
                            options,
                            &lookup,
                            &mut *diagonals,
                            xdrop,
                            oid,
                            x_containment,
                        );
                        super::x_timing::batch(XBatch::ScanInline, x_t);
                        computed?
                    };
                    sink.subject_start(oid, batch, &params, computed.prepared_link)?;
                    for (h, contained) in &computed.containment {
                        sink.containment(h, *contained);
                    }
                    let SubjectComputation { raw, purged, .. } = computed;
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-899
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:550-552
                    if options.gapped && !raw.is_empty() {
                        sink.preliminary_purge(oid, &raw, &purged)?;
                    }
                    sink.subject(oid, &purged, batch, &params)?;
                }
                super::x_timing::batch(XBatch::Replay, x_t_replay);
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:268-271
            let x_t = super::x_timing::lap();
            sink.chunk_end(&chunk, is_split)?;
            super::x_timing::batch(XBatch::ChunkEnd, x_t);
        }
        // full preparation retains full DNA lengths, masked/nomask bytes and ordinal.
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:292-301
        let x_t = super::x_timing::lap();
        sink.batch_end(&full, &full_params)?;
        super::x_timing::batch(XBatch::BatchEnd, x_t);
        super::x_timing::batch(XBatch::Sum, x_t_batch);
        super::x_timing::batch_count();
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:172-179
// Only independent seed/extension state belongs to workers. Existing LOSAT
// utils::threading owns pool lifetime; all sink state is replayed by input OID.
// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:243-260
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1413
// x_compute_chunk is everything the serial chunk loop does for one chunk before it writes to the
// sink: the chunk's scoring parameters, lookup table, diagonal table, subject parameters and the
// subject loop. The code is the same; it returns the results instead of writing them to the sink.
/// Everything the chunk loop computes before it talks to the sink, for one chunk; `None` when no
/// context of the chunk has valid Karlin parameters (the loop then skips the chunk). `x_par`: the
/// chunk's subjects as tasks of that many subjects on the pool (`LOSAT_X_BXUNITCAP`'s waves);
/// `x_parts`: the parts of its lookup table's build (`LOSAT_X_BXBUILDPAR`).
#[cfg(feature = "parallel")]
#[allow(clippy::type_complexity)]
fn x_compute_chunk(
    chunk: &QueryChunk,
    subjects: Subjects<'_>,
    options: &ResolvedOptions,
    fixed_spaces: &[i64],
    numbers: super::parameters::SubjectNumbers,
    x_containment: bool,
    x_par: Option<usize>,
    x_parts: usize,
) -> Result<Option<(Vec<ContextParameters>, Vec<(usize, SubjectComputation)>)>> {
    let batch = &chunk.prepared;
    let mut params = score_block_for(batch, options.scoring);
    if !params.iter().any(|p| p.valid) {
        return Ok(None);
    }
    let lookup = BxLookup::x_new_in_parts(batch, &params, options, x_parts)?;
    let last = batch.contexts.last().expect("BLASTX contexts");
    if let Some(per_task) = x_par {
        let xdrop = subject_parameters(batch, &mut params, options, numbers, Some(fixed_spaces))?;
        let results = x_bx_mem::chunk_subjects_on_pool(
            per_task,
            subjects,
            batch,
            &params,
            options,
            &lookup,
            xdrop,
            x_containment,
        )?;
        return Ok(Some((params, results)));
    }
    let mut diagonals = Diagonals::new(last.offset + last.length, options.window_size)?;
    let xdrop = subject_parameters(batch, &mut params, options, numbers, Some(fixed_spaces))?;
    // Only the subjects with something for the sink are kept, by oid: the chunks of a wave are
    // all in memory at once, and a database holds hundreds of thousands of subjects, most
    // without a hit. The sink gets the others as empty computations (`SubjectComputation::
    // empty`), so it sees the same calls as from the serial chunk loop.
    let mut results = Vec::new();
    // The blocks of the chunk's pass.
    let mut x_pass = subjects.x_pass();
    for oid in 0..subjects.len() {
        let computed = compute_subject_any(
            x_pass.view(oid)?,
            batch,
            &params,
            options,
            &lookup,
            &mut diagonals,
            xdrop,
            oid,
            x_containment,
        )?;
        if !computed.is_empty() {
            results.push((oid, computed));
        }
    }
    Ok(Some((params, results)))
}

// `-db_soft_mask` and `-db_hard_mask` (series NP, M4b): the search of a masked database sequence.
#[path = "search_db_masked.rs"]
mod db_masked;

struct SubjectComputation {
    raw: Vec<PreliminaryHsp>,
    purged: Vec<PreliminaryHsp>,
    containment: Vec<(PreliminaryHsp, bool)>,
    prepared_link: Option<crate::core::link_hsps::LinkCutoffs>,
}

/// A piece of a database sequence that a search with `-db_soft_mask` or `-db_hard_mask` scans in
/// place of the sequence (`search_db_masked.rs`): its letters with the sentinel after them (and
/// before them too for a batch with `x_ungap_sent`), and its scan ranges (`None`: every word).
#[derive(Clone, Copy)]
struct SubjectPiece<'a> {
    buffer: &'a [u8],
    ranges: Option<&'a [(i32, i32)]>,
}

#[cfg(feature = "parallel")]
impl SubjectComputation {
    /// A subject without anything for the sink (the chunks searched side by side keep only the
    /// others).
    fn empty() -> Self {
        Self {
            raw: Vec::new(),
            purged: Vec::new(),
            containment: Vec::new(),
            prepared_link: None,
        }
    }

    /// Whether the sink would get nothing from this subject but its oid.
    fn is_empty(&self) -> bool {
        self.raw.is_empty()
            && self.purged.is_empty()
            && self.containment.is_empty()
            && self.prepared_link.is_none()
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:172-179
// Only independent seed/extension state belongs to workers. Existing LOSAT
// utils::threading owns pool lifetime; all sink state is replayed by input OID.
#[allow(clippy::too_many_arguments)]
fn compute_subject(
    subjects: Subjects<'_>,
    batch: &PreparedQueryBatch,
    params: &[ContextParameters],
    options: &ResolvedOptions,
    lookup: &BxLookup,
    diagonals: &mut Diagonals,
    xdrop: i32,
    oid: usize,
    x_containment: bool,
    piece: Option<SubjectPiece<'_>>,
) -> Result<SubjectComputation> {
    compute_subject_with(
        super::seed::Blosum62Bx,
        subjects,
        batch,
        params,
        options,
        lookup,
        diagonals,
        xdrop,
        oid,
        x_containment,
        piece,
    )
}

/// [`compute_subject`] of a search with another matrix than BLOSUM62 (item SC6 of line NP), its
/// own instance of the stages; the callers choose it once per subject by the search's matrix. Not
/// inlined, so that BLOSUM62's keeps its code.
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn compute_subject_of(
    subjects: Subjects<'_>,
    batch: &PreparedQueryBatch,
    params: &[ContextParameters],
    options: &ResolvedOptions,
    lookup: &BxLookup,
    diagonals: &mut Diagonals,
    xdrop: i32,
    oid: usize,
    x_containment: bool,
    piece: Option<SubjectPiece<'_>>,
) -> Result<SubjectComputation> {
    compute_subject_with(
        super::seed::MatrixBx::of(options.scoring),
        subjects,
        batch,
        params,
        options,
        lookup,
        diagonals,
        xdrop,
        oid,
        x_containment,
        piece,
    )
}

/// [`compute_subject`] or [`compute_subject_of`] by the search's matrix, with the chunk's
/// parameters `params`, or, where NCBI computes the statistics of each subject (the batch's
/// `subject_option_spaces`, a database length of 0), with the subject's own
/// (`compute_own_subject`): one test of the batch per subject on the default path.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1428-1443 (cited above
// `search_preliminary`): the subject's own statistics before its scan.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn compute_subject_any(
    subjects: Subjects<'_>,
    batch: &PreparedQueryBatch,
    params: &[ContextParameters],
    options: &ResolvedOptions,
    lookup: &BxLookup,
    diagonals: &mut Diagonals,
    xdrop: i32,
    oid: usize,
    x_containment: bool,
) -> Result<SubjectComputation> {
    if let Some(option_spaces) = batch.subject_option_spaces.as_deref() {
        return compute_own_subject(
            subjects,
            batch,
            params,
            option_spaces,
            options,
            lookup,
            diagonals,
            xdrop,
            oid,
            x_containment,
        );
    }
    compute_subject_scored(
        subjects,
        batch,
        params,
        options,
        lookup,
        diagonals,
        xdrop,
        oid,
        x_containment,
    )
}

/// [`compute_subject`] or [`compute_subject_of`] by the search's matrix.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn compute_subject_scored(
    subjects: Subjects<'_>,
    batch: &PreparedQueryBatch,
    params: &[ContextParameters],
    options: &ResolvedOptions,
    lookup: &BxLookup,
    diagonals: &mut Diagonals,
    xdrop: i32,
    oid: usize,
    x_containment: bool,
) -> Result<SubjectComputation> {
    // `-db_soft_mask` and `-db_hard_mask` (`search_db_masked.rs`): one test of the options per
    // subject on the default path.
    if let Some(mask) =
        crate::algorithm::common::db_mask_search::SearchMask::of(options.db_mask.as_ref())
    {
        return db_masked::compute_masked_subject(
            subjects,
            batch,
            params,
            options,
            lookup,
            diagonals,
            xdrop,
            oid,
            x_containment,
            mask,
        );
    }
    if options.scoring == crate::config::ScoringMatrix::Blosum62 {
        compute_subject(
            subjects,
            batch,
            params,
            options,
            lookup,
            diagonals,
            xdrop,
            oid,
            x_containment,
            None,
        )
    } else {
        compute_subject_of(
            subjects,
            batch,
            params,
            options,
            lookup,
            diagonals,
            xdrop,
            oid,
            x_containment,
            None,
        )
    }
}

/// [`compute_subject_scored`] with the subject's own parameters where NCBI computes the
/// statistics of each subject (`parameters::subject_update`, with the whole query's options'
/// spaces `option_spaces`): a subject with residues, the one NCBI searches (blast_engine.c:1428-
/// 1432); an empty one keeps the chunk's. The sink computes them again for the e-values of the
/// subject's HSPs (`Runtime::subject`). Apart from the default path.
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn compute_own_subject(
    subjects: Subjects<'_>,
    batch: &PreparedQueryBatch,
    params: &[ContextParameters],
    option_spaces: &[i64],
    options: &ResolvedOptions,
    lookup: &BxLookup,
    diagonals: &mut Diagonals,
    xdrop: i32,
    oid: usize,
    x_containment: bool,
) -> Result<SubjectComputation> {
    let length = subjects.length(oid);
    let own = if length > 0 {
        Some(super::parameters::subject_update(
            batch,
            params,
            options,
            subjects.numbers()?,
            length,
            option_spaces,
        )?)
    } else {
        None
    };
    let params = own.as_deref().unwrap_or(params);
    let mut computed = compute_subject_scored(
        subjects,
        batch,
        params,
        options,
        lookup,
        diagonals,
        xdrop,
        oid,
        x_containment,
    )?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1448-1455 (cited at
    // `compute_subject_with`): a reader batch searched in units (`reader_batch.rs`) links with
    // `cutoff_score_min` of the word parameters, the minimum over the batch's contexts
    // (blast_parameters.c:401-416, `parameters::subject_word_cutoff_min`) at the subject's length.
    if let (Some(shape), Some(_)) = (batch.reader_shape.as_deref(), computed.prepared_link) {
        if length > 0 {
            computed.prepared_link = Some(super::linking::link_cutoffs(
                batch,
                params,
                length as i32,
                subjects.total_length()?,
                super::parameters::subject_word_cutoff_min(
                    &shape.batch_contexts,
                    options,
                    subjects.numbers()?,
                    length,
                )?,
            ));
        }
    }
    Ok(computed)
}

/// The body of [`compute_subject`] with the matrix of `scoring`: BLOSUM62's stages are today's
/// functions (`BxScoring::BLOSUM62`), another matrix's their instances for it. With `piece`, the
/// search of a piece of a masked database sequence (`search_db_masked.rs`: the whole sequence over
/// its soft ranges, or a hard chunk, whose link cutoffs of the chunk's length are not read); the
/// masked search calls these bodies (no code placement in NCBI), so that the stages keep one call
/// site each and stay inlined here.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn compute_subject_with<S: super::seed::BxScoring>(
    scoring: S,
    subjects: Subjects<'_>,
    batch: &PreparedQueryBatch,
    params: &[ContextParameters],
    options: &ResolvedOptions,
    lookup: &BxLookup,
    diagonals: &mut Diagonals,
    xdrop: i32,
    oid: usize,
    x_containment: bool,
    piece: Option<SubjectPiece<'_>>,
) -> Result<SubjectComputation> {
    #[cfg(feature = "blastx-worker-probe")]
    worker_probe("start", oid);
    let mut containment = Vec::new();
    // EXPERIMENT (LOSAT_X_BXUNGAPSENT): for a batch on the path of the seed loops' sentinel walks,
    // the subject with the sentinel before its residues too (`blastx/x_ungap_sent.rs`);
    // `encoded` is the same bytes either way. A masked sequence's piece brings its own.
    let x_lead = batch.x_ungap_sent;
    let x_buffer = match piece {
        Some(piece) => std::borrow::Cow::Borrowed(piece.buffer),
        None if x_lead => subjects.x_encoded_with_sentinels(oid),
        None => subjects.encoded_with_sentinel(oid),
    };
    let encoded: &[u8] = &x_buffer[usize::from(x_lead)..];
    // `subjects.length(oid)` for every kind of subjects: the residues before the sentinel byte (a
    // database read in blocks then looks the subject up once).
    let subject_length = encoded.len() - 1;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1451-1455
    // Worker-local cutoff preparation precedes WordFinder as in NCBI; only
    // adoption/observer notification is replayed on the serial OID collector.
    // `db_length` is `BlastSeqSrcGetTotLen` (blast_engine.c:1407), not the statistics' total:
    // `Subjects::total_length`.
    // With the statistics of each subject, `compute_own_subject` makes them again for a reader
    // batch searched in units.
    let prepared_link = if super::linking::computes_link_cutoffs(options) {
        Some(super::linking::link_cutoffs(
            batch,
            params,
            subject_length as i32,
            subjects.total_length()?,
            None,
        ))
    } else {
        None
    };

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:500-515
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1450-1458
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:502-511
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:494-500
    // The seed loop reads the scanner's per-thread buffer in place (x_seed_scan.rs;
    // LOSAT_X_BXDIAG, chosen inside `x_word_finder_scan`, reads the lookup cells in place instead,
    // x_seed_diag.rs). A `-subject` with lower-case letters (`-lcase_masking`) is scanned over its
    // scan ranges by the reference word finder (x_soft_subject.rs). The tests compare both with
    // the reference word finder (`seed::word_finder_bx`).
    // A masked sequence's piece: its own scan ranges (a database sequence has no lower case).
    let scan_ranges = super::x_soft_subject::x_subject_scan_ranges(subjects, oid);
    let ranges = match piece {
        Some(piece) => piece.ranges,
        None => scan_ranges.as_deref(),
    };
    let (raw, purged) = subject_hits_with(
        scoring,
        batch,
        params,
        options,
        lookup,
        diagonals,
        xdrop,
        &x_buffer,
        x_lead,
        ranges,
        x_containment,
        &mut containment,
    )?;
    #[cfg(feature = "blastx-worker-probe")]
    worker_probe("finish", oid);
    Ok(SubjectComputation {
        raw,
        purged,
        containment,
        prepared_link,
    })
}

/// The preliminary HSPs of a subject or of a piece of a masked database sequence
/// (`compute_subject_with`, its one caller): its word finder over the whole subject or over `ranges`, its
/// gapped or ungapped extensions, and the list before (`raw`) and after the purge of the HSPs
/// with common endpoints and the sort by score (`purged`). `x_buffer` is the subject's residues
/// and the sentinel after them, and with `x_lead` the sentinel before them too.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn subject_hits_with<S: super::seed::BxScoring>(
    scoring: S,
    batch: &PreparedQueryBatch,
    params: &[ContextParameters],
    options: &ResolvedOptions,
    lookup: &BxLookup,
    diagonals: &mut Diagonals,
    xdrop: i32,
    x_buffer: &[u8],
    x_lead: bool,
    ranges: Option<&[(i32, i32)]>,
    x_containment: bool,
    containment: &mut Vec<(PreliminaryHsp, bool)>,
) -> Result<(Vec<PreliminaryHsp>, Vec<PreliminaryHsp>)> {
    let encoded: &[u8] = &x_buffer[usize::from(x_lead)..];
    let initial = if let Some(ranges) = ranges {
        super::x_soft_subject::x_word_finder(
            scoring, batch, params, options, encoded, lookup, diagonals, ranges, None,
        )?
    } else {
        let x_sent = if x_lead {
            super::seed::x_ungap_sent::XSent::new(scoring, &batch.sequence_start, x_buffer)
        } else {
            None
        };
        super::seed::x_seed_scan::x_word_finder_scan(
            scoring, batch, params, options, encoded, lookup, diagonals, x_sent,
        )?
    };
    let subject_length = (encoded.len() - 1) as i32;
    match hits_of_initial(
        scoring,
        batch,
        params,
        options,
        encoded,
        subject_length,
        &initial,
        xdrop,
        x_containment,
        containment,
    )? {
        Some(hits) => Ok(hits),
        // The HSPs of a subject's letters lie inside the tree's bounds and are in order.
        None => anyhow::bail!("BLASTX: an interval tree insertion of a whole subject never ends"),
    }
}

/// The preliminary HSPs of the initial hits `initial` of a subject (`subject_hits_with`) or of a
/// hard chunk shorter than the word (`search_db_masked.rs`): `encoded` is its letters and the
/// sentinel after them, `subject_length` its length (NCBI's `subject->length`, the chunk's
/// `Int4`), and the gapped or ungapped HSPs, before (`raw`) and after the purge of the HSPs with
/// common endpoints and the sort by score (`purged`); `None` where NCBI's insertion of an HSP into
/// the interval tree never ends (`preliminary::gapped_observed`).
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn hits_of_initial<S: super::seed::BxScoring>(
    scoring: S,
    batch: &PreparedQueryBatch,
    params: &[ContextParameters],
    options: &ResolvedOptions,
    encoded: &[u8],
    subject_length: i32,
    initial: &[super::seed::InitHsp],
    xdrop: i32,
    x_containment: bool,
    containment: &mut Vec<(PreliminaryHsp, bool)>,
) -> Result<Option<(Vec<PreliminaryHsp>, Vec<PreliminaryHsp>)>> {
    let raw = if options.gapped {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3916-3921
        let mut observe = |h: &super::preliminary::PreliminaryHsp, contained: bool| {
            // No NCBI counterpart: the copy is skipped when nothing reads it; it does not
            // change any value NCBI computes.
            if x_containment {
                containment.push((h.clone(), contained))
            }
        };
        let gapped = if S::BLOSUM62 {
            super::preliminary::gapped_observed(
                batch,
                params,
                options,
                encoded,
                subject_length,
                initial,
                xdrop,
                &mut observe,
            )?
        } else {
            super::preliminary::gapped_observed_of(
                scoring,
                batch,
                params,
                options,
                encoded,
                subject_length,
                initial,
                xdrop,
                &mut observe,
            )?
        };
        let Some(raw) = gapped else {
            return Ok(None);
        };
        raw
    } else {
        ungapped(batch, initial)
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:552-563
    let mut purged = raw.clone();
    if options.gapped {
        purge_endpoints(&mut purged);
    } else {
        sort_by_score(&mut purged);
    }
    Ok(Some((raw, purged)))
}

/// [`hits_of_initial`] of a hard chunk shorter than the word by the search's matrix
/// (`search_db_masked.rs`): the seeds `initial` of its one word, on its stored letters `letters`
/// (the word's and the sentinel after them) with the chunk's length `length`.
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn one_word_hits(
    batch: &PreparedQueryBatch,
    params: &[ContextParameters],
    options: &ResolvedOptions,
    letters: &[u8],
    length: i32,
    initial: &[super::seed::InitHsp],
    xdrop: i32,
    x_containment: bool,
    containment: &mut Vec<(PreliminaryHsp, bool)>,
) -> Result<Option<(Vec<PreliminaryHsp>, Vec<PreliminaryHsp>)>> {
    if options.scoring == crate::config::ScoringMatrix::Blosum62 {
        hits_of_initial(
            super::seed::Blosum62Bx,
            batch,
            params,
            options,
            letters,
            length,
            initial,
            xdrop,
            x_containment,
            containment,
        )
    } else {
        hits_of_initial(
            super::seed::MatrixBx::of(options.scoring),
            batch,
            params,
            options,
            letters,
            length,
            initial,
            xdrop,
            x_containment,
            containment,
        )
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1659-1664
// Diagnostic-only host clock/worker identity; never read by search or reporting.
#[cfg(feature = "blastx-worker-probe")]
fn worker_probe(event: &str, oid: usize) {
    let ns = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .expect("host clock after epoch")
        .as_nanos();
    crate::utils::threading::probe_write(format!(
        "[blastx-worker-probe] event={event} oid={oid} thread={:?} ns={ns}\n",
        std::thread::current().id()
    ));
}
