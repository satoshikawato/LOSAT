//! psiblast's rounds of FASTA queries (`-num_iterations` other than 1; series PS, item PS5a) and of
//! a PSSM query (`-in_pssm` with `-num_iterations` other than 1 or a PSSM file, item PS5b): the
//! round driver (`search_rounds`), which `search` calls for such a run, and the hook of one
//! round's search (`round_outcome`), which `run_resolved_in_pool` calls beside `psiblast_pssms`.
//! Both are cold: blastp's default path only tests `ResolvedBlastpArgs::psi_rounds` and
//! `BlastpReports::psi_round`.
//!
//! The driver runs `run_resolved_in_pool` once per round per query, as NCBI's `CPsiBlast` searches
//! one query per batch and round (psiblast_app.cpp:468, `x_RunLocalPsiBlastIterations`; cited at
//! `algorithm::psiblast::rounds`), on one search pool and the database opened once: round 1 is
//! blastp's search of the FASTA query (with `-in_pssm` the PSSM search of the input PSSM), a later
//! round the PSSM search (`algorithm::pssm_search`) of the PSSM that the engine made from the round
//! before (`PsiSearch::from_engine`), whose query stays the round-1 query. The outfmt 0
//! prolog comes with the first call; every call writes no epilog, which the driver writes after
//! the last query (outfmt 7's count is that of the rounds). The PSSM files are written between
//! the calls, where NCBI writes them, after the report of the round before (with `-in_pssm` the
//! input PSSM before round 1, `algorithm::psiblast::in_pssm_files`).

use std::io::Write as _;
use std::sync::Arc;

use anyhow::{bail, Result};

use super::{
    blastp_database_name, run_resolved_in_pool, write_pairwise_prologs,
    write_reports_without_queries, BlastpRanges, BlastpReports, BlastpSubjectSource, ShownLetters,
    ShownSubjects, SubjectSet,
};
use crate::algorithm::blastp::args::{BlastpCompositionMode, ResolvedBlastpArgs};
use crate::algorithm::blastp::hsp::BlastpHitList;
use crate::algorithm::psiblast::pssm::{
    ncbistdaa_letters, query_pssm, PssmAlignment, PssmAncillary, PssmQueryRecord, PssmScoring,
    QueryPssm,
};
use crate::algorithm::psiblast::rounds::{
    HandleIndex, InclusionSet, IterationState, PsiRoundIo, PsiRoundPlan, RoundFiles, RoundOutcome,
    SetKey,
};
use crate::algorithm::pssm_search::PsiSearch;
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::api::local_blast::ReportOutputs;
use crate::blastinput::app::CompositionMode;
use crate::blastinput::app::ErrorAfterBatchReports;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::query_batch::QueryBatches;
use crate::blastinput::seq_range;
use crate::core::psi::KarlinBlk;
use crate::report::db_subject::ParsedSubject;
use crate::stats::KarlinParams;

/// The run of a psiblast search of FASTA queries with `-num_iterations` other than 1, after the
/// query batches are read (`search`): the rounds of each searched batch's query, the epilog, then
/// the stop of the batch after them.
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(super) fn search_rounds(
    args: &ResolvedBlastpArgs,
    plan: &PsiRoundPlan,
    query_records: &[FastaRecord],
    searched_queries: &[FastaRecord],
    query_input: &seq_range::QueryInput,
    batches: QueryBatches,
    subjects: &SubjectSet<'_>,
    subject_source: BlastpSubjectSource<'_>,
    subject_label: &str,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    let epilog = batches.stop.is_none();
    // The search strategy's requests (`-export_search_strategy`): one per query, written before
    // its first round (`run_queries`), and the one of the batch that stops the run
    // (`QueryBatches::stop_batch`). NCBI's posts of the bytes replaced in a query's request follow
    // the reader's messages of its batch: round 1's call places them before the query's report
    // (`CallRequests::take_posts`), and an error of the batch after those messages
    // (`write_request_posts`).
    let requests = args.strategy.as_ref().map(|_| {
        crate::algorithm::blastp::archive::requests(
            subject_source
                .database()
                .is_none()
                .then(|| subject_source.fasta_records()),
            args.lcase_masking,
            args.psi_search.as_ref().and_then(|psi| psi.request_pssm()),
        )
    });
    let write_stop_request = || -> Result<()> {
        match (
            args.strategy.as_ref(),
            requests.as_ref(),
            batches.stop_batch.clone(),
        ) {
            (Some(strategy), Some(requests), Some(span)) => requests.write(
                strategy,
                span.start,
                &searched_queries[span.clone()],
                query_input.placements.get(span.start),
            ),
            _ => Ok(()),
        }
    };
    if batches.searched_done == 0 {
        // As `search`: the first batch's stop after the prolog, or the reports of no query.
        if let Some((_, error)) = batches.stop.as_ref() {
            // The prologs of the first report, then the request of a batch that NCBI made a query
            // factory of; NCBI's stream that throws while it writes the request
            // (`encoder::NullElement`) ends the run after the prolog (psiblast_app.cpp:506, 340).
            write_pairwise_prologs(args, subject_source, subject_label, outputs, false)?;
            write_stop_request()?;
            // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:246-248 (cited at
            // `in_pssm_files`: the loop writes the PSSM of each round at its start, the input PSSM
            // of `-in_pssm` in round 1)
            // An error of the PSSM query's search itself (NCBI's error, not LOSAT's refusal: a
            // `kbp_psi` that is not valid, a database without sequences or whose sequence file
            // cannot be opened, `psiblast_stop_at_invalid_query`) comes after round 1 has written
            // the input PSSM's files; NCBI's check of the PSSM (`validation_error`) stops before
            // (NCBI 2.17.0, probes `y.*`, `y2.*` of the PS5c review; line PS DECISIONS 31), and so
            // do the search strategy's request without a database and the PSSM engine's scores
            // (`PsiSearch::stop_before_search`: `SaveSearchStrategy` and `CPsiBlast` before the
            // rounds), but not the query set-up that `CPsiBlast::Run` starts with (the query's
            // invalid residues, `PssmQueryStop::after_round_files`).
            // The object manager's exception (`PsiSearch::annot_abort`) stops before round 1 too.
            let search_error = error.downcast_ref::<crate::cli::NativeError>().is_some()
                && args.psi_search.as_deref().is_some_and(|psi| {
                    psi.annot.abort.is_none()
                        && psi
                            .stop_before_search(
                                args.strategy.is_some(),
                                subject_source.database().is_some(),
                            )
                            .is_none_or(|stop| stop.after_round_files())
                });
            if let Some(input) = plan.input.as_ref().filter(|_| search_error) {
                let files = input.files(&plan.request)?;
                write_file_posts(args, &files.posts, outputs)?;
                RoundFiles::new(
                    &plan.request,
                    args.out.as_deref().filter(|path| path.as_os_str() != "-"),
                )
                .write(files.asn.as_deref(), files.ascii.as_deref())?;
                // NCBI's throw in the ASN.1 file ends the run there, before the search.
                if let Some(thrown) = files.thrown {
                    return Err(thrown.into());
                }
            }
        } else {
            write_reports_without_queries(args, subject_source, subject_label, outputs)?;
        }
    } else {
        crate::utils::threading::with_search_pool(args.num_threads, "blastp", |pool| {
            run_queries(
                args,
                plan,
                query_records,
                searched_queries,
                query_input,
                &batches,
                subjects,
                subject_source,
                subject_label,
                outputs,
                pool,
                epilog,
                requests.as_ref(),
            )
        })?;
        // The request of the batch that stops the run follows the queries' rounds (NCBI wrote it
        // before the search that fails).
        write_stop_request()?;
    }
    // The stopping batch's reader messages, NCBI's posts of the bytes replaced in its request
    // (written above, after the rounds of the queries before), its error.
    if let Some((warnings, error)) = batches.stop {
        outputs.diagnostics.write_all(&warnings)?;
        super::write_request_posts(args, &mut *outputs.diagnostics)?;
        return Err(error);
    }
    Ok(())
}

/// `SavePssmToFile` of a round's PSSM `current` (cited at `algorithm::psiblast::pssm::PssmFiles`):
/// NCBI's posts of the bytes that its stream replaces in the ASN.1 text come as it writes the file
/// (`pssm::file_posts`), after what the run wrote before; then the files (`RoundFiles::write`).
fn write_round_pssm(
    args: &ResolvedBlastpArgs,
    files: &mut RoundFiles,
    current: &QueryPssm,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    write_file_posts(args, &current.asn_posts, outputs)?;
    files.write(current.asn.as_deref(), current.ascii.as_deref())
}

/// NCBI's posts of the bytes that its ASN.1 text stream replaces in a PSSM file
/// (`pssm::file_posts`), made as the file is written.
fn write_file_posts(
    args: &ResolvedBlastpArgs,
    replaced: &[crate::report::serial::encoder::ReplacedChar],
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    if !replaced.is_empty() {
        crate::algorithm::psiblast::pssm::file_posts(args.app.ncbi_name(), replaced)
            .write(&mut *outputs.diagnostics, args.archive.as_ref())?;
        outputs.diagnostics.flush()?;
    }
    Ok(())
}

/// The composition-based mode of the options' reader (`in_pssm::setup_messages`).
fn app_mode(mode: BlastpCompositionMode) -> CompositionMode {
    match mode {
        BlastpCompositionMode::NoCompositionBasedStats => CompositionMode::NoCompositionBasedStats,
        BlastpCompositionMode::CompositionBasedStats => CompositionMode::CompositionBasedStats,
        BlastpCompositionMode::CompositionMatrixAdjust => CompositionMode::CompositionMatrixAdjust,
        BlastpCompositionMode::ForceFullMatrixAdjust => CompositionMode::CompoForceFullMatrixAdjust,
    }
}

/// The rounds of every searched batch's query on `pool`, then the epilog (`epilog`).
// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:508-532
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:196-205 (cited at
// `in_pssm::setup_messages`: a PSSM search of mode 2 or 3 sets the options to 1, so a later
// round of the query neither runs 2 or 3 nor warns again)
#[allow(clippy::too_many_arguments)]
fn run_queries(
    args: &ResolvedBlastpArgs,
    plan: &PsiRoundPlan,
    query_records: &[FastaRecord],
    searched_queries: &[FastaRecord],
    query_input: &seq_range::QueryInput,
    batches: &QueryBatches,
    subjects: &SubjectSet<'_>,
    subject_source: BlastpSubjectSource<'_>,
    subject_label: &str,
    outputs: &mut ReportOutputs<'_>,
    pool: &crate::utils::threading::SearchPool<'_>,
    epilog: bool,
    requests: Option<&crate::report::archive::requests::RunRequests<'_>>,
) -> Result<()> {
    let mut files = RoundFiles::new(
        &plan.request,
        args.out.as_deref().filter(|path| path.as_os_str() != "-"),
    );
    // With `-in_pssm` the options' mode is already the PSSM search's (`check_options`), so a later
    // round warns no more.
    let original_mode = args.comp_based_stats.mode;
    let pssm_mode =
        crate::algorithm::psiblast::in_pssm::position_based_mode(app_mode(original_mode));
    let pssm_round_mode = match pssm_mode {
        CompositionMode::NoCompositionBasedStats => BlastpCompositionMode::NoCompositionBasedStats,
        _ => BlastpCompositionMode::CompositionBasedStats,
    };
    // The rounds formatted so far (`m_QueriesFormatted`: each round's `PrintOneResultSet`).
    let mut calls = 0usize;
    for (batch, &(input_start, query)) in batches.starts.iter().enumerate() {
        let records = std::slice::from_ref(&searched_queries[query]);
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:144-155
        // (cited at `seq_range::cut_queries`: a record whose range starts past its end is
        // dropped from the batch after its reader added it to the scope)
        // NCBI's scope at each round's prefetch (`BlastpReports::prefetch_scope`), made only where
        // it stops: the input records that NCBI's reader read for the query's batch (those that
        // `-query_loc` drops before the query included; each with its `Query_N` of the run,
        // `FastaRecord::local_id`), which every round of the query keeps; with `-in_pssm` the
        // PSSM's query (`PsiSearch::query_scope`; a later round's PSSM has no ids of its own).
        let input_end = batches
            .starts
            .get(batch + 1)
            .map_or(batches.input_done, |&(input, _)| input);
        let batch_records = &query_records[input_start..input_end];
        let pssm_query = args.psi_search.as_deref();
        let prefetch_scope = || match pssm_query {
            Some(psi) => psi.query_scope(),
            None => crate::blastinput::query_batch::QueryScope::of_shown_records(batch_records),
        };
        // The query as round 1 searches it, in every round: with `-query_loc` its range of the
        // record (`placement`; the PSSM is made of the whole record, `next_pssm`, and searched
        // at the range's offsets, `PsiSearch::from_engine`). A psiblast batch holds one searched
        // query, after the records that `-query_loc` skips (a range that starts past the end of
        // its record; NCBI's batch reader drops them and numbers the query by its place in the
        // input, `Query_<n>`): the query's record is its input ordinal's, not the batch's first.
        let ordinal = query_input.ordinal(query);
        let record = &query_records[ordinal];
        let placement = query_input.placements.get(query);
        let query_input = seq_range::QueryInput {
            placements: seq_range::Placements(placement.into_iter().collect()),
            ordinals: vec![0],
            input_lengths: vec![record.seq().len()],
            skipped: vec![false],
            range_length: query_input.range_length,
        };
        let ranges = BlastpRanges {
            input_records: std::slice::from_ref(record),
            query_input: &query_input,
            subject_placements: &subjects.placements,
            x_ordinal_base: ordinal,
            searched_base: query,
            split_search_spaces: None,
        };
        let mut state = IterationState::new(plan.num_iterations);
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:118-130
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_impl.cpp:115-119
        // Every round of a FASTA query searches the query factory of its batch, which NCBI splits
        // as blastp's (`split_protein_batch`); a PSSM query's (`-in_pssm`) is never split.
        // Sizes that NCBI cannot convert stopped the first batch before any (`search_rounds`).
        let split = args.psi_search.is_none()
            && crate::algorithm::blastp::protein_query_split::split_protein_batch(
                &[records[0].seq().len()],
                crate::algorithm::blastp::protein_query_split::protein_query_split_sizes(10000)?,
            )
            .is_some();
        // The query as the round before searched it, from round 2 on (`RoundQuery`).
        let mut round_query: Option<RoundQuery> = None;
        // The PSSM that the next round searches (`pssm`, made after the round before); before
        // round 1 of `-in_pssm` the input PSSM (`args.psi_search`), whose files are `input_files`.
        let mut pssm: Option<QueryPssm> = None;
        let mut input_files = match &plan.input {
            Some(input) => Some(input.files(&plan.request)?),
            None => None,
        };
        // `kbp_psi[0]` and `kbp_gap_psi[0]` of the input PSSM's search, after a round 1 without
        // alignments (`in_pssm_files::InPssmRounds::written_back`).
        let mut searched_input = None;
        // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:337-340
        // The query's search strategy request, before its first round (also when round 1 fails
        // with a header file NCBI cannot read). Round 1's call notes the query's searched space
        // (`CallRequests`, which finds the request written) and writes NCBI's posts of the bytes
        // replaced in it after the reader's messages of the query's batch, before its report
        // (`CallRequests::take_posts`); the later rounds' calls find all three done. NCBI's
        // handle keeps the last round's space, which is round 1's: every round's effective
        // lengths read the same `kbp_gap_std` (blast_setup.c:766-768).
        if let (Some(strategy), Some(requests)) = (args.strategy.as_ref(), requests) {
            requests.write(strategy, query, records, placement)?;
        }
        while state.go() {
            let round = state.number();
            match (&pssm, input_files.as_mut()) {
                (Some(current), _) => write_round_pssm(args, &mut files, current, outputs)?,
                (None, Some(input)) if round == 1 => {
                    // NCBI's posts of the bytes replaced in the search strategy request (written
                    // above, before round 1) come before those of the PSSM file that round 1
                    // writes.
                    super::write_request_posts(args, &mut *outputs.diagnostics)?;
                    write_file_posts(args, &input.posts, outputs)?;
                    // NCBI's prolog comes before the rounds (`PrintProlog`), so before a file
                    // that cannot be written, or NCBI's throw in it, stops the run here.
                    let written = files
                        .write(input.asn.as_deref(), input.ascii.as_deref())
                        .and_then(|()| match input.thrown.take() {
                            Some(thrown) => Err(thrown.into()),
                            None => Ok(()),
                        });
                    if let Err(error) = written {
                        if calls == 0 {
                            write_pairwise_prologs(
                                args,
                                subject_source,
                                subject_label,
                                outputs,
                                false,
                            )?;
                        }
                        return Err(error);
                    }
                    // Only the ASCII PSSM is written again (after a round 1 without alignments;
                    // the ASN.1 text is then made with the write-back).
                    input.asn = None;
                }
                _ => {}
            }
            let mut round_args = args.clone();
            if let Some(current) = &pssm {
                // The options' mode when this round's score block is set up: the given one in
                // the query's first PSSM round, 1 (or 0) after it.
                let setup_mode = if round == 2 {
                    app_mode(original_mode)
                } else {
                    pssm_mode
                };
                // `PsiBlastSetupScoreBlock`'s warnings look at every column of the PSSM, also
                // past a range that the search reads (`pssm_setup_messages`).
                let setup_messages = crate::algorithm::psiblast::in_pssm::pssm_setup_messages(
                    current
                        .engine
                        .diagnostics
                        .frequency_ratios
                        .as_deref()
                        .unwrap_or(&[]),
                    setup_mode,
                    &args.matrix_name,
                );
                let mut search =
                    PsiSearch::from_engine(&current.engine, records[0].clone(), setup_messages);
                // The PSSM that the round's archive holds (made only for a run that writes one,
                // `query_pssm`'s `archive`).
                search.archive_pssm = args.archive.as_ref().and(current.object.clone());
                // The query source is fixed for the run: a PSSM query's warning label keeps its
                // title rule in every round (`PsiSearch::warning_title`).
                search.warning_title = args
                    .psi_search
                    .as_deref()
                    .and_then(|psi| psi.warning_title.clone());
                // So are the letters of its Bioseq that its reports read (`PsiSearch::
                // bioseq_letters`): NCBI's reports of every round read the query's Bioseq.
                search.bioseq_letters = args
                    .psi_search
                    .as_deref()
                    .and_then(|psi| psi.bioseq_letters.clone());
                // And its molecule: a query of the molecule `other` is unreadable for the
                // reports of every round (`PsiSearch::query_other_coding`).
                search.query_other_coding = args
                    .psi_search
                    .as_deref()
                    .and_then(|psi| psi.query_other_coding);
                // A query that NCBI splits: round 1's query, which the round searches, and the query
                // of its set-up (`RoundQuery::split`).
                let (query_now, split_setup) = if split {
                    let (query, setup) = RoundQuery::split(args, &records[0]);
                    (query, Some(setup))
                } else {
                    let query = match round_query.take() {
                        Some(before) => before.next(args),
                        None => RoundQuery::first(args, &records[0]).next(args),
                    };
                    (query, None)
                };
                // A later round whose letters give no Karlin-Altschul parameters (SEG again on the
                // letters that the rounds before masked can mask every letter left): NCBI's
                // setup crashes there, after this round's PSSM files, as for a PSSM query
                // (`INVALID_PSSM_QUERY_REFUSAL`).
                //
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2782-2799
                // The set-up's query (`kbp_std`), which a split query's round reads too.
                let frames = [vec![split_setup
                    .as_ref()
                    .unwrap_or(&query_now.frame)
                    .clone()]];
                if !super::pssm_hooks::psi_contexts(args, &frames)[0].is_valid {
                    bail!(super::INVALID_PSSM_QUERY_REFUSAL);
                }
                search.round_query = Some(query_now.frame.clone());
                // A split query's next round starts again from the record (`RoundQuery::split`).
                if split_setup.is_none() {
                    round_query = Some(query_now);
                }
                search.split_setup = split_setup;
                round_args.psi_search = Some(Arc::new(search));
                round_args.comp_based_stats.mode = pssm_round_mode;
            }
            let io = PsiRoundIo {
                state: state.clone(),
                save_last: plan.save_last,
                request: plan.request.clone(),
                query_id: args
                    .psi_search
                    .as_deref()
                    .and_then(|psi| psi.query_seq_loc_id().cloned()),
                outcome: std::sync::Mutex::new(None),
            };
            // The library's records of the round (`FormatObserver::round_begin`): the call's hit
            // list and HSP numbers are the round's.
            if let Some(observer) = outputs.observer.as_deref_mut() {
                observer.round_begin(query, round);
            }
            // The reader's messages of the batch come with its first round.
            let reader_warnings: &[Vec<u8>] = if round == 1 {
                &batches.warnings[query..query + 1]
            } else {
                &[]
            };
            let result = run_resolved_in_pool(
                round_args,
                records,
                subject_source,
                "",
                subject_label,
                outputs,
                pool,
                None,
                &ranges,
                &BlastpReports {
                    reader_warnings,
                    identity_warned: 0,
                    epilog: false,
                    x_prolog: calls == 0,
                    // The result sets formatted before (`m_QueriesFormatted`): BLAST XML numbers
                    // each round's `Iteration` after them and writes its header with the first.
                    x_queries_before: calls,
                    x_xml2_carry: None,
                    last_call: false,
                    psi_round: Some(&io),
                    stopped_pssms: None,
                    query_scopes: batches
                        .query_scopes
                        .as_ref()
                        .map(|scopes| scopes.numbered_from(query)),
                    prefetch_scope: Some(&prefetch_scope),
                    strategy_next: None,
                    call_requests: None,
                    x_plan_blocks: false,
                },
                super::x_lookup_ahead::XAheadBlock::none(),
            );
            // A header file that NCBI cannot read fails the round whose results reach it first,
            // at its prefetch, after the reports of the rounds before (`search`'s rule for the
            // batches): the outfmt 0 prolog if nothing was written, and the reader's messages
            // of a batch whose first round fails.
            if let Err(error) = result {
                if error.is::<ErrorAfterBatchReports>() {
                    let failure = error
                        .downcast::<ErrorAfterBatchReports>()
                        .expect("the error was checked to be an ErrorAfterBatchReports");
                    if calls == 0 {
                        write_pairwise_prologs(
                            args,
                            subject_source,
                            subject_label,
                            outputs,
                            false,
                        )?;
                    }
                    if round == 1 {
                        outputs.diagnostics.write_all(&batches.warnings[query])?;
                        // NCBI's posts of the bytes replaced in the query's request, which
                        // round 1's report did not place.
                        super::write_request_posts(args, &mut *outputs.diagnostics)?;
                    }
                    // NCBI's warning of its first read of a UTF8String header, where the round's
                    // Seq-aligns made it before the failure (once per run: a round before posted
                    // its own).
                    if let Some(posts) = subject_source.database().and_then(|database| {
                        failure
                            .utf8
                            .as_ref()?
                            .post(database.seqdb(), args.app.ncbi_name())
                    }) {
                        posts.write(&mut *outputs.diagnostics, args.archive.as_ref())?;
                    }
                    return Err(failure.error);
                }
                return Err(error);
            }
            calls += 1;
            let outcome = io
                .outcome
                .lock()
                .expect("the round's outcome")
                .take()
                .expect("every round's search reaches its outcome hook");
            if !outcome.has_alignments {
                if round == 1 {
                    searched_input = outcome.searched;
                }
                break;
            }
            state.advance(outcome.set);
            if state.go() || plan.save_last {
                // The PSSM's warnings and the engine's error come here, after the round's report
                // (no epilog).
                outputs
                    .diagnostics
                    .write_all(&outcome.next_warnings.bytes())?;
                pssm = Some(
                    outcome
                        .next
                        .expect("the hook makes the PSSM where NCBI makes it")?,
                );
            }
        }
        // After a round without alignments NCBI writes the PSSM it searched last, after that
        // round's K/A write-back (psiblast_impl.cpp:150-168): an engine PSSM's own values (the
        // same bytes), or the input PSSM's after round 1 (`InPssmRounds::written_back`), whose
        // ASCII PSSM is round 1's (`m_AncillaryData` is still the input PSSM's).
        if plan.save_last {
            match (&pssm, &plan.input, &input_files, searched_input) {
                (Some(current), _, _, _) => write_round_pssm(args, &mut files, current, outputs)?,
                (None, Some(input), Some(round_1), Some((ungapped, gapped))) => {
                    let (asn, posts) = input.written_back(&plan.request, ungapped, gapped)?;
                    write_file_posts(args, &posts, outputs)?;
                    files.write(asn.as_deref(), round_1.ascii.as_deref())?;
                }
                _ => {}
            }
        }
        if state.has_converged() {
            if let Some(observer) = outputs.observer.as_deref_mut() {
                observer.rounds_converged(query);
            }
            write_converged(outputs)?;
        }
        // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:513
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2295-2314
        if let Some(database) = subject_source.database() {
            database.reset_scope();
        }
    }
    if epilog {
        write_round_epilogs(args, subject_source, subject_label, outputs, calls)?;
    }
    Ok(())
}

/// The query of one round of a psiblast query, as NCBI's search sets it up. NCBI keeps the
/// query's factory for all its rounds (`CPsiBlast` holds it; `SetPssm` changes the PSSM only), and
/// the factory keeps its local query data, so every round's setup works on the same sequence
/// block: the hard masks of a round (`BlastSetUp_MaskQuery`, not with `-soft_masking true`)
/// write X over the letters that the next round filters again, the lower-case mask is taken
/// out of the block by round 1's filtering (`seqloc_array[context]` set to NULL), and the unmasked
/// copy that the identities read is made again only by a round that masks letters. A query that
/// NCBI splits into chunks is the exception: its block is made again from the record in every
/// round (`RoundQuery::split`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/query_data.cpp:52-58
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:206-210
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1241-1249
/// (`s_GetFilteringLocationsForOneContext`)
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1375-1382
/// (`BlastSetUp_MaskQuery`)
struct RoundQuery {
    /// The letters that the round filtered: the record's, with X where the rounds before
    /// masked them.
    letters: Vec<u8>,
    /// The round's query frame: its masks, its letters masked, and the letters that its
    /// identities read (`aa_seq_nomask`).
    frame: crate::algorithm::tblastx::translation::QueryFrame,
}

impl RoundQuery {
    /// Round 1's query: the record's letters, as the search encodes them
    /// (`blastp_query_frames`: `-seg`, `-lcase_masking`, `-soft_masking`).
    fn first(args: &ResolvedBlastpArgs, record: &FastaRecord) -> Self {
        let frame = super::blastp_query_frames(args, std::slice::from_ref(record))
            .pop()
            .and_then(|mut frames| frames.pop())
            .expect("one frame of one protein query");
        Self {
            letters: record.seq().to_vec(),
            frame,
        }
    }

    /// A later round's query of a query that NCBI splits into chunks, and the query of the
    /// round's set-up. NCBI's set-up of a split query frees its sequence block after the score
    /// block (`kbp_std`) is made, and the preliminary stage makes it again from the record after
    /// the chunks, for the traceback: so every round searches round 1's query (the record's
    /// letters, filtered with the lower-case mask), and the set-up of a later round reads the
    /// block that the round before made again, filtered again (`next`) — not the rounds before
    /// it, as a query that NCBI does not split carries them (Q4b).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:249-256
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:292-301
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:206-219 (a freed
    /// block is set up again from the query's Seq-loc: `SafeSetupQueries`, with the lower-case
    /// mask of the query vector)
    fn split(
        args: &ResolvedBlastpArgs,
        record: &FastaRecord,
    ) -> (Self, crate::algorithm::tblastx::translation::QueryFrame) {
        let first = Self::first(args, record);
        let setup = first.next(args).frame;
        (first, setup)
    }

    /// The next round's query: this round's letters with X over its masks (none with
    /// `-soft_masking true`), filtered by SEG again, without the lower-case mask. Its identities
    /// read its letters when it masks some, else what this round's read.
    fn next(&self, args: &ResolvedBlastpArgs) -> Self {
        let mut letters = self.letters.clone();
        if !args.soft_masking {
            for &(start, end) in &self.frame.seg_masks {
                letters[start..end].fill(b'X');
            }
        }
        let seg = args.seg.params();
        let mut frame = if args.soft_masking {
            crate::algorithm::blastp::encoding::encode_protein_query_frame_soft_masked(
                &letters,
                seg.as_ref(),
                false,
            )
        } else {
            crate::algorithm::blastp::encoding::encode_protein_query_frame_with_seg(
                &letters,
                seg.as_ref(),
            )
        };
        if frame.aa_seq_nomask.is_none() {
            frame.aa_seq_nomask = self.frame.aa_seq_nomask.clone();
        }
        Self { letters, frame }
    }
}

/// The query frames of a later round's search (`PsiSearch::round_query`, set by `run_queries`).
#[cold]
#[inline(never)]
pub(super) fn round_query_frames(
    args: &ResolvedBlastpArgs,
) -> Vec<Vec<crate::algorithm::tblastx::translation::QueryFrame>> {
    let frame = args
        .psi_search
        .as_deref()
        .and_then(|psi| psi.round_query.clone())
        .expect("a later round's query");
    vec![vec![frame]]
}

/// Whether the search read other letters than its query record's: a later round of psiblast
/// (`RoundQuery`), whose HSPs the reports count again on the record's letters
/// (`recount_shown_identities`).
pub(super) fn searched_other_letters(args: &ResolvedBlastpArgs) -> bool {
    args.psi_search
        .as_deref()
        .is_some_and(|psi| psi.round_query.is_some())
}

/// The query contexts that the reports of a later round read: the search's, with the record's
/// letters as the unmasked copy. NCBI's formatters read the query from the scope (the record's
/// Bioseq): the rows and identities of the pairwise report, and the tabular fields, which count
/// the identities again from the rows. The round's own unmasked copy, made of the letters that the
/// rounds before masked (`RoundQuery`), is read only by the engine's identity count
/// (`num_ident`), which a tabular report prints where it does not count them again
/// (`report::tabular`: `nident` alone with an engine count of 0). Before round 2 (and in every
/// other search) the contexts are the search's.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:985-1007
pub(super) fn report_contexts<'a>(
    args: &ResolvedBlastpArgs,
    contexts: &'a [QueryContext],
    query_records: &[FastaRecord],
) -> std::borrow::Cow<'a, [QueryContext]> {
    if !searched_other_letters(args) {
        return std::borrow::Cow::Borrowed(contexts);
    }
    later_round_report_contexts(contexts, query_records)
}

#[cold]
#[inline(never)]
fn later_round_report_contexts<'a>(
    contexts: &[QueryContext],
    query_records: &[FastaRecord],
) -> std::borrow::Cow<'a, [QueryContext]> {
    let mut contexts = contexts.to_vec();
    for context in &mut contexts {
        let record = &query_records[context.q_idx as usize];
        context.aa_seq_nomask =
            Some(crate::algorithm::blastp::encoding::encode_protein_sequence(record.seq()).aa_seq);
    }
    std::borrow::Cow::Owned(contexts)
}

// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:525-528
// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:556-558
/// The line `Search has CONVERGED!` between two `NcbiEndl`s, for every report that is not
/// structured (psiblast's outfmt 0, 6, 7, 10 and 20), on the report's one stream: the
/// application writes it after the query's rounds, outside the formatter's regions (`Site::Report`:
/// no `CBlastTabularInfo` lives there, so a failed write there is the application's "BLAST failed
/// to write output", exit 6), each `NcbiEndl` a line feed and a flush (`ostream::ncbi_endl`). The
/// next query's report starts in the format's region again (`Site::of_format`).
fn write_converged(outputs: &mut ReportOutputs<'_>) -> Result<()> {
    use crate::blastinput::app::ReportFormat;
    use crate::report::ostream::{ncbi_endl, ReportWrite, Site};
    for format in outputs.formats.iter_mut() {
        let number = crate::blastinput::app::parse_formatting_string(format.outfmt)?.number;
        if matches!(
            ReportFormat::from_number(number),
            Some(
                ReportFormat::Pairwise
                    | ReportFormat::Tabular
                    | ReportFormat::TabularWithComments
                    | ReportFormat::CommaSeparatedValues
                    | ReportFormat::CommaSeparatedValuesWithHeader
            )
        ) {
            let mut writer = crate::cli::NcbiOstream::sink(&mut format.sink)?;
            writer.site(Site::Report)?;
            ncbi_endl(&mut writer)?;
            writer.write_all(b"Search has CONVERGED!")?;
            ncbi_endl(&mut writer)?;
            writer.site(Site::of_format(number))?;
            writer.flush()?;
        }
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2233-2238 (cited at
// `tabular::write_tabular_without_queries`: outfmt 7's count of the result sets formatted)
/// The epilog of each report after the last round: outfmt 0's, outfmt 7's `# BLAST processed R
/// queries` (R the rounds formatted, `m_QueriesFormatted`), BLAST XML's end, nothing for outfmt
/// 6, 8 to 11 and 20.
fn write_round_epilogs(
    args: &ResolvedBlastpArgs,
    subject_source: BlastpSubjectSource<'_>,
    subject_label: &str,
    outputs: &mut ReportOutputs<'_>,
    rounds: usize,
) -> Result<()> {
    let (num_sequences, total_letters) = subject_source.shown_totals();
    let pairwise_database = subject_source
        .database()
        .map(super::database::BlastpDatabase::pairwise_database);
    let taxdb_found = match subject_source.database() {
        Some(database) => database.taxdb_found(),
        None => args.subject_taxdb_found,
    };
    for format in outputs.formats.iter_mut() {
        if let Some(spec) = super::super::tabular::custom_tabular_spec(format.outfmt, args.app) {
            let mut writer = crate::cli::NcbiOstream::sink(&mut format.sink)?;
            super::super::tabular::write_tabular_round_epilog(
                spec,
                taxdb_found,
                args.app.ncbi_name(),
                &mut writer,
                &mut *outputs.diagnostics,
                rounds,
            )?;
            writer.flush()?;
            continue;
        }
        match crate::blastinput::app::parse_formatting_string(format.outfmt)?.number {
            0 => {}
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2242-2246
            // BLAST XML's end, after the rounds' `Iteration`s (`BlastXmlWriter::finish`).
            5 => {
                let mut bytes = Vec::new();
                crate::report::blast_xml::StructuredXml::finish(
                    &mut crate::report::blast_xml::BlastXmlWriter::new(rounds),
                    &mut bytes,
                )?;
                let mut writer = format.sink.open()?;
                writer.write_all(&bytes)?;
                writer.flush()?;
                continue;
            }
            7 => {
                // The default outfmt 7 (`write_blastp_tabular_hit_lists`' epilog).
                let mut writer = format.sink.open()?;
                writeln!(writer, "# BLAST processed {rounds} queries")?;
                writer.flush()?;
                continue;
            }
            _ => continue,
        }
        let report = crate::report::BlastpPairwiseReport {
            version: crate::wording::REPORT_VERSION.to_string(),
            database_name: blastp_database_name(args, subject_label).into(),
            database_num_sequences: num_sequences,
            database_total_letters: total_letters,
            matrix_name: Some(args.matrix_name.clone()),
            gap_open: args.scoring.gap_open,
            gap_extend: args.scoring.gap_extend,
            word_threshold: args.threshold,
            window_size: args.window_size as i32,
            gapped_karlin: crate::stats::lookup_protein_params(&args.scoring),
            gumbel: None,
            gapped: true,
            show_linked_set_size: false,
            num_descriptions: args.num_descriptions,
            num_alignments: args.num_alignments,
        };
        let mut writer = crate::cli::NcbiOstream::sink(&mut format.sink)?;
        crate::report::pairwise::psiblast_round::write_psiblast_epilog(
            &mut writer,
            &report,
            pairwise_database.as_ref(),
        )?;
    }
    Ok(())
}

/// A Karlin-Altschul block of the PSSM engine (`core::psi`).
fn engine_block(params: &KarlinParams) -> KarlinBlk {
    KarlinBlk {
        lambda: params.lambda,
        k: params.k,
        h: params.h,
    }
}

/// The outcome of one round's search (`PsiRoundIo::outcome`), from the round's hit list of its
/// one query before the formatter's cut (`psiblast_pssms`' Seq-aligns): whether it has
/// alignments, its inclusion set, and where NCBI makes one the next round's PSSM (an error of
/// the engine is kept for after the round's report). Returns, when the previous round's set is
/// not empty, whether each subject of the round (by subject index, a database's by rank) is in
/// it, for the two-part description table (`psiblast_round::ReportRound`, `FoundBefore`).
/// Refuses keys of the inclusion sets that NCBI's handles can take for one (`HandleIndex`).
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:529-566
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:158
// The Seq-id of a subject is its key (`rounds::SetKey`, `SubjectKeys`): a `-subject` record's
// index (`lcl|Subject_N`), a database sequence's OID (`gnl|BL_ORD_ID|oid`; the hit lists name it
// by its rank among the listed ones, `SubjectKeys` maps the rank to the OID), or its parsed Seq-align id. The PSSM's multiple
// alignment counts its rows by the same keys and joins following Seq-aligns whose ids `Match`
// (`pssm::input_msa`). Keys that differ but that NCBI's handles can take for one are refused
// where NCBI compares them (`rounds::HandleIndex`): in the inclusion sets before the round's
// report, a shown subject against the previous set before the outfmt 0 table (`FoundBefore`).
// NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:483-509 (cited at
// `core::psi::PsiScoreBlk::set_ungapped_statistical_params`: the round's PSI blocks)
/// Whether the PSSM takes a Seq-align of the e-value `e_value` (`GetLowestEvalue` of its scores,
/// below `-inclusion_ethresh`, `x_ExtractAlignmentData`): the round's inclusion set and the
/// subjects whose letters the PSSM reads (`ListedTitles::pssm_reads`).
pub(super) fn included(e_value: f64, inclusion_ethresh: f64) -> bool {
    const BLAST_EXPECT_VALUE: f64 = 10.0;
    crate::report::serial::seqalign::seqalign_evalue(e_value).min(BLAST_EXPECT_VALUE)
        < inclusion_ethresh
}

#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(super) fn round_outcome(
    args: &ResolvedBlastpArgs,
    io: &PsiRoundIo,
    hit_lists: &[Option<BlastpHitList>],
    ranges: &BlastpRanges<'_>,
    shown_subjects: ShownSubjects<'_>,
    contexts: &[QueryContext],
    gapped: &KarlinParams,
    keys: SubjectKeys<'_>,
) -> Result<Option<FoundBefore>> {
    let key = |s_idx: u32| keys.key(s_idx);
    let hsps: Vec<&crate::algorithm::blastp::hsp::BlastpHsp> = hit_lists
        .first()
        .and_then(Option::as_ref)
        .map(|hit_list| {
            let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
            hit_list.hsplist_array[..count]
                .iter()
                .flat_map(|hsp_list| &hsp_list.hsps)
                .collect()
        })
        .unwrap_or_default();
    let has_alignments = !hsps.is_empty();
    let previous = io.state.previously_found();
    let set: InclusionSet = hsps
        .iter()
        .filter(|hsp| included(hsp.e_value, args.inclusion_ethresh))
        .map(|hsp| key(hsp.s_idx))
        .collect();
    // The handles that NCBI compares: the set's keys with the previous set's (convergence,
    // `CPsiBlastIterationState::HasConverged`), the PSSM's rows among themselves (below, where the
    // round makes the next PSSM), and the shown subjects with the previous set in the two-part
    // table (`FoundBefore::unsure`, refused where it shows them). The previous set is one that
    // passed the rows' check: a round follows only one that made the next PSSM.
    let earlier = HandleIndex::new(previous)?;
    if let Some((other, key)) = set
        .iter()
        .find_map(|key| earlier.unsure(key).map(|other| (other, key)))
    {
        return Err(SetKey::refuse_pair(other, key));
    }
    let found_before = if previous.is_empty() {
        None
    } else {
        let subjects = hsps
            .iter()
            .map(|hsp| hsp.s_idx as usize + 1)
            .max()
            .unwrap_or(0);
        let mut found = vec![false; subjects];
        let mut unsure = std::collections::BTreeMap::new();
        for hsp in &hsps {
            let subject = key(hsp.s_idx);
            found[hsp.s_idx as usize] = previous.contains(&subject);
            if let Some(other) = earlier.unsure(&subject) {
                unsure.insert(hsp.s_idx, (other.clone(), subject));
            }
        }
        Some(FoundBefore { found, unsure })
    };
    // NCBI makes the next PSSM after `Advance` when another round follows (`itr`) or the last PSSM
    // is saved (`GetSaveLastPssm`), and after a round with alignments only.
    let makes_next = if has_alignments {
        let mut after = io.state.clone();
        after.advance(set.clone());
        after.go() || io.save_last
    } else {
        false
    };
    if makes_next {
        HandleIndex::new(&set)?;
    }
    let mut next_warnings = crate::algorithm::psiblast::pssm::PssmReadMessages::default();
    let next = makes_next.then(|| {
        next_pssm(
            args,
            io,
            &hsps,
            ranges,
            shown_subjects,
            contexts,
            gapped,
            keys,
            &mut next_warnings,
        )
    });
    // The blocks of a PSSM round's search, which NCBI writes back into the PSSM it searched
    // (`in_pssm_files::InPssmRounds::written_back`).
    let searched = args.psi_search.as_deref().map(|psi| {
        let karlin = super::pssm_hooks::psi_karlin(psi, contexts, gapped);
        (engine_block(&karlin.ungapped), engine_block(&karlin.gapped))
    });
    *io.outcome.lock().expect("the round's outcome") = Some(RoundOutcome {
        has_alignments,
        set,
        next,
        next_warnings,
        searched,
    });
    Ok(found_before)
}

/// What the two-part description table of a round after one with an inclusion set reads of the
/// previous set (`round_outcome`).
pub(super) struct FoundBefore {
    /// Whether each subject (by subject index) is in the previous round's set.
    pub(super) found: Vec<bool>,
    /// The subjects whose keys a key of the previous set differs from but can be one handle with
    /// (`HandleIndex::unsure`), with the two keys: refused where the table shows them
    /// (`refuse_unsure_shown`).
    pub(super) unsure: std::collections::BTreeMap<u32, (SetKey, SetKey)>,
}

impl FoundBefore {
    /// LOSAT's refusal of a subject of `shown` (subject indices of the outfmt 0 report) that the
    /// two-part table cannot place as NCBI does (`unsure`). `shown` is every reported subject, a
    /// superset of the first `-num_descriptions` Seq-aligns that NCBI's `x_SplitSeqAlign` looks up
    /// (blast_format.cpp:541-553); the checks of the shown titles refuse every such subject first
    /// in the cases built so far (review 2 of line DBID pssmrows, L1), so this stays a backstop.
    pub(super) fn refuse_unsure_shown(&self, shown: impl IntoIterator<Item = usize>) -> Result<()> {
        for s_idx in shown {
            if let Some((earlier, subject)) = u32::try_from(s_idx)
                .ok()
                .and_then(|s_idx| self.unsure.get(&s_idx))
            {
                return Err(SetKey::refuse_pair(earlier, subject));
            }
        }
        Ok(())
    }
}

/// The keys of the subjects of a search's hit lists (`rounds::SetKey`, by subject index; a
/// database's subjects by rank): the inclusion sets of psiblast's rounds and the rows of its
/// PSSMs (`round_outcome`, `next_pssm`, `psiblast_pssms`).
#[derive(Clone, Copy)]
pub(super) struct SubjectKeys<'a> {
    /// The ordinal ids of a database's listed subjects with the database's OID of each
    /// (`None` for `-subject` records).
    pub(super) listed_ordinals: Option<(&'a [u32], crate::seqdb::OidMap<'a>)>,
    /// The parsed Seq-ids of a database's subjects by rank (empty without them).
    pub(super) parsed: &'a [Option<Arc<ParsedSubject>>],
    /// What the PSSM reads by those Seq-align ids where it is not the subject's own letters
    /// (`ListedTitles::pssm_reads`; empty without parsed Seq-ids).
    pub(super) reads: &'a super::database::PssmReads,
}

impl SubjectKeys<'_> {
    /// The letters that the PSSM reads by the Seq-align id of subject `s_idx` (`PssmReads`), in
    /// NCBIstdaa: its own are `shown_subjects`' (`None` for no sequence). A subject whose own
    /// letters are missing or empty (no database or `-subject` record has none) is LOSAT's
    /// internal error, not NCBI's warning (review 2 of line DBID pssmread, L2).
    pub(super) fn letters(
        &self,
        shown_subjects: ShownSubjects<'_>,
        s_idx: u32,
    ) -> Result<Option<Vec<u8>>> {
        self.reads.letters(s_idx, || {
            let letters = match shown_subjects.letters(s_idx as usize) {
                // A subject's letters are never a query's Bioseq letters; they would read as these.
                Some(ShownLetters::Ascii(letters) | ShownLetters::Bioseq(letters)) => {
                    ncbistdaa_letters(letters)
                }
                Some(ShownLetters::Stored(residues)) => residues.to_vec(),
                None => Vec::new(),
            };
            if letters.is_empty() {
                bail!("internal error of LOSAT's PSIBLAST: subject {s_idx} of the PSSM has no letters");
            }
            Ok(letters)
        })
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:445-450
    /// NCBI's warning for each Seq-align whose letters the PSSM cannot read, by subject
    /// (`query_pssm`'s `unread`): the Seq-align id as a FASTA string, each with its subject. Only
    /// a database's subjects can fail to be read: one with parsed Seq-ids by that id, one without
    /// them by its ordinal id `gnl|BL_ORD_ID|<oid>` where its load throws (`PssmReads`).
    pub(super) fn unread_warnings(
        &self,
        program: &str,
        unread: &[u32],
        out: &mut crate::algorithm::psiblast::pssm::PssmReadMessages,
    ) {
        for &s_idx in unread {
            let id = match self.parsed.get(s_idx as usize).and_then(Option::as_deref) {
                Some(subject) => subject.alignment_fasta(),
                None => match self.listed_ordinals {
                    Some((ordinals, oids)) => {
                        format!("gnl|BL_ORD_ID|{}", oids.oid(ordinals[s_idx as usize]).get())
                            .into_bytes()
                    }
                    None => continue,
                },
            };
            let mut text = b"Failed to retrieve sequence ".to_vec();
            text.extend_from_slice(&id);
            // A post of NCBI's main thread (`ncbi_exception::post_line`).
            let line = crate::ncbi_exception::post_line(
                crate::ncbi_exception::Severity::Warning,
                None,
                Some(program),
                &text,
            );
            out.push(s_idx, line);
        }
    }
    pub(super) fn key(&self, s_idx: u32) -> SetKey {
        match self.listed_ordinals {
            None => SetKey::Index(s_idx),
            Some((ordinals, oids)) => {
                match self.parsed.get(s_idx as usize).and_then(Option::as_deref) {
                    Some(subject) => SetKey::Id(subject.alignment_id().clone()),
                    // The Seq-align id `gnl|BL_ORD_ID|<oid>` (`SetKey::ordinal_id`).
                    None => SetKey::Index(oids.oid(ordinals[s_idx as usize]).get()),
                }
            }
        }
    }
}

/// `ComputePssmForNextIteration` of the round's Seq-aligns `hsps` (`psiblast_pssms`' PSSM of the
/// FASTA query, its Seq-align id and title, which every round keeps) with the round's ancillary
/// data: round 1's standard blocks, or the PSSM's `kbp_psi` and `kbp_gap_psi` (and `kbp_gap`,
/// which is `kbp_gap_psi`) after a PSSM round.
// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:167-179
#[allow(clippy::too_many_arguments)]
fn next_pssm(
    args: &ResolvedBlastpArgs,
    io: &PsiRoundIo,
    hsps: &[&crate::algorithm::blastp::hsp::BlastpHsp],
    ranges: &BlastpRanges<'_>,
    shown_subjects: ShownSubjects<'_>,
    contexts: &[QueryContext],
    gapped: &KarlinParams,
    keys: SubjectKeys<'_>,
    warnings: &mut crate::algorithm::psiblast::pssm::PssmReadMessages,
) -> Result<QueryPssm> {
    let subject = |s_idx: u32| keys.letters(shown_subjects, s_idx);
    let scoring = PssmScoring {
        matrix: args.scoring.matrix,
        matrix_name: &args.matrix_name,
        gap_open: args.scoring.gap_open,
        gap_extend: args.scoring.gap_extend,
    };
    let q_shift = ranges.query_input.placements.offset(0);
    let alignments: Vec<PssmAlignment<'_>> = hsps
        .iter()
        .map(|hsp| PssmAlignment::of_hsp(hsp, q_shift))
        .collect();
    let ordinal = ranges.query_input.ordinal(0);
    let record = &ranges.input_records[ordinal];
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:148-158
    // The Seq-align's query id: a FASTA query's (`fasta_query_id`: `Query_<n>`, or the id that
    // `-parse_deflines` read), the PSSM query's own id with `-in_pssm` (a `local id` stays an
    // integer). The title is the query's (`s_GetQueryBioseq` gives the current PSSM's Bioseq, whose
    // title is the input PSSM's, `""` without one), and its letters the round-1 query's, O as X
    // (`encode_protein_sequence` in `query_pssm`).
    let id = match &io.query_id {
        Some(id) => id.clone(),
        None => crate::algorithm::psiblast::pssm::fasta_query_id(
            record,
            args.first_query_ordinal + ranges.x_ordinal_base + ordinal + 1,
        ),
    };
    let standard_ungapped = engine_block(&contexts[0].karlin_params);
    let table_gapped = engine_block(gapped);
    let ancillary = match args.psi_search.as_deref() {
        None => PssmAncillary {
            ungapped: standard_ungapped,
            gapped: table_gapped,
            psi_ungapped: standard_ungapped,
            psi_gapped: table_gapped,
        },
        Some(psi) => {
            let karlin = super::pssm_hooks::psi_karlin(psi, contexts, gapped);
            let psi_gapped = engine_block(&karlin.gapped);
            PssmAncillary {
                ungapped: standard_ungapped,
                gapped: psi_gapped,
                psi_ungapped: engine_block(&karlin.ungapped),
                psi_gapped,
            }
        }
    };
    let mut unread = Vec::new();
    let pssm = query_pssm(
        &io.request,
        scoring,
        table_gapped,
        args.inclusion_ethresh,
        PssmQueryRecord {
            letters: record.seq(),
            id: &id,
            title: record.title(),
        },
        &alignments,
        &subject,
        &|s_idx| keys.key(s_idx),
        ancillary,
        &mut unread,
        args.archive.is_some(),
    );
    keys.unread_warnings(args.app.ncbi_name(), &unread, warnings);
    pssm
}

#[cfg(test)]
#[path = "psi_round_tests.rs"]
mod tests;
