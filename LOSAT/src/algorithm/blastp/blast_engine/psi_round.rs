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
//! before (`PsiSearch::from_engine`), whose query stays the round-1 query. deltablast's round 1 of
//! a batch runs the phases of its calls one at a time for all the batch's queries (their searches,
//! their prefetch, their reports: `search_delta_batch`, as `CDeltaBlast::Run` searches the whole
//! batch before its reports). The outfmt 0
//! prolog comes with the first call; every call writes no epilog, which the driver writes after
//! the last query (outfmt 7's count is that of the rounds). The PSSM files are written between
//! the calls, where NCBI writes them, after the report of the round before (with `-in_pssm` the
//! input PSSM before round 1, `algorithm::psiblast::in_pssm_files`).

use std::io::Write as _;
use std::sync::Arc;

use anyhow::{bail, Result};

use super::{
    blastp_database_name, prefetch_call, report_call, run_resolved_in_pool, search_call,
    write_pairwise_prologs, write_reports_without_queries, BlastpRanges, BlastpReports,
    BlastpSubjectSource, ShownLetters, ShownSubjects, SubjectSet,
};
use crate::algorithm::blastp::args::{BlastpApp, BlastpCompositionMode, ResolvedBlastpArgs};
use crate::algorithm::blastp::hsp::BlastpHitList;
use crate::algorithm::psiblast::pssm::{
    ncbistdaa_letters, query_pssm, seqalign_pssm_evalue, PssmAlignment, PssmAncillary,
    PssmQueryRecord, PssmScoring, QueryPssm,
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
        if let Some((warnings, error)) = batches.stop.as_ref() {
            // NCBI reference (598d8ae6): c++/src/app/blast/deltablast_app.cpp:288-329
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/deltablast.cpp:83-90
            // When NCBI already built a query factory for deltablast's first batch (`stop_batch`
            // is `Some`), `CDeltaBlast::Run` opens `-rpsdb` and checks the domain search options
            // before the query splitter (`CHUNK_SIZE`) and before `BlastSetup_Validate`
            // (`Sequence contains no data` / `Bad argument to function detected`).
            if let Some(delta) = args
                .delta
                .as_deref()
                .filter(|_| batches.stop_batch.is_some())
            {
                refuse_before_output(args, subject_source)?;
                outputs.diagnostics.write_all(warnings)?;
                if let Err(error) = delta.open_domain_side(&mut *outputs.diagnostics) {
                    if is_ncbi_stop(&error) {
                        write_pairwise_prologs(args, subject_source, subject_label, outputs)?;
                    }
                    return Err(error);
                }
                let stop = batches.stop.expect("stop").1;
                if is_ncbi_stop(&stop) {
                    write_pairwise_prologs(args, subject_source, subject_label, outputs)?;
                }
                return Err(stop);
            }
            // The request of a batch that NCBI made a query factory of, before any output; NCBI's
            // stream that throws while it writes the request (`encoder::NullElement`) ends the
            // run after the prolog (psiblast_app.cpp:506, 340).
            if let Err(error) = write_stop_request() {
                if error.is::<crate::cli::NativeError>() {
                    write_pairwise_prologs(args, subject_source, subject_label, outputs)?;
                }
                return Err(error);
            }
            write_pairwise_prologs(args, subject_source, subject_label, outputs)?;
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
                    args.app.ncbi_name(),
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

/// deltablast's refusals that the inputs decide before any search, made before the run writes
/// anything: a subject set (`-subject`) of more records than `-max_target_seqs` (the round-1
/// search, `RoundResults`), and a domain database of several volumes
/// (`DeltaRun::refuse_several_volumes`).
fn refuse_before_output(
    args: &ResolvedBlastpArgs,
    subject_source: BlastpSubjectSource<'_>,
) -> Result<()> {
    let subjects = subject_source.fasta_records().len();
    if subject_source.database().is_none() && subjects > 1 && subjects > args.max_target_seqs {
        bail!(
            "-subject with more records than the hit-list limit (-max_target_seqs, -num_descriptions, or -num_alignments: {}) is not supported by LOSAT's DELTABLAST (NCBI BLAST+ keeps the best subjects of the set in its first search)",
            args.max_target_seqs
        );
    }
    if let Some(delta) = args.delta.as_deref() {
        delta.refuse_several_volumes()?;
    }
    Ok(())
}

/// The composition-based modes of a run's rounds (`pssm_round_args`): the options' mode
/// (`original`), the mode that a PSSM search sets the options to (`pssm`,
/// `in_pssm::position_based_mode`), the mode of a round that searches a PSSM (`pssm_round`), and
/// the first round that searches one (`first_pssm_round`: deltablast's round 1, psiblast's
/// round 2).
#[derive(Clone, Copy)]
struct RoundModes {
    original: BlastpCompositionMode,
    pssm: CompositionMode,
    pssm_round: BlastpCompositionMode,
    first_pssm_round: u32,
}

/// What `search_delta_batch` reads of `run_queries`: each field is `run_queries`' value of the same
/// name.
#[derive(Clone, Copy)]
struct RoundsRun<'r> {
    args: &'r ResolvedBlastpArgs,
    plan: &'r PsiRoundPlan,
    query_records: &'r [FastaRecord],
    searched_queries: &'r [FastaRecord],
    query_input: &'r seq_range::QueryInput,
    batches: &'r QueryBatches,
    subjects: &'r SubjectSet<'r>,
    subject_source: BlastpSubjectSource<'r>,
    subject_label: &'r str,
    pool: &'r crate::utils::threading::SearchPool<'r>,
    modes: RoundModes,
}

/// The query input of the rounds of the searched query `query` of the run's `query_input`, whose
/// record is `record`: the query as round 1 searches it, in every round: with `-query_loc` its
/// range of the record (its placement; the PSSM is made of the whole record, `next_pssm`, and
/// searched at the range's offsets, `PsiSearch::from_engine`). A psiblast batch holds one searched
/// query, after the records that `-query_loc` skips (a range that starts past the end of its
/// record; NCBI's batch reader drops them and numbers the query by its place in the input,
/// `Query_<n>`): the query's record is its input ordinal's, not the batch's first.
fn round_query_input(
    query_input: &seq_range::QueryInput,
    query: usize,
    record: &FastaRecord,
) -> seq_range::QueryInput {
    let placement = query_input.placements.get(query);
    seq_range::QueryInput {
        placements: seq_range::Placements(placement.into_iter().collect()),
        ordinals: vec![0],
        input_lengths: vec![record.seq().len()],
        skipped: vec![false],
        range_length: query_input.range_length,
    }
}

/// The arguments of round `round` of a query whose record is `record`: the run's, and where the
/// round searches a PSSM (`pssm`) the search of that PSSM (`PsiSearch`), with its set-up's
/// warnings, the round's query (`round_query` holds the query of the round before, and takes this
/// round's where the next round starts from it) and the mode of a PSSM round. `calls`: the rounds
/// searched before in the run; `split`: whether NCBI splits the query (`RoundQuery::split`);
/// `cdd_query`: deltablast's query of its PSSM's Bioseq (`CddPssm::query`).
#[allow(clippy::too_many_arguments)]
fn pssm_round_args(
    args: &ResolvedBlastpArgs,
    modes: RoundModes,
    pssm: Option<&QueryPssm>,
    record: &FastaRecord,
    split: bool,
    round: u32,
    calls: usize,
    round_query: &mut Option<RoundQuery>,
    cdd_query: Option<&FastaRecord>,
) -> Result<ResolvedBlastpArgs> {
    let delta = args.delta.as_deref();
    let mut round_args = args.clone();
    let Some(current) = pssm else {
        return Ok(round_args);
    };

    // The options' mode when this round's score block is set up: the given one in
    // the query's first PSSM round, 1 (or 0) after it.
    let setup_mode = if delta.is_some() {
        // The options object is the run's: its first PSSM search resets the mode.
        if calls == 0 && round == 1 {
            app_mode(modes.original)
        } else {
            modes.pssm
        }
    } else if round == modes.first_pssm_round {
        app_mode(modes.original)
    } else {
        modes.pssm
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
        delta.is_some(),
    );
    let mut search = PsiSearch::from_engine(&current.engine, record.clone(), setup_messages);
    // The PSSM that the round's archive holds (made only for a run that writes one,
    // `query_pssm`'s `archive`).
    search.archive_pssm = args.archive.as_ref().and(current.object.clone());
    // The query source is fixed for the run: a PSSM query's warning label keeps its
    // title rule in every round (`PsiSearch::warning_title`).
    search.warning_title = args
        .psi_search
        .as_deref()
        .and_then(|psi| psi.warning_title.clone())
        // deltablast's warning has no query title.
        .or_else(|| delta.map(|_| Vec::new()));
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
    // And the letter past its Seq-inst's length (`PsiSearch::query_tail`).
    search.query_tail = args
        .psi_search
        .as_deref()
        .map_or_else(Vec::new, |psi| psi.query_tail.clone());
    // A query that NCBI splits: round 1's query, which the round searches, and the query
    // of its set-up (`RoundQuery::split`).
    let (query_now, split_setup) = if split {
        let (query, setup) = RoundQuery::split(args, record);
        (query, Some(setup))
    } else {
        let query = match (round_query.take(), cdd_query) {
            // NCBI reference (598d8ae6): c++/src/app/blast/deltablast_app.cpp:526-531
            // deltablast's rounds 1 and 2 are new searches (`CPsiBlast`) of the query
            // of the PSSM's Bioseq; from round 3 on `SetPssm` keeps round 2's query.
            (before, Some(cdd_query)) if round <= 2 || before.is_none() => {
                RoundQuery::first(args, cdd_query)
            }
            (Some(before), _) => before.next(args),
            (None, _) => RoundQuery::first(args, record).next(args),
        };
        (query, None)
    };
    // A later round whose letters give no Karlin-Altschul parameters (SEG again on the
    // letters that the rounds before masked can mask every letter left): NCBI's
    // setup crashes there, after this round's PSSM files, as for a PSSM query
    // (`invalid_pssm_query_refusal`).
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2782-2799
    // The set-up's query (`kbp_std`), which a split query's round reads too.
    let frames = [vec![split_setup
        .as_ref()
        .unwrap_or(&query_now.frame)
        .clone()]];
    if !super::pssm_hooks::psi_contexts(args, &frames)[0].is_valid {
        return Err(super::invalid_pssm_query_refusal(args.app));
    }
    search.round_query = Some(query_now.frame.clone());
    // A split query's next round starts again from the record (`RoundQuery::split`).
    if split_setup.is_none() {
        *round_query = Some(query_now);
    }
    search.split_setup = split_setup;
    round_args.psi_search = Some(Arc::new(search));
    round_args.comp_based_stats.mode = modes.pssm_round;
    Ok(round_args)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/deltablast.cpp:152-172
/// The results of a round, each searched by a call of its own (`RoundResults::call`): one, but in
/// deltablast's round 1 of a subject set (`-subject`). That round is `CPsiBl2Seq`'s search, whose
/// results are one for each subject, each reported as a query's report; the search reads the
/// statistics of the whole set. LOSAT searches each subject in a call of its own that reads those
/// statistics (`ResolvedBlastpArgs::whole_subject_set`). The preliminary stage of NCBI's one
/// search keeps the best `-max_target_seqs` subjects of the set only, which separate searches do
/// not: a larger set is LOSAT's refusal, made before the run writes anything
/// (`refuse_before_output`).
struct RoundResults<'s> {
    count: usize,
    whole_subject_set: Option<(i64, i32)>,
    subject_records: &'s [FastaRecord],
}

impl<'s> RoundResults<'s> {
    /// The results of round `round` against `subject_source`.
    fn of(args: &ResolvedBlastpArgs, round: u32, subject_source: BlastpSubjectSource<'s>) -> Self {
        let subject_records = subject_source.fasta_records();
        let count = if args.delta.is_some() && round == 1 && subject_source.database().is_none() {
            subject_records.len().max(1)
        } else {
            1
        };
        let whole_subject_set = (count > 1).then(|| {
            (
                subject_records
                    .iter()
                    .map(|record| record.seq().len() as i64)
                    .sum::<i64>(),
                i32::try_from(subject_records.len()).unwrap_or(i32::MAX),
            )
        });
        Self {
            count,
            whole_subject_set,
            subject_records,
        }
    }

    /// The call of result `index`: its arguments (the round's, `round_args`), its subjects and
    /// their placements (of `subjects`, the run's).
    fn call(
        &self,
        index: usize,
        round_args: &ResolvedBlastpArgs,
        subject_source: BlastpSubjectSource<'s>,
        subjects: &SubjectSet<'_>,
    ) -> (
        ResolvedBlastpArgs,
        BlastpSubjectSource<'s>,
        seq_range::Placements,
    ) {
        let mut call_args = round_args.clone();
        if self.count == 1 {
            return (call_args, subject_source, subjects.placements.clone());
        }
        call_args.whole_subject_set = self.whole_subject_set;
        (
            call_args,
            BlastpSubjectSource::Fasta(&self.subject_records[index..index + 1]),
            seq_range::Placements(subjects.placements.get(index).into_iter().collect()),
        )
    }
}

/// What the round driver hands a round's call (`PsiRoundIo`), whose state before the round is
/// `state`: `group_first` whether its query is its batch's first, `extra_result` whether the call
/// is not the round's first result of the batch (`PsiRoundIo::extra_result`), `cdd_query`
/// deltablast's query of its PSSM's Bioseq, `first_round_ancillary` round 1's ancillary data from
/// round 2 on.
fn round_io(
    args: &ResolvedBlastpArgs,
    plan: &PsiRoundPlan,
    state: &IterationState,
    group_first: bool,
    extra_result: bool,
    cdd_query: Option<&FastaRecord>,
    first_round_ancillary: Option<PssmAncillary>,
) -> PsiRoundIo {
    PsiRoundIo {
        state: state.clone(),
        save_last: plan.save_last && group_first,
        request: plan.request.clone(),
        query_id: args
            .psi_search
            .as_deref()
            .and_then(|psi| psi.query_seq_loc_id().cloned()),
        silent: args.delta.is_some() && plan.num_iterations == 1,
        extra_result,
        pssm_letters: cdd_query.map(|query| query.seq().to_vec()),
        first_round_ancillary,
        outcome: std::sync::Mutex::new(None),
    }
}

/// What a round's call of the searched query `query` reports besides its search
/// (`BlastpReports`): the reader's messages before its report (`reader_warnings`), the result sets
/// formatted before it (`calls`, NCBI's `m_QueriesFormatted`: BLAST XML numbers each round's
/// `Iteration` after them and writes its header with the first), whether its outfmt 0 report starts
/// with the prolog (`prolog`), its round (`io`), and NCBI's scope at its prefetch
/// (`prefetch_scope`).
#[allow(clippy::too_many_arguments)]
fn round_reports<'r>(
    reader_warnings: &'r [Vec<u8>],
    calls: usize,
    prolog: bool,
    io: &'r PsiRoundIo,
    batches: &'r QueryBatches,
    query: usize,
    searched_queries: &'r [FastaRecord],
    prefetch_scope: &'r (dyn Fn() -> crate::blastinput::query_batch::QueryScope + Sync),
) -> BlastpReports<'r> {
    BlastpReports {
        reader_warnings,
        identity_warned: 0,
        epilog: false,
        x_prolog: prolog,
        x_queries_before: calls,
        x_xml2_carry: None,
        last_call: false,
        psi_round: Some(io),
        stopped_pssms: None,
        query_scopes: batches
            .query_scopes
            .as_ref()
            .map(|scopes| scopes.numbered_from(query)),
        all_searched_queries: Some(searched_queries),
        prefetch_scope: Some(prefetch_scope),
        strategy_next: None,
        call_requests: None,
        x_plan_blocks: false,
    }
}

/// Whether `error` is an NCBI error that stops a deltablast batch after NCBI's `PrintProlog`
/// (`deltablast_app.cpp:288`), rather than a LOSAT refusal before any report is written.
fn is_ncbi_stop(error: &anyhow::Error) -> bool {
    error.is::<crate::cli::NativeError>()
        || error.is::<crate::algorithm::common::database::ErrorAfterPosts>()
}

/// The error of a failed round's call (`error`), with what NCBI writes before it. A header file
/// that NCBI cannot read fails the round whose results reach it first, at its prefetch, after the
/// reports of the rounds before (`search`'s rule for the batches): psiblast writes the outfmt 0
/// prolog if nothing was written (`calls` 0), and the reader's messages of a batch whose first
/// round fails (`reader_warnings`, with NCBI's posts of the bytes replaced in the query's request,
/// which round 1's report did not place); deltablast wrote both at the start of its batch. Then
/// NCBI's warning of its first read of a UTF8String header, where the round's Seq-aligns made it
/// before the failure (once per run: a round before posted its own).
#[allow(clippy::too_many_arguments)]
fn round_call_error(
    error: anyhow::Error,
    args: &ResolvedBlastpArgs,
    calls: usize,
    round: u32,
    reader_warnings: &[u8],
    subject_source: BlastpSubjectSource<'_>,
    subject_label: &str,
    outputs: &mut ReportOutputs<'_>,
) -> Result<anyhow::Error> {
    let failure = match error.downcast::<ErrorAfterBatchReports>() {
        Ok(failure) => failure,
        Err(error) => return Ok(error),
    };
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:262, 572-576
    // NCBI's warning of its first read of a UTF8String header, where the round's Seq-aligns made
    // it before the failure (once per run: a round before posted its own), then the posts of the
    // titles that its Seq-aligns read before the failure. Checked before the round writes its
    // prolog or warnings so a refusal of the failed round's title posts (e.g. both 803.10 and
    // 803.7 in one round) stops cleanly when nothing has been written yet.
    let database = subject_source.database();
    let utf8 = database.and_then(|database| {
        failure
            .utf8
            .as_ref()?
            .post(database.seqdb(), args.app.ncbi_name())
    });
    let titled = match database {
        Some(database) => database.failed_title_posts(&failure, utf8.is_some())?,
        None => Default::default(),
    };
    let delta = args.delta.is_some();
    if calls == 0 {
        write_pairwise_prologs(args, subject_source, subject_label, outputs)?;
    }
    if round == 1 && !delta {
        outputs.diagnostics.write_all(reader_warnings)?;
        super::write_request_posts(args, &mut *outputs.diagnostics)?;
    }
    for posts in utf8.iter().chain([&titled]) {
        posts.write(&mut *outputs.diagnostics, args.archive.as_ref())?;
    }
    crate::algorithm::common::database::ErrorAfterPosts::written(
        failure.error,
        &mut *outputs.diagnostics,
        args.archive.as_ref(),
    )
}

/// Round 1 of a query of a deltablast batch (`search_delta_batch`), for the query's later rounds
/// (`run_queries`): the outcome of its first result, and the query of its PSSM's Bioseq
/// (`CddPssm::query`).
struct RoundOne {
    outcome: RoundOutcome,
    cdd_query: FastaRecord,
}

/// A call of deltablast's round 1 of a batch (`search_delta_batch`), held from its search to its
/// report: what its phases read besides what the search hands on (`SearchedCall`). Each field is
/// `run_queries`' value of the same name for that call.
struct DeltaCall<'s> {
    query: usize,
    ordinal: usize,
    io: PsiRoundIo,
    call_source: BlastpSubjectSource<'s>,
    call_placements: seq_range::Placements,
    query_input: seq_range::QueryInput,
}

impl DeltaCall<'_> {
    /// The call's ranges (`BlastpRanges`), of its query's record in the run's `query_records`.
    fn ranges<'r>(&'r self, query_records: &'r [FastaRecord]) -> BlastpRanges<'r> {
        BlastpRanges {
            input_records: std::slice::from_ref(&query_records[self.ordinal]),
            query_input: &self.query_input,
            subject_placements: &self.call_placements,
            x_ordinal_base: self.ordinal,
            searched_base: self.query,
            split_search_spaces: None,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/deltablast.cpp:83-174
// NCBI reference (598d8ae6): c++/src/app/blast/deltablast_app.cpp:326-375
/// Round 1 of every query of a deltablast batch, the queries `batch_queries` (`CDeltaBlast::Run`,
/// then the application's steps after it): the domain search of every query, then for each query
/// its PSSM of its domains and its search, through the reads of its results' Seq-align ids
/// (`search_call`; NCBI's errors and warnings of those reads come there, before the next query's
/// PSSM); then the batch's first PSSM's files (`SavePssmToFile`), the batch's prefetch
/// (`prefetch_call`, each call's in the queries' order) and its reports (`report_call`). The round
/// is reported with round labels where rounds follow (`PsiRoundIo::silent`); its outcomes are each
/// query's first result's (`RoundOne`). `calls` counts the result sets formatted in the run.
#[allow(clippy::too_many_arguments)]
fn search_delta_batch(
    run: &RoundsRun<'_>,
    delta: &super::DeltaRun,
    batch_queries: std::ops::Range<usize>,
    prefetch_scope: &(dyn Fn() -> crate::blastinput::query_batch::QueryScope + Sync),
    files: &mut RoundFiles,
    carried_ascii: &mut Option<Vec<(&'static str, KarlinBlk)>>,
    calls: &mut usize,
    outputs: &mut ReportOutputs<'_>,
) -> Result<Vec<RoundOne>> {
    let RoundsRun {
        args,
        plan,
        query_records,
        searched_queries,
        query_input,
        batches,
        subjects,
        subject_source,
        subject_label,
        pool,
        modes,
    } = *run;
    // The query of each PSSM (`CddQuery`): its record and its Seq-align id.
    let alone = batch_queries.len() == 1;
    let ids: Vec<_> = batch_queries
        .clone()
        .map(|query| {
            let ordinal = query_input.ordinal(query);
            crate::algorithm::psiblast::pssm::fasta_query_id(
                &query_records[ordinal],
                args.first_query_ordinal + ordinal + 1,
            )
        })
        .collect();
    let cdd_queries: Vec<_> = (0..ids.len())
        .map(|index| super::delta_hooks::CddQuery {
            record: &query_records[query_input.ordinal(batch_queries.start + index)],
            id: &ids[index],
            alone,
        })
        .collect();
    // NCBI reference (598d8ae6): c++/src/app/blast/deltablast_app.cpp:288
    // Write the pairwise prolog before an NCBI error in batch 0; LOSAT refusals before the
    // batch's reports leave stdout empty.
    let with_prolog =
        |error: anyhow::Error, outputs: &mut ReportOutputs<'_>| -> Result<anyhow::Error> {
            if *calls == 0 && is_ncbi_stop(&error) {
                write_pairwise_prologs(args, subject_source, subject_label, outputs)?;
            }
            Ok(error)
        };
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/deltablast.cpp:89-90
    // One domain search of every query of the batch (`x_FindDomainHits`) comes before the first
    // query's PSSM.
    let domains = match delta.cdd_domains(&cdd_queries, pool, &mut *outputs.diagnostics) {
        Ok(domains) => domains,
        Err(error) => return Err(with_prolog(error, outputs)?),
    };
    let results = RoundResults::of(args, 1, subject_source);
    let state = IterationState::new(plan.num_iterations);
    let mut held = Vec::with_capacity(ids.len() * results.count);
    let mut cdd_queries_out = Vec::with_capacity(ids.len());
    let mut first_pssm = None;
    for (index, domains) in domains.into_iter().enumerate() {
        let query = batch_queries.start + index;
        let cdd = match delta.cdd_pssm(args, &plan.request, &cdd_queries[index], domains) {
            Ok(cdd) => cdd,
            Err(error) => return Err(with_prolog(error, outputs)?),
        };
        let mut current = cdd.pssm;
        // NCBI reference (598d8ae6): c++/src/app/blast/deltablast_app.cpp:328-353
        // NCBI saves the batch's first PSSM only (`SavePssmToFile` of the first query's), its
        // ASCII PSSM with the application's ancillary data (`m_AncillaryData`).
        if index > 0 {
            current.asn = None;
            current.ascii = None;
            current.asn_posts.clear();
        } else if let Some(object) = current.object.clone().filter(|_| current.ascii.is_some()) {
            let blocks = carried_ascii.get_or_insert_with(|| {
                crate::algorithm::psiblast::in_pssm_files::ascii_ancillary(&object).to_vec()
            });
            current.ascii = Some(crate::algorithm::psiblast::in_pssm_files::ascii_pssm_with(
                &object, blocks,
            ));
        }
        let record = &searched_queries[query];
        // Round 1's query is that of the PSSM's Bioseq; round 2 starts from it again
        // (`pssm_round_args`), so the round's query is not kept.
        let round_args = match pssm_round_args(
            args,
            modes,
            Some(&current),
            record,
            false,
            1,
            *calls + held.len(),
            &mut None,
            Some(&cdd.query),
        ) {
            Ok(round_args) => round_args,
            Err(error) => return Err(with_prolog(error, outputs)?),
        };
        for result_index in 0..results.count {
            let (call_args, call_source, call_placements) =
                results.call(result_index, &round_args, subject_source, subjects);
            let ordinal = query_input.ordinal(query);
            let call = DeltaCall {
                query,
                ordinal,
                io: round_io(
                    args,
                    plan,
                    &state,
                    index == 0,
                    index > 0 || result_index > 0,
                    Some(&cdd.query),
                    None,
                ),
                call_source,
                call_placements,
                query_input: round_query_input(query_input, query, &query_records[ordinal]),
            };
            // The library's records of the round (`FormatObserver::round_begin`).
            if let Some(observer) = outputs.observer.as_deref_mut() {
                observer.round_begin(query, 1);
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/deltablast.cpp:150-156
            // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:425
            // `CPsiBlast::Run` meets an empty `-db` or unopenable `.psq` here, after the domain
            // search and CDD PSSM, and before `SavePssmToFile` or report formatting.
            if let Some(failure) = call.call_source.database().and_then(|db| {
                crate::algorithm::common::database::PsiSearchReadFailure::of(db.seqdb())
            }) {
                return Err(with_prolog(failure.error(args.app.name()), outputs)?);
            }
            let formatted = *calls + held.len();
            let searched = search_call(
                call_args,
                std::slice::from_ref(record),
                call.call_source,
                subject_label,
                outputs,
                pool,
                None,
                &call.ranges(query_records),
                &round_reports(
                    &[],
                    formatted,
                    false,
                    &call.io,
                    batches,
                    query,
                    searched_queries,
                    prefetch_scope,
                ),
                super::x_lookup_ahead::XAheadBlock::none(),
            );
            let searched = match searched {
                Ok(searched) => {
                    searched.expect("a round's call is not of the preliminary stage only")
                }
                Err(error) => {
                    if error.is::<ErrorAfterBatchReports>() {
                        return Err(round_call_error(
                            error,
                            args,
                            *calls,
                            1,
                            &[],
                            subject_source,
                            subject_label,
                            outputs,
                        )?);
                    }
                    return Err(with_prolog(error, outputs)?);
                }
            };
            held.push((call, searched));
        }
        if index == 0 {
            first_pssm = Some(current);
        }
        cdd_queries_out.push(cdd.query);
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/deltablast_app.cpp:353,405
    // `SavePssmToFile` of the batch's first PSSM: after the searches of the batch
    // (`CDeltaBlast::Run` has returned), before its prefetch.
    if let Some(first) = first_pssm.as_ref() {
        if let Err(error) = write_round_pssm(args, files, first, outputs) {
            return Err(with_prolog(error, outputs)?);
        }
    }
    let mut prefetched = Vec::with_capacity(held.len());
    for (formatted, (call, searched)) in held.into_iter().enumerate() {
        let formatted = *calls + formatted;
        let record = std::slice::from_ref(&searched_queries[call.query]);
        let done = prefetch_call(
            searched,
            record,
            call.call_source,
            subject_label,
            outputs,
            pool,
            &call.ranges(query_records),
            &round_reports(
                &[],
                formatted,
                false,
                &call.io,
                batches,
                call.query,
                searched_queries,
                prefetch_scope,
            ),
        );
        match done {
            Ok((searched, done)) => prefetched.push((call, searched, done)),
            Err(error) => {
                if error.is::<ErrorAfterBatchReports>() {
                    return Err(round_call_error(
                        error,
                        args,
                        *calls,
                        1,
                        &[],
                        subject_source,
                        subject_label,
                        outputs,
                    )?);
                }
                return Err(with_prolog(error, outputs)?);
            }
        }
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/deltablast_app.cpp:288
    // Write the pairwise prolog once all searches and prefetches of batch 0 have passed
    // LOSAT's refusals; NCBI errors before this point write the prolog above.
    if *calls == 0 {
        write_pairwise_prologs(args, subject_source, subject_label, outputs)?;
    }
    let mut outcomes = Vec::with_capacity(cdd_queries_out.len());
    for (formatted, (call, searched, done)) in prefetched.into_iter().enumerate() {
        let formatted = *calls + formatted;
        let record = std::slice::from_ref(&searched_queries[call.query]);
        let reported = report_call(
            searched,
            done,
            record,
            call.call_source,
            outputs,
            pool,
            &call.ranges(query_records),
            &round_reports(
                &[],
                formatted,
                false,
                &call.io,
                batches,
                call.query,
                searched_queries,
                prefetch_scope,
            ),
        );
        if let Err(error) = reported {
            return Err(round_call_error(
                error,
                args,
                formatted.max(1),
                1,
                &[],
                subject_source,
                subject_label,
                outputs,
            )?);
        }
        let outcome = call
            .io
            .outcome
            .lock()
            .expect("the round's outcome")
            .take()
            .expect("every round's search reaches its outcome hook");
        // A query's first result is its round's outcome (`(*results)[0]` of its search).
        if call.query - batch_queries.start == outcomes.len() {
            outcomes.push(outcome);
        }
    }
    *calls += ids.len() * results.count;
    Ok(outcomes
        .into_iter()
        .zip(cdd_queries_out)
        .map(|(outcome, cdd_query)| RoundOne { outcome, cdd_query })
        .collect())
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
        args.app.ncbi_name(),
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
    // deltablast: the first PSSM search of the run is the one that warns of a mode that PSSM
    // searches reset (`delta_hooks`); its batches are NCBI's (`read_deltablast_query_batches`).
    let delta = args.delta.as_deref();
    let modes = RoundModes {
        original: original_mode,
        pssm: pssm_mode,
        pssm_round: pssm_round_mode,
        first_pssm_round: if delta.is_some() { 1 } else { 2 },
    };
    if delta.is_some() {
        refuse_before_output(args, subject_source)?;
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/deltablast_app.cpp:127-146
    // deltablast's `m_AncillaryData`: the Karlin-Altschul blocks that the ASCII PSSM of the first
    // query of a batch prints, taken from that PSSM where it is empty, set by each PSSM made from
    // a search's alignments, and emptied after each batch of a run with rounds only: a one-round
    // run prints a later batch's with the first batch's data (or those of the last PSSM that the
    // batch before saved).
    let mut carried_ascii: Option<Vec<(&'static str, KarlinBlk)>> = None;
    // The rounds formatted so far (`m_QueriesFormatted`: each round's `PrintOneResultSet`).
    let mut calls = 0usize;
    for (batch, &(input_start, first_query)) in batches.starts.iter().enumerate() {
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
        // NCBI reference (598d8ae6): c++/src/app/blast/deltablast_app.cpp:291-447
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/deltablast.cpp:83-174
        // A deltablast batch (several queries with one round, or `-num_iterations 0`): its prolog
        // and its reader's messages, then round 1 of every query of the batch, searched and
        // reported together (`search_delta_batch`); the later rounds are those of the batch's first
        // query (below). The first query's last PSSM of a one-round run comes after the batch's
        // reports (`batch_last`).
        let batch_queries = match delta {
            Some(_) => batches.searched_span(first_query),
            None => first_query..first_query + 1,
        };
        let round_ones = match delta {
            Some(delta) => {
                outputs
                    .diagnostics
                    .write_all(&batches.warnings[first_query])?;
                if plan.num_iterations != 1 {
                    carried_ascii = None;
                }
                let run = RoundsRun {
                    args,
                    plan,
                    query_records,
                    searched_queries,
                    query_input,
                    batches,
                    subjects,
                    subject_source,
                    subject_label,
                    pool,
                    modes,
                };
                search_delta_batch(
                    &run,
                    delta,
                    batch_queries.clone(),
                    &prefetch_scope,
                    &mut files,
                    &mut carried_ascii,
                    &mut calls,
                    outputs,
                )?
            }
            None => Vec::new(),
        };
        let mut round_ones = round_ones.into_iter();
        // A one-round deltablast batch's last PSSM (`-save_pssm_after_last_round`, of its first
        // query's round 1): its messages, the PSSM and the ancillary data it was made with.
        let mut batch_last = None;
        for query in batch_queries.clone() {
            let records = std::slice::from_ref(&searched_queries[query]);
            let ordinal = query_input.ordinal(query);
            let record = &query_records[ordinal];
            let placement = query_input.placements.get(query);
            let query_input = round_query_input(query_input, query, record);
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
            // deltablast's rounds search the PSSM's query, never split.
            let split = args.psi_search.is_none()
                && delta.is_none()
                && crate::algorithm::blastp::protein_query_split::split_protein_batch(
                    &[records[0].seq().len()],
                    crate::algorithm::blastp::protein_query_split::protein_query_split_sizes(
                        10000,
                    )?,
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
            // deltablast: round 1 was the batch's (`search_delta_batch`): its outcome, and the
            // query of its PSSM's Bioseq, that of the later PSSMs and of round 2
            // (`CddPssm::query`).
            let group_first = query == first_query;
            let (mut round_one, cdd_query) = match round_ones.next() {
                Some(RoundOne { outcome, cdd_query }) => (Some(outcome), Some(cdd_query)),
                None => (None, None),
            };
            // deltablast: round 1's ancillary data, whether round 1 had alignments, and the last
            // PSSM made for a round without alignments (`RoundOutcome::last`).
            let mut first_round_ancillary = None;
            let mut first_round_alignments = false;
            let mut last_candidate = None;
            while state.go() {
                let round = state.number();
                let outcome = match round_one.take() {
                    // deltablast's round 1, searched and reported with its batch.
                    Some(outcome) => outcome,
                    None => {
                        match (&pssm, input_files.as_mut()) {
                            (Some(current), _) => {
                                write_round_pssm(args, &mut files, current, outputs)?
                            }
                            (None, Some(input)) if round == 1 => {
                                // NCBI's posts of the bytes replaced in the search strategy request
                                // (written above, before round 1) come before those of the PSSM
                                // file that round 1 writes.
                                super::write_request_posts(args, &mut *outputs.diagnostics)?;
                                write_file_posts(args, &input.posts, outputs)?;
                                // NCBI's prolog comes before the rounds (`PrintProlog`), so before
                                // a file that cannot be written, or NCBI's throw in it, stops the
                                // run here.
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
                                        )?;
                                    }
                                    return Err(error);
                                }
                                // Only the ASCII PSSM is written again (after a round 1 without
                                // alignments; the ASN.1 text is then made with the write-back).
                                input.asn = None;
                            }
                            _ => {}
                        }
                        let round_args = pssm_round_args(
                            args,
                            modes,
                            pssm.as_ref(),
                            &records[0],
                            split,
                            round,
                            calls,
                            &mut round_query,
                            cdd_query.as_ref(),
                        )?;
                        let results = RoundResults::of(args, round, subject_source);
                        let mut first_outcome: Option<RoundOutcome> = None;
                        for result_index in 0..results.count {
                            let io = round_io(
                                args,
                                plan,
                                &state,
                                group_first,
                                result_index > 0,
                                cdd_query.as_ref(),
                                first_round_ancillary.filter(|_| round >= 2),
                            );
                            // The library's records of the round (`FormatObserver::round_begin`):
                            // the call's hit list and HSP numbers are the round's.
                            if let Some(observer) = outputs.observer.as_deref_mut() {
                                observer.round_begin(query, round);
                            }
                            // The reader's messages of the batch come with its first round
                            // (deltablast's before its domain search).
                            let reader_warnings: &[Vec<u8>] = if round == 1 && delta.is_none() {
                                &batches.warnings[query..query + 1]
                            } else {
                                &[]
                            };
                            let (call_args, call_source, call_placements) =
                                results.call(result_index, &round_args, subject_source, subjects);
                            let call_ranges = BlastpRanges {
                                subject_placements: &call_placements,
                                ..ranges
                            };
                            let result = run_resolved_in_pool(
                                call_args,
                                records,
                                call_source,
                                "",
                                subject_label,
                                outputs,
                                pool,
                                None,
                                &call_ranges,
                                &round_reports(
                                    reader_warnings,
                                    calls,
                                    calls == 0 && delta.is_none(),
                                    &io,
                                    batches,
                                    query,
                                    searched_queries,
                                    &prefetch_scope,
                                ),
                                super::x_lookup_ahead::XAheadBlock::none(),
                            );
                            if let Err(error) = result {
                                return Err(round_call_error(
                                    error,
                                    args,
                                    calls,
                                    round,
                                    &batches.warnings[query],
                                    subject_source,
                                    subject_label,
                                    outputs,
                                )?);
                            }
                            calls += 1;
                            let outcome = io
                                .outcome
                                .lock()
                                .expect("the round's outcome")
                                .take()
                                .expect("every round's search reaches its outcome hook");
                            if result_index == 0 {
                                first_outcome = Some(outcome);
                            }
                        }
                        first_outcome.expect("a round searches at least one result")
                    }
                };
                // NCBI reference (598d8ae6): c++/src/app/blast/deltablast_app.cpp:507-516
                // deltablast's later rounds are the batch's first query's (`(*results)[0]`): the
                // other queries of the batch end with round 1.
                if !group_first {
                    break;
                }
                if round == 1 {
                    first_round_ancillary = outcome.ancillary;
                }
                if !outcome.has_alignments {
                    if round == 1 {
                        searched_input = outcome.searched;
                    }
                    break;
                }
                first_round_alignments |= round == 1;
                last_candidate = outcome.last;
                state.advance(outcome.set);
                if delta.is_some() && plan.num_iterations == 1 {
                    // NCBI reference (598d8ae6): c++/src/app/blast/deltablast_app.cpp:384-396
                    // One round: NCBI makes the last PSSM after the batch's reports.
                    if group_first {
                        batch_last = Some((outcome.next_warnings, outcome.next, outcome.ancillary));
                    }
                } else if state.go() || plan.save_last {
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
            // NCBI reference (598d8ae6): c++/src/app/blast/deltablast_app.cpp:507-510,569-577
            // deltablast with rounds saves no last PSSM after a round 1 without alignments; after a
            // later round without them it saves the PSSM of the alignments before with round 1's
            // ancillary data (`RoundOutcome::last`; from round 1's alignments that is the PSSM that
            // round searched). One round saves it after the batch (`batch_last`).
            if plan.save_last && group_first && delta.is_some() {
                if plan.num_iterations != 1 && first_round_alignments {
                    match (last_candidate, &pssm) {
                        (Some((last, warnings)), _) => {
                            outputs.diagnostics.write_all(&warnings.bytes())?;
                            write_round_pssm(args, &mut files, &last?, outputs)?;
                        }
                        (None, Some(current)) => {
                            write_round_pssm(args, &mut files, current, outputs)?
                        }
                        (None, None) => {}
                    }
                }
            } else if plan.save_last && group_first {
                match (&pssm, &plan.input, &input_files, searched_input) {
                    (Some(current), _, _, _) => {
                        write_round_pssm(args, &mut files, current, outputs)?
                    }
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
        }
        // NCBI reference (598d8ae6): c++/src/app/blast/deltablast_app.cpp:384-396
        // A one-round deltablast batch's last PSSM, from its first query's round 1 after the
        // batch's reports, with that round's ancillary data, which `m_AncillaryData` keeps.
        if delta.is_some() && plan.num_iterations == 1 && plan.save_last {
            let Some((warnings, Some(last), ancillary)) = batch_last else {
                bail!(
                    "-save_pssm_after_last_round with one round where the first of the queries read together (up to 10000 letters) has no alignments (NCBI BLAST+ stops with an exception for a null pointer) is not supported by LOSAT's DELTABLAST"
                );
            };
            outputs.diagnostics.write_all(&warnings.bytes())?;
            let last = last?;
            if let Some(ancillary) = ancillary {
                carried_ascii = Some(crate::algorithm::psiblast::pssm::ancillary_blocks(
                    &ancillary,
                ));
            }
            write_round_pssm(args, &mut files, &last, outputs)?;
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
// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_args.hpp:1135-1145
/// The line `Search has CONVERGED!` between two `NcbiEndl`s, for every report that is not
/// structured (`HasStructuredOutputFormat`: BLAST XML, the Seq-annots, the XML2 and JSON documents
/// and SAM) and not the archive (psiblast's outfmt 0 to 4, 6, 7, 10, 18 and 20), on the report's
/// one stream: the application writes it after the query's rounds, outside the formatter's regions
/// (`Site::Report`: no `CBlastTabularInfo` lives there, so a failed write there is the
/// application's "BLAST failed to write output", exit 6), each `NcbiEndl` a line feed and a flush
/// (`ostream::ncbi_endl`). The next query's report starts in the format's region again
/// (`Site::of_format`).
fn write_converged(outputs: &mut ReportOutputs<'_>) -> Result<()> {
    use crate::blastinput::app::ReportFormat;
    use crate::report::ostream::{ncbi_endl, ReportWrite, Site};
    for format in outputs.formats.iter_mut() {
        let number = crate::blastinput::app::parse_formatting_string(format.outfmt)?.number;
        if !matches!(
            ReportFormat::from_number(number),
            None | Some(
                ReportFormat::Xml
                    | ReportFormat::AsnText
                    | ReportFormat::AsnBinary
                    | ReportFormat::Xml2
                    | ReportFormat::Json
                    | ReportFormat::Xml2S
                    | ReportFormat::JsonS
                    | ReportFormat::JsonSeqalign
                    | ReportFormat::Sam
                    | ReportFormat::ArchiveFormat
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
/// The epilog of each report after the last round: outfmt 0's (also that of the query-anchored
/// reports 1 to 4), outfmt 7's `# BLAST processed R queries` (R the rounds formatted,
/// `m_QueriesFormatted`), BLAST XML's end, nothing for outfmt 6, 8 to 12, 18 and 20.
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2237-2238
            // The formats below the tabular ones end with the pairwise footer; the tabular ones and
            // those after them write no footer.
            0..=4 => {}
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
            database_name: blastp_database_name(args, subject_label),
            database_num_sequences: num_sequences,
            database_total_letters: total_letters,
            matrix_name: Some(args.matrix_name.clone()),
            gap_open: args.scoring.gap_open,
            gap_extend: args.scoring.gap_extend,
            word_threshold: args.threshold,
            window_size: args.window_size,
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
            args.html.is_some(),
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
    seqalign_pssm_evalue(e_value) < inclusion_ethresh
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
    let earlier = HandleIndex::new(previous, args.app)?;
    if let Some((other, key)) = set
        .iter()
        .find_map(|key| earlier.unsure(key).map(|other| (other, key)))
    {
        return Err(SetKey::refuse_pair(other, key, args.app));
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
    // NCBI reference (598d8ae6): c++/src/app/blast/deltablast_app.cpp:466-580
    // A later subject's result of deltablast's round 1 of `-subject` makes none: the round's
    // results are one for each subject and its PSSM is made of the first's (`extra_result`).
    let goes_on = has_alignments && !io.extra_result && {
        let mut after = io.state.clone();
        after.advance(set.clone());
        after.go()
    };
    let makes_next = has_alignments && !io.extra_result && (goes_on || io.save_last);
    if makes_next {
        HandleIndex::new(&set, args.app)?;
    }
    let ancillary = round_ancillary(args, contexts, gapped);
    // NCBI reference (598d8ae6): c++/src/app/blast/deltablast_app.cpp:569-577
    // deltablast's last saved PSSM after round 2 or later is made with round 1's ancillary data:
    // here where no round follows, and, where one follows, as the PSSM saved if that round has
    // no alignments (`RoundOutcome::last`).
    let next_ancillary = match io.first_round_ancillary {
        Some(first) if !goes_on => first,
        _ => ancillary,
    };
    let mut next_warnings = crate::algorithm::psiblast::pssm::PssmReadMessages::default();
    let pssm_of =
        |ancillary: PssmAncillary,
         warnings: &mut crate::algorithm::psiblast::pssm::PssmReadMessages| {
            next_pssm(
                args,
                io,
                &hsps,
                ranges,
                shown_subjects,
                gapped,
                keys,
                ancillary,
                warnings,
            )
        };
    let next = makes_next.then(|| pssm_of(next_ancillary, &mut next_warnings));
    let last = io
        .first_round_ancillary
        .filter(|_| makes_next && goes_on && io.save_last)
        .map(|first| {
            let mut warnings = crate::algorithm::psiblast::pssm::PssmReadMessages::default();
            let pssm = pssm_of(first, &mut warnings);
            (pssm, warnings)
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
        ancillary: args.delta.is_some().then_some(ancillary),
        last,
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
    /// `app`: the run's program, which the refusal names.
    pub(super) fn refuse_unsure_shown(
        &self,
        shown: impl IntoIterator<Item = usize>,
        app: BlastpApp,
    ) -> Result<()> {
        for s_idx in shown {
            if let Some((earlier, subject)) = u32::try_from(s_idx)
                .ok()
                .and_then(|s_idx| self.unsure.get(&s_idx))
            {
                return Err(SetKey::refuse_pair(earlier, subject, app));
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
                bail!("LOSAT internal error: subject {s_idx} of the PSSM has no letters");
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

/// The round's ancillary data (`CBlastAncillaryData` of its results): round 1's standard blocks,
/// or the PSSM's `kbp_psi` and `kbp_gap_psi` (and `kbp_gap`, which is `kbp_gap_psi`) after a PSSM
/// round.
fn round_ancillary(
    args: &ResolvedBlastpArgs,
    contexts: &[QueryContext],
    gapped: &KarlinParams,
) -> PssmAncillary {
    let standard_ungapped = engine_block(&contexts[0].karlin_params);
    let table_gapped = engine_block(gapped);
    match args.psi_search.as_deref() {
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
    }
}

/// `ComputePssmForNextIteration` of the round's Seq-aligns `hsps` (`psiblast_pssms`' PSSM of the
/// FASTA query, its Seq-align id and title, which every round keeps; deltablast's letters are the
/// domains' PSSM's, `PsiRoundIo::pssm_letters`) with the ancillary data `ancillary`
/// (`round_ancillary`).
// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:167-179
#[allow(clippy::too_many_arguments)]
fn next_pssm(
    args: &ResolvedBlastpArgs,
    io: &PsiRoundIo,
    hsps: &[&crate::algorithm::blastp::hsp::BlastpHsp],
    ranges: &BlastpRanges<'_>,
    shown_subjects: ShownSubjects<'_>,
    gapped: &KarlinParams,
    keys: SubjectKeys<'_>,
    ancillary: PssmAncillary,
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
    let table_gapped = engine_block(gapped);
    let mut unread = Vec::new();
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:240-261
    // The PSSM engine reads the query's whole Seq-data, also the letter past a Seq-inst `length`
    // that is one less than the letters (`PsiSearch::query_tail`).
    let query_letters = [
        record.seq(),
        args.psi_search
            .as_deref()
            .map_or(&[][..], |psi| psi.query_tail.as_slice()),
    ]
    .concat();
    let pssm = query_pssm(
        &io.request,
        scoring,
        table_gapped,
        args.inclusion_ethresh,
        PssmQueryRecord {
            letters: io.pssm_letters.as_deref().unwrap_or(&query_letters),
            id: &id,
            title: record.title(),
        },
        &alignments,
        &subject,
        &|s_idx| keys.key(s_idx),
        ancillary,
        &mut unread,
        args.archive.is_some(),
        args.app,
    );
    keys.unread_warnings(args.app.ncbi_name(), &unread, warnings);
    pssm
}

#[cfg(test)]
#[path = "psi_round_tests.rs"]
mod tests;
