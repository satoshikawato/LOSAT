//! Three parts of a BLASTX run (`native::run_inner`), each called once at its place: the set-up of
//! the run's options, with the search strategy file and the requests of the BLAST archive and of
//! the search strategy (`RunSetup::open`), the read of the `-subject` file (`read_subjects`), and
//! the BLAST archive of one query batch (`ArchiveBatch::write`, outfmt 11).
//!
//! The NCBI reference for the set-up of the options (`CBlastAppArgs::SetOptions`, whose argument
//! handlers process the options and read the `-subject` file) is at `RunSetup::open`.
//!
//! This module is a child of `native` for the functions of that file it calls (`reported_error`,
//! `archive_error`, `batch_titles`, `ReportedAsks`). Each
//! function is the statements of `run_inner` at its place, in their order.

use super::{archive_error, batch_titles, reported_error, ReportedAsks};
use crate::algorithm::blastx::args::{
    BlastxArgs, ProgramName, ResolvedOptions, SemanticOptionsError,
};
use crate::algorithm::blastx::database::{BlastxDatabase, ReportedSubjects};
use crate::algorithm::blastx::input::FastaRecord;
use crate::algorithm::blastx::input_source::from_upper;
use crate::algorithm::blastx::runtime::BatchResults;
use crate::algorithm::blastx::subjects::Subjects;
use crate::algorithm::blastx::{archive, seqalign};
use crate::blastinput::app::{self, EarlyFormat, Program};
use crate::blastinput::fasta_reader;
use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::blastinput::query_batch::QueryScopes;
use crate::blastinput::seq_range::SequenceRange;
use crate::report::archive::requests::RunRequests;
use crate::report::archive::strategy::StrategyExport;
use crate::report::archive::ArchiveHandle;
use crate::report::pairwise::PairwiseDatabase;
use crate::seqid::id_scope::IdScope;
use anyhow::Result;
use std::{cell::RefCell, io::Write, sync::Arc};

/// What the set-up of the run keeps for the rest of `run_inner`: the `-outfmt` choice, the
/// requests of the BLAST archive (outfmt 11) and of the search strategy
/// (`-export_search_strategy`), which the caller's slots hold too, and the run's options
/// (`BlastxArgs::resolve`; a semantic error is raised later, where NCBI's handlers raise it).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3629-3630
pub(super) struct RunSetup {
    pub(super) choice: EarlyFormat,
    pub(super) archive: Option<ArchiveHandle>,
    pub(super) strategy: Option<StrategyExport>,
    pub(super) resolution: Result<ResolvedOptions>,
}

impl RunSetup {
    /// Open the search strategy file, parse `-outfmt`, make the requests of the BLAST archive and
    /// of the search strategy (written to `archive_slot` and `strategy_slot` at once), resolve the
    /// options, look for an ungapped search's matrix file and check `-num_threads`. An error of
    /// `resolve` that is not a semantic one is returned here.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3584-3585
    pub(super) fn open(
        args: &BlastxArgs,
        settings: &ApplicationSettings,
        program: ProgramName,
        archive_slot: &mut Option<ArchiveHandle>,
        strategy_slot: &mut Option<StrategyExport>,
    ) -> Result<Self> {
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3631-3639
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3624-3627
        // `ArchiveFormatRequested` parses `-outfmt` (blast_args.cpp:2745-2748) before the
        // option handlers run, with NCBI's errors; a format that NCBI runs and LOSAT does not
        // write is rejected there too. The custom tabular fields are read as NCBI reads them
        // (`ResolvedOptions::tabular`).
        // `GetExportSearchStrategyStream` comes first: it creates the search strategy file
        // (`strategy::open_strategy_file`), which an error after it leaves empty; with it a value
        // that does not parse is left to the formatting handler (`app::early_format_choice`).
        let strategy_file = crate::report::archive::strategy::open_strategy_file(
            args.export_search_strategy.as_deref(),
            program.upper,
        )?;
        let outfmt_bytes = crate::utils::os_bytes::os_str_bytes(&args.outfmt);
        let choice = app::early_format_choice(&outfmt_bytes, strategy_file.is_some())?;
        let format = choice.report_format(Program::Blastx, args.out.as_deref())?;
        // The BLAST archive (outfmt 11): its request, before anything is read or written.
        // The error archives follow the `-outfmt` text (`writes_error_archive`).
        // The run's options handle, which its archive and its search strategy share.
        let options_handle = crate::report::archive::OptionsHandle::default();
        let archive = archive::handle_for(
            args,
            format == Some(app::ReportFormat::ArchiveFormat),
            &options_handle,
        )
        .map(|handle| handle.with_outfmt_text(&outfmt_bytes));
        archive_slot.clone_from(&archive);
        // The search strategy's request, with the archive's options handle.
        let strategy = strategy_file.map(|file| {
            crate::report::archive::strategy::StrategyExport::for_args(
                args,
                file,
                choice.number(),
                &options_handle,
            )
        });
        strategy_slot.clone_from(&strategy);
        // Preserve explicit product refusals; defer supported semantic failures
        // until subject and standard query/output streams have been extracted.
        // A deferred `-outfmt` value is read again by the formatting handler (`at_formatting_handler`),
        // which fails with NCBI's error after the errors of the arguments that come before it.
        let mut resolution = match &choice {
            EarlyFormat::Deferred => Err(SemanticOptionsError::from(String::new()).into()),
            EarlyFormat::Parsed(_) => match args.resolve() {
                Err(error) if error.downcast_ref::<SemanticOptionsError>().is_none() => {
                    return Err(error)
                }
                resolution => resolution,
            },
        };
        // An ungapped search's `-matrix` that NCBI's tables do not have is looked for as a file when
        // NCBI sets up a query batch (`ApplicationSettings::matrix_file`): the first batch is refused
        // then (`query_setup::matrix_error_message`), and a run without a batch reads no file.
        if let Ok(options) = &mut resolution {
            options.find_matrix_file(&settings.matrix_search)?;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3157-3165
        // Existing LOSAT search-scoped pool validates the explicit thread contract
        // before input/output side effects; serial WASI remains explicitly serial.
        crate::utils::threading::validate_threads(args.num_threads as usize)?;
        Ok(RunSetup {
            choice,
            archive,
            strategy,
            resolution,
        })
    }
}

/// What the read of the `-subject` file reads: the run's arguments and settings (the data
/// loaders), the id scope of `-parse_deflines`, the open subject file (`None` with `-db`), the
/// requests of the BLAST archive and of the search strategy (which keep each record's definition
/// line), the run's standard error (the reader's warnings) and the function that writes an error's
/// archive (`early_archive_error`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2525-2527
pub(super) struct SubjectRead<'a> {
    pub(super) args: &'a BlastxArgs,
    pub(super) settings: &'a ApplicationSettings,
    pub(super) id_scope: &'a Option<Arc<IdScope>>,
    pub(super) subject_file: Option<crate::blastinput::input_files::InputFile>,
    pub(super) archive: &'a Option<ArchiveHandle>,
    pub(super) strategy: &'a Option<StrategyExport>,
    pub(super) stderr: &'a mut dyn Write,
    pub(super) early_archive_error: &'a dyn Fn(anyhow::Error) -> anyhow::Error,
}

/// Parse `-subject_loc` and read the records of the `-subject` file (none with `-db`); returns the
/// interval and the records, not yet cut to it.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2463-2465
pub(super) fn read_subjects(
    input: SubjectRead<'_>,
) -> Result<(Option<SequenceRange>, Vec<FastaRecord>)> {
    let SubjectRead {
        args,
        settings,
        id_scope,
        subject_file,
        archive,
        strategy,
        stderr,
        early_archive_error,
    } = input;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2540-2556
    // The handler of the database arguments reads `-subject_loc` after it opened the subject
    // file, and then the subjects, one record at a time: a record whose interval starts more
    // than one letter past its end stops the reading there (`seq_range::record_cut`), so the
    // reader hands over each record as it closes (a batch of 0 letters) and no warning of a
    // later record is written.
    let subject_range = crate::blastinput::seq_range::parse_optional_range(
        args.subject_loc.as_deref(),
        crate::blastinput::seq_range::RangeRole::Subject,
    )
    .map_err(early_archive_error)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:221-247
    // The subjects are read as the other programs' are (`fasta_reader::read_subjects_upper`), each
    // record's messages as it is read, and each is given its BLASTX form as it is read.
    let subjects = match subject_file {
        Some(subject_file) => {
            let mut source = fasta_reader::FastaInputSource::from_file(
                subject_file,
                fasta_reader::ReaderConfig::subject("BLASTX", true, settings.data_loaders)
                    .with_defline(archive.is_some() || strategy.is_some())
                    .with_parse_deflines(args.parse_deflines),
            );
            source.set_id_scope(id_scope.clone());
            fasta_reader::read_subjects_upper(
                &mut source,
                subject_range.as_ref(),
                args.lcase_masking,
                &mut |message| stderr.write_all(message),
                |record, masks| Ok(from_upper(record, masks)),
            )
            .map_err(early_archive_error)?
        }
        None => Vec::new(),
    };
    Ok((subject_range, subjects))
}

/// What the BLAST archive of a query batch (outfmt 11) reads of the run: the archive's handle and
/// requests, the database of `-db` and its pairwise description, the subjects the batch searched,
/// the options, the `-subject` records and the scopes of the queries (`-parse_deflines`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1067-1069
pub(super) struct ArchiveBatch<'a> {
    pub(super) handle: &'a ArchiveHandle,
    pub(super) archive_requests: &'a RunRequests<'a>,
    pub(super) database: &'a Option<BlastxDatabase>,
    pub(super) pairwise_database: &'a Option<PairwiseDatabase>,
    pub(super) search_subjects: Subjects<'a>,
    pub(super) options: &'a ResolvedOptions,
    pub(super) subjects: &'a Vec<FastaRecord>,
    pub(super) query_scopes: &'a RefCell<Option<QueryScopes>>,
}

impl ArchiveBatch<'_> {
    /// Write the BLAST archive of the query batch `queries` (after `queries_before` queries of the
    /// run) with its results `results` to `writer`: the hit lists sorted by e-value, with `-db`
    /// the titles of the reported sequences read first.
    ///
    /// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:283-285
    pub(super) fn write(
        self,
        writer: &mut crate::cli::NcbiOstream<'static>,
        diagnostics: &mut dyn Write,
        queries: &[FastaRecord],
        results: &mut Vec<BatchResults>,
        queries_before: usize,
    ) -> Result<()> {
        let ArchiveBatch {
            handle,
            archive_requests,
            database,
            pairwise_database,
            search_subjects,
            options,
            subjects,
            query_scopes,
        } = self;
        // The BLAST archive (outfmt 11): one document per query batch, written from the
        // data of outfmt 8's Seq-aligns, in the reports' order (`report::sort_lists`).
        crate::algorithm::blastx::report::sort_lists(
            results,
            crate::algorithm::blastx::empty_rows::SubjectLoc::of_run(database.is_some()),
        );
        let mut reported;
        let (subject_names, call_subjects): (_, &[FastaRecord]) =
            match (&database, &pairwise_database) {
                (Some(database), Some(pairwise)) => {
                    reported = match batch_titles(
                        database,
                        search_subjects,
                        results,
                        options.gapped,
                        // The archive reads no title besides the Seq-aligns' (no prefetch;
                        // RPSTBLASTN still stays `Reports`).
                        ReportedAsks::of_archive(options),
                        (&mut *diagnostics, Some(handle)),
                    )? {
                        None => ReportedSubjects::default(),
                        Some(titles) => ReportedSubjects::for_batch(
                            database,
                            search_subjects,
                            results,
                            titles,
                            options.gapped,
                            false,
                            false,
                            false,
                            false,
                            false,
                            false,
                            false,
                            crate::report::tabular::PrintedIds::default(),
                            false,
                            crate::report::db_subject::ReportRead {
                                one_by_one: false,
                                batch: Some(queries_before),
                                prefetch_fails: false,
                            },
                            usize::MAX,
                            (crate::report::db_subject::FormatThrow::Unread, None),
                        )
                        .map_err(reported_error)?,
                    };
                    if let Some(scopes) = query_scopes.borrow().as_ref() {
                        reported.refuse_ids_of_queries(
                            results,
                            options.gapped,
                            (
                                queries_before,
                                &crate::blastinput::query_batch::BatchQueries::Shown(scopes),
                            ),
                            &crate::report::query_bioseq::SeenRun::REFUSED,
                        )?;
                    }
                    (
                        seqalign::SeqalignSubjects::Database {
                            reported: &reported,
                            title: &pairwise.definition,
                        },
                        reported.records(),
                    )
                }
                _ => (seqalign::SeqalignSubjects::Records, &subjects[..]),
            };
        // The batch's search and the read of its definition lines have not failed: NCBI
        // writes its archive.
        archive::write_blastx_archive(
            &archive::ArchiveCall {
                handle,
                batch: archive::CallBatch {
                    options,
                    queries,
                    requests: archive_requests,
                    batches: results,
                },
                subjects: call_subjects,
                names: subject_names,
            },
            writer,
            diagnostics,
        )
        .map_err(archive_error)
    }
}
