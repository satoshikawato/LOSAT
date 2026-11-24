//! TBLASTN's application (NCBI c++/src/app/blast/tblastn_app.cpp:114-361, `Run` and
//! `x_RunMTBySplitDB`): the command line's `TblastnArgs::run`, web ABI v2's `run_local`, and the
//! search over the query batches that both reach. The options are in `args.rs`.
use anyhow::{Context, Result};
use std::io::Write;

use super::args::{search_settings, ResolvedTblastnArgs, SearchSettings, TblastnArgs};
use super::database::TblastnDatabase;
use super::stage_d_pipeline::{
    run_local_for_report_threads, LocalStageDProfile, SubjectSetNumbers,
};
use super::stage_d_pipeline::{ReadFails, SetupEngineError};
use super::stage_d_stats::LocalSubjectParameters;
use super::stage_e_report::{render, ShownSubjects};
use super::x_subject_store::SubjectLetters;
use crate::algorithm::blastn::blast_engine::QueryReading;
use crate::algorithm::common::database::PreliminaryReadError;
use crate::algorithm::common::x_db_stream::XReadError;
use crate::api::local_blast::{OutputSink, ReportOutputs};
use crate::api::QueryEnd;
use crate::blastinput::app::ErrorAfterBatchReports;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::core::composition_adjustment::redo_alignment::BlastCompoAdjustMode;
use crate::report::tabular::SubjectFetch;

/// Whether `error` may fail one batch, after the reports of the batches before, rather than the
/// run: the errors `batch_failure` reads (a read of the database that failed, `XReadError`; a
/// batch whose setup failed, `SetupEngineError`; or one whose search LOSAT refuses or fails after
/// the reports of the batches before, `ErrorAfterBatchReports`: a step where NCBI's search crashes,
/// does not finish or stops, and the failed or endless Smith-Waterman redo of `-use_sw_tback`). The
/// batches searched ahead (`stage_d_pipeline::x_search_batches`) are then searched again one at a
/// time, so that the earliest batch that fails does in its place, after the reports of the batches
/// before it, and nothing of it or of the batches after it is written (`batch_failure` decides).
fn fails_one_batch(error: &anyhow::Error) -> bool {
    error.is::<XReadError>()
        || error.is::<SetupEngineError>()
        || error.is::<ErrorAfterBatchReports>()
}

/// The failure of a batch whose search failed where NCBI's fails it with its engine error: a
/// read of the database (`common::x_db_stream::XReadError`), as the check's error fails the first
/// batch with a valid context, when the reader's error is NCBI's own, and the batch's setup
/// (`stage_d_pipeline::SetupEngineError`: a PSSM query whose `kbp_psi` is not valid, the average
/// subject length); and where LOSAT refuses a step at which NCBI's search of the batch crashes,
/// does not finish or stops, after it has printed the reports of the batches before
/// (`ErrorAfterBatchReports`: line RG item X-short-bat, and the setup's statistics count of 0,
/// `stage_d_pipeline::preliminary_setup_checks`, D1a-5 of line RG), or the batch's search fails
/// after those reports (the failed or endless Smith-Waterman redo of `-use_sw_tback`,
/// `match_redo`): that error, as the batch's failure. Any other refusal of LOSAT's (and any other
/// error) is returned as `Err` to stop the run before its first byte.
fn batch_failure(error: anyhow::Error) -> std::result::Result<anyhow::Error, anyhow::Error> {
    if let Some(SetupEngineError(text)) = error.downcast_ref::<SetupEngineError>() {
        return Ok(crate::blastinput::app::engine_error(text));
    }
    let error = match error.downcast::<ErrorAfterBatchReports>() {
        Ok(ErrorAfterBatchReports(refusal)) => return Ok(refusal),
        Err(error) => error,
    };
    let error = match error.downcast::<XReadError>() {
        Ok(XReadError(read)) => read,
        Err(error) => return Err(error),
    };
    match error.downcast::<PreliminaryReadError>() {
        Ok(read_error) => Ok(crate::blastinput::app::engine_error(read_error.0)),
        Err(error) => Err(error),
    }
}

/// The first read of the database `database` (`-db`), where NCBI reads it, in the preliminary
/// search of the first batch with a valid context (`search`): its stream, or the read's error (a
/// database without sequences, a Gumbel block's shortest length that leaves no residue, a sequence
/// file that cannot be read; with `-in_pssm`, SeqDB's own error), and whether that error comes
/// before the batch's setup checks (`PsiSearchReadFailure::precedes_setup_checks`). The error
/// comes with the step of NCBI's search of a split batch that meets it (`ReadFails`, decided once,
/// here): the setup's read of the length of OID 0 (`Setup`, before the chunk's Karlin-Altschul
/// check), else the search's first read of a sequence (`Search`, after it; also the subject-length
/// error, which the chunk's own setup checks meet first). One owner for the search and the set-up
/// of a refused batch (`refused_batch_stop`).
fn first_database_read<'d>(
    database: &'d TblastnDatabase,
    resolved: &ResolvedTblastnArgs,
    subject_set: &SubjectSetNumbers,
) -> (
    std::result::Result<
        crate::algorithm::common::x_db_stream::XDbStream<'d>,
        (anyhow::Error, ReadFails),
    >,
    bool,
) {
    // `-in_pssm` (PSI-TBLASTN): `CPsiBlast` reads the database without `CLocalBlast`'s
    // catch, so its first read gives SeqDB's own error (`OID not found`, a sequence file
    // not found; `shared::PsiSearchReadFailure`, psiblast's rule), where a FASTA query
    // gets `Database memory map file error`.
    let pssm_failure = resolved
        .pssm
        .as_ref()
        .and_then(|_| database.pssm_search_read_failure());
    let before_setup = pssm_failure
        .as_ref()
        .is_some_and(|failure| failure.precedes_setup_checks());
    if database.first_oid_read_fails() {
        let error = match pssm_failure {
            Some(failure) => failure.error("TBLASTN"),
            None => PreliminaryReadError("Database memory map file error").into(),
        };
        return (Err((error, ReadFails::Setup)), before_setup);
    }
    let read = if !subject_set.old_fsc
        && crate::stats::setup::setup_subject_length(
            0,
            1,
            Some((database.numbers().min_length, true)),
        )
        .is_none()
    {
        Err(PreliminaryReadError(crate::stats::setup::SUBJECT_LENGTH_TOO_SHORT).into())
    } else if let Some(failure) = pssm_failure {
        Err(failure.error("TBLASTN"))
    } else {
        database.x_stream(resolved.num_threads, resolved.db_mask.as_ref())
    };
    (
        read.map_err(|error| (error, ReadFails::Search)),
        before_setup,
    )
}

/// NCBI's set-up of the records of a refused batch that are first in their scope (`first`, by
/// input record from `input_done`; the batch's searched queries start at `searched_done`): its
/// error where NCBI's set-up of the batch fails for the subjects
/// (`stage_d_pipeline::set_up_batch`, then the first read of the database where the batch's search
/// gets to it, as the search of a batch fails in `search`), in the one form of every program
/// (`query_batch::refused_batch_error`; line WQA D10 round 5, 2026-11-20); `Searched` where it
/// passes with a valid context, so that the refusal comes unless the run's XML failure takes it;
/// `Unknown` without a valid context or where LOSAT's error is not NCBI's. The batch that NCBI
/// splits holds the repeated ids too, by the lengths of their own records (line RG item
/// X-short-bat, review 4 M1). NCBI's errors of the batches before come first: their searches fail
/// before the stop is acted on. `queries`: the input records, the searched ones and how they map
/// (`QueryInput`).
fn refused_batch_setup_error(
    first: &[bool],
    (input_done, searched_done): (usize, usize),
    queries: (
        &[FastaRecord],
        &[FastaRecord],
        &crate::blastinput::seq_range::QueryInput,
    ),
    resolved: &ResolvedTblastnArgs,
    subject_set: SubjectSetNumbers,
    engine_positives: bool,
    database: Option<&TblastnDatabase>,
    subject_letters: SubjectLetters<'_>,
    (compo_adjust_mode, scoring): (
        crate::core::composition_adjustment::redo_alignment::BlastCompoAdjustMode,
        super::stage_d_pipeline::LocalStageDScoring,
    ),
) -> crate::blastinput::query_batch::BatchSetUp {
    use crate::blastinput::query_batch::BatchSetUp;
    let (input_records, searched_queries, query_input) = queries;
    // NCBI's batch (`stage_d_pipeline::BatchRecord`): each searched query by its own length, with
    // the index of its letters among those set up where it is first in its scope.
    let mut queries: Vec<usize> = Vec::new();
    let mut records: Vec<super::stage_d_pipeline::BatchRecord> = Vec::new();
    for query in (searched_done..searched_queries.len())
        .take_while(|&query| query_input.ordinal(query) < input_done + first.len())
    {
        let own = first[query_input.ordinal(query) - input_done];
        records.push((
            searched_queries[query].seq().len(),
            own.then_some(queries.len()),
        ));
        if own {
            queries.push(query);
        }
    }
    let letters: Vec<Vec<u8>> = queries
        .iter()
        .map(|&query| searched_queries[query].seq().to_vec())
        .collect();
    let whole: Vec<(&[u8], usize)> = queries
        .iter()
        .map(|&query| {
            (
                input_records[query_input.ordinal(query)].seq(),
                query_input.placements.offset(query),
            )
        })
        .collect();
    // NCBI's first read of the database: where the batch's search gets to it (`set_up_batch`, whose
    // chunks meet a read that fails where NCBI's do), its error.
    let read = database.map(|database| first_database_read(database, resolved, &subject_set).0);
    let subject_set = SubjectSetNumbers {
        read_fails: match &read {
            Some(Err((_, class))) => Some(*class),
            _ => None,
        },
        ..subject_set
    };
    let seg = resolved.seg.params();
    let profile = LocalStageDProfile {
        seg: seg.as_ref(),
        soft_masking: resolved.soft_masking,
        mask_lowercase: resolved.lcase_masking,
        genetic_code: resolved.db_gencode,
        expect_value: resolved.evalue,
        max_target_seqs: resolved.hitlist_size,
        query_records: &whole,
        subject_set,
        engine_positives,
        culling_limit: resolved.culling_limit,
        best_hit: resolved.best_hit,
        subject_besthit: resolved.subject_besthit,
        pssm: None,
        fetch_excluded: None,
        max_intron_length: resolved.max_intron_length,
        query_cov_hsp_perc: resolved.query_cov_hsp_perc,
    };
    let setup = super::stage_d_pipeline::set_up_batch(
        &letters,
        &records,
        match database {
            Some(_) => SubjectLetters::Packed(&[]),
            None => subject_letters,
        },
        profile,
        compo_adjust_mode,
        resolved.sum_stats,
        scoring,
    );
    match setup {
        Err(error) => batch_failure(error).ok().into(),
        Ok(false) => BatchSetUp::Unknown,
        Ok(true) => match read {
            None | Some(Ok(_)) => BatchSetUp::Searched,
            Some(Err((error, _))) => database_read_failure(error).ok().into(),
        },
    }
}

/// What NCBI's engine error is of a failed first read of the database (`first_database_read`)
/// in a batch with a valid context: the preliminary search's text, or SeqDB's own error with
/// `-in_pssm` (`TblastnDatabase::pssm_search_read_failure`); `Err`: LOSAT's refusal.
fn database_read_failure(error: anyhow::Error) -> Result<anyhow::Error> {
    match error.downcast::<PreliminaryReadError>() {
        Ok(read_error) => Ok(crate::blastinput::app::engine_error(read_error.0)),
        Err(error) if error.is::<crate::cli::NativeError>() => Ok(error),
        Err(error) => Err(error),
    }
}

impl TblastnArgs {
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:288-301:
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1411-1458:
    // PrintOneResultSet dispatches the complete result to outfmt 0/6/7.
    /// `settings` are the NCBI application settings that LOSAT reproduces
    /// (`ncbi_environment::ApplicationSettings`): the readers of the subjects
    /// and of the queries (`fasta_reader`) use them.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_scope_src.cpp:67-72
    pub fn run(
        self,
        settings: crate::blastinput::ncbi_environment::ApplicationSettings,
    ) -> Result<()> {
        self.run_io(settings, RunIo::Files)
    }

    /// `run` over the inputs of web ABI v1 (`web_api::v1_tblastn`): the query and the subjects
    /// are the call's bytes, read as the regular files with the same bytes, and the report is
    /// returned instead of written to standard output. The `query` and `subject` of `self` are
    /// only the names the report shows for the inputs; `out` is `None`. The diagnostics go to
    /// standard error, as the command line's.
    #[cfg(any(target_arch = "wasm32", test))]
    pub(crate) fn run_in_memory(
        self,
        settings: crate::blastinput::ncbi_environment::ApplicationSettings,
        query: &[u8],
        subject: &[u8],
    ) -> Result<Vec<u8>> {
        let report = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        let ran = self.run_io(
            settings,
            RunIo::Memory {
                query,
                subject,
                report: report.clone(),
            },
        );
        // `CATCH_ALL`'s silent end returns 0 with what the run wrote (the error archive of the
        // saved messages), as the library's call does (`api::run::run_local`).
        if let Err(error) = ran {
            if !crate::cli::is_silent_end(&error) {
                return Err(error);
            }
        }
        let bytes = std::mem::take(&mut *report.lock().expect("report mutex poisoned"));
        Ok(bytes)
    }

    fn run_io(
        self,
        settings: crate::blastinput::ncbi_environment::ApplicationSettings,
        io: RunIo<'_>,
    ) -> Result<()> {
        use crate::blastinput::app;
        use crate::blastinput::fasta_reader::{FastaInputSource, ReaderConfig};
        use crate::blastinput::input_files;
        use crate::blastinput::seq_range;
        let mut args = self;
        args.id_scope = args.parse_deflines.then(Default::default);
        // The run's `SERIAL_VERIFY_DATA_WRITE`, which `-in_pssm`'s reader reads
        // (`TblastnArgs::write_verification`).
        args.write_verification = settings.write_verification.clone();
        // NCBI's score block looks for the file of an ungapped search's `-matrix` that its tables
        // do not have when it sets up a batch (`TblastnArgs::find_matrix_file`); the search
        // refuses it there (`args::check_search_limits`).
        args.find_matrix_file(&settings.matrix_search)?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3624-3627
        // `ArchiveFormatRequested` parses `-outfmt` (blast_args.cpp:2745-2748) before the
        // option handlers run; a format that NCBI runs and LOSAT does not write is rejected
        // there too. `GetExportSearchStrategyStream` comes first: it creates the search strategy
        // file (`strategy::open_strategy_file`), which an error after it leaves empty; with it a
        // value that does not parse is left to the formatting handler
        // (`app::early_format_choice`).
        let strategy_file = crate::report::archive::strategy::open_strategy_file(
            args.export_search_strategy.as_deref(),
            args.out.as_deref(),
            "TBLASTN",
        )?;
        let choice = app::early_format_choice(&args.outfmt, strategy_file.is_some())?;
        let format = choice.report_format(app::Program::Tblastn, args.out.as_deref())?;
        // The BLAST archive (outfmt 11): its request, before anything is read or written.
        // The error archives follow the `-outfmt` text (`writes_error_archive`).
        // The run's options handle, which its archive and its search strategy share.
        let options = crate::report::archive::OptionsHandle::default();
        args.archive = super::archive::handle_for(
            &args,
            format == Some(app::ReportFormat::ArchiveFormat),
            &options,
        )
        .map(|handle| handle.with_outfmt_text(&args.outfmt));
        // The search strategy's request, with the archive's options handle.
        args.strategy = strategy_file.map(|file| {
            crate::report::archive::strategy::StrategyExport::for_args(
                &args,
                file,
                choice.number(),
                &options,
            )
        });
        let archive = args.archive.clone();
        // LOSAT's thread capability is checked before any input or output, as for BLASTN.
        crate::utils::threading::validate_threads(args.num_threads)?;
        // The diagnostics stream: with the archive it also saves what it writes as the
        // messages of the next archive.
        let mut stderr = crate::cli::CliStderr::process();
        let mut archive_tee;
        let diagnostics: &mut (dyn Write + Send) = match &archive {
            Some(handle) => {
                archive_tee = crate::report::archive::DiagTee::new(&mut stderr, handle.clone());
                &mut archive_tee
            }
            None => &mut stderr,
        };
        // The taxid or seqid list (line DI-8): NCBI's handler of the database arguments
        // (`s_GetTaxIDList`, `CSeqDBFileGiList`) comes before the other option handlers
        // (`m_BlastDbArgs` is among the first of `CBlastpAppArgs`' `m_Args`, blastp_args.cpp:61-65)
        // and before `CStdCmdLineArgs`, which opens the query and the output: its warning before
        // theirs (probe F04 of line DI-8), and a list error leaves no `-out`, only the error archive
        // that NCBI's `PrintErrorArchive` writes for outfmt 11 (probe E of line DI-8).
        let list = match args.database_name() {
            Some(_) => match args
                .db_filter
                .resolve(&settings, "tblastn", &mut *diagnostics)
            {
                Ok(list) => list,
                Err(error) => {
                    write_run_error_archive(&args, &io, archive.as_ref(), &error);
                    return Err(error);
                }
            },
            None => None,
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2553-2562
        // The handler of the database arguments reads the subjects (an empty subject set
        // fails there), before the query and the output are opened.
        // `-db`: the database is opened with the options (`search_cli`), as NCBI's
        // `InitializeSubject` opens it after `SetOptions`.
        let subjects_read = (|| -> Result<(
            Vec<FastaRecord>,
            seq_range::Placements,
            Option<Vec<crate::sequence::packed_subject::PackedNucleotide>>,
        )> {
            Ok(if args.db.is_some() {
            (Vec::new(), seq_range::Placements::default(), None)
        } else {
            let Some(subject_path) = args.subject.clone() else {
                return Err(app::missing_subject_error());
            };
            // Web ABI v1's subjects are the call's bytes: no file to name or open.
            let subject_input = match &io {
                RunIo::Files => {
                    input_files::check_utf8_file_name(&subject_path, "subject", "TBLASTN")?;
                    Input::File(input_files::open_input(&subject_path, "subject", "TBLASTN")?)
                }
                RunIo::Memory { subject, .. } => Input::Bytes(subject),
            };
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2537-2545
            // The subject range is read after the subject file is opened, before it is read.
            let subject_range = seq_range::parse_optional_range(
                args.subject_loc.as_deref(),
                seq_range::RangeRole::Subject,
            )?;
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:230-246
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2466-2470
            // TBLASTN's subjects are nucleotide (`read_proteins` = `IsProtein()` is false, so the
            // reader assumes nucleotides, `fAssumeNuc`). NCBI reads the records one at a time, writes
            // each one's messages as it reads it, and checks each one's range after reading it: a
            // range that starts past the end of a record stops the reading there
            // (`read_subjects`). Records without residues (and intervals without letters) are
            // read without a message and reported when the formatter sets up the subjects
            // (`before_prolog`).
            let config = ReaderConfig::subject("TBLASTN", false, settings.data_loaders)
                .with_defline(archive.is_some() || args.strategy.is_some())
                .with_parse_deflines(args.parse_deflines);
            // EXPERIMENT (LOSAT_X_TBNSTORE): each record's letters packed as it is read
            // (`x_read_subject_stores`), with a subject range after the record is cut to its
            // interval (the record `cut_subjects` makes).
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:515-520
            // Dispatch point: NCBI keeps the subjects in its scope and makes each one's search forms
            // from it. With the switch the records keep their deflines and the letters are held packed,
            // made again when a subject is searched or reported; without it the records keep them.
            // Diagnostic (LOSAT_TIMING): the subject FASTA for the common line (utils/xtiming.rs);
            // with LOSAT_X_TBNSTORE the packing of each record is in this interval too. No NCBI
            // counterpart.
            let t_read_fasta = crate::utils::xtiming::clock();
            // With the switch the records come back cut to the range (`Some` placements).
            let (read_subject_records, x_subject_stores, x_placements) = {
                let mut warn = |message: &[u8]| std::io::Write::write_all(&mut *diagnostics, message);
                // Both sides of the stores are recorded in x_subject_store.rs. The archive of
                // `-subject` writes every subject's whole letters (and its deflines) in every
                // batch, which the packed stores do not keep (they hold the interval of a
                // `-subject_loc` and no letters), so a run that writes it (or the search
                // strategy's requests, which hold the same Bioseqs) reads the records as the
                // switch's `=0` side does.
                let x_packed = super::x_subject_store::x_tbn_store()
                    && archive.is_none()
                    && args.strategy.is_none();
                match subject_input {
                    Input::File(file) => read_subject_input(
                        FastaInputSource::from_file(file, config),
                        &args,
                        subject_range.as_ref(),
                        x_packed,
                        &mut warn,
                    )?,
                    Input::Bytes(bytes) => read_subject_input(
                        FastaInputSource::from_bytes(bytes, config),
                        &args,
                        subject_range.as_ref(),
                        x_packed,
                        &mut warn,
                    )?,
                }
            };
            crate::utils::xtiming::add(crate::utils::xtiming::Phase::ReadFasta, t_read_fasta);
            let (subjects, subject_placements) = match x_placements {
                Some(placements) => (read_subject_records, placements),
                None => {
                    match seq_range::cut_subjects(&read_subject_records, subject_range.as_ref())? {
                        Some((cut, placements)) => (cut, placements),
                        None => (read_subject_records, seq_range::Placements::default()),
                    }
                }
            };
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
            // `Reset` of `m_Subjects` with a new `blast::CObjMgr_QueryFactory` of `subjects` (blast_args.cpp:2557)
            // fails for a subject input without records.
            if subjects.is_empty() {
                return Err(app::empty_subjects_error());
            }
            (subjects, subject_placements, x_subject_stores)
            })
        })();
        let (subjects, subject_placements, x_subject_stores) = match subjects_read {
            Ok(read) => read,
            Err(error) => {
                // NCBI reads the `-subject` file inside `Run`, after the diagnostics handler
                // that saves the messages is set: an error archive follows the error.
                write_run_error_archive(&args, &io, archive.as_ref(), &error);
                return Err(error);
            }
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3456-3481
        // The query is opened and the output file created before the options are checked;
        // `-` is standard input and standard output. Web ABI v1's query is the call's bytes.
        let query_input = match &io {
            // `-in_pssm`: the query is the PSSM, which the options' handler reads
            // (`check_options`); NCBI opens the default `-query -` and never reads it.
            RunIo::Files if args.in_pssm.is_some() => Input::Bytes(&[]),
            RunIo::Files => {
                input_files::check_utf8_file_name(&args.query, "query", "TBLASTN")?;
                // NCBI opens the query file inside `Run` too: its error is saved and an error
                // archive follows.
                match input_files::open_input(&args.query, "query", "TBLASTN") {
                    Ok(file) => Input::File(file),
                    Err(error) => {
                        write_run_error_archive(&args, &io, archive.as_ref(), &error);
                        return Err(error);
                    }
                }
            }
            RunIo::Memory { query, .. } => Input::Bytes(query),
        };
        // The organism report's `-out` file is written as NCBI's buffer writes it.
        let tax_report = format == Some(app::ReportFormat::TaxFormat);
        let mut stream = report_stream(&args, &io, tax_report)?;
        // Outfmt 0 to 4: the pairwise report (`report::pairwise::PairwiseView`).
        let pairwise = format
            .and_then(crate::report::pairwise::PairwiseView::of)
            .is_some();
        // The format as the reports read it: the number and the custom fields, or the value as
        // given when it has a delimiter, which `normalized` drops and the custom tabular report
        // reads (`tabular::tabular_spec`, which parses it again as NCBI does).
        // A value not read yet (`EarlyFormat::Deferred`) fails at the options' handlers below.
        let outfmt = match choice.parsed() {
            Some(parsed) if parsed.delimiter.is_empty() => parsed.normalized(),
            _ => args.outfmt.clone(),
        };
        let strategy = args.strategy.clone();
        let result = {
            let mut outputs =
                ReportOutputs::single(&outfmt, OutputSink::Writer(&mut stream), &mut *diagnostics);
            search_cli(
                args,
                list,
                query_input,
                &subjects,
                &subject_placements,
                x_subject_stores.as_deref(),
                &mut outputs,
                choice,
                settings,
            )
        };
        // A post of the search strategy left unwritten is LOSAT's internal error (`end_run`).
        let result = crate::report::archive::strategy::end_run(strategy.as_ref(), result);
        // NCBI writes an error archive after the archives of the batches before when the
        // application stops with one of its own errors, and when the run ends with messages that
        // no archive carried (`PrintErrorArchive` after its `CATCH_ALL`, whatever the status).
        if let Some(handle) = &archive {
            match &result {
                Err(error) => {
                    // The `-in_pssm` file, which NCBI opens inside `Run` too.
                    let in_run = crate::algorithm::pssm_search::in_pssm::file_error_in_run(error);
                    if let Some(lines) = super::archive::error_lines(error, in_run) {
                        let _ = write_error_archive(&io, handle, error, &lines, &mut stream);
                    } else if let Some(lines) =
                        crate::report::archive::silent_end_lines(handle, error)
                    {
                        // An exception that `CATCH_ALL` ends without a post: the archive of the
                        // messages saved before it, with no text before it.
                        let _ = crate::report::archive::write_error_archive(
                            handle,
                            &lines,
                            &mut stream,
                        );
                    }
                }
                Ok(()) => {
                    if handle.has_messages() {
                        let _ =
                            crate::report::archive::write_error_archive(handle, &[], &mut stream);
                    }
                }
            }
        }
        // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:252-255
        // The outfmt 0 formatter's stream throws when a write fails; for outfmt 6/7 NCBI
        // aborts and LOSAT reports the error (`cli::ReportStream::finish`).
        stream.finish(result, pairwise || tax_report)
    }
}

/// Where `run` reads its inputs and writes its report.
enum RunIo<'a> {
    /// The command line: the files of `-query` and `-subject` (`-`: standard input), and the
    /// `-out` file or standard output.
    Files,
    /// Web ABI v1 (`run_in_memory`): the call's bytes, and the report kept in memory.
    #[cfg_attr(
        not(any(target_arch = "wasm32", test)),
        expect(
            dead_code,
            reason = "made only by `run_in_memory`, for the web ABI (wasm32) and tests"
        )
    )]
    Memory {
        query: &'a [u8],
        subject: &'a [u8],
        report: std::sync::Arc<std::sync::Mutex<Vec<u8>>>,
    },
}

/// An input of `run`: the opened file of `-query` or `-subject`, or web ABI v1's bytes.
enum Input<'a> {
    File(crate::blastinput::input_files::InputFile),
    Bytes(&'a [u8]),
}

/// The reader of `search_cli`'s queries over an `Input`.
enum QuerySource<'a> {
    File(
        crate::blastinput::fasta_reader::FastaInputSource<
            crate::blastinput::input_files::InputFile,
        >,
    ),
    Bytes(crate::blastinput::fasta_reader::FastaInputSource<&'a [u8]>),
}

/// The report kept in memory (`RunIo::Memory`) as the report's stream.
struct MemoryReport(std::sync::Arc<std::sync::Mutex<Vec<u8>>>);

impl Write for MemoryReport {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.0
            .lock()
            .expect("report mutex poisoned")
            .extend_from_slice(buf);
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// The subjects of `source` as `run` reads them (`read_subjects`), or with LOSAT_X_TBNSTORE's
/// packed stores (`packed`, see `run`), with the reader's messages to `warn`.
fn read_subject_input<R: std::io::Read>(
    mut source: crate::blastinput::fasta_reader::FastaInputSource<R>,
    args: &TblastnArgs,
    range: Option<&crate::blastinput::seq_range::SequenceRange>,
    packed: bool,
    warn: &mut dyn FnMut(&[u8]) -> std::io::Result<()>,
) -> Result<(
    Vec<FastaRecord>,
    Option<Vec<crate::sequence::packed_subject::PackedNucleotide>>,
    Option<crate::blastinput::seq_range::Placements>,
)> {
    source.set_id_scope(args.id_scope.clone());
    if packed {
        let (records, stores, placements) = super::x_subject_store::x_read_subject_stores(
            &mut source,
            range,
            warn,
            args.num_threads,
        )?;
        Ok((records, stores, Some(placements)))
    } else {
        Ok((
            crate::blastinput::fasta_reader::read_subjects(&mut source, range, warn)?,
            None,
            None,
        ))
    }
}

/// The error archive of a run that stopped in `Run` with `error` before its search wrote
/// anything, when NCBI writes one (`PrintErrorArchive`): into `-out`, or nothing when that
/// cannot be opened (NCBI's call swallows that).
fn write_run_error_archive(
    args: &TblastnArgs,
    io: &RunIo<'_>,
    archive: Option<&crate::report::archive::ArchiveHandle>,
    error: &anyhow::Error,
) {
    let Some(handle) = archive else {
        return;
    };
    // `-outfmt "11 qseqid"`: NCBI's `PrintErrorArchive` stops at the number before it opens `-out`.
    if !handle.error_archive_enabled() {
        return;
    }
    let Some(lines) = super::archive::error_lines(error, true) else {
        return;
    };
    if let Ok(mut report) = report_stream(args, io, false) {
        let _ = write_error_archive(io, handle, error, &lines, &mut report);
        let _ = std::io::Write::flush(&mut report);
    }
}

/// The error archive of a run that stopped with `error` (`lines`), into `report`: on the command
/// line after the error's text (`write_error_archive_after_post`); web ABI v1 (`RunIo::Memory`)
/// reports the error itself, so there only the archive.
fn write_error_archive(
    io: &RunIo<'_>,
    handle: &crate::report::archive::ArchiveHandle,
    error: &anyhow::Error,
    lines: &[crate::report::diag_log::Message],
    report: &mut crate::cli::ReportStream,
) -> Result<()> {
    match io {
        RunIo::Files => {
            crate::report::archive::write_error_archive_after_post(handle, error, lines, report)
        }
        RunIo::Memory { .. } => crate::report::archive::write_error_archive(handle, lines, report),
    }
}

/// The report stream of a run: the `-out` file (created) or standard output, or the report
/// in memory (`RunIo::Memory`).
fn report_stream(
    args: &TblastnArgs,
    io: &RunIo<'_>,
    tax_report: bool,
) -> Result<crate::cli::ReportStream> {
    Ok(match io {
        RunIo::Files => crate::cli::ReportStream::of_cli(
            open_out_file(args)?,
            tax_report,
            args.strategy.as_ref(),
        ),
        RunIo::Memory { report, .. } => crate::cli::ReportStream {
            inner: Box::new(MemoryReport(report.clone())),
            failed: false,
        },
    })
}

/// The `-out` file, created (`-`: standard output).
fn open_out_file(args: &TblastnArgs) -> Result<Option<std::io::BufWriter<std::fs::File>>> {
    args.out
        .as_deref()
        .filter(|path| path.as_os_str() != "-")
        .map(|path| crate::cli::create_out_file(path, Some("TBLASTN")).map(std::io::BufWriter::new))
        .transpose()
}

/// What `search_cli` and `run_local_database` do after the options are processed and before the
/// query: the opening of the database of `-db` (`list`: its taxid or seqid list, `None` for a
/// search without `-db`), and what the requested reports and requests ask of it, which a run whose
/// PSSM query never joins NCBI's scope does not make (`joins_scope`, `PsiSearch::joins_scope`:
/// NCBI's tblastn makes its formatter and its search strategy after `AddTopLevelSeqEntry`,
/// tblastn_app.cpp:226, 241-270, 307; NCBI BLAST+ 2.17.0 writes the object manager's stop with
/// `-dbsize` and a list that keeps no sequence, with a database list whose lengths overflow, with
/// a database named `n/a`: line RG D3-1 fix 4).
fn open_search_database(
    args: &mut TblastnArgs,
    list: Option<&crate::seqdb::DbList>,
    settings: &crate::blastinput::ncbi_environment::ApplicationSettings,
    joins_scope: bool,
    outputs: &mut ReportOutputs<'_>,
) -> Result<Option<TblastnDatabase>> {
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:189-193
    // The database of `-db` is opened after the options are processed, before
    // `Query is Empty!` (as BLASTN's).
    let database = match args.database_name() {
        Some(name) => {
            let database = super::database::TblastnDatabase::open_restricted(
                name,
                settings,
                list,
                args.db_mask.option().as_ref(),
                outputs.diagnostics,
            )?
            .showing_gis(args.show_gis)
            .after_calls(args.first_query_ordinal > 0);
            if !joins_scope {
                return Ok(Some(database));
            }
            crate::algorithm::common::database::note_archive_database_filter(
                args.archive.as_ref(),
                args.strategy.as_ref(),
                &args.db_filter,
                database.list_taxids(),
                &database.statistics(),
                true,
                args.searchsp.is_some_and(|space| space > 0),
                args.dbsize.is_some_and(|size| size > 0),
                "TBLASTN",
            )?;
            crate::algorithm::common::database::check_formats_report_sums(
                || database.pairwise_database(),
                outputs
                    .formats
                    .iter_mut()
                    .map(|format| (format.outfmt, &mut format.refusal)),
                "TBLASTN",
            )?;
            Some(database)
        }
        None => None,
    };
    if !joins_scope {
        return Ok(database);
    }
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:126-133
    // A tabular field list that looks every taxid up (the plural names of the taxonomy database,
    // the organism report) asks whether `taxdb.bti` is on the search path: the plural names'
    // formatters warn without it (`tabular`), as BLASTN. A `-db` search asks its database, a
    // `-subject` search the search path of NCBI's databases (`subject_taxonomy`).
    match database.as_ref() {
        Some(database) => {
            if outputs
                .formats
                .iter()
                .filter_map(|format| super::tabular::tabular_spec(format.outfmt))
                .any(|spec| spec.looks_up_every_taxid())
            {
                args.tabular_taxdb = Some(database.taxdb_found());
            }
        }
        None => subject_taxonomy(
            args,
            Some(&crate::algorithm::common::database::TaxdbSetting::of(
                settings,
            )),
            outputs,
        )?,
    }
    // LOSAT's refusal of a database name for the Seq-annot formats, before the search
    // (`report::serial::seqalign::check_database`).
    if let Some(database) = database.as_ref() {
        for format in &outputs.formats {
            let choice = crate::blastinput::app::parse_formatting_string(format.outfmt)?;
            if let Some(serial) = crate::blastinput::app::ReportFormat::from_number(choice.number)
                .and_then(crate::report::serial::SerialFormat::of)
            {
                crate::report::serial::seqalign::check_database(
                    serial,
                    database.name(),
                    crate::blastinput::app::Program::Tblastn,
                )?;
            }
        }
    }
    Ok(database)
}

/// What the reports of `outputs` of a `-subject` search ask of the taxonomy database, through
/// `setting` (`common::database::subject_taxonomy`): asked by the command line
/// (`open_search_database`) and by a search of records (`run_local_formats`) at the same point,
/// after the options are processed.
fn subject_taxonomy(
    args: &mut TblastnArgs,
    setting: Option<&crate::algorithm::common::database::TaxdbSetting>,
    outputs: &ReportOutputs<'_>,
) -> Result<()> {
    let specs: Vec<_> = outputs
        .formats
        .iter()
        .filter_map(|format| super::tabular::tabular_spec(format.outfmt))
        .collect();
    let taxonomy =
        crate::algorithm::common::database::subject_taxonomy(setting, "TBLASTN", &specs)?;
    args.subject_taxdb = taxonomy.taxdb;
    args.tabular_taxdb = taxonomy.found;
    Ok(())
}

/// The part of `run` after the files are opened: NCBI's processing of the options,
/// `Query is Empty!`, the query batch size and the subjects without letters, the reading of
/// the queries and the search.
#[allow(clippy::too_many_arguments)]
fn search_cli(
    mut args: TblastnArgs,
    list: Option<crate::seqdb::DbList>,
    query_input: Input<'_>,
    subjects: &[FastaRecord],
    subject_placements: &crate::blastinput::seq_range::Placements,
    // EXPERIMENT (LOSAT_X_TBNSTORE): the subjects' letters, packed (the records keep none).
    x_subject_stores: Option<&[crate::sequence::packed_subject::PackedNucleotide]>,
    outputs: &mut ReportOutputs<'_>,
    choice: crate::blastinput::app::EarlyFormat,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
) -> Result<()> {
    use crate::blastinput::fasta_reader::{
        read_queries, FastaInputSource, QueryEnd, QueryRecords, ReaderConfig,
    };
    let resolved = args.check_options(outputs.diagnostics)?;
    // The options' handlers have read a value that the early reading left to them (`Deferred`).
    let choice = &choice.at_formatting_handler(&args.outfmt)?;
    let joins_scope = resolved
        .pssm
        .as_deref()
        .is_none_or(crate::algorithm::pssm_search::PsiSearch::joins_scope);
    let database = open_search_database(&mut args, list.as_ref(), &settings, joins_scope, outputs)?;
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:194-238
    // `-in_pssm` (PSI-TBLASTN): the PSSM's Bioseq is the one query, read with the options; no
    // FASTA input, no `Query is Empty!` and no query batches (one search of the PSSM,
    // tblastn_app.cpp:305-327).
    if let Some(psi) = resolved.pssm.clone() {
        drop(query_input);
        return search_pssm_query(
            &args,
            resolved,
            &psi,
            &[choice.number],
            args.id_scope.as_deref(),
            subjects,
            subject_placements,
            x_subject_stores,
            database.as_ref(),
            outputs,
        );
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:200-221
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:51,106
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:328-330
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-860
    // The query of tblastn is a protein (`QueryIsProtein()`, `fAssumeProt`) and its data
    // loaders are those of the subjects (`SDataLoaderConfig`). The position is taken on the
    // opened file, before it is read (`stream_is_empty`). When the subjects were read from
    // standard input too, `cin` has reached its end (a failed stream), so NCBI gets no
    // position for the query either, and its reader is at its end: no query batch is read.
    let config = ReaderConfig::query("TBLASTN", true, settings.data_loaders)
        .with_defline(args.archive.is_some() || args.strategy.is_some())
        .with_parse_deflines(args.parse_deflines);
    let (mut query_source, query_empty) = match query_input {
        Input::File(query_file) => {
            let shared_standard_input = args.query.as_os_str() == "-"
                && args
                    .subject
                    .as_deref()
                    .is_some_and(|path| path.as_os_str() == "-");
            let mut source = FastaInputSource::from_file(query_file, config);
            source.set_id_scope(args.id_scope.clone());
            let empty = !shared_standard_input && source.stream_is_empty();
            (QuerySource::File(source), empty)
        }
        // Web ABI v1's bytes, read as the regular file with the same bytes: `stream_is_empty`
        // finds such a file empty when it holds white space only.
        Input::Bytes(bytes) => {
            let mut source = FastaInputSource::from_bytes(bytes, config);
            source.set_id_scope(args.id_scope.clone());
            let empty = crate::ncbi_str::is_blank(bytes);
            (QuerySource::Bytes(source), empty)
        }
    };
    if query_empty {
        outputs
            .diagnostics
            .write_all(b"Warning: [tblastn] Query is Empty!\n")?;
        // NCBI returns from `Run` here, before its `PrintErrorArchive`: the message is not an
        // archive's.
        if let Some(handle) = &args.archive {
            handle.take_messages();
        }
        return Ok(());
    }
    let batch_size = before_prolog(subjects, x_subject_stores, outputs.diagnostics)?;
    let mut resolved = resolved;
    resolve_db_mask(&args, &mut resolved, database.as_ref(), outputs.diagnostics)?;
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:267-270
    // The XInclude formats written to standard output fail after the formatter is made.
    crate::blastinput::app::xinclude_check_for(choice, args.out.as_deref())?;
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:275-277
    // The records are read here, and their messages are written batch by batch
    // (`read_query_batches`). What ends the reading (the end of the input, an `eEOF`, a
    // reader error or LOSAT's rejection of a Seq-id line) comes in the batch that NCBI reads
    // it in, after the reports of the batches before.
    // Diagnostic (LOSAT_TIMING): the query FASTA for the common line (utils/xtiming.rs). No NCBI
    // counterpart.
    let t_read_fasta = crate::utils::xtiming::clock();
    let QueryRecords { records, end } = match &mut query_source {
        QuerySource::File(source) => read_queries(source),
        QuerySource::Bytes(source) => read_queries(source),
    };
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::ReadFasta, t_read_fasta);
    drop(query_source);
    let reading = match end {
        QueryEnd::Input => QueryReading::End,
        QueryEnd::EofError => QueryReading::BlankLines,
        QueryEnd::Error { error, warnings } => QueryReading::Error { error, warnings },
    };
    search(
        &args,
        resolved,
        batch_size,
        &records,
        reading,
        subjects,
        subject_placements,
        x_subject_stores,
        database.as_ref(),
        outputs,
    )
}

/// tblastn `-in_pssm` (PSI-TBLASTN): the search of the PSSM's query (the file read with the
/// options, `check_options`), from where the entries would read their FASTA queries: the command
/// line's (`search_cli`) and those of the records (`run_local_formats`, `search_database`, the
/// library API's `Inputs::pssm`). No FASTA input, no `Query is Empty!` and no query batches (one
/// search of the PSSM, tblastn_app.cpp:305-327). A request whose PSSM file has members that the
/// reader skipped is refused here, once the database is opened (NCBI's error of a database it
/// cannot open comes first) and before any output, and so are LOSAT's limits of the PSSM search for
/// each requested format (`outfmts`, `in_pssm::check_limits`), and the query record's refusals;
/// then the messages of NCBI's object manager about the query's annotations and its exception
/// (`PsiSearch::write_annot_posts`, `annot_abort`), the formatter's set-up of the subjects (a
/// subject without letters) and the XInclude formats without a file, as for a FASTA query
/// (`before_prolog`, `ReportOutputs::xinclude_check`).
///
/// The query is the PSSM's Bioseq as the reports show it (`PsiSearch::query_record`); with
/// `-parse_deflines` its ids join `scope`, the scope of the `-subject` records (`IdScope::
/// add_record`; NCBI BLAST+ 2.17.0 fails on a subject with the query's id), and the database's
/// sequences are checked as for a FASTA query (`QueryScopes`). The library's observer receives it
/// (`FormatObserver::pssm_query`).
#[allow(clippy::too_many_arguments)]
fn search_pssm_query(
    args: &TblastnArgs,
    resolved: ResolvedTblastnArgs,
    psi: &crate::algorithm::pssm_search::PsiSearch,
    outfmts: &[i32],
    scope: Option<&crate::seqid::id_scope::IdScope>,
    subjects: &[FastaRecord],
    subject_placements: &crate::blastinput::seq_range::Placements,
    // EXPERIMENT (LOSAT_X_TBNSTORE): the subjects' letters, packed (the records keep none).
    x_subject_stores: Option<&[crate::sequence::packed_subject::PackedNucleotide]>,
    database: Option<&TblastnDatabase>,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // An XInclude format without a file name ends the run with NCBI's error before the search
    // (`ReportOutputs::xinclude_check` below): LOSAT's refusals of the search's reports wait for it (NCBI
    // 2.17.0, line RG F/d36 probes tbnp.o13, tbnp.o14: `Please provide a file name for outfmt
    // 14.`).
    if outputs.xinclude_check(args.out.as_deref()).is_ok() {
        for &outfmt in outfmts {
            super::in_pssm::check_limits(args, &resolved, psi, outfmt)?;
        }
    }
    // The query record's refusal comes before any output, in a run that searches.
    let searches = psi.searches(args.strategy.is_some(), database.is_some());
    let records = vec![psi.query_record_in(args.parse_deflines, scope, "TBLASTN", searches)?];
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:222-226, 241-270 (the query's
    // Bioseq joins the scope, cited at `PsiSearch::write_annot_posts`, before the formatter is
    // made and its prolog written)
    // The messages of NCBI's object manager about the query's annotations, after the database
    // and the subjects, before the formatter's warnings and anything on the standard output;
    // its exception ends the run there (`PsiSearch::write_annot_posts`, `annot_abort`).
    psi.write_annot_posts(outputs.diagnostics, args.archive.as_ref())?;
    if let Some(abort) = psi.annot_abort() {
        return Err(abort);
    }
    crate::blastinput::input_files::write_empty_subject_warnings_by_length(
        subjects,
        |index| x_subject_stores.map_or(subjects[index].seq().len(), |stores| stores[index].len()),
        "tblastn",
        outputs.diagnostics,
    )?;
    let mut resolved = resolved;
    resolve_db_mask(args, &mut resolved, database, outputs.diagnostics)?;
    outputs.xinclude_check(args.out.as_deref())?;
    if let Some(observer) = outputs.observer.as_deref_mut() {
        observer.pssm_query(&records[0]);
    }
    search(
        args,
        resolved,
        1,
        &records,
        QueryReading::End,
        subjects,
        subject_placements,
        x_subject_stores,
        database,
        outputs,
    )
}

/// NCBI's steps of tblastn between `Query is Empty!` and the outfmt 0 prolog that concern
/// the subjects: the query batch size, then the formatter, which sets up the subjects and
/// warns about each one without letters (a record without residues, or an interval that
/// starts just past its record's end). Such a subject stays in the database statistics
/// with no letters and is never searched (`x_scan_plan` skips subjects
/// shorter than a word). Returns the query batch size.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:217-221
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:129-136
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:773-789
/// NCBI's gapped tblastn with the Gumbel block never stops for the average subject length,
/// also when every subject is without letters: its shortest subject is raised to
/// `BLAST_SEQSRC_MINLENGTH` before the division by 3 (`run_local_search_with_pool` ports
/// blast_setup.c:969-980 and seqsrc_multiseq.cpp:226-240). Without it (`-searchsp`, `OLD_FSC`
/// or an ungapped search) the average is used, and below one letter the batch fails
/// (`stage_d_pipeline::SetupEngineError::SUBJECT_LENGTH_INVALID`).
fn before_prolog(
    subjects: &[FastaRecord],
    x_subject_stores: Option<&[crate::sequence::packed_subject::PackedNucleotide]>,
    diagnostics: &mut dyn Write,
) -> Result<u32> {
    let batch_size = crate::blastinput::app::query_batch_size(20000)?;
    // EXPERIMENT (LOSAT_X_TBNSTORE): a packed subject's length is its store's.
    crate::blastinput::input_files::write_empty_subject_warnings_by_length(
        subjects,
        |index| x_subject_stores.map_or(subjects[index].seq().len(), |stores| stores[index].len()),
        "tblastn",
        diagnostics,
    )?;
    Ok(batch_size)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:140-153 (cited at
// `blastinput::db_mask::DbMaskOption::resolve`)
/// The formatter's resolution of the database masking (`-db_soft_mask`, `-db_hard_mask`; series
/// NP, M4b PR 2) into `resolved.db_mask`, where NCBI makes its formatter: after `Query is Empty!`
/// and the formatter's set-up of the subjects, before the XInclude checks. Its warning (one
/// message of the archive) or its database error, and the program options of the archive and of
/// the search strategy. Nothing without a database or without the options.
fn resolve_db_mask(
    args: &TblastnArgs,
    resolved: &mut ResolvedTblastnArgs,
    database: Option<&TblastnDatabase>,
    diagnostics: &mut dyn Write,
) -> Result<()> {
    if let (Some(database), Some(option)) = (database, args.db_mask.option()) {
        resolved.db_mask = database.resolve_mask_option(
            &option,
            args.archive.as_ref(),
            args.strategy.as_ref(),
            diagnostics,
        )?;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:288-301:
// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:429-467:
// NCBI formats one result set without searching again; several requested formats
// are several CBlastFormat printers over the same result set.
/// Runs one TBLASTN search over already parsed records and writes every requested
/// output format from the same result (the entry of web ABI v2; the command line reaches
/// `search` through `search_cli`).
///
/// Each requested `-outfmt` is read as on the command line, and NCBI's checks of the
/// options and of the subjects come before the search. The `-query` and `-subject` values
/// of `args` are used only as display names.
///
/// The records are those of NCBI's reader (`fasta_reader`): the subjects as NCBI's
/// nucleotide reader reads them, the queries as its protein reader reads them. Queries
/// without a record are NCBI's `Query is Empty!`.
pub fn run_local(
    mut args: TblastnArgs,
    query_records: &[FastaRecord],
    query_end: &QueryEnd,
    subject_records: &[FastaRecord],
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // NCBI's options factory checks the task first (`TaskParser::check`); the command line's
    // parser has checked it, a caller of `run_local` has not.
    super::args::TBLASTN_TASK.check(&mut args.task)?;
    // `run_local` searches subject records; a `-db` search runs through `run_local_database`.
    if args.db.is_some() {
        anyhow::bail!(crate::api::local_blast::RUN_LOCAL_DB_REFUSAL);
    }
    report_formats(outputs, &args)?;
    // The BLAST archive (outfmt 11) among the formats: its request and a diagnostics stream that
    // saves the messages. The records' deflines are
    // those the caller's reader kept; records without one get `>` and their title.
    let archive_requested = outputs
        .formats
        .iter()
        .any(|format| super::stage_e_report::format_number(format.outfmt) == 11);
    // The archive shares the options handle of the run's search strategy (LOSAT Web's
    // `-export_search_strategy`, `web_adapter::give_strategy`), as on the command line (`run_io`).
    if archive_requested {
        args.archive = super::archive::handle_for(
            &args,
            true,
            &crate::report::archive::strategy::options_of(args.strategy.as_ref()),
        )
        .map(|handle| outputs.with_archive_outfmt(handle));
    }
    let archive = args.archive.clone();
    let strategy = args.strategy.clone();
    let no_query = query_records.is_empty() && query_end.is_input() && args.in_pssm.is_none();
    outputs.run_records(strategy.as_ref(), archive.as_ref(), no_query, |outputs| {
        run_local_formats(args, query_records, query_end, subject_records, outputs)
    })
}

/// Each requested format as NCBI parses `-outfmt` before its option handlers run, and LOSAT's
/// rejection of a format that its TBLASTN does not write (`ReportOutputs::early_choices`): the
/// command line's checks before it opens `-out` (`run_io`). The formatting handler's come after
/// (`read_formats`), and so does the error of a value left to it with a search strategy file.
fn report_formats(outputs: &ReportOutputs<'_>, args: &TblastnArgs) -> Result<()> {
    outputs
        .early_choices(
            crate::blastinput::app::Program::Tblastn,
            args.out.as_deref(),
            args.export_search_strategy.is_some(),
        )
        .map(drop)
}

/// The numbers of the requested formats, each read as on the command line. `out`: the file of
/// `api::run_to_file` (a format with the caller's receiver of side files takes that receiver's,
/// `FormatOutput::side_out`). `strategy`: the argv names a search strategy file, with which a value
/// that does not parse is left to the formatting handler of the options
/// (`TblastnArgs::check_options_for_formats`, which reads the arguments' `-outfmt`: LOSAT Web gives
/// such a value as the run's one format and as the arguments' `-outfmt`); it has no number here.
fn read_formats(
    outputs: &ReportOutputs<'_>,
    out: Option<&std::path::Path>,
    strategy: bool,
) -> Result<Vec<i32>> {
    use crate::blastinput::app;
    let mut outfmts = Vec::with_capacity(outputs.formats.len());
    for format in &outputs.formats {
        let Some(choice) = app::early_format_choice(format.outfmt, strategy)?
            .parsed()
            .cloned()
        else {
            continue;
        };
        // A format without a report: NCBI's formatting handler refuses SAM, AIRR and FASTA here;
        // the XInclude formats without a file come after `Query is Empty!`
        // (`ReportOutputs::xinclude_check`).
        if app::report_format(&choice, app::Program::Tblastn, format.side_out(out))?.is_none() {
            app::formatting_handler_check(&choice, app::Program::Tblastn)?;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2887
        // The formatting handler checks the custom delimiter of each format (the command line's
        // in `check_options`).
        app::formatting_handler_check(&choice, app::Program::Tblastn)?;
        outfmts.push(choice.number);
    }
    Ok(outfmts)
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:189-193
/// Runs one TBLASTN search of the records `query_records` against the database of `args` and
/// writes every requested output format from the same result (the library API's entry; the
/// command line's is `run`). `settings` are the NCBI application settings, which find the
/// database. The steps are `run`'s and `search_cli`'s, in their order; what differs is what they
/// read: the records are read already, and the outputs are the caller's. An archive (outfmt 11)
/// ends as the command line's does (`ReportOutputs::finish_archive`).
pub fn run_local_database(
    mut args: TblastnArgs,
    query_records: &[FastaRecord],
    query_end: &QueryEnd,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // NCBI's options factory checks the task first (`TaskParser::check`); the command line's
    // parser has checked it, a caller of `run_local` has not.
    super::args::TBLASTN_TASK.check(&mut args.task)?;
    // The file of an unknown `-matrix`, as the command line's `run_io` looks it up.
    args.find_matrix_file(&settings.matrix_search)?;
    args.write_verification = settings.write_verification.clone();
    report_formats(outputs, &args)?;
    let archive_requested = outputs
        .formats
        .iter()
        .any(|format| super::stage_e_report::format_number(format.outfmt) == 11);
    // The archive shares the options handle of the run's search strategy (LOSAT Web's
    // `-export_search_strategy`, `web_adapter::give_strategy`), as on the command line (`run_io`).
    args.archive = super::archive::handle_for(
        &args,
        archive_requested,
        &crate::report::archive::strategy::options_of(args.strategy.as_ref()),
    )
    .map(|handle| outputs.with_archive_outfmt(handle));
    crate::utils::threading::validate_threads(args.num_threads)?;
    let archive = args.archive.clone();
    let strategy = args.strategy.clone();
    let no_query = query_records.is_empty() && query_end.is_input() && args.in_pssm.is_none();
    outputs.run_records(strategy.as_ref(), archive.as_ref(), no_query, |outputs| {
        search_database(args, query_records, query_end, settings, outputs)
    })
}

/// `run_local_database` from the taxid or seqid list on.
fn search_database(
    mut args: TblastnArgs,
    query_records: &[FastaRecord],
    query_end: &QueryEnd,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // The handler of the database arguments reads the list before the other handlers (`run`).
    let list = args
        .db_filter
        .resolve(&settings, "tblastn", &mut *outputs.diagnostics)?;
    // `-out` (`run_to_file`) is made where the command line opens it (`run_io`); the formatting
    // handler's checks come after it.
    outputs.open_out()?;
    let outfmts = read_formats(
        outputs,
        args.out.as_deref(),
        args.export_search_strategy.is_some(),
    )?;
    let resolved = args.check_options_for_formats(&outfmts, outputs.diagnostics)?;
    let joins_scope = resolved
        .pssm
        .as_deref()
        .is_none_or(crate::algorithm::pssm_search::PsiSearch::joins_scope);
    let database = open_search_database(&mut args, list.as_ref(), &settings, joins_scope, outputs)?
        .ok_or_else(crate::blastinput::app::missing_subject_error)?;
    // `-in_pssm` (the library API's `Inputs::pssm_database`): the query is the PSSM's. No subject
    // record holds an id (the database's sequences are checked against the query's,
    // `QueryScopes`).
    if let Some(psi) = resolved.pssm.clone() {
        return search_pssm_query(
            &args,
            resolved,
            &psi,
            &outfmts,
            None,
            &[],
            &crate::blastinput::seq_range::Placements::default(),
            None,
            Some(&database),
            outputs,
        );
    }
    // Records of an input read to its end; an end that reads past them (blank lines, a reader
    // error) is met in the first batch (`QueryEnd`).
    if query_records.is_empty() && query_end.is_input() {
        outputs
            .diagnostics
            .write_all(b"Warning: [tblastn] Query is Empty!\n")?;
        return Ok(());
    }
    let batch_size = before_prolog(&[], None, outputs.diagnostics)?;
    let mut resolved = resolved;
    resolve_db_mask(&args, &mut resolved, Some(&database), outputs.diagnostics)?;
    outputs.xinclude_check(args.out.as_deref())?;
    search(
        &args,
        resolved,
        batch_size,
        query_records,
        query_end.reading(),
        &[],
        &crate::blastinput::seq_range::Placements::default(),
        None,
        Some(&database),
        outputs,
    )
}

fn run_local_formats(
    mut args: TblastnArgs,
    query_records: &[FastaRecord],
    query_end: &QueryEnd,
    subject_records: &[FastaRecord],
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    use crate::blastinput::app;
    use crate::blastinput::seq_range;
    // The thread count is checked before the inputs are read, as on the command line (`run_io`).
    crate::utils::threading::validate_threads(args.num_threads)?;
    // The subject range and its record checks come where NCBI reads the subjects (`run`):
    // the messages of each record read, up to a record whose range starts past its end.
    let subject_range = seq_range::parse_optional_range(
        args.subject_loc.as_deref(),
        seq_range::RangeRole::Subject,
    )?;
    for record in
        &subject_records[..seq_range::subjects_read(subject_records, subject_range.as_ref())]
    {
        outputs.diagnostics.write_all(&record.warnings)?;
    }
    let ranged_subjects = seq_range::cut_subjects(subject_records, subject_range.as_ref())?;
    if subject_records.is_empty() {
        return Err(app::empty_subjects_error());
    }
    // `-out` (`run_to_file`) is made where the command line opens it (`run_io`); the formatting
    // handler's checks come after it.
    outputs.open_out()?;
    let outfmts = read_formats(
        outputs,
        args.out.as_deref(),
        args.export_search_strategy.is_some(),
    )?;
    let resolved = args.check_options_for_formats(&outfmts, outputs.diagnostics)?;
    // The taxonomy database of the reports, where the command line asks (`open_search_database`).
    let setting = args.taxdb_setting.take();
    subject_taxonomy(&mut args, setting.as_ref(), outputs)?;
    let whole_subjects = seq_range::Placements::default();
    let (searched_subjects, subject_placements) = match &ranged_subjects {
        Some((cut, placements)) => (cut.as_slice(), placements),
        None => (subject_records, &whole_subjects),
    };
    // `-in_pssm` (the library API's `Inputs::pssm`): the query is the PSSM's, whose ids join the
    // scope of the subject records read (`pssm_search::scope_of_subjects`).
    if let Some(psi) = resolved.pssm.clone() {
        let id_scope = crate::algorithm::pssm_search::scope_of_subjects(
            args.parse_deflines,
            subject_records,
            "TBLASTN",
        )?;
        return search_pssm_query(
            &args,
            resolved,
            &psi,
            &outfmts,
            id_scope.as_ref(),
            searched_subjects,
            subject_placements,
            None,
            None,
            outputs,
        );
    }
    // Records of an input read to its end; an end that reads past them (blank lines, a reader
    // error) is met in the first batch (`QueryEnd`).
    if query_records.is_empty() && query_end.is_input() {
        outputs
            .diagnostics
            .write_all(b"Warning: [tblastn] Query is Empty!\n")?;
        return Ok(());
    }
    // The caller owns the records: their letters are read in place.
    let batch_size = before_prolog(searched_subjects, None, outputs.diagnostics)?;
    outputs.xinclude_check(args.out.as_deref())?;
    search(
        &args,
        resolved,
        batch_size,
        query_records,
        query_end.reading(),
        searched_subjects,
        subject_placements,
        None,
        None,
        outputs,
    )
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:270-280
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:118-119
/// The prologs of the formats (`formatter.PrintProlog()`). `before_search`: the outfmt 0 to 4
/// prolog of every sink (a library `File` sink is kept for the run,
/// `api::local_blast::KeptFile`), where NCBI writes it: before it reads the first query batch and
/// writes its search strategy request, so that the prolog's error (a malformed `./.ncbirc` under
/// `-html`) and a failed write end the run there, before LOSAT's refusal of a strategy file on
/// the report stream; flushed as every program's (`pairwise::flushes_after_prolog`). Otherwise,
/// for a run that fails in its first query batch, the prologs that come with the first report,
/// written and flushed: the start of the XML2 and JSON documents, and outfmt 20's line of the
/// field names with its formatter's warning. `subject_path` names a FASTA subject set of
/// `num_subjects` records (`subject_letters`); a database search shows its `Database:` lines
/// (`database`).
#[allow(clippy::too_many_arguments)]
pub(super) fn write_first_batch_prologs(
    outputs: &mut ReportOutputs<'_>,
    subject_path: Option<&std::path::Path>,
    num_subjects: usize,
    subject_letters: SubjectLetters<'_>,
    database: Option<&crate::report::pairwise::PairwiseDatabase>,
    taxdb_found: bool,
    // The program line of outfmt 0 (`TBLASTN`; `PSITBLASTN` with `-in_pssm`).
    program: &str,
    // The run's `-html` data (`ResolvedTblastnArgs::html`; `None`: the plain report).
    html: Option<&crate::report::html_ids::HtmlLinks>,
    before_search: bool,
) -> Result<()> {
    let name = subject_path.map_or_else(Vec::new, crate::report::pairwise::subject_set_name);
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_multiseq.cpp:175-181
    // The subject set's total as NCBI's `Int8`.
    let total: i64 = match database {
        Some(_) => 0,
        None => subject_letters.lengths().map(|length| length as i64).sum(),
    };
    for format in outputs.formats.iter_mut() {
        if (0..=4).contains(&crate::blastinput::app::parse_formatting_string(format.outfmt)?.number)
        {
            if !before_search {
                continue;
            }
            let flush = crate::report::pairwise::flushes_after_prolog(
                html.is_some(),
                !format.sink.is_writer(),
            );
            let mut writer = format.sink.open()?;
            crate::report::pairwise::write_tblastn_pairwise_prolog(
                &mut writer,
                program,
                crate::wording::REPORT_VERSION,
                &name,
                num_subjects as i64,
                total,
                database,
                html.map(|links| {
                    let posts = crate::report::query_warnings::html_posts_sink(
                        outputs.html_posts.as_deref_mut(),
                        &mut *outputs.diagnostics,
                    );
                    (links, posts as &mut dyn Write)
                }),
            )?;
            if flush {
                writer.flush()?;
            }
        } else if before_search {
            // The other formats' prologs come with the first report.
        } else if let Some(kind) = super::stage_e_report::xml2_kind(format.outfmt) {
            // The XML2 and JSON documents start with their prolog too.
            let mut writer = format.sink.open()?;
            writer.write_all(&crate::report::blast_xml2::framing::prolog_bytes(kind))?;
            writer.flush()?;
        } else if let Some(spec) = super::tabular::tabular_spec(format.outfmt) {
            // Outfmt 20's line of the field names and its formatter's warning (nothing for the
            // other tabular kinds).
            let mut writer = format.sink.open()?;
            super::tabular::write_tabular_prolog(
                &spec,
                taxdb_found,
                &mut writer,
                &mut *outputs.diagnostics,
            )?;
            writer.flush()?;
        }
    }
    Ok(())
}

/// The traceback's filters of a batch's final hit lists, one per query of `queries`.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1763-1777
// `-max_hsps`, `-qcov_hsp_perc` and `-subject_besthit` filter the batch's final hit
// lists (each query's context is its searched letters), which the search has sorted by
// e-value where NCBI sorts the results (a positive `BlastSeqSrcGetTotLen`); the sort is
// a total order (best e-value, first score, OID), so sorting again after the filter
// gives NCBI's order, and where NCBI does not sort the lists keep their order.
// `-culling_limit` runs NCBI's culling pipe first on each query's hit list
// (`BlastHSPStreamTBackClose`, before the filters of `s_FilterBlastResults`;
// hspfilter_culling.c:699-733); its lists leave in the order of the tree, and the
// sort by e-value here is NCBI's after the pipes, where it sorts. `-best_hit_overhang` /
// `-best_hit_score_edge` run NCBI's best hit pipe there instead (hspfilter_besthit.c:520-539),
// with the preliminary hit list size of the search (`compo_adjust_mode`: the batch's composition
// mode, as the preliminary stage's).
fn filter_final_hitlists(
    batch_results: &mut [super::stage_d_results::KappaResultHitList],
    queries: &[Vec<u8>],
    resolved: &ResolvedTblastnArgs,
    compo_adjust_mode: BlastCompoAdjustMode,
    sorts_results: bool,
) -> Result<()> {
    if resolved.max_hsps > 0
        || resolved.query_cov_hsp_perc > 0.0
        || resolved.culling_limit > 0
        || resolved.best_hit.is_some()
        || resolved.subject_besthit
    {
        let best_hit = resolved.best_hit.map(|params| {
            (
                params,
                super::stage_d_pipeline::preliminary_hitlist_size(
                    resolved.hitlist_size,
                    compo_adjust_mode,
                    resolved.gapped,
                ),
            )
        });
        for (hitlist, query) in batch_results.iter_mut().zip(queries) {
            if resolved.culling_limit > 0 {
                hitlist.cull(resolved.culling_limit, query.len())?;
            }
            if let Some((params, prelim_size)) = best_hit {
                hitlist.best_hit(params, prelim_size)?;
            }
            if resolved.max_hsps > 0 || resolved.query_cov_hsp_perc > 0.0 {
                hitlist.filter_results(resolved.max_hsps, resolved.query_cov_hsp_perc, query.len());
            }
            if resolved.subject_besthit {
                hitlist.subject_besthit();
            }
            if sorts_results {
                hitlist.sort_by_evalue_for_positive_totlen();
            }
        }
    }
    Ok(())
}

/// The search of `run_local` and of the CLI, after NCBI's checks of the options, the
/// reading of the subjects and the steps before the prolog (`before_prolog`): LOSAT's
/// limits and the environment, then the query batches (`read_query_batches`), the search of
/// the batches that NCBI searches and their reports.
#[allow(clippy::too_many_arguments)]
/// How the pairwise reports of `formats` read each query's subjects, as the first pairwise
/// format's writer reads them (`PairwiseView::reading` of `-sorthits`, `-sorthsps`,
/// `-num_descriptions` and `-num_alignments`; the formats of a run whose loads they make share
/// one view's order, `FormatThrow::and`), for the loads of a subject whose lookup throws and the
/// texts of its reads (line DBID item pfoff-3); `None` without one.
fn pairwise_reading(
    resolved: &ResolvedTblastnArgs,
    formats: &[crate::api::local_blast::FormatOutput<'_>],
) -> Option<crate::report::pairwise::PairwiseReading> {
    let view = formats.iter().find_map(|format| {
        format
            .report_format()
            .and_then(crate::report::pairwise::PairwiseView::of)
    })?;
    Some(view.reading(
        resolved.formatting.report_sort(),
        crate::report::db_subject::Shown {
            descriptions: resolved.num_descriptions,
            alignments: resolved.num_alignments,
        },
    ))
}

fn search(
    args: &TblastnArgs,
    resolved: ResolvedTblastnArgs,
    batch_size: u32,
    query_records: &[FastaRecord],
    reading: QueryReading,
    subject_records: &[FastaRecord],
    subject_placements: &crate::blastinput::seq_range::Placements,
    // EXPERIMENT (LOSAT_X_TBNSTORE): the subjects' letters, packed (the records keep none).
    x_subject_stores: Option<&[crate::sequence::packed_subject::PackedNucleotide]>,
    // `-db`: the database (the subject records are then none).
    database: Option<&TblastnDatabase>,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    use crate::blastinput::app;
    use crate::blastinput::seq_range;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:253-259
    // What the run writes before its query batches ends here.
    if let Some(observer) = outputs.observer.as_deref_mut() {
        observer.queries_begin();
    }
    // The records are NCBI's reader's records: the subjects of its nucleotide reader (a
    // subject without letters stays in the database, its warning written before the prolog,
    // `before_prolog`) and the queries of its protein reader. The query range applies to
    // every query record as NCBI's batch reader reads it (`seq_range::cut_queries`): a record
    // is searched cut, and one whose interval starts just past its end is a query without
    // data (`Sequence contains no data`, as a record without residues).
    let input_records = query_records;
    let (searched_queries, query_input) = match resolved.query_range.as_ref() {
        Some(range) => {
            let ranged = seq_range::cut_queries(query_records, range);
            (std::borrow::Cow::Owned(ranged.records), ranged.input)
        }
        None => (
            std::borrow::Cow::Borrowed(query_records),
            seq_range::QueryInput::whole(query_records),
        ),
    };
    let SearchSettings {
        compo_adjust_mode,
        identity_cbs_reset,
        scoring,
    } = search_settings(&resolved)?;
    // `-in_pssm` (PSI-TBLASTN): the query's PSSM, which the batch's search and reports read.
    let pssm = resolved.pssm.clone();
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_memento_priv.hpp:64-74
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:288
    // IDENTITY's reset (`search_settings`) writes through the memento's pointers into the options
    // that every query batch reuses: the first batch resets the mode and warns its queries, and
    // the later batches start in mode 0 without the warning. These are the searched queries of
    // the first batch.
    let identity_warned = if identity_cbs_reset {
        query_input
            .batches(batch_size)
            .first()
            .map_or(0, |batch| batch.searched.end)
    } else {
        0
    };
    // LOSAT's refusals of the masked search (`-db_hard_mask`; series NP, M4b PR 2) come where NCBI
    // starts the search, after its checks.
    if let (Some(database), Some(mask)) = (database, resolved.db_mask.as_ref()) {
        database.check_mask_limits(mask, usize::try_from(resolved.word_size)?)?;
    }
    app::check_search_environment("TBLASTN", database.is_none())?;
    let subjects = subject_records;
    // The search reads the residues of the subject records in place (NCBI's subjects stay in the
    // reader's records, `CBlastQueryVector`); they are not copied for it. EXPERIMENT
    // (LOSAT_X_TBNSTORE): or from their packed form.
    let subject_seqs: Vec<&[u8]> = match x_subject_stores {
        Some(_) => Vec::new(),
        None => subjects.iter().map(FastaRecord::seq).collect(),
    };
    let subject_letters = match x_subject_stores {
        Some(stores) => SubjectLetters::Packed(stores),
        None => SubjectLetters::Records(&subject_seqs),
    };
    // What the requested formats need beyond outfmt 0, 6 and 7, decided once for the search
    // (`ReportNeeds`): with the Seq-aligns (outfmt 8, 12, 11), the engine positives of the
    // composition-adjusted matches (`LocalStageDProfile::engine_positives`).
    // The archive's Seq-aligns are those of outfmt 8.
    let engine_positives =
        crate::report::record::ReportNeeds::of(outputs.formats.iter().filter_map(|format| {
            app::ReportFormat::from_number(super::stage_e_report::format_number(format.outfmt))
        }))
        .seqalign
            || args.archive.is_some();
    // A database's numbers and `-dbsize` (`SubjectSetNumbers`), and the `Database:` lines of
    // a database search.
    let subject_set = SubjectSetNumbers {
        database: database.map(TblastnDatabase::numbers),
        dbsize: args.dbsize,
        searchsp: args.searchsp,
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:922-923
        // `-searchsp` sets OLD_FSC (blast_args.cpp:303-307).
        old_fsc: args.searchsp.is_some() || crate::blastinput::app::old_fsc_environment(),
        read_fails: None,
    };
    // Without the Gumbel block NCBI's shortest subject length is the statistics' average
    // (blast_setup.c:975), which divides by their number of sequences: a count of 0 is LOSAT's
    // refusal in the setup of the first batch with a valid context
    // (`stage_d_pipeline::run_local_search_with_pool`).
    let pairwise_database =
        database.map(|database| database.pairwise_database_masked(resolved.db_mask.as_ref()));
    let subject_path = match database {
        Some(_) => None,
        None => Some(
            args.subject
                .as_ref()
                .context("TBLASTN subject file missing")?,
        ),
    };
    // The program line of outfmt 0 (`PSITBLASTN` with `-in_pssm`).
    let program = if pssm.is_some() {
        super::in_pssm::PROGRAM
    } else {
        "TBLASTN"
    };
    // The outfmt 0 to 4 prolog, before the first batch is read and its request written
    // (`write_first_batch_prologs`), as in every program.
    write_first_batch_prologs(
        outputs,
        subject_path.map(|path| path.as_path()),
        subjects.len(),
        subject_letters,
        pairwise_database.as_ref(),
        args.tabular_taxdb.unwrap_or(false),
        program,
        resolved.html.as_deref(),
        true,
    )?;
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:275-278
    // NCBI reads the query batches one at a time after the outfmt 0 prolog. The run stops at
    // the first batch that NCBI fails (a reader error, a batch without a query, a batch of
    // queries without letters), after the reports of the batches before (without the
    // epilog) and the messages of that batch's records (`read_query_batches`).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_impl.cpp:59-68
    // `-in_pssm`: NCBI's checks of the PSSM (`CPsiBlastValidate::Pssm`) and the PSSM engine's
    // scores of a PSSM of frequency ratios only (`PsiBlastComputePssmScores`) fail when the search
    // object is made, after the outfmt 0 prolog (tblastn_app.cpp:270, 322) and before any
    // subject is read.
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:307,324
    // With `-export_search_strategy` the request is made first: with `-subject` it stops the run
    // ("Remote PSI-BL2SEQ is not supported", exit status 1, the strategy file empty; NCBI 2.17.0
    // `tbn.pssm.sub*` of line RG F1d probe 7) before the PSSM's check, which a request with `-db`
    // makes (no request then); the PSSM engine's error comes after the request is written
    // (`PsiSearch::stop_before_search`, psiblast's rule too).
    let pssm_stop = pssm
        .as_deref()
        .and_then(|psi| psi.stop_before_search(args.strategy.is_some(), database.is_some()));
    if let Some(stop) = pssm_stop {
        write_first_batch_prologs(
            outputs,
            subject_path.map(|path| path.as_path()),
            subjects.len(),
            subject_letters,
            pairwise_database.as_ref(),
            args.tabular_taxdb.unwrap_or(false),
            program,
            resolved.html.as_deref(),
            false,
        )?;
        // The request of the one PSSM query, which NCBI wrote before its constructor failed.
        if let (true, Some(strategy)) = (stop.after_request(), args.strategy.as_ref()) {
            super::archive::requests(
                database.is_none().then_some(subject_records),
                resolved.lcase_masking,
            )
            .with_pssm_query(
                pssm.as_deref()
                    .and_then(crate::algorithm::pssm_search::PsiSearch::request_pssm),
            )
            .write(strategy, 0, input_records, None)?;
            // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:270,307
            // NCBI's posts of the bytes replaced in that request (the `-db` name; NCBI 2.17.0
            // probe RG K/k4b/agent7/probe4 `tbn_dbe_m45_f0`: the prolog, the posts, the error),
            // saved in the error archive with outfmt 11 (`tbn_dbe_m45_f11`). The PSSM's own strings
            // hold none: NCBI's reader replaced them before the prolog with its own messages, and
            // LOSAT refuses such a PSSM (`in_pssm::psi_search_and_object`).
            strategy.write_posts(&mut *outputs.diagnostics, args.archive.as_ref())?;
        }
        return Err(stop.error());
    }
    // Under `-parse_deflines`, the scopes of the batches' queries (`QueryScopes`); with
    // `-in_pssm`, of the PSSM's query when the reports show its ids (`PsiSearch::query_record`).
    let believed_queries = crate::blastinput::query_batch::believed_queries(
        args.parse_deflines,
        resolved.pssm.is_some(),
        input_records,
    );
    // LOSAT's refusal of two queries of one scope with one id is the stop of their batch
    // (`QueryBatches::scope_refusal`), after the reports of the batches before unless NCBI stops at
    // one of them first (line WQA D10 round 4, 2026-11-20): with `-export_search_strategy` the
    // refusal (no request of the batch: `stop_batch` stays `None`), else NCBI's error of the
    // batch's set-up where it fails (`refused_batch_setup_error`), else the refusal.
    let mut batches = crate::blastinput::query_batch::read_query_batches(
        input_records,
        &searched_queries,
        &query_input,
        reading,
        batch_size,
        true,
        crate::blastinput::query_batch::QueryScopes::of_search(
            believed_queries,
            database.is_none(),
            believed_queries
                && crate::blastinput::query_batch::outfmts_have_xml(
                    outputs.formats.iter().map(|format| format.outfmt),
                ),
            "TBLASTN",
        )
        .map(|scopes| {
            scopes
                .with_refusal(args.scope_refusal.clone())
                .with_xml_matrix_failure(super::stage_e_report::xml_matrix_failure(
                    &outputs.formats,
                    &resolved.matrix_name,
                    database.is_some(),
                    pssm.is_some(),
                ))
        }),
    );
    // Where the run's XML failure for its matrix name takes a refusal (`ScopeRefusal::xml_takes`),
    // a batch whose set-up passes is searched (`QueryBatches::resolve_refusals`).
    batches.resolve_refusals(args.strategy.is_some(), |start, first| {
        refused_batch_setup_error(
            first,
            start,
            (input_records, &searched_queries, &query_input),
            &resolved,
            subject_set,
            engine_positives,
            database,
            subject_letters,
            (compo_adjust_mode, scoring),
        )
    });
    // NCBI crashes with IDENTITY's compressed lookup table once it searches a query batch, and
    // LOSAT's limits of the search apply where it searches (`check_search_limits`): a run that
    // stops at its first batch (or before it, `pssm_stop` above) gives NCBI's error.
    if batches.searched_done > 0 {
        super::args::check_identity_compressed_lookup_of(&resolved)?;
        super::args::check_search_limits(&resolved)?;
    }
    // The search strategy's requests of the batches (`-export_search_strategy`), with every
    // `-subject` record: written in the loop below, and for the batch that stops the run (every
    // query without letters) after the batches before it.
    // With `-in_pssm` every request holds the PSSM as read (`RunRequests::with_pssm_query`).
    let requests = args.strategy.as_ref().map(|_| {
        super::archive::requests(
            database.is_none().then_some(subject_records),
            resolved.lcase_masking,
        )
        .with_pssm_query(
            pssm.as_deref()
                .and_then(crate::algorithm::pssm_search::PsiSearch::request_pssm),
        )
    });
    if batches.searched_done == 0 && batches.stop.is_some() {
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
        // The first batch fails after the prologs (`formatter.PrintProlog()`; those of the first
        // report here). A first batch whose set-up fails has its search strategy request (NCBI's
        // posts of the bytes replaced in it follow the batch's reader messages, before its error,
        // below).
        write_first_batch_prologs(
            outputs,
            subject_path.map(|path| path.as_path()),
            subjects.len(),
            subject_letters,
            pairwise_database.as_ref(),
            args.tabular_taxdb.unwrap_or(false),
            program,
            resolved.html.as_deref(),
            false,
        )?;
        if let (Some(strategy), Some(requests), Some(span)) = (
            args.strategy.as_ref(),
            requests.as_ref(),
            batches.stop_batch.clone(),
        ) {
            requests.write(
                strategy,
                span.start,
                &searched_queries[span.clone()],
                query_input.placements.get(span.start),
            )?;
        }
    } else {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1426
        // `-db`: the database's sequences, checked once for every batch (`XDbStream::new`)
        // where NCBI reads them first, in the preliminary search of the first batch with a
        // valid context (a batch without one is not searched, `CheckInternalData`). A database
        // without sequences (BLASTN's `Database memory map file error`,
        // blastn/blast_engine/run.rs) or a sequence file that cannot be read fails that batch,
        // after the reports of the batches before, which are searched without sequences; a run
        // without such a batch never reads the file. The
        // database is checked as NCBI's first read would fail, with the same error, and each
        // sequence is read where it is searched (`common/x_db_stream.rs`).
        let mut x_stream = None;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-980
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_message.c:215-217
        // With the Gumbel block a database's `GetMinLength` below 3 (an alias file's
        // `MIN_SEQ_LENGTH` 1 or 2, a list's shortest sequence) leaves no residue
        // (`min_subject_length` is NCBI's `Uint4`): the setup of the first batch with a valid
        // context fails after the database's first OID is found and before any sequence is read.
        // NCBI 2.17.0 (probe of line RG D1a, `db.ak.tbn.min1.6`): `BLAST engine error: The
        // average subject length is too short`, exit status 3.
        // With `-in_pssm`, whether the read's failure comes before the batch's setup checks
        // (`PsiSearchReadFailure::precedes_setup_checks`: the empty database).
        let mut read_before_setup = false;
        let mut read_fails = None;
        let (database_stores, mut stores_error) = match database {
            Some(database) => {
                // `-in_pssm` (PSI-TBLASTN): `CPsiBlast` reads the database without `CLocalBlast`'s
                // catch, so its first read gives SeqDB's own error (`OID not found`, a sequence file
                // not found; `shared::PsiSearchReadFailure`, psiblast's rule), where a FASTA query
                // gets `Database memory map file error`.
                let (read, before_setup) = first_database_read(database, &resolved, &subject_set);
                read_before_setup = before_setup;
                match read {
                    Ok(stream) => {
                        x_stream = Some(stream);
                        (Vec::new(), None)
                    }
                    Err((error, class)) => {
                        read_fails = Some(class);
                        (Vec::new(), Some(error))
                    }
                }
            }
            None => (Vec::new(), None),
        };
        let search_letters = match (database, &x_stream) {
            (Some(_), Some(stream)) => SubjectLetters::Stream(stream),
            (Some(_), None) => SubjectLetters::Packed(&database_stores),
            (None, _) => subject_letters,
        };
        // Where NCBI's search first fails to read the database (the chunks of a split batch): the
        // class was decided with the error (`first_database_read`).
        let subject_set = SubjectSetNumbers {
            read_fails,
            ..subject_set
        };
        // The batches before the stop are searched and reported; a reader at the end of its
        // input reads no batch (its `End()` is true), and the reports have no query (an empty
        // pipe, or `-query -` after `-subject -`).
        let searched_done = batches.searched_done;
        let queries = &searched_queries[..searched_done];
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xmlformat.cpp:213-218 (the BLAST
        // XML report data of the first result NCBI formats: the matrix of its formatting table,
        // found by the name as typed)
        // NCBI fails when it formats a batch in BLAST XML with a matrix name its XML data does not
        // know (`stage_e_report::xml_matrix_failure`, the decision of `stage_e_report::render`):
        // the batch of the run's first searched query, or where the report is held (`-subject`,
        // `-in_pssm`; `XmlMatrixFailure::SecondQueryOrEnd`) of NCBI's next query of another
        // Seq-id (`QueryScopes::xml_failing_query` of the batches' scopes);
        // the search strategy's file holds the requests up to that
        // batch's, whatever the number of batches after it (probe 7 of RG-F1, probe K4f of line
        // RG). A held report that fails at the run's end, or a run that stops before, writes every
        // request it makes.
        if let Some(strategy) = args.strategy.as_ref() {
            if let Some(end) = super::stage_e_report::xml_matrix_failure(
                &outputs.formats,
                &resolved.matrix_name,
                database.is_some(),
                pssm.is_some(),
            )
            .and_then(|failure| batches.xml_matrix_request_end(failure))
            {
                strategy.stop_after(end);
            }
        }
        // Each searched query's whole record and where its letters start in it (`-query_loc`).
        let whole_records: Vec<(&[u8], usize)> = (0..queries.len())
            .map(|index| {
                (
                    input_records[query_input.ordinal(index)].seq(),
                    query_input.placements.offset(index),
                )
            })
            .collect();
        let query_seqs: Vec<_> = queries.iter().map(|record| record.seq().to_vec()).collect();
        let seg = resolved.seg.params();
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:172-188:
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1411-1417:
        // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:251-252,275-301:
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:70-127:
        // Search each accumulated batch independently; preliminary linking
        // and sum-statistics depend on its concatenated query contexts.
        let mut results = Vec::with_capacity(queries.len());
        let mut lengths: Option<LocalSubjectParameters> = None;
        let mut ungapped_karlin = Vec::with_capacity(queries.len());
        let mut query_validity = Vec::with_capacity(queries.len());
        let mut query_batch_skipped = Vec::with_capacity(queries.len());
        // `-db`: the stored titles read so far, and the batch that a header file or the
        // sequence file failed with its error.
        let mut titles = super::database::DatabaseTitles::new();
        let mut failed: Option<(std::ops::Range<usize>, anyhow::Error)> = None;
        // NCBI's fetches with `check_oid_exclusion` of the redo and the traceback, one check for
        // the search's batches (`database::fetch_excluded`), and their reads of headers.
        // The first format, whose prefetch each batch asks (`database::prefetch_format`).
        let prefetch_format = super::database::prefetch_format(&outputs.formats);
        let reads = super::database::CallReads::of(
            database,
            || query_input.batches(batch_size),
            prefetch_format,
        );
        // Graduated 2026-10-22 (was LOSAT_X_TBNSSIDE, LOSAT_X_TBNBATCH): the batches of the
        // loop below searched ahead (`x_search_batches`), each with its range; `None` when there
        // is none. A batch whose range does not match is searched by the loop.
        // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:275-291
        // Dispatch point: NCBI searches the query batches one after another (one CLocalBlast per
        // batch). x_search_batches searches them ahead (on one shared pool, or side by side) and
        // returns the reports in batch order; the loop below takes them from it.
        // A database whose sequences could not be read searches its batches one by one below.
        let x_batches = match stores_error {
            Some(_) => Vec::new(),
            None => query_input
                .batches(batch_size)
                .into_iter()
                .map(|batch| batch.searched)
                .filter(|searched| !searched.is_empty() && searched.end <= searched_done)
                .collect(),
        };
        // An error that may fail one batch while the batches were searched ahead
        // (`fails_one_batch`): they are searched again one at a time below, so that the batch
        // that fails does after the reports of the batches before.
        let x_searched_ahead = super::stage_d_pipeline::x_search_batches(
            &query_seqs,
            search_letters,
            x_batches,
            |range| LocalStageDProfile {
                max_intron_length: resolved.max_intron_length,
                seg: seg.as_ref(),
                soft_masking: resolved.soft_masking,
                mask_lowercase: resolved.lcase_masking,
                genetic_code: resolved.db_gencode,
                expect_value: resolved.evalue,
                max_target_seqs: resolved.hitlist_size,
                query_records: &whole_records[range.clone()],
                subject_set,
                engine_positives,
                culling_limit: resolved.culling_limit,
                best_hit: resolved.best_hit,
                subject_besthit: resolved.subject_besthit,
                pssm: pssm.as_deref(),
                fetch_excluded: reads.of_batch(&range),
                query_cov_hsp_perc: resolved.query_cov_hsp_perc,
            },
            compo_adjust_mode,
            resolved.sum_stats,
            scoring,
            resolved.num_threads,
        );
        let mut x_searched = match x_searched_ahead {
            Ok(searched) => searched.map(Vec::into_iter),
            Err(error) if fails_one_batch(&error) => None,
            Err(error) => return Err(error),
        };
        for range in query_input
            .batches(batch_size)
            .into_iter()
            .map(|batch| batch.searched)
            .filter(|searched| !searched.is_empty() && searched.end <= searched_done)
        {
            // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:277-280
            // The batch's search strategy request (`-export_search_strategy`): once its query
            // factory is made, before any error of the batch is acted on (`x_search_batches` may
            // have searched it ahead; only the order of the file and of what the run writes
            // matters, and the file is apart from the reports). Its `EffectiveSearchSpace` option
            // is the one the batch before left (`note_searched` below). NCBI's posts of the bytes
            // replaced in it wait for the batch's start (`query_posts` below, or before the
            // batch's error).
            if let (Some(strategy), Some(requests)) = (args.strategy.as_ref(), requests.as_ref()) {
                requests.write(
                    strategy,
                    range.start,
                    &queries[range.clone()],
                    query_input.placements.get(range.start),
                )?;
            }
            let (mut batch_results, batch_lengths, batch_karlin, batch_validity) = {
                // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:277-291
                // Dispatch point: a report searched ahead for this very range is used; otherwise
                // run_local_for_report_threads searches the batch here, which is the port of
                // CLocalBlast::Run for one batch.
                let x_searched_here = match x_searched.as_mut().and_then(Iterator::next) {
                    Some((x_range, report)) if x_range == range => Ok(report),
                    _ => run_local_for_report_threads(
                        &query_seqs[range.clone()],
                        search_letters,
                        LocalStageDProfile {
                            max_intron_length: resolved.max_intron_length,
                            seg: seg.as_ref(),
                            soft_masking: resolved.soft_masking,
                            mask_lowercase: resolved.lcase_masking,
                            genetic_code: resolved.db_gencode,
                            expect_value: resolved.evalue,
                            max_target_seqs: resolved.hitlist_size,
                            query_records: &whole_records[range.clone()],
                            subject_set,
                            engine_positives,
                            culling_limit: resolved.culling_limit,
                            best_hit: resolved.best_hit,
                            subject_besthit: resolved.subject_besthit,
                            pssm: pssm.as_deref(),
                            fetch_excluded: reads.of_batch(&range),
                            query_cov_hsp_perc: resolved.query_cov_hsp_perc,
                        },
                        compo_adjust_mode,
                        resolved.sum_stats,
                        scoring,
                        resolved.num_threads,
                    ),
                };
                // A read of the database that failed in this
                // batch fails it as the check's error would (below); so does NCBI's engine error
                // of its setup.
                match x_searched_here {
                    Ok(report) => {
                        // The batch's set-up has written its `EffectiveSearchSpace` option, which
                        // the next batch's search strategy request holds.
                        if let Some(strategy) = args.strategy.as_ref() {
                            let (_, lengths, _, validity) = &report;
                            strategy.note_searched(
                                range.clone(),
                                super::archive::batch_searched_space(
                                    validity.iter().any(|&valid| valid),
                                    &lengths.lengths[0],
                                    &queries[range.start],
                                    validity[0],
                                ),
                            )?;
                        }
                        report
                    }
                    // NCBI's setup reads OID 0 before its checks (`BLAST_GapAlignSetUp`, quoted
                    // at `PsiSearchReadFailure::precedes_setup_checks`): with `-in_pssm` an empty
                    // database fails the batch before the PSSM's Karlin-Altschul check or the
                    // average subject length would, as psiblast's
                    // (`psiblast_stop_at_invalid_query`).
                    Err(error)
                        if read_before_setup
                            && stores_error.is_some()
                            && error.is::<super::stage_d_pipeline::SetupEngineError>() =>
                    {
                        failed = Some((
                            range,
                            stores_error
                                .take()
                                .expect("the database's read failure before the setup"),
                        ));
                        break;
                    }
                    // `-db`: the first batch whose search reads the sequences (one with a valid
                    // context searched whole, a split one whose chunks get there) meets the read's
                    // failure; a reader error of NCBI's own fails it (the preliminary search's
                    // text, or with `-in_pssm` SeqDB's own error,
                    // `TblastnDatabase::pssm_search_read_failure`), LOSAT's refusal stops the run.
                    // (The search gives the marker only where the read fails: `stores_error`.)
                    Err(error) if error.is::<super::stage_d_pipeline::DatabaseReadReached>() => {
                        let error = stores_error.take().expect("the database's read failure");
                        failed = Some((range, database_read_failure(error)?));
                        break;
                    }
                    Err(error) => match batch_failure(error) {
                        Ok(failure) => {
                            reads.fetched_before_failure(&range);
                            failed = Some((range, failure));
                            break;
                        }
                        Err(error) => return Err(error),
                    },
                }
            };
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:283-293
            filter_final_hitlists(
                &mut batch_results,
                &query_seqs[range.clone()],
                &resolved,
                compo_adjust_mode,
                subject_set.sorts_results(),
            )?;
            // `-db`: the definition lines of the batch's final hit lists are read where NCBI
            // makes the batch's Seq-aligns (`database::read_titles`); a header file that NCBI
            // cannot read fails this batch, after the reports of the batches before.
            if database.is_some() {
                if let Some(error) = reads.titles(&mut batch_results, &range, &mut titles)? {
                    failed = Some((range, error));
                    break;
                }
            }
            // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:294-299
            // The prefetch before the batch's reports reads `PRE_FETCH_SEQS_LIMIT`, after the
            // batch's search and the read of its titles: a value that NCBI cannot convert fails
            // the first batch here, before its reports (`app::check_prefetch_environment`).
            if let Err(error) = app::check_prefetch_environment(outputs.only_archive()) {
                failed = Some((range, error));
                break;
            }
            // LOSAT's refusal of a reported sequence with the id of a query of its scope, after
            // NCBI's own errors of the batch: the reports would show it (the prefetch's error
            // above comes first).
            if let (Some(database), Some(scopes)) = (database, batches.query_scopes.as_ref()) {
                crate::algorithm::common::database::refuse_reported_ids_of_queries(
                    &titles,
                    batch_results
                        .iter()
                        .enumerate()
                        .flat_map(|(index, hitlist)| {
                            hitlist.lists().iter().filter_map(move |list| {
                                Some((range.start + index, u32::try_from(list.oid).ok()?))
                            })
                        }),
                    database.oid_map(),
                    scopes,
                    "TBLASTN",
                )?;
            }
            results.append(&mut batch_results);
            if let Some(all_lengths) = &mut lengths {
                // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:288-301 formats each batch's own
                // query-context lengths; retain those values in input order.
                all_lengths.lengths.extend(batch_lengths.lengths);
                all_lengths.cutoffs.extend(batch_lengths.cutoffs);
            } else {
                lengths = Some(batch_lengths);
            }
            ungapped_karlin.extend(batch_karlin);
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-224:
            // When `CheckInternalData` returns a nonzero status, each result in this `Run` batch
            // receives an unsearched ancillary state and null align set.
            let skipped = batch_validity.iter().all(|&valid| !valid);
            query_batch_skipped.extend(std::iter::repeat_n(skipped, batch_validity.len()));
            query_validity.extend(batch_validity);
        }
        // A run without a searched query has no parameters of a query set.
        let mut lengths = lengths.unwrap_or(LocalSubjectParameters {
            lengths: Vec::new(),
            cutoffs: Vec::new(),
            link: None,
            hit_cutoff_min: 0,
            word_cutoff_min: 0,
            prelim_evalue: 0.0,
        });
        // Whether the first format prefetches the subjects of each batch (outfmt 10 and 20 as the
        // tabular 6 and 7, below).
        let first_number = outputs
            .formats
            .first()
            .map(|format| super::stage_e_report::format_number(format.outfmt));
        let first_prefetches = first_number.is_some()
            && crate::algorithm::common::database::prefetches_subjects(Some(prefetch_format));
        // The XML2 and JSON formats (13 to 16), whose prefetch NCBI runs too but which write the
        // taxonomy texts at each `Hit`: of the prefetch only its stop
        // (`shared::held_batch_prefetch_stop`).
        let first_xml2 = first_number.is_some_and(|number| matches!(number, 13..=16));
        // The query records that NCBI's scope holds at the prefetch of the batch of the searched
        // queries `span` (`shared::PrefetchQueries`, made only where the batch stops): the PSSM's
        // query of `-in_pssm` (`PsiSearch::query_scope`), else the input records that NCBI's
        // reader added for the batch (`RangedBatch::input`: the records that `-query_loc` drops
        // included, numbered in the run, `FastaRecord::local_id`), and in outfmt 5, whose scope
        // the formatter keeps (`ResetScopeHistory`, cited at `shared::stop_lookup_unsure`), those
        // of the batches before.
        let xml_scope = crate::blastinput::query_batch::outfmts_have_xml(
            outputs.formats.iter().map(|format| format.outfmt),
        );
        // The batches are made here only (review 3 of line DBID loadord, M1).
        let batch_scope = |span: &std::ops::Range<usize>| {
            if let Some(psi) = pssm.as_deref() {
                return psi.query_scope();
            }
            let input = query_input
                .batches(batch_size)
                .into_iter()
                .find(|batch| batch.searched == *span)
                .map_or(0..0, |batch| batch.input);
            let input = if xml_scope { 0..input.end } else { input };
            crate::blastinput::query_batch::QueryScope::of_shown_records(&input_records[input])
        };
        // A batch whose prefetch loads a subject whose lookup fails on a file stops the run in
        // its prefetch (`database::prefetch_stop`): as a batch that a database file failed, it is
        // not reported, nor are the batches after it; its error is made after the prefetches of
        // the batches before (below), which look the taxonomy up first.
        let prefetch_stop = if database.is_some() && (first_prefetches || first_xml2) {
            let spans: Vec<std::ops::Range<usize>> = query_input
                .batches(batch_size)
                .into_iter()
                .map(|batch| batch.searched)
                .filter(|searched| !searched.is_empty() && searched.end <= results.len())
                .collect();
            super::database::prefetch_stop(
                &titles,
                &results,
                &spans,
                (prefetch_format, resolved.gapped),
            )
            .filter(|(span, _)| {
                failed
                    .as_ref()
                    .is_none_or(|(range, _)| span.start < range.start)
            })
        } else {
            None
        };
        if let Some((span, _)) = &prefetch_stop {
            failed = Some((
                span.clone(),
                anyhow::anyhow!("the prefetch of the batch stops (its error is made below)"),
            ));
        }
        // A batch that a database file failed is not reported, nor are the batches after it.
        let queries = match &failed {
            Some((range, _)) => {
                results.truncate(range.start);
                lengths.lengths.truncate(range.start);
                lengths.cutoffs.truncate(range.start);
                ungapped_karlin.truncate(range.start);
                query_validity.truncate(range.start);
                query_batch_skipped.truncate(range.start);
                &queries[..range.start]
            }
            None => queries,
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1443-1451:
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:632-640
        // The warnings of each query are written before its report, once (`QueryWarnings`):
        // the reader's messages of the query's batch (before its first query,
        // `read_query_batches`), then the query's own, on one line in NCBI's order
        // (`query_warning`): the O characters read as X, and the message of its set-up (a
        // query without letters) or of its Karlin-Altschul parameters (an invalid query
        // with letters).
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1431-1452
        // A structured format (outfmt 8, 12) writes none of a query's own warnings
        // (`ReportFormat::writes_query_warnings`); the warnings go with the first format, or, where
        // it is structured, they are kept by query for the later formats that write them
        // (`later_lines`, `Posts::later`, `ReportOutputs::later_warnings`).
        let result_warnings = outputs.first_writes_query_warnings();
        let later = outputs.keeps_later_query_warnings();
        let mut later_lines: Vec<Vec<u8>> = Vec::new();
        let filtered_empty = database.is_some_and(|database| database.filtered_empty());
        // The texts of the subjects that the reports of the queries whose batch's prefetch is
        // off read one by one, by query, where the first format places them (`places_reads`).
        let mut subject_reads = crate::report::query_warnings::SubjectReads::default();
        // Outfmt 5 first: whether NCBI's XML stream cuts a query's report before its hits
        // (`shared::XmlCuts`), the one answer for the record of the batches' loads before a
        // stopping batch (`LoadedBefore`: a cut report loads nothing) and for the texts that the
        // cut reports left (`shared::XmlReportsBefore`, from the reads of the queries' reports
        // modelled here).
        let xml_first = outputs
            .formats
            .first()
            .and_then(crate::api::local_blast::FormatOutput::report_format)
            == Some(crate::blastinput::app::ReportFormat::Xml);
        let xml_cut = |query: usize| {
            queries.get(query).is_some_and(|record| {
                crate::report::blast_xml::names::iteration_cut_before_hits(
                    record.local_id.as_bytes(),
                    record.warning_id(),
                    record.believed(),
                    record.title(),
                    true,
                )
            })
        };
        let xml_cuts = crate::algorithm::common::database::XmlCuts::of(xml_first, &xml_cut);
        // The posts of each batch's prefetch, with its first query (`query_warnings::query_posts`).
        let mut prefetch_posts = Vec::new();
        let mut query_warning_lines: Vec<Vec<u8>> = queries
            .iter()
            .zip(&query_validity)
            .enumerate()
            .map(|(index, (query, valid))| {
                let mut lines = batches.warnings[index].clone();
                if !result_warnings && !later {
                    return lines;
                }
                let messages = crate::report::query_warnings::query_messages(
                    query.seq(),
                    *valid,
                    index < identity_warned,
                    filtered_empty && !query_batch_skipped[index],
                );
                let text = match pssm.as_deref() {
                    // A PSSM query's label is its Seq-id and the title of NCBI's query source of
                    // the PSSM's Bioseq (`PsiSearch::warning_title`), as psiblast's
                    // (`blastp::blast_engine::pssm_hooks::pssm_query_warning`).
                    Some(psi) => crate::report::query_warnings::query_warning_with_seq_id(
                        "tblastn",
                        psi.query.warning_id(),
                        psi.warning_title.as_deref().unwrap_or(&psi.query.title),
                        &messages,
                    ),
                    None => crate::report::query_warnings::query_warning(
                        "tblastn",
                        args.first_query_ordinal + query_input.ordinal(index),
                        query,
                        &messages,
                    ),
                };
                match result_warnings {
                    true => lines.extend(text),
                    false => later_lines.push(text),
                }
                lines
            })
            .collect();
        // The query batches that were read and searched, by their searched queries (the batches
        // before a database failure or a reader stop): each is reported.
        let batch_spans: Vec<std::ops::Range<usize>> = query_input
            .batches(batch_size)
            .into_iter()
            .map(|batch| batch.searched)
            .filter(|searched| !searched.is_empty() && searched.end <= queries.len())
            .collect();
        // `-db_hard_mask` (series NP, M4b PR 2): the search read `B` over the masked bases, and the
        // reports count again on the stored ones (`stage_e_report::recounts_shown_counts`).
        let hard_masked = resolved
            .db_mask
            .as_ref()
            .is_some_and(|mask| mask.kind == crate::blastinput::db_mask::MaskKind::Hard);
        // `-db`: the reported sequences by rank, with their records and lengths, and their packed
        // subjects where a report reads their letters; the database's other sequences are dropped
        // before the reports.
        let reported = match database {
            Some(database) => Some(
                super::database::reported_subjects(
                    database,
                    &titles,
                    &mut results,
                    search_letters,
                    super::stage_e_report::reads_subject_letters(
                        outputs,
                        pssm.as_deref(),
                        hard_masked,
                    ),
                    &outputs.formats,
                    (outputs.other_outfmts, resolved.html.is_some()),
                    pairwise_reading(&resolved, &outputs.formats),
                    &batch_spans,
                    resolved.gapped,
                )
                // A read of a reported sequence that
                // failed: the run's error, as the check's error would be.
                .map_err(|error| batch_failure(error).unwrap_or_else(|error| error))?,
            ),
            None => None,
        };
        // The library's records name a database's sequences through this event.
        if let (Some(reported), Some(database), true, Some(observer)) = (
            &reported,
            database,
            outputs.hits.is_some(),
            outputs.observer.as_deref_mut(),
        ) {
            observer.database_subjects(reported.library_subjects(database));
        }
        drop(database_stores);
        // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:185-203 (cited at
        // `db_subject::RunLoads`)
        // A report whose load ends the run ends NCBI's application there: no batch after it is
        // prefetched, so a stop of a later batch's prefetch is never reached
        // (`HitLoads::reaches`; review 3 of line DBID pfoff-3, N1). The writer ends the run at
        // that report with its error, before `failed` (the stop's, above) is asked.
        let reaches = |query: usize| {
            reported
                .as_ref()
                .is_none_or(|reported| reported.loads.reaches(query))
        };
        let prefetch_stop = prefetch_stop.filter(|(span, _)| reaches(span.start));
        // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:298-302
        // Before the reports of a batch NCBI fetches every subject that the batch reports, and
        // with `taxdb.bti` on the search path the taxonomy database writes its messages there
        // (`TblastnDatabase::prefetch_messages`: the texts of opening it at the run's first
        // lookup, a volume's lookup of a taxid once per run): after the reader's messages of the
        // batch and before the warnings of its first query. They go with the warnings, where the
        // first format is outfmt 0, 5, 6 or 7, or with `PRE_FETCH_SEQS_LIMIT=2147483647` 8, 9 or
        // 12, which write the reader's messages only (`database::prefetches_subjects`); XML2 and JSON
        // write theirs at each `Hit` (`TaxdbPrinter`). They are read here, after the refusals of
        // the reported subjects and before any format is written, so a refusal (a taxonomy
        // database or a taxid that LOSAT cannot read as NCBI does) comes before any byte, as
        // those refusals do; without `taxdb.bti` nothing is read. A failed batch is not reported
        // and prefetches nothing (`failed`, a header file or the sequence file that NCBI cannot
        // read during the batch's search).
        // Which rows of the first format read the subjects' Bioseqs where a batch's prefetch is off
        // (`TabularSpec::subject_fetch`; every row of a format that is not a field list): what the
        // reports of such a batch load (`shared::batch_loads`).
        let first_rows = outputs
            .formats
            .first()
            .and_then(|format| super::tabular::tabular_spec(format.outfmt))
            .map_or(SubjectFetch::Every, |spec| {
                spec.subject_fetch(crate::blastinput::app::Program::Tblastn)
            });
        if let (Some(reported), Some(database)) = (&reported, database) {
            // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:716-725
            // Outfmt 10 and 20 prefetch as the tabular 6 and 7 do (a name of the taxonomy
            // database that the prefetch would look up is LOSAT's refusal, `tabular`).
            // Where the prefetch fails at a subject with parsed Seq-ids, NCBI writes its warning
            // after those messages (`TblastnDatabase::batch_prefetch`), with or without
            // `taxdb.bti`.
            let taxdb_found = database.taxdb_found();
            let xml2_reads = (reported, &mut subject_reads);
            let after = super::database::after_throw(
                &outputs.formats,
                (outputs.other_outfmts, resolved.html.is_some()),
            );
            prefetch_posts.extend(reads.xml2(
                first_xml2,
                xml2_reads,
                &results,
                &batch_spans,
                after.run(),
            )?);
            if first_prefetches && (taxdb_found || database.parsed_seqids()) {
                // A custom tabular report with a plural name looks every taxid of its rows up,
                // taxid 0 too, which the prefetch does not: the rows of a batch with hits may open
                // the taxonomy database before a later batch's prefetch
                // (`TaxdbRun::open_for_rows`).
                let first_spec = outputs
                    .formats
                    .first()
                    .and_then(|format| super::tabular::tabular_spec(format.outfmt));
                // Whether the first format's writer writes the texts of the subjects that a
                // query's report reads one by one where it reads each
                // (`ReportFormat::places_subject_reads`).
                let first_report = outputs
                    .formats
                    .first()
                    .and_then(crate::api::local_blast::FormatOutput::report_format);
                let places_reads = first_report
                    .is_some_and(crate::blastinput::app::ReportFormat::places_subject_reads);
                // The pairwise report (outfmt 0 to 4) reads its subjects in the order of its
                // `-sorthits`.
                let pairwise = first_report
                    .and_then(crate::report::pairwise::PairwiseView::of)
                    .and(pairwise_reading(&resolved, &outputs.formats));
                let plural_names = taxdb_found
                    && first_spec
                        .as_ref()
                        .is_some_and(crate::report::tabular::TabularSpec::looks_up_every_taxid);
                let rows_open = |hitlists: &[super::stage_d_results::KappaResultHitList]| {
                    if plural_names && reported.loads_a_row(hitlists) {
                        database.open_for_tabular_rows()?;
                    }
                    Ok::<(), anyhow::Error>(())
                };
                for span in &batch_spans {
                    // No batch after the run's end is prefetched, nor read one by one.
                    if !reaches(span.start) {
                        break;
                    }
                    reads.reach(span);
                    let hitlists = &results[span.clone()];
                    // The batch's Seq-aligns, one for each HSP of a gapped search, one for each
                    // subject of an ungapped one (`database::seq_aligns_of`).
                    let seq_aligns = super::database::seq_aligns_of(hitlists, resolved.gapped);
                    // Where `PRE_FETCH_SEQS_LIMIT` turns the batch's prefetch off
                    // (`shared::prefetch_off`), each query's report fetches the subjects it
                    // shows, after the query's own warnings, where the first format's rows read
                    // their Bioseqs (`TabularSpec::subject_fetch`; for `nident` alone the rows
                    // with identities, `database::reads_subject`). NCBI looks the taxonomy up
                    // only where it fetches a subject: in such a report, or in the prefetch of a
                    // batch with Seq-aligns (one without them loads no id, below). Every batch's
                    // lookups are made here before any report is written, so the run's refusal
                    // of NCBI's taxonomy cache (`shared::TaxdbRun::cache_collision`) is asked
                    // only there: a later batch that fetches still refuses before any byte.
                    if crate::algorithm::common::database::prefetch_off(prefetch_format, seq_aligns)
                    {
                        // No prefetch; no fetch writes without `taxdb.bti` or such a lookup.
                        if !taxdb_found && !database.scope_utf8_seen() {
                            continue;
                        }
                        let rows = first_rows;
                        for (query, hitlist) in span.clone().zip(hitlists) {
                            if !reaches(query) {
                                break;
                            }
                            let shown = std::slice::from_ref(hitlist);
                            if hitlist
                                .lists()
                                .iter()
                                .any(|list| super::database::reads_subject(list, rows))
                            {
                                // The first format writes the texts where it reads each subject
                                // (`places_subject_reads`), else before the report.
                                if places_reads {
                                    subject_reads.set(
                                        query,
                                        database.subject_read_texts(
                                            reported,
                                            shown,
                                            rows,
                                            plural_names,
                                            pairwise,
                                        )?,
                                    );
                                } else {
                                    let messages =
                                        database.prefetch_messages(reported, shown, rows)?;
                                    query_warning_lines[query].extend(messages);
                                }
                            }
                            rows_open(shown)?;
                        }
                        continue;
                    }
                    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:768-772
                    // A batch without Seq-aligns prefetches no sequence.
                    if seq_aligns == 0 {
                        continue;
                    }
                    let prefetched = |query: usize| {
                        prefetch_posts
                            .iter()
                            .find(|(first, _)| *first == query)
                            .map(|(_, posts)| posts)
                    };
                    let xml_before = xml_cuts.reports_before(subject_reads.queries(), &prefetched);
                    let scope = || batch_scope(span);
                    // The batches before (`results`): what NCBI's data loader loaded for them
                    // (`super::database::LoadedBefore`), read only where this batch stops; a call
                    // of the library's batches after the first (`api::run_batches`) has those of
                    // the calls before unrecorded.
                    let messages = database.batch_prefetch(
                        reported,
                        super::database::LoadedBefore {
                            hits: &results[..span.start],
                            spans: &batch_spans,
                            format: prefetch_format,
                            rows: first_rows,
                            cuts: xml_cuts,
                        },
                        hitlists,
                        &crate::algorithm::common::database::PrefetchQueries::new(
                            span.start, &scope,
                        )
                        .after_xml_reports(xml_before.as_ref())
                        .where_reports_unchanged_by_load_throws(
                            crate::report::db_subject::unchanged_by_load_throws_outfmts(
                                outputs.formats.iter().map(|format| format.outfmt),
                                crate::blastinput::app::Program::Tblastn,
                                crate::report::archive::options::RequestKind::Tblastn,
                            ),
                        ),
                    )?;
                    rows_open(hitlists)?;
                    // After the reader's messages of the batch, in its first query
                    // (`query_warnings::query_posts`).
                    if !messages.is_empty() {
                        prefetch_posts.push((span.start, messages));
                    }
                }
            }
        }
        // The stopping batch's prefetch, after those of the batches before: NCBI's error, or
        // LOSAT's refusal before any byte.
        if let (Some((span, subjects)), Some(database)) = (prefetch_stop, database) {
            reads.reach(&span);
            // The batches before are those of `results`, which the stop truncated (above).
            let prefetched = |query: usize| {
                prefetch_posts
                    .iter()
                    .find(|(first, _)| *first == query)
                    .map(|(_, posts)| posts)
            };
            let xml_before = xml_cuts.reports_before(subject_reads.queries(), &prefetched);
            let scope = || batch_scope(&span);
            let error = database.prefetch_stop_error(
                reported.as_ref(),
                super::database::LoadedBefore {
                    hits: &results,
                    spans: &batch_spans,
                    format: prefetch_format,
                    rows: first_rows,
                    cuts: xml_cuts,
                },
                &subjects,
                &crate::algorithm::common::database::PrefetchQueries::new(span.start, &scope)
                    .after_xml_reports(xml_before.as_ref()),
                first_xml2,
            );
            match error.downcast::<crate::algorithm::common::database::HeaderReadError>() {
                Ok(failure) => failed = Some((span, failure.error)),
                Err(error) => return Err(error),
            }
        }
        // `-db`, XML2 and JSON: the taxid of each reported sequence, and what the taxonomy
        // database says about it (the refusals of both come before any byte).
        let xml_taxa = match (&reported, database) {
            (Some(reported), Some(database))
                if outputs
                    .formats
                    .iter()
                    .any(|format| super::stage_e_report::xml2_kind(format.outfmt).is_some()) =>
            {
                database.read_taxids(reported, &results)?
            }
            _ => Default::default(),
        };
        // The custom tabular reports (`tabular`): whether `taxdb.bti` was found, the taxids of the
        // reported database sequences where a field shows them and the taxonomy database's names
        // that the rows look up (read here, before any byte, so that their refusals come before
        // it), and the `-query_loc` range's length (`qlen`, `qcovs`).
        let tabular_specs: Vec<_> = outputs
            .formats
            .iter()
            .filter_map(|format| super::tabular::tabular_spec(format.outfmt))
            .collect();
        let taxonomy = match (&reported, database) {
            (Some(reported), Some(database))
                if tabular_specs.iter().any(|spec| spec.needs().taxonomy) =>
            {
                database.tabular_taxonomy(reported, &results)?
            }
            _ => Vec::new(),
        };
        // A `-subject` search's organism report (outfmt 18) looks names up in the run's
        // database.
        let names = match (tabular_specs.first(), database) {
            (Some(spec), Some(database)) => database.tabular_names(spec, &taxonomy)?,
            (Some(_), None) => args
                .subject_taxdb
                .as_ref()
                .map(|taxdb| taxdb.report_names("TBLASTN")),
            _ => None,
        };
        let tabular_inputs = super::tabular::TabularInputs {
            taxdb_found: args.tabular_taxdb.unwrap_or(false),
            taxonomy,
            names,
            range_length: resolved
                .query_range
                .as_ref()
                .map(crate::blastinput::seq_range::SequenceRange::length),
            ungapped: args.ungapped,
            html: resolved.html.clone(),
            is_html: args.html,
        };
        // The records of a subject set whose ids `-parse_deflines` read, by subject index.
        let records_parsed = if args.parse_deflines && reported.is_none() {
            crate::report::db_subject::parsed_records(
                subjects
                    .iter()
                    .map(|record| (record.believed(), record.title.as_slice())),
                args.show_gis,
            )
        } else {
            Vec::new()
        };
        // LOSAT's refusals of these records in outfmt 0 to 4, each such format's
        // (`FormatRefusal`).
        if !records_parsed.is_empty() {
            for format in outputs.formats.iter_mut() {
                let number = crate::blastinput::app::parse_formatting_string(format.outfmt)?.number;
                if (0..=4).contains(&number) {
                    format
                        .refusal
                        .check(crate::report::db_subject::check_pairwise_records(
                            &records_parsed,
                            false,
                            "TBLASTN",
                        ))?;
                }
            }
        }
        let (report_records, report_letters, shown) = match (&reported, database) {
            (Some(reported), Some(database)) => (
                reported.records.as_slice(),
                reported.letters(),
                ShownSubjects::Database {
                    name: database.name(),
                    pairwise: pairwise_database
                        .as_ref()
                        .expect("a database search has its Database: lines"),
                    reported,
                    oids: database.oid_map(),
                    hard_masked,
                },
            ),
            _ => (
                subjects,
                subject_letters,
                ShownSubjects::Fasta {
                    path: subject_path.expect("a search of -subject has its path"),
                    parsed: &records_parsed,
                },
            ),
        };
        // The BLAST archive (outfmt 11): one document per query batch that was read and searched
        // (`batch_spans`).
        let archive_render = args
            .archive
            .as_ref()
            .map(|handle| super::archive::ArchiveRender {
                handle,
                spans: &batch_spans,
            });
        // The posts before each query's report: NCBI's posts of the bytes replaced in the search
        // strategy request of each reported batch follow the reader's messages of the batch, in
        // its first query (`query_warnings::query_posts`).
        let query_posts = crate::report::query_warnings::query_posts(
            query_warning_lines,
            args.strategy
                .as_ref()
                .map_or_else(Vec::new, |strategy| strategy.take_posts(queries.len())),
            prefetch_posts,
            |query| batches.warnings[query].len(),
            reads.header_read(queries.len()),
        );
        // NCBI reference (598d8ae6): blast_format.cpp:1450-1452 (at `QueryWarnings::with_later`).
        let query_posts = crate::report::diag_log::Posts::with_later(query_posts, &later_lines);
        // Diagnostic (LOSAT_TIMING): the report for the common line (utils/xtiming.rs); it is
        // written after all batches are searched. No NCBI counterpart.
        let t_output = crate::utils::xtiming::clock();
        render(
            outputs,
            queries,
            report_records,
            report_letters,
            &shown,
            subject_set,
            &mut results,
            &lengths,
            &ungapped_karlin,
            &query_validity,
            &query_batch_skipped,
            scoring,
            &resolved.matrix_name,
            crate::blastinput::report_counts::ReportCounts {
                hitlist_size: resolved.hitlist_size,
                num_descriptions: resolved.num_descriptions,
                num_alignments: resolved.num_alignments,
            },
            resolved.formatting,
            resolved.db_gencode,
            seg.as_ref(),
            resolved.seg.filter_options(),
            resolved.lcase_masking,
            &query_posts,
            subject_reads.queries(),
            &super::stage_e_report::ReportRanges {
                queries: &query_input.placements,
                subjects: subject_placements,
                // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:342
                // A run that stops at a query batch writes no epilog.
                epilog: batches.stop.is_none() && failed.is_none(),
                query_scopes: batches.query_scopes.as_ref(),
            },
            super::blast_xml::XmlValues {
                expect: resolved.evalue,
                soft_masking: resolved.soft_masking,
                gapped: resolved.gapped,
                // IDENTITY's reset is NCBI's option of the search (`search_settings`).
                cbs: if identity_cbs_reset {
                    0
                } else {
                    resolved.composition_mode as i32
                },
                db_gencode: i32::from(resolved.db_gencode),
                identity_warned,
                filtered_empty,
                pssm: pssm.is_some(),
            },
            xml_taxa,
            args.out.as_deref(),
            archive_render.as_ref(),
            &tabular_inputs,
            pssm.as_deref(),
            args.show_gis,
            &batches
                .starts
                .iter()
                .map(|&(_, searched)| searched)
                .collect::<Vec<_>>(),
        )?;
        crate::utils::xtiming::add(crate::utils::xtiming::Phase::Output, t_output);
        // The failed batch's reading messages, NCBI's posts of the bytes replaced in its search
        // strategy request (written before its search), then its error.
        if let Some((range, error)) = failed {
            for warnings in &batches.warnings[range.clone()] {
                outputs.diagnostics.write_all(warnings)?;
            }
            if let Some(strategy) = args.strategy.as_ref() {
                strategy.write_posts(&mut *outputs.diagnostics, args.archive.as_ref())?;
            }
            reads.write_failed(&range, &mut *outputs.diagnostics, args.archive.as_ref())?;
            return Err(error);
        }
    }
    // A batch whose set-up fails (every query without letters) has its search strategy request:
    // NCBI made its query factory before the search (`QueryBatches::stop_batch`), after the
    // batches above (a first batch's was written after the prologs).
    if let (Some(strategy), Some(requests), Some(span)) = (
        args.strategy.as_ref(),
        requests.as_ref(),
        batches.stop_batch.clone(),
    ) {
        requests.write(
            strategy,
            span.start,
            &searched_queries[span.clone()],
            query_input.placements.get(span.start),
        )?;
    }
    // The batch's reader messages, NCBI's posts of the bytes replaced in its request, its error.
    if let Some((warnings, error)) = batches.stop {
        outputs.diagnostics.write_all(&warnings)?;
        if let Some(strategy) = args.strategy.as_ref() {
            strategy.write_posts(&mut *outputs.diagnostics, args.archive.as_ref())?;
        }
        return Err(error);
    }
    Ok(())
}

#[cfg(test)]
#[path = "tblastn_app_tests.rs"]
mod tests;
