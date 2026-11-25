//! Native serial caller, tested through the internal probe before public dispatch.
#[cfg(feature = "parallel")]
use super::reader_batch::{
    attach, prepare_lone_unit, prepare_reader_batch, prepare_unit, ReaderBatchShape, UnitPrep,
};
use super::{
    archive,
    args::BlastxArgs,
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
    input::{
        read_fasta_stream_batches, stream_is_empty, FastaParseError, InputWarning, NoDeflineError,
        QueryStreamSettings,
    },
    query_setup::BATCH_SIZE,
    report,
    runtime::search_internal,
};
use super::{
    database::{BlastxDatabase, ReportedSubjects},
    input::FastaRecord,
    report::ShownSubjects,
    runtime::BatchResults,
    subjects::{Subjects, UnreadableSubjects},
};
use crate::algorithm::common::database::{
    prefetches_subjects, HeaderReadError, PreliminaryReadError,
};
use crate::report::blast_xml2::{FsSideFiles, NoSideFiles, RunSideFiles, Xml2Kind};
use crate::wording::{DATABASE_ERROR, ENGINE_ERROR};
use anyhow::{bail, Result};
use std::{
    fmt,
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
    io::{self, Write},
    path::Path,
};

// The NCBI-style exit error is shared with the other programs' command lines.
pub use crate::cli::{exit_on_native_error, NativeError};

// Parts of `run_inner` (the set-up of the options, the read of the `-subject` file, the archive of
// a query batch): a child module, so that it calls the private functions of this file.
mod run_parts;
use run_parts::{read_subjects, ArchiveBatch, RunSetup, SubjectRead};

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:225-227
pub(crate) fn engine(message: impl fmt::Display) -> anyhow::Error {
    NativeError {
        exit: 3,
        message: format!("{ENGINE_ERROR} {message}\n").into(),
    }
    .into()
}

/// The error of the database's sequences (`BlastxDatabase::x_check_stores`, or a read of
/// `x_db_stream.rs`): NCBI's engine error of a reader error, or LOSAT's refusal as it is.
fn database_error(kept: anyhow::Error) -> anyhow::Error {
    match kept.downcast::<PreliminaryReadError>() {
        Ok(PreliminaryReadError(text)) => engine(text),
        Err(refusal) => refusal,
    }
}
/// A reader error of the input as the application reports it (`CObjReaderParseException`:
/// `blastinput::app::query_error`, the one producer of that line); any other error as it is.
pub(super) fn input_error(error: anyhow::Error) -> anyhow::Error {
    let text = if let Some(e) = error.downcast_ref::<FastaParseError>() {
        e.to_string()
    } else if let Some(e) = error.downcast_ref::<NoDeflineError>() {
        e.to_string()
    } else {
        return error;
    };
    crate::blastinput::app::query_error(&text)
}

/// The error that ends the read of the `-subject` file (`ReadSequencesToBlast`, NCBI's
/// `CBlastInput::GetAllSeqs`): the reader's `eNoDefline` is `GetAllSeqs`'s input error
/// (`fasta_reader::no_defline_error`), every other error is `input_error`'s.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:205-214
pub(super) fn subject_input_error(error: anyhow::Error) -> anyhow::Error {
    match error.downcast_ref::<NoDeflineError>() {
        Some(NoDeflineError(line)) => crate::blastinput::fasta_reader::no_defline_error(*line),
        None => input_error(error),
    }
}
/// The error of writing one batch's archive: the encoder's refusal of a byte NCBI writes as `#`
/// (`SerialError`; a document with such a byte is written through `Traced`, `write_archive`, so
/// only a string that its check does not see reaches it) is LOSAT's refusal, not an error of
/// NCBI's engine, so it gets no `BLAST engine error:` and no error archive (`archive::error_lines`
/// takes `NativeError`s only), as `XmlRefusal` gets none from `report_error`.
pub(super) fn archive_error(error: anyhow::Error) -> anyhow::Error {
    if error.is::<crate::report::serial::encoder::SerialError>() {
        error
    } else {
        report_error(error)
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:225-227
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:251-255
pub(super) fn report_error(error: anyhow::Error) -> anyhow::Error {
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:160-163
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:118-119
    if error.is::<super::blast_xml::XmlRefusal>()
        || error.is::<crate::algorithm::common::database::ReportSumRefusal>()
        || crate::report::pairwise::is_total_score_refusal(&error)
        || error.is::<crate::report::tax_report::TaxReportRefusal>()
        || error.is::<crate::report::tabular::FailedLookupRefusal>()
        || error.is::<crate::algorithm::common::database::ParsedRefusal>()
        || error.is::<super::empty_rows::EmptyRowsRefusal>()
        || error.is::<crate::algorithm::common::database::FetchRefusal>()
    {
        // LOSAT's refusal of a title NCBI's XML crashes on, of the database counts a report
        // would print where NCBI's sums overflow (`ReportSumRefusal`), of a `-sorthits` table whose
        // last total NCBI pads without end (`TotalScoreRefusal`, an `io::Error` of the report's
        // writes), of an organism report (outfmt 18, made before its first byte), of a tabular
        // row at which NCBI stops (`report::tabular::check_failed_lookups`), of a database
        // sequence that a query-anchored report (outfmt 1 to 4) would show as a query
        // (`ParsedRefusal`, `report::write_pairwise`), or of a report of an ungapped HSP of
        // length 0 in a form LOSAT does not write (`empty_rows::EmptyRowsRefusal`): not an error
        // of NCBI's engine or of its output stream.
        error
    } else if error.downcast_ref::<NativeError>().is_some_and(|e| {
        (e.exit == 1 && e.message.starts_with(b"Command line argument error"))
            || (e.exit == 6 && e.message == b"BLAST failed to write output\n")
            || (e.exit == 3 && e.message.starts_with(ENGINE_ERROR.as_bytes()))
            || (e.exit == 2 && e.message.starts_with(DATABASE_ERROR.as_bytes()))
            || e.exit == 255
            || e.exit == 0
    }) {
        // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:185-185
        // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:199-202
        // An outfmt 13 or 14 side file NCBI cannot open (`CArgException`, exit status 1), a
        // write to one that fails (`framing::side_file_write_error`; on the command line NCBI's
        // write fault, which `CATCH_ALL` swallows: exit status 0, `ostream::write_side_file`),
        // NCBI's engine error made
        // already (`blast_xml::matrix::unsupported_score_matrix`), or any exception of NCBI's
        // that a report raises and that ends the program through `CATCH_ALL`'s `CException`
        // branch (`NativeError::exception`, exit status 255, blast_app_util.hpp:247-250; such as
        // an exception of NCBI's XML stream, `blast_xml::writer::run_ending_exception`, or the
        // formatters' exceptions on an ungapped HSP of length 0, `empty_rows`) or through its
        // `CSeqDBException` branch, which comes first (exit status 2: the lookup of a reported
        // subject that throws in an outfmt 13 or 14 report, `Xml2Document::failed`), which is
        // not an engine error.
        error
    } else if error
        .downcast_ref::<io::Error>()
        .is_some_and(crate::report::ostream::is_report_write_failure)
    {
        // The report stream's failure: NCBI's end is the stream's (`ReportStream::finish`).
        error
    } else if let Some(native) = NativeError::of_writer_exception(&error) {
        // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:248-251
        // An exception of NCBI's that a report's writer raised (`NcbiException::into_io_error`:
        // the `-html` registries, `formatter_registry::FormatterNcbirc` in the prolog and
        // `InitConfig` at the first id link; a subject's load that ends the run) is `CATCH_ALL`'s
        // `CException` report, exit status 255 (as every program's `NativeError::of_writer`), not
        // a failed write.
        native.into()
    } else if error.downcast_ref::<io::Error>().is_some() {
        // Any other error carried as an `io::Error` is LOSAT's own error (an internal error of a
        // formatter, `report::query_anchored::internal`): only the report stream's failure (a
        // closed pipe's too) is NCBI's "BLAST failed to write output" (line RG T3c review 3 N4).
        error
    } else {
        engine(error)
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:62-70
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2537-2537
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2554-2557
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3468-3468
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3479-3479
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:198-200
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:213-217
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:254-254
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:260-297
fn run_inner(
    args: &BlastxArgs,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
    output: &mut Option<crate::cli::NcbiOstream<'static>>,
    archive_slot: &mut Option<crate::report::archive::ArchiveHandle>,
    strategy_slot: &mut Option<crate::report::archive::strategy::StrategyExport>,
) -> Result<()> {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1031-1044
    // NCBI's application layer reads its environment and registry files before blastx's own
    // code runs; LOSAT rejects the settings that change the output first, before any option or
    // file of the run, in `main` as for its other programs: `settings` are those of the run's
    // application start. The registry's `BLASTDB` finds `-db`, and the data loaders decide
    // whether the reader tries the first line of the query or the `-subject` file as a Seq-id
    // (`FastaStream::seq_id_role`).
    // An RPSTBLASTN run (`BlastxArgs::rpstblastn`) is the NCBI application rpstblastn: its
    // settings, messages and query batches are that program's (`args::ProgramName`).
    let program = args.program();
    let id_scope: Option<std::sync::Arc<crate::seqid::id_scope::IdScope>> =
        args.parse_deflines.then(Default::default);
    // Diagnostic (LOSAT_TIMING): the clock of the run (x_timing.rs). No NCBI counterpart.
    let x_t_run = super::x_timing::lap();
    use super::x_timing::Run as XRun;
    use crate::blastinput::app::{self, Program};
    let RunSetup {
        choice,
        archive,
        strategy,
        resolution,
    } = RunSetup::open(args, &settings, program, archive_slot, strategy_slot)?;
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:727-735
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
    // `-db`: the database opens after the options are processed (below).
    let subject_path = args.subject.as_deref();
    // NCBI opens the `-subject` and `-query` files inside `Run`, after the diagnostics handler
    // that saves the messages is set: an error archive follows their errors.
    let early_archive_error = |error: anyhow::Error| -> anyhow::Error {
        write_run_error_archive(args, archive.as_ref(), &error);
        error
    };
    // `-subject -` is `cin` (ncbiargs.cpp:717-721, cited at `input_files::open_input`), the
    // stream that `-query -` names too: one owner opens both (`open_input`).
    let subject_file = match (&args.db, subject_path) {
        (Some(_), _) => None,
        (None, Some(path)) => Some(
            crate::blastinput::input_files::open_input(path, "subject", program.upper)
                .map_err(early_archive_error)?,
        ),
        (None, None) => return Err(crate::blastinput::app::missing_subject_error()),
    };
    // With 2+ threads batches are searched on other threads while this function runs. A diagnostic
    // or a panic message written there must not wait for a lock this thread keeps for the whole
    // run, so the handle is used unlocked: the same bytes, the lock taken for each write.
    // No NCBI counterpart: only the stderr lock is taken per write instead of once for the run, so
    // that other threads can write diagnostics (the bytes and their order are unchanged); it does
    // not change any value NCBI computes.
    // With the BLAST archive the diagnostics also go to its message log (`DiagTee`): the
    // messages of each archive are the diagnostics posted since the previous one.
    let mut archive_stderr = crate::cli::CliStderr::process();
    let mut archive_tee;
    let mut stderr: Box<dyn Write + '_> = match &archive {
        Some(handle) => {
            archive_tee = crate::report::archive::DiagTee::new(&mut archive_stderr, handle.clone());
            Box::new(&mut archive_tee)
        }
        None => Box::new(crate::cli::CliStderr::process()),
    };
    // The taxid or seqid list (line DI-8): NCBI's handler of the database arguments
    // (`s_GetTaxIDList`, `CSeqDBFileGiList`) comes before the formatting handler and before
    // `CStdCmdLineArgs`, which opens the query and the output: its warning before theirs (probe
    // F04 of line DI-8), and a list error leaves no `-out`, only the error archive that NCBI's
    // `PrintErrorArchive` writes for outfmt 11 (probe E of line DI-8). The registry's `BLASTDB`
    // finds `-db` (`settings`, read at the start of the run).
    let list = match args.database_name() {
        Some(_) => args
            .db_filter
            .resolve(&settings, "blastx", &mut stderr)
            .map_err(early_archive_error)?,
        None => None,
    };
    // Diagnostic (LOSAT_TIMING): the subject FASTA for the common line (utils/xtiming.rs). No NCBI
    // counterpart.
    let x_t_read_fasta = crate::utils::xtiming::clock();
    let (subject_range, mut subjects) = read_subjects(SubjectRead {
        args,
        settings: &settings,
        id_scope: &id_scope,
        subject_file,
        archive: &archive,
        strategy: &strategy,
        stderr: &mut stderr,
        early_archive_error: &early_archive_error,
    })?;
    // The subjects as NCBI searches them: each cut to its interval (an interval without letters
    // leaves a record without letters, `write_empty_subject_warnings`). The Bioseqs of the archive
    // and of the search strategy and the lower case of outfmt 0's Sbjct rows (`x_soft_subject::x_set_subject_row_masks`, records
    // with lower case only) read the whole record.
    if let Some(range) = &subject_range {
        // Outfmt 0 to 4: the rows of outfmt 0 and the query-anchored sections (1 to 4).
        let shows_lower_case = args.lcase_masking
            && choice
                .number()
                .is_some_and(|number| (0..=4).contains(&number));
        for record in &mut subjects {
            let (from, to_exclusive) =
                crate::blastinput::seq_range::record_cut(range, record.sequence.len())?;
            let keep_whole = archive.is_some()
                || strategy.is_some()
                || (shows_lower_case && !record.lowercase_masks.is_empty());
            record.cut_in_place(from, to_exclusive, keep_whole);
        }
    }
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::ReadFasta, x_t_read_fasta);
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
    if args.db.is_none() && subjects.is_empty() {
        return Err(early_archive_error(engine("Empty CBlastQueryVector")));
    }
    let mut query_file =
        crate::blastinput::input_files::open_input(&args.query, "query", program.upper)
            .map_err(early_archive_error)?;
    // BLAST XML (outfmt 5) and the XML2 and JSON documents (outfmt 13 to 16): the output can
    // hold the whole document in memory (`blast_xml::HeldOutput`), decided below before the
    // search.
    // Outfmt 13 and 14 write side files next to `-out` and the master file at `-out`: with no
    // file name they stop at NCBI's `xinclude_check` below (`Please provide a file name`).
    let out_file = args.out.as_deref().filter(|p| *p != Path::new("-"));
    let structured_xml = match choice.number() {
        Some(5 | 15 | 16) => true,
        Some(13 | 14) => out_file.is_some(),
        _ => false,
    };
    // The hold of the side files of outfmt 13 and 14 follows the master file's.
    let side_files = out_file.and_then(|out| {
        let kind = app::ReportFormat::from_number(choice.number()?).and_then(Xml2Kind::of)?;
        kind.is_files()
            .then(|| super::blast_xml::HeldSideFiles::new(Box::new(FsSideFiles::new(out, kind))))
    });
    // The `-out` file, where NCBI opens it (`cli::create_out_file`, one owner with the library).
    let create = |path: &Path| crate::cli::create_out_file(path, None);
    let mut held = None;
    // The report's stream (`cli::NcbiOstream`). The organism report's: NCBI's buffer of the file or
    // of standard output. Every other format's: LOSAT's buffer before the file or standard output,
    // behind the held document.
    let device = crate::cli::ReportStream::of_cli(
        out_file.map(create).transpose()?,
        program.lower,
        &args.outfmt,
        strategy.as_ref(),
    );
    let stream = match structured_xml {
        true => {
            let report = device.handle();
            crate::cli::NcbiOstream::before_report(
                held.insert(
                    super::blast_xml::HeldOutput::new(Box::new(device), side_files.clone())
                        .with_report(report.clone()),
                )
                .clone(),
                report,
            )
        }
        false => crate::cli::NcbiOstream::report(device),
    };
    // The run's standard error goes through the held output (`HeldOutput::stderr`): a held
    // document gets the texts written in the meantime where NCBI writes them. The structured
    // formats write no BLAST archive, so the run's standard error is the process's, tied to the
    // report on standard output (`cli::CliStderr`).
    if let Some(held) = &held {
        debug_assert!(archive.is_none());
        stderr = Box::new(held.stderr(Box::new(crate::cli::CliStderr::process())));
    }
    // The run keeps its output in the caller's slot: an error of NCBI's own is followed by the
    // error archive in the same output (`run`).
    let writer = output.insert(stream);
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2975-2977
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2874-2886
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:110-126
    // The formatting handler raises NCBI's errors for the formats of other programs (SAM,
    // AIRR rearrangement, FASTA) before its warning. The deferred errors of `resolve` are
    // those of the composition-based statistics handler and of `Validate`, which come later.
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1995-1999
    // The query options handler comes before the genetic code and formatting handlers
    // (blastx_args.cpp:104-115). LOSAT's own refusals of option values come earlier (NCBI's
    // argument reader, `cli::reader`, and `resolve`).
    // The filtering handler comes before the query options handler (blastx_args.cpp:87-107):
    // NCBI's error of `-seg` first (`local::check_filtering`).
    super::local::check_filtering(args)?;
    let query_range = crate::blastinput::seq_range::parse_optional_range(
        args.query_loc.as_deref(),
        crate::blastinput::seq_range::RangeRole::Query,
    )?;
    // A value that the early reading left to this handler (`EarlyFormat::Deferred`) fails here.
    let choice = choice.at_formatting_handler(&args.outfmt)?;
    app::formatting_handler_check(&choice, Program::Blastx)?;
    let html = app::html_links(
        args.html,
        Program::Blastx,
        args.database_name(),
        &[choice.number],
    )?;
    // The warnings of the hit list size and of an ignored -num_descriptions
    // (`blastinput::report_counts`).
    stderr.write_all(&warnings_of(args, choice.number))?;
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:172-175
    let mut options = resolution.map_err(options_error)?;
    options.html = html;
    options.query_range = query_range;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:195-200
    // The database of `-db` opens after the arguments (the query and output files) and the
    // options are processed and before the query is read, as NCBI's `InitializeSubject`.
    let database = match args.database_name() {
        Some(name) => {
            let (keeps_scope, separate_xml2) =
                crate::report::db_subject::FormatThrow::scope_flags(Some(&args.outfmt));
            let database = BlastxDatabase::open_restricted(
                name,
                &settings,
                list.as_ref(),
                args.db_mask.option().as_ref(),
                &mut stderr,
                program.upper,
            )?
            .showing_gis(args.show_gis)
            .after_calls(args.first_query_ordinal > 0)
            .keeping_scope(keeps_scope, separate_xml2);
            // RPSTBLASTN: each batch's `EffectiveSearchSpace` option is its last search unit's
            // own (`BatchResults::request_space`), and its units read `-dbsize` themselves
            // (`rps_search::volumes::RpsUnits`): no program of RPS-BLAST takes a list.
            let rps = options.rpstblastn.is_some();
            crate::algorithm::common::database::note_archive_database_filter(
                archive.as_ref(),
                strategy.as_ref(),
                &args.db_filter,
                database.list_taxids(),
                &database.statistics(),
                false,
                rps || args.searchsp.is_some_and(|space| space > 0),
                !rps && args.dbsize.is_some_and(|size| size > 0),
                program.upper,
            )?;
            crate::algorithm::common::database::check_formats_report_sums(
                || database.pairwise_database(),
                // The command line's one format fails the call.
                [(
                    args.outfmt.as_str(),
                    &mut crate::api::local_blast::FormatRefusal::FailsCall,
                )],
                program.upper,
            )?;
            Some(database)
        }
        None => None,
    };
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:126-135
    options.filtered_empty = database
        .as_ref()
        .is_some_and(|database| database.filtered_empty());
    // A custom tabular field list that looks every taxid up (the plural names of the taxonomy
    // database, the organism report) asks whether `taxdb.bti` is on the search path: the plural
    // names' formatters warn without it. A `-db` search asks its database; a `-subject` search
    // the search path of NCBI's databases (the environment and the registry), as BLASTP's and
    // BLASTN's (`local::subject_taxonomy`).
    match &database {
        Some(database) => {
            if options
                .tabular
                .is_some_and(|spec| spec.looks_up_every_taxid())
            {
                options.taxdb_found = database.taxdb_found();
            }
        }
        None => super::local::subject_taxonomy(
            std::slice::from_mut(&mut options),
            Some(&crate::algorithm::common::database::TaxdbSetting::of(
                &settings,
            )),
        ),
    }
    // LOSAT's refusal of a database name for the Seq-annot formats, before the search
    // (`report::serial::seqalign::check_database`).
    if let (Some(database), Some(serial)) = (
        database.as_ref(),
        app::ReportFormat::from_number(choice.number)
            .and_then(crate::report::serial::SerialFormat::of),
    ) {
        crate::report::serial::seqalign::check_database(serial, database.name(), Program::Blastx)?;
    }
    // How the search reads the query file (`-parse_deflines`, its ids, `-strand`), which the
    // reads of the file before the search share (`QueryStreamSettings`).
    let query_settings = QueryStreamSettings {
        program,
        parse_deflines: args.parse_deflines,
        id_scope: id_scope.clone(),
        strand: crate::report::archive::bioseq::SeqLocStrand::of_option(&args.strand),
        data_loaders: settings.data_loaders,
        batch_count: crate::blastinput::query_batch::BatchCount::of_search(
            args.parse_deflines,
            database.is_none(),
            choice.number == 5,
            program.upper,
        ),
    };
    // The letters of a reader batch (below, where NCBI reads `BATCH_SIZE`, after the check of
    // an empty query), which the reads before the search take too.
    let default_batch_size = match &args.rpstblastn {
        Some(_) => super::query_setup::RPSTBLASTN_BATCH_SIZE,
        None => BATCH_SIZE,
    };
    let batch_size = crate::blastinput::app::query_batch_size(default_batch_size as u32)
        .map(|size| size as usize);
    // LOSAT's other refusals of BLAST XML (a title NCBI crashes on, `names::Refusal`)
    // come before the first byte, though the query file is read batch by batch: the query
    // titles are checked by a read of the file before the search, and a run in which a
    // `-subject` title, the message of a query or a subject of a database with parsed Seq-ids
    // could still be refused later holds the document in memory until its epilog (`blast_xml::HeldOutput`; outfmt 13 and 14: the side files
    // too). Only with outfmt 5 and 13 to 16; the other formats read the file once.
    // The run stops at a batch after this many queries (`blast_xml::QueryScan::stops_after`).
    let mut stops_after = None;
    // Under `-parse_deflines`, the scopes of the batches' queries (`QueryScopes`): with outfmt 5
    // and 13 to 16 the read before the search makes them (`QueryScopes::read_ahead`), else the
    // batches add theirs where the loop takes them (`take`, `add_batch_scope`). Either way LOSAT's
    // refusal of two queries with one id comes where the loop takes the batch of the second, after
    // the reports of the batches before (line WQA D10 round 4, 2026-11-20), before the batch's
    // request and search (round 5, `add_batch_scope`).
    let query_scopes_of = || {
        crate::blastinput::query_batch::QueryScopes::of_search(
            args.parse_deflines,
            database.is_none(),
            choice.number == 5,
            program.upper,
        )
    };
    let mut read_ahead_scopes = None;
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xmlformat.cpp:213-218 (cited at
    // `blast_xml::write_blastx_xml`)
    // Where the run's XML format fails for a matrix name that NCBI's BLAST XML report data does not
    // know (`blast_xml::xml_matrix_failure`, decided once for the run): before any report is
    // written, so NCBI writes no title. Such a run reads the queries before the search only for
    // their scopes (`-parse_deflines`), checks no title and holds nothing.
    let xml_matrix =
        super::blast_xml::xml_matrix_failure(std::slice::from_ref(&options), database.is_none());
    let reads_ahead = held
        .as_ref()
        .filter(|_| xml_matrix.is_none() || query_scopes_of().is_some());
    if let Some(held) = reads_ahead {
        // A query that is not a regular file (a pipe, `<(...)`, `/dev/stdin` from a pipe, a
        // FIFO) or `-query -` cannot be read a second time: the read before the search would
        // take every byte (DECISIONS 22). Refused before the search.
        if args.query.as_os_str() == "-"
            || !query_file.metadata().is_ok_and(|m| m.is_file())
            || !crate::blastinput::input_files::has_stream_position(&query_file)
        {
            // RPSTBLASTN reads such a query into an unnamed temporary file first and reads it
            // from there (twice); BLASTX keeps DECISIONS 22's refusal. The copy is read as
            // the stream is read (`InputFile::into_seekable_copy`).
            let copy = args
                .rpstblastn
                .is_some()
                .then(|| query_file.into_seekable_copy())
                .flatten();
            query_file = match copy {
                Some(copy) => copy,
                None => bail!(
                    "unsupported {} query input for -outfmt {}: use a regular file",
                    program.upper,
                    choice.number
                ),
            };
        }
        // The read before the search reads the open file itself (a `.gz` query is an unnamed
        // inflated file, `input_files::open_file_input`, which no path reaches) and gives it back
        // at the position it had (`prescan_queries`): no second handle, which the WASI targets
        // cannot make.
        // A `BATCH_SIZE` that the run refuses (below, before any query is read) leaves no title
        // to check, as a batch size of 0 does.
        // A `-subject` XML2 or JSON document holds each report back until the next query's
        // arrives (`blast_xml2::reports_held`).
        let holds_back = super::blast_xml::xml2_kind(options.outfmt).is_some()
            && crate::report::blast_xml2::reports_held(database.is_none(), false);
        let (scan, file) = super::blast_xml::prescan_queries(
            query_file,
            &args.query,
            &query_settings,
            &options,
            args.lcase_masking,
            *batch_size.as_ref().unwrap_or(&0),
            query_range.as_ref(),
            holds_back,
            xml_matrix.is_none(),
            query_scopes_of().map(|scopes| scopes.with_xml_matrix_failure(xml_matrix)),
        )?;
        query_file = file;
        stops_after = scan.stops_after;
        read_ahead_scopes = scan.query_scopes;
        // A `-db` search in which a later batch can still refuse the taxonomy database or a
        // subject with parsed Seq-ids (outfmt 5 looks the batch's subjects up in its prefetch,
        // `ReportedSubjects::batch_prefetch`, DBID item T, and asks no scientific name), or whose
        // database counts the reports would refuse at the first report (`write_blastx_xml`,
        // after NCBI's own errors of the search, which stop the run after the document's start,
        // written with the prolog below), holds its document too, as does a `-subject` title
        // (`blast_xml::documents_may_refuse`, the library's decision too).
        // A run whose XML format fails for its matrix name holds nothing (`xml_matrix`: no title is
        // checked, the failure comes before any report and any start).
        if xml_matrix.is_none()
            && super::blast_xml::documents_may_refuse(
                &[choice.number],
                database.as_ref(),
                &subjects,
            )
        {
            held.hold();
        }
    }
    // The query whose batch the failure stops, of the scopes read ahead
    // (`QueryScopes::xml_failing_query`).
    let xml_failing = xml_matrix.map(|failure| {
        crate::blastinput::query_batch::QueryScopes::xml_failing_query(
            read_ahead_scopes.as_ref(),
            failure,
        )
    });
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-873
    let mut query_stream = query_settings.stream(query_file);
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:146-155 (cited at
    // `input::emit_fasta_batch`: the batch drops a record whose range starts past its end, after
    // the reader added it to its scope)
    // The records that `-query_loc` skips are kept for NCBI's scope of a batch's prefetch, which
    // only a stop of the prefetch asks (`BatchScopes::of_batch`): where a prefetch of the run can
    // stop (review 5 of line DBID loadord, L2; the gate opens the database's lookup files once);
    // and under `-parse_deflines` for the scopes of the queries (`add_batch_scope`,
    // `QueryStreamSettings::stream`).
    query_stream.keep_skipped |= query_range.is_some()
        && database
            .as_ref()
            .is_some_and(|database| database.prefetch_can_stop());
    // With `-subject`, the letters around each `-query_loc` interval that the Seq-loc mapper can
    // read into a query row (`empty_rows::query_letters_around`): NCBI reads the subjects first.
    if args.db.is_none() {
        query_stream.around = super::empty_rows::query_letters_around(&options, &subjects);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:242-245
    // With `-query -` and `-subject -` the subjects were read from `cin` to its end above: the
    // stream has failed, its `tellg()` is negative, and the test above answers "not empty" even
    // for a regular file as the standard input. No query batch is read then: outfmt 7 writes
    // `# BLAST processed 0 queries`, never `Query is Empty!` (blastp's `shared_standard_input`,
    // the same rule).
    let shared_standard_input =
        args.query.as_os_str() == "-" && subject_path.is_some_and(|path| path.as_os_str() == "-");
    if !shared_standard_input && stream_is_empty(&mut query_stream) {
        stderr.write_all(format!("Warning: [{}] Query is Empty!\n", program.lower).as_bytes())?;
        // NCBI returns before the formatter exists: no archive carries the message.
        if let Some(handle) = &archive {
            handle.take_messages();
        }
        return Ok(());
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:217-218
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:85-91
    // The letters of a reader / report batch: `BATCH_SIZE` when it is set, else 10002. The
    // search itself keeps units that end at the first query that brings them to 10002 letters
    // (`BATCH_SIZE` the constant), however large a reader batch is: a reader batch is split into
    // units that are searched as the batches of a run without `BATCH_SIZE` are, and its report is
    // written when its last unit is searched; where a unit cannot give the batch's rows alone it
    // is the whole reader batch (`query_setup::unit_ranges`). (`batch_size` is read above, for the
    // reads before the search.)
    let batch_size = batch_size?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:140-153 (cited at
    // `blastinput::db_mask::DbMaskOption::resolve`)
    // The formatter resolves the database masking (`-db_soft_mask`, `-db_hard_mask`; series NP,
    // M4b) before its other work (`BlastxDatabase::resolve_mask_option`).
    if let (Some(database), Some(option)) = (&database, args.db_mask.option()) {
        options.db_mask = database.resolve_mask_option(
            &option,
            archive.as_ref(),
            strategy.as_ref(),
            &mut *stderr,
        )?;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:773-789
    write_empty_subject_warnings(&subjects, &mut *stderr)?;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:251-254
    // The XInclude formats written to standard output fail after the formatter is made (a
    // format 13 or 14 to a file is LOSAT's rejection, above).
    app::xinclude_check_for(&choice, args.out.as_deref())?;
    // `BL2SEQ_LEGACY` makes NCBI search each `-subject` sequence on its own and write the
    // bl2seq report; LOSAT refuses it before the search, where an empty query has already ended
    // the run as NCBI ends it (`Query is Empty!`, exit status 0).
    // LOSAT's refusals of the masked search (`ResolvedMask::check_losat_limits`) come where NCBI
    // starts the search, after its checks.
    if let (Some(database), Some(mask)) = (&database, options.db_mask.as_ref()) {
        database.check_mask_limits(mask, &options)?;
    }
    app::check_search_environment(program.upper, database.is_none())?;
    // RPSTBLASTN's refusals of the database's files and of its options, before the first byte
    // of the reports, as rpsblast refuses them (`args::RpsTblastnSearch::for_database`, once per
    // run).
    if let Some(search) = options.rpstblastn.clone() {
        let handle = archive
            .as_ref()
            .map(crate::report::archive::ArchiveHandle::options)
            .or_else(|| strategy.as_ref().map(|strategy| strategy.options()));
        let search = search.for_database(&options, database.as_ref(), handle.is_some())?;
        // What each batch's RPS set-up writes into the options handle, which the requests of the
        // archive and of the search strategy hold (`OptionsHandle::pass`): the search for the
        // database knows its volumes' costs and its units' database length.
        if let Some(handle) = handle {
            handle.set_search_sets(Vec::new(), search.handle_sets())?;
        }
        options.rpstblastn = Some(search);
    }
    // LOSAT_X_DBIDXREAD acts on a database's index files only; the new side is marked where a
    // volume is opened (`seqdb::Volume::open`).
    if database.is_none() && crate::seqdb::x_dbidxread() != crate::seqdb::XDbIdxRead::Off {
        crate::utils::xreach::x_reach!("LOSAT_X_DBIDXREAD", false);
    }
    let pairwise_database = database
        .as_ref()
        .map(|database| database.pairwise_database_masked(options.db_mask.as_ref()));
    let no_reported = ReportedSubjects::default();
    let records_parsed = if args.parse_deflines {
        report::RecordsParsed::of(&subjects, args.show_gis, |parsed| {
            if options.pairwise_view().is_some() {
                crate::report::db_subject::check_pairwise_records(parsed, true, "BLASTX")
            } else {
                Ok(())
            }
        })?
    } else {
        report::RecordsParsed::default()
    };
    let shown_run = match (&database, &pairwise_database) {
        (Some(database), Some(pairwise)) => ShownSubjects::Database {
            name: database.name(),
            pairwise,
            reported: &no_reported,
        },
        _ => ShownSubjects::Fasta {
            path: subject_path.expect("a subject file without -db"),
            parsed: &records_parsed,
        },
    };
    super::x_timing::run(XRun::Startup, x_t_run);
    let x_t = super::x_timing::lap();
    // The writer of the XML2 or JSON document (outfmt 15 and 16), which the run keeps across the
    // query batches; its start is written with the prolog, as NCBI writes `PrintProlog`'s.
    let xml2 = std::cell::RefCell::new(super::blast_xml::xml2_kind(options.outfmt).map(|kind| {
        let files: RunSideFiles = match &side_files {
            Some(side_files) => Box::new(side_files.clone()),
            None => Box::new(NoSideFiles),
        };
        super::blast_xml::xml2_writer(kind, &options, &subjects, &shown_run, files)
    }));
    report::write_prolog(
        &mut *writer,
        &mut *stderr,
        None,
        &options,
        &subjects,
        &shown_run,
    )
    .map_err(report_error)?;
    if let Some(xml2) = xml2.borrow_mut().as_mut() {
        let mut start = Vec::new();
        xml2.document_mut().start(&mut start);
        writer.write_all(&start)?;
    }
    // NCBI writes the prolog out before it searches the first batch, so the warnings of the
    // batch come after it when both streams go to one place.
    writer.flush()?;
    super::x_timing::run(XRun::Prolog, x_t);
    let x_t = super::x_timing::lap();
    // The database's sequences, read in blocks where they are searched and checked where they are
    // first read (`x_db_stream.rs`); the empty database's error belongs to the first batch with a
    // valid context (`Subjects::Unreadable`, `search_error`). No NCBI counterpart for the blocks:
    // NCBI maps the volume files.
    let (x_longest, unreadable) = match &database {
        None => (None, None),
        Some(database) => match database.x_open_stream() {
            Ok(longest) => (Some(longest), None),
            Err(error) => (None, Some(error)),
        },
    };
    let x_stream = match (&database, x_longest) {
        (Some(database), Some(longest)) => Some(super::x_db_stream::XStreamDb::new(
            database.x_seqdb(),
            longest,
        )),
        _ => None,
    };
    super::x_timing::run(XRun::ReadStores, x_t);
    let unreadable = std::cell::RefCell::new(unreadable);
    let search_subjects = match (&database, &x_stream) {
        (None, _) => Subjects::Fasta(&subjects),
        (Some(_), Some(stream)) => Subjects::Stream(super::x_db_stream::XStream {
            db: stream,
            block: None,
        }),
        (Some(_), None) => Subjects::Unreadable,
    };
    let search_error = |error: anyhow::Error| -> anyhow::Error {
        database_search_error(error, x_stream.as_ref(), &unreadable)
    };
    // The reports of one query batch: with `-db`, the titles of its reported sequences are
    // read first and its hit lists renamed to ranks (`ReportedSubjects::for_batch`); NCBI's
    // reader error there stops the run at this batch, with no epilog.
    // A refusal found while a batch is reported drops what the run holds
    // (`blast_xml::drops_held_document`: a title, the database counts, a reported sequence with
    // parsed Seq-ids such as RPSTBLASTN's profile with a query's generated id, PS15a-2).
    let refused = |error: anyhow::Error| -> anyhow::Error {
        let error = report_error(error);
        if let (true, Some(held)) = (super::blast_xml::drops_held_document(&error), &held) {
            held.discard();
        }
        error
    };
    let query_scopes = std::cell::RefCell::new(read_ahead_scopes.or_else(query_scopes_of));
    // The requests of the archives (outfmt 11), with every `-subject` record, made once for every
    // batch.
    let archive_requests = archive.as_ref().map(|_| {
        archive::requests(match (&database, &pairwise_database) {
            (Some(_), Some(_)) => None,
            _ => Some(&subjects[..]),
        })
    });
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:261-264
    // The search strategy's request of a reader batch (`-export_search_strategy`), which NCBI
    // writes once its query factory is made, before its search: the loop takes each batch here
    // before it looks at the batch's search (`run_query_batches`), so a batch whose search fails
    // (a query without letters in a batch without a valid context, `Warning: Sequence contains no
    // data`; the sequence file of `-db`) has one, and the file is apart from the reports. An empty
    // batch has none: its `Empty CBlastQueryVector` comes before it. Each unit of a batch the
    // search split takes it once (`StrategyExport::takes`).
    let batch_take = BatchTake::new(
        &query_scopes,
        strategy.as_ref(),
        archive.as_ref(),
        &options,
        search_subjects,
        &search_error,
        database.is_none().then_some(&subjects[..]),
    )
    .with_held(held.as_ref());
    // NCBI's scope of the query records at each batch's prefetch (`BatchScopes`).
    let batch_scopes = BatchScopes::new(options.outfmt == 5, database.as_ref());
    let render_batch = |writer: &mut crate::cli::NcbiOstream<'static>,
                        diagnostics: &mut dyn Write,
                        queries: &[FastaRecord],
                        results: &mut Vec<BatchResults>,
                        queries_before: usize,
                        pool: Option<&crate::utils::threading::SearchPool<'_>>|
     -> Result<()> {
        // The batch's set-up has written its `EffectiveSearchSpace` option, which the next batch's
        // search strategy request holds (`take` below).
        batch_take.note_searched(queries, results, queries_before)?;
        if let (Some(handle), Some(archive_requests)) = (&archive, &archive_requests) {
            return ArchiveBatch {
                handle,
                archive_requests,
                database: &database,
                pairwise_database: &pairwise_database,
                search_subjects,
                options: &options,
                subjects: &subjects,
                query_scopes: &query_scopes,
            }
            .write(writer, diagnostics, queries, results, queries_before);
        }
        // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:283-288
        // The prefetch before the batch's reports reads `PRE_FETCH_SEQS_LIMIT`, after the batch's
        // search and, with `-db`, the read of its reported sequences' titles (NCBI's traceback,
        // `batch_titles`): a value that NCBI cannot convert stops the first batch here, after the
        // outfmt 0 prolog and the messages of its records, before its query warnings
        // (`app::check_prefetch_environment`).
        let drop_held = |error: &anyhow::Error| {
            if let (true, Some(held)) = (super::blast_xml::drops_held_document(error), &held) {
                held.discard();
            }
        };
        let asks = ReportedAsks::of(&options);
        let titles = match &database {
            Some(database) => batch_titles(
                database,
                search_subjects,
                results,
                options.gapped,
                asks,
                &mut *diagnostics,
            )
            .inspect_err(drop_held)?,
            None => None,
        };
        app::check_prefetch_environment(false)?;
        // The run stops after this batch (`stops_after`).
        let stops = stops_after == Some(queries_before + queries.len());
        match (&database, &pairwise_database) {
            (Some(database), Some(pairwise)) => {
                let scope = || batch_scopes.of_batch(queries);
                let reported = reported_batch(
                    database,
                    search_subjects,
                    titles,
                    &options,
                    asks,
                    results,
                    (queries, queries_before),
                    query_scopes
                        .borrow()
                        .as_ref()
                        .map(|scopes| scopes.numbered_from(queries_before)),
                    &crate::algorithm::common::database::PrefetchQueries::new(0, &scope),
                    diagnostics,
                )
                .inspect_err(drop_held)?;
                batch_scopes.reported(queries);
                report::render_queries(
                    writer,
                    diagnostics,
                    &options,
                    queries,
                    reported.records(),
                    &ShownSubjects::Database {
                        name: database.name(),
                        pairwise,
                        reported: &reported,
                    },
                    results,
                    queries_before,
                    stops,
                    xml2.borrow_mut().as_mut(),
                    xml_failing,
                    pool,
                )
                .map_err(refused)
            }
            _ => report::render_queries(
                writer,
                diagnostics,
                &options,
                queries,
                &subjects,
                &shown_run,
                results,
                queries_before,
                stops,
                xml2.borrow_mut().as_mut(),
                xml_failing,
                pool,
            )
            .map_err(refused),
        }
    };
    // NCBI: a batch size of 0 makes the first batch empty, and its query factory refuses it.
    if batch_size == 0 && !query_stream.at_eof() {
        return Err(engine("Empty CBlastQueryVector"));
    }
    let mut read = |deliver: &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
                    warn: &mut dyn FnMut(&[u8]) -> Result<()>|
     -> Result<()> {
        read_fasta_stream_batches(
            &mut query_stream,
            false,
            args.lcase_masking,
            batch_size,
            archive.is_some() || strategy.is_some(),
            query_range.as_ref(),
            deliver,
            &mut |warning| {
                let mut bytes = Vec::new();
                reader_warning(warning, &mut bytes)?;
                warn(&bytes)
            },
        )
    };
    let mut render = |diagnostics: &mut dyn Write,
                      queries: &[FastaRecord],
                      results: &mut Vec<BatchResults>,
                      queries_before: usize,
                      pool: Option<&crate::utils::threading::SearchPool<'_>>|
     -> Result<()> {
        render_batch(
            &mut *writer,
            diagnostics,
            queries,
            results,
            queries_before,
            pool,
        )?;
        writer.flush().map_err(Into::into)
    };
    let mut take = |diagnostics: &mut dyn Write,
                    queries: &[FastaRecord],
                    queries_before: usize|
     -> Result<()> { batch_take.take(diagnostics, queries, queries_before) };
    // RPSTBLASTN's search of the query batches takes the place of BLASTX's, once per run
    // (`args::RpsTblastnSearch`, DECISIONS 38 D1).
    let processed = match &options.rpstblastn {
        Some(search) => search.run_query_batches(
            &options,
            database.as_ref(),
            &mut read,
            &mut *stderr,
            &search_error,
            &mut render,
            &mut take,
        ),
        None => run_query_batches(
            &options,
            search_subjects,
            args.parse_deflines,
            &mut read,
            &mut *stderr,
            &search_error,
            &mut render,
            &mut take,
        ),
    }
    .map_err(input_error)?;
    let x_t = super::x_timing::lap();
    report::write_epilog(
        &mut *writer,
        &mut *stderr,
        &options,
        &subjects,
        &shown_run,
        processed,
        xml2.borrow_mut().as_mut(),
        xml_failing,
    )
    .map_err(report_error)?;
    writer.flush()?;
    // The document a hold kept: written at the end of the run, with the epilog.
    if let Some(held) = &held {
        held.release().map_err(report_error)?;
    }
    super::x_timing::run(XRun::Epilog, x_t);
    super::x_timing::run(XRun::Total, x_t_run);
    // Diagnostic (LOSAT_TIMING): the run's share of the common line (x_timing.rs).
    super::x_timing::common_run();
    super::x_timing::write(&mut *stderr, options.num_threads as usize)?;
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:286-291
/// What the reports of a query batch read of its database sequences
/// (`ReportedSubjects::for_batch`): the arguments of one format (`of`), or of several (`or`).
#[derive(Clone, Copy, Default)]
pub(super) struct ReportedAsks {
    tabular: bool,
    tabular_reads: bool,
    pairwise: bool,
    structured: bool,
    taxids: bool,
    joined: bool,
    taxonomy: bool,
    words: crate::report::tabular::PrintedIds,
    prefetched: bool,
    /// Whether the batches' prefetch runs where it is on (`prefetched`, or the XML2 and JSON
    /// formats 13 to 16, which write the taxonomy texts at each `Hit`): a subject whose lookup
    /// fails on a file stops the run in it (`shared::batch_prefetch`, `batch_prefetch_stop`).
    prefetch_runs: bool,
    /// Whether the reports may read the reported sequences one by one
    /// (`common::database::prefetch_may_be_off` of the first format), for the reader of the
    /// headers (`HeaderReader::new`).
    one_by_one: bool,
    /// What the reports do with a subject whose lookup throws where a batch's prefetch is off
    /// (`FormatThrow::and` of the reports, for `HitLoads::of`).
    throw: crate::report::db_subject::FormatThrow,
    /// Whether a report writes otherwise where every load of a subject with a positive taxid
    /// throws (not `report::db_subject::unchanged_by_load_throws`), for the batches' prefetch.
    changed_by_load_throws: bool,
}

impl ReportedAsks {
    /// What the reports of `options.outfmt` ask.
    pub(super) fn of(options: &super::args::ResolvedOptions) -> Self {
        // The organism report (18) makes no tabular subject ids (`report::tabular::writes_rows`).
        let tabular = matches!(options.outfmt, 6 | 7 | 10 | 18 | 20)
            && crate::report::tabular::writes_rows(options.tabular);
        Self {
            tabular,
            tabular_reads: tabular
                && crate::report::tabular::reads_subject(
                    options.tabular,
                    crate::blastinput::app::Program::Blastx,
                ),
            pairwise: options.pairwise_view().is_some(),
            structured: matches!(options.outfmt, 5 | 13..=16),
            taxids: matches!(options.outfmt, 13..=16),
            joined: matches!(options.outfmt, 0..=5 | 13..=16)
                || options
                    .tabular
                    .is_some_and(|spec| spec.needs().joins_subjects()),
            taxonomy: matches!(options.outfmt, 13..=16)
                || options.tabular.is_some_and(|spec| spec.needs().taxonomy),
            words: if tabular {
                crate::report::tabular::PrintedIds::of(options.tabular)
            } else {
                crate::report::tabular::PrintedIds::default()
            },
            // Outfmt 10 and 20 prefetch as 6 does (`ResolvedOptions::prefetch_format`).
            prefetched: prefetches_subjects(Some(options.prefetch_format()))
                && !crate::algorithm::common::database::prefetch_always_off(
                    options.prefetch_format(),
                ),
            prefetch_runs: (prefetches_subjects(Some(options.prefetch_format()))
                || matches!(options.outfmt, 13..=16))
                && !crate::algorithm::common::database::prefetch_always_off(
                    options.prefetch_format(),
                ),
            one_by_one: crate::algorithm::common::database::prefetch_may_be_off(
                options.prefetch_format(),
            ),
            throw: crate::blastinput::app::ReportFormat::from_number(i32::from(options.outfmt))
                .map_or(
                    crate::report::db_subject::FormatThrow::Refused(
                        crate::report::db_subject::ThrowRefusal::Format,
                    ),
                    |format| {
                        // With `-html` where the run has its links (the writer's
                        // `PairwiseConfig::html`).
                        crate::report::db_subject::FormatThrow::of(
                            format,
                            options.tabular,
                            crate::blastinput::app::Program::Blastx,
                            options.html.is_some(),
                        )
                    },
                ),
            changed_by_load_throws: !crate::blastinput::app::ReportFormat::from_number(i32::from(
                options.outfmt,
            ))
            .is_some_and(|format| {
                crate::report::db_subject::unchanged_by_load_throws(
                    format,
                    options.tabular,
                    crate::blastinput::app::Program::Blastx,
                )
            }),
        }
    }

    /// What the reports `options` of one call ask: each of them, with the batches' prefetch of
    /// the first, which NCBI's app runs for its one format
    /// (`BlastFormatter_PreFetchSequenceData` with `GetFormattedOutputChoice()`), and the call's
    /// refusal of a subject whose lookup throws over the `-outfmt` values of its other searches,
    /// `others`, too (`FormatThrow::of_search`, `ReportOutputs::other_outfmts`).
    pub(super) fn of_all(options: &[super::args::ResolvedOptions], others: &[&str]) -> Self {
        // The first search's asks start the call's: `FormatThrow::and` of a format that reads
        // nothing (`Self::default()`'s) with one that stops is a refusal (`FormatThrow::of_formats`).
        let mut asks = options
            .iter()
            .map(Self::of)
            .reduce(Self::or)
            .unwrap_or_default();
        asks.throw = asks.throw.of_search(
            others,
            crate::blastinput::app::Program::Blastx,
            crate::report::archive::options::RequestKind::Blastx,
            options.iter().any(|options| options.html.is_some()),
        );
        asks.prefetched = options
            .first()
            .is_some_and(|options| Self::of(options).prefetched);
        asks.prefetch_runs = options
            .first()
            .is_some_and(|options| Self::of(options).prefetch_runs);
        asks.one_by_one = options
            .first()
            .is_some_and(|options| Self::of(options).one_by_one);
        asks
    }

    /// What two reports of one batch ask.
    pub(super) fn or(self, other: Self) -> Self {
        Self {
            tabular: self.tabular || other.tabular,
            tabular_reads: self.tabular_reads || other.tabular_reads,
            pairwise: self.pairwise || other.pairwise,
            structured: self.structured || other.structured,
            taxids: self.taxids || other.taxids,
            joined: self.joined || other.joined,
            taxonomy: self.taxonomy || other.taxonomy,
            words: self.words.or(other.words),
            prefetched: self.prefetched || other.prefetched,
            prefetch_runs: self.prefetch_runs || other.prefetch_runs,
            one_by_one: self.one_by_one || other.one_by_one,
            throw: self.throw.and(other.throw),
            changed_by_load_throws: self.changed_by_load_throws || other.changed_by_load_throws,
        }
    }
}

/// The read error of a query batch's reported sequences as the run reports it: NCBI's reader
/// error of a header file (`HeaderReadError`) or of a record's residues (`XReadError`,
/// `database_error`), or LOSAT's refusal as it is.
fn reported_error(error: anyhow::Error) -> anyhow::Error {
    match error.downcast::<HeaderReadError>() {
        Ok(header) => header.error,
        // A failed read of a record's residues.
        Err(error) => match error.downcast::<super::x_db_stream::XReadError>() {
            Ok(super::x_db_stream::XReadError(kept)) => database_error(kept),
            Err(error) => error,
        },
    }
}

/// The titles of the reported sequences of the query batch `results` of a `-db` search
/// (`ReportedSubjects::read_titles`, for the reports that `asks` names), read where
/// NCBI's traceback makes the batch's Seq-aligns, before the batch's prefetch, after the hit lists
/// of no Seq-align leave `results` (`ReportedSubjects::drop_unaligned`): NCBI's reader error
/// stops the run at this batch, with no epilog. `None` where the database's sequences could not be
/// read (the batch's search failed with that error).
/// NCBI's warning of its first read of a UTF8String header (`BatchTitles::utf8`, or that of the
/// failed read) is written to `diagnostics` there: after the batch's reader messages and request
/// posts, before its prefetch.
pub(super) fn batch_titles(
    database: &BlastxDatabase,
    search_subjects: Subjects<'_>,
    results: &mut [BatchResults],
    gapped: bool,
    asks: ReportedAsks,
    diagnostics: &mut dyn Write,
) -> Result<Option<super::database::BatchTitles>> {
    use super::x_timing::Run as XRun;
    if matches!(search_subjects, Subjects::Unreadable) {
        return Ok(None);
    }
    use crate::algorithm::common::database::{EarlierUtf8, Utf8Read};
    let first_query = results.first().map_or(0, |results| results.query_ordinal);
    let fetched = results
        .iter()
        .filter_map(|results| results.fetched_utf8)
        .min()
        .and_then(|ordinal| database.db.location_of(ordinal));
    let aligned = ReportedSubjects::drop_unaligned(database, results, gapped);
    let earlier = EarlierUtf8::of_batch(first_query, fetched, aligned);
    let x_t = super::x_timing::lap();
    let titles = ReportedSubjects::read_titles(
        database,
        results,
        gapped,
        asks.pairwise,
        (earlier, asks.one_by_one),
    );
    super::x_timing::run(XRun::Titles, x_t);
    if let Some(posts) = Utf8Read::post_of(
        Utf8Read::of_read(&titles, |titles| titles.utf8.as_ref()),
        &database.db,
        &database.program.to_ascii_lowercase(),
    ) {
        posts.write(diagnostics, None)?;
    }
    titles.map(Some).map_err(reported_error)
}

/// What NCBI's prefetch (`BlastFormatter_PreFetchSequenceData`, cited at `reported_batch`) does
/// with a query batch.
pub(super) struct BatchPrefetch {
    /// The batch's Seq-aligns, which `PRE_FETCH_SEQS_LIMIT` counts.
    seq_aligns: usize,
    /// The format, whose prefetch the batch asks (`ResolvedOptions::prefetch_format`).
    format: crate::algorithm::common::database::PrefetchFormat,
    /// `PRE_FETCH_SEQS_LIMIT` turns it off (`common::database::prefetch_off`).
    pub(super) off: bool,
    /// It throws at a subject row whose range it cannot take before it loads anything and writes
    /// its warning (`empty_rows::prefetch_reads_unranged_row`).
    pub(super) fails: bool,
}

/// What NCBI's prefetch does with the query batch `results` (`BatchPrefetch`) for the format of
/// `options`, of a report of `loc`'s subjects (`-db`'s, `reported_batch`, or `-subject`'s,
/// `report::render_queries_with`). The batch's Seq-aligns: one for each HSP that the reports show
/// (`has_seqalign`), one for each subject of an ungapped search; where `PRE_FETCH_SEQS_LIMIT`
/// turns the batch's prefetch off (`shared::prefetch_off`), it reads nothing of a subject that no
/// report reads (`asks.prefetched`). Where it reads the Seq-aligns and one has a subject row whose
/// range it cannot take (an empty one: a subject whose every HSP has length 0, as the mapper leaves
/// it; one that starts with a gap: RPSTBLASTN's Smith-Waterman redo), it throws before it loads
/// anything and writes its warning: the batch's subjects are then read as with the prefetch off.
pub(super) fn batch_prefetch(
    results: &[BatchResults],
    options: &super::args::ResolvedOptions,
    loc: super::empty_rows::SubjectLoc,
) -> BatchPrefetch {
    // (RPSTBLASTN: with the HSPs of the profiles that left the merged lists past the formatter's
    // cut, `BatchResults::unlisted_hsps`; its search is gapped.)
    let seq_aligns: usize = results
        .iter()
        .flat_map(|batch| batch.queries.iter().flatten())
        .map(|list| {
            crate::algorithm::common::database::subject_seq_aligns(
                list.hsps
                    .iter()
                    .filter(|hsp| super::report::has_seqalign(hsp, options.gapped))
                    .count(),
                options.gapped,
            )
        })
        .sum::<usize>()
        + results
            .iter()
            .map(|batch| batch.unlisted.iter().sum::<usize>())
            .sum::<usize>();
    let format = options.prefetch_format();
    let off = crate::algorithm::common::database::prefetch_off(format, seq_aligns);
    // A prefetch that runs reads the Seq-aligns (`PrefetchFormat::batch_prefetches`).
    let fails =
        !off && super::empty_rows::prefetch_reads_unranged_row(results, options.gapped, loc);
    BatchPrefetch {
        seq_aligns,
        format,
        off,
        fails,
    }
}

impl BatchPrefetch {
    /// What NCBI's BLAST database data loader loaded for the batch, the run's record
    /// (`ReportedSubjects::note_loaded`): `common::database::batch_loads` of its Seq-aligns (every
    /// subject of its Seq-aligns where its prefetch is on, else what the first format's rows read,
    /// `rows`), and those rows where its prefetch failed before it loaded anything (`fails`, a
    /// subject row whose range it cannot take: BLASTX's ungapped HSPs of length 0, RPSTBLASTN's
    /// alignments that start with a gap; the reports then read the subjects as with the prefetch
    /// off).
    fn loads(
        &self,
        rows: crate::report::tabular::SubjectFetch,
    ) -> crate::algorithm::common::database::BatchLoads {
        if self.fails {
            crate::algorithm::common::database::BatchLoads::Read(rows)
        } else {
            crate::algorithm::common::database::batch_loads(self.format, self.seq_aligns, rows)
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2310-2313 (cited at
// `common::database::stop_lookup_unsure`: the formatter keeps its scope after a batch in outfmt
// 5 only)
/// NCBI's scope of the query records at the prefetch of each query batch of a BLASTX run
/// (`common::database::PrefetchQueries`), whose reader drops each batch once it is reported: the
/// batch's records (those that the `-query_loc` range skipped included, `input::shown_ids`;
/// `query_batch::QueryScope::of_shown`), and in outfmt 5, whose scope the formatter keeps, those
/// of the batches before, which the run keeps as they are reported (`reported`), against a
/// database whose prefetch can stop only (`BlastxDatabase::prefetch_can_stop`, asked at the first
/// batch reported; review 3 of line DBID loadord, M1). `query_batch::QueryScopes`, the owner of the
/// run's scopes, holds the ids that `-parse_deflines` reads only (it refuses two records of one
/// id), where `refuse_id_of_a_query` refuses every reported sequence with an id of the run's
/// scope before the prefetch; the scope kept here is that of the generated ids too, so it stays
/// apart.
pub(super) struct BatchScopes<'a> {
    /// The run's database where its reports include outfmt 5.
    database: Option<&'a BlastxDatabase>,
    before: std::cell::RefCell<Option<crate::blastinput::query_batch::QueryScope>>,
    /// Whether the database's prefetch can stop, asked once.
    keeps: std::cell::OnceCell<bool>,
}

impl<'a> BatchScopes<'a> {
    /// The scopes of a run whose reports include outfmt 5 (`xml`) against `database` (none with
    /// `-subject`).
    pub(super) fn new(xml: bool, database: Option<&'a BlastxDatabase>) -> Self {
        Self {
            database: database.filter(|database| xml && database.parsed_seqids()),
            before: Default::default(),
            keeps: Default::default(),
        }
    }

    /// The scope at the prefetch of the batch of the records `queries`.
    pub(super) fn of_batch(
        &self,
        queries: &[FastaRecord],
    ) -> crate::blastinput::query_batch::QueryScope {
        let mut scope =
            crate::blastinput::query_batch::QueryScope::of_shown(super::input::shown_ids(queries));
        if let Some(before) = self.before.borrow().as_ref() {
            scope.extend(before);
        }
        scope
    }

    /// The batch of the records `queries` was reported: in outfmt 5 its records stay in the scope
    /// of the batches after it, where a prefetch of the run can stop.
    pub(super) fn reported(&self, queries: &[FastaRecord]) {
        let Some(database) = self.database else {
            return;
        };
        if !*self.keeps.get_or_init(|| database.prefetch_can_stop()) {
            return;
        }
        let batch =
            crate::blastinput::query_batch::QueryScope::of_shown(super::input::shown_ids(queries));
        self.before
            .borrow_mut()
            .get_or_insert_with(Default::default)
            .extend(&batch);
    }
}

/// The reported sequences of the query batch `results` of a `-db` search, with its hit lists
/// renamed to ranks (`ReportedSubjects::for_batch`, `asks`), from their titles (`titles`,
/// `batch_titles`; `None`: no sequence could be read), and what NCBI does before the reports of a
/// batch, for the format of `options`: its prefetch and the values of a custom tabular field
/// list. NCBI's reader error stops the run at this batch, with no epilog. A run that holds its
/// document drops it on a refusal from here (`blast_xml::drops_held_document`: of the XML2
/// document, of a subject with parsed Seq-ids, or of the taxonomy database). The stop of the
/// batch's prefetch asks `prefetch_queries` (its scope, `BatchScopes`; its error names query 0:
/// each batch's reports are written before the next is searched, `reported_error`); a batch whose
/// prefetch does not stop is recorded in the run (`ReportedSubjects::note_loaded`), without the
/// loads of an outfmt 5 report that NCBI's XML stream cuts before its hits (`queries`, the reader
/// batch's records: `XmlCuts`), as the loads of the batch numbered `queries_before` (the queries
/// of the run before it), in whose scope its reports read (`ReportedSubjects::for_batch`).
pub(super) fn reported_batch(
    database: &BlastxDatabase,
    search_subjects: Subjects<'_>,
    titles: Option<super::database::BatchTitles>,
    options: &super::args::ResolvedOptions,
    asks: ReportedAsks,
    results: &mut [BatchResults],
    (queries, queries_before): (&[FastaRecord], usize),
    query_scopes: Option<crate::blastinput::query_batch::QueryScopesFrom<'_>>,
    prefetch_queries: &crate::algorithm::common::database::PrefetchQueries<'_>,
    diagnostics: &mut dyn Write,
) -> Result<ReportedSubjects> {
    use super::x_timing::Run as XRun;
    let Some(titles) = titles else {
        return Ok(ReportedSubjects::default());
    };
    let prefetch = batch_prefetch(results, options, super::empty_rows::SubjectLoc::Whole);
    let prefetch_fails = prefetch.fails;
    let unfetched = prefetch.off || prefetch_fails;
    let x_t = super::x_timing::lap();
    let reported = ReportedSubjects::for_batch(
        database,
        search_subjects,
        results,
        titles,
        options.gapped,
        asks.tabular,
        asks.tabular_reads,
        asks.pairwise,
        asks.structured,
        asks.taxids,
        asks.joined,
        asks.taxonomy,
        asks.words,
        asks.prefetched && !unfetched,
        asks.prefetch_runs && !unfetched,
        unfetched,
        options.shown_lists(),
        (
            asks.throw,
            // The pairwise report reads the subjects that it shows as its writer reads them.
            options.pairwise_reading(),
            queries_before,
        ),
    );
    super::x_timing::run(XRun::Titles, x_t);
    let mut reported = reported.map_err(reported_error)?;
    // `-db_soft_mask` and `-db_hard_mask`: the rows of outfmt 0 to 4 show the stored masks in
    // lower case (`ReportedSubjects::shown_db_masks`).
    if let Some(mask) = options.db_mask.as_ref().filter(|_| asks.pairwise) {
        reported.read_db_masks(database, mask.id)?;
    }
    if let Some(scopes) = query_scopes {
        reported.refuse_ids_of_queries(results, options.gapped, scopes)?;
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:286-291
    // Before the reports of a batch NCBI fetches every subject that the batch
    // reports, and with `taxdb.bti` on the search path the taxonomy database
    // writes its messages there (`ReportedSubjects::batch_prefetch`: the
    // texts of opening it at the run's first lookup, a volume's lookup of a
    // taxid once per run): after the reader's messages of the reader batch,
    // which every batch loop writes before it calls this, and before the
    // warnings of its first query (`report::render_queries`). Outfmt 0, 5, 6
    // and 7, and with `PRE_FETCH_SEQS_LIMIT=2147483647` 8, 9 and 12
    // (`database::prefetches_subjects`); XML2 and JSON write theirs at
    // each `Hit` (`TaxdbPrinter`). Without `taxdb.bti` nothing is read. A
    // refusal comes after the reports of the batches before, as the refusals
    // of the batch's titles above (line DBID, DECISIONS 9 (9)); a refusal of
    // the taxonomy database (`TaxonomyRefusal`) drops what an outfmt 5 run
    // holds, and the run holds where one could come
    // (`BlastxDatabase::later_batch_may_refuse`), as XML2's does.
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:716-725
    // Outfmt 10 and 20 read the subjects as 6 and 7 do; the rows' own
    // lookups of the taxonomy names come with the rows (`tabular`).
    let prefetch_format = options.prefetch_format();
    let mut reported = reported;
    if prefetch_fails {
        diagnostics.write_all(&database.prefetch_warning())?;
    }
    // Where the prefetch fails at a subject with parsed Seq-ids, NCBI writes its warning after
    // the taxonomy texts (`ReportedSubjects::batch_prefetch`), with or without `taxdb.bti`.
    let taxdb_found = database.taxdb_found();
    if prefetches_subjects(Some(prefetch_format)) && (taxdb_found || database.parsed_seqids()) {
        // Where `PRE_FETCH_SEQS_LIMIT` turns the batch's prefetch off
        // (`prefetch_off`, above), each query's report fetches the
        // subjects it shows, after the query's own warnings
        // (`ReportedSubjects::fetch_by_query`), where the format's rows read
        // their Bioseqs (`TabularSpec::subject_fetch`).
        // No prefetch; without `taxdb.bti` the reports' fetches look nothing up nor, but where a
        // lookup reads a UTF8String header (`load_read_texts`), post a warning.
        if unfetched && (taxdb_found || database.db.scope_utf8_seen()) {
            let rows = options
                .tabular
                .map_or(crate::report::tabular::SubjectFetch::Every, |spec| {
                    spec.subject_fetch(crate::blastinput::app::Program::Blastx)
                });
            // The writers of outfmt 0 (description lines and alignment headings), 5 (`Hit`s), the
            // tabular rows and the organism report write the texts where they read each subject
            // (`ReportFormat::places_subject_reads`).
            use crate::blastinput::app::ReportFormat;
            let first_report = ReportFormat::from_number(i32::from(options.outfmt));
            let places = first_report.is_some_and(ReportFormat::places_subject_reads);
            // A field list with a plural name and the organism report look every taxid up, taxid
            // 0 too: the first row or report line opens the taxonomy database after the first
            // read (`common::database::subject_read_texts`).
            let rows_open = options
                .tabular
                .is_some_and(|spec| spec.looks_up_every_taxid());
            // The pairwise report (outfmt 0 to 4) reads the subjects that it shows as its writer
            // reads them.
            let pairwise = options.pairwise_reading();
            reported.fetch_by_query(
                database,
                results,
                options.gapped,
                (rows, prefetch_format),
                places,
                rows_open,
                pairwise,
            )?;
            if options.outfmt == 5 {
                diagnostics.write_all(&reported.take_query_fetches())?;
            }
        } else if !unfetched {
            // A lookup that fails on a file stops the run here, after the reports of the
            // batches before and this batch's reader messages (`shared::batch_prefetch`).
            let messages = reported
                .batch_prefetch(
                    database,
                    results,
                    options.gapped,
                    &prefetch_queries
                        .where_reports_unchanged_by_load_throws(!asks.changed_by_load_throws),
                )
                .map_err(reported_error)?;
            messages.write(diagnostics, None)?;
        }
    } else if asks.prefetch_runs && !unfetched && database.parsed_seqids() {
        // The XML2 and JSON formats write the taxonomy texts at each `Hit`: of the batch's
        // prefetch only its stop, or the warning of a UTF8String header that a load's lookup reads
        // (`ReportedSubjects::batch_prefetch_stop`).
        let warning = reported
            .batch_prefetch_stop(database, results, options.gapped, prefetch_queries)
            .map_err(reported_error)?;
        diagnostics.write_all(&warning)?;
    } else if matches!(options.outfmt, 13..=16) && database.parsed_seqids() {
        // The XML2 and JSON formats, a batch whose prefetch is off: the `Hit`s load their subjects
        // one by one.
        reported.scope_utf8_hit_reads(
            database,
            results,
            options.gapped,
            crate::report::db_subject::FormatThrow::of_run([asks.throw]),
        );
    }
    // The batch's prefetch (or its reports' reads one by one) did not stop: the run's record of
    // what NCBI's data loader loaded for it, for a later batch's stop
    // (`ReportedSubjects::note_loaded`, `BatchPrefetch::loads`): the subjects of its Seq-aligns
    // where its prefetch read them, else those that the first format's rows read, of the reports
    // that NCBI writes (outfmt 5: not one that its XML stream cuts before its hits, which loads
    // nothing; `XmlCuts`, review 6 of line DBID loadord, H2).
    let loads = prefetch.loads(
        options
            .tabular
            .map_or(crate::report::tabular::SubjectFetch::Every, |spec| {
                spec.subject_fetch(crate::blastinput::app::Program::Blastx)
            }),
    );
    let xml_cut = |query: usize| {
        queries.get(query).is_some_and(|record| {
            crate::report::blast_xml::names::iteration_cut_before_hits(
                record.internal_id.as_bytes(),
                record.warning_id(),
                record.believed(),
                &record.title,
                false,
            )
        })
    };
    reported.note_loaded(
        (database, queries_before),
        results,
        options.gapped,
        loads,
        crate::algorithm::common::database::XmlCuts::of(options.outfmt == 5, &xml_cut),
    );
    // The values of a custom tabular field list that reads more of the
    // definition lines than the default report (`sall*`) or the taxonomy
    // (`staxid`, `staxids`, the names), read before the batch's first byte so
    // that their refusals come before it, after the prefetch (which opens the
    // taxonomy database first, as NCBI's does).
    if let Some(spec) = options.tabular {
        reported.read_tabular(database, spec, results, options.gapped)?;
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:259
    // NCBI reference (598d8ae6): c++/src/app/blast/rpstblastn_app.cpp:191
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2295-2314
    database.reset_scope();
    Ok(reported)
}

/// The error of a batch whose search fails after its reads of header records (`HeaderReadError`
/// with `utf8`: BLASTX's Smith-Waterman redo that fails the batch after fetches that read a
/// UTF8String where a VisibleString is expected, `kappa_batch::sw_failed_after`): NCBI's warning
/// of that read is written to `diagnostics`, where the run has not posted it (`Utf8Read::post`),
/// after the batch's reader messages and request posts, before the error, which is returned.
/// Any other error as it is.
fn failed_batch_reads(
    error: anyhow::Error,
    options: &super::args::ResolvedOptions,
    subjects: Subjects<'_>,
    diagnostics: &mut dyn Write,
) -> Result<anyhow::Error> {
    let header = match error.downcast::<HeaderReadError>() {
        Ok(header) if header.utf8.is_some() => header,
        Ok(header) => return Ok(header.into()),
        Err(error) => return Ok(error),
    };
    if let (Some(read), Some(db)) = (&header.utf8, subjects.seqdb()) {
        let program = if options.rpstblastn.is_some() {
            "rpstblastn"
        } else {
            "blastx"
        };
        if let Some(posts) = read.post(db, program) {
            posts.write(diagnostics, None)?;
        }
    }
    Ok(header.error)
}

/// Whether an error of a search ends the run as it is: LOSAT's refusal, not an error of NCBI's
/// engine, an error that has NCBI's report already (the exception of a query split size that
/// NCBI cannot convert, `split::split_sizes`), or an error of a query batch's search after the
/// reports of the batches before (`app::ErrorAfterBatchReports`: the failed or endless
/// Smith-Waterman redo of `-use_sw_tback`, `kappa`), which `search_error_held` unwraps. The one
/// list of both callers (`database_search_error`, `local::search_records`).
pub(super) fn search_error_as_is(error: &anyhow::Error) -> bool {
    error.is::<super::parameters::SearchSpaceOverflow>()
        || error.is::<super::split::QuerySplitRefusal>()
        || error.is::<crate::algorithm::common::database::FetchRefusal>()
        || error.is::<super::database::LongSubjectRefusal>()
        || error.is::<crate::algorithm::common::hard_chunk_scan::ShortChunkSeedRefusal>()
        || error.is::<NativeError>()
        || error.is::<crate::blastinput::app::ErrorAfterBatchReports>()
}

/// An error that ends the run as it is (`search_error_as_is`) as the run reports it: the error
/// that an `app::ErrorAfterBatchReports` holds (NCBI's engine error `Traceback failed`, exit
/// status 3, or LOSAT's refusal), any other as it is. The flow writes the reports of the
/// batches before when it meets the batch's error, so the marker has done its work here.
pub(super) fn search_error_held(error: anyhow::Error) -> anyhow::Error {
    match error.downcast::<crate::blastinput::app::ErrorAfterBatchReports>() {
        Ok(crate::blastinput::app::ErrorAfterBatchReports(held)) => held,
        Err(error) => error,
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:226-227
/// The error of the search of a query batch as a `-db` run reports it, and a `-subject` run
/// (no `x_stream`, nothing `unreadable`).
pub(super) fn database_search_error(
    error: anyhow::Error,
    x_stream: Option<&super::x_db_stream::XStreamDb<'_>>,
    unreadable: &std::cell::RefCell<Option<anyhow::Error>>,
) -> anyhow::Error {
    // An error the search raised after the point where a database refused before the search
    // fails it: the database's first error comes before it (`search::search_core`).
    let error = match error.downcast::<super::x_db_stream::XPastHead>() {
        Ok(super::x_db_stream::XPastHead(error)) => {
            if let Some(stream) = x_stream.filter(|stream| !stream.x_checked()) {
                if let Some(first) = stream.x_first_error() {
                    return database_error(first);
                }
            }
            error
        }
        Err(error) => error,
    };
    if search_error_as_is(&error) {
        return search_error_held(error);
    }
    if error.is::<UnreadableSubjects>() {
        let kept = unreadable
            .borrow_mut()
            .take()
            .expect("the database's error, raised once");
        return database_error(kept);
    }
    // A read or a check of the database that failed during the search, reported as the
    // check's error while the reads have not checked every sequence (the first sequence in
    // ordinal order that fails, whichever the search's threads read first; the search's own
    // where the check finds none), and as its own after that (a file changed during the run).
    match error.downcast::<super::x_db_stream::XReadError>() {
        Ok(super::x_db_stream::XReadError(kept)) => database_error(match x_stream {
            Some(stream) if !stream.x_checked() => stream.x_first_error().unwrap_or(kept),
            _ => kept,
        }),
        Err(error) => engine(error),
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:172-175
/// NCBI's error of the options that `BlastxArgs::resolve` defers (`SemanticOptionsError`).
/// Exception (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`, exception 12): uses `USER_MANUAL_FOOTER`.
/// The message of a `SemanticOptionsError` keeps its bytes (probe R9 of line LIB-B).
pub(super) fn options_error(error: anyhow::Error) -> anyhow::Error {
    match error.downcast_ref::<super::args::SemanticOptionsError>() {
        Some(semantic) => crate::blastinput::app::options_error_bytes(&semantic.0),
        None => crate::blastinput::app::options_error(&error.to_string()),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:773-789
/// The warning of each `-subject` record without letters, as every program writes it
/// (`input_files::write_empty_subject_warnings_of`: the id, then the title's bytes as read).
pub(super) fn write_empty_subject_warnings(
    subjects: &[FastaRecord],
    stderr: &mut dyn Write,
) -> Result<()> {
    crate::blastinput::input_files::write_empty_subject_warnings_of(
        subjects.iter().map(|s| (s.names(), s.sequence.len())),
        "blastx",
        stderr,
    )?;
    Ok(())
}

/// A reader of a run's query batches (`run_query_batches`): it delivers each reader batch and
/// writes the bytes of each reader warning, in the reader's order, and keeps no record (the search
/// holds the batches it searches). The command line reads the query file
/// (`read_fasta_stream_batches`); LOSAT Web's `run_local` delivers the records its host registered.
pub(crate) type QueryBatchReader<'r> = dyn FnMut(
        &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
        &mut dyn FnMut(&[u8]) -> Result<()>,
    ) -> Result<()>
    + 'r;

/// The render of a query batch's reports (`run_query_batches`): the diagnostics, the batch's
/// records, its results, the number of queries reported before it, and the run's pool where no
/// search runs on it while the batch is written (its threads are free to format the batch's
/// reports, `x_pairwise_par::write`), else `None`.
pub(crate) type BatchRender<'r> = dyn FnMut(
        &mut dyn Write,
        &[FastaRecord],
        &mut Vec<BatchResults>,
        usize,
        Option<&crate::utils::threading::SearchPool<'_>>,
    ) -> Result<()>
    + 'r;

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:261-264 (cited at `run_inner`'s
// `take`: the batch's query factory, then its request, before its search)
/// The add of the reader batch `queries` (its first query `queries_before` of the run) to the
/// scopes of the queries (`query_scopes`), where the loop takes the batch: once per batch
/// (`taken`: the batch taken last, whose units the loop takes one by one), after the batch's reader
/// messages and the reports of the batches before, before the batch's request and search. LOSAT's
/// refusal of two queries of one scope with one id stops the run there, in the one form of every
/// program (`query_batch::refused_batch_error`; line WQA D10 round 5, 2026-11-20): with
/// `-export_search_strategy` (`export`) no request of the batch is written; else NCBI's error of
/// the set-up of the batch's records that NCBI's set-up reads with their own letters
/// (`reader_batch::set_up_error`, as `search_error` makes the search's), else the batch's search
/// where that set-up passes and the run's XML failure for its matrix name takes the refusal
/// (`ScopeRefusal::xml_takes`), else the refusal.
/// RPSTBLASTN's set-up reads its database's files first: refused. (Round 4 added the batch after
/// its search, and wrote its request.) The reader's batches are NCBI's (`input::emit_fasta_batch`
/// counts as `query_batch::BatchCount`; line WQA D10 round 6, 2026-11-20).
#[allow(clippy::too_many_arguments)]
pub(super) fn add_batch_scope(
    query_scopes: &std::cell::RefCell<Option<crate::blastinput::query_batch::QueryScopes>>,
    taken: &std::cell::Cell<Option<usize>>,
    queries: &[FastaRecord],
    queries_before: usize,
    export: bool,
    options: &super::args::ResolvedOptions,
    subjects: Subjects<'_>,
    search_error: &dyn Fn(anyhow::Error) -> anyhow::Error,
) -> Result<()> {
    if taken.replace(Some(queries_before)) == Some(queries_before) {
        return Ok(());
    }
    // Every record read, those that the range skipped among them (`input::scope_records`).
    let ordered = super::input::scope_records(queries);
    let Some(scopes) = query_scopes.borrow_mut().as_mut().map(|scopes| {
        scopes.add_batch_records(
            ordered
                .iter()
                .map(|&(believed, length, place)| (believed, length, place.is_some())),
            queries_before,
        )
    }) else {
        return Ok(());
    };
    scopes.or_else(|refusal| {
        use crate::blastinput::query_batch::BatchSetUp;
        crate::blastinput::query_batch::refused_batch_error(refusal, export, |first| {
            if options.rpstblastn.is_some() {
                return BatchSetUp::Unknown;
            }
            let records: Vec<FastaRecord> = ordered
                .iter()
                .zip(first)
                .filter_map(|(&(_, _, place), &first)| place.filter(|_| first))
                .map(|place| queries[place].clone())
                .collect();
            match super::reader_batch::set_up_error(&records, options, subjects) {
                BatchSetUp::Fails(error) => BatchSetUp::Fails(search_error(error)),
                other => other,
            }
        })
    })
}

/// What the query loop does where it takes a reader batch (`run_query_batches`' `take`) and where
/// it reports one, for the command line's run (`run_inner`) and LOSAT Web's search of records
/// (`local::render`), so that both keep one order of a batch's search strategy request: the batch
/// joins the scopes of the queries first (`add_batch_scope`); with `-export_search_strategy` its
/// request is written (`RunRequests::write`), and NCBI's posts of the bytes that its stream
/// replaced in it go to the run's diagnostics, as messages of the archive
/// (`StrategyExport::write_posts`); where the batch is reported, the search space that its set-up
/// wrote is noted for the next request (`note_searched`).
pub(super) struct BatchTake<'a> {
    query_scopes: &'a std::cell::RefCell<Option<crate::blastinput::query_batch::QueryScopes>>,
    scope_taken: std::cell::Cell<Option<usize>>,
    strategy: Option<&'a crate::report::archive::strategy::StrategyExport>,
    requests: Option<crate::report::archive::requests::RunRequests<'a>>,
    archive: Option<&'a crate::report::archive::ArchiveHandle>,
    held: Option<&'a super::blast_xml::HeldOutput>,
    options: &'a super::args::ResolvedOptions,
    subjects: Subjects<'a>,
    search_error: &'a dyn Fn(anyhow::Error) -> anyhow::Error,
}

impl<'a> BatchTake<'a> {
    /// The loop's takes of a run with the export `strategy`, the archive `archive` and, with
    /// `-subject`, its records `request_subjects` (the requests' subjects).
    pub(super) fn new(
        query_scopes: &'a std::cell::RefCell<Option<crate::blastinput::query_batch::QueryScopes>>,
        strategy: Option<&'a crate::report::archive::strategy::StrategyExport>,
        archive: Option<&'a crate::report::archive::ArchiveHandle>,
        options: &'a super::args::ResolvedOptions,
        subjects: Subjects<'a>,
        search_error: &'a dyn Fn(anyhow::Error) -> anyhow::Error,
        request_subjects: Option<&'a [FastaRecord]>,
    ) -> Self {
        Self {
            query_scopes,
            scope_taken: std::cell::Cell::new(None),
            strategy,
            requests: strategy.map(|_| archive::requests(request_subjects)),
            archive,
            held: None,
            options,
            subjects,
            search_error,
        }
    }

    // No NCBI counterpart: attaches LOSAT's `HeldOutput` so strategy request releases interleave
    // with held BLASTX/RPSTBLASTN XML/XML2/JSON reports (`HeldOutput::push_posts`).
    /// Attaches the run's `HeldOutput` so strategy request releases are queued with held
    /// single-stream XML/XML2/JSON reports (`HeldOutput::push_posts`).
    pub(super) fn with_held(mut self, held: Option<&'a super::blast_xml::HeldOutput>) -> Self {
        self.held = held;
        self
    }

    /// Takes the reader batch `queries`, after the `queries_before` queries of the batches
    /// before it.
    pub(super) fn take(
        &self,
        diagnostics: &mut dyn Write,
        queries: &[FastaRecord],
        queries_before: usize,
    ) -> Result<()> {
        add_batch_scope(
            self.query_scopes,
            &self.scope_taken,
            queries,
            queries_before,
            self.strategy.is_some(),
            self.options,
            self.subjects,
            self.search_error,
        )?;
        let (Some(strategy), Some(requests)) = (self.strategy, &self.requests) else {
            return Ok(());
        };
        requests.write(
            strategy,
            queries_before,
            queries,
            queries.first().map(FastaRecord::placement),
        )?;
        let posts = strategy.take_all_posts();
        if self.held.is_some_and(|held| held.push_posts(posts.clone())) {
            return Ok(());
        }
        Ok(posts.write(diagnostics, self.archive)?)
    }

    /// The batch `queries` (after `queries_before` queries) is reported with the results
    /// `results`: its set-up's `EffectiveSearchSpace`, which the next batch's request holds
    /// (`StrategyExport::note_searched`).
    pub(super) fn note_searched(
        &self,
        queries: &[FastaRecord],
        results: &[BatchResults],
        queries_before: usize,
    ) -> Result<()> {
        if let Some(strategy) = self.strategy {
            strategy.note_searched(
                queries_before..queries_before + queries.len(),
                archive::batch_searched_space(results, self.options),
            )?;
        }
        Ok(())
    }
}

/// The query batches of a run: each reader batch from `read` is searched with `options` against
/// `search_subjects` (with 2+ threads, on one pool for the whole run) and handed to `render` with
/// the number of queries reported before it, in the reader's order; the reader's warnings go to
/// `diagnostics` in that order. `serial_batches` searches and reports each batch before the next
/// is read (`-parse_deflines`). Returns the number of queries reported. `take` takes each reader
/// batch (each unit of a split one), with the number of queries reported before it, before its
/// search's result is looked at (the command line's search strategy request of that batch,
/// `native::run_inner`).
#[allow(clippy::too_many_arguments)]
pub(super) fn run_query_batches(
    options: &super::args::ResolvedOptions,
    search_subjects: Subjects<'_>,
    serial_batches: bool,
    read: &mut QueryBatchReader<'_>,
    diagnostics: &mut dyn Write,
    search_error: &dyn Fn(anyhow::Error) -> anyhow::Error,
    render: &mut BatchRender<'_>,
    take: &mut dyn FnMut(&mut dyn Write, &[FastaRecord], usize) -> Result<()>,
) -> Result<usize> {
    use super::x_timing::Run as XRun;
    #[cfg(not(feature = "parallel"))]
    let _ = serial_batches;
    let mut processed = 0usize;
    let x_t_loop = super::x_timing::lap();
    let x_batches = {
        // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:260-261
        // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:289-290
        // The reader and formatter use the same diagnostic stream at distinct
        // times. RefCell only serializes Rust borrows; it does not defer events.
        let diagnostics = std::cell::RefCell::new(diagnostics);
        // The batch loop as a function of the pool its searches run on, so that one pool can serve
        // all batches.
        // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:259-291
        // NCBI runs this loop: take the next batch of queries, search it, print its results.
        // x_run_batches is the same loop (the code below the batched branch is the original loop
        // body) as a closure over the pool its searches run on, so that one pool can serve all
        // batches.
        let mut x_run_batches = |x_pool: Option<&crate::utils::threading::SearchPool<'_>>| {
            // With a pool of 2+ threads: search several query batches at a time.
            //
            // A batch (`BATCH_SIZE` bases of queries) is searched on its own:
            // nothing is carried from one batch to the next except the order of
            // what is written. A file of many short queries is therefore a long
            // row of small searches whose serial parts (query masking, lookup
            // table, hand-over between stages) leave the other threads idle. Here
            // the reader's events are queued instead of acted on: a warning as the
            // bytes it would have written, a batch as a copy of its records. When
            // enough batches are queued they are searched side by side on the
            // pool, and the queue is then replayed in order: warnings are written,
            // each batch is reported exactly where the unbatched loop reported it.
            // A batch that fails ends the replay there, as it ended the loop; what
            // the reader delivered before it failed itself is replayed first.
            // A batch with a long query is split into chunks that use the whole
            // pool already; it is flushed at once, never queued with others.
            // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:345-357
            // NCBI has a similar mode (-mt_mode 1, x_RunMTBySplitQuery): several query batches are
            // searched at the same time as nodes of a master node. Here the batches are queued and
            // replayed in input order; the results printed are those of the one-batch-at-a-time
            // loop above.
            #[cfg(feature = "parallel")]
            {
                // Without a pool of 2+ threads the batches are searched one at a time;
                // LOSAT_X_BXFLOW acts only with the queued batches.
                let x_pooled = x_pool.is_some_and(|pool| pool.enabled());
                if !x_pooled {
                    crate::utils::xreach::x_reach!("LOSAT_X_BXBATCH", false);
                }
                if super::runtime::x_bx_flow() {
                    crate::utils::xreach::x_reach!("LOSAT_X_BXFLOW", x_pooled);
                }
            }
            #[cfg(feature = "parallel")]
            // With `-parse_deflines` a batch is searched and reported before the next is read:
            // the batches' scopes of query ids (`query_batch::QueryScopes`, added where each batch
            // is reported) and the reader's refusals keep the order of the one-batch loop. The
            // queued batches below have not been checked with them (line RG, Q1c).
            if let Some(pool) = x_pool.filter(|pool| pool.enabled() && !serial_batches) {
                use rayon::prelude::*;
                enum XEvent {
                    Warn(Vec<u8>),
                    Batch(XUnit),
                }
                let x_events = std::cell::RefCell::new(Vec::<XEvent>::new());
                let x_failed = std::cell::Cell::new(false);
                let mut pending = ReaderBatchResults::default();
                let x_wave = pool.threads() * 4;
                // A unit of more than two search units' letters (one long query, a reader batch
                // whose last query is long, or a reader batch searched whole) is searched with the
                // pool inside, never queued with others, and so is a unit of several places
                // (`XUnit::places`). A query is split into chunks that use the pool when it is the
                // only query of its reader batch; in a unit of several queries (NCBI splits none)
                // its subjects use the pool. Measured (2 queries of 280 kb, `BATCH_SIZE=1000000`,
                // 4 and 35 subjects, 8 and 16 threads): 0.10-0.16 s, against 0.17-0.26 s for the
                // same queries split.
                let x_large = |queries: &[super::input::FastaRecord]| {
                    queries.iter().map(|q| q.sequence.len()).sum::<usize>() > 2 * BATCH_SIZE
                };
                // A reader batch becomes its units; a batch of several is prepared first (phase
                // 1, `reader_batch.rs`): by the waves at once, on the pool, by the flow as tasks.
                let x_split = |queries: &[super::input::FastaRecord]| {
                    let units = XUnit::split(queries, options, search_subjects, false);
                    XUnit::prepare_now(&units, &|batch, ranges| {
                        prepare_reader_batch(batch, ranges, options, search_subjects, pool)
                    });
                    units
                };
                // EXPERIMENT (LOSAT_X_BXFLOW): the queued batches as a flow instead of waves.
                //
                // The waves below search `x_wave` batches and then wait for the last of them
                // before any is written or the next wave is read; at the end of the input the
                // remaining batches are searched one after another. Here every batch becomes a
                // task of its own once the reader has delivered as many batches after it as there
                // are threads, and this thread (worker 0 of the pool) writes the queued warnings
                // and reports in the reader's order as the results come in. The last batches of
                // the input (at most one per thread, the waves' tail) are each searched with the
                // pool inside, as the waves search their tail, but side by side. At most `x_wave`
                // batches are read and not yet written; at that bound, and at the end of the
                // input, this thread runs pending tasks itself and otherwise waits for a result.
                // Each batch is searched by a function the waves use, so the reports, the
                // warnings, their order, and the stop at a batch that fails are those of the
                // waves; only the time at which a batch is searched changes.
                // No NCBI counterpart: NCBI's -mt_mode 1 (blastx_app.cpp x_RunMTBySplitQuery, cited
                // above x_flush) also searches several query batches at once and prints them in
                // order; the flow changes only LOSAT's scheduling of the searches above.
                if super::runtime::x_bx_flow() {
                    // Dispatch point (LOSAT_X_BXUNITCAP=K, in the strict set): each batch searched
                    // with the pool inside (the subject flow) instead of one batch per thread, and
                    // at most K at a time: the flow spawns a batch once one batch follows it and
                    // while fewer than K are searched, and reads on while fewer than
                    // `flow_bound(K)` are unwritten (blastx/x_bx_mem.rs).
                    let x_cap =
                        super::search::x_bx_mem::flow_unit_cap(pool, options, search_subjects);
                    // The places (`XUnit::places`) scanning at a time: K with the cap, else the
                    // threads, at one place a thread.
                    let (x_bound, x_lookahead, x_max_live, x_max_places) = match x_cap {
                        Some(k) => (super::search::x_bx_mem::flow_bound(k), 1, k, k),
                        None => (x_wave, pool.threads(), usize::MAX, pool.threads()),
                    };
                    let x_ahead =
                        x_cap.is_some() && super::search::x_bx_mem::prepares_batches_ahead();
                    return x_flow(
                    x_bound,
                    x_lookahead,
                    x_max_live,
                    x_max_places,
                    &mut *read,
                    // Dispatch point (LOSAT_X_BXCAPFILL, in the strict set): in the capped
                    // flow a reader batch of one unit is prepared by a task of its own before
                    // its search, outside the K searched at a time (blastx/x_bx_mem.rs).
                    &|queries: &[super::input::FastaRecord]| {
                        XUnit::split(queries, options, search_subjects, x_ahead)
                    },
                    &|queries: &[super::input::FastaRecord],
                      layout: Option<&std::sync::Arc<ReaderBatchShape>>| match layout {
                        Some(layout) => prepare_unit(queries, layout, options, search_subjects),
                        None => prepare_lone_unit(queries, options, false),
                    },
                    &|unit: &XUnit, alone: bool, scan_end: &dyn Fn()| {
                        let queries = unit.queries();
                        let large = x_cap.is_some() || alone || unit.places > 1 || x_large(queries);
                        crate::utils::xreach::x_reach!("LOSAT_X_BXBATCH", !large);
                        if large {
                            let _x_live = x_cap
                                .map(|_| super::search::x_bx_mem::XCappedBatch::enter(pool));
                            super::runtime::x_search_unit_in_pool(
                                queries,
                                search_subjects,
                                options,
                                pool,
                                unit.take_prep(options, search_subjects)?,
                                Some(scan_end),
                            )
                        } else {
                            super::runtime::x_search_unit_serial(
                                queries,
                                search_subjects,
                                options,
                                unit.take_prep(options, search_subjects)?,
                                Some(scan_end),
                            )
                        }
                    },
                    &mut |reader_batch: &[super::input::FastaRecord],
                          start: usize,
                          last: bool,
                          results: Result<Vec<BatchResults>>|
                     -> Result<()> {
                        take(&mut **diagnostics.borrow_mut(), reader_batch, processed)?;
                        let results = match results {
                            Ok(results) => results,
                            Err(error) => {
    return Err(search_error(failed_batch_reads(
        error,
        options,
        search_subjects,
        &mut **diagnostics.borrow_mut(),
    )?))
}
                        };
                        if let Some(mut results) = pending.add(start, last, results) {
                            // The flow's searches of the batches after it run on the pool while
                            // this one is written: its reports are formatted on this thread
                            // (`BatchRender`).
                            render(
                                &mut **diagnostics.borrow_mut(),
                                reader_batch,
                                &mut results,
                                processed,
                                None,
                            )?;
                            processed += reader_batch.len();
                        }
                        Ok(())
                    },
                    &mut |bytes: &[u8]| -> Result<()> {
                        diagnostics.borrow_mut().write_all(bytes)?;
                        Ok(())
                    },
                );
                }
                // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:259-291
                // x_flush searches the queued batches (side by side, or one after another with the
                // whole pool inside each) and then writes warnings and results in the order the
                // reader produced them, which is the order of the loop above.
                let mut x_flush = |events: Vec<XEvent>| -> Result<()> {
                    let batches: Vec<&XUnit> = events
                        .iter()
                        .filter_map(|event| match event {
                            XEvent::Batch(unit) => Some(unit),
                            XEvent::Warn(_) => None,
                        })
                        .collect();
                    // Fewer batches than threads (the end of the input): each one
                    // uses the pool inside, as the unbatched loop does.
                    let x_t = super::x_timing::lap();
                    let x_tail = batches.len() < pool.threads();
                    crate::utils::xreach::x_reach!("LOSAT_X_BXBATCH", !x_tail);
                    let searched = if x_tail {
                        batches
                            .iter()
                            .map(|unit| {
                                super::runtime::x_search_unit_in_pool(
                                    unit.queries(),
                                    search_subjects,
                                    options,
                                    pool,
                                    unit.take_prep(options, search_subjects)?,
                                    None,
                                )
                            })
                            .collect::<Vec<_>>()
                    } else {
                        pool.install(|| {
                            batches
                                .par_iter()
                                .with_max_len(1)
                                .map(|unit| {
                                    let x_t = super::x_timing::lap();
                                    // A long unit uses the pool inside (a single long query
                                    // is split into chunks, otherwise its subjects use the
                                    // pool); a short one is one serial job.
                                    let queries = unit.queries();
                                    let searched = if unit.places > 1 || x_large(queries) {
                                        super::runtime::x_search_unit_in_pool(
                                            queries,
                                            search_subjects,
                                            options,
                                            pool,
                                            unit.take_prep(options, search_subjects)?,
                                            None,
                                        )
                                    } else {
                                        super::runtime::x_search_unit_serial(
                                            queries,
                                            search_subjects,
                                            options,
                                            unit.take_prep(options, search_subjects)?,
                                            None,
                                        )
                                    };
                                    super::x_timing::run(XRun::WaveBusy, x_t);
                                    searched
                                })
                                .collect::<Vec<_>>()
                        })
                    };
                    // Diagnostic (LOSAT_TIMING): the wall time of the searches. No NCBI
                    // counterpart.
                    let x_wall = super::x_timing::run(XRun::Search, x_t);
                    let (x_count, x_wall_stage) = if x_tail {
                        (XRun::TailBatches, XRun::TailWall)
                    } else {
                        super::x_timing::run_add(XRun::Waves, 1);
                        (XRun::WaveBatches, XRun::WaveWall)
                    };
                    super::x_timing::run_add(x_count, batches.len() as u64);
                    super::x_timing::run_add(x_wall_stage, x_wall);
                    let x_t = super::x_timing::lap();
                    let mut searched = searched.into_iter();
                    for event in &events {
                        match event {
                            XEvent::Warn(bytes) => diagnostics.borrow_mut().write_all(bytes)?,
                            XEvent::Batch(unit) => {
                                take(&mut **diagnostics.borrow_mut(), &unit.batch, processed)?;
                                let results =
                                    match searched.next().expect("one search per queued batch") {
                                        Ok(results) => results,
                                        Err(error) => {
                                            return Err(search_error(failed_batch_reads(
                                                error,
                                                options,
                                                search_subjects,
                                                &mut **diagnostics.borrow_mut(),
                                            )?))
                                        }
                                    };
                                if let Some(mut results) =
                                    pending.add(unit.range.start, unit.last, results)
                                {
                                    // The queued batches are searched: the pool is free.
                                    render(
                                        &mut **diagnostics.borrow_mut(),
                                        &unit.batch,
                                        &mut results,
                                        processed,
                                        Some(pool),
                                    )?;
                                    processed += unit.batch.len();
                                }
                            }
                        }
                    }
                    super::x_timing::run(XRun::Render, x_t);
                    Ok(())
                };
                let parsed = read(
                    &mut |queries| {
                        if queries.is_empty() {
                            return Err(engine("Empty CBlastQueryVector"));
                        }
                        // The reader batch is queued as its units (`XUnit::split`), which are
                        // searched as the batches of a run without
                        // `BATCH_SIZE` are; its report is written at its last unit.
                        for unit in x_split(queries) {
                            let large = unit.places > 1 || x_large(unit.queries());
                            let events = {
                                let mut events = x_events.borrow_mut();
                                events.push(XEvent::Batch(unit));
                                let queued = events
                                    .iter()
                                    .filter(|event| matches!(event, XEvent::Batch(_)))
                                    .count();
                                if !large && queued < x_wave {
                                    continue;
                                }
                                std::mem::take(&mut *events)
                            };
                            let flushed = x_flush(events);
                            if flushed.is_err() {
                                x_failed.set(true);
                                return flushed;
                            }
                        }
                        Ok(())
                    },
                    &mut |bytes: &[u8]| {
                        x_events.borrow_mut().push(XEvent::Warn(bytes.to_vec()));
                        Ok(())
                    },
                );
                if !x_failed.get() {
                    x_flush(x_events.take())?;
                }
                return parsed;
            }
            read(
                &mut |queries| {
                    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:259-262
                    // NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
                    if queries.is_empty() {
                        return Err(engine("Empty CBlastQueryVector"));
                    }
                    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
                    // Dispatch point: with a kept pool the batch is searched on it; without one,
                    // search_internal runs as before. Both searches compute the same result (the
                    // pool only changes who runs the work).
                    take(&mut **diagnostics.borrow_mut(), queries, processed)?;
                    let x_t = super::x_timing::lap();
                    let searched = match x_pool {
                        Some(pool) => super::runtime::x_search_internal_in_pool(
                            queries,
                            search_subjects,
                            options,
                            pool,
                        ),
                        None => search_internal(queries, search_subjects, options),
                    };
                    super::x_timing::run(XRun::Search, x_t);
                    let x_t = super::x_timing::lap();
                    let mut results = match searched {
                        Ok(results) => results,
                        Err(error) => {
                            return Err(search_error(failed_batch_reads(
                                error,
                                options,
                                search_subjects,
                                &mut **diagnostics.borrow_mut(),
                            )?))
                        }
                    };
                    reader_batch_skipped(&mut results);
                    render(
                        &mut **diagnostics.borrow_mut(),
                        queries,
                        &mut results,
                        processed,
                        x_pool,
                    )?;
                    super::x_timing::run(XRun::Render, x_t);
                    processed += queries.len();
                    Ok(())
                },
                &mut |bytes: &[u8]| Ok(diagnostics.borrow_mut().write_all(bytes)?),
            )
        };
        // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
        // Dispatch point: NCBI builds a CLocalBlast, and so its worker threads, for every batch.
        // With 2+ threads one pool is built for the whole run and every batch is searched on
        // it. Reuse of threads only.
        crate::utils::xreach::x_reach!("LOSAT_X_BXPOOL", options.num_threads > 1);
        if options.num_threads > 1 {
            // A pool that cannot be built is the engine error it was when each
            // batch built its own.
            match crate::utils::threading::x_with_search_pool_local(
                options.num_threads as usize,
                "BLASTX",
                |pool| {
                    // The run's counts of LOSAT_X_BXUNITCAP (scheduling only), on its pool.
                    #[cfg(feature = "parallel")]
                    let x_caps = super::search::x_bx_mem::XCapRun::default();
                    #[cfg(feature = "parallel")]
                    let pool = &x_caps.attach(pool);
                    Ok(x_run_batches(Some(pool)))
                },
            ) {
                Ok(batches) => batches,
                Err(error) => Err(engine(error)),
            }
        } else {
            x_run_batches(None)
        }
    };
    x_batches?;
    super::x_timing::run(XRun::Loop, x_t_loop);
    Ok(processed)
}

/// Runs the search; NCBI writes an error archive after the archives of the batches before when
/// the application stops with one of its own errors, and when the run ends with messages that no
/// archive carried (`PrintErrorArchive` after its `CATCH_ALL`, whatever the status), into the
/// output the run opened.
pub fn run(
    args: &BlastxArgs,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
) -> Result<()> {
    let mut output: Option<crate::cli::NcbiOstream<'static>> = None;
    let mut archive: Option<crate::report::archive::ArchiveHandle> = None;
    let mut strategy = None;
    let result = run_inner(args, settings, &mut output, &mut archive, &mut strategy);
    // A post of the search strategy left unwritten is LOSAT's internal error (`end_run`).
    let result = crate::report::archive::strategy::end_run(strategy.as_ref(), result);
    // `~CBlastFormat` and `CATCH_ALL`: a failed report ends the run as NCBI's does.
    let result = match output.as_ref().and_then(|output| output.report_handle()) {
        Some(report) => report.finish(result),
        None => result,
    };
    if let (Some(handle), Some(writer)) = (&archive, output.as_mut()) {
        match &result {
            Err(error) => {
                if let Some(lines) = archive::error_lines(error, false) {
                    let _ = crate::report::archive::write_error_archive_after_post(
                        handle, error, &lines, writer,
                    );
                } else if let Some(lines) = crate::report::archive::silent_end_lines(handle, error)
                {
                    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:307-310
                    // NCBI reference (598d8ae6): c++/src/app/blast/rpstblastn_app.cpp:232-235
                    // A run that `CATCH_ALL` ended without a post (`cli::NativeError::silent_end`,
                    // such as a strategy stream `CIOException::eWrite` fault): write the error
                    // archive of the messages saved before it.
                    let _ = crate::report::archive::write_error_archive(handle, &lines, writer);
                }
            }
            Ok(()) => {
                if handle.has_messages() {
                    let _ = crate::report::archive::write_error_archive(handle, &[], writer);
                }
            }
        }
        let _ = writer.flush();
    }
    result
}

/// The error archive of a run that stopped in `Run` before its output was opened (the taxid or
/// seqid list, the files of `-subject` and `-query`, the subjects): into the run's report stream
/// (`-out`, or standard output), as the other programs write it; nothing when `-out` cannot be
/// opened (NCBI's call swallows that).
fn write_run_error_archive(
    args: &BlastxArgs,
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
    let Some(lines) = archive::error_lines(error, true) else {
        return;
    };
    let out_file = match args.out.as_deref().filter(|p| *p != Path::new("-")) {
        Some(path) => match crate::cli::create_out_file(path, None) {
            Ok(file) => Some(file),
            Err(_) => return,
        },
        None => None,
    };
    let mut report =
        crate::cli::ReportStream::of_cli(out_file, args.program().lower, &args.outfmt, None);
    let _ =
        crate::report::archive::write_error_archive_after_post(handle, error, &lines, &mut report);
    let _ = report.finish(Ok(()));
}
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:989-1013
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2196-2217
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_exception.cpp:106-137
// NCBI reference (598d8ae6): c++/include/objtools/readers/fasta_exception.hpp:90-92
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2117-2125
// NCBI reference (598d8ae6): c++/include/corelib/ncbidiag.hpp:226-230
// NCBI reference (598d8ae6): c++/include/corelib/ncbidiag.hpp:747-747
// ParseDataLine precedes AssignMolType (fasta.cpp:409,1360), so fresh
// BLAST FASTA records have no molecule label when this warning is posted.
// LOG_POST_X uses the log flags without severity/prefix. Positions in the
// existing B reader state have already been converted to one-based indexes.
pub(crate) fn reader_warning(warning: &InputWarning, stderr: &mut dyn Write) -> Result<()> {
    // A warning of `ParseGapLine` (`fasta_reader::reader::gap_line_size`) is the bytes that the
    // reader of the other programs writes: a modifier's value as it is in the input.
    if warning.kind == "gap_line" {
        let mut line = Vec::with_capacity(warning.text.len() + 1);
        line.extend_from_slice(&warning.text);
        line.push(b'\n');
        stderr.write_all(&line)?;
        return Ok(());
    }
    let message=match warning.kind {
            "hyphens"=>format!("CFastaReader: Hyphens are invalid and will be ignored around line {}",warning.line),
            "title_nucleotide_residues"=>"FASTA-Reader: Title ends with at least 20 valid nucleotide characters.  Was the sequence accidentally put in the title line?".into(),
            "title_amino_acids"=>"FASTA-Reader: Title ends with at least 50 valid amino acid characters.  Was the sequence accidentally put in the title line?".into(),
            "invalid_residues"=>{
                let mut ranges:Vec<(usize,usize)>=Vec::new();
                for &p in &warning.positions {if let Some(last)=ranges.last_mut().filter(|r|r.1+1==p){last.1=p;}else{if ranges.len()>=1000{break;}ranges.push((p,p));}}
                let text=ranges.iter().map(|&(a,b)|if a==b{a.to_string()}else{format!("{a}-{b}")}).collect::<Vec<_>>().join(", ");
                format!("FASTA-Reader: Ignoring invalid {}residues at position(s): On line {}: {text}",String::from_utf8_lossy(&warning.text),warning.line)
            },
            "defline_id"=>String::from_utf8_lossy(&warning.text).into_owned(),
            _=>unreachable!("source-registered FASTA warning"),
        };
    writeln!(stderr, "{message}")?;
    Ok(())
}

/// NCBI searches a reader batch as a whole: it is not searched only when none of its contexts is
/// valid. The units of a batch are searched on their own, so a unit without a valid context is
/// marked skipped; `search_skipped` of all of them is set from the whole batch.
fn reader_batch_skipped(results: &mut [BatchResults]) {
    let skipped = results.iter().all(|results| results.search_skipped);
    for results in results {
        results.search_skipped = skipped;
    }
}

// No NCBI counterpart: how LOSAT schedules the search of a reader batch (NCBI's query batch,
// `BATCH_SIZE` letters). The search keeps units that end at the first query that brings them to
// 10002 letters (`batch_ranges` with the constant), whatever `BATCH_SIZE` is (the reader batch is
// one unit with `-use_sw_tback` and
// where the `-db_hard_mask` walk keeps a database sequence apart, `query_setup::unit_ranges`:
// there a unit cannot give the batch's rows alone): the flow and the waves schedule the units as
// they scheduled the batches of a run without `BATCH_SIZE`, and the report of a reader batch is
// written when its last unit is searched. The values that NCBI takes of the whole query batch
// (the split decision, the ungapped link cutoffs, the sum-statistics average length, the
// composition redo's minimum cutoff) are taken of the reader batch, not of the unit: phase 1
// prepares every unit of the batch and reduces them (`reader_batch.rs`), phase 2 searches the
// units with the result. The other values NCBI takes of the batch (the first valid context, the
// search-skipped flag) are taken by `ReaderBatchResults` below.
/// One unit of the search: `range` of a reader batch, which `last` ends.
#[cfg(feature = "parallel")]
#[derive(Clone)]
struct XUnit {
    batch: std::sync::Arc<Vec<FastaRecord>>,
    range: std::ops::Range<usize>,
    last: bool,
    /// The places the unit takes in the flow (`x_flow`, `query_setup::unit_places`): the search
    /// units that `batch_ranges` makes of its letters, more than one only for a reader batch of
    /// more than one such unit searched whole (`query_setup::searches_batch_whole`), whose scan
    /// holds a lookup table and the workers' diagonal tables of all its letters (its batch end
    /// holds neither: the flow counts the places until the scan ends); 1 for every other unit.
    places: usize,
    /// Phase 1 of the unit's reader batch when the batch has several units (or, in the capped batch
    /// flow, `LOSAT_X_BXCAPFILL`, is the unit itself), and the unit's place in it. `None`: the unit
    /// is the whole reader batch and prepares itself when it is searched.
    prep: Option<(std::sync::Arc<XReaderPrep>, usize)>,
}

/// Phase 1 of a reader batch of several units (`reader_batch.rs`), or of one unit prepared ahead
/// (`LOSAT_X_BXCAPFILL`): the units' preparations, and (the flow) the searches that wait for them.
#[cfg(feature = "parallel")]
struct XReaderPrep {
    /// The batch's values that follow from its records; each unit is prepared with them. `None`:
    /// a reader batch of one unit prepared ahead of its search (`LOSAT_X_BXCAPFILL`), which has no
    /// shape (`reader_batch::prepare_lone_unit`).
    layout: Option<std::sync::Arc<ReaderBatchShape>>,
    /// The units' preparations, from the end of phase 1 until each unit's search takes its own.
    /// All `None` if phase 1 failed: each unit is then prepared again when it is searched
    /// (`XUnit::take_prep`), and the first one whose preparation fails raises the error there (a
    /// reader batch with a search error writes nothing, so the values of the units before it
    /// are not used).
    slots: Vec<std::sync::Mutex<Option<UnitPrep>>>,
    /// The flow's tasks of phase 1 not yet ended, and whether one failed.
    remaining: std::sync::atomic::AtomicUsize,
    failed: std::sync::atomic::AtomicBool,
    /// Whether phase 1 has ended, and the units whose search was due before it ended (the flow:
    /// number, unit, searched alone).
    gate: std::sync::Mutex<(bool, Vec<(usize, XUnit, bool)>)>,
}

#[cfg(feature = "parallel")]
impl XReaderPrep {
    /// The end of phase 1: the batch's shape is reduced and attached to every preparation (or,
    /// if a unit failed, no unit keeps one), and the searches that waited are handed back.
    fn finish(&self) -> Vec<(usize, XUnit, bool)> {
        use std::sync::atomic::Ordering;
        let mut slots: Vec<_> = self
            .slots
            .iter()
            .map(|slot| slot.lock().expect("a unit's preparation"))
            .collect();
        if self.failed.load(Ordering::Acquire) {
            slots.iter_mut().for_each(|slot| **slot = None);
        } else if let Some(layout) = &self.layout {
            let shape = ReaderBatchShape::reduce(
                layout,
                slots
                    .iter()
                    .map(|slot| slot.as_ref().expect("a prepared unit")),
            );
            attach(
                slots
                    .iter_mut()
                    .map(|slot| slot.as_mut().expect("a prepared unit")),
                shape,
            );
        }
        drop(slots);
        let mut gate = self.gate.lock().expect("the reader batch's gate");
        gate.0 = true;
        std::mem::take(&mut gate.1)
    }
}

#[cfg(feature = "parallel")]
impl XUnit {
    /// The units of a reader batch, in order (the last one marked); a copy of the batch is
    /// shared by them. A batch of several units gets its phase 1 (`prepare_now`, or the flow's
    /// tasks), which its units' searches take; with `ahead` (the capped batch flow,
    /// `LOSAT_X_BXCAPFILL`) a batch of one unit of at most two search units' letters gets one too
    /// (the flow's task, without a shape). A longer batch (one long query, or a reader batch
    /// searched whole) prepares itself when it is searched, as before: preparations held ahead are
    /// bounded by letters, not only by the flow's count of unwritten batches.
    fn split(
        queries: &[FastaRecord],
        options: &super::args::ResolvedOptions,
        subjects: Subjects<'_>,
        ahead: bool,
    ) -> Vec<XUnit> {
        let batch = std::sync::Arc::new(queries.to_vec());
        let db = subjects.seqdb();
        let ranges = super::query_setup::unit_ranges(&batch, options, db);
        let count = ranges.len();
        let places = super::query_setup::unit_places(&batch, options, db);
        let ahead = ahead
            && count == 1
            && batch.iter().map(|q| q.sequence.len()).sum::<usize>() <= 2 * BATCH_SIZE;
        let prep = (count > 1 || ahead).then(|| {
            std::sync::Arc::new(XReaderPrep {
                layout: (count > 1)
                    .then(|| std::sync::Arc::new(ReaderBatchShape::layout(&batch, options))),
                slots: (0..count).map(|_| std::sync::Mutex::new(None)).collect(),
                remaining: std::sync::atomic::AtomicUsize::new(count),
                failed: std::sync::atomic::AtomicBool::new(false),
                gate: std::sync::Mutex::new((false, Vec::new())),
            })
        });
        ranges
            .into_iter()
            .enumerate()
            .map(|(index, range)| XUnit {
                batch: std::sync::Arc::clone(&batch),
                range,
                last: index + 1 == count,
                places,
                prep: prep
                    .as_ref()
                    .map(|prep| (std::sync::Arc::clone(prep), index)),
            })
            .collect()
    }

    /// Phase 1 of the units' reader batch (`units` from `split`), all at once by `prepare`.
    fn prepare_now(
        units: &[XUnit],
        prepare: &dyn Fn(&[FastaRecord], &[std::ops::Range<usize>]) -> Result<Vec<UnitPrep>>,
    ) {
        let Some((prep, _)) = units.first().and_then(|unit| unit.prep.as_ref()) else {
            return;
        };
        let ranges: Vec<_> = units.iter().map(|unit| unit.range.clone()).collect();
        if let Ok(preps) = prepare(&units[0].batch, &ranges) {
            for (slot, unit) in prep.slots.iter().zip(preps) {
                *slot.lock().expect("a unit's preparation") = Some(unit);
            }
        }
        prep.gate.lock().expect("the reader batch's gate").0 = true;
    }

    /// The flow's task of phase 1 for this unit: prepares it by `prepare` (with the batch's layout,
    /// or none for a batch of one unit prepared ahead); the last task of the reader batch ends
    /// phase 1 and returns the searches that waited for it.
    fn prepare_task(
        &self,
        prepare: &(dyn Fn(&[FastaRecord], Option<&std::sync::Arc<ReaderBatchShape>>) -> Result<UnitPrep>
              + Sync),
    ) -> Vec<(usize, XUnit, bool)> {
        use std::sync::atomic::Ordering;
        let Some((prep, position)) = self.prep.as_ref() else {
            return Vec::new();
        };
        // A panic is a failure of phase 1: the unit prepares itself again when it is searched,
        // and its search raises the panic where the unit is written. (With `panic = "abort"`, the
        // release profile, the process ends at the first panic instead, which may be before the
        // batches ahead of this one are written.)
        match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            prepare(self.queries(), prep.layout.as_ref())
        })) {
            Ok(Ok(unit)) => {
                *prep.slots[*position].lock().expect("a unit's preparation") = Some(unit)
            }
            _ => prep.failed.store(true, Ordering::Release),
        }
        if prep.remaining.fetch_sub(1, Ordering::AcqRel) == 1 {
            prep.finish()
        } else {
            Vec::new()
        }
    }

    /// The unit, if its search can start now: a unit of a reader batch whose phase 1 has not
    /// ended waits for it (number `index`, searched `alone`), and its search is started when
    /// phase 1 ends.
    fn ready(self, index: usize, alone: bool) -> Option<XUnit> {
        let Some((prep, _)) = self.prep.as_ref() else {
            return Some(self);
        };
        let prep = std::sync::Arc::clone(prep);
        let mut gate = prep.gate.lock().expect("the reader batch's gate");
        if gate.0 {
            return Some(self);
        }
        gate.1.push((index, self, alone));
        None
    }

    /// The unit's phase 1 preparation, once. If phase 1 failed, the unit of a batch of several is
    /// prepared again as phase 1 prepares it, after the split sizes that NCBI reads before the
    /// setup (with the reader batch's query count): the error is the reader batch's, as with one
    /// thread (`search::search_core`).
    fn take_prep(
        &self,
        options: &super::args::ResolvedOptions,
        subjects: Subjects<'_>,
    ) -> Result<Option<UnitPrep>> {
        let Some((prep, position)) = self.prep.as_ref() else {
            return Ok(None);
        };
        let taken = prep.slots[*position]
            .lock()
            .expect("a unit's preparation")
            .take();
        match (taken, &prep.layout) {
            (Some(unit), _) => Ok(Some(unit)),
            (None, Some(layout)) => {
                super::split::split_sizes(layout.query_count)?;
                prepare_unit(self.queries(), layout, options, subjects).map(Some)
            }
            (None, None) => Ok(None),
        }
    }

    fn queries(&self) -> &[FastaRecord] {
        &self.batch[self.range.clone()]
    }
}

/// The results of the units of a reader batch searched so far. At the batch's last unit they are
/// handed over as the results of the reader batch: `query_ordinal` counts from the reader
/// batch's first query, and `search_skipped` is set for all of them when no unit had a valid
/// context (NCBI searches the batch, or does not, as a whole).
#[cfg(feature = "parallel")]
#[derive(Default)]
struct ReaderBatchResults {
    results: Vec<BatchResults>,
}

#[cfg(feature = "parallel")]
impl ReaderBatchResults {
    fn add(
        &mut self,
        start: usize,
        last: bool,
        mut unit: Vec<BatchResults>,
    ) -> Option<Vec<BatchResults>> {
        for results in &mut unit {
            results.query_ordinal += start;
            // The report reads the contexts of the prepared queries and nothing of their
            // residues, which are not needed once the unit is searched; the results of a unit
            // wait here for the last unit of the batch.
            results.prepared.sequence_start = Vec::new();
            results.prepared.sequence_start_nomask = Vec::new();
            results.prepared.lookup_segments = Vec::new();
        }
        self.results.append(&mut unit);
        if !last {
            return None;
        }
        let mut results = std::mem::take(&mut self.results);
        reader_batch_skipped(&mut results);
        Some(results)
    }
}

/// A unit's search in the flow, or the panic it ended with (raised again where it is written).
#[cfg(feature = "parallel")]
type XSearched = std::thread::Result<Result<Vec<BatchResults>>>;

// EXPERIMENT (LOSAT_X_BXFLOW): the reading loop of `run` with every query batch searched as a task
// of its own (see the call in `run`). `search` searches one batch (the waves' function; `alone`:
// with the pool inside), and calls its third argument when the batch's scan has ended (before its
// redo and traceback, `runtime::x_search_unit_in_pool`), `write_batch` writes one batch's reports
// (or ends the run at its error),
// `write_warning` writes the bytes of one reader warning. A batch is spawned once `lookahead`
// batches follow it and while fewer than `max_live` are spawned and not yet searched (by the
// reading thread or by a search that ends); the batches still held when the input ends are
// searched alone. At most `bound` batches are read and not yet written. `split` makes a reader
// batch its units (the flow's batches); a reader batch of several units (in the capped flow with
// `LOSAT_X_BXCAPFILL`, every reader batch) is prepared by `prepare` (phase 1, one task per unit)
// before any of its units is searched. Returns the reader's error, if any (the reading stops at
// the first batch whose search or write fails).
// No NCBI counterpart: scheduling of the batch searches only (NCBI's -mt_mode 1 runs batches as
// nodes and prints them in order, blastx_app.cpp:345-357, cited in `run`); every value written
// is computed by the same calls as without the switch, in the same order of writing.
#[cfg(feature = "parallel")]
#[allow(clippy::too_many_arguments)]
fn x_flow(
    bound: usize,
    lookahead: usize,
    max_live: usize,
    max_places: usize,
    read: &mut QueryBatchReader<'_>,
    split: &dyn Fn(&[FastaRecord]) -> Vec<XUnit>,
    prepare: &(dyn Fn(&[FastaRecord], Option<&std::sync::Arc<ReaderBatchShape>>) -> Result<UnitPrep>
          + Sync),
    search: &(dyn Fn(&XUnit, bool, &dyn Fn()) -> Result<Vec<BatchResults>> + Sync),
    write_unit: &mut dyn FnMut(
        &[FastaRecord],
        usize,
        bool,
        Result<Vec<BatchResults>>,
    ) -> Result<()>,
    write_warning: &mut dyn FnMut(&[u8]) -> Result<()>,
) -> Result<()> {
    use super::x_timing::Run as XRun;
    use std::cell::{Cell, RefCell};
    use std::collections::{BTreeMap, VecDeque};
    use std::sync::mpsc;
    // A reader event, written in the reader's order: a warning's bytes, or a unit by its number.
    enum XEvent {
        Warn(Vec<u8>),
        Batch(usize, XUnit),
    }
    struct XState {
        events: VecDeque<XEvent>,
        searched: BTreeMap<usize, XSearched>,
        read: usize,
        unwritten: usize,
    }
    // Units read and not yet spawned (the newest `lookahead` until the end of the input), in
    // reader order, the number spawned and not yet searched (at most `max_live`; a unit that
    // waits for its reader batch's phase 1 is spawned), and the places (`XUnit::places`) of those
    // whose scan has not ended (a unit gives back its places at its batch end, `spawn_search`).
    struct XHeld {
        batches: VecDeque<(usize, XUnit)>,
        live: usize,
        places: usize,
        ended: bool,
        // The writing ended (at the end, at an error, or at a panic it raises again): spawn no more.
        stopped: bool,
    }
    struct XLaunch<'f> {
        held: std::sync::Mutex<XHeld>,
        lookahead: usize,
        max_live: usize,
        max_places: usize,
        search: &'f (dyn Fn(&XUnit, bool, &dyn Fn()) -> Result<Vec<BatchResults>> + Sync),
    }
    impl XLaunch<'_> {
        // The held units; a panic elsewhere cannot leave them unusable (nothing panics while
        // they are locked, and a search's result must still be sent).
        fn held(&self) -> std::sync::MutexGuard<'_, XHeld> {
            self.held
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner())
        }

        // Spawns the held units beyond the newest `lookahead` (all of them at the end of the
        // input, then searched alone) while fewer than `max_live` are spawned and unsearched,
        // and a unit of several places (`XUnit::places`) while its places and those of the
        // spawned units still scanning are at most `max_places`, or none is scanning: a reader
        // batch searched whole takes the places of the units the split would make of it, so
        // that the letters scanned at a time (each scan holds a lookup table and diagonal tables
        // of its letters) stay those of the split's units (one place each), and one of more
        // places than `max_places` is scanned with no other unit. Its batch end (the redo and
        // the traceback) holds none of those tables, and leaves the threads it does not use to
        // the next unit's scan (`spawn_search`). A unit of a reader batch whose phase 1 has not
        // ended is searched when it ends.
        fn release<'s>(
            &'s self,
            scope: &rayon::Scope<'s>,
            sender: &mpsc::Sender<(usize, XSearched)>,
        ) {
            loop {
                let next = {
                    let mut held = self.held();
                    let keep = if held.ended { 0 } else { self.lookahead };
                    let fits = held.batches.front().is_some_and(|(_, unit)| {
                        unit.places <= 1
                            || held.places == 0
                            || held.places + unit.places <= self.max_places
                    });
                    if !held.stopped
                        && held.batches.len() > keep
                        && held.live < self.max_live
                        && fits
                    {
                        held.live += 1;
                        let places = held.batches.front().map_or(1, |(_, unit)| unit.places);
                        held.places += places;
                        let alone = held.ended;
                        held.batches
                            .pop_front()
                            .map(|(index, unit)| (index, unit, alone))
                    } else {
                        None
                    }
                };
                let Some((index, unit, alone)) = next else {
                    break;
                };
                if alone {
                    super::x_timing::run_add(XRun::FlowAlone, 1);
                }
                if let Some(unit) = unit.ready(index, alone) {
                    self.spawn_search(scope, sender, index, unit, alone);
                }
            }
        }

        // Spawns the search of unit `index`, spawned by `release`. The search gives back the
        // unit's places when its scan ends and spawns again (the units that wait for places
        // start their scans beside its batch end), sends its result (or panic) once and then
        // spawns again, so a place freed by a search is filled without waiting for the writing
        // thread. A unit whose phase 1 ends after the writing ended is not searched.
        fn spawn_search<'s>(
            &'s self,
            scope: &rayon::Scope<'s>,
            sender: &mpsc::Sender<(usize, XSearched)>,
            index: usize,
            unit: XUnit,
            alone: bool,
        ) {
            {
                let mut held = self.held();
                if held.stopped {
                    held.live -= 1;
                    held.places -= unit.places;
                    return;
                }
            }
            let sender = sender.clone();
            scope.spawn(move |scope| {
                // The places are given back once: at the end of the scan, or with the search
                // (a batch without a valid context, or an error before its batch end).
                let given = std::cell::Cell::new(false);
                let scan_end = || {
                    if !given.replace(true) {
                        self.held().places -= unit.places;
                        self.release(scope, &sender);
                    }
                };
                let searched = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    (self.search)(&unit, alone, &scan_end)
                }));
                {
                    let mut held = self.held();
                    held.live -= 1;
                    if !given.get() {
                        held.places -= unit.places;
                    }
                }
                // The receiver is gone only when the run already ended at an error.
                let _ = sender.send((index, searched));
                self.release(scope, &sender);
            });
        }
    }
    debug_assert!(rayon::current_thread_index().is_some());
    debug_assert!(lookahead < bound);
    let launch = XLaunch {
        held: std::sync::Mutex::new(XHeld {
            batches: VecDeque::new(),
            live: 0,
            places: 0,
            ended: false,
            stopped: false,
        }),
        lookahead,
        max_live,
        max_places,
        search,
    };
    // Marks the writing ended when the scope's body returns or unwinds, before the scope waits
    // for the searches it spawned: the searches that end then spawn no more.
    struct XStop<'a, 'f>(&'a XLaunch<'f>);
    impl Drop for XStop<'_, '_> {
        fn drop(&mut self) {
            self.0.held().stopped = true;
        }
    }
    rayon::in_place_scope(|scope| {
        let _stop = XStop(&launch);
        let (sender, receiver) = mpsc::channel::<(usize, XSearched)>();
        let state = RefCell::new(XState {
            events: VecDeque::new(),
            searched: BTreeMap::new(),
            read: 0,
            unwritten: 0,
        });
        let failed = Cell::new(false);
        let release = || launch.release(scope, &sender);
        // Writes the events at the front whose batches are searched; then, until fewer than
        // `bound` batches are unwritten (or, with `all`, until every event is written), runs a
        // pending task or waits for a result, and writes again.
        let mut write = |all: bool| -> Result<()> {
            loop {
                {
                    let mut state = state.borrow_mut();
                    while let Ok((index, searched)) = receiver.try_recv() {
                        state.searched.insert(index, searched);
                    }
                }
                loop {
                    let event = {
                        let mut state = state.borrow_mut();
                        let ready = match state.events.front() {
                            Some(XEvent::Warn(_)) => true,
                            Some(XEvent::Batch(index, _)) => state.searched.contains_key(index),
                            None => false,
                        };
                        if ready {
                            state.events.pop_front()
                        } else {
                            None
                        }
                    };
                    let x_t = super::x_timing::lap();
                    match event {
                        None => break,
                        Some(XEvent::Warn(bytes)) => write_warning(&bytes)?,
                        Some(XEvent::Batch(index, unit)) => {
                            let searched = {
                                let mut state = state.borrow_mut();
                                state.unwritten -= 1;
                                state.searched.remove(&index).expect("a searched batch")
                            };
                            match searched {
                                Ok(results) => {
                                    write_unit(&unit.batch, unit.range.start, unit.last, results)?
                                }
                                Err(panic) => std::panic::resume_unwind(panic),
                            }
                        }
                    }
                    super::x_timing::run(XRun::Render, x_t);
                }
                {
                    let state = state.borrow();
                    if (all && state.events.is_empty()) || (!all && state.unwritten < bound) {
                        return Ok(());
                    }
                }
                let x_t = super::x_timing::lap();
                if let Some(rayon::Yield::Executed) = rayon::yield_now() {
                    let x_ns = super::x_timing::run(XRun::FlowJobTime, x_t);
                    super::x_timing::run_add(XRun::FlowJobs, 1);
                    super::x_timing::run_add(XRun::Search, x_ns);
                } else {
                    // No task is pending now: the unwritten batches are searched on other
                    // threads (each sends its result once). Wait for a result, but look for
                    // pending work again soon: a batch searched alone hands out parts of its
                    // search as it goes.
                    match receiver.recv_timeout(std::time::Duration::from_millis(1)) {
                        Ok((index, searched)) => {
                            state.borrow_mut().searched.insert(index, searched);
                        }
                        Err(mpsc::RecvTimeoutError::Timeout) => {}
                        Err(mpsc::RecvTimeoutError::Disconnected) => {
                            unreachable!("the reading thread keeps a sender")
                        }
                    }
                    let x_ns = super::x_timing::run(XRun::FlowBlocked, x_t);
                    super::x_timing::run_add(XRun::Search, x_ns);
                }
            }
        };
        let parsed = read(
            &mut |queries| {
                if queries.is_empty() {
                    return Err(engine("Empty CBlastQueryVector"));
                }
                // The reader batch enters the flow as its units (`XUnit::split`), each as a
                // batch of a run without `BATCH_SIZE` does; the bound counts units. A
                // reader batch of several units (in the capped flow with LOSAT_X_BXCAPFILL, every
                // reader batch) is prepared first (phase 1) by one task per unit,
                // beside the searches of the batches before it, while this thread reads on: the
                // last task ends phase 1 and starts the searches that were due before it ended
                // (no thread waits for phase 1).
                let units = split(queries);
                if units.first().is_some_and(|unit| unit.prep.is_some()) {
                    for unit in &units {
                        let (unit, sender, launch) = (unit.clone(), sender.clone(), &launch);
                        scope.spawn(move |scope| {
                            for (index, unit, alone) in unit.prepare_task(prepare) {
                                launch.spawn_search(scope, &sender, index, unit, alone);
                            }
                        });
                    }
                }
                for unit in units {
                    {
                        let mut state = state.borrow_mut();
                        let index = state.read;
                        state.read += 1;
                        state.unwritten += 1;
                        state.events.push_back(XEvent::Batch(index, unit.clone()));
                        let mut held = launch.held();
                        held.batches.push_back((index, unit));
                    }
                    super::x_timing::run_add(XRun::FlowBatches, 1);
                    release();
                    let written = write(false);
                    if written.is_err() {
                        failed.set(true);
                        return written;
                    }
                }
                Ok(())
            },
            &mut |bytes: &[u8]| {
                state
                    .borrow_mut()
                    .events
                    .push_back(XEvent::Warn(bytes.to_vec()));
                Ok(())
            },
        );
        if !failed.get() {
            // The end of the input: the held batches (the last of the input) are searched alone,
            // side by side.
            launch.held().ended = true;
            release();
            write(true)?;
        }
        parsed
    })
}

/// The formatting handler's warnings of a BLASTX run with output format `outfmt`
/// (`ReportCounts::warnings`).
pub(crate) fn warnings_of(args: &BlastxArgs, outfmt: i32) -> Vec<u8> {
    crate::blastinput::report_counts::report_counts(
        outfmt,
        super::args::max_target_seqs_count(args.max_target_seqs),
        args.num_descriptions,
        args.num_alignments,
    )
    .warnings_with(
        args.program().lower,
        outfmt,
        args.num_descriptions.is_some(),
        &args.formatting,
    )
}

#[cfg(test)]
#[path = "native_archive_error_tests.rs"]
mod archive_error_tests;

#[cfg(all(test, feature = "parallel"))]
#[path = "native_reader_batch_tests.rs"]
mod reader_batch_tests;
