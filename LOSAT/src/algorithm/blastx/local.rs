//! BLASTX over records in memory: one search over the records that LOSAT Web registered or the
//! library API read (`run_local`), or of the library's records against a database
//! (`run_local_database`), written in every requested output format. The query batches are
//! searched and reported by the command line's loop (`native::run_query_batches`).

use super::{
    archive,
    args::{BlastxArgs, ResolvedOptions, SemanticOptionsError},
    database::{BlastxDatabase, ReportedSubjects},
    empty_rows::SubjectLoc,
    input::FastaRecord,
    input_source::{from_shared_ref, BatchCutter, CutSettings},
    native::{
        archive_error, database_search_error, engine, options_error, report_error, reported_batch,
        run_query_batches, write_empty_subject_warnings, ReportedAsks,
    },
    query_setup::BATCH_SIZE,
    report::{self, LocalReport, RecordsParsed, ShownSubjects},
    runtime::BatchResults,
    subjects::Subjects,
};
use crate::algorithm::common::database::library_subject;
use crate::api::local_blast::{FormatProbe, ReportOutputs, ReportedSubject};
use crate::api::QueryEnd;
use crate::blastinput::app::{self, Program};
use crate::blastinput::fasta_reader;
use crate::blastinput::seq_range::{parse_optional_range, record_cut, RangeRole, SequenceRange};
use anyhow::{bail, ensure, Result};
use std::io::{self, Write};

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:429-467
// NCBI formats one result set without searching again; several requested formats are several
// CBlastFormat printers over the same result set.
/// Runs one BLASTX search over already parsed records and writes every requested output format
/// from the same result (the entry of web ABI v2; the command line is `native::run`), and gives
/// the caller the final hit list (`ReportOutputs::hits`).
///
/// Each requested `-outfmt` is read as on the command line, and NCBI's checks of the options and
/// of the subjects come in the command line's order, before the search. The `-query` and
/// `-subject` values of `args` are used only as display names. Each format is the command line's
/// writer of it (`report::render_queries`); the side files of outfmt 13 and 14 and `-db` are
/// refused. Under `-parse_deflines` the records are those that the reader read with the ids of
/// their deflines (`api::input::read_search_inputs`, which LOSAT Web reads them with too), as for
/// the other programs. LOSAT Web asks for outfmt 0, 6 and
/// 7, 6 and 7 with a custom field list or delimiter too. `-subject_loc` and `-query_loc` cut the records as the command line cuts them: each
/// subject as it is read, each query as its batch is read (`native::run_inner`).
///
/// The records are those of NCBI's reader (`fasta_reader`): the queries as its nucleotide reader
/// reads them, the subjects as its protein reader reads them (`input_source::from_shared_ref`). Their messages are
/// written where BLASTX's reader writes them: a subject's as the subjects are read, a query's
/// before the report of the batch that holds it. Queries without a record are NCBI's
/// `Query is Empty!`.
pub fn run_local(
    mut args: BlastxArgs,
    query_records: &[fasta_reader::FastaRecord],
    query_end: &QueryEnd,
    subject_records: &[fasta_reader::FastaRecord],
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // NCBI's options factory checks the task first (`TaskParser::check`); the command line's
    // parser has checked it, a caller of `run_local` has not.
    super::args::BLASTX_TASK.check(&mut args.task)?;
    // A search of records reads no database: `run_local_database` searches the database of `-db`.
    // `check_options` (LOSAT Web's `validate`) checks `-db` as the command line does (line WQA D12a).
    if args.db.is_some() {
        bail!(crate::api::local_blast::RUN_LOCAL_DB_REFUSAL);
    }
    let requested = requested_formats(outputs);
    let (choices, resolutions) = read_options(&args, &requested)?;
    crate::utils::threading::validate_threads(args.num_threads as usize)?;
    let archive_requested = choices.iter().any(|choice| choice.number() == Some(11));
    // The archive shares the options handle of the run's search strategy (LOSAT Web's
    // `-export_search_strategy`, `web_adapter::give_strategy`), as on the command line
    // (`native::RunSetup::open`).
    let archive = archive::handle_for(
        &args,
        archive_requested,
        &crate::report::archive::strategy::options_of(args.strategy.as_ref()),
    )
    .map(|handle| outputs.with_archive_outfmt(handle));
    let no_query = query_records.is_empty() && query_end.is_input();
    outputs.run_records(
        args.strategy.as_ref(),
        archive.as_ref(),
        no_query,
        |outputs| {
            search_records(
                &args,
                query_records,
                query_end,
                subject_records,
                (choices, resolutions),
                archive.as_ref(),
                outputs,
            )
        },
    )
}

/// `run_local` from the reading of the subjects on.
fn search_records(
    args: &BlastxArgs,
    query_records: &[fasta_reader::FastaRecord],
    query_end: &QueryEnd,
    subject_records: &[fasta_reader::FastaRecord],
    (choices, resolutions): ReadOptions,
    archive: Option<&crate::report::archive::ArchiveHandle>,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // An imported strategy's subjects stand for `-subject`, with no name (`native::run_inner`).
    let subject_path = match (
        args.subject.as_deref(),
        args.import
            .as_ref()
            .is_some_and(|run| run.subjects.is_some()),
    ) {
        (Some(path), _) => path,
        (None, true) => std::path::Path::new(""),
        (None, false) => return Err(app::missing_subject_error()),
    };
    let diagnostics = &mut outputs.diagnostics;
    let subject_range = parse_subject_range(args)?;
    // NCBI reads the subjects before the queries, each record's messages as it is read, and cuts
    // each to its `-subject_loc` interval: a record whose interval starts more than one letter
    // past its end stops the reading after its messages (`native::run_inner`).
    let mut subjects = Vec::with_capacity(subject_records.len());
    for record in subject_records {
        let mut subject = from_shared_ref(record, args.lcase_masking)?;
        diagnostics.write_all(&record.warnings)?;
        if let Some(range) = &subject_range {
            let (from, to_exclusive) = record_cut(range, subject.sequence.len())?;
            // The lower case of the whole record, which outfmt 0's Sbjct rows and the hits show
            // (`x_soft_subject::x_set_subject_row_masks`), and the archive's Bioseqs.
            let keep_whole = archive.is_some()
                || args.strategy.is_some()
                || (args.lcase_masking && !subject.lowercase_masks.is_empty());
            subject.cut_in_place(from, to_exclusive, keep_whole);
        }
        subjects.push(subject);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
    if subjects.is_empty() {
        return Err(engine("Empty CBlastQueryVector"));
    }
    // `-out` (`run_to_file`) is made where the command line opens it (`native::run`).
    crate::api::local_blast::open_out(&mut outputs.formats)?;
    check_filtering(args)?;
    let query_range = parse_query_range(args)?;
    let (choices, mut resolved) =
        check_resolutions(args, &choices, resolutions, &mut **diagnostics)?;
    for options in &mut resolved {
        options.query_range = query_range;
    }
    // The taxonomy database of the reports, where the command line asks (`native::run`).
    subject_taxonomy(&mut resolved, args.taxdb_setting.as_ref());
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:213-217
    // Records of an input read to its end; an end that reads past them (blank lines, a reader
    // error) is met in the first batch (`reader_batches`).
    if query_records.is_empty() && query_end.is_input() {
        diagnostics
            .write_all(format!("Warning: [{}] Query is Empty!\n", args.app_name()).as_bytes())?;
        return Ok(());
    }
    // The letters of a reader batch (`native::run_inner`).
    let batch_size = app::query_batch_size(BATCH_SIZE as u32)? as usize;
    let xml_matrix = super::blast_xml::xml_matrix_failure(&resolved, true);
    let xml_scan = scan_xml_titles(
        &resolved,
        query_records,
        query_end,
        args.lcase_masking,
        query_range,
        batch_size,
        archive.is_some() || args.strategy.is_some(),
        true,
        xml_matrix,
        args.parse_deflines,
        args.scope_refusal.as_ref(),
    )?;
    write_empty_subject_warnings(&subjects, args.app_name(), &mut **diagnostics)?;
    // NCBI's error of an XInclude format without a file name, where the command line gives it
    // (`native::run_inner`, blastx_app.cpp:251).
    outputs.xinclude_check(args.out.as_deref())?;
    // The environment of a search of subject records, after the subjects' warnings
    // (`native::run_inner`).
    app::check_search_environment("BLASTX", true)?;
    let records_parsed = if args.parse_deflines {
        // LOSAT's refusals of the records in outfmt 0 to 4, each such format's (`FormatRefusal`).
        RecordsParsed::of(&subjects, args.show_gis, |parsed| {
            for (format, options) in outputs.formats.iter_mut().zip(&resolved) {
                if options.pairwise_view().is_some() {
                    format
                        .refusal
                        .check(crate::report::db_subject::check_pairwise_records(
                            parsed, true, "BLASTX",
                        ))?;
                }
            }
            Ok(())
        })?
    } else {
        RecordsParsed::default()
    };
    let shown = ShownSubjects::Fasta {
        path: subject_path,
        parsed: &records_parsed,
    };
    // The errors of a `-subject` search (`native::run_inner`).
    let search_error = |error: anyhow::Error| -> anyhow::Error {
        // LOSAT's refusal, not an error of NCBI's engine, and an error that has NCBI's report
        // already (`native::search_error_as_is`).
        if super::native::search_error_as_is(&error) {
            super::native::search_error_held(error)
        } else {
            engine(error)
        }
    };
    report_batches(
        RunReports {
            lowercase: args.lcase_masking,
            query_records,
            query_end,
            query_range,
            batch_size,
            choices: &choices,
            resolved: &resolved,
            subjects: &subjects,
            shown,
            search_subjects: Subjects::Fasta(&subjects),
            serial_batches: args.parse_deflines,
            parse_deflines: args.parse_deflines,
            scope_refusal: args.scope_refusal.as_ref(),
            database: None,
            archive,
            strategy: args.strategy.as_ref(),
            out: args.out.as_deref(),
            search_error: &search_error,
            xml_scan,
            xml_matrix,
        },
        outputs,
    )
}

/// What the reports `resolved` of a `-subject` search ask of the taxonomy database, through
/// `setting` (`common::database::subject_taxonomy`): asked by the command line (`native::run`)
/// and by a search of records (`search_records`) at the same point, after the options are
/// resolved. The reports share the organism report's database, as the other programs' formats.
pub(super) fn subject_taxonomy(
    resolved: &mut [ResolvedOptions],
    setting: Option<&crate::algorithm::common::database::TaxdbSetting>,
) {
    let taxonomy = crate::algorithm::common::database::subject_taxonomy(
        setting,
        resolved.iter().filter_map(|options| options.tabular),
    );
    for options in resolved {
        options.taxdb_found = taxonomy.found.unwrap_or(false);
        options.subject_taxdb = taxonomy.taxdb.clone();
    }
}

/// A `-db` search's part of `report_batches`.
struct DatabaseReports<'a> {
    database: &'a BlastxDatabase,
    pairwise: &'a crate::report::pairwise::PairwiseDatabase,
}

/// What `report_batches` reports: the checks of the options and of the inputs are done.
struct RunReports<'a> {
    /// `-lcase_masking`.
    lowercase: bool,
    query_records: &'a [fasta_reader::FastaRecord],
    /// How the reading of the records ended (`reader_batches`).
    query_end: &'a QueryEnd,
    query_range: Option<SequenceRange>,
    /// The letters of a reader batch.
    batch_size: usize,
    choices: &'a [app::FormatChoice],
    resolved: &'a [ResolvedOptions],
    /// The subject records; none for a database.
    subjects: &'a [FastaRecord],
    /// The subjects as the prolog and the epilog show them.
    shown: ShownSubjects<'a>,
    search_subjects: Subjects<'a>,
    /// `run_query_batches`.
    serial_batches: bool,
    /// `-parse_deflines` (`query_batch::QueryScopes`).
    parse_deflines: bool,
    /// `BlastxArgs::scope_refusal`: the library's refusal of the call's batch.
    scope_refusal: Option<&'a crate::blastinput::query_batch::ScopeRefusal>,
    database: Option<DatabaseReports<'a>>,
    /// The BLAST archive (outfmt 11), when it is requested.
    archive: Option<&'a crate::report::archive::ArchiveHandle>,
    /// The search strategy file (`-export_search_strategy`, LOSAT Web's), when it is requested.
    strategy: Option<&'a crate::report::archive::strategy::StrategyExport>,
    /// The file of `api::run_to_file`, next to which outfmt 13 and 14 write their side files.
    out: Option<&'a std::path::Path>,
    search_error: &'a dyn Fn(anyhow::Error) -> anyhow::Error,
    /// The read of the queries before the search of BLAST XML and the XML2 and JSON documents
    /// (`scan_xml_titles`).
    xml_scan: super::blast_xml::QueryScan,
    /// Where the run's XML formats fail for a matrix name that NCBI's XML does not know
    /// (`blast_xml::xml_matrix_failure`, decided once for the run).
    xml_matrix: Option<crate::report::blast_xml::matrix::XmlMatrixFailure>,
}

/// BLAST XML and the XML2 and JSON documents refuse a query title that NCBI crashes on
/// (`names::Refusal`) before the first byte, in the reports that they write; a `-subject` XML2 or
/// JSON document never writes
/// the last query of a run that stops at a later batch (`blast_xml::scan_queries`, on the records
/// in memory, read into the search's reader batches as the command line's `prescan_queries` reads
/// the file, and at the same place: before the subjects' warnings). With outfmt 5 among the
/// formats every report is written. `keep_whole`: the search keeps the whole letters of a cut
/// record (the archive's). Under `-parse_deflines` (`parse_deflines`) the read adds each batch
/// that the reader delivers to the scopes of the queries (`QueryScan::query_scopes`, read ahead),
/// as the command line's read does (`blast_xml::scan_query_stream`): LOSAT's refusal of two
/// queries of one scope with one id stops the run at its batch (`QueryScopes::read_ahead`; line WQA
/// D10 round 4, 2026-11-20), and the batch that a reader error ends adds nothing (review M3,
/// 2026-11-20).
#[allow(clippy::too_many_arguments)]
fn scan_xml_titles(
    resolved: &[ResolvedOptions],
    query_records: &[fasta_reader::FastaRecord],
    query_end: &QueryEnd,
    lowercase: bool,
    query_range: Option<SequenceRange>,
    batch_size: usize,
    keep_whole: bool,
    subject_set: bool,
    xml_matrix: Option<crate::report::blast_xml::matrix::XmlMatrixFailure>,
    parse_deflines: bool,
    scope_refusal: Option<&crate::blastinput::query_batch::ScopeRefusal>,
) -> Result<super::blast_xml::QueryScan> {
    let structured: Vec<u8> = resolved
        .iter()
        .map(|options| options.outfmt)
        .filter(|outfmt| matches!(outfmt, 5 | 13..=16))
        .collect();
    let query_scopes = crate::blastinput::query_batch::QueryScopes::of_search(
        parse_deflines,
        subject_set,
        structured.contains(&5),
        "BLASTX",
    )
    .map(|scopes| {
        scopes
            .with_refusal(scope_refusal.cloned())
            .with_xml_matrix_failure(xml_matrix)
    });
    let count = crate::blastinput::query_batch::BatchCount::of(query_scopes.as_ref());
    // A matrix name that NCBI's BLAST XML report data does not know fails the report data before
    // any report is written (`xml_matrix`, the run's `blast_xml::xml_matrix_failure`): no title is
    // checked, and the queries are read for their scopes alone.
    if structured.is_empty() || batch_size == 0 || (xml_matrix.is_some() && query_scopes.is_none())
    {
        return Ok(super::blast_xml::QueryScan::default());
    }
    let holds_back =
        crate::report::blast_xml2::reports_held(subject_set, false) && !structured.contains(&5);
    super::blast_xml::scan_queries(
        &resolved[0],
        holds_back,
        xml_matrix.is_none(),
        query_scopes,
        &mut |deliver| {
            reader_batches(
                query_records,
                query_end,
                lowercase,
                query_range,
                batch_size,
                keep_whole,
                0,
                parse_deflines,
                count.clone(),
                deliver,
                &mut |_| Ok(()),
            )
        },
    )
    .map_err(report_error)
}

/// The subjects as the reports of a query batch show them: `shown`, or a database's sequences of
/// the batch.
fn batch_shown<'a>(
    database: Option<&DatabaseReports<'a>>,
    reported: Option<&'a ReportedSubjects>,
    shown: ShownSubjects<'a>,
) -> ShownSubjects<'a> {
    match (database, reported) {
        (Some(run), Some(reported)) => ShownSubjects::Database {
            name: run.database.name(),
            pairwise: run.pairwise,
            reported,
        },
        _ => shown,
    }
}

/// The reader batches of `query_records` as BLASTX's reader makes them (`input_source::BatchCutter`,
/// the command line's too), given to `deliver` in order: a record's messages (`warn`), then the
/// record. `query_range` (`-query_loc`) cuts each record to its interval (its whole letters kept
/// with `keep_whole`, the archive's); the skipped records are kept for the batch's scope
/// (`FastaRecord::skipped`) with `keep_skipped` (`input_source::QueryStream::keep_skipped`). Under
/// `-parse_deflines` (`count`) each record counts as NCBI's reader counts it
/// (`query_batch::BatchCount`). After the records comes `end`, as BLASTX's reader meets it after
/// the last record it read (`input_source::QueryStream::read_batches`): a reader error stops the
/// reading in the batch being read, after the messages of its records and of the record the
/// reader failed on, so that batch is never delivered (line WQA's D10, 2026-11-19); blank lines
/// alone, an input without a record, are a batch without a query (NCBI's `Empty
/// CBlastQueryVector`).
fn reader_batches(
    query_records: &[fasta_reader::FastaRecord],
    end: &QueryEnd,
    lowercase: bool,
    query_range: Option<SequenceRange>,
    batch_size: usize,
    keep_whole: bool,
    around: usize,
    keep_skipped: bool,
    count: Option<crate::blastinput::query_batch::BatchCount>,
    deliver: &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
    warn: &mut dyn FnMut(&[u8]) -> Result<()>,
) -> Result<()> {
    let mut cutter = BatchCutter::new(CutSettings {
        range: query_range.as_ref(),
        keep_whole,
        around,
        keep_skipped,
        batch_size,
        count,
        blank_batch: true,
    });
    for record in query_records {
        let query = from_shared_ref(record, lowercase)?;
        if !record.warnings.is_empty() {
            warn(&record.warnings)?;
        }
        cutter.add(query, deliver)?;
    }
    if let QueryEnd::Error(error) = end {
        if !error.warnings().is_empty() {
            warn(error.warnings())?;
        }
        return Err(error.app_error());
    }
    // An input without a line gives no batch; blank lines alone give one without a query.
    if query_records.is_empty() && !matches!(end, QueryEnd::BlankLines) {
        return Ok(());
    }
    cutter.finish(deliver)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:254-254
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:297-297
/// The prologs, the query batches (`run_query_batches`) with the reports of each in every format,
/// the epilogs and the hit list of a search of records.
fn report_batches(mut run: RunReports<'_>, outputs: &mut ReportOutputs<'_>) -> Result<()> {
    // The scopes that the read before the search made (`scan_xml_titles`), which the search's
    // batches then keep.
    let read_ahead_scopes = run.xml_scan.query_scopes.take();
    let run = &run;
    let &RunReports {
        resolved,
        subjects,
        shown,
        ref database,
        out,
        ..
    } = run;
    // Whether every format is the archive, which NCBI writes without the prefetch (`render`).
    let only_archive = outputs.only_archive();
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1450-1452
    // Whether the queries' warnings that a structured first format does not write go to the
    // caller's receivers of the later formats that write them, each format's to its own
    // (`ReportOutputs::later_warnings`).
    let keeps_later = outputs.keeps_later_query_warnings();
    let ReportOutputs {
        other_outfmts,
        formats,
        diagnostics,
        later_warnings,
        hits: hits_sink,
        observer,
        html_posts,
    } = outputs;
    let mut later = later_warnings.as_deref_mut().filter(|_| keeps_later);
    // The caller's receivers of the side files of outfmt 13 and 14.
    let side_sinks: Vec<_> = formats
        .iter()
        .map(|format| format.sink.side_files().cloned())
        .collect();
    // What a refusal of each format's checks does (`FormatRefusal`), beside its writer.
    let mut refusals: Vec<&mut crate::api::local_blast::FormatRefusal> = Vec::new();
    let writers = formats
        .iter_mut()
        .map(|format| {
            refusals.push(&mut format.refusal);
            crate::cli::NcbiOstream::sink(&mut format.sink)
        })
        .collect::<io::Result<Vec<_>>>()?;
    // A run whose XML documents could still refuse after their first byte (a `-subject` title, a
    // later batch's taxonomy database or subject with parsed Seq-ids, the database counts of the
    // reports: `blast_xml::documents_may_refuse`, the command line's decision too) holds every
    // format in order (`blast_xml::HeldFormats`), as the command line holds its document
    // (`blast_xml::HeldOutput`), so that the refusal follows no byte.
    let outfmts: Vec<i32> = resolved
        .iter()
        .map(|options| i32::from(options.outfmt))
        .collect();
    let xml_may_refuse = super::blast_xml::documents_may_refuse(
        &outfmts,
        database.as_ref().map(|run| run.database),
        subjects,
    );
    let hold = xml_may_refuse.then(super::blast_xml::HeldFormats::new);
    // The hold comes between a format's stream and the caller's writer, which has no stream of
    // NCBI's (`cli::NcbiOstream::held`): it keeps the writes and flushes in order.
    let mut writers: Vec<crate::cli::NcbiOstream<'_>> = match &hold {
        Some(hold) => {
            hold.hold();
            writers
                .into_iter()
                .map(|writer| crate::cli::NcbiOstream::held(hold.writer(writer)))
                .collect()
        }
        None => writers,
    };
    // The writer of each XML2 or JSON document (outfmt 13 to 16), which the run keeps across the
    // query batches (`native::run_inner`). The side files of 13 and 14 go next to `out` or to the
    // caller's receiver.
    let mut documents: Vec<_> = resolved
        .iter()
        .zip(&side_sinks)
        .map(|(options, sink)| {
            super::blast_xml::xml2_kind(options.outfmt).map(|kind| {
                let files = crate::report::blast_xml2::run_side_files(kind, out, sink.as_ref());
                // A run that holds its formats holds the side files of 13 and 14 too.
                let files = match &hold {
                    Some(hold) => hold.side_files(files),
                    None => files,
                };
                super::blast_xml::xml2_writer(kind, options, subjects, &shown, files)
            })
        })
        .collect();
    let wants_hits = hits_sink.is_some();
    let written = write_run(
        run,
        read_ahead_scopes,
        &mut writers,
        &mut refusals,
        &mut documents,
        wants_hits,
        only_archive,
        diagnostics,
        html_posts,
        &mut later,
        observer,
        other_outfmts,
    );
    // What the hold kept. A run that ends: every format, in order. A run that stops with an
    // error: nothing after one of the documents' refusals and the bytes of the batches before
    // after any other error (`HeldFormats::end_document_with`), as the command line. A write that
    // fails there is the run's error in place of the run's own, as without the hold. A format
    // that failed alone (`FormatRefusal`) gets nothing of what was held: its bytes are dropped
    // with the hold, so that its writer never had them.
    let refused = |sink: usize| refusals[sink].refused();
    // The observer's events of a held run come before its bytes reach their sinks, so a position
    // read from a sink at an event is not that of the event (`FormatObserver`): the same for every
    // hold of this function, a `-subject` run's too.
    match (&hold, &written) {
        (Some(hold), Ok(_)) => hold.release(&refused)?,
        (Some(hold), Err(error)) => hold.end_document_with(error, &refused)?,
        (None, _) => {}
    }
    let (run_hits, run_subjects) = written?;
    if let Some(sink) = hits_sink {
        sink(std::sync::Arc::new(run_hits));
        if let (true, Some(observer)) = (database.is_some(), observer.as_deref_mut()) {
            observer.database_subjects(run_subjects);
        }
    }
    Ok(())
}

/// The run's diagnostics as the writer of a query batch's shared warnings for one format
/// (`write_run`): it tells whether the format wrote any text there, which marks the batch as
/// warned whatever became of the format's report afterwards. No NCBI counterpart: NCBI's
/// command line writes one format. It counts any text, not only the queries' warnings: the
/// formats that share the diagnostics in LOSAT Web (outfmt 6 and 7 of the other search) run the
/// same checks and stop at the same row (review r3 of 2026-11-20, L1).
struct SharedWarnings<'a> {
    diagnostics: &'a mut dyn Write,
    wrote: bool,
}

impl Write for SharedWarnings<'_> {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        let written = self.diagnostics.write(bytes)?;
        self.wrote |= written > 0;
        Ok(written)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.diagnostics.flush()
    }
}

/// `report_batches` from the prologs to the epilogs, into `writers` (each format's, behind the
/// hold of a run that holds its formats) and `documents` (the writer of each XML2 or JSON
/// document): whatever ends it, `report_batches` then ends the hold. `refusals` say what a
/// refusal of each format's checks does (`FormatRefusal`): a format that fails alone writes
/// nothing after it, and the first format that is still written writes the warnings that the
/// formats share. `only_archive`: every format is the archive (`ReportOutputs::only_archive`).
/// `later`: by format, the receivers of the queries' warnings that the first format does not write
/// (a structured format): each later format that writes them writes its own there, as the
/// command line's run of it writes them to standard error (`ReportOutputs::later_warnings`).
/// `html_posts`: where the `-html` reports' posts go apart from `diagnostics`
/// (`ReportOutputs::html_posts`). Returns the run's hit list (`wants_hits`) and the database
/// sequences it names.
#[allow(clippy::too_many_arguments)]
fn write_run(
    run: &RunReports<'_>,
    read_ahead_scopes: Option<crate::blastinput::query_batch::QueryScopes>,
    writers: &mut [crate::cli::NcbiOstream<'_>],
    refusals: &mut [&mut crate::api::local_blast::FormatRefusal],
    documents: &mut [Option<super::blast_xml::RunXml2>],
    wants_hits: bool,
    only_archive: bool,
    diagnostics: &mut &mut (dyn Write + Send),
    html_posts: &mut Option<&mut (dyn Write + Send)>,
    later: &mut Option<&mut [Vec<u8>]>,
    observer: &mut Option<&mut (dyn crate::api::local_blast::FormatObserver + Send)>,
    other_outfmts: &[&str],
) -> Result<(
    Vec<crate::report::pairwise::PairwiseHit>,
    Vec<ReportedSubject>,
)> {
    let &RunReports {
        lowercase,
        query_records,
        query_end,
        query_range,
        batch_size,
        choices,
        resolved,
        subjects,
        shown,
        search_subjects,
        serial_batches,
        parse_deflines,
        scope_refusal,
        ref database,
        archive,
        strategy,
        search_error,
        ref xml_scan,
        xml_matrix,
        ..
    } = run;
    // The query whose batch the run's XML failure stops, of the scopes read ahead
    // (`QueryScopes::xml_failing_query`).
    let xml_matrix = xml_matrix.map(|failure| {
        crate::blastinput::query_batch::QueryScopes::xml_failing_query(
            read_ahead_scopes.as_ref(),
            failure,
        )
    });
    // Under `-parse_deflines`, the scopes of the batches' queries (`QueryScopes`): with outfmt 5 and
    // 13 to 16 the read before the search made them (`scan_xml_titles`, `QueryScopes::read_ahead`),
    // whose refusal comes where the loop takes its batch; else the batches add theirs where the loop
    // takes them (`native::add_batch_scope`), as the command line's (`native::run_inner`). Either
    // way a batch that a reader error ends adds none (review M3, 2026-11-20), and LOSAT's refusal of
    // two queries with one id stops the run at the batch of the second, after the reports of the
    // batches before (line WQA D10 round 4, 2026-11-20), in the form of every program (round 5).
    let xml = choices.iter().any(|choice| choice.number == 5);
    let query_scopes = std::cell::RefCell::new(read_ahead_scopes.or_else(|| {
        crate::blastinput::query_batch::QueryScopes::of_search(
            parse_deflines,
            database.is_none(),
            xml,
            "BLASTX",
        )
        .map(|scopes| scopes.with_refusal(scope_refusal.cloned()))
    }));
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:254-254
    for (index, (writer, options)) in writers.iter_mut().zip(resolved).enumerate() {
        let mut quiet = io::sink();
        // A later format's `-html` prolog posts its references' reads as a run of it alone
        // (`report::write_prolog`'s `html_posts`; line RG T4a2 fix round 2, review c C3).
        let (warnings, posts): (&mut dyn Write, Option<&mut dyn Write>) = if index == 0 {
            (
                &mut **diagnostics,
                html_posts
                    .as_deref_mut()
                    .map(|posts| posts as &mut dyn Write),
            )
        } else {
            let posts = crate::report::query_warnings::html_posts_sink(
                html_posts
                    .as_deref_mut()
                    .map(|posts| posts as &mut dyn Write),
                &mut **diagnostics,
            );
            (&mut quiet, Some(posts))
        };
        report::write_prolog(writer, warnings, posts, options, subjects, &shown)
            .map_err(report_error)?;
        if let Some(document) = &mut documents[index] {
            let mut start = Vec::new();
            document.document_mut().start(&mut start);
            writer.write_all(&start)?;
        }
        writer.flush()?;
    }
    // NCBI: a batch size of 0 makes the first batch empty, and its query factory refuses it.
    if batch_size == 0 {
        return Err(engine("Empty CBlastQueryVector"));
    }
    // The search's options are those of every format; the formats differ only in the counts of
    // the pairwise report (`run_report_counts`).
    let search_options = &resolved[0];
    // The Bioseq of a cut record in the archive and the search strategy's request holds its whole
    // letters (`FastaRecord::cut_in_place`), as the command line's reader keeps them
    // (`native::run_inner`).
    let keep_whole = archive.is_some() || strategy.is_some();
    // With `-subject`, the letters around each `-query_loc` interval that the Seq-loc mapper can
    // read into a query row (`empty_rows::query_letters_around`), as the command line's reader.
    let around = match database {
        Some(_) => 0,
        None => super::empty_rows::query_letters_around(search_options, subjects),
    };
    // The records that `-query_loc` skips, for NCBI's scope of a batch's prefetch, where a
    // prefetch of the run can stop, and under `-parse_deflines` for the scopes of the queries
    // (as the command line's reader, `native`).
    let keep_skipped = query_range.is_some()
        && (parse_deflines
            || database
                .as_ref()
                .is_some_and(|run| run.database.prefetch_can_stop()));
    let mut read = |deliver: &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
                    warn: &mut dyn FnMut(&[u8]) -> Result<()>|
     -> Result<()> {
        reader_batches(
            query_records,
            query_end,
            lowercase,
            query_range,
            batch_size,
            keep_whole,
            around,
            keep_skipped,
            crate::blastinput::query_batch::BatchCount::of_search(
                parse_deflines,
                database.is_none(),
                xml,
                "BLASTX",
            ),
            deliver,
            warn,
        )
    };
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:253-259
    // What the run writes before its query batches ends here.
    if let Some(observer) = observer.as_deref_mut() {
        observer.queries_begin();
    }
    let made_hits = wants_hits || choices.iter().any(|choice| choice.number == 0);
    let mut run_hits = Vec::new();
    // A database's sequences that the run's hit list names (`FormatObserver::database_subjects`),
    // and the place in them of each ordinal id: a batch names them by its own ranks.
    let mut run_subjects: Vec<ReportedSubject> = Vec::new();
    let mut subject_places = std::collections::HashMap::new();
    let mut first_hsp = 0;
    let asks = ReportedAsks::of_all(resolved, other_outfmts);
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:261-264 (cited at
    // `native::run_inner`'s `batch_take`)
    // The loop's takes of a reader batch, as the command line's (`native::BatchTake`): the batch
    // joins the scopes of the queries (line WQA D10 round 5, 2026-11-20) and, with LOSAT Web's
    // `-export_search_strategy` (the library refuses the option, `api::search::file_options`), its
    // request is written before the loop looks at the batch's search.
    let batch_take = super::native::BatchTake::new(
        &query_scopes,
        strategy,
        archive,
        search_options,
        search_subjects,
        search_error,
        database.is_none().then_some(subjects),
    );
    // NCBI's scope of the query records at each batch's prefetch (`native::BatchScopes`).
    let batch_scopes =
        super::native::BatchScopes::new(xml, database.as_ref().map(|run| run.database));
    let mut render = |diagnostics: &mut dyn Write,
                      queries: &[FastaRecord],
                      results: &mut Vec<BatchResults>,
                      queries_before: usize,
                      pool: Option<&crate::utils::threading::SearchPool<'_>>|
     -> Result<()> {
        // The batch's set-up has written its `EffectiveSearchSpace` option, which the next batch's
        // search strategy request holds (`take` below).
        batch_take.note_searched(queries, results, queries_before)?;
        // Every format and the hit list read the lists in the order of the reports.
        report::sort_lists(results, SubjectLoc::of_run(database.is_some()));
        // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:283-288 (cited at
        // `native::run_inner`'s batches)
        // A database: the titles of the batch's reported sequences, read in NCBI's traceback
        // (`native::batch_titles`). Then the prefetch before the batch's reports reads
        // `PRE_FETCH_SEQS_LIMIT` (not with the archive alone): a value that NCBI cannot convert
        // stops the first batch here (`app::check_prefetch_environment`).
        let titles = match database {
            Some(run) => super::native::batch_titles(
                run.database,
                search_subjects,
                results,
                search_options.gapped,
                asks,
                (&mut *diagnostics, archive),
            )?,
            None => None,
        };
        app::check_prefetch_environment(only_archive)?;
        // A database: the batch's reported sequences, and its hit lists named by their ranks.
        let scope = || batch_scopes.of_batch(queries);
        let mut reported = match database {
            Some(run) => Some(reported_batch(
                run.database,
                search_subjects,
                titles,
                search_options,
                asks,
                results,
                (queries, queries_before),
                query_scopes
                    .borrow()
                    .as_ref()
                    .map(|scopes| scopes.numbered_from(queries_before)),
                &crate::algorithm::common::database::PrefetchQueries::new(0, &scope),
                &mut *diagnostics,
            )?),
            None => None,
        };
        batch_scopes.reported(queries);
        let hits = match made_hits {
            true => Some(
                report::pairwise_hits(
                    search_options,
                    queries,
                    reported
                        .as_ref()
                        .map_or(subjects, ReportedSubjects::records),
                    &batch_shown(database.as_ref(), reported.as_ref(), shown),
                    results,
                    None,
                )
                .map_err(report_error)?,
            ),
            false => None,
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:82-89
        // The statistics of each query of the batch for the observer.
        if let Some(observer) = observer.as_deref_mut() {
            super::blast_xml::observe_query_stats(
                search_options,
                results,
                queries_before,
                &mut *observer,
            );
        }
        // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:289-290
        // Whether a format of the batch wrote the warnings that the formats share. This is
        // LOSAT's own: NCBI's command line writes one format. The first format that is written
        // writes them, and the diagnostics hold that format's standard error of the batch. The
        // batch is warned once a format wrote its report of the batch whole or wrote any text to
        // the diagnostics, also where its report then stopped (it failed alone after a query's
        // warning: the formats after it do not write that warning again). A format that its
        // checks refused before its first text there leaves the warnings to the next
        // (`SharedWarnings`; line WQA, 2026-11-20: the flag was read off the refusal, so a format
        // that stopped after a warning left the next one to write the same warning again).
        let mut warned = false;
        for (index, ((writer, options), refusal)) in writers
            .iter_mut()
            .zip(resolved)
            .zip(refusals.iter_mut())
            .enumerate()
        {
            // A format that failed alone writes nothing more.
            if refusal.refused() {
                continue;
            }
            // The warnings of the batch's queries are written once, with the first format that
            // is written; NCBI's posts of the bytes that a Seq-annot's stream replaces (outfmt 8,
            // 9; the only text that their writer writes to it, `seqalign::write_blastx_seqalign`)
            // at every format, as the other programs write them.
            let mut quiet = io::sink();
            let report_format = app::ReportFormat::from_number(i32::from(options.outfmt));
            let seq_annot = report_format
                .and_then(crate::report::serial::SerialFormat::of)
                .is_some();
            let writes = report_format.is_none_or(app::ReportFormat::writes_query_warnings);
            let warns = !warned;
            // A later format that writes the queries' warnings where the first, a structured
            // one, does not, writes its standard error of the batch's reports to its own
            // receiver (`later`): one whose report stops has given what it wrote until then, and
            // the caller keeps each line of the formats' texts once (line WQA's sweep, D3 review
            // r3, 2026-11-19).
            let own_later = later
                .as_deref_mut()
                .and_then(|later| later.get_mut(index))
                .filter(|_| writes);
            let mut shared = SharedWarnings {
                diagnostics: &mut *diagnostics,
                wrote: false,
            };
            // A `-html` report's posts of its reads of NCBI's formatter registries are what a run
            // of its format alone posts, at any place among the formats, as in the other four
            // programs (`query_warnings::report_posts`): the diagnostics also where a later
            // format's standard error goes to its receiver or nowhere (line RG T4a2 fix round 2,
            // review c C3 = review d D5), or the call's receiver of them
            // (`ReportOutputs::html_posts`).
            let mut caller_posts = html_posts
                .as_deref_mut()
                .map(|posts| posts as &mut dyn Write);
            let (warnings, report_posts): (&mut dyn Write, Option<&mut dyn Write>) =
                match (warns || seq_annot, own_later) {
                    (true, _) => (&mut shared, caller_posts.take()),
                    (false, Some(later)) => {
                        (later, Some(caller_posts.take().unwrap_or(&mut shared)))
                    }
                    _ => (&mut quiet, Some(caller_posts.take().unwrap_or(&mut shared))),
                };
            // The values of a custom tabular field list are each format's (`reported_batch`
            // read the first format's).
            if let (Some(run), Some(reported), Some(spec), true) =
                (database, reported.as_mut(), options.tabular, index > 0)
            {
                reported.read_tabular(run.database, spec, results, options.gapped)?;
            }
            // The BLAST archive (outfmt 11): one document per query batch, written from the
            // data of outfmt 8's Seq-aligns.
            if let (Some(handle), 11) = (archive, options.outfmt) {
                let (subject_records, names, call_subjects): (_, _, &[FastaRecord]) =
                    match (database, reported.as_ref()) {
                        (Some(run), Some(reported)) => (
                            None,
                            super::seqalign::SeqalignSubjects::Database {
                                reported,
                                title: &run.pairwise.definition,
                            },
                            reported.records(),
                        ),
                        _ => (
                            Some(subjects),
                            super::seqalign::SeqalignSubjects::Records,
                            subjects,
                        ),
                    };
                let archived = archive::write_blastx_archive(
                    &archive::ArchiveCall {
                        handle,
                        batch: archive::CallBatch {
                            options,
                            queries,
                            requests: &archive::requests(subject_records),
                            batches: results,
                        },
                        subjects: call_subjects,
                        names,
                    },
                    writer,
                    &mut shared,
                )
                .map_err(archive_error)
                .and_then(|()| writer.flush().map_err(Into::into));
                // An error of the format's report of the batch fails the call, or this format
                // alone (`FormatRefusal`), which writes nothing in the later batches. The
                // archive, as the first format, takes the place of the warnings' writer.
                let written = refusal.check(archived)?;
                warned |= warns && (written || shared.wrote);
                continue;
            }
            let mut probe = observer
                .as_deref_mut()
                .map(|observer| FormatProbe::new(observer, index));
            let rendered = report::render_queries_with(
                writer,
                warnings,
                options,
                queries,
                reported
                    .as_ref()
                    .map_or(subjects, ReportedSubjects::records),
                &batch_shown(database.as_ref(), reported.as_ref(), shown),
                results,
                queries_before,
                xml_scan.stops_after == Some(queries_before + queries.len()),
                documents[index].as_mut(),
                xml_matrix,
                query_scopes
                    .borrow()
                    .as_ref()
                    .map(|scopes| scopes.numbered_from(queries_before)),
                pool,
                Some(LocalReport {
                    probe: probe.as_mut(),
                    first_hsp,
                    hits: hits.as_deref(),
                    refusal,
                    report: index,
                    html_posts: report_posts,
                }),
            )
            .map_err(report_error)
            .and_then(|()| writer.flush().map_err(Into::into));
            // An error of the format's report of the batch fails the call, or this format alone
            // (`FormatRefusal`), which writes nothing in the later batches.
            let written = refusal.check(rendered)?;
            // A format refused by its checks wrote no warning (they come before its warnings and
            // rows, `report::render_queries_with`); one whose report stopped after a warning
            // wrote it.
            warned |= warns && (written || shared.wrote);
        }
        first_hsp += match &hits {
            Some(hits) => hits.len(),
            None => report::reported_hsps(results, search_options.gapped),
        };
        // The run's hit list numbers the queries of the run, and a database's sequences in the
        // order that it first names them.
        if wants_hits {
            for mut hit in hits.into_iter().flatten() {
                hit.hit.q_idx = u32::try_from(queries_before + hit.hit.q_idx as usize)?;
                if let (Some(run), Some(reported)) = (database, reported.as_ref()) {
                    let rank = hit.hit.s_idx as usize;
                    let ordinal = reported.ordinal(rank);
                    let place = *subject_places.entry(ordinal).or_insert_with(|| {
                        run_subjects.push(library_subject(
                            run.database.x_seqdb(),
                            reported.stored_titles()[rank].as_ref(),
                            crate::report::pairwise::stored_source(reported.sources(), rank as u32),
                            ordinal,
                            reported.parsed_subject(rank),
                            reported.records()[rank].sequence.len(),
                            true,
                            "BLASTX",
                        ));
                        run_subjects.len() - 1
                    });
                    hit.hit.s_idx = u32::try_from(place)?;
                }
                run_hits.push(hit);
            }
        }
        Ok(())
    };
    let mut take = |diagnostics: &mut dyn Write,
                    queries: &[FastaRecord],
                    queries_before: usize|
     -> Result<()> { batch_take.take(diagnostics, queries, queries_before) };
    let processed = run_query_batches(
        search_options,
        search_subjects,
        serial_batches,
        &mut read,
        &mut **diagnostics,
        search_error,
        &mut render,
        &mut take,
    )?;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:297-297
    // The first format that is still written writes the warnings of the epilogs.
    let first = refusals.iter().position(|refusal| !refusal.refused());
    for (index, ((writer, options), refusal)) in writers
        .iter_mut()
        .zip(resolved)
        .zip(refusals.iter_mut())
        .enumerate()
    {
        if refusal.refused() {
            continue;
        }
        let mut quiet = io::sink();
        let warnings: &mut dyn Write = if Some(index) == first {
            &mut **diagnostics
        } else {
            &mut quiet
        };
        let ended = report::write_epilog(
            writer,
            warnings,
            options,
            subjects,
            &shown,
            processed,
            documents[index].as_mut(),
            xml_matrix,
        )
        .map_err(report_error)
        .and_then(|()| writer.flush().map_err(Into::into));
        // An error of the format's epilog fails the call, or this format alone (`FormatRefusal`).
        refusal.check(ended)?;
    }
    Ok((run_hits, run_subjects))
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:195-200
/// Runs one BLASTX search of the records `query_records` against the database of `args` and
/// writes every requested output format from the same result (the library API's entry; the
/// command line's is `native::run`). `settings` are the NCBI application settings, which find the
/// database. The steps are `native::run_inner`'s, in its order; what differs is what they read:
/// the records are read already, and the outputs are the caller's: where the command line holds
/// its output for a refusal that a later batch can still make, every format of the call is held
/// in order (`report_batches`). An archive (outfmt 11) ends as the command line's does
/// (`ReportOutputs::finish_archive`).
pub fn run_local_database(
    mut args: BlastxArgs,
    query_records: &[fasta_reader::FastaRecord],
    query_end: &QueryEnd,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // NCBI's options factory checks the task first (`TaskParser::check`); the command line's
    // parser has checked it, a caller of `run_local` has not.
    super::args::BLASTX_TASK.check(&mut args.task)?;
    args.matrix_search = Some(settings.matrix_search.clone());
    let requested = requested_formats(outputs);
    let (choices, resolutions) = read_format_options(&args, &requested)?;
    let archive_requested = choices.iter().any(|choice| choice.number() == Some(11));
    // The archive shares the options handle of the run's search strategy (LOSAT Web's
    // `-export_search_strategy`, `web_adapter::give_strategy`), as on the command line
    // (`native::RunSetup::open`).
    let archive = archive::handle_for(
        &args,
        archive_requested,
        &crate::report::archive::strategy::options_of(args.strategy.as_ref()),
    )
    .map(|handle| outputs.with_archive_outfmt(handle));
    crate::utils::threading::validate_threads(args.num_threads as usize)?;
    let no_query = query_records.is_empty() && query_end.is_input();
    outputs.run_records(
        args.strategy.as_ref(),
        archive.as_ref(),
        no_query,
        |outputs| {
            search_database(
                &args,
                query_records,
                query_end,
                settings,
                (choices, resolutions),
                archive.as_ref(),
                outputs,
            )
        },
    )
}

/// `run_local_database` from the taxid or seqid list on.
fn search_database(
    args: &BlastxArgs,
    query_records: &[fasta_reader::FastaRecord],
    query_end: &QueryEnd,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
    (choices, resolutions): ReadOptions,
    archive: Option<&crate::report::archive::ArchiveHandle>,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    let name = args
        .database_name()
        .ok_or_else(app::missing_subject_error)?;
    // The handler of the database arguments reads the list before the other handlers
    // (`native::run_inner`).
    let list = args
        .db_filter
        .resolve(&settings, args.app_name(), &mut *outputs.diagnostics)?;
    // `-out` (`run_to_file`) is made where the command line opens it (`native::run`).
    outputs.open_out()?;
    check_filtering(args)?;
    let query_range = parse_query_range(args)?;
    let (choices, mut resolved) =
        check_resolutions(args, &choices, resolutions, &mut *outputs.diagnostics)?;
    // The database opens after the options are processed and before the query is read, as NCBI's
    // `InitializeSubject`.
    let (keeps_scope, separate_xml2) = crate::report::db_subject::FormatThrow::scope_flags(
        outputs.formats.first().map(|format| format.outfmt),
    );
    let database = BlastxDatabase::open_restricted(
        &name,
        &settings,
        list.as_ref(),
        args.db_mask.option().as_ref(),
        &mut *outputs.diagnostics,
        super::args::ProgramName::BLASTX.upper,
        args.app_name(),
    )?
    .showing_gis(args.show_gis)
    .after_calls(args.first_query_ordinal > 0)
    .keeping_scope(keeps_scope, separate_xml2);
    crate::algorithm::common::database::note_archive_database_filter(
        archive,
        args.strategy.as_ref(),
        &args.db_filter,
        database.list_gis(),
        database.list_taxids(),
        &database.statistics(),
        false,
        args.searchsp.is_some_and(|space| space != 0),
        args.dbsize.is_some_and(|size| size > 0),
        "BLASTX",
    )?;
    crate::algorithm::common::database::check_formats_report_sums(
        || database.pairwise_database(),
        outputs
            .formats
            .iter_mut()
            .map(|format| (format.outfmt, &mut format.refusal)),
        "BLASTX",
    )?;
    for options in &mut resolved {
        options.query_range = query_range;
        options.filtered_empty = database.filtered_empty();
        if options
            .tabular
            .is_some_and(|spec| spec.looks_up_every_taxid())
        {
            options.taxdb_found = database.taxdb_found();
        }
    }
    // LOSAT's refusals of a database for the Seq-annot formats and of a name that BLAST XML cannot
    // write, before the search.
    for choice in &choices {
        if let Some(serial) = app::ReportFormat::from_number(choice.number)
            .and_then(crate::report::serial::SerialFormat::of)
        {
            crate::report::serial::seqalign::check_database(
                serial,
                database.name(),
                Program::Blastx,
            )?;
        }
    }
    if query_records.is_empty() && query_end.is_input() {
        outputs
            .diagnostics
            .write_all(format!("Warning: [{}] Query is Empty!\n", args.app_name()).as_bytes())?;
        return Ok(());
    }
    // The letters of a reader batch (`native::run_inner`).
    let batch_size = app::query_batch_size(BATCH_SIZE as u32)? as usize;
    // The formatter resolves the database masking before its other work (`native::run_inner`).
    if let Some(option) = args.db_mask.option() {
        let mask = database.resolve_mask_option(
            &option,
            args.app_name(),
            archive,
            args.strategy.as_ref(),
            &mut *outputs.diagnostics,
        )?;
        if let Some(mask) = &mask {
            // Once per word and window of the formats' options (they are the same search's).
            let mut checked: Vec<(i32, i32)> = Vec::new();
            for options in &resolved {
                let scan = (options.word_size, options.window_size);
                if !checked.contains(&scan) {
                    database.check_mask_limits(mask, options)?;
                    checked.push(scan);
                }
            }
        }
        for options in &mut resolved {
            options.db_mask = mask.clone();
        }
    }
    // NCBI's error of an XInclude format without a file name, where the command line gives it
    // (`native::run_inner`, blastx_app.cpp:251).
    outputs.xinclude_check(args.out.as_deref())?;
    app::check_search_environment("BLASTX", false)?;
    let pairwise = database.pairwise_database_masked(
        resolved
            .first()
            .and_then(|options| options.db_mask.as_ref()),
    );
    let no_reported = ReportedSubjects::default();
    let shown = ShownSubjects::Database {
        name: database.name(),
        pairwise: &pairwise,
        reported: &no_reported,
    };
    // The database's sequences, read in blocks where they are searched; the empty database's error
    // belongs to the first batch with a valid context (`native::run_inner`).
    let allow_long = resolved
        .first()
        .and_then(|options| options.db_mask.as_ref())
        .is_some_and(|mask| mask.is_hard());
    let (longest, unreadable) = match database.x_open_stream_with_long(allow_long) {
        Ok(longest) => (Some(longest), None),
        Err(error) => (None, Some(error)),
    };
    let stream = longest.map(|longest| {
        super::x_db_stream::XStreamDb::new_with_long(database.x_seqdb(), longest, allow_long)
    });
    let unreadable = std::cell::RefCell::new(unreadable);
    let search_subjects = match &stream {
        Some(stream) => Subjects::Stream(super::x_db_stream::XStream {
            db: stream,
            block: None,
        }),
        None => Subjects::Unreadable,
    };
    let search_error = |error: anyhow::Error| -> anyhow::Error {
        database_search_error(error, stream.as_ref(), &unreadable)
    };
    let xml_matrix = super::blast_xml::xml_matrix_failure(&resolved, false);
    let xml_scan = scan_xml_titles(
        &resolved,
        query_records,
        query_end,
        args.lcase_masking,
        query_range,
        batch_size,
        archive.is_some() || args.strategy.is_some(),
        false,
        xml_matrix,
        args.parse_deflines,
        args.scope_refusal.as_ref(),
    )?;
    report_batches(
        RunReports {
            lowercase: args.lcase_masking,
            query_records,
            query_end,
            query_range,
            batch_size,
            choices: &choices,
            resolved: &resolved,
            subjects: &[],
            shown,
            search_subjects,
            serial_batches: args.parse_deflines,
            parse_deflines: args.parse_deflines,
            scope_refusal: args.scope_refusal.as_ref(),
            database: Some(DatabaseReports {
                database: &database,
                pairwise: &pairwise,
            }),
            archive,
            // The search strategy file of LOSAT Web's run (`web_adapter::give_strategy`); the
            // library writes none (`api::search::file_options`).
            strategy: args.strategy.as_ref(),
            out: args.out.as_deref(),
            search_error: &search_error,
            xml_scan,
            xml_matrix,
        },
        outputs,
    )
}

/// The checks of the options alone, without the records (the `validate` of web ABI v2), for the
/// output formats `outfmts` of one `run_local` call: `run_local`'s refusals and LOSAT's, the
/// thread count, the syntax of `-subject_loc` and `-query_loc`, then NCBI's errors of the
/// options, with the messages of a run. A run reports NCBI's errors after the subjects are read,
/// so for an argv and records that both fail the run reports the records' error.
///
/// Outfmt 13 and 14 are read as in a run that gives them the caller's receiver of their side files
/// (`OutputSink::Files`, which LOSAT Web gives them): a run without one gives NCBI's error after the
/// records are read (`ReportOutputs::xinclude_check`).
pub fn check_options(args: &BlastxArgs, outfmts: &[&str]) -> Result<()> {
    let outfmt_bytes: Vec<&[u8]> = outfmts.iter().map(|s| s.as_bytes()).collect();
    let (choices, resolutions) = read_options(args, &outfmt_bytes)?;
    crate::utils::threading::validate_threads(args.num_threads as usize)?;
    parse_subject_range(args)?;
    check_filtering(args)?;
    parse_query_range(args)?;
    check_resolutions(args, &choices, resolutions, &mut io::sink()).map(|_| ())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2540-2545
/// `-subject_loc`, which the command line reads after it opened the subject file and before it
/// reads the subjects (`native::run_inner`).
fn parse_subject_range(args: &BlastxArgs) -> Result<Option<SequenceRange>> {
    parse_optional_range(args.subject_loc.as_deref(), RangeRole::Subject)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1995-1999
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:87-88,106-107
/// The filtering handler (`CFilteringArgs`), which comes after the subjects are read and before
/// the query options and formatting handlers: NCBI's error of `-seg` (`BlastxArgs::seg_spec`;
/// `resolve_with_format` keeps it with the other deferred errors).
pub(super) fn check_filtering(args: &BlastxArgs) -> Result<()> {
    match args.seg_spec() {
        Err(message) => Err(options_error(SemanticOptionsError::from(message).into())),
        Ok(_) => Ok(()),
    }
}

/// `-query_loc`, which the command line reads after the subjects, with the query options and
/// before the formatting handler (`native::run_inner`).
fn parse_query_range(args: &BlastxArgs) -> Result<Option<SequenceRange>> {
    // An imported strategy's query range stands for `-query_loc`'s (`import::apply`).
    match args.import.as_ref().and_then(|import| import.query_range) {
        Some(range) => Ok(Some(range)),
        None => parse_optional_range(args.query_loc.as_deref(), RangeRole::Query),
    }
}

/// The output formats `outfmts` as `run_local` reads them and the options of each, before the
/// records are read: LOSAT Web's refusals and LOSAT's fail here; NCBI's errors of the options
/// (`SemanticOptionsError`) stay in the resolutions, which the command line reports after the
/// subjects are read (`native::run_inner`).
fn read_options(args: &BlastxArgs, outfmts: &[&[u8]]) -> Result<ReadOptions> {
    // `-parse_deflines` is the records' (`run_local`): LOSAT Web reads them with the ids of their
    // deflines since line WQA's sweep, D7, 2026-11-19 (`register_inputs`), so its refusal of
    // records that kept none, which also refused two empty inputs that the command line ends with
    // NCBI's `Empty CBlastQueryVector`, is gone.
    read_format_options(args, outfmts)
}

/// The `-outfmt` value of each output of `outputs`.
fn requested_formats<'o>(outputs: &'o ReportOutputs<'_>) -> Vec<&'o [u8]> {
    outputs.formats.iter().map(|format| format.outfmt).collect()
}

/// The output formats of a search of records as the command line's early reading reads them
/// (`app::EarlyFormat`) and the options of each (`read_format_options`), whose NCBI errors
/// `check_resolutions` reports where the command line does.
type ReadOptions = (Vec<app::EarlyFormat>, Vec<Result<ResolvedOptions>>);

/// `read_options` without its refusals of what a search of records in memory cannot do. Outfmt 13
/// and 14 whose side files have no file name to go by get NCBI's error later, where the command
/// line gives it (`ReportOutputs::xinclude_check`).
fn read_format_options(args: &BlastxArgs, outfmts: &[&[u8]]) -> Result<ReadOptions> {
    // NCBI parses -outfmt before its option handlers read the subjects (`native::run_inner`);
    // where the argv names a search strategy file a value that does not parse is left to the
    // formatting handler (`app::early_format_choice`, `check_resolutions`).
    let mut choices = Vec::with_capacity(outfmts.len());
    for &outfmt in outfmts {
        let choice = app::early_format_choice(outfmt, args.export_search_strategy.is_some())?;
        choice.report_format(Program::Blastx, args.out.as_deref())?;
        choices.push(choice);
    }
    ensure!(!choices.is_empty(), "no output format was requested");
    // The options of each format: LOSAT's refusals now, NCBI's errors of the options after the
    // subjects are read (`native::run_inner`). A value left to the formatting handler has none
    // (its error comes there first), as the command line's (`native::RunSetup::open`).
    let mut resolutions = Vec::with_capacity(choices.len());
    for choice in &choices {
        let Some(choice) = choice.parsed() else {
            resolutions.push(Err(SemanticOptionsError::from(String::new()).into()));
            continue;
        };
        match args.resolve_with_format(choice) {
            Err(error) if error.downcast_ref::<SemanticOptionsError>().is_none() => {
                return Err(error)
            }
            // The file of an unknown matrix, looked up where the command line looks it up
            // (`native::run`), with the settings of the library's start.
            Ok(mut options) => {
                if let Some(search) = &args.matrix_search {
                    options.find_matrix_file(search)?;
                }
                resolutions.push(Ok(options));
            }
            resolution => resolutions.push(resolution),
        }
    }
    Ok((choices, resolutions))
}

/// The checks of the options that the command line makes after the subjects are read, in its
/// order: the formatting handler's (a value that the early reading left to it, which LOSAT Web
/// gives as the run's one format and as the arguments' `-outfmt`, is read again and fails, as the
/// command line's `native::run_inner`), the hit list size (its warnings and those of an ignored
/// -num_descriptions, the first format's, to `diagnostics`; formats that need different hit lists
/// are refused, `run_report_counts`), then NCBI's errors of each format's options. The formats'
/// choices and their options.
fn check_resolutions(
    args: &BlastxArgs,
    early: &[app::EarlyFormat],
    resolutions: Vec<Result<ResolvedOptions>>,
    diagnostics: &mut dyn Write,
) -> Result<(Vec<app::FormatChoice>, Vec<ResolvedOptions>)> {
    let mut choices = Vec::with_capacity(early.len());
    let outfmt_bytes = crate::utils::os_bytes::os_str_bytes(&args.outfmt);
    for choice in early {
        let choice = choice.at_formatting_handler(&outfmt_bytes)?;
        app::formatting_handler_check_of_run(&choice, args.app(), Program::Blastx)?;
        choices.push(choice);
    }
    let outfmts: Vec<i32> = choices.iter().map(|choice| choice.number).collect();
    let html = app::html_links(
        args.html,
        args.app_request_kind(),
        args.request_kind(),
        args.database_name().as_deref(),
        &outfmts,
    )?;
    let (_, warnings) = crate::blastinput::report_counts::run_report_counts_with(
        args.app_name(),
        &outfmts,
        super::args::max_target_seqs_count(args.max_target_seqs),
        args.num_descriptions,
        args.num_alignments,
        &args.formatting,
    )?;
    diagnostics.write_all(&warnings)?;
    // An imported strategy's warnings for the options it replaces
    // (`x_IssueWarningsForIgnoredOptions`, after the formatting handler's).
    if let Some(import) = &args.import {
        diagnostics.write_all(&import.option_warning_lines())?;
    }
    let resolved = resolutions
        .into_iter()
        .map(|resolution| {
            let mut options = resolution.map_err(options_error)?;
            options.html = html.clone();
            Ok(options)
        })
        .collect::<Result<Vec<_>>>()?;
    Ok((choices, resolved))
}
