//! BLASTX over records in memory: one search over the records that LOSAT Web registered or the
//! library API read (`run_local`), or of the library's records against a database
//! (`run_local_database`), written in every requested output format. The query batches are
//! searched and reported by the command line's loop (`native::run_query_batches`).

use super::{
    archive,
    args::{BlastxArgs, ResolvedOptions, SemanticOptionsError},
    database::{BlastxDatabase, ReportedSubjects},
    input::FastaRecord,
    native::{
        archive_error, database_search_error, engine, input_error, options_error, output_error,
        report_error, reported_batch, run_query_batches, write_empty_subject_warnings,
        ReportedAsks,
    },
    query_setup::BATCH_SIZE,
    report::{self, LocalReport, RecordsParsed, ShownSubjects},
    runtime::BatchResults,
    subjects::Subjects,
};
use crate::algorithm::common::database::library_subject;
use crate::api::local_blast::{FormatProbe, ReportOutputs, ReportedSubject};
use crate::blastinput::app::{self, Program};
use crate::blastinput::fasta_reader;
use crate::blastinput::seq_range::{
    parse_optional_range, record_cut, record_interval, RangeRole, RecordInterval, SequenceRange,
};
use anyhow::{bail, ensure, Context, Result};
use std::io::{self, Write};

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:429-467
// ```c++
// CRef<CSearchResultSet> results = m_RmtBlast->GetResultSet();
// formatter.PrintProlog();
// ...
// ITERATE(CSearchResultSet, result, *results) {
//     ...
//         formatter.PrintOneResultSet(**result, queries);
//     ...
// }
// ```
// NCBI formats one result set without searching again; several requested formats are several
// CBlastFormat printers over the same result set.
/// Runs one BLASTX search over already parsed records and writes every requested output format
/// from the same result (the entry of web ABI v2; the command line is `native::run`), and gives
/// the caller the final hit list (`ReportOutputs::hits`).
///
/// Each requested `-outfmt` is read as on the command line, and NCBI's checks of the options and
/// of the subjects come in the command line's order, before the search. The `-query` and
/// `-subject` values of `args` are used only as display names. Each format is the command line's
/// writer of it (`report::render_queries`); the side files of outfmt 13 and 14, `-db`, and
/// `-parse_deflines` for records without their ids are refused. LOSAT Web asks for outfmt 0, 6 and
/// 7, 6 and 7 with a custom field list or delimiter too. `-subject_loc` and `-query_loc` cut the records as the command line cuts them: each
/// subject as it is read, each query as its batch is read (`native::run_inner`).
///
/// The records are those of NCBI's reader (`fasta_reader`): the queries as its nucleotide reader
/// reads them, the subjects as its protein reader reads them (`blastx_record`). Their messages are
/// written where BLASTX's reader writes them: a subject's as the subjects are read, a query's
/// before the report of the batch that holds it. Queries without a record are NCBI's
/// `Query is Empty!`.
pub fn run_local(
    mut args: BlastxArgs,
    query_records: &[fasta_reader::FastaRecord],
    subject_records: &[fasta_reader::FastaRecord],
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // NCBI's options factory checks the task first (`TaskParser::check`); the command line's
    // parser has checked it, a caller of `run_local` has not.
    super::args::BLASTX_TASK.check(&mut args.task)?;
    let requested = requested_formats(&args, outputs);
    // LOSAT Web registers records before it knows the options: they keep no ids of their
    // deflines. The library API's `read_inputs` reads them.
    let ids_read = query_records
        .iter()
        .chain(subject_records)
        .any(|record| record.believed().is_some());
    let (choices, resolutions) = read_options(&args, &requested, ids_read)?;
    crate::utils::threading::validate_threads(args.num_threads as usize)?;
    let archive_requested = choices.iter().any(|choice| choice.number == 11);
    let archive = archive::handle_for(
        &args,
        archive_requested,
        &crate::report::archive::OptionsHandle::default(),
    )?
    .map(|handle| outputs.with_archive_outfmt(handle));
    match archive {
        Some(handle) => {
            let result = outputs.with_archive_diagnostics(handle.clone(), |view| {
                search_records(
                    &args,
                    query_records,
                    subject_records,
                    (choices, resolutions),
                    Some(&handle),
                    view,
                )
            });
            outputs.finish_archive(&handle, query_records.is_empty(), result)
        }
        None => search_records(
            &args,
            query_records,
            subject_records,
            (choices, resolutions),
            None,
            outputs,
        ),
    }
}

/// `run_local` from the reading of the subjects on.
fn search_records(
    args: &BlastxArgs,
    query_records: &[fasta_reader::FastaRecord],
    subject_records: &[fasta_reader::FastaRecord],
    (choices, resolutions): (Vec<app::FormatChoice>, Vec<Result<ResolvedOptions>>),
    archive: Option<&crate::report::archive::ArchiveHandle>,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    let Some(subject_path) = args.subject.as_deref() else {
        return Err(app::missing_subject_error());
    };
    let diagnostics = &mut outputs.diagnostics;
    let subject_range = parse_subject_range(args)?;
    // NCBI reads the subjects before the queries, each record's messages as it is read, and cuts
    // each to its `-subject_loc` interval: a record whose interval starts more than one letter
    // past its end stops the reading after its messages (`native::run_inner`).
    let mut subjects = Vec::with_capacity(subject_records.len());
    for record in subject_records {
        let mut subject = blastx_record(record, args.lcase_masking).map_err(input_error)?;
        diagnostics.write_all(&record.warnings)?;
        if let Some(range) = &subject_range {
            let (from, to_exclusive) = record_cut(range, subject.sequence.len())?;
            // The lower case of the whole record, which outfmt 0's Sbjct rows and the hits show
            // (`x_soft_subject::x_masked_subject_rows`), and the archive's Bioseqs.
            let keep_whole =
                archive.is_some() || (args.lcase_masking && !subject.lowercase_masks.is_empty());
            subject.cut_in_place(from, to_exclusive, keep_whole);
        }
        subjects.push(subject);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
    // ```c++
    // CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
    //     : m_QueryVector(& queries)
    // {
    //     if (queries.Empty()) {
    //         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
    //     }
    // ```
    if subjects.is_empty() {
        return Err(engine("Empty CBlastQueryVector"));
    }
    check_filtering(args)?;
    let query_range = parse_query_range(args)?;
    let mut resolved = check_resolutions(args, &choices, resolutions, &mut **diagnostics)?;
    for options in &mut resolved {
        options.query_range = query_range;
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:213-217
    // ```c++
    //         if(IsIStreamEmpty(m_CmdLineArgs->GetInputStream())){
    //            	ERR_POST(Warning << "Query is Empty!");
    //            	return BLAST_EXIT_SUCCESS;
    //         }
    // ```
    if query_records.is_empty() {
        diagnostics.write_all(b"Warning: [blastx] Query is Empty!\n")?;
        return Ok(());
    }
    // The letters of a reader batch (`native::run_inner`).
    let batch_size = app::query_batch_size("BLASTX", BATCH_SIZE as u32)? as usize;
    let xml_scan = scan_xml_titles(
        &resolved,
        query_records,
        args.lcase_masking,
        query_range,
        batch_size,
        archive.is_some(),
        true,
    )?;
    write_empty_subject_warnings(&subjects, &mut **diagnostics)?;
    // The environment of a search of subject records, after the subjects' warnings
    // (`native::run_inner`).
    app::check_search_environment("BLASTX", true, outputs.only_archive())?;
    let records_parsed = if args.parse_deflines {
        let pairwise = resolved.iter().any(|options| options.outfmt == 0);
        RecordsParsed::of(&subjects, args.show_gis, pairwise)?
    } else {
        RecordsParsed::default()
    };
    let shown = ShownSubjects::Fasta {
        path: subject_path,
        parsed: &records_parsed,
    };
    // The errors of a `-subject` search (`native::run_inner`).
    let search_error = |error: anyhow::Error| -> anyhow::Error {
        // LOSAT's refusal, not an error of NCBI's engine.
        if error.is::<super::parameters::SearchSpaceOverflow>()
            || error.is::<super::split::QuerySplitRefusal>()
        {
            error
        } else {
            engine(error)
        }
    };
    report_batches(
        &RunReports {
            lowercase: args.lcase_masking,
            query_records,
            query_range,
            batch_size,
            choices: &choices,
            resolved: &resolved,
            subjects: &subjects,
            shown,
            search_subjects: Subjects::Fasta(&subjects),
            serial_batches: args.parse_deflines,
            parse_deflines: args.parse_deflines,
            database: None,
            archive,
            out: args.out.as_deref(),
            search_error: &search_error,
            xml_scan,
        },
        outputs,
    )
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
    database: Option<DatabaseReports<'a>>,
    /// The BLAST archive (outfmt 11), when it is requested.
    archive: Option<&'a crate::report::archive::ArchiveHandle>,
    /// The file of `api::run_to_file`, next to which outfmt 13 and 14 write their side files.
    out: Option<&'a std::path::Path>,
    search_error: &'a dyn Fn(anyhow::Error) -> anyhow::Error,
    /// The read of the queries before the search of BLAST XML and the XML2 and JSON documents
    /// (`scan_xml_titles`).
    xml_scan: super::blast_xml::QueryScan,
}

/// BLAST XML and the XML2 and JSON documents refuse a query title that NCBI cannot write before
/// the first byte, in the reports that they write; a `-subject` XML2 or JSON document never writes
/// the last query of a run that stops at a later batch (`blast_xml::scan_queries`, on the records
/// in memory, read into the search's reader batches as the command line's `prescan_queries` reads
/// the file, and at the same place: before the subjects' warnings). With outfmt 5 among the
/// formats every report is written. `keep_whole`: the search keeps the whole letters of a cut
/// record (the archive's).
#[allow(clippy::too_many_arguments)]
fn scan_xml_titles(
    resolved: &[ResolvedOptions],
    query_records: &[fasta_reader::FastaRecord],
    lowercase: bool,
    query_range: Option<SequenceRange>,
    batch_size: usize,
    keep_whole: bool,
    subject_set: bool,
) -> Result<super::blast_xml::QueryScan> {
    let structured: Vec<u8> = resolved
        .iter()
        .map(|options| options.outfmt)
        .filter(|outfmt| matches!(outfmt, 5 | 13..=16))
        .collect();
    if structured.is_empty() || batch_size == 0 {
        return Ok(super::blast_xml::QueryScan::default());
    }
    let holds_back = subject_set && !structured.contains(&5);
    super::blast_xml::scan_queries(&resolved[0], holds_back, &mut |deliver| {
        reader_batches(
            query_records,
            lowercase,
            query_range,
            batch_size,
            keep_whole,
            deliver,
            &mut |_| Ok(()),
        )
    })
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

/// The reader batches of `query_records` as BLASTX's reader makes them
/// (`input::emit_fasta_batch`), given to `deliver` in order: a record's messages (`warn`), then
/// the record, and the batch once it holds `batch_size` letters or more; the rest at the end. With
/// `query_range` (`-query_loc`) a record is cut to its interval (its whole letters kept with
/// `keep_whole`, the archive's) and the batch counts the whole record; a record whose interval
/// starts more than one letter past its end is not searched and adds nothing, so the last batch
/// can hold none (NCBI's `Empty CBlastQueryVector`). Only the batch being read is kept.
fn reader_batches(
    query_records: &[fasta_reader::FastaRecord],
    lowercase: bool,
    query_range: Option<SequenceRange>,
    batch_size: usize,
    keep_whole: bool,
    deliver: &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
    warn: &mut dyn FnMut(&[u8]) -> Result<()>,
) -> Result<()> {
    let mut batch = Vec::new();
    let (mut letters, mut records_read) = (0, 0);
    for record in query_records {
        let mut query = blastx_record(record, lowercase)?;
        if !record.warnings.is_empty() {
            warn(&record.warnings)?;
        }
        let length = query.sequence.len();
        records_read += 1;
        match query_range.map(|range| record_interval(&range, length)) {
            None => {}
            Some(RecordInterval::Letters { from, to }) => {
                // The whole letters of a record whose letter before the interval is lower
                // case, for the chunks of a split query (`split::x_chunk_masks`).
                let keep_whole = keep_whole || query.lower_case_before(from);
                query.cut_in_place(from, to + 1, keep_whole)
            }
            Some(RecordInterval::Empty) => query.cut_in_place(length, length, keep_whole),
            Some(RecordInterval::PastEnd) => continue,
        }
        letters += length;
        batch.push(query);
        if letters >= batch_size {
            deliver(&batch)?;
            batch.clear();
            (letters, records_read) = (0, 0);
        }
    }
    if records_read > 0 {
        deliver(&batch)?;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:254-254
// ```c++
//         formatter.PrintProlog();
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:297-297
// ```c++
//         formatter.PrintEpilog(opt);
// ```
/// The prologs, the query batches (`run_query_batches`) with the reports of each in every format,
/// the epilogs and the hit list of a search of records.
fn report_batches(run: &RunReports<'_>, outputs: &mut ReportOutputs<'_>) -> Result<()> {
    let &RunReports {
        lowercase,
        query_records,
        query_range,
        batch_size,
        resolved,
        subjects,
        shown,
        ref database,
        archive,
        out,
        ref xml_scan,
        ..
    } = run;
    let ReportOutputs {
        formats,
        diagnostics,
        hits: hits_sink,
        observer,
    } = outputs;
    // The caller's receivers of the side files of outfmt 13 and 14.
    let side_sinks: Vec<_> = formats
        .iter()
        .map(|format| format.sink.side_files().cloned())
        .collect();
    let writers = formats
        .iter_mut()
        .map(|format| format.sink.ostream())
        .collect::<io::Result<Vec<_>>>()?;
    // The archive (outfmt 11) refuses a title where it writes the batch, after its search
    // (`archive::check_batch_titles`, below): when a batch that the search may write has a query
    // title it cannot write (`archive::batches_may_refuse`, on the records in memory) or the
    // subjects have one, every format is held in order from its prolog on (`archive::ArchiveHold`),
    // so that the refusal follows no byte, as on the command line. Otherwise no batch's check can
    // refuse, and the batches are not checked again (`titles_may_refuse`).
    let titles_may_refuse = archive.is_some()
        && batch_size > 0
        && (archive::batches_may_refuse(&resolved[0], true, false, &mut |deliver| {
            reader_batches(
                query_records,
                lowercase,
                query_range,
                batch_size,
                true,
                deliver,
                &mut |_| Ok(()),
            )
        }) || (database.is_none()
            && archive::check_titles(subjects, "subject", resolved[0].program().upper).is_err()));
    // A run in which a `-subject` title (refused with its first hit) or the message of a query
    // (refused where NCBI does not search its batch) could still be refused by those documents
    // holds every format in order, as the archive's titles do and as the command line holds its
    // document (`blast_xml::HeldOutput`), so that the refusal follows no byte. So does a `-db` run
    // in which a later batch could still refuse the taxonomy database or a subject with parsed
    // Seq-ids (`BlastxDatabase::later_batch_may_refuse`; the scientific names are those of the
    // XML2 and JSON documents).
    let writes_names = resolved
        .iter()
        .any(|options| matches!(options.outfmt, 13..=16));
    let xml_may_refuse = resolved
        .iter()
        .any(|options| matches!(options.outfmt, 5 | 13..=16))
        && (xml_scan.messages_may_be_refused
            || match database {
                Some(run) => run.database.later_batch_may_refuse(writes_names),
                None => super::blast_xml::subjects_may_be_refused(subjects),
            });
    let hold = (titles_may_refuse || xml_may_refuse).then(archive::ArchiveHold::new);
    // The hold comes between a format's stream and the caller's writer, which has no stream of
    // NCBI's (`cli::NcbiOstream`): it keeps the bytes in order.
    let mut writers: Vec<crate::cli::NcbiOstream<'_>> = match &hold {
        Some(hold) => {
            hold.hold();
            writers
                .into_iter()
                .map(|writer| crate::cli::NcbiOstream::new(Box::new(hold.writer(writer))))
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
        &mut writers,
        &mut documents,
        titles_may_refuse,
        wants_hits,
        diagnostics,
        observer,
    );
    // What the hold kept. A run that ends: every format, in order. A run that stops with an
    // error: when the archive's titles are what it held for (`titles_may_refuse`), NCBI's bytes
    // if the error is NCBI's own (its error archive follows where NCBI writes one,
    // `ReportOutputs::finish_archive`) and nothing after LOSAT's refusal (`ArchiveHold::end_with`);
    // when only the documents are, nothing after one of their refusals and the bytes of the
    // batches before after any other error (`ArchiveHold::end_document_with`), as the command line.
    // The observer's events of a held run come before its bytes reach their sinks, so a position
    // read from a sink at an event is not that of the event (`FormatObserver`): the same for every
    // hold of this function, a `-subject` run's too.
    match (&hold, &written) {
        (Some(hold), Ok(_)) => hold.release()?,
        (Some(hold), Err(error)) if titles_may_refuse => hold.end_with(error),
        (Some(hold), Err(error)) => hold.end_document_with(error),
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

/// `report_batches` from the prologs to the epilogs, into `writers` (each format's, behind the
/// hold of a run that holds its formats) and `documents` (the writer of each XML2 or JSON
/// document): whatever ends it, `report_batches` then ends the hold. `titles_may_refuse`: a batch
/// that the search may write has a title the archive cannot write. Returns the run's hit list
/// (`wants_hits`) and the database sequences it names.
fn write_run(
    run: &RunReports<'_>,
    writers: &mut [crate::cli::NcbiOstream<'_>],
    documents: &mut [Option<super::blast_xml::RunXml2>],
    titles_may_refuse: bool,
    wants_hits: bool,
    diagnostics: &mut &mut (dyn Write + Send),
    observer: &mut Option<&mut (dyn crate::api::local_blast::FormatObserver + Send)>,
) -> Result<(
    Vec<crate::report::pairwise::PairwiseHit>,
    Vec<ReportedSubject>,
)> {
    let &RunReports {
        lowercase,
        query_records,
        query_range,
        batch_size,
        choices,
        resolved,
        subjects,
        shown,
        search_subjects,
        serial_batches,
        parse_deflines,
        ref database,
        archive,
        search_error,
        ref xml_scan,
        ..
    } = run;
    // Under `-parse_deflines`, the scopes of the batches' queries (`QueryScopes`): with `-subject`
    // or outfmt 5 one scope of every record in memory, made before the search so that its refusal
    // comes before any byte (`QueryScopes::read_ahead`); else the batches add theirs as they come.
    let xml = choices.iter().any(|choice| choice.number == 5);
    let query_scopes = std::cell::RefCell::new(
        match crate::blastinput::query_batch::QueryScopes::of_search(
            parse_deflines,
            database.is_none(),
            xml,
            "BLASTX",
        ) {
            Some(scopes)
                if scopes.span() == crate::blastinput::query_batch::QueryScopeSpan::Run =>
            {
                Some(
                    crate::blastinput::query_batch::QueryScopes::of_records(
                        crate::blastinput::query_batch::QueryScopeSpan::Run,
                        "BLASTX",
                        query_records,
                    )?
                    .read_ahead(),
                )
            }
            scopes => scopes,
        },
    );
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:254-254
    // ```c++
    //         formatter.PrintProlog();
    // ```
    for (index, (writer, options)) in writers.iter_mut().zip(resolved).enumerate() {
        let mut quiet = io::sink();
        let warnings: &mut dyn Write = if index == 0 {
            &mut **diagnostics
        } else {
            &mut quiet
        };
        report::write_prolog(writer, warnings, options, subjects, &shown).map_err(report_error)?;
        if let Some(document) = &mut documents[index] {
            let mut start = Vec::new();
            document.document_mut().start(&mut start);
            writer.write_all(&start).map_err(output_error)?;
        }
        writer.flush().map_err(output_error)?;
    }
    // NCBI: a batch size of 0 makes the first batch empty, and its query factory refuses it.
    if batch_size == 0 {
        return Err(engine("Empty CBlastQueryVector"));
    }
    // The search's options are those of every format; the formats differ only in the counts of
    // the pairwise report (`run_report_counts`).
    let search_options = &resolved[0];
    // The archive's Bioseq of a cut record holds its whole letters (`FastaRecord::cut_in_place`).
    let keep_whole = archive.is_some();
    let mut read = |deliver: &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
                    warn: &mut dyn FnMut(&[u8]) -> Result<()>|
     -> Result<()> {
        reader_batches(
            query_records,
            lowercase,
            query_range,
            batch_size,
            keep_whole,
            deliver,
            warn,
        )
    };
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:253-259
    // ```c++
    //         formatter.PrintProlog();
    //         ...
    //         /*** Process the input ***/
    //         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
    //         ...
    //             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
    // ```
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
    let asks = ReportedAsks::of_all(resolved);
    let mut render = |diagnostics: &mut dyn Write,
                      queries: &[FastaRecord],
                      results: &mut Vec<BatchResults>,
                      queries_before: usize|
     -> Result<()> {
        if let Some(scopes) = query_scopes.borrow_mut().as_mut() {
            scopes.add_batch_ids(
                queries.iter().filter_map(|record| record.believed()),
                queries_before,
            )?;
        }
        // Every format and the hit list read the lists in the order of the reports.
        report::sort_lists(results);
        // A database: the batch's reported sequences, and its hit lists named by their ranks.
        let mut reported = match database {
            Some(run) => Some(reported_batch(
                run.database,
                search_subjects,
                search_options,
                asks,
                results,
                query_scopes
                    .borrow()
                    .as_ref()
                    .map(|scopes| scopes.numbered_from(queries_before)),
                &mut *diagnostics,
            )?),
            None => None,
        };
        // The archive (outfmt 11) of this batch, which NCBI writes once its search and the read of
        // its definition lines have not failed: the titles it cannot write are refused before any
        // format of the batch (`archive::check_batch_titles`; the run holds its formats, so the
        // refusal follows no byte of the batches before); only a run whose read before the search
        // found such a title can refuse (`titles_may_refuse`).
        if titles_may_refuse {
            archive::check_batch_titles(
                queries,
                database.is_none().then_some(subjects),
                queries_before == 0,
                search_options.program().upper,
            )?;
        }
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
                )
                .map_err(report_error)?,
            ),
            false => None,
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:82-89
        // ```c++
        //     // find the first valid context corresponding to this query
        //     for (i = 0; i < context_per_query; i++) {
        //         BlastContextInfo *ctx = query_info->contexts +
        //                                 query_number * context_per_query + i;
        //         if (ctx->is_valid) {
        //             m_SearchSpace = ctx->eff_searchsp;
        //             m_LengthAdjustment = ctx->length_adjustment;
        //             break;
        // ```
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
        // ```c++
        //             	ITERATE(CSearchResultSet, result, *results) {
        //                	    formatter.PrintOneResultSet(**result, query_batch);
        // ```
        for (index, (writer, options)) in writers.iter_mut().zip(resolved).enumerate() {
            // The warnings of the batch's queries are written once, with the first format.
            let mut quiet = io::sink();
            let warnings: &mut dyn Write = if index == 0 {
                &mut *diagnostics
            } else {
                &mut quiet
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
                archive::write_blastx_archive(
                    &archive::ArchiveCall {
                        handle,
                        batch: archive::CallBatch {
                            options,
                            queries,
                            requests: &archive::requests(subject_records, "BLASTX"),
                            batches: results,
                        },
                        subjects: call_subjects,
                        names,
                    },
                    writer,
                )
                .map_err(archive_error)?;
                writer.flush().map_err(output_error)?;
                continue;
            }
            let mut probe = observer
                .as_deref_mut()
                .map(|observer| FormatProbe::new(observer, index));
            report::render_queries_with(
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
                Some(LocalReport {
                    probe: probe.as_mut(),
                    first_hsp,
                    hits: hits.as_deref(),
                }),
            )
            .map_err(report_error)?;
            writer.flush().map_err(output_error)?;
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
    let processed = run_query_batches(
        search_options,
        search_subjects,
        serial_batches,
        &mut read,
        &mut **diagnostics,
        search_error,
        &mut render,
        // The library and LOSAT Web write no search strategy (`api::search::file_options`).
        &mut |_, _| Ok(()),
    )
    .map_err(input_error)?;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:297-297
    // ```c++
    //         formatter.PrintEpilog(opt);
    // ```
    for (index, (writer, options)) in writers.iter_mut().zip(resolved).enumerate() {
        let mut quiet = io::sink();
        let warnings: &mut dyn Write = if index == 0 {
            &mut **diagnostics
        } else {
            &mut quiet
        };
        report::write_epilog(
            writer,
            warnings,
            options,
            subjects,
            &shown,
            processed,
            documents[index].as_mut(),
        )
        .map_err(report_error)?;
        writer.flush().map_err(output_error)?;
    }
    Ok((run_hits, run_subjects))
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:195-200
// ```c++
//         /*** Initialize the database/subject ***/
//         CRef<CBlastDatabaseArgs> db_args(m_CmdLineArgs->GetBlastDatabaseArgs());
//         CRef<CLocalDbAdapter> db_adapter;
//         CRef<CScope> scope;
//         InitializeSubject(db_args, m_OptsHndl, m_CmdLineArgs->ExecuteRemotely(),
//                          db_adapter, scope);
// ```
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
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // NCBI's options factory checks the task first (`TaskParser::check`); the command line's
    // parser has checked it, a caller of `run_local` has not.
    super::args::BLASTX_TASK.check(&mut args.task)?;
    args.matrix_search = Some(settings.matrix_search.clone());
    let requested = requested_formats(&args, outputs);
    let (choices, resolutions) = read_format_options(&args, &requested)?;
    let archive_requested = choices.iter().any(|choice| choice.number == 11);
    let archive = archive::handle_for(
        &args,
        archive_requested,
        &crate::report::archive::OptionsHandle::default(),
    )?
    .map(|handle| outputs.with_archive_outfmt(handle));
    crate::utils::threading::validate_threads(args.num_threads as usize)?;
    match archive {
        Some(handle) => {
            let result = outputs.with_archive_diagnostics(handle.clone(), |view| {
                search_database(
                    &args,
                    query_records,
                    settings,
                    (choices, resolutions),
                    Some(&handle),
                    view,
                )
            });
            outputs.finish_archive(&handle, query_records.is_empty(), result)
        }
        None => search_database(
            &args,
            query_records,
            settings,
            (choices, resolutions),
            None,
            outputs,
        ),
    }
}

/// `run_local_database` from the taxid or seqid list on.
fn search_database(
    args: &BlastxArgs,
    query_records: &[fasta_reader::FastaRecord],
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
    (choices, resolutions): (Vec<app::FormatChoice>, Vec<Result<ResolvedOptions>>),
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
        .resolve(&settings, &mut *outputs.diagnostics)?;
    check_filtering(args)?;
    let query_range = parse_query_range(args)?;
    let mut resolved = check_resolutions(args, &choices, resolutions, &mut *outputs.diagnostics)?;
    // The database opens after the options are processed and before the query is read, as NCBI's
    // `InitializeSubject`.
    let database = BlastxDatabase::open_restricted(
        name,
        &settings,
        list.as_ref(),
        &mut *outputs.diagnostics,
        super::args::ProgramName::BLASTX.upper,
    )?
    .showing_gis(args.show_gis);
    crate::algorithm::common::database::note_archive_database_filter(
        archive,
        None,
        &args.db_filter,
        database.list_taxids(),
        database.statistics_length_is_empty(),
        args.searchsp.is_some_and(|space| space > 0),
        args.dbsize.is_some_and(|size| size > 0),
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
                &database.pairwise_database().definition,
                Program::Blastx,
            )?;
        }
        if matches!(choice.number, 5 | 15 | 16) {
            super::blast_xml::check_database_name(database.name())?;
        }
    }
    if query_records.is_empty() {
        outputs
            .diagnostics
            .write_all(b"Warning: [blastx] Query is Empty!\n")?;
        return Ok(());
    }
    // The letters of a reader batch (`native::run_inner`).
    let batch_size = app::query_batch_size("BLASTX", BATCH_SIZE as u32)? as usize;
    app::check_search_environment("BLASTX", false, outputs.only_archive())?;
    let pairwise = database.pairwise_database();
    let no_reported = ReportedSubjects::default();
    let shown = ShownSubjects::Database {
        name: database.name(),
        pairwise: &pairwise,
        reported: &no_reported,
    };
    // The database's sequences, read in blocks where they are searched; the empty database's error
    // belongs to the first batch with a valid context (`native::run_inner`).
    let (longest, unreadable) = match database.x_open_stream() {
        Ok(longest) => (Some(longest), None),
        Err(error) => (None, Some(error)),
    };
    let stream =
        longest.map(|longest| super::x_db_stream::XStreamDb::new(database.x_seqdb(), longest));
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
    let xml_scan = scan_xml_titles(
        &resolved,
        query_records,
        args.lcase_masking,
        query_range,
        batch_size,
        archive.is_some(),
        false,
    )?;
    report_batches(
        &RunReports {
            lowercase: args.lcase_masking,
            query_records,
            query_range,
            batch_size,
            choices: &choices,
            resolved: &resolved,
            subjects: &[],
            shown,
            search_subjects,
            serial_batches: args.parse_deflines,
            parse_deflines: args.parse_deflines,
            database: Some(DatabaseReports {
                database: &database,
                pairwise: &pairwise,
            }),
            archive,
            out: args.out.as_deref(),
            search_error: &search_error,
            xml_scan,
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
/// (`OutputSink::Files`, which LOSAT Web gives them).
pub fn check_options(args: &BlastxArgs, outfmts: &[&str]) -> Result<()> {
    let requested: Vec<(&str, bool)> = outfmts.iter().map(|&outfmt| (outfmt, true)).collect();
    let (choices, resolutions) = read_options(args, &requested, false)?;
    crate::utils::threading::validate_threads(args.num_threads as usize)?;
    parse_subject_range(args)?;
    check_filtering(args)?;
    parse_query_range(args)?;
    check_resolutions(args, &choices, resolutions, &mut io::sink()).map(|_| ())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2540-2545
// ```c++
//         TSeqRange subj_range;
//         if (args.Exist(kArgSubjectLocation) && args[kArgSubjectLocation]) {
//             subj_range =
//                 ParseSequenceRange(args[kArgSubjectLocation].AsString(),
//                             "Invalid specification of subject location");
//         }
// ```
/// `-subject_loc`, which the command line reads after it opened the subject file and before it
/// reads the subjects (`native::run_inner`).
fn parse_subject_range(args: &BlastxArgs) -> Result<Option<SequenceRange>> {
    parse_optional_range(args.subject_loc.as_deref(), RangeRole::Subject, "BLASTX")
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1995-1999
// ```c++
//     // set the sequence range
//     if (args.Exist(kArgQueryLocation) && args[kArgQueryLocation]) {
//         m_Range = ParseSequenceRange(args[kArgQueryLocation].AsString(),
//                                      "Invalid specification of query location");
//     }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:87-88,106-107
// ```c++
//     arg.Reset(new CFilteringArgs( !kQueryIsProtein ));
//     m_Args.push_back(arg);
// ...
//     m_QueryOptsArgs.Reset(new CQueryOptionsArgs(kQueryIsProtein));
//     arg.Reset(m_QueryOptsArgs);
// ```
/// The filtering handler (`CFilteringArgs`), which comes after the subjects are read and before
/// the query options and formatting handlers: NCBI's error of `-seg` (`BlastxArgs::seg_spec`;
/// `resolve_with_format` keeps it with the other deferred errors, LOSAT's refusal of a value that
/// NCBI reads comes from `resolve_with_format`).
pub(super) fn check_filtering(args: &BlastxArgs) -> Result<()> {
    match args.seg_spec() {
        Err(app::SegOptionError::Ncbi(message)) => {
            Err(options_error(SemanticOptionsError::from(message).into()))
        }
        Ok(_) | Err(app::SegOptionError::Losat(_)) => Ok(()),
    }
}

/// `-query_loc`, which the command line reads after the subjects, with the query options and
/// before the formatting handler (`native::run_inner`).
fn parse_query_range(args: &BlastxArgs) -> Result<Option<SequenceRange>> {
    parse_optional_range(args.query_loc.as_deref(), RangeRole::Query, "BLASTX")
}

/// The output formats `outfmts` as `run_local` reads them and the options of each, before the
/// records are read: LOSAT Web's refusals and LOSAT's fail here; NCBI's errors of the options
/// (`SemanticOptionsError`) stay in the resolutions, which the command line reports after the
/// subjects are read (`native::run_inner`).
fn read_options(
    args: &BlastxArgs,
    outfmts: &[(&str, bool)],
    ids_read: bool,
) -> Result<(Vec<app::FormatChoice>, Vec<Result<ResolvedOptions>>)> {
    // The host has no database files; LOSAT Web searches the records it gives.
    if args.db.is_some() {
        bail!("-db (a database search) is not supported in the browser yet");
    }
    // The records keep no ids read from their deflines (the host registers them before it knows
    // the options).
    if args.parse_deflines && !ids_read {
        bail!("-parse_deflines is not supported in the browser yet");
    }
    read_format_options(args, outfmts)
}

/// The `-outfmt` value of each output of `outputs`, and whether the side files of outfmt 13 and 14
/// have somewhere to go: a file next to which they are written (`api::run_to_file`), or the
/// caller's receiver of them (`FormatOutput::side_out`).
fn requested_formats<'o>(
    args: &BlastxArgs,
    outputs: &'o ReportOutputs<'_>,
) -> Vec<(&'o str, bool)> {
    outputs
        .formats
        .iter()
        .map(|format| {
            (
                format.outfmt,
                format.side_out(args.out.as_deref()).is_some(),
            )
        })
        .collect()
}

/// `read_options` without its refusals of what a search of records in memory cannot do: each
/// format comes with whether its side files have somewhere to go (`requested_formats`).
fn read_format_options(
    args: &BlastxArgs,
    outfmts: &[(&str, bool)],
) -> Result<(Vec<app::FormatChoice>, Vec<Result<ResolvedOptions>>)> {
    // NCBI parses -outfmt before its option handlers read the subjects (`native::run_inner`).
    let mut choices = Vec::with_capacity(outfmts.len());
    for &(outfmt, side_files) in outfmts {
        let choice = app::parse_formatting_string(outfmt)?;
        app::report_format(&choice, Program::Blastx, args.out.as_deref())?;
        // No NCBI counterpart: the side files of outfmt 13 and 14 have nowhere to go without a file
        // or a receiver of them.
        if !side_files && matches!(choice.number, 13 | 14) {
            bail!(
                "output format {} is not supported in memory for BLASTX",
                choice.number
            );
        }
        choices.push(choice);
    }
    ensure!(!choices.is_empty(), "no output format was requested");
    // The options of each format: LOSAT's refusals now, NCBI's errors of the options after the
    // subjects are read (`native::run_inner`).
    let mut resolutions = Vec::with_capacity(choices.len());
    for choice in &choices {
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
/// order: the formatting handler's, the hit list size (its warnings and those of an ignored
/// -num_descriptions, the first format's, to `diagnostics`; formats that need different hit lists
/// are refused, `run_report_counts`), then NCBI's errors of each format's options.
fn check_resolutions(
    args: &BlastxArgs,
    choices: &[app::FormatChoice],
    resolutions: Vec<Result<ResolvedOptions>>,
    diagnostics: &mut dyn Write,
) -> Result<Vec<ResolvedOptions>> {
    for choice in choices {
        app::formatting_handler_check(choice, Program::Blastx)?;
    }
    let outfmts: Vec<i32> = choices.iter().map(|choice| choice.number).collect();
    let (_, warnings) = crate::blastinput::report_counts::run_report_counts_with(
        "blastx",
        &outfmts,
        super::args::max_target_seqs_count(args.max_target_seqs),
        args.num_descriptions,
        args.num_alignments,
        &args.formatting,
    )?;
    diagnostics.write_all(&warnings)?;
    resolutions
        .into_iter()
        .map(|resolution| resolution.map_err(options_error))
        .collect()
}

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1082-1090
// ```c++
//     m_MaskRangeStart = GetCurrentPos(ePosWithGapsAndSegs);
// }
//
// void CFastaReader::x_CloseMask(void)
// {
//     _ASSERT(m_MaskRangeStart != kInvalidSeqPos);
//     m_CurrentMask->SetPacked_int().AddInterval
//         (GetBestID(), m_MaskRangeStart, GetCurrentPos(ePosWithGapsAndSegs) - 1,
//          eNa_strand_plus);
// ```
/// A record of NCBI's reader (`fasta_reader`, the records that LOSAT Web registers) as BLASTX's
/// reader (`input.rs`) reads the same input: the residues in upper case and, with
/// `-lcase_masking` (`lowercase`), the lowercase mask's intervals, which `fasta_reader` keeps as
/// the case of the residues (a `>?` gap inside a masked interval, too); a title that is not UTF-8
/// is BLASTX's reader's error. The reader's messages stay with the caller's record.
fn blastx_record(record: &fasta_reader::FastaRecord, lowercase: bool) -> Result<FastaRecord> {
    let title = String::from_utf8(record.title.clone()).context("invalid UTF-8 FASTA title")?;
    if record.sequence.len() > i32::MAX as usize {
        bail!("BLASTX sequence exceeds Int4 length");
    }
    let mut sequence = record.sequence.clone();
    let mut lowercase_masks: Vec<(i32, i32)> = Vec::new();
    for (offset, residue) in sequence.iter_mut().enumerate() {
        if !residue.is_ascii_lowercase() {
            continue;
        }
        residue.make_ascii_uppercase();
        if lowercase {
            // The length is at most `i32::MAX` (above).
            let offset = offset as i32;
            match lowercase_masks.last_mut().filter(|r| r.1 + 1 == offset) {
                Some(last) => last.1 = offset,
                None => lowercase_masks.push((offset, offset)),
            }
        }
    }
    Ok(FastaRecord {
        internal_id: record.local_id.clone(),
        title,
        sequence,
        lowercase_masks,
        warnings: Vec::new(),
        source: record.source.clone(),
    })
}
