//! TBLASTX's application: the command line, the input and the loop over the query batches.
//!
//! This file is a child module of `blast_engine/run_impl.rs` (declared there with `#[path]`, as
//! `x_seed_release.rs` is), so that it reaches the run's types and the search of a query batch
//! without changing them.
//!
//! NCBI reference: `c++/src/app/blast/tblastx_app.cpp:106-209` (`CTblastxApp::Run`), the option
//! checks of `c++/src/algo/blast/blastinput/blast_args.cpp:1995-3636`, and
//! `c++/src/algo/blast/api/local_blast.cpp:197-293` (`CLocalBlast::Run`).

use super::blast_format::{report_order, write_pairwise_prologs, BatchReports, SearchedBatches};
use super::{
    search_query_batch, x_batch_pipe, ErrorAfterBatchReports, QueryRuns, RunSubjects,
    ShownSubjects, SubjectSource,
};
use crate::algorithm::blastn::blast_engine::{QueryReading, NO_DATA_MESSAGE};
use crate::algorithm::tblastx::args::TblastxArgs;
use crate::algorithm::tblastx::blast_engine::{report, OutputSink, ReportOutputs, TblastxBatch};
use crate::algorithm::tblastx::database::DatabaseTitles;
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
use crate::algorithm::tblastx::diagnostics::diagnostics_enabled;
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
use crate::algorithm::tblastx::stage_dump;
use crate::algorithm::tblastx::subjects::{SubjectSetNumbers, Subjects, SubjectsUnreadable};
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
use crate::algorithm::tblastx::tracing::trace_hsp_target;
use crate::blastinput::fasta_reader::{
    read_queries, read_subjects, FastaInputSource, FastaRecord, QueryEnd, QueryRecords,
    ReaderConfig,
};
use crate::blastinput::report_counts::ReportCounts;
use crate::blastinput::seq_range::{
    cut_queries, cut_subjects, next_ranged_batch, parse_optional_range, subjects_read, Placements,
    QueryInput, RangeRole, SequenceRange,
};
use anyhow::{Context, Result};
use std::sync::Arc;

/// `settings` are the NCBI application settings that LOSAT reproduces
/// (`ncbi_environment::check_ncbi_application_settings`): the input readers use them once
/// the program reads with `fasta_reader` (steps S3-S8 of the port plan).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_scope_src.cpp:67-72
/// ```c++
///     CNcbiApplication* app = CNcbiApplication::Instance();
///     if (app) {
///         const CNcbiRegistry& registry = app->GetConfig();
///         x_LoadDataLoadersConfig(registry);
///         x_LoadBlastDbDataLoaderConfig(registry);
///     }
/// ```
pub fn run(
    args: TblastxArgs,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
) -> Result<()> {
    // NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:106-111
    // ```c
    // if(RecoverSearchStrategy(args, m_CmdLineArgs)) {
    // 	opts_hndl.Reset(&*m_CmdLineArgs->SetOptionsForSavedStrategy(args));
    // }
    // else {
    // 	opts_hndl.Reset(&*m_CmdLineArgs->SetOptions(args));
    // }
    // ```
    // NCBI's option handlers read the subjects first, then open the query and the output,
    // then process the formatting options (blast_args.cpp:3631-3636). The thread count is
    // checked before the inputs are read, as before. LOSAT's own debug and timing output
    // (LOSAT_WASI_THREADS_DEBUG, LOSAT_TIMING, LOSAT_DEBUG_SCAN_SOFF) follows the reads and
    // does not count their time.
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3624-3627
    // ```c
    //     if (GetExportSearchStrategyStream(args) ||
    //            m_FormattingArgs->ArchiveFormatRequested(args)) {
    //         locality = CBlastOptions::eBoth;
    //     }
    // ```
    // `ArchiveFormatRequested` parses `-outfmt` (blast_args.cpp:2745-2748) before the
    // option handlers run; a format that NCBI runs and LOSAT does not write is rejected
    // there too.
    let choice = crate::blastinput::app::parse_formatting_string(&args.outfmt)?;
    let format = crate::blastinput::app::report_format(
        &choice,
        crate::blastinput::app::Program::Tblastx,
        args.out.as_deref(),
    )?;
    crate::utils::threading::validate_threads(args.num_threads)?;
    use crate::blastinput::input_files;
    // `-db`: no subject file; the database is opened after the options are processed
    // (`search_cli`), where NCBI's `InitializeSubject` opens it (tblastx_app.cpp:113-119).
    let (subjects, subject_placements) = if args.db.is_some() {
        (Vec::new(), Placements::default())
    } else {
        read_subject_set(&args, settings.clone())?
    };
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3456-3481
    // ```c
    //     if (args.Exist(kArgQuery) && args[kArgQuery].HasValue() &&
    //         m_InputStream == NULL) {
    // ...
    //         else {
    //             m_InputStream = &args[kArgQuery].AsInputFile();
    //         }
    //     }
    // ...
    //     else {
    //         m_OutputStream = &args[kArgOutput].AsOutputFile();
    //     }
    // ```
    // The query is opened, then the output file is created (`-` is standard output); the
    // query is read after both, so an output file that is the query file reads empty.
    input_files::check_utf8_file_name(&args.query, "query", "TBLASTX")?;
    let query_file = input_files::open_input(&args.query, "query", "TBLASTX")?;
    if let Some(path) = args.out.as_deref() {
        input_files::check_utf8_file_name(path, "out", "TBLASTX")?;
    }
    let out_file = match args.out.as_deref().filter(|path| path.as_os_str() != "-") {
        Some(path) => Some(std::io::BufWriter::new(
            std::fs::File::create(path).map_err(|_| crate::cli::inaccessible("out", path))?,
        )),
        None => None,
    };
    let pairwise = format == Some(crate::blastinput::app::ReportFormat::Pairwise);
    // LOSAT's reports read the format number that NCBI reads.
    let outfmt = choice.normalized();
    let mut stderr = std::io::stderr();
    let mut stream = crate::cli::ReportStream {
        inner: match out_file {
            Some(file) => Box::new(file) as Box<dyn std::io::Write + Send>,
            None => crate::cli::report_standard_output(),
        },
        failed: false,
    };
    let result = {
        let mut outputs =
            ReportOutputs::single(&outfmt, OutputSink::Writer(&mut stream), &mut stderr);
        search_cli(
            args,
            query_file,
            &subjects,
            &subject_placements,
            &mut outputs,
            settings,
        )
    };
    // What was written before an error (such as the outfmt 0 prolog) stays in the file.
    let flushed = std::io::Write::flush(&mut stream);
    // NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.hpp:252-255
    // ```c
    //     catch (const std::ios::failure&) {                                      \
    //         LOG_POST(Error << "BLAST failed to write output");                  \
    //         exit_code = BLAST_OUTPUT_ERROR;                                     \
    //     }                                                                       \
    // ```
    // The outfmt 0 formatter's stream throws when a write fails (blast_format.cpp:118-119).
    // For outfmt 6/7 NCBI aborts instead, and LOSAT reports the error
    // (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES).
    if stream.failed && pairwise {
        return Err(crate::cli::NativeError {
            exit: 6,
            message: "BLAST failed to write output\n".to_string(),
        }
        .into());
    }
    result?;
    flushed.context("failed to write the output")
}

/// The subjects of `-subject` (after `-subject_loc`) as NCBI's option handlers read them, with
/// their messages.
fn read_subject_set(
    args: &TblastxArgs,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
) -> Result<(Vec<FastaRecord>, Placements)> {
    use crate::blastinput::input_files;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2553-2557
    // ```c
    //         CRef<blast::CBlastQueryVector> subjects;
    //         m_Scope = ReadSequencesToBlast(*subj_input_stream, IsProtein(),
    //                                        subj_range, parse_deflines,
    //                                        use_lcase_masks, subjects, m_IsMapper);
    //         m_Subjects.Reset(new blast::CObjMgr_QueryFactory(*subjects));
    // ```
    // The first handler opens and reads the subjects (an empty subject set fails there),
    // as BLASTN's `run`.
    let Some(subject_path) = args.subject.clone() else {
        return Err(crate::blastinput::app::missing_subject_error());
    };
    input_files::check_utf8_file_name(&subject_path, "subject", "TBLASTX")?;
    let subject_file = input_files::open_input(&subject_path, "subject", "TBLASTX")?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2537-2545
    // ```c++
    //             subj_input_stream = &args[kArgSubject].AsInputFile();
    //         }
    //
    //         TSeqRange subj_range;
    //         if (args.Exist(kArgSubjectLocation) && args[kArgSubjectLocation]) {
    //             subj_range =
    //                 ParseSequenceRange(args[kArgSubjectLocation].AsString(),
    //                             "Invalid specification of subject location");
    //         }
    // ```
    // The subject range is read after the subject file is opened, before it is read.
    let subject_range =
        parse_optional_range(args.subject_loc.as_deref(), RangeRole::Subject, "TBLASTX")?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:242-246
    // ```c++
    //     CRef<CBlastFastaInputSource> fasta(new CBlastFastaInputSource(in, iconfig));
    //     CRef<CBlastInput> input(new CBlastInput(fasta));
    //     CRef<CScope> scope(new CScope(*CObjectManager::GetInstance()));
    //     sequences = input->GetAllSeqs(*scope);
    //     return scope;
    // ```
    // NCBI reads the records one at a time, writes each one's messages as it reads it, and
    // checks each one's range after reading it: a range that starts past the end of a
    // record stops the reading there (`read_subjects`). Records without residues (and
    // intervals without letters) are read without a message and reported when the search
    // is set up (`search`).
    let mut subject_source = FastaInputSource::from_argument(
        &subject_path,
        subject_file,
        ReaderConfig::subject("TBLASTX", false, settings.data_loaders),
    );
    // Diagnostic (LOSAT_TIMING): the subject FASTA for the common line (utils/xtiming.rs). No NCBI
    // counterpart.
    let t_read_fasta = crate::utils::xtiming::clock();
    let read_subject_records = {
        let mut stderr = std::io::stderr();
        read_subjects(
            &mut subject_source,
            subject_range.as_ref(),
            &mut |message: &[u8]| std::io::Write::write_all(&mut stderr, message),
        )?
    };
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::ReadFasta, t_read_fasta);
    drop(subject_source);
    let (subjects, subject_placements) =
        match cut_subjects(&read_subject_records, subject_range.as_ref())? {
            Some((cut, placements)) => (cut, placements),
            None => (read_subject_records, Placements::default()),
        };
    check_subjects_not_empty(&subjects)?;
    Ok((subjects, subject_placements))
}

/// NCBI's error for a subject file without records, raised when it reads the subjects.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
/// ```c
/// CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
///     : m_QueryVector(& queries)
/// {
///     if (queries.Empty()) {
///         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
///     }
/// ```
/// NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.hpp:225-227
/// ```c
///             LOG_POST(Error << "BLAST engine error: " << e.GetMsg());        \
///             exit_code = BLAST_ENGINE_ERROR;                                 \
/// ```
fn check_subjects_not_empty(subjects: &[FastaRecord]) -> Result<()> {
    if subjects.is_empty() {
        return Err(crate::cli::NativeError {
            exit: 3,
            message: "BLAST engine error: Empty CBlastQueryVector\n".to_string(),
        }
        .into());
    }
    Ok(())
}

/// NCBI's processing of the options of tblastx after the files are opened: the handlers in
/// the order of `CTblastxAppArgs` (the -seg value; the formats of other programs and the
/// warning for a hit list size below 5), then `Validate` (the threshold and the word size
/// of the lookup table, then the e-value), with NCBI's errors.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2975-2977
/// ```c
///     if(hitlist_size < 5){
///    		ERR_POST(Warning << "Examining 5 or more matches is recommended");
///     }
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_options.c:1303-1310,1358-1364
/// ```c
///     if (program_number != eBlastTypeBlastn &&
///         program_number != eBlastTypeMapping &&
///         (!Blast_ProgramIsRpsBlast(program_number)) &&
///         options->threshold <= 0)
///     {
///         Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
///                          "Non-zero threshold required");
/// ...
///         else {
///             Blast_MessageWrite(blast_msg, eBlastSevError,
///                                kBlastMessageNoContext,
///                                "Word-size must be less "
///                                "than 6 for protein comparison");
///             return BLASTERR_OPTION_VALUE_INVALID;
///         }
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_options.c:1518-1523
/// ```c
///     if (options->expect_value <= 0.0 && options->cutoff_score <= 0)
///     {
///         Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
///          "expect value or cutoff score must be greater than zero");
///         return BLASTERR_OPTION_VALUE_INVALID;
///     }
/// ```
/// TBLASTX has no cutoff score option, so an `-evalue` of 0 (or one that reads as 0) fails.
fn check_ncbi_options(
    args: &TblastxArgs,
    outfmts: &[i32],
    diagnostics: &mut dyn std::io::Write,
) -> Result<(Option<SequenceRange>, ReportCounts)> {
    use crate::blastinput::app::{
        formatting_handler_check, options_error, parse_formatting_string, Program,
    };
    args.seg_spec()?;
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1995-1999
    // ```c++
    //     // set the sequence range
    //     if (args.Exist(kArgQueryLocation) && args[kArgQueryLocation]) {
    //         m_Range = ParseSequenceRange(args[kArgQueryLocation].AsString(),
    //                                      "Invalid specification of query location");
    //     }
    // ```
    // The query options handler comes after the filtering handler and before the
    // formatting handler (tblastx_args.cpp:44-125).
    let query_range = crate::blastinput::seq_range::parse_optional_range(
        args.query_loc.as_deref(),
        crate::blastinput::seq_range::RangeRole::Query,
        "TBLASTX",
    )?;
    formatting_handler_check(&parse_formatting_string(&args.outfmt)?, Program::Tblastx)?;
    // The hit list size and the counts of the pairwise report
    // (`blastinput::report_counts`), with the handler's warnings.
    let (counts, warnings) = crate::blastinput::report_counts::run_report_counts(
        "tblastx",
        outfmts,
        args.max_target_seqs,
        args.num_descriptions,
        args.num_alignments,
    )?;
    diagnostics.write_all(&warnings)?;
    // (`threshold <= 0` is false for NaN, as in C.)
    if args.threshold <= 0.0 {
        return Err(options_error("Non-zero threshold required"));
    }
    if args.word_size > 4 {
        return Err(options_error(
            "Word-size must be less than 6 for protein comparison",
        ));
    }
    // The hit saving options' check of the hit list size comes before the e-value's.
    counts.validate()?;
    if args.evalue <= 0.0 {
        return Err(options_error(
            "expect value or cutoff score must be greater than zero",
        ));
    }
    Ok((query_range, counts))
}

/// The number of `-outfmt` of the command line (`check_ncbi_options`).
fn cli_outfmts(args: &TblastxArgs) -> Result<[i32; 1]> {
    Ok([crate::blastinput::app::parse_formatting_string(&args.outfmt)?.number])
}

/// The options that LOSAT's TBLASTX rejects, where NCBI would start the search.
pub(crate) fn check_losat_limits(args: &TblastxArgs) -> Result<()> {
    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:214-224
    // ```c
    //     if (ewp->diag_table->multiple_hits) {
    //         status = s_BlastAaWordFinder_TwoHit(subject, query,
    //     ...
    //     } else {
    //         status = s_BlastAaWordFinder_OneHit(subject, query,
    // ```
    // LOSAT's TBLASTX scans with the two-hit word finder only; the one-hit word finder of
    // a window of 0 is not ported.
    if args.window_size == 0 {
        anyhow::bail!(
            "-window_size 0 (the one-hit word finder) is not supported by LOSAT's TBLASTX"
        );
    }
    // NCBI reference: c++/src/algo/blast/core/blast_aalookup.c:237
    // ```c
    //     lookup->word_length = opt->word_size;
    // ```
    // LOSAT's TBLASTX builds and scans 3-residue words only.
    if args.word_size != 3 {
        anyhow::bail!(
            "-word_size {} is not supported by LOSAT's TBLASTX (it implements word size 3)",
            args.word_size
        );
    }
    Ok(())
}

/// The checks of the options alone, without the inputs (the `validate` of web ABI v2):
/// NCBI's processing of the options and LOSAT's limits.
pub fn check_options(args: &TblastxArgs) -> Result<()> {
    // The host has no database files (`run_local_with`).
    if args.db.is_some() {
        anyhow::bail!("-db (a database search) is not supported by LOSAT Web yet");
    }
    check_ncbi_options(args, &cli_outfmts(args)?, &mut std::io::sink())?;
    check_losat_limits(args)?;
    // The host validates the argv it runs, with its `-num_threads` (docs/web/abi_v2.md).
    crate::utils::threading::validate_threads(args.num_threads)
}

/// The part of `run` after the output is opened: the options, `Query is Empty!`, the
/// reading of the queries and the search.
fn search_cli(
    mut args: TblastxArgs,
    query_file: std::fs::File,
    subjects: &[FastaRecord],
    subject_placements: &Placements,
    outputs: &mut ReportOutputs<'_>,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
) -> Result<()> {
    let (query_range, counts) =
        check_ncbi_options(&args, &cli_outfmts(&args)?, outputs.diagnostics)?;
    args.report_counts = counts;
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:113-119
    // ```c++
    //         /*** Initialize the database/subject ***/
    //         CRef<CBlastDatabaseArgs> db_args(m_CmdLineArgs->GetBlastDatabaseArgs());
    //         CRef<CLocalDbAdapter> db_adapter;
    //         CRef<CScope> scope;
    //         InitializeSubject(db_args, opts_hndl, m_CmdLineArgs->ExecuteRemotely(),
    //                          db_adapter, scope);
    // ```
    // `-db`: the database opens after the options are processed, before `Query is Empty!`
    // (`TblastxDatabase::open`).
    let database = match args.database_name() {
        Some(name) => Some(crate::algorithm::tblastx::database::TblastxDatabase::open(
            name, &settings,
        )?),
        None => None,
    };
    // LOSAT's refusals of a database for the Seq-annot formats, before the search
    // (`report::serial::seqalign::check_database`).
    if let Some(database) = database.as_ref() {
        for format in &outputs.formats {
            if let report::TblastxOutputFormat::Seqalign(serial) =
                report::output_format(format.outfmt)
            {
                crate::report::serial::seqalign::check_database(
                    serial,
                    database.name(),
                    &database.pairwise_database().definition,
                    crate::blastinput::app::Program::Tblastx,
                )?;
            }
        }
    }
    // The BLAST XML (outfmt 5, and the XML2 and JSON documents of 16 and 15) cannot write a
    // `-db` name outside US-ASCII: refused before the search (`names::check_database_name`).
    if let Some(database) = database.as_ref() {
        if outputs
            .formats
            .iter()
            .any(|format| report::output_format(format.outfmt).is_structured_xml())
        {
            crate::report::blast_xml::names::check_database_name("TBLASTX", database.name())?;
        }
    }
    let source = match &database {
        Some(database) => SubjectSource::Database(database),
        None => SubjectSource::Fasta {
            records: subjects,
            placements: subject_placements,
        },
    };
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:125-135
    // ```c++
    //         SDataLoaderConfig dlconfig =
    //             InitializeQueryDataLoaderConfiguration(query_opts->QueryIsProtein(),
    //                                                    db_adapter);
    //         CBlastInputSourceConfig iconfig(dlconfig, query_opts->GetStrand(),
    //                                      query_opts->UseLowercaseMasks(),
    //                                      query_opts->GetParseDeflines(),
    //                                      query_opts->GetRange());
    //         if(IsIStreamEmpty(m_CmdLineArgs->GetInputStream())){
    //            	ERR_POST(Warning << "Query is Empty!");
    //            	return BLAST_EXIT_SUCCESS;
    //         }
    // ```
    //
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-860
    // ```c++
    // 	char c;
    // 	CNcbiStreampos orig_p = in.tellg();
    // 	// Piped input
    // 	if(orig_p < 0)
    // 		return false;
    // ```
    // The query's data loaders are those of the subjects (`SDataLoaderConfig`). The position
    // is taken on the opened file, before it is read (`stream_is_empty`), as BLASTN's
    // `search_cli`. When the subjects were read from standard input too, `cin` has reached
    // its end (a failed stream), so NCBI gets no position for the query either, and its
    // reader is at its end: no query batch is read.
    let shared_standard_input = args.query.as_os_str() == "-"
        && args
            .subject
            .as_deref()
            .is_some_and(|path| path.as_os_str() == "-");
    let mut query_source = FastaInputSource::from_argument(
        &args.query,
        query_file,
        ReaderConfig::query("TBLASTX", false, settings.data_loaders),
    );
    if !shared_standard_input && query_source.stream_is_empty() {
        outputs
            .diagnostics
            .write_all(b"Warning: [tblastx] Query is Empty!\n")?;
        return Ok(());
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:176-178
    // ```c++
    //         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
    //
    //             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
    // ```
    // The records are read here, and their messages are written batch by batch
    // (`run_in_pool`). What ends the reading (the end of the input, an `eEOF`, a reader
    // error or LOSAT's rejection of a Seq-id line) comes in the batch that NCBI reads it in,
    // after the reports of the batches before (`run_in_pool`).
    // Diagnostic (LOSAT_TIMING): the query FASTA for the common line (utils/xtiming.rs). No NCBI
    // counterpart.
    let t_read_fasta = crate::utils::xtiming::clock();
    let QueryRecords { records, end } = read_queries(&mut query_source);
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::ReadFasta, t_read_fasta);
    drop(query_source);
    let reading = match end {
        QueryEnd::Input => QueryReading::End,
        QueryEnd::EofError => QueryReading::BlankLines,
        QueryEnd::Error { error, warnings } => QueryReading::Error { error, warnings },
    };
    search(args, &records, reading, query_range, source, outputs, None)
}

// NCBI reference: ncbi-blast/c++/src/app/blast/blast_formatter.cpp:429-467
// ```c
// CRef<CSearchResultSet> results = m_RmtBlast->GetResultSet();
// formatter.PrintProlog();
// ...
// ITERATE(CSearchResultSet, result, *results) {
//     ...
//         formatter.PrintOneResultSet(**result, queries);
//     ...
// }
// ```
// NCBI formats one result set without searching again; several requested formats
// are several CBlastFormat printers over the same result set.
/// Runs one TBLASTX search over already parsed records and writes every requested
/// output format from the same result (the shared entry of the CLI, web ABI v1 and
/// v2), and gives the caller the final hit list (`ReportOutputs::hits`).
///
/// Each requested `-outfmt` is validated as on the command line, before the search
/// starts. The `-query` and `-subject` values of `args` are used only as display names.
/// The records are those of NCBI's reader (`blastinput::fasta_reader`; their messages
/// are written where NCBI reads them): queries without a record are NCBI's
/// `Query is Empty!`.
pub fn run_local(
    args: TblastxArgs,
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    run_local_with(args, query_records, subject_records, outputs, None)
}

/// `run_local`, with the checks that ABI v1 makes of its `bio` records (`v1`,
/// `web_api::v1_tblastx::run_web_pair`).
pub(crate) fn run_local_with(
    mut args: TblastxArgs,
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    outputs: &mut ReportOutputs<'_>,
    v1: Option<&dyn V1Checks>,
) -> Result<()> {
    // The host has no database files; LOSAT Web searches the records it gives.
    if args.db.is_some() {
        anyhow::bail!("-db (a database search) is not supported by LOSAT Web yet");
    }
    // NCBI parses -outfmt before its option handlers read the subjects (BLASTN's
    // `run_local`); `search` checks the formats again.
    check_report_formats(outputs, None)?;
    // No NCBI counterpart: ABI v1 is frozen (plan TD-1; line OF, outfmt-20261019 DECISIONS
    // 10): it adds no format, so it refuses the Seq-annot formats (outfmt 8, 9, 12) with the
    // message that `check_report_formats` gave them before TBLASTX wrote them.
    if v1.is_some() {
        for format in &outputs.formats {
            let number = crate::blastinput::app::parse_formatting_string(format.outfmt)?.number;
            if matches!(number, 8 | 9 | 12) {
                anyhow::bail!("output format {number} is not supported by LOSAT's TBLASTX");
            }
        }
    }
    // The subject range and its record checks come where NCBI reads the subjects (`run`):
    // the messages of each record read, up to a record whose range starts past its end.
    let subject_range =
        parse_optional_range(args.subject_loc.as_deref(), RangeRole::Subject, "TBLASTX")?;
    for record in &subject_records[..subjects_read(subject_records, subject_range.as_ref())] {
        outputs.diagnostics.write_all(&record.warnings)?;
    }
    let ranged_subjects = cut_subjects(subject_records, subject_range.as_ref())?;
    check_subjects_not_empty(subject_records)?;
    let mut outfmts = Vec::with_capacity(outputs.formats.len());
    for format in &outputs.formats {
        outfmts.push(crate::blastinput::app::parse_formatting_string(format.outfmt)?.number);
    }
    let (query_range, counts) = check_ncbi_options(&args, &outfmts, outputs.diagnostics)?;
    args.report_counts = counts;
    let whole_subjects = Placements::default();
    let (searched_subjects, subject_placements) = match &ranged_subjects {
        Some((cut, placements)) => (cut.as_slice(), placements),
        None => (subject_records, &whole_subjects),
    };
    search(
        args,
        query_records,
        QueryReading::Records,
        query_range,
        SubjectSource::Fasta {
            records: searched_subjects,
            placements: subject_placements,
        },
        outputs,
        v1,
    )
}

/// Each requested format as NCBI reads `-outfmt`: NCBI's errors, and a rejection of a
/// format that LOSAT's TBLASTX does not write.
///
/// `out` is the run's `-out`: outfmt 13 and 14 write side files next to it, and fail on standard
/// output (`xinclude_check_for`). `run_local` has no output file.
fn check_report_formats(outputs: &ReportOutputs<'_>, out: Option<&std::path::Path>) -> Result<()> {
    use crate::blastinput::app;
    for format in &outputs.formats {
        let choice = app::parse_formatting_string(format.outfmt)?;
        if app::report_format(&choice, app::Program::Tblastx, out)?.is_none() {
            app::formatting_handler_check(&choice, app::Program::Tblastx)?;
            app::xinclude_check_for(&choice, out)?;
        }
    }
    Ok(())
}

/// The search of `run_local` and of the CLI, after NCBI's checks of the options and the
/// reading of the subjects: `Query is Empty!` for records without a query, the query batch
/// size, the subjects without letters, and the query batches.
#[allow(clippy::too_many_arguments)]
fn search(
    args: TblastxArgs,
    query_records: &[FastaRecord],
    reading: QueryReading,
    query_range: Option<SequenceRange>,
    source: SubjectSource<'_>,
    outputs: &mut ReportOutputs<'_>,
    v1: Option<&dyn V1Checks>,
) -> Result<()> {
    // NCBI reads the subjects before the queries (`run`); the records are those of NCBI's
    // reader (`fasta_reader`; `run_local`'s callers read them so). The CLI has already
    // reported an empty query file and read its queries; the records of the other callers
    // get NCBI's warning for an input without records.
    if query_records.is_empty() && matches!(reading, QueryReading::Records) {
        outputs
            .diagnostics
            .write_all(b"Warning: [tblastx] Query is Empty!\n")?;
        return Ok(());
    }
    // NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:136-137
    // ```c
    //         CBlastFastaInputSource fasta(m_CmdLineArgs->GetInputStream(), iconfig);
    //         CBlastInput input(&fasta, m_CmdLineArgs->GetQueryBatchSize());
    // ```
    // The batch size is read before the formatter is made and before the queries, which
    // NCBI reads one batch at a time after the outfmt 0 prolog (`run_in_pool`).
    let batch_size = crate::blastinput::app::query_batch_size("TBLASTX", 10002)?;
    // ABI v1's checks of its `bio` records come where they came before the records of
    // NCBI's reader (plan TD-1).
    if let Some(v1) = v1 {
        v1.check_search(&args, query_range.as_ref(), outputs)?;
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:129-136
    // ```c
    //     if(m_IsDbScan) {
    // 	int num_seqs=0;
    //         int total_length=0;
    // 	if (!is_remote_search)
    //         {
    //                 BlastSeqSrc* seqsrc = db_adapter.MakeSeqSrc();
    //                 num_seqs=BlastSeqSrcGetNumSeqs(seqsrc);
    //                 total_length=static_cast<int>(BlastSeqSrcGetTotLen(seqsrc));
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:773-789
    // ```c
    //         catch(CBlastException & e ) {
    //         	// Skip bad subject sequence
    //         	if(e.GetErrCode() == CBlastException::eInvalidArgument) {
    //         		seqblk_vec->push_back(subj);
    //         ...
    //         		warning += "Subject sequence contains no data";
    //         		ERR_POST(Warning << warning);
    //         		continue;
    // ```
    // The formatter, made after the query batch size and before the outfmt 0 prolog, sets
    // up the subjects (as BLASTN's `search`): a subject without letters (a record without
    // residues, or an interval that starts just past its record's end) gets its warning,
    // stays in the database statistics with no letters, and is never searched
    // (`process_subject` skips subjects shorter than a word).
    // A database has no such warnings (its sequences are not set up by the formatter).
    if let SubjectSource::Fasta { records, .. } = &source {
        crate::blastinput::input_files::write_empty_subject_warnings(
            records,
            "tblastx",
            outputs.diagnostics,
        )?;
    }
    // NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:170-173
    // ```c
    // 	if(UseXInclude(*fmt_args, args[kArgOutput].AsString())) {
    //         	formatter.SetBaseFile(args[kArgOutput].AsString());
    //         }
    //         formatter.PrintProlog();
    // ```
    // The XInclude formats written to standard output fail after the formatter is made
    // (`check_report_formats`: `xinclude_check`), as the other format checks of each
    // requested format.
    check_report_formats(outputs, args.out.as_deref())?;
    // LOSAT's limits come where NCBI starts the search, after its checks.
    check_losat_limits(&args)?;
    crate::blastinput::app::check_unsupported_environment("TBLASTX")?;
    // The query range applies to every query record as NCBI's batch reader reads it
    // (`seq_range::cut_queries`): a record is searched cut, and one whose interval starts
    // just past its end is a query without data (`run_in_pool`).
    let input_records = query_records;
    let (query_records, query_input) = match query_range.as_ref() {
        Some(range) => {
            let ranged = cut_queries(query_records, range);
            (std::borrow::Cow::Owned(ranged.records), ranged.input)
        }
        None => (
            std::borrow::Cow::Borrowed(query_records),
            QueryInput::whole(query_records),
        ),
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1426
    // ```c
    //     itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
    //
    //     /* iterate over all subject sequences */
    //     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
    //            != BLAST_SEQSRC_EOF) {
    //        ...
    //        if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
    // ```
    // `-db`: the database's sequences, read once for every batch (`TblastxDatabase::read_stores`)
    // where NCBI reads them first, in the preliminary search of the first batch with a valid
    // context (a batch without one is not searched, `CheckInternalData`). A database without
    // sequences (BLASTN's `Database memory map file error`) or a sequence file that cannot be
    // read fails that batch, after the reports of the batches before (`run_in_pool`), and
    // LOSAT's refusal of a damaged file stops the run there; a run without such a batch never
    // reads the file.
    let whole_subjects = Placements::default();
    let (stores, stores_error) = match &source {
        SubjectSource::Database(database) => {
            let read = if database.numbers(None).count == 0 {
                Err(crate::algorithm::common::database::PreliminaryReadError(
                    "Database memory map file error",
                )
                .into())
            } else {
                database.read_stores()
            };
            match read {
                Ok(stores) => (stores, None),
                Err(error) => (Vec::new(), Some(error)),
            }
        }
        SubjectSource::Fasta { .. } => (Vec::new(), None),
    };
    let (run_subjects, subject_placements) = match source {
        SubjectSource::Fasta {
            records,
            placements,
        } => (
            RunSubjects {
                subjects: Subjects::Fasta(records),
                numbers: SubjectSetNumbers::of_records(records, args.dbsize),
                shown: ShownSubjects::Fasta(records),
            },
            placements,
        ),
        SubjectSource::Database(database) => (
            RunSubjects {
                subjects: if stores_error.is_some() {
                    Subjects::Unreadable
                } else {
                    Subjects::Stored(&stores)
                },
                numbers: database.numbers(args.dbsize),
                shown: ShownSubjects::Database {
                    database,
                    pairwise: database.pairwise_database(),
                    stores: &stores,
                },
            },
            &whole_subjects,
        ),
    };
    // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
    // TBlastThreads the_threads(GetNumberOfThreads());
    // (*thread)->Run(); (*thread)->Join(&result);
    crate::utils::threading::with_search_pool(args.num_threads, "tblastx", |pool| {
        run_in_pool(
            args,
            &QueryRuns {
                input_records,
                records: &query_records,
                input: &query_input,
                subject_placements,
            },
            reading,
            &run_subjects,
            stores_error,
            outputs,
            pool,
            batch_size,
            v1.is_some(),
        )
    })
}

/// ABI v1's checks of its `bio` records (plan TD-1), which the v1 layer makes on its own
/// records (`web_api::v1_tblastx`): the search calls them where ABI v1 made them when it
/// searched those records (and `check_shown_subject_titles` before the reports), and
/// searches the records of NCBI's reader made from them, which those checks guarantee to be
/// read alike.
pub(crate) trait V1Checks: Sync {
    /// ABI v1's checks where the search starts, in their order: the outfmt 0 and 7 titles,
    /// LOSAT's limits and environment, the residues, the records without residues and the
    /// intervals without letters (the subjects cut to `-subject_loc`, the queries to
    /// `-query_loc`).
    fn check_search(
        &self,
        args: &TblastxArgs,
        query_range: Option<&SequenceRange>,
        outputs: &ReportOutputs<'_>,
    ) -> Result<()>;
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_input_aux.cpp:85-91
// ```c
//     // used for experimentation purposes
//     char* batch_sz_str = getenv("BATCH_SIZE");
//     if (batch_sz_str) {
//         retval = NStr::StringToInt(batch_sz_str);
//         _TRACE("DEBUG: Using query batch size " << retval);
//         return retval;
//     }
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_input_aux.cpp:130-134
// ```c
//     case eTblastx:
//         // N.B.: the splitting is done on the nucleotide query sequences, then
//         // each of these chunks is translated
//         retval = 10002;
//         break;
// ```

// NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:136-137
// ```c
//         CBlastFastaInputSource fasta(m_CmdLineArgs->GetInputStream(), iconfig);
//         CBlastInput input(&fasta, m_CmdLineArgs->GetQueryBatchSize());
// ```
// NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:176-207
// ```c
//         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//
//             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//             CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
//     ...
//                 CLocalBlast lcl_blast(queries, opts_hndl, db_adapter);
//                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
//                 results = lcl_blast.Run();
//     ...
//                 ITERATE(CSearchResultSet, result, *results) {
//                     formatter.PrintOneResultSet(**result, query_batch);
//                 }
//         }
// ```
// Each query batch is searched on its own: the linking cutoffs use the average length
// and the smallest Lambda of the batch's contexts (blast_parameters.c:1023-1026), so the
// batches change the HSPs of a run with several queries. NCBI writes the reports of a
// batch after its search; LOSAT searches the batches one at a time and writes the reports
// of all of them once (`write_tblastx_outputs`), so that what stops the run at a batch
// comes after the reports of the batches before (without the epilog).
#[allow(clippy::too_many_arguments)]
fn run_in_pool(
    args: TblastxArgs,
    query_runs: &QueryRuns<'_>,
    reading: QueryReading,
    subjects: &RunSubjects<'_>,
    stores_error: Option<anyhow::Error>,
    outputs: &mut ReportOutputs<'_>,
    parallel_pool: &crate::utils::threading::SearchPool<'_>,
    batch_size: u32,
    v1_titles: bool,
) -> Result<()> {
    check_report_formats(outputs, args.out.as_deref())?;
    let query_records = query_runs.records;
    let input = query_runs.input;
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:474-504
    // (`s_ReplaceLocalId`)
    // The query ID of the tabular formats is the title's first word, or the local ID
    // (`FastaRecord::shown_id`), as BLASTN's.
    let query_ids: Vec<Arc<[u8]>> = query_records
        .iter()
        .map(|record| Arc::from(record.shown_id()))
        .collect();
    let reports = BatchReports {
        args: &args,
        query_runs,
        query_ids: &query_ids,
        subjects: &subjects.shown,
        numbers: &subjects.numbers,
        titles: std::cell::RefCell::new(DatabaseTitles::new()),
        v1_titles,
        x_search_pool: parallel_pool,
    };
    // NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:173
    // ```c
    //         formatter.PrintProlog();
    // ```
    write_pairwise_prologs(outputs, &args, &subjects.shown, true)?;
    // NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:176-178
    // ```c
    //         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
    //
    //             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
    // ```
    // A reader at the end of its input reads no batch (`!input.End()`): the report of no
    // query, the prolog and the epilog (an empty pipe, or `-query -` after `-subject -`).
    if query_runs.input_records.is_empty() && matches!(reading, QueryReading::End) {
        return reports.write(outputs, SearchedBatches::default(), &[], true);
    }
    // EXPERIMENT (LOSAT_X_TBXPIPE), dispatch: the batches the loop below will search, planned
    // before it (`x_pipe_plan`), when the pipe may run (`x_pipe_wanted`); `None` keeps the loop's
    // own searches.
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let x_plan = x_pipe_wanted(parallel_pool)
        .then(|| x_pipe_plan(input, &reading, batch_size, query_records));
    // NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:193-195
    // ```c
    //                 results = lcl_blast.Run();
    // ```
    // The search of one batch, the queries `range`; with LOSAT_X_TBXPIPE its seed-stage end
    // (`seed_end`) may start the next batch (tblastx/x_batch_pipe.rs).
    let search_subjects = subjects.subjects;
    let numbers = &subjects.numbers;
    let search: &x_batch_pipe::XSearch<'_, TblastxBatch> = &|range, seed_end| {
        search_query_batch(
            &args,
            &query_records[range],
            search_subjects,
            numbers,
            parallel_pool,
            seed_end,
        )
    };
    // The loop's result of the batch of queries `start..end`: the search above, or with
    // LOSAT_X_TBXPIPE the pipe's (`XPipe::take`).
    let batch_loop = |take: &mut dyn FnMut(std::ops::Range<usize>) -> Result<TblastxBatch>| {
        run_batch_loop(
            &args,
            query_runs,
            reading,
            subjects,
            stores_error,
            outputs,
            batch_size,
            &reports,
            take,
        )
    };
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    if let Some(plan) = x_plan {
        if let Some(pipe) = x_batch_pipe::XPipe::new(plan, search_subjects.count(), search) {
            let result = rayon::in_place_scope(|scope| {
                let _closer = pipe.closer();
                batch_loop(&mut |range| pipe.take(scope, range))
            });
            for line in pipe.lines() {
                eprintln!("{line}");
            }
            return result;
        }
    }
    batch_loop(&mut |range| search(range, None))
}

// EXPERIMENT (LOSAT_X_TBXPIPE): whether the pipe may run in this search (tblastx/x_batch_pipe.rs,
// "Off"). No NCBI counterpart.
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
fn x_pipe_wanted(parallel_pool: &crate::utils::threading::SearchPool<'_>) -> bool {
    if !x_batch_pipe::on() {
        return false;
    }
    let wanted = parallel_pool.enabled()
        && parallel_pool.threads() > 1
        && rayon::current_thread_index().is_some()
        && !std::env::vars_os().any(|(name, _)| {
            name.to_str()
                .is_some_and(|name| name.starts_with("LOSAT_DEBUG"))
        })
        && trace_hsp_target().is_none()
        && !diagnostics_enabled()
        && !crate::utils::threading::diagnostics_enabled()
        && !stage_dump::enabled();
    // The batches that start early are marked where they start (`XPipe::seed_ended`).
    if !wanted {
        crate::utils::xreach::x_reach!("LOSAT_X_TBXPIPE", false);
    }
    wanted
}

// EXPERIMENT (LOSAT_X_TBXPIPE): the batches `run_batch_loop` will search unless a search fails, in
// its order: its batch boundaries and the checks of a batch that need no search (the reader's
// error where the reading ends, a batch without a query, a batch whose queries all lack a codon
// with one without letters). The check of the query splitter's environment concerns the first
// batch only, which the loop starts itself after it. No NCBI counterpart (scheduling only; the
// statements are the loop's, see there).
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
fn x_pipe_plan(
    input: &QueryInput,
    reading: &QueryReading,
    batch_size: u32,
    query_records: &[FastaRecord],
) -> Vec<std::ops::Range<usize>> {
    let mut plan = Vec::new();
    let reads_past_records = reading.reads_past_records();
    let reading_fails = matches!(reading, QueryReading::Error { .. });
    let mut reading_ended = false;
    let mut input_start = 0;
    while input_start < input.input_lengths.len() || (reads_past_records && !reading_ended) {
        let (input_end, reached_size) = next_ranged_batch(
            &input.input_lengths,
            &input.skipped,
            input_start,
            batch_size,
        );
        let std::ops::Range { start, end } = input.searched_in(input_start..input_end);
        if input_end == input.input_lengths.len() && !reached_size && reads_past_records {
            reading_ended = true;
            if reading_fails {
                break;
            }
        }
        if start == end {
            break;
        }
        let batch_queries = &query_records[start..end];
        if batch_queries.iter().any(|record| record.seq().is_empty())
            && batch_queries
                .iter()
                .all(|record| record.seq().len() < CODON_LENGTH)
        {
            break;
        }
        plan.push(start..end);
        input_start = input_end;
    }
    plan
}

/// The batch loop of `run_in_pool`: reads the batches, checks them, takes each one's search
/// result from `take`, and writes the reports.
#[allow(clippy::too_many_arguments)]
fn run_batch_loop(
    args: &TblastxArgs,
    query_runs: &QueryRuns<'_>,
    reading: QueryReading,
    subjects: &RunSubjects<'_>,
    mut stores_error: Option<anyhow::Error>,
    outputs: &mut ReportOutputs<'_>,
    batch_size: u32,
    reports: &BatchReports<'_>,
    take: &mut dyn FnMut(std::ops::Range<usize>) -> Result<TblastxBatch>,
) -> Result<()> {
    let query_records = query_runs.records;
    let input = query_runs.input;
    let mut done = SearchedBatches::default();
    // The warnings written before each query's report (`QueryWarnings`).
    let mut warnings: Vec<Vec<u8>> = vec![Vec::new(); query_records.len()];
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1431-1452
    // ```c++
    //     if (m_FormatType == CFormattingArgs::eAsnText
    //       ...
    //       || m_FormatType == CFormattingArgs::eJsonSeqalign
    //       || m_FormatType == CFormattingArgs::eSAM)
    //     {
    //         x_PrintStructuredReport(results, queries);
    //         return;
    //     }
    //     ...
    //     if (results.HasWarnings()) {
    //         ERR_POST(Warning << results.GetWarningStrings());
    //     }
    // ```
    // A structured format (outfmt 5, 8, 9, 12) writes none of a query's own warnings
    // (`ReportFormat::writes_query_warnings`; outfmt 5 puts them in the Iteration's message,
    // `blast_xml`); the warnings go with the first format.
    let result_warnings = outputs
        .formats
        .first()
        .is_none_or(|format| report::output_format(format.outfmt).writes_query_warnings());
    // A batch is read from the input records (`input_start..input_end`); its searched
    // queries are `start..end` (all of them without `-query_loc`). A batch size counts the
    // whole records, and a record whose range starts past its end is skipped
    // (`next_ranged_batch`).
    let mut reading = reading;
    // Whether what ended the reading after the records (`reading`) has come.
    let mut reading_ended = false;
    let mut input_start = 0;
    // The reader is at the end of its input after the last record, unless an `eEOF` or an
    // error ends the reading there; that comes in the batch being read when the records
    // run out before the batch reaches its size, and otherwise in a batch of its own.
    while input_start < input.input_lengths.len()
        || (reading.reads_past_records() && !reading_ended)
    {
        let (input_end, reached_size) = next_ranged_batch(
            &input.input_lengths,
            &input.skipped,
            input_start,
            batch_size,
        );
        let std::ops::Range { start, end } = input.searched_in(input_start..input_end);
        let reading_ends_here =
            input_end == input.input_lengths.len() && !reached_size && reading.reads_past_records();
        // `CFastaReader` writes its messages about a batch's queries (their lines and
        // titles) when it reads them, after the report of the batch before: before the
        // report of the batch's first query (`QueryWarnings`). Skipped records are read too.
        let batch_warnings: Vec<u8> = query_runs.input_records[input_start..input_end]
            .iter()
            .flat_map(|record| record.warnings.iter().copied())
            .collect();
        // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_input.cpp:146-152
        // ```c
        //         try { q.Reset(m_Source->GetNextSequence(scope)); }
        //         catch (const CObjReaderParseException& e) {
        //             if (e.GetErrCode() == CObjReaderParseException::eEOF) {
        //                 break;
        //             }
        //             throw;
        //         }
        // ```
        // An `eEOF` ends the batch with the records read so far. A reader error (or LOSAT's
        // rejection of a Seq-id line) stops the run while the batch is read: after the
        // reports of the batches before (without the epilog) and the messages of the
        // batch's records, which are neither searched nor reported.
        if reading_ends_here {
            reading_ended = true;
            if let QueryReading::Error {
                error,
                warnings: error_warnings,
            } = std::mem::replace(&mut reading, QueryReading::End)
            {
                reports.write_before_failed_batch(outputs, done, &warnings)?;
                outputs.diagnostics.write_all(&batch_warnings)?;
                outputs.diagnostics.write_all(&error_warnings)?;
                return Err(error.into_app_error());
            }
        }
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/objmgr_query_data.cpp:378-380
        // ```c
        //     if (queries.Empty()) {
        //         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
        //     }
        // ```
        // A batch without a query (a batch size of 0, a batch of skipped records only with
        // `-query_loc`, or a batch that an `eEOF` ends before its first record) fails after
        // the reports of the batches before (without the epilog), when NCBI has read that
        // batch (its records' messages).
        if start == end {
            reports.write_before_failed_batch(outputs, done, &warnings)?;
            outputs.diagnostics.write_all(&batch_warnings)?;
            return Err(crate::cli::NativeError {
                exit: 3,
                message: "BLAST engine error: Empty CBlastQueryVector\n".to_string(),
            }
            .into());
        }
        // The query splitter reads its environment variables for every batch, before the
        // queries are set up (`check_query_split_environment`); the first batch meets an
        // error there.
        if input_start == 0 {
            if let Err(error) =
                crate::blastinput::app::check_query_split_environment("TBLASTX", true)?
            {
                reports.write_before_failed_batch(outputs, done, &warnings)?;
                outputs.diagnostics.write_all(&batch_warnings)?;
                return Err(error);
            }
        }
        let batch_queries = &query_records[start..end];
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:632-652
        // ```c
        //         } catch (const CException& e) {
        //             ...
        //             CRef<CSearchMessage> m
        //                 (new CSearchMessage(eBlastSevWarning, index, e.GetMsg()));
        //             messages[index].push_back(m);
        //             s_InvalidateQueryContexts(qinfo, index);
        //         }
        //     ...
        //     // Validate that at least one query context is valid
        //     if (BlastSetup_Validate(qinfo, NULL) != 0 && messages.HasMessages()) {
        //         NCBI_THROW(CBlastException, eSetup, messages.ToString());
        //     }
        // ```
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_aux.cpp:1013-1025
        // ```c
        // TSearchMessages::ToString() const
        // {
        //     string retval;
        //     ITERATE(vector<TQueryMessages>, qm, *this) {
        //         if (qm->empty()) {
        //             continue;
        //         }
        //         ITERATE(TQueryMessages, msg, *qm) {
        //             retval += (*msg)->GetMessage() + " ";
        //         }
        //     }
        //     return retval;
        // }
        // ```
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:923-929
        // ```c
        // BLAST_GetTranslatedProteinLength(size_t nucleotide_length, unsigned int context)
        // {
        //     if (nucleotide_length == 0 || nucleotide_length <= context % CODON_LENGTH) {
        //         return 0;
        //     }
        //     return (nucleotide_length - context % CODON_LENGTH) / CODON_LENGTH;
        // }
        // ```
        // A query's contexts are its six frames, of length 0 (invalid,
        // `s_QueryInfo_SetContext`) below one codon. A query without letters
        // (`Sequence contains no data`, blast_setup.hpp:190-198) gets a message and its
        // contexts are invalidated; a query of one or two letters gets none. A batch whose
        // queries all lack a valid context stops NCBI only when one of them has a message:
        // one message per query without letters (`GetMessage` puts the severity first),
        // after the reports of the batches before (`CATCH_ALL`: `BLAST engine error: `,
        // exit 3). Any other batch is set up.
        let no_data_queries = batch_queries
            .iter()
            .filter(|record| record.seq().is_empty())
            .count();
        if no_data_queries > 0
            && batch_queries
                .iter()
                .all(|record| record.seq().len() < CODON_LENGTH)
        {
            reports.write_before_failed_batch(outputs, done, &warnings)?;
            outputs.diagnostics.write_all(&batch_warnings)?;
            let mut message = String::from("BLAST engine error: ");
            for _ in 0..no_data_queries {
                message.push_str("Warning: ");
                message.push_str(NO_DATA_MESSAGE);
                message.push(' ');
            }
            message.push('\n');
            return Err(crate::cli::NativeError { exit: 3, message }.into());
        }
        warnings[start].extend_from_slice(&batch_warnings);
        let batch = match take(start..end) {
            Ok(batch) => {
                // LOSAT_TIMING: the batch's block, kept in its result by a search of the pipe
                // (LOSAT_X_TBXPIPE), in the batch order.
                for line in &batch.x_timing {
                    eprintln!("{line}");
                }
                batch
            }
            Err(error) if error.is::<SubjectsUnreadable>() => {
                // The first batch with a valid context of a `-db` search whose sequences could
                // not be read (`search`): NCBI's read error is an engine error after the
                // reports of the batches before; LOSAT's refusal (a truncated sequence file,
                // unreadable ambiguity data) stops the run, as any other refusal does.
                let error = stores_error
                    .take()
                    .expect("unreadable subjects keep their read error");
                let read_error =
                    error.downcast::<crate::algorithm::common::database::PreliminaryReadError>()?;
                reports.write_before_failed_batch(outputs, done, &warnings)?;
                outputs.diagnostics.write_all(&batch_warnings)?;
                return Err(crate::blastinput::app::engine_error(read_error.0));
            }
            Err(error) => match error.downcast::<ErrorAfterBatchReports>() {
                // NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:203-209
                // ```c
                //                 ITERATE(CSearchResultSet, result, *results) {
                //                     formatter.PrintOneResultSet(**result, query_batch);
                //                 }
                //             }
                //         }
                //
                //         formatter.PrintEpilog(opt);
                // ```
                // The reports of the batches before are written (for the first batch, the
                // outfmt 0 prolog), then the messages of the failing batch's records (read
                // before its search), then the error, which skips `PrintEpilog`.
                Ok(ErrorAfterBatchReports(error)) => {
                    reports.write_before_failed_batch(outputs, done, &warnings)?;
                    outputs.diagnostics.write_all(&batch_warnings)?;
                    return Err(error);
                }
                Err(error) => return Err(error),
            },
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:283-293
        // ```c++
        //     m_TbackSearch.Reset(new CBlastTracebackSearch(m_QueryFactory,
        //                                                   m_InternalData,
        //                                                   m_Opts,
        //                                                   seqinfo_src,
        //                                                   search_msgs));
        // ```
        // `-db`: the definition lines of the batch's reported sequences are read where NCBI
        // makes the batch's Seq-aligns (`database::read_titles`, in report order); a header file
        // that NCBI cannot read fails this batch, after the reports of the batches before.
        if let ShownSubjects::Database { database, .. } = &subjects.shown {
            let ordered = report_order(&args, batch_queries, batch.hits.clone());
            let read = crate::algorithm::tblastx::database::read_titles(
                database,
                &ordered,
                start,
                &mut reports.titles.borrow_mut(),
            );
            if let Err(error) = read {
                match error.downcast::<crate::algorithm::common::database::HeaderReadError>() {
                    Ok(failure) => {
                        reports.write_before_failed_batch(outputs, done, &warnings)?;
                        outputs.diagnostics.write_all(&batch_warnings)?;
                        return Err(failure.error);
                    }
                    Err(error) => return Err(error),
                }
            }
        }
        // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:1450-1452
        // ```c
        //     if (results.HasWarnings()) {
        //         ERR_POST(Warning << results.GetWarningStrings());
        //     }
        // ```
        // NCBI reference: c++/src/algo/blast/core/blast_stat.c:2815-2821
        // ```c
        //    if (valid_context == FALSE)
        //    {   /* No valid contexts were found. */
        //        /* Message for non-translated search issued above. */
        //        if (Blast_QueryIsTranslated(program) ) {
        //             Blast_MessageWrite(blast_message, eBlastSevWarning, kBlastMessageNoContext,
        //             kBlastErrMsg_CantCalculateUngappedKAParams);
        //        }
        // ```
        // NCBI reference: c++/src/algo/blast/api/local_blast.cpp:197-201,221
        // ```c
        //          for (index=0; index<local_query_data->GetNumQueries(); index++)
        //          {
        //               CConstRef<objects::CSeq_id> query_id(local_query_data->GetSeq_loc(index)->GetId());
        //               TQueryMessages q_msg;
        //               local_query_data->GetQueryMessages(index, q_msg);
        //     ...
        //          msg_vec.Combine(m_PrelimSearch->GetSearchMessages());
        // ```
        // The warnings of a query come with its report (`PrintOneResultSet`), before its
        // preamble (`QueryWarnings`): a query without letters has the message of its set-up
        // (above); a translated query has no Karlin-Altschul message of its own, and the
        // queries of a batch that is not searched (no valid context) all get the message of
        // the batch.
        for (offset, query) in batch_queries.iter().enumerate().filter(|_| result_warnings) {
            let index = start + offset;
            let mut messages = Vec::new();
            if !batch.searched {
                messages.push(crate::report::query_warnings::INVALID_QUERY_MESSAGE.to_string());
            }
            if query.seq().is_empty() {
                messages.push(NO_DATA_MESSAGE.to_string());
            }
            warnings[index].extend(crate::report::query_warnings::query_warning(
                "tblastx",
                input.ordinal(index),
                query,
                &messages,
            ));
        }
        // The batch numbers its queries from 0.
        done.hits.extend(batch.hits.into_iter().map(|mut hsp| {
            hsp.hit.q_idx += start as u32;
            hsp
        }));
        done.query_stats.extend(batch.queries);
        done.searched
            .extend(std::iter::repeat_n(batch.searched, end - start));
        input_start = input_end;
    }
    reports.write(outputs, done, &warnings, true)
}

/// NCBI's `CODON_LENGTH`: a query shorter than one codon has no frame of positive length.
///
/// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_def.h:62-64
/// ```c
/// #ifndef CODON_LENGTH
/// #define CODON_LENGTH 3
/// #endif
/// ```
const CODON_LENGTH: usize = 3;
